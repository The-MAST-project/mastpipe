"""The ``identify`` command: the attended bootstrap of a configuration's first reference
solution, and the seeded re-solve of a full arc behind ``solve``.

A reference solution is made once per configuration and layout preset, from a full arc
that lights every fiber. The chain from the frame to its measured lines is shared
(:func:`_arc_lines`): the arc is loaded and calibrated with the bias and background
masters found in ``out_dir/calib`` (nothing but the variance is applied when there are
none), the traces come from a flat master there when one exists and from the arc itself
otherwise (percentile collapse), the layout is matched, the spectra are extracted and
continuum-fitted and the lines are found. Then:

- :func:`run_identify` **without anchors** draws the anchor panels
  (:func:`mastpipe.io.qa.plot_anchor_panels`) and writes the numbered candidates, so that
  the user can read a few template-measured pairs off the figure and run again with them;
  **with anchors** (given, or the instrument's stored pairs for the configuration and
  preset) it builds the anchor model of the anchor trace, solves every fiber from it and
  writes the reference product with the fiber-to-fiber agreement table;
- :func:`run_solve` is the same chain seeded by a previous
  :class:`~mastpipe.core.model.ReferenceSolution` (the caller finds it).

The template of the panels is the Murphy et al. atlas when ``atlas_dir`` is given and the
synthetic atlas of the shipped line table otherwise. This is pipeline-layer code: it reads
files, holds the :class:`~mastpipe.config.models.EffectiveConfig` and writes products.
"""

from __future__ import annotations

import argparse
import glob
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from mastpipe.__version__ import __version__
from mastpipe.config.loader import build_effective, config_hash, load_yaml, parse_override
from mastpipe.config.models import EffectiveConfig
from mastpipe.core.model import (
    CalibrationSet,
    Frame,
    FrameRef,
    MasterCalib,
    ReferenceSolution,
    TraceSet,
)
from mastpipe.instruments.registry import (
    load_anchor_trace,
    load_anchors,
    load_header_map,
    load_instrument,
    load_layout,
    load_line_table,
)
from mastpipe.io import products
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.io.qa import plot_anchor_panels
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.atlas import (
    catalogue_arrays,
    convolve_to_resolving_power,
    flag_usable_lines,
    load_murphy_atlas,
    remove_pedestal,
    synthetic_atlas,
)
from mastpipe.steps.calib import calibrate_frame, calibration_steps
from mastpipe.steps.extract import continuum_subtract_all, extract_aperture
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.lines import find_lines
from mastpipe.steps.trace import find_traces, layout_indexed, match_layout, trace_source
from mastpipe.steps.wavecal import build_reference, fiber_agreement, solve_full_arc
from mastpipe.steps.wavemodel import anchor_model

CALIB_FOLDER = "calib"
"""Sub-folder of the output folder holding the masters read and the reference written."""
PANELS_FILE = "anchor_panels.png"
CANDIDATES_FILE = "anchor_candidates.csv"
AGREEMENT_FILE = "fiber_agreement.csv"
ATLAS_DIR_VARIABLE = "MASTPIPE_ATLAS_DIR"
PAD_NM = 0.5
"""Band padding of the atlas spectrum, so a line at the edge keeps its wings."""


# ------------------------------------------------------------------------ frame to lines


@dataclass(eq=False)
class ArcLines:
    """One full arc taken to its measured lines by :func:`_arc_lines`.

    ``frame`` is the calibrated arc, ``traceset`` the matched traces, ``fluxes`` the
    continuum-subtracted flux of every trace (``flux - continuum``), ``tables`` the line
    tables of :func:`~mastpipe.steps.lines.find_lines`, ``calibrations`` the file names of
    the masters used and ``steps`` the operations applied, for provenance.
    """

    frame: Frame
    traceset: TraceSet
    fluxes: list[np.ndarray]
    tables: list[pd.DataFrame]
    calibrations: list[str]
    steps: list[str]


MASTER_KINDS = ("bias", "background", "flat")
"""Master kinds :func:`_arc_lines` looks for in ``out_dir/calib``."""
NO_TIME = "notime"
"""The stamp :func:`~mastpipe.io.products.product_name` writes for a frame without a mid-time."""


def _latest_master(calib_dir: Path, effective: EffectiveConfig, kind: str) -> Path | None:
    """The latest ``{instrument}_{configuration}_{camera}_{stamp}_{kind}.fits`` of ``calib_dir``.

    Only the masters of ``effective``'s instrument, configuration and camera qualify; the
    latest is the one with the largest ``YYYYmmddTHHMMSS`` stamp (compared as a stamp, not
    as a file name), and a ``notime`` master is taken only when no dated one exists.
    ``None`` without a candidate.
    """
    prefix = f"{effective.instrument}_{effective.configuration}_{effective.camera}_"
    suffix = f"_{kind}.fits"
    dated: list[tuple[str, Path]] = []
    undated: list[Path] = []
    for path in calib_dir.glob(f"{glob.escape(prefix)}*{glob.escape(suffix)}"):
        if not path.is_file():
            continue
        stamp = path.name[len(prefix) : -len(suffix)]
        if stamp == NO_TIME:
            undated.append(path)
        elif (
            stamp.isascii()
            and len(stamp) == 15
            and stamp[8] == "T"
            and stamp.replace("T", "").isdigit()
        ):
            dated.append((stamp, path))
    if dated:
        return max(dated, key=lambda item: (item[0], item[1].name))[1]
    return sorted(undated)[-1] if undated else None


def _master_name(kind: str, master: MasterCalib) -> str:
    """Identifier of an injected master: its provenance ``name`` when the caller set one,
    else ``{kind}:{validity mid_time | notime}``."""
    name = master.provenance.get("name")
    if name:
        return str(name)
    return f"{kind}:{master.validity.get('mid_time') or NO_TIME}"


def _frame_ref(full_arc: Path, instrument: str):
    pattern = glob.escape(full_arc.name)
    refs = scan_folder(full_arc.parent, instrument, load_header_map(instrument), pattern=pattern)
    if not refs:
        raise ValueError(f"{full_arc}: no image HDU to reduce")
    return refs[0]


def _arc_lines(
    full_arc: Path | FrameRef,
    effective: EffectiveConfig,
    layout: pd.DataFrame,
    out_dir: Path,
    masters: dict[str, MasterCalib] | None = None,
    traceset: TraceSet | None = None,
) -> ArcLines:
    """Load, calibrate, trace, extract, continuum-fit and find the lines of a full arc.

    ``full_arc`` is the frame's path (scanned on its own through the instrument's header
    map) or a :class:`~mastpipe.core.model.FrameRef` a runner has already scanned, whose
    metadata (a manifest's overrides included) is then used as given. Masters come from
    ``masters`` (keyed by kind, as a runner holding its library's choice passes them) or,
    when that is ``None``, from ``out_dir/calib``: the latest
    ``*_bias.fits``, ``*_background.fits`` and ``*_flat.fits`` of the effective instrument,
    configuration and camera (:func:`_latest_master`). The bias and background calibrate the
    frame, the background only when ``params.calib.subtract_background`` is on (otherwise it
    serves the saturation mask alone, which the steps record and ``calibrations`` does not);
    without them the frame only gains its variance and mask. ``traceset`` is used as given
    (a runner's night traces); otherwise the traces are found on the flat master when there
    is one and on the arc itself with ``collapse="percentile"`` when not, matched to the
    layout with ``params.trace.layout_tolerance`` and keyed by the layout's trace id
    (:func:`~mastpipe.steps.trace.layout_indexed`: detections matching no layout row are
    dropped, so ``traceset`` holds layout fibers only). Returns an :class:`ArcLines`; its
    ``calibrations`` are the masters that were applied (file names when read from
    ``out_dir/calib``, :func:`_master_name` when injected).
    """
    out_dir = Path(out_dir)
    params = effective.params
    detector = effective.detector
    if isinstance(full_arc, FrameRef):
        ref = full_arc
    else:
        ref = _frame_ref(Path(full_arc), effective.instrument)
    frame = load_frame(ref, detector)
    shape = (int(frame.shape[0]), int(frame.shape[1]))

    found: dict[str, tuple[MasterCalib, str]] = {}
    if masters is not None:
        found = {kind: (master, _master_name(kind, master)) for kind, master in masters.items()}
    else:
        calib_dir = out_dir / CALIB_FOLDER
        if calib_dir.is_dir():
            for kind in MASTER_KINDS:
                path = _latest_master(calib_dir, effective, kind)
                if path is not None:
                    found[kind] = (products.read_master(path)[0], path.name)
    bias = found.get("bias")
    background = found.get("background")
    flat = found.get("flat")
    calibrations: list[str] = []
    if bias is not None:
        calibrations.append(bias[1])
    t_background = None if background is None else background[0].validity.get("exposure_time")
    calib = CalibrationSet(
        bias=None if bias is None else bias[0].data,
        dark_rate=0.0,
        background=None if background is None else background[0].data,
        t_background=None if t_background is None else float(t_background),
        gain=detector.gain,
        read_noise=detector.read_noise,
        shape=shape,
    )
    subtract_background = background is not None and params.calib.subtract_background
    arc = calibrate_frame(frame, calib, subtract_background, saturation=detector.saturation)
    steps = calibration_steps(calib, float(frame.meta.exposure_time or 0.0), subtract_background)
    if background is not None:
        if subtract_background:
            calibrations.append(background[1])
        else:
            steps.append(
                f"background {background[1]} present but not subtracted "
                "(calib.subtract_background off): saturation mask only"
            )

    if traceset is not None:
        source = "given"
    else:
        flat_frames: list[Frame] = []
        if flat is not None:
            flat_frames.append(Frame(data=flat[0].data, meta=frame.meta, variance=flat[0].variance))
            calibrations.append(flat[1])
        source, trace_frame = trace_source({"flat": flat_frames, "arc_full": [arc]})
        trace_params = params.trace
        if source == "arc_full":
            trace_params = trace_params.model_copy(update={"collapse": "percentile"})
        assert trace_frame is not None
        matched = match_layout(
            find_traces(trace_frame.data, trace_params), layout, params.trace.layout_tolerance
        )
        traceset = layout_indexed(matched)
    if len(traceset) == 0:
        raise ValueError(f"{ref.meta.source}: no trace found on the {source.replace('_', ' ')}")
    dropped = f", {len(traceset.extra)} extra dropped" if traceset.extra else ""
    steps.append(
        f"traces {source} ({len(traceset)} found, {len(traceset.missing)} missing{dropped})"
    )

    spectra = extract_aperture(arc, traceset, params.extract)
    spectra = continuum_subtract_all(spectra, params.extract)
    fluxes = [spectrum.flux - spectrum.continuum for spectrum in spectra]
    tables, _ = find_lines(arc.data, arc.variance, traceset, params.lines)
    steps += ["extract", "continuum", "lines"]
    return ArcLines(
        frame=arc,
        traceset=traceset,
        fluxes=fluxes,
        tables=tables,
        calibrations=calibrations,
        steps=steps,
    )


# -------------------------------------------------------------------------- the template


def _template(effective: EffectiveConfig, atlas_dir: Path | None):
    """``(wave, flux_conv, usable)`` of the anchor panels' template.

    The Murphy et al. atlas of ``atlas_dir`` over the band, pedestal removed and convolved
    to the configuration's resolving power, with
    :func:`~mastpipe.steps.atlas.flag_usable_lines` at ``params.wavecal.max_shift_fwhm`` /
    ``min_snr_atlas``; or, without an atlas folder, the synthetic atlas of the shipped
    line table convolved the same way, with the table's own ``usable`` flags.
    """
    band = effective.band_nm
    R = effective.resolving_power
    if atlas_dir is not None:
        spectrum, catalogue = load_murphy_atlas(Path(atlas_dir), band, pad=PAD_NM)
        wave = spectrum["wavelength"].to_numpy(dtype=float)
        flux = remove_pedestal(spectrum["flux"].to_numpy(dtype=float))
        conv, fwhm = convolve_to_resolving_power(wave, flux, R)
        wavecal = effective.params.wavecal
        usable = flag_usable_lines(
            wave, conv, catalogue, fwhm, wavecal.max_shift_fwhm, wavecal.min_snr_atlas
        )
        return wave, conv, usable
    table = load_line_table(effective.instrument, effective.configuration)
    catalogue = pd.DataFrame({"wavelength": table["wavelength_nm"], "logI": table["logI"]})
    wave, flux = synthetic_atlas(catalogue, (band[0] - PAD_NM, band[1] + PAD_NM), R)
    conv, _ = convolve_to_resolving_power(wave, flux, R)
    usable = pd.DataFrame(
        {
            "wavelength": table["wavelength_nm"],
            "species": table["species"],
            "usable": table["usable"].to_numpy(dtype=bool),
        }
    )
    return wave, conv, usable


def _write_candidates(path: Path, template: pd.DataFrame, measured: pd.DataFrame) -> None:
    """``id, wavelength_nm, x_center``: ``template_<id>`` rows carry the wavelength,
    ``measured_<id>`` rows the pixel; the other cell is empty."""
    rows = pd.concat(
        [
            pd.DataFrame(
                {
                    "id": [f"template_{k}" for k in template["id"]],
                    "wavelength_nm": template["wavelength"].to_numpy(dtype=float),
                    "x_center": np.nan,
                }
            ),
            pd.DataFrame(
                {
                    "id": [f"measured_{k}" for k in measured["id"]],
                    "wavelength_nm": np.nan,
                    "x_center": measured["x_center"].to_numpy(dtype=float),
                }
            ),
        ],
        ignore_index=True,
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        rows.to_csv(handle, index=False)


# ----------------------------------------------------------------------- the reference


def _provenance(
    effective: EffectiveConfig, lines: ArcLines, steps: list[str], layout: str | None = None
) -> Provenance:
    return Provenance(
        pipeline_version=__version__,
        config_hash=config_hash(effective),
        inputs=[lines.frame.meta.checksum],
        calibrations=list(lines.calibrations),
        wavelength_status="solved",
        steps=[*lines.steps, *steps],
        layout=layout,
        timezone=lines.frame.meta.timezone,
    )


def _solve_and_package(
    lines: ArcLines,
    effective: EffectiveConfig,
    preset: str,
    start: Any,
    ref_fiber: int,
) -> ReferenceSolution:
    """:func:`~mastpipe.steps.wavecal.solve_full_arc` on every trace, packaged as a
    reference created at the arc's mid-time; ``ValueError`` when no fiber solves or the
    frame has no mid-time."""
    geom = grating_geometry(effective.geometry, effective.detector.pixel_size_um)
    table = load_line_table(effective.instrument, effective.configuration)
    cat_wl, cat_sp = catalogue_arrays(table, effective.params.wavecal.usable_only)
    solutions = solve_full_arc(
        lines.tables,
        lines.fluxes,
        lines.traceset,
        start,
        cat_wl,
        cat_sp,
        effective.params.wavecal,
        ref_fiber=ref_fiber,
        geom=geom,
    )
    meta = lines.frame.meta
    if not solutions:
        raise ValueError(f"{meta.source}: no fiber solved from the starting model")
    if meta.mid_time is None:
        raise ValueError(f"{meta.source}: the frame has no mid-time to date the reference")
    return build_reference(
        effective.instrument,
        effective.configuration,
        preset,
        solutions,
        created=meta.mid_time,
        source_checksums=(meta.checksum,),
    )


def _write_reference_product(
    out_dir: Path, ref: ReferenceSolution, prov: Provenance, write: bool, camera: str
) -> Path:
    """The reference at ``out_dir/calib/<identifier>_wref.fits``, registered in the
    :class:`~mastpipe.pipeline.callib.CalibLibrary` of that folder under ``camera``, with
    ``fiber_agreement.csv`` next to it; only the path when ``write`` is false."""
    path = Path(out_dir) / CALIB_FOLDER / products.reference_name(ref)
    if not write:
        return path
    CalibLibrary(path.parent).add_reference(ref, prov, camera)
    with (path.parent / AGREEMENT_FILE).open("w", encoding="utf-8", newline="") as handle:
        fiber_agreement(ref.solutions).to_csv(handle, index=False)
    return path


# ----------------------------------------------------------------------------- commands


def run_identify(
    full_arc: Path,
    effective: EffectiveConfig,
    layout: str,
    out_dir: Path,
    anchors: list[tuple[float, float]] | None = None,
    anchor_trace: int | None = None,
    atlas_dir: Path | None = None,
    write_reference: bool = True,
) -> tuple[ReferenceSolution | None, Path]:
    """Bootstrap the reference solution of ``effective``'s configuration from a full arc.

    ``layout`` is the fiber layout preset name (the table is loaded from the instrument
    definition; the name also selects the stored anchors and names the reference's
    preset). ``anchors`` are ``(wavelength nm, pixel)`` pairs measured on trace
    ``anchor_trace``, the layout's trace id of that fiber (``0`` by default); when ``None``
    the instrument's stored pairs for the configuration and preset are used, and
    ``anchor_trace`` defaults to the trace they were measured on.

    Without any anchors the two-panel figure is written to ``out_dir/anchor_panels.png``
    (template from the atlas in ``atlas_dir`` when given, else synthetic) and the
    numbered candidates to ``out_dir/anchor_candidates.csv``; returns ``(None, png path)``.

    With anchors the anchor model seeds :func:`~mastpipe.steps.wavecal.solve_full_arc`,
    and the reference is written to ``out_dir/calib/<identifier>_wref.fits``, registered in
    that folder's calibration library, with ``fiber_agreement.csv`` beside it
    (``write_reference=False`` skips the writing);
    returns ``(reference, reference path)``. ``ValueError`` when ``anchor_trace`` is not
    among the matched traces or no fiber solves.
    """
    full_arc = Path(full_arc)
    out_dir = Path(out_dir)
    instrument, configuration = effective.instrument, effective.configuration
    if anchors is None:
        anchors = load_anchors(instrument, configuration, layout)
        if anchor_trace is None:
            anchor_trace = load_anchor_trace(instrument, configuration, layout)
    trace = 0 if anchor_trace is None else int(anchor_trace)
    lines = _arc_lines(full_arc, effective, load_layout(instrument, layout), out_dir)
    positions = {int(t.index): k for k, t in enumerate(lines.traceset.traces)}
    if trace not in positions:
        raise ValueError(
            f"anchor trace {trace} is not among the matched traces {sorted(positions)}"
        )
    at = positions[trace]

    if not anchors:
        wave, conv, usable = _template(effective, atlas_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        png = out_dir / PANELS_FILE
        template, measured = plot_anchor_panels(
            wave,
            conv,
            usable,
            lines.fluxes[at],
            lines.tables[at],
            effective.band_nm,
            flip_x=effective.geometry.dispersion_sign < 0,
            out=png,
        )
        _write_candidates(out_dir / CANDIDATES_FILE, template, measured)
        return None, png

    geom = grating_geometry(effective.geometry, effective.detector.pixel_size_um)
    start = anchor_model(list(anchors), geom, int(lines.frame.shape[1]))
    ref = _solve_and_package(lines, effective, layout, start, trace)
    pairs = ", ".join(f"{wl:.4f} nm at {px:.2f} px" for wl, px in anchors)
    prov = _provenance(
        effective, lines, [f"wavecal from anchors on trace {trace}: {pairs}"], layout=layout
    )
    return ref, _write_reference_product(out_dir, ref, prov, write_reference, effective.camera)


def run_solve(
    full_arc: Path,
    effective: EffectiveConfig,
    layout: str,
    out_dir: Path,
    previous: ReferenceSolution,
) -> ReferenceSolution:
    """Solve every fiber of a full arc seeded by ``previous`` and write the new reference.

    The same chain as :func:`run_identify`, with :func:`~mastpipe.steps.wavecal.solve_full_arc`
    starting every fiber from its own earlier solution. ``previous`` must be a reference of
    the same instrument, configuration and layout preset (``ValueError`` otherwise); its
    identifier and hash are recorded as the new product's ``reference_id`` /
    ``reference_hash`` provenance, naming the seed. Returns the new reference, written to
    ``out_dir/calib/<identifier>_wref.fits`` (registered in that folder's calibration
    library) with ``fiber_agreement.csv``.
    """
    full_arc = Path(full_arc)
    out_dir = Path(out_dir)
    expected = (effective.instrument, effective.configuration, layout)
    actual = (previous.instrument, previous.configuration, previous.fiber_preset)
    if actual != expected:
        raise ValueError(
            f"reference {previous.identifier} is for instrument/configuration/preset "
            f"{actual}, the arc is reduced as {expected}"
        )
    lines = _arc_lines(full_arc, effective, load_layout(effective.instrument, layout), out_dir)
    ref = _solve_and_package(lines, effective, layout, previous, 0)
    prov = _provenance(
        effective, lines, [f"wavecal seeded by reference {previous.identifier}"], layout=layout
    )
    prov.reference_id = previous.identifier
    prov.reference_hash = previous.hash or None
    _write_reference_product(out_dir, ref, prov, True, effective.camera)
    return ref


# ------------------------------------------------------------------------ command line


def parse_anchor_pairs(text: str) -> list[tuple[float, float]]:
    """``"515.8605:267.22,525.8360:1423.22"`` to ``[(515.8605, 267.22), (525.836, 1423.22)]``.

    Pairs are comma-separated ``wavelength_nm:pixel``; blanks around them are ignored.
    ``ValueError`` for an empty string or a malformed pair.
    """
    pairs: list[tuple[float, float]] = []
    for item in text.split(","):
        parts = item.strip().split(":")
        if len(parts) != 2:
            raise ValueError(f"anchor pair {item.strip()!r} must be wavelength_nm:pixel")
        try:
            pairs.append((float(parts[0]), float(parts[1])))
        except ValueError as error:
            raise ValueError(f"anchor pair {item.strip()!r} must be wavelength_nm:pixel") from error
    return pairs


def read_anchors_file(path: Path) -> tuple[list[tuple[float, float]], int | None]:
    """``(pairs, trace)`` of a YAML anchors file: ``pairs: [[nm, px], ...]`` and an optional
    ``trace:``; ``trace`` is ``None`` when the file does not name one."""
    data = load_yaml(path)
    raw = data.get("pairs") or []
    pairs: list[tuple[float, float]] = []
    for pair in raw:
        if not isinstance(pair, list | tuple) or len(pair) != 2:
            raise ValueError(f"{Path(path).name}: each anchor pair must be [wavelength_nm, pixel]")
        pairs.append((float(pair[0]), float(pair[1])))
    trace = data.get("trace")
    return pairs, None if trace is None else int(trace)


def register_identify(subparsers: argparse._SubParsersAction) -> None:
    """Add the ``identify`` sub-command."""
    parser = subparsers.add_parser(
        "identify",
        help="bootstrap a reference wavelength solution from a full arc",
        description=(
            "Reduce one full arc to its lines and either draw the anchor panels (no anchors) "
            "or solve every fiber from the anchor pairs and write the reference solution."
        ),
    )
    parser.add_argument("full_arc", type=Path, help="the full-arc FITS file")
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    parser.add_argument("--instrument", required=True, help="instrument name, e.g. highspec")
    parser.add_argument("--configuration", required=True, help="configuration, e.g. Mg")
    parser.add_argument("--layout", help="fiber layout preset (default: the instrument's)")
    parser.add_argument("--camera", help="camera name (default: the instrument's only camera)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--anchors",
        metavar="NM:PX,NM:PX",
        help="anchor pairs wavelength_nm:pixel, comma-separated: 515.8605:267.22,525.836:1423.22",
    )
    group.add_argument(
        "--anchors-file",
        type=Path,
        metavar="FILE",
        help="YAML with pairs: [[nm, px], ...] and trace:",
    )
    group.add_argument(
        "--stored-anchors",
        action="store_true",
        help="use the anchor pairs stored with the instrument for this configuration and layout",
    )
    parser.add_argument(
        "--trace",
        type=int,
        help="layout trace id the anchors were measured on (default: the file's, else 0)",
    )
    parser.add_argument(
        "--atlas-dir",
        type=Path,
        default=_atlas_dir_default(),
        help=f"folder of the Th-Ar atlas for the panels (default: ${ATLAS_DIR_VARIABLE})",
    )
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        metavar="KEY.PATH=VALUE",
        help="override a step parameter, e.g. trace.trace_halfwidth=4 (repeatable)",
    )
    parser.set_defaults(func=identify_command)


def _atlas_dir_default() -> Path | None:
    value = os.environ.get(ATLAS_DIR_VARIABLE)
    return Path(value) if value else None


def identify_command(args: argparse.Namespace) -> int:
    """Run ``identify`` from parsed arguments; prints the panel or reference path."""
    instrument = load_instrument(args.instrument)
    camera = args.camera
    if camera is None:
        if len(instrument.cameras) != 1:
            raise SystemExit(
                f"{args.instrument} has cameras {sorted(instrument.cameras)}; pass --camera"
            )
        camera = next(iter(instrument.cameras))
    layout = args.layout or instrument.default_layout
    overrides = [parse_override(text) for text in args.set]
    effective = build_effective(instrument, args.configuration, camera, overrides, layout=layout)

    anchors: list[tuple[float, float]] | None = []
    trace = args.trace
    if args.anchors:
        anchors = parse_anchor_pairs(args.anchors)
    elif args.anchors_file is not None:
        anchors, file_trace = read_anchors_file(args.anchors_file)
        if trace is None:
            trace = file_trace
    elif args.stored_anchors:
        anchors = None  # run_identify looks them up
    ref, path = run_identify(
        args.full_arc,
        effective,
        layout,
        args.out,
        anchors=anchors,
        anchor_trace=trace,
        atlas_dir=args.atlas_dir,
    )
    if ref is None:
        print(f"anchor panels: {path}")
        print(f"candidates: {Path(args.out) / CANDIDATES_FILE}")
    else:
        print(f"reference {ref.identifier}: {len(ref.solutions)} fibers solved")
        print(f"written: {path}")
    return 0
