"""The quick-look recipe: a night folder to masters, traces, references, drift models,
wavelength-calibrated and sky-subtracted spectra and per-target combinations, with QC.

:func:`run_quicklook` scans the folder and its subfolders (manifest applied; frames of
another instrument and manifest rows matching no frame go into ``skipped``), drops the
frames whose type is ``unknown`` into ``skipped`` with their classification rule, and
reduces every ``(camera, configuration)`` group in turn:

1. **masters** — bias, dark and background frames are stacked (dark and background one
   exposure time at a time, the largest batch serving the night), the flats calibrated
   and stacked, and every master is added to the calibration library (``out/calib`` or
   ``options.library``) — unless the library already holds one stacked from the same
   frames under the same stacking parameters and detector (``calib_hash``; and, for the
   flat, calibrated with the same masters), which is reused rather than rebuilt and
   overwritten; ``options.force`` rebuilds. The frames are loaded one at a time: one that
   cannot be read, is not binned as the group's science frames (the
   most common binning among them) or has another shape than the detector's (after crop
   and orientation) is recorded in ``failed`` and left out. A kind the night lacks, or
   none of whose frames survive, is taken from the library when it holds one that
   applies;

   1b. **light frames** — a frame whose header type is the provisional ``light`` (what the
   control software writes for arcs and science frames alike) is calibrated with the
   masters (no cosmic-ray flagging) and resolved from its content
   (:func:`~mastpipe.steps.classify.classify_light`) into an arc or a science frame, its
   ``classification_rule`` set to ``content:arc_like`` or ``content:continuum``, or
   skipped with ``content:no_signal``. This comes before the traces because on a night
   without a flat the full arc is itself a ``light`` frame and is the trace source. A
   manifest row's ``exposure_type`` wins: such a frame is never ``light``. The rule is in
   the QC of every arc and science frame and, when it came from the content, in the
   product's provenance steps;

2. **traces** — found on the flat master, or on a full arc with the percentile collapse
   when there is no flat, matched to the layout preset (``options.layout``, else the
   frames' ``fiber_preset`` from a manifest, else the instrument's default) and keyed by
   the layout's trace id (detections matching no layout row are dropped, logged and
   listed in the trace QC as ``extra``), then written with a QC file; the relative
   throughput comes from the flat (``{}`` without one). Every per-fiber key downstream
   (spectra, reference solutions, drift models, throughput) is that layout id, so a fiber
   missing on one night is simply absent rather than shifting its neighbours' identities;
3. **arcs** — classified by their lit traces when the header does not say which kind they
   are; a **full arc** is re-solved seeded by the latest reference of the configuration
   and preset and added to the library (skipped with "run identify" when there is none,
   and skipped when the library already holds a reference solved from that very frame);
   a **calibration-fiber arc** is taken to its lines on the arc fibers and fitted for its
   drift against the reference valid at its time, kept in memory (tagged with that
   reference, so that only arcs fitted against a frame's own reference are ever paired
   with it) and written as a QC file;
4. **science frames** — calibrated with cosmic rays flagged, extracted, given wavelengths
   by the time policy (``transferred`` through the paired drift arcs, ``reference`` from
   the reference alone, ``nominal`` from the configuration's geometry when there is no
   reference), sky-subtracted, passed through the flux and telluric stages, written with
   their provenance and QC, and combined per target. A fiber whose trace centre sits
   farther than ``trace.layout_tolerance`` from the reference solution's is not given a
   wavelength (a warning, and ``wavelength_skipped`` in the QC); a reference on another
   pixel grid than the frame's is refused with ``ValueError``.

A science product is identified by its idempotence key — the frame's checksum, the
configuration hash (which covers every override, ``options.time_policy`` included), the
reference solution's hash (``MPWREFH``), the sorted drift-arc checksums, the layout preset
and the hash of its CSV file (``MPLAYOUT``, ``MPLAYH``), the hash of the masters' file
checksums (``MPMKH``) and the frame's time zone (``MPTZ``); a product on disk with the same
key is skipped unless ``options.force``. A product whose key changed (a ``reference``
product whose exposure now has a drift arc, a ``nominal`` one once a reference exists, an
edited layout, a reference re-solved under the same identifier, a rebuilt master, another
time policy) is reduced again and listed in ``reduced_again``; one written before a part
of the key existed is reduced again once.

Errors follow ``options.on_error``: with ``"continue"`` a frame's exception is logged and
recorded in ``failed`` and the run goes on, and a failure of the masters or traces aborts
its group with a recorded failure; with ``"raise"`` it propagates. A night folder that is
not a directory raises ``FileNotFoundError`` before anything is written; one without a
readable frame records a failure for the night. The ``mastpipe`` logger writes
``out/logs/quicklook_<UTC stamp>.log`` and the console for the run; both handlers are
removed when it ends.
"""

from __future__ import annotations

import hashlib
import logging
from collections import Counter, defaultdict
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from functools import partial
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.__version__ import __version__
from mastpipe.config.loader import build_effective, calib_hash, config_hash
from mastpipe.config.models import DriftParams, EffectiveConfig, InstrumentConfig, StepParams
from mastpipe.core.model import (
    PER_FIBER_COLUMNS,
    CalibrationSet,
    DriftModel,
    ExposureSpectra,
    Frame,
    FrameMeta,
    FrameRef,
    MasterCalib,
    QCReport,
    ReferenceSolution,
    Spectrum,
    TraceSet,
)
from mastpipe.core.times import to_utc, utc_iso
from mastpipe.instruments.registry import (
    layout_hash,
    load_header_map,
    load_instrument,
    load_layout,
)
from mastpipe.io import products
from mastpipe.io.ingest import (
    ASSUMED_UTC,
    apply_manifest,
    expected_shape,
    file_checksum,
    load_frame,
    read_manifest,
    scan_folder,
)
from mastpipe.pipeline.arcs import pair_arcs, status_for
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.identify import CALIB_FOLDER, _arc_lines, _provenance, _solve_and_package
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.calib import (
    build_calibration_set,
    build_master,
    calibrate_frame,
    calibration_steps,
)
from mastpipe.steps.classify import ARC_LIKE, CONTINUUM, classify_light
from mastpipe.steps.combine import combine_spectra, group_by_target
from mastpipe.steps.cosmics import flag_cosmic_rays
from mastpipe.steps.drift import apply_drift, drift_unusable_reason, fit_drift
from mastpipe.steps.extract import extract_aperture
from mastpipe.steps.fluxcal import apply_response, flux_status
from mastpipe.steps.geometry import GratingGeometry, grating_geometry
from mastpipe.steps.sky import subtract_sky_all
from mastpipe.steps.telluric import correct_telluric, telluric_status
from mastpipe.steps.trace import (
    classify_arc,
    find_traces,
    layout_indexed,
    match_layout,
    relative_throughput,
    trace_source,
)
from mastpipe.steps.wavemodel import model_from_config
from mastpipe.vocab import GRATING_NAMES

log = logging.getLogger("mastpipe")

SPEC_FOLDER = "spec"
QC_FOLDER = "qa"
LOG_FOLDER = "logs"
MASTER_KINDS = ("bias", "dark", "background")
"""The master kinds stacked from raw frames, in the order they are built."""
TIME_POLICIES = ("nearest_before", "nearest", "interpolate")
ALREADY_REDUCED = "already reduced"
"""Start of the ``skipped`` reason of a science product whose idempotence key is unchanged."""
NO_FRAMES = "no readable frame in the night folder"
"""The failure recorded for a night folder that holds nothing to reduce."""


# ------------------------------------------------------------------------------ options


@dataclass
class RunOptions:
    """What a quick-look run is told from outside.

    ``overrides`` are nested parameter mappings layered on the instrument configuration
    (see :func:`~mastpipe.config.loader.parse_override`); ``layout`` names the fiber
    preset (default: the frames' ``fiber_preset`` from a manifest, else the instrument's
    default);
    ``reference_id`` pins the reference solution instead of the one valid at each frame's
    time; ``time_policy`` overrides ``params.drift.time_policy`` (as the last override, so
    the configuration hash covers it); ``force`` rebuilds masters and science products that
    are up to date; ``library`` is the calibration library folder (default
    ``out_dir/calib``).
    """

    instrument: str
    configuration: str | None = None
    camera: str | None = None
    layout: str | None = None
    overrides: Sequence[dict[str, Any]] = ()
    on_error: Literal["continue", "raise"] = "continue"
    force: bool = False
    manifest: Path | None = None
    reference_id: str | None = None
    time_policy: str | None = None
    library: Path | None = None


@dataclass
class RunSummary:
    """What a run produced: product and QC paths, and the frames it skipped, failed on or
    reduced again (``skipped`` holds ``"<source>: <reason>"`` strings, ``failed`` pairs of
    source and error, ``reduced_again`` the sources of rewritten science products)."""

    products: list[Path] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    failed: list[tuple[str, str]] = field(default_factory=list)
    qc: list[Path] = field(default_factory=list)
    reduced_again: list[str] = field(default_factory=list)


@dataclass(eq=False)
class _Group:
    """The state of one ``(camera, configuration)`` group while it is reduced."""

    effective: EffectiveConfig
    instrument: InstrumentConfig
    preset: str
    layout: pd.DataFrame
    library: CalibLibrary
    out_dir: Path
    options: RunOptions
    summary: RunSummary
    cfg_hash: str
    binning: tuple[int, int]
    readout_key: str
    """:attr:`~mastpipe.core.model.FrameMeta.readout_key` of the group (:func:`_group_readout`):
    the readout mode and amplifier library masters must share."""
    amplifier: str | None
    """The group's output amplifier, whose gain and read noise calibrate its frames."""
    temperature: float | None
    near_time: str | None
    timezone: str | None = None
    """How the group's frame times were read (:func:`_timezone` of their metadata)."""
    layout_hash: str | None = None
    """:func:`~mastpipe.instruments.registry.layout_hash` of the group's preset."""
    calib_hash: str | None = None
    """:func:`~mastpipe.config.loader.calib_hash` of ``effective``: masters are reused under it."""
    masters: dict[str, tuple[MasterCalib, Path]] = field(default_factory=dict)
    masters_hash: str | None = None
    """:func:`_masters_hash` of ``masters``, set once they are chosen."""
    calib: CalibrationSet | None = None
    traceset: TraceSet | None = None
    traces_path: Path | None = None
    throughput: dict[int, float] = field(default_factory=dict)
    arcs: list[tuple[FrameMeta, DriftModel, str]] = field(default_factory=list)
    """The night's drift models with the identifier of the reference each was fitted against."""

    @property
    def params(self) -> StepParams:
        return self.effective.params

    @property
    def geom(self) -> GratingGeometry:
        return grating_geometry(self.effective.geometry, self.effective.detector.pixel_size_um)

    @property
    def drift_params(self) -> DriftParams:
        """``params.drift``; ``options.time_policy`` is already in it (and so in the
        configuration hash) as an override of the effective configuration."""
        return self.params.drift

    @property
    def spec_dir(self) -> Path:
        return self.out_dir / SPEC_FOLDER

    @property
    def qc_dir(self) -> Path:
        return self.out_dir / QC_FOLDER

    def traces(self) -> TraceSet:
        assert self.traceset is not None
        return self.traceset

    def calibration_names(self) -> list[str]:
        """File names of the masters and the trace product a science frame was reduced with."""
        names = [path.name for _, path in self.masters.values()]
        if self.traces_path is not None:
            names.append(self.traces_path.name)
        return names

    def injected_masters(self) -> dict[str, MasterCalib]:
        """The bias and background masters for :func:`~mastpipe.pipeline.identify._arc_lines`,
        each named after its file so the provenance lists the file."""
        return {
            kind: replace(master, provenance={**master.provenance, "name": path.name})
            for kind, (master, path) in self.masters.items()
            if kind in ("bias", "background")
        }


# ------------------------------------------------------------------------------ helpers


def _describe(error: BaseException) -> str:
    return f"{type(error).__name__}: {error}"


def _time_key(ref: FrameRef) -> tuple[bool, str, str]:
    return (ref.meta.mid_time is None, ref.meta.mid_time or "", ref.meta.source)


def _skip(group: _Group, source: str, reason: str) -> None:
    group.summary.skipped.append(f"{source}: {reason}")
    log.info("skipped %s: %s", source, reason)


def _per_frame(group: _Group, ref: FrameRef, action: Callable[[_Group, FrameRef], None]) -> None:
    """Run ``action`` on one frame under the error policy."""
    try:
        action(group, ref)
    except Exception as error:
        if group.options.on_error == "raise":
            raise
        group.summary.failed.append((ref.meta.source, _describe(error)))
        log.exception("failed %s", ref.meta.source)


def _by_exposure(refs: list[FrameRef]) -> list[list[FrameRef]]:
    """``refs`` batched by exposure time, shortest first (frames without one last)."""
    batches: dict[float | None, list[FrameRef]] = {}
    for ref in refs:
        batches.setdefault(ref.meta.exposure_time, []).append(ref)
    return [batches[t] for t in sorted(batches, key=lambda t: (t is None, t or 0.0))]


def _median_snr(spectrum: Spectrum) -> float | None:
    good = (
        ~spectrum.mask
        & np.isfinite(spectrum.flux)
        & np.isfinite(spectrum.error)
        & (spectrum.error > 0)
    )
    if not good.any():
        return None
    return float(np.median(spectrum.flux[good] / spectrum.error[good]))


def _write_qc(group: _Group, name: str, stage: str, values: dict[str, Any]) -> Path:
    path = products.write_qc(group.qc_dir / name, QCReport(stage, values))
    group.summary.qc.append(path)
    return path


def _master_meta(group: _Group, master: MasterCalib) -> FrameMeta:
    """A :class:`~mastpipe.core.model.FrameMeta` standing for a master (from the library,
    whose frames are not at hand), for naming and for the frames built from it."""
    validity = master.validity
    binning = validity.get("binning") or group.binning
    return FrameMeta(
        instrument=group.effective.instrument,
        camera=group.effective.camera,
        exposure_type=master.kind,
        configuration=group.effective.configuration,
        exposure_time=validity.get("exposure_time"),
        mid_time=validity.get("mid_time"),
        binning=(int(binning[0]), int(binning[1])),
        readout_mode=validity.get("readout_mode"),
        amplifier=validity.get("amplifier"),
        detector_temperature=validity.get("detector_temperature"),
        fiber_preset=group.preset,
        source=f"master {master.kind}",
    )


def _calibrated(group: _Group, ref: FrameRef, cosmics: bool = False) -> tuple[Frame, list[str]]:
    """``ref`` loaded and calibrated with the group's masters (background subtracted when
    ``params.calib.subtract_background`` and a background master exists), cosmic rays
    flagged into the mask when ``cosmics``; with the steps applied."""
    detector, params = group.effective.detector, group.params
    calib = group.calib
    assert calib is not None
    frame = load_frame(ref, detector)
    subtract = params.calib.subtract_background and calib.background is not None
    out = calibrate_frame(
        frame, calib, subtract_background=subtract, saturation=detector.saturation
    )
    steps = calibration_steps(calib, float(ref.meta.exposure_time or 0.0), subtract)
    if cosmics:
        assert out.variance is not None
        before = 0 if out.mask is None else int(out.mask.sum())
        mask = flag_cosmic_rays(
            out.data, out.variance, out.mask, params.calib.cosmic_sigma, params.calib.cosmic_size
        )
        out = Frame(data=out.data, meta=out.meta, variance=out.variance, mask=mask)
        steps.append(
            f"cosmic rays: {int(mask.sum()) - before} px flagged "
            f"(sigma {params.calib.cosmic_sigma:g}, size {params.calib.cosmic_size})"
        )
    return out, steps


# ------------------------------------------------------------------------------ masters


def _timezone(metas: Iterable[FrameMeta]) -> str | None:
    """The ``FrameMeta.timezone`` labels of ``metas`` for ``Provenance.timezone``: the one
    label they share, the distinct labels joined by commas when they differ, ``None`` when
    no frame has a time."""
    zones = sorted({meta.timezone for meta in metas if meta.timezone is not None})
    return ",".join(zones) if zones else None


def _master_provenance(
    group: _Group, frames: list[Frame], calibrations: list[str], steps: list[str]
) -> Provenance:
    return Provenance(
        pipeline_version=__version__,
        config_hash=group.cfg_hash,
        inputs=[frame.meta.checksum for frame in frames],
        calibrations=calibrations,
        steps=steps,
        timezone=_timezone(frame.meta for frame in frames),
    )


def _science_or_all(refs: list[FrameRef]) -> list[FrameRef]:
    """The frames of ``refs`` that may be science (``science`` and the not yet resolved
    ``light`` ones), else all of ``refs``: the frames whose readout a group follows."""
    candidates = [ref for ref in refs if ref.meta.exposure_type in ("science", "light")]
    return candidates or refs


def _group_binning(refs: list[FrameRef]) -> tuple[int, int]:
    """The binning of a group: the most common among the frames that may be science
    (``science`` and the not yet resolved ``light`` ones), else among all ``refs``."""
    return Counter(ref.meta.binning for ref in _science_or_all(refs)).most_common(1)[0][0]


def _group_readout(refs: list[FrameRef]) -> tuple[str, str | None]:
    """The readout key and amplifier of a group: the most common readout key among the
    frames that may be science, else among all ``refs`` (as :func:`_group_binning`), with
    the amplifier of the frames that carry it."""
    candidates = _science_or_all(refs)
    key = Counter(ref.meta.readout_key for ref in candidates).most_common(1)[0][0]
    amplifier = next(ref.meta.amplifier for ref in candidates if ref.meta.readout_key == key)
    return key, amplifier


def _check_binning(group: _Group, ref: FrameRef) -> None:
    """``ValueError`` when ``ref`` is not binned as the group."""
    if ref.meta.binning != group.binning:
        (x, y), (gx, gy) = ref.meta.binning, group.binning
        raise ValueError(f"binning {x}x{y} is not the group's {gx}x{gy}")


def _unlisted_amplifiers(refs: Iterable[FrameRef], instrument: InstrumentConfig) -> list[str]:
    """``"<camera>: <amplifier>"`` for every amplifier named by ``refs`` that a camera with
    a non-empty ``readouts`` table does not list (sorted): such frames take the camera's
    default gain and read noise, which for a misspelt amplifier is a silent mistake."""
    unlisted = set()
    for ref in refs:
        detector = instrument.cameras.get(ref.meta.camera)
        amplifier = ref.meta.amplifier
        if detector is None or not detector.readouts or amplifier is None:
            continue
        if amplifier not in detector.readouts:
            unlisted.add(f"{ref.meta.camera}: {amplifier}")
    return sorted(unlisted)


def _check_readout(group: _Group, ref: FrameRef) -> None:
    """``ValueError`` when ``ref`` was not read out as the group (another readout mode or
    output amplifier: :attr:`~mastpipe.core.model.FrameMeta.readout_key`)."""
    if ref.meta.readout_key != group.readout_key:
        raise ValueError(f"readout {ref.meta.readout_key} is not the group's {group.readout_key}")


def _usable_frames(group: _Group, refs: list[FrameRef], calibrate: bool) -> list[Frame]:
    """The frames of ``refs`` that can go into one master, loaded one at a time (and, when
    ``calibrate``, calibrated as a flat is: background subtracted when
    ``params.calib.flat_subtract_background`` and the set has one). A frame that fails to
    load or calibrate, whose binning is not ``group.binning``, whose readout key is not
    ``group.readout_key`` (:func:`_check_readout`) or whose loaded shape is not
    the detector's at that binning (:func:`~mastpipe.io.ingest.expected_shape`) is dropped
    under the error policy of :func:`_per_frame`: recorded in ``failed`` with the reason,
    or raised under ``on_error='raise'``."""
    detector = group.effective.detector
    shape = expected_shape(detector, group.binning)
    frames: list[Frame] = []

    def use(group: _Group, ref: FrameRef) -> None:
        _check_binning(group, ref)
        _check_readout(group, ref)
        frame = load_frame(ref, detector)
        if frame.data.shape != shape:
            raise ValueError(f"shape {frame.data.shape} is not the detector's {shape}")
        if calibrate:
            calib = group.calib
            assert calib is not None
            subtract = group.params.calib.flat_subtract_background and calib.background is not None
            frame = calibrate_frame(
                frame, calib, subtract_background=subtract, saturation=detector.saturation
            )
        frames.append(frame)

    for ref in refs:
        _per_frame(group, ref, use)
    return frames


MASTERS_HASH_LENGTH = 12
"""Hex characters kept of the SHA-256 in :func:`_masters_hash`."""


def _masters_hash(masters: dict[str, tuple[MasterCalib, Path]]) -> str | None:
    """SHA-256 over the sorted file checksums of ``masters``' files (first
    :data:`MASTERS_HASH_LENGTH` hex characters), part of a science product's idempotence
    key; ``None`` without masters."""
    if not masters:
        return None
    checksums = sorted(file_checksum(path) for _, path in masters.values())
    return hashlib.sha256("\n".join(checksums).encode()).hexdigest()[:MASTERS_HASH_LENGTH]


def _library_master(
    group: _Group, kind: str, frames: list[Frame], calibrations: list[str]
) -> tuple[MasterCalib, Path] | None:
    """The library's master of ``kind`` stacked from exactly ``frames`` (by checksum) under
    the group's ``calib_hash`` with exactly ``calibrations`` applied, to reuse instead of
    building it again; ``None`` under ``options.force`` or without one."""
    if group.options.force:
        return None
    checksums = [frame.meta.checksum for frame in frames]
    assert group.calib_hash is not None
    found = group.library.find_by_inputs(kind, checksums, group.effective.camera, group.calib_hash)
    if found is None:
        return None
    stored = Provenance.from_header(fits.getheader(found[1])).calibrations
    if stored != calibrations:
        log.info("master %s in the library was calibrated with %s; building it again", kind, stored)
        return None
    log.info("master %s reused from the library: %s", kind, found[1].name)
    return found


def _masters(group: _Group, by_type: dict[str, list[FrameRef]]) -> None:
    """Stage 1: stack the night's masters into the library and build the calibration set.
    Every frame goes through :func:`_usable_frames`; a master the library already holds,
    stacked from the same surviving frames under the same ``calib_hash`` (and, for the
    flat, calibrated with the same masters, none of them built in this run), is reused
    rather than built and overwritten unless ``options.force``; a kind none of whose
    frames survive is sought in the library.
    Sets ``group.masters`` and ``group.masters_hash``."""
    detector, params, library = group.effective.detector, group.params, group.library
    chosen: dict[str, tuple[MasterCalib, Path]] = {}
    built = False
    for kind in MASTER_KINDS:
        refs = by_type.get(kind, [])
        batches = [refs] if kind == "bias" and refs else _by_exposure(refs)
        best: tuple[tuple[int, float], MasterCalib, Path] | None = None
        for batch in batches:
            frames = _usable_frames(group, batch, calibrate=False)
            if not frames:
                continue
            found = _library_master(group, kind, frames, [])
            if found is not None:
                master, path = found
            else:
                master = build_master(frames, kind, params.calib)
                steps = [f"stack {params.calib.stack_method} of {len(frames)} {kind} frames"]
                prov = _master_provenance(group, frames, [], steps)
                path = library.add(master, prov, frames[0].meta, group.calib_hash)
                group.summary.products.append(path)
                built = True
                log.info(
                    "master %s: %d frames at %s s -> %s",
                    kind,
                    len(frames),
                    master.validity.get("exposure_time"),
                    path.name,
                )
            exposure = master.validity.get("exposure_time")
            score = (len(frames), float(exposure or 0.0))
            if best is None or score > best[0]:
                best = (score, master, path)
        if best is not None:
            chosen[kind] = (best[1], best[2])
            continue
        found = library.find(
            kind,
            group.effective.camera,
            group.binning,
            group.readout_key,
            temperature=group.temperature,
            near_time=group.near_time,
        )
        if found is not None:
            chosen[kind] = found
            log.info("master %s from the library: %s", kind, found[1].name)
    if "bias" not in chosen:
        raise RuntimeError("no bias frames in the night and no bias master in the library")
    bias = chosen["bias"][0]
    dark = chosen["dark"][0] if "dark" in chosen else None
    background = chosen["background"][0] if "background" in chosen else None
    group.calib = build_calibration_set(bias, dark, background, detector, group.amplifier)
    applied = [path.name for _, path in chosen.values()]

    flats = by_type.get("flat", [])
    frames = []
    if flats:
        batch = max(_by_exposure(flats), key=lambda b: (len(b), b[0].meta.exposure_time or 0.0))
        frames = _usable_frames(group, batch, calibrate=True)
    found = None if not frames or built else _library_master(group, "flat", frames, applied)
    if found is not None:
        chosen["flat"] = found
    elif frames:
        subtract = params.calib.flat_subtract_background and background is not None
        master = build_master(frames, "flat", params.calib)
        steps = calibration_steps(group.calib, float(frames[0].meta.exposure_time or 0.0), subtract)
        steps.append(f"stack {params.calib.stack_method} of {len(frames)} flat frames")
        prov = _master_provenance(group, frames, applied, steps)
        path = library.add(master, prov, frames[0].meta, group.calib_hash)
        group.summary.products.append(path)
        chosen["flat"] = (master, path)
        log.info("master flat: %d of %d frames -> %s", len(frames), len(flats), path.name)
    else:
        found = library.find(
            "flat",
            group.effective.camera,
            group.binning,
            group.readout_key,
            configuration=group.effective.configuration,
            temperature=group.temperature,
            near_time=group.near_time,
        )
        if found is not None:
            chosen["flat"] = found
            log.info("master flat from the library: %s", found[1].name)
    group.masters = chosen
    group.masters_hash = _masters_hash(chosen)


# ------------------------------------------------------------------------------- traces


def _resolve_light(group: _Group, by_type: dict[str, list[FrameRef]]) -> None:
    """Stage 1b: every ``light`` frame, in time order and under the error policy, calibrated
    with the group's masters (no cosmic-ray flagging) and classified from its content: moved
    into ``by_type["arc"]`` or ``by_type["science"]`` with its ``classification_rule``, or
    skipped with its reason; a frame not read out as the group fails (:func:`_check_readout`).
    ``by_type`` loses its ``light`` entry."""

    def resolve(group: _Group, ref: FrameRef) -> None:
        _check_readout(group, ref)
        frame, _ = _calibrated(group, ref)
        kind, rule = classify_light(frame, group.params.classify)
        log.info("light %s: %s (%s)", ref.meta.source, kind, rule)
        if kind not in ("arc", "science"):
            _skip(group, ref.meta.source, f"exposure type {kind} ({rule})")
            return
        arc_kind = ref.meta.arc_kind if kind == "arc" else None
        meta = replace(ref.meta, exposure_type=kind, classification_rule=rule, arc_kind=arc_kind)
        by_type.setdefault(kind, []).append(replace(ref, meta=meta))

    for ref in sorted(by_type.pop("light", []), key=_time_key):
        _per_frame(group, ref, resolve)


def _content_steps(meta: FrameMeta) -> list[str]:
    """The provenance step recording a type resolved from the frame's content, if it was."""
    if meta.classification_rule in (ARC_LIKE, CONTINUUM):
        return [
            f"exposure type {meta.exposure_type} from the frame's content "
            f"({meta.classification_rule})"
        ]
    return []


def _full_arc_candidates(group: _Group, arc_refs: list[FrameRef]) -> list[Frame]:
    """Calibrated arcs that may serve as the trace source, the one lighting the most traces
    first: the arcs the header calls full, else the arcs of unknown kind. Each candidate
    runs under the error policy of :func:`_per_frame`: one not binned or read out as the
    group or that cannot be calibrated is recorded in ``failed`` and removed from
    ``arc_refs``, so that the arc stage does not record it again; the others still
    compete."""
    known = [ref for ref in arc_refs if ref.meta.arc_kind == "full"]
    candidates = known or [ref for ref in arc_refs if ref.meta.arc_kind is None]
    trace_params = group.params.trace.model_copy(update={"collapse": "percentile"})
    scored: list[tuple[int, int, Frame]] = []

    def score(group: _Group, ref: FrameRef, n: int) -> None:
        _check_binning(group, ref)
        _check_readout(group, ref)
        frame, _ = _calibrated(group, ref)
        scored.append((len(find_traces(frame.data, trace_params)), -n, frame))

    for n, ref in enumerate(candidates):
        before = len(scored)
        _per_frame(group, ref, partial(score, n=n))
        if len(scored) == before:
            arc_refs.remove(ref)
    scored.sort(key=lambda item: (item[0], item[1]), reverse=True)
    return [frame for _, _, frame in scored]


def _traces(group: _Group, by_type: dict[str, list[FrameRef]]) -> None:
    """Stage 2: find and match the traces, write the trace product and its QC."""
    params = group.params
    flat_frame: Frame | None = None
    flat_name: list[str] = []
    inputs: list[str] = []
    if "flat" in group.masters:
        master, path = group.masters["flat"]
        flat_frame = Frame(
            data=master.data, meta=_master_meta(group, master), variance=master.variance
        )
        flat_name = [path.name]
        inputs = [str(c) for c in master.provenance.get("inputs", [])]
    arcs = [] if flat_frame is not None else _full_arc_candidates(group, by_type.get("arc", []))
    source, frame = trace_source({"flat": [flat_frame] if flat_frame else [], "arc_full": arcs})
    if frame is None:
        raise RuntimeError("no flat master and no full arc: nothing to find traces on")
    trace_params = params.trace
    if source == "arc_full":
        trace_params = trace_params.model_copy(update={"collapse": "percentile"})
        inputs = [frame.meta.checksum]
    matched = match_layout(
        find_traces(frame.data, trace_params), group.layout, params.trace.layout_tolerance
    )
    extras = [trace for trace in matched if trace.role is None]
    traceset = layout_indexed(matched)
    where = source.replace("_", " ")
    if len(traceset) == 0:
        raise RuntimeError(
            f"no trace of layout {group.preset} found on the {where} ({len(matched)} detected)"
        )
    if extras:
        log.warning(
            "%d detection(s) on the %s match no row of layout %s and are dropped: y = %s",
            len(extras),
            where,
            group.preset,
            ", ".join(f"{trace.y_center:.1f}" for trace in extras),
        )
    path = group.out_dir / CALIB_FOLDER / products.product_name(frame.meta, "traces")
    prov = Provenance(
        pipeline_version=__version__,
        config_hash=group.cfg_hash,
        inputs=inputs,
        calibrations=flat_name,
        steps=[
            f"traces from {where} ({len(matched)} detected, {len(traceset)} matched to layout "
            f"{group.preset}, {len(traceset.missing)} missing, {len(extras)} extra dropped)"
        ],
        layout=group.preset,
        timezone=frame.meta.timezone if source == "arc_full" else group.timezone,
    )
    products.write_traceset(path, traceset, prov)
    group.summary.products.append(path)
    _write_qc(
        group,
        path.with_suffix(".json").name,
        "trace",
        {
            "source": source,
            "n_detected": len(matched),
            "n_traces": len(traceset),
            "layout": group.preset,
            "layout_shift": traceset.layout_shift,
            "missing": list(traceset.missing),
            "extra": list(traceset.extra),
            "extra_y_center": [t.y_center for t in extras],
            "extra_snr": [t.snr for t in extras],
            "index": [t.index for t in traceset],
            "y_center": [t.y_center for t in traceset],
            "rms": [t.rms for t in traceset],
            "snr": [t.snr for t in traceset],
            "roles": [t.role for t in traceset],
            "telescopes": [t.telescope for t in traceset],
        },
    )
    group.traceset = traceset
    group.traces_path = path
    group.throughput = (
        relative_throughput(flat_frame, traceset, params.extract.halfwidth) if flat_frame else {}
    )
    log.info(
        "traces from the %s: %d of %d detections matched to layout %s, shift %s px, "
        "missing %s, %d extra dropped -> %s",
        where,
        len(traceset),
        len(matched),
        group.preset,
        "none" if traceset.layout_shift is None else f"{traceset.layout_shift:.2f}",
        list(traceset.missing),
        len(extras),
        path.name,
    )


# --------------------------------------------------------------------------------- arcs


def _reference_for(group: _Group, at_time: str | None) -> ReferenceSolution | None:
    """The reference a frame at ``at_time`` uses: the pinned ``options.reference_id``
    (``KeyError`` when the library lacks it), else the one valid at that time (the latest
    overall for a frame without a time)."""
    effective, library = group.effective, group.library
    keys = (effective.instrument, effective.configuration, group.preset)
    if group.options.reference_id is not None:
        ref = library.find_reference(*keys, identifier=group.options.reference_id)
        if ref is None:
            raise KeyError(
                f"reference {group.options.reference_id!r} is not in the library {library.root}"
            )
        return ref
    return library.find_reference(*keys, at_time=at_time)


def _arcs(group: _Group, arc_refs: list[FrameRef]) -> None:
    """Stage 3: every arc in time order, under the error policy."""
    for ref in sorted(arc_refs, key=_time_key):
        _per_frame(group, ref, _arc)


def _arc(group: _Group, ref: FrameRef) -> None:
    _check_readout(group, ref)
    meta = ref.meta
    kind = meta.arc_kind
    detail = "from the header"
    if kind is None:
        if group.instrument.arc_classification == "lit_traces":
            frame, _ = _calibrated(group, ref)
            kind, flags = classify_arc(frame, group.traces(), group.params.extract.halfwidth)
            lit = sum(1 for value in flags.values() if value)
            detail = f"lights {lit} of {len(flags)} traces"
        else:
            kind, detail = "unknown", "no arc kind in the header"
    log.info("arc %s: %s (%s)", meta.source, kind, detail)
    if kind == "full":
        _full_arc(group, ref)
    elif kind == "calfiber":
        _cal_arc(group, ref)
    else:
        _skip(group, meta.source, f"arc kind unknown ({detail})")


def _full_arc(group: _Group, ref: FrameRef) -> None:
    """A full arc re-solved from the latest reference of its configuration and preset."""
    meta, effective, library = ref.meta, group.effective, group.library
    keys = (effective.instrument, effective.configuration, group.preset)
    existing = library.reference_for_source(*keys, meta.checksum)
    if existing is not None and not group.options.force:
        _skip(group, meta.source, f"already solved into reference {existing['identifier']}")
        return
    previous = library.find_reference(*keys)
    if previous is None:
        _skip(group, meta.source, f"no reference for {'/'.join(keys)}; run identify")
        return
    lines = _arc_lines(
        ref,
        effective,
        group.layout,
        group.out_dir,
        masters=group.injected_masters(),
        traceset=group.traces(),
    )
    solution = _solve_and_package(lines, effective, group.preset, previous, 0)
    prov = _provenance(
        effective,
        lines,
        [*_content_steps(meta), f"wavecal seeded by reference {previous.identifier}"],
        layout=group.preset,
    )
    prov.reference_id = previous.identifier
    prov.reference_hash = previous.hash or None
    path = library.add_reference(solution, prov, effective.camera)
    group.summary.products.append(path)
    _write_qc(
        group,
        products.product_name(meta, "wref", ".json"),
        "reference",
        {
            "source": meta.source,
            "classification_rule": meta.classification_rule,
            "identifier": solution.identifier,
            "valid_from": solution.valid_from,
            "seeded_by": previous.identifier,
            "n_fibers": len(solution.solutions),
            "n_traces": len(lines.traceset),
            "fibers": sorted(solution.solutions),
            "rms_px": [solution.solutions[f].rms_px for f in sorted(solution.solutions)],
            "n_lines": [solution.solutions[f].n_lines for f in sorted(solution.solutions)],
        },
    )
    log.info(
        "full arc %s: reference %s (%d of %d fibers solved, seeded by %s) -> %s",
        meta.source,
        solution.identifier,
        len(solution.solutions),
        len(lines.traceset),
        previous.identifier,
        path.name,
    )


def _cal_arc(group: _Group, ref: FrameRef) -> None:
    """A calibration-fiber arc fitted for its drift against the reference valid at its time."""
    meta, effective = ref.meta, group.effective
    reference = _reference_for(group, meta.mid_time)
    if reference is None:
        _skip(group, meta.source, f"no reference valid at {meta.mid_time}; run identify")
        return
    arc_traces = group.traces().by_role("arc")
    if not arc_traces:
        _skip(group, meta.source, f"layout {group.preset} has no arc fiber")
        return
    lines = _arc_lines(
        ref,
        effective,
        group.layout,
        group.out_dir,
        masters=group.injected_masters(),
        traceset=TraceSet(traces=arc_traces, shape=group.traces().shape),
    )
    tables = {trace.index: table for trace, table in zip(arc_traces, lines.tables, strict=True)}
    drift = fit_drift(tables, reference, group.traces(), meta, group.drift_params)
    per_fiber = drift.per_fiber
    shifts = per_fiber["shift_px"].to_numpy(dtype=float)
    worst = float(np.nanmax(np.abs(shifts))) if len(shifts) and np.isfinite(shifts).any() else 0.0
    unusable = drift_unusable_reason(drift, group.drift_params)
    if unusable is not None:
        log.warning(
            "drift arc %s: model not usable: %s (%d of %d fibers usable, scatter %.3f px)",
            meta.source,
            unusable,
            drift.n_usable,
            len(per_fiber),
            drift.scatter_px,
        )
    elif worst > group.drift_params.warn_shift_px:
        log.warning(
            "drift arc %s: shift %.2f px exceeds warn_shift_px %.2f",
            meta.source,
            worst,
            group.drift_params.warn_shift_px,
        )
    group.arcs.append((meta, drift, reference.identifier))
    _write_qc(
        group,
        products.product_name(meta, "drift", ".json"),
        "drift",
        {
            "source": meta.source,
            "classification_rule": meta.classification_rule,
            "arc_mid_time": drift.arc_mid_time,
            "reference_id": reference.identifier,
            "usable": bool(drift.usable),
            "unusable_reason": unusable,
            "n_usable": int(drift.n_usable or 0),
            "min_fibers": int(group.drift_params.min_fibers),
            "scatter_px": float(drift.scatter_px),
            "max_abs_shift_px": worst,
            "shift_coeffs": [float(c) for c in drift.shift_coeffs],
            "stretch_coeffs": [float(c) for c in drift.stretch_coeffs],
            "per_fiber": per_fiber.to_dict(orient="records"),
        },
    )
    log.info(
        "drift arc %s: shift %.3f px at the detector centre, scatter %.3f px, usable %s",
        meta.source,
        float(drift.shift_at(group.traces().shape[0] / 2.0)),
        drift.scatter_px,
        drift.usable,
    )


# ------------------------------------------------------------------------------ science


def _arcs_for(
    arcs: list[tuple[FrameMeta, DriftModel, str]], reference: ReferenceSolution
) -> list[tuple[FrameMeta, DriftModel]]:
    """The drift models fitted against ``reference``: a shift measured against another
    reference (a full arc solved mid-night) must never be applied to this one's map."""
    return [(meta, drift) for meta, drift, identifier in arcs if identifier == reference.identifier]


def _arc_weights(pairs: list[tuple[DriftModel, float]], checksums: list[str]) -> list[float]:
    """The pairing weight of every checksum of ``checksums`` (the arcs' union), for QC."""
    weights = {c: 0.0 for c in checksums}
    for model, weight in pairs:
        for checksum in model.arc_checksums:
            weights[checksum] = weights.get(checksum, 0.0) + float(weight)
    return [weights[c] for c in checksums]


def _blend_drifts(
    pairs: list[tuple[DriftModel, float]],
) -> tuple[DriftModel, dict[int, float] | None]:
    """One :class:`~mastpipe.core.model.DriftModel` for the product: the model itself for a
    single arc (weights ``None``: every fiber 1.0); for several, the weighted blend of their
    polynomials (right-aligned when the degrees differ) and per-fiber values (``n_lines`` the
    smallest, ``usable`` when all are), the sorted union of their checksums, the weighted
    mid-time, the largest scatter, the smallest ``n_usable`` and ``usable`` when every model
    is. The returned weights feed the product's ``WEIGHT`` column: per fiber, the summed
    normalised weight of the models that hold the fiber (1.0 when all do), which says how
    much of the blend that row represents; the per-arc weights themselves are recorded in
    the science QC JSON (``drift_weights``, aligned with ``drift_arcs``)."""
    if len(pairs) == 1:
        return pairs[0][0], None
    total = sum(weight for _, weight in pairs)
    pairs = [(model, weight / total) for model, weight in pairs]
    n_shift = max(len(model.shift_coeffs) for model, _ in pairs)
    n_stretch = max(len(model.stretch_coeffs) for model, _ in pairs)
    shift = np.zeros(n_shift)
    stretch = np.zeros(n_stretch)
    for model, weight in pairs:
        shift[n_shift - len(model.shift_coeffs) :] += weight * np.asarray(model.shift_coeffs)
        stretch[n_stretch - len(model.stretch_coeffs) :] += weight * np.asarray(
            model.stretch_coeffs
        )
    tables = [(model.per_fiber.set_index("fiber"), weight) for model, weight in pairs]
    fibers = sorted({int(f) for table, _ in tables for f in table.index})
    rows: list[dict[str, Any]] = []
    fiber_weights: dict[int, float] = {}
    for fiber in fibers:
        present = [(table.loc[fiber], weight) for table, weight in tables if fiber in table.index]
        weight_sum = sum(weight for _, weight in present)
        fiber_weights[fiber] = float(weight_sum)
        rows.append(
            {
                "fiber": fiber,
                "y_center": float(present[0][0]["y_center"]),
                "shift_px": sum(w * float(r["shift_px"]) for r, w in present) / weight_sum,
                "stretch": sum(w * float(r["stretch"]) for r, w in present) / weight_sum,
                "n_lines": int(min(int(r["n_lines"]) for r, _ in present)),
                "rms_px": sum(w * float(r["rms_px"]) for r, w in present) / weight_sum,
                "usable": bool(all(bool(r["usable"]) for r, _ in present)),
            }
        )
    origin = to_utc(pairs[0][0].arc_mid_time)
    offset = sum(
        ((to_utc(model.arc_mid_time) - origin).total_seconds() * weight for model, weight in pairs),
        0.0,
    )
    mid_time = utc_iso(origin + pd.Timedelta(seconds=offset))
    blended = DriftModel(
        arc_checksums=tuple(sorted({c for model, _ in pairs for c in model.arc_checksums})),
        arc_mid_time=mid_time,
        per_fiber=pd.DataFrame(rows, columns=list(PER_FIBER_COLUMNS)),
        shift_coeffs=shift,
        stretch_coeffs=stretch,
        scatter_px=max(float(model.scatter_px) for model, _ in pairs),
        usable=all(model.usable for model, _ in pairs),
        n_pix=pairs[0][0].n_pix,
        n_usable=min(int(model.n_usable or 0) for model, _ in pairs),
    )
    return blended, fiber_weights


def _wavelengths(
    group: _Group,
    reference: ReferenceSolution | None,
    pairs: list[tuple[DriftModel, float]],
    n_pix: int,
) -> tuple[str, dict[int, np.ndarray], DriftModel | None, dict[int, float] | None, dict[int, str]]:
    """``(status, wavelength by fiber, drift model, weights, skipped)`` of a science frame.

    ``skipped`` maps the fibers that keep no wavelength to the reason (a trace centre
    beyond ``trace.layout_tolerance`` from the reference solution's, see
    :func:`~mastpipe.steps.drift.apply_drift`). ``ValueError`` when the reference was
    solved on another pixel grid than the frame's ``n_pix`` columns.
    """
    if reference is None:
        model = model_from_config(group.effective.geometry.model, group.geom, n_pix)
        lam = np.asarray(model.wavelength(np.arange(n_pix, dtype=float)), dtype=float)
        return "nominal", {trace.index: lam.copy() for trace in group.traces()}, None, None, {}
    grids = sorted({int(solution.n_pix) for solution in reference.solutions.values()})
    if grids and grids != [int(n_pix)]:
        raise ValueError(
            f"reference {reference.identifier} is solved on {grids} pixels, the frame has "
            f"{n_pix} columns: the reference belongs to another detector crop"
        )
    wavelengths, skipped = apply_drift(
        reference,
        pairs,
        group.traces(),
        group.geom,
        tolerance=group.params.trace.layout_tolerance,
        n_pix=n_pix,
    )
    if not pairs:
        return status_for(pairs), wavelengths, None, None, skipped
    drift, weights = _blend_drifts(pairs)
    return status_for(pairs), wavelengths, drift, weights, skipped


def _stored_key(path: Path) -> tuple[tuple[Any, ...], str] | None:
    """The idempotence key and wavelength status of the product at ``path``; ``None`` when
    there is none (an unreadable product counts as none and is rewritten)."""
    if not path.is_file():
        return None
    try:
        prov = Provenance.from_header(fits.getheader(path))
    except (OSError, KeyError, ValueError) as error:
        log.warning("%s cannot be read as a product (%s); it will be rewritten", path.name, error)
        return None
    key = (
        prov.inputs[0] if prov.inputs else None,
        prov.config_hash,
        prov.reference_hash,
        tuple(sorted(prov.drift_arcs)),
        prov.layout,
        prov.layout_hash,
        prov.masters_hash,
        prov.timezone,
    )
    return key, prov.wavelength_status


def _science_frame(group: _Group, ref: FrameRef) -> None:
    """Stage 4 for one science frame: wavelengths by the time policy, sky, products, QC.
    A frame not read out as the group (:func:`_check_readout`) fails rather than being
    calibrated with masters and noise of another amplifier."""
    _check_readout(group, ref)
    meta, params = ref.meta, group.params
    reference = _reference_for(group, meta.mid_time)
    arcs = [] if reference is None else _arcs_for(group.arcs, reference)
    pairs = [] if reference is None else pair_arcs(meta, arcs, reference, group.drift_params)
    drift_arcs = sorted({c for model, _ in pairs for c in model.arc_checksums})
    reference_hash = None if reference is None else (reference.hash or None)
    key = (
        meta.checksum,
        group.cfg_hash,
        reference_hash,
        tuple(drift_arcs),
        group.preset,
        group.layout_hash,
        group.masters_hash,
        meta.timezone,
    )
    path = group.spec_dir / products.product_name(meta, "spec")
    stored = _stored_key(path)
    if stored is not None:
        if stored[0] == key and not group.options.force:
            _skip(group, meta.source, f"{ALREADY_REDUCED} as {path.name}")
            return
        why = "forced" if stored[0] == key else f"its key changed (was {stored[1]})"
        log.info("reducing %s again: %s", meta.source, why)

    frame, steps = _calibrated(group, ref, cosmics=True)
    steps = [*_content_steps(meta), *steps]
    traceset = group.traces()
    spectra = extract_aperture(frame, traceset, params.extract)
    steps.append(
        f"extract {params.extract.aperture} aperture of halfwidth {params.extract.halfwidth:g} "
        f"({params.extract.method})"
    )
    status, wavelengths, drift, weights, skipped = _wavelengths(
        group, reference, pairs, frame.shape[1]
    )
    for spectrum in spectra:
        if spectrum.fiber in wavelengths:
            spectrum.wavelength = wavelengths[spectrum.fiber]
    withheld = [
        {"fiber": int(fiber), "reason": reason}
        for fiber, reason in sorted(skipped.items())
        if any(spectrum.fiber == fiber for spectrum in spectra)
    ]
    for item in withheld:
        log.warning(
            "%s: fiber %d is given no wavelength: %s", meta.source, item["fiber"], item["reason"]
        )
    if reference is None:
        steps.append("wavelength nominal from the configuration geometry")
    elif drift is None:
        steps.append(f"wavelength from reference {reference.identifier} (no drift arc paired)")
    else:
        steps.append(
            f"wavelength transferred from reference {reference.identifier} through "
            f"{len(pairs)} drift arc(s) ({group.drift_params.time_policy})"
        )
    if withheld:
        steps.append(
            f"wavelength withheld from {len(withheld)} fiber(s) whose trace is off the "
            "reference's centre"
        )
    spectra, sky_reasons = subtract_sky_all(spectra, traceset, group.throughput, params.sky)
    n_sky = sum(1 for spectrum in spectra if "skysub" in spectrum.versions)
    steps.append(f"skysub: {n_sky} science fibers subtracted, {len(sky_reasons)} left raw")
    spectra = [correct_telluric(apply_response(spectrum, None), None) for spectrum in spectra]
    es = ExposureSpectra(
        meta=meta,
        spectra=spectra,
        wavelength_status=status,
        flux_status=flux_status(spectra),
        telluric_status=telluric_status(spectra),
        reference_id=None if reference is None else reference.identifier,
        reference_hash=reference_hash,
        drift_arcs=tuple(drift_arcs),
    )
    prov = Provenance(
        pipeline_version=__version__,
        config_hash=group.cfg_hash,
        inputs=[meta.checksum],
        calibrations=group.calibration_names(),
        wavelength_status=es.wavelength_status,
        flux_status=es.flux_status,
        telluric_status=es.telluric_status,
        reference_id=es.reference_id,
        reference_hash=es.reference_hash,
        drift_arcs=list(drift_arcs),
        steps=steps,
        layout=group.preset,
        timezone=meta.timezone,
        layout_hash=group.layout_hash,
        masters_hash=group.masters_hash,
    )
    products.write_spectra(path, es, prov, drift=drift, weights=weights)
    group.summary.products.append(path)
    if stored is not None:
        group.summary.reduced_again.append(meta.source)
    traces = {trace.index: trace for trace in traceset}
    _write_qc(
        group,
        path.with_suffix(".json").name,
        "science",
        {
            "source": meta.source,
            "classification_rule": meta.classification_rule,
            "mid_time": meta.mid_time,
            "targets": dict(meta.targets),
            "wavelength_status": status,
            "reference_id": es.reference_id,
            "drift_arcs": list(drift_arcs),
            "drift_weights": _arc_weights(pairs, drift_arcs),
            "drift_scatter_px": None if drift is None else float(drift.scatter_px),
            "wavelength_skipped": withheld,
            "fibers": [s.fiber for s in spectra],
            "roles": [s.role for s in spectra],
            "telescopes": [s.telescope for s in spectra],
            "median_snr": [_median_snr(s) for s in spectra],
            "masked_fraction": [float(np.mean(s.mask)) for s in spectra],
            "drift_shift_px": [
                None if drift is None else float(drift.shift_at(traces[s.fiber].y_center))
                for s in spectra
            ],
            "drift_stretch": [
                None if drift is None else float(drift.stretch_at(traces[s.fiber].y_center))
                for s in spectra
            ],
            "sky": sky_reasons,
        },
    )
    log.info(
        "science %s: %d spectra, wavelength %s, %d sky-subtracted -> %s",
        meta.source,
        len(spectra),
        status,
        n_sky,
        path.name,
    )
    _combine_targets(group, es, prov)


def _combine_targets(group: _Group, es: ExposureSpectra, prov: Provenance) -> None:
    """Stage 5: one combined product per target of the exposure (fibers without a target
    are not combined)."""
    meta = es.meta
    for target, members in group_by_target(es.spectra).items():
        if any(spectrum.target is None for spectrum in members):
            log.info(
                "%s: %d science fibers without a target are not combined", meta.source, len(members)
            )
            continue
        try:
            cs = combine_spectra(es.spectra, target, group.params.combine, version=None)
            path = group.spec_dir / products.combined_name(meta, target)
            steps = [*prov.steps, f"combine {len(cs.contributors)} science fibers of {target}"]
            products.write_combined(path, cs, meta, replace(prov, steps=steps))
            group.summary.products.append(path)
            log.info(
                "combined %s: %d fibers, %d excluded -> %s",
                target,
                len(cs.contributors),
                len(cs.excluded),
                path.name,
            )
        except Exception as error:
            if group.options.on_error == "raise":
                raise
            group.summary.failed.append((meta.source, f"combine {target}: {_describe(error)}"))
            log.exception("failed combining %s of %s", target, meta.source)


# ---------------------------------------------------------------------------- the run


def _start_logging(out_dir: Path) -> tuple[list[logging.Handler], int]:
    folder = out_dir / LOG_FOLDER
    folder.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S")
    file_handler = logging.FileHandler(folder / f"quicklook_{stamp}.log", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    console.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    previous = log.level
    log.setLevel(logging.INFO)
    log.addHandler(file_handler)
    log.addHandler(console)
    return [file_handler, console], previous


def _run_group(
    refs: list[FrameRef],
    instrument: InstrumentConfig,
    camera: str,
    configuration: str,
    library: CalibLibrary,
    out_dir: Path,
    options: RunOptions,
    summary: RunSummary,
) -> None:
    presets = {ref.meta.fiber_preset for ref in refs if ref.meta.fiber_preset}
    if options.layout is not None:
        preset = options.layout
    elif len(presets) == 1:
        preset = presets.pop()
    else:
        preset = instrument.default_layout
        if presets:
            log.warning("frames carry presets %s; using the default %s", sorted(presets), preset)
    layout = load_layout(instrument.name, preset)
    overrides = list(options.overrides)
    if options.time_policy is not None:
        # an override like any --set, so that the configuration hash covers it
        overrides.append({"drift": {"time_policy": options.time_policy}})
    effective = build_effective(instrument, configuration, camera, overrides, layout=preset)
    by_type: dict[str, list[FrameRef]] = defaultdict(list)
    for ref in refs:
        by_type[ref.meta.exposure_type].append(ref)
    # binning and the time near which library masters are sought: from the frames that
    # may be science (light ones included, before their content is known)
    candidates = sorted(by_type.get("science", []) + by_type.get("light", []), key=_time_key)
    first = candidates[0].meta if candidates else refs[0].meta
    temperatures = [
        m.detector_temperature for m in (r.meta for r in refs) if m.detector_temperature is not None
    ]
    readout_key, amplifier = _group_readout(refs)
    group = _Group(
        effective=effective,
        instrument=instrument,
        preset=preset,
        layout=layout,
        library=library,
        out_dir=out_dir,
        options=options,
        summary=summary,
        cfg_hash=config_hash(effective),
        binning=_group_binning(refs),
        readout_key=readout_key,
        amplifier=amplifier,
        temperature=float(np.mean(temperatures)) if temperatures else None,
        near_time=first.mid_time,
        timezone=_timezone(ref.meta for ref in refs),
        layout_hash=layout_hash(instrument.name, preset),
        calib_hash=calib_hash(effective),
    )
    log.info(
        "group %s/%s: %s, layout %s, config %s",
        camera,
        configuration,
        ", ".join(f"{len(v)} {k}" for k, v in sorted(by_type.items())),
        preset,
        group.cfg_hash,
    )
    _masters(group, by_type)
    _resolve_light(group, by_type)
    science = sorted(by_type.get("science", []), key=_time_key)
    _traces(group, by_type)
    _arcs(group, by_type.get("arc", []))
    for ref in science:
        _per_frame(group, ref, _science_frame)


def run_quicklook(night_dir: Path | str, out_dir: Path | str, options: RunOptions) -> RunSummary:
    """Reduce the night in ``night_dir`` into ``out_dir`` (see the module docstring)."""
    night_dir, out_dir = Path(night_dir), Path(out_dir)
    if options.time_policy is not None and options.time_policy not in TIME_POLICIES:
        raise ValueError(f"time_policy must be one of {TIME_POLICIES}, got {options.time_policy!r}")
    if options.on_error not in ("continue", "raise"):
        raise ValueError(f"on_error must be 'continue' or 'raise', got {options.on_error!r}")
    if not night_dir.is_dir():
        raise FileNotFoundError(f"night folder {night_dir} is not a directory")
    out_dir.mkdir(parents=True, exist_ok=True)
    handlers, level = _start_logging(out_dir)
    summary = RunSummary()
    try:
        log.info("quicklook of %s into %s (mastpipe %s)", night_dir, out_dir, __version__)
        instrument = load_instrument(options.instrument)
        manifest = read_manifest(options.manifest) if options.manifest is not None else None
        unreadable: list[Path] = []
        other: list[tuple[str, str]] = []
        refs = scan_folder(
            night_dir,
            options.instrument,
            load_header_map(options.instrument),
            unreadable=unreadable,
            skipped=other,
        )
        for path in unreadable:
            source = path.relative_to(night_dir).as_posix()
            summary.failed.append((source, "unreadable FITS file"))
            log.error("unreadable FITS file skipped: %s", source)
        for source, reason in other:
            summary.skipped.append(f"{source}: {reason}")
            log.warning("skipped %s: %s", source, reason)
        if manifest is not None:
            unmatched: list[str] = []
            refs = apply_manifest(refs, manifest, unmatched)
            for row in unmatched:
                summary.skipped.append(f"{row}: matches no frame")
                log.warning("%s: matches no frame", row)
        log.info("scanned %d frames in %s", len(refs), night_dir)
        if any(ref.meta.timezone == ASSUMED_UTC for ref in refs):
            log.warning("frame times carry no time zone; taken as UTC (header map 'assume-utc')")
        unlisted = _unlisted_amplifiers(refs, instrument)
        if unlisted:
            log.warning(
                "no gain and read noise listed for amplifier %s; the camera defaults apply",
                ", ".join(unlisted),
            )
        if not refs:
            summary.failed.append((f"night {night_dir.name}", NO_FRAMES))
            log.error("%s: %s", night_dir, NO_FRAMES)
        library = CalibLibrary(options.library or out_dir / CALIB_FOLDER)
        default_conf = options.configuration
        groups: dict[tuple[str, str | None], list[FrameRef]] = defaultdict(list)
        for ref in refs:
            meta = ref.meta
            if meta.exposure_type == "unknown":
                reason = f"exposure type unknown ({meta.classification_rule})"
                summary.skipped.append(f"{meta.source}: {reason}")
                log.warning("skipped %s: %s", meta.source, reason)
                continue
            groups[(meta.camera, meta.configuration or default_conf)].append(ref)
        for (camera, configuration), group_refs in sorted(groups.items(), key=str):
            label = f"group {camera}/{configuration}"
            if options.camera is not None and camera != options.camera:
                summary.skipped.append(f"{label}: not camera {options.camera}")
                log.info("skipped %s: not camera %s", label, options.camera)
                continue
            if configuration is None:
                summary.failed.append(
                    (label, f"configuration unknown; pass --configuration one of {GRATING_NAMES}")
                )
                log.error("%s: configuration unknown, group skipped", label)
                continue
            try:
                _run_group(
                    group_refs,
                    instrument,
                    camera,
                    configuration,
                    library,
                    out_dir,
                    options,
                    summary,
                )
            except Exception as error:
                if options.on_error == "raise":
                    raise
                summary.failed.append((label, _describe(error)))
                log.exception("%s failed", label)
        log.info(
            "done: %d products, %d QC files, %d skipped, %d reduced again, %d failed",
            len(summary.products),
            len(summary.qc),
            len(summary.skipped),
            len(summary.reduced_again),
            len(summary.failed),
        )
    finally:
        for handler in handlers:
            log.removeHandler(handler)
            handler.close()
        log.setLevel(level)
    return summary
