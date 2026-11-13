"""The ``quicklook``, ``solve`` and ``show`` sub-commands.

Each registrar adds one sub-parser to :data:`mastpipe.cli.COMMANDS`; the command functions
turn the parsed arguments into calls of :func:`~mastpipe.pipeline.runner.run_quicklook`,
:func:`~mastpipe.pipeline.identify.run_solve` (seeded by the calibration library) and a
reader of any product file, and print what happened.

``quicklook`` exits 0 when the run is ``complete`` or ``empty``, 1 on an unexpected error,
2 on a usage error (the argument parser's), 3 when it is ``partial`` (some frame or group
failed) and 4 when it ``failed`` (a :class:`~mastpipe.core.errors.RunError` stopped it); see
:mod:`mastpipe.pipeline.summary`.
"""

from __future__ import annotations

import argparse
import os
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits

from mastpipe.config.loader import build_effective, override_argument, parse_override
from mastpipe.config.models import InstrumentConfig
from mastpipe.instruments.registry import load_instrument
from mastpipe.io import products
from mastpipe.pipeline import summary as result
from mastpipe.pipeline.callib import CalibLibrary, registry_geometry
from mastpipe.pipeline.group import stopped_by_run_error
from mastpipe.pipeline.identify import CALIB_FOLDER, run_solve
from mastpipe.pipeline.provenance import Provenance
from mastpipe.pipeline.runner import TIME_POLICIES, RunOptions, RunSummary, run_quicklook

RUN_IDENTIFY_FIRST = "run identify first"


def _resolve_camera(instrument: InstrumentConfig, camera: str | None) -> str:
    if camera is not None:
        return camera
    if len(instrument.cameras) != 1:
        raise SystemExit(f"{instrument.name} has cameras {sorted(instrument.cameras)}; pass --camera")
    return next(iter(instrument.cameras))


def _add_set(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--set",
        action="append",
        default=[],
        type=override_argument,
        metavar="KEY.PATH=VALUE",
        help="override a step parameter, e.g. extract.halfwidth=4 (repeatable)",
    )


# ------------------------------------------------------------------------------ quicklook


def register_quicklook(subparsers: argparse._SubParsersAction) -> None:
    """Add the ``quicklook`` sub-command."""
    parser = subparsers.add_parser(
        "quicklook",
        help="reduce a night folder non-interactively",
        description=(
            "Reduce every frame of a night folder: masters into the calibration library, "
            "traces, full arcs re-solved from the library's reference, drift from the "
            "calibration-fiber arcs, and wavelength-calibrated, sky-subtracted science spectra "
            "combined per target, with provenance and QC files."
        ),
    )
    parser.add_argument("night_dir", type=Path, help="folder of raw FITS frames")
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    parser.add_argument("--instrument", required=True, help="instrument name, e.g. highspec")
    parser.add_argument("--configuration", help="configuration of frames whose header has none, e.g. Mg")
    parser.add_argument("--camera", help="reduce this camera only")
    parser.add_argument(
        "--layout",
        help="fiber layout preset (default: a manifest's fiber_preset, else the instrument's)",
    )
    parser.add_argument("--manifest", type=Path, help="CSV of metadata overrides by path, pattern or file name")
    parser.add_argument(
        "--on-error",
        choices=["continue", "raise"],
        default="continue",
        help="record a frame's failure and go on (default), or stop at the first",
    )
    parser.add_argument("--force", action="store_true", help="reduce products that exist and rebuild masters")
    parser.add_argument("--reference", metavar="ID", help="use this reference solution instead of the one valid")
    parser.add_argument("--time-policy", choices=list(TIME_POLICIES), help="how drift arcs are paired in time")
    parser.add_argument("--library", type=Path, help="calibration library (default: OUT/calib)")
    parser.add_argument(
        "--summary-json",
        type=Path,
        metavar="PATH",
        help="write the run's result as JSON to PATH (any outcome but a usage error)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="keep progress lines off the console; warnings and errors still show",
    )
    _add_set(parser)
    parser.set_defaults(func=quicklook_command)


def quicklook_command(args: argparse.Namespace) -> int:
    """Run ``quicklook`` from parsed arguments and return its exit code (see
    :func:`mastpipe.pipeline.summary.outcome`): 0 ``complete`` or ``empty``, 1 ``error``,
    3 ``partial``, 4 ``failed``. A :class:`~mastpipe.core.errors.RunError` raised before the
    calibration library was open is one ``FAILED:`` line on stderr, any other exception
    (a run error re-raised from a group included) its traceback; the printed summary
    follows only a run that ended normally. It is left out silently when stdout is closed;
    a stdout that cannot take it (an encoding that cannot spell a file name, a full device)
    is one line on stderr, its descriptor is then pointed at the null device
    (:func:`_discard_stdout`), and the exit code and the result file are left alone. With
    ``--summary-json`` the result document is built and written last; a failure to do so
    is one line on stderr and changes nothing else."""
    options = RunOptions(
        instrument=args.instrument,
        configuration=args.configuration,
        camera=args.camera,
        layout=args.layout,
        overrides=[parse_override(text) for text in args.set],
        on_error=args.on_error,
        force=args.force,
        manifest=args.manifest,
        reference_id=args.reference,
        time_policy=args.time_policy,
        library=args.library,
        quiet=args.quiet,
    )
    summary = RunSummary()
    error: Exception | None = None
    try:
        run_quicklook(args.night_dir, args.out, options, summary)
    except Exception as caught:
        error = caught
        if stopped_by_run_error(summary, caught):
            print(f"FAILED: {caught}", file=sys.stderr)
        else:
            traceback.print_exception(caught, file=sys.stderr)
    if error is None and sys.stdout is not None:  # None: file descriptor 1 is closed
        try:
            _print_summary(summary, options, args.night_dir)
            sys.stdout.flush()
        except (OSError, UnicodeError) as failure:
            print(f"could not print the summary: {failure}", file=sys.stderr)
            _discard_stdout()
    if args.summary_json is not None:
        try:
            doc = result.document(summary, options, args.night_dir, args.out, args.set, error)
            result.write(args.summary_json, doc)
        except Exception as failure:
            print(f"could not write the result file {args.summary_json}: {failure}", file=sys.stderr)
    return result.outcome(summary, error)[1]


def _discard_stdout() -> None:
    """Point file descriptor 1 at the null device, so that the text a failed print left
    in stdout's buffer is flushed there at shutdown instead of failing once more (which
    would end the process with the interpreter's status 120, not the run's). A stdout
    without a real descriptor (one replaced in-process) is left alone."""
    try:
        target = sys.stdout.fileno()
    except (OSError, ValueError, AttributeError):  # io.UnsupportedOperation is both of the first two
        return
    null = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(null, target)
    finally:
        os.close(null)


def _print_summary(summary: RunSummary, options: RunOptions, night_dir: Path) -> None:
    print(
        f"{len(summary.products)} products, {len(summary.qc)} QC files, "
        f"{len(summary.skipped)} skipped, {len(summary.reduced_again)} reduced again, "
        f"{len(summary.failed)} failed"
    )
    if summary.frames_scanned == 0:
        print(f"no frames of {options.instrument} in {night_dir}")
    for item in summary.skipped:
        print(f"  skipped {item.label}: {item.reason}")
    for source in summary.reduced_again:
        print(f"  reduced again {source}")
    for failure in summary.failed:
        cause = failure.message if failure.error_type is None else f"{failure.error_type}: {failure.message}"
        print(f"  FAILED {failure.label}: {cause}")


# ---------------------------------------------------------------------------------- solve


def register_solve(subparsers: argparse._SubParsersAction) -> None:
    """Add the ``solve`` sub-command."""
    parser = subparsers.add_parser(
        "solve",
        help="re-solve a full arc seeded by the library's reference solution",
        description=(
            "Solve every fiber of a full arc starting from the latest reference solution of "
            "the configuration and layout in the calibration library, and add the result."
        ),
    )
    parser.add_argument("full_arc", type=Path, help="the full-arc FITS file")
    parser.add_argument("--out", type=Path, required=True, help="output folder")
    parser.add_argument("--instrument", required=True, help="instrument name, e.g. highspec")
    parser.add_argument("--configuration", required=True, help="configuration, e.g. Mg")
    parser.add_argument("--layout", help="fiber layout preset (default: the instrument's)")
    parser.add_argument("--camera", help="camera name (default: the instrument's only camera)")
    parser.add_argument("--library", type=Path, help="calibration library (default: OUT/calib)")
    _add_set(parser)
    parser.set_defaults(func=solve_command)


def solve_command(args: argparse.Namespace) -> int:
    """Run ``solve`` from parsed arguments; exit code 2 without a reference to seed from."""
    instrument = load_instrument(args.instrument)
    camera = _resolve_camera(instrument, args.camera)
    layout = args.layout or instrument.default_layout
    overrides = [parse_override(text) for text in args.set]
    effective = build_effective(instrument, args.configuration, camera, overrides, layout=layout)
    root = Path(args.library) if args.library is not None else Path(args.out) / CALIB_FOLDER
    library = CalibLibrary(root)
    keys = (effective.instrument, effective.configuration, layout)
    previous = library.find_reference(*keys)
    if previous is None:
        print(f"no reference for {'/'.join(keys)} in {root}; {RUN_IDENTIFY_FIRST}")
        return 2
    ref = run_solve(args.full_arc, effective, layout, args.out, previous, library=root)
    path = Path(args.out) / CALIB_FOLDER / products.reference_name(ref)
    if root.resolve() != path.parent.resolve():
        _, prov = products.read_reference(path, registry_geometry(camera))
        library.add_reference(ref, prov, camera)
    print(f"reference {ref.identifier}: {len(ref.solutions)} fibers solved, seeded by {previous.identifier}")
    print(f"written: {path}")
    return 0


# ----------------------------------------------------------------------------------- show


def register_show(subparsers: argparse._SubParsersAction) -> None:
    """Add the ``show`` sub-command."""
    parser = subparsers.add_parser(
        "show",
        help="print what a product file holds",
        description="Print the kind, statuses, provenance and per-fiber summary of a product.",
    )
    parser.add_argument("product", type=Path, help="a mastpipe FITS product")
    parser.set_defaults(func=show_command)


def _product_kind(hdul: fits.HDUList) -> str:
    names = {hdu.name for hdu in hdul}
    if "MPKIND" in hdul[0].header:
        return f"master {hdul[0].header['MPKIND']}"
    for table, kind in (
        ("SPECTRA", "spectra"),
        ("SOLUTIONS", "reference"),
        ("COMBINED", "combined"),
        ("TRACES", "traces"),
    ):
        if table in names:
            return kind
    return "unknown product"


def _fmt(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def _print_provenance(prov: Provenance) -> None:
    print(f"mastpipe {prov.pipeline_version}, configuration {prov.config_hash}")
    print(f"wavelength {prov.wavelength_status}, flux {prov.flux_status}, telluric {prov.telluric_status}")
    if prov.reference_id is not None:
        print(f"reference {prov.reference_id} ({prov.reference_hash})")
    print(f"inputs: {', '.join(prov.inputs) or '-'}")
    print(f"calibrations: {', '.join(prov.calibrations) or '-'}")
    if prov.drift_arcs:
        print(f"drift arcs: {', '.join(prov.drift_arcs)}")
    for step in prov.steps:
        print(f"  step: {step}")


def _targets_text(targets: tuple[tuple[str, str], ...]) -> str:
    """``mast01=HD 1; mast02=HD 2`` (``*`` for a target of every telescope), ``-`` if none."""
    return "; ".join(f"{telescope}={target}" for telescope, target in targets) or "-"


def _show_spectra(path: Path) -> None:
    es, _, drift = products.read_spectra(path)
    meta = es.meta
    print(f"frame {meta.source} ({meta.exposure_type}), mid-time {meta.mid_time}, targets {_targets_text(meta.targets)}")
    if drift is not None:
        print(f"drift arc {drift.arc_mid_time}: scatter {drift.scatter_px:.3f} px, usable {drift.usable}")
    print(f"{'fiber':>5} {'role':<8} {'tel':<5} {'target':<16} {'S/N':>8} {'masked':>7} versions")
    traces = {} if drift is None else {int(f): True for f in drift.per_fiber["fiber"]}
    for s in es.spectra:
        good = ~s.mask & np.isfinite(s.flux) & np.isfinite(s.error) & (s.error > 0)
        snr = float(np.median(s.flux[good] / s.error[good])) if good.any() else None
        tag = " drift" if s.fiber in traces else ""
        print(
            f"{s.fiber:>5} {s.role or '-':<8} {s.telescope or '-':<5} {s.target or '-':<16} "
            f"{_fmt(snr):>8} {float(np.mean(s.mask)):>7.3f} "
            f"{','.join(s.versions) or '-'}{tag}"
        )


def _show_reference(hdul: fits.HDUList) -> None:
    header = hdul[0].header
    print(
        f"reference {header['MPREFID']} of {header['MPINST']}/{header['MPCONF']}/"
        f"{header['MPPRESET']}, valid from {header['MPVALID']}, hash {header['MPREFH']}"
    )
    rows = hdul["SOLUTIONS"].data
    print(f"{'fiber':>5} {'lines':>5} {'rms px':>8} {'from nm':>9} {'to nm':>9} {'y':>8}")
    for i in range(len(rows)):
        print(
            f"{int(rows['FIBER'][i]):>5} {int(rows['NLINES'][i]):>5} "
            f"{float(rows['RMSPX'][i]):>8.4f} {float(rows['ANCHNM0'][i]):>9.3f} "
            f"{float(rows['ANCHNM1'][i]):>9.3f} {float(rows['YCENTER'][i]):>8.2f}"
        )


def _show_traces(hdul: fits.HDUList) -> None:
    header = hdul[0].header
    rows = hdul["TRACES"].data
    print(
        f"{len(rows)} traces on {header['MPNY']} x {header['MPNX']}, layout shift "
        f"{header['MPSHIFT']}, missing [{header['MPMISS']}], extra [{header['MPEXTRA']}]"
    )
    print(f"{'index':>5} {'y':>8} {'role':<8} {'tel':<5} {'snr':>8} {'rms':>6}")
    for i in range(len(rows)):
        print(
            f"{int(rows['INDEX'][i]):>5} {float(rows['YCENTER'][i]):>8.2f} "
            f"{str(rows['ROLE'][i]) or '-':<8} {str(rows['TELESCOPE'][i]) or '-':<5} "
            f"{float(rows['SNR'][i]):>8.1f} {float(rows['RMS'][i]):>6.3f}"
        )


def _show_combined(hdul: fits.HDUList) -> None:
    header = hdul[0].header
    columns = hdul["COMBINED"].columns.names
    versions = [name[5:].lower() for name in columns if name.startswith("FLUX_")]
    print(
        f"target {header['MPTARGET']} ({header['MPGROUP']}): {len(hdul['COMBINED'].data)} pixels, "
        f"versions {', '.join(versions)}"
    )
    contributors = hdul["CONTRIB"].data
    print(
        f"{len(contributors)} contributors: "
        + ", ".join(f"fiber {int(row['FIBER'])} ({row['TELESCOPE']}, w {float(row['WEIGHT']):.2f})" for row in contributors)
    )
    for row in hdul["EXCLUDED"].data:
        print(f"excluded fiber {int(row['FIBER'])} ({row['TELESCOPE']}): {row['REASON']}")


def _show_master(hdul: fits.HDUList) -> None:
    header = hdul[0].header
    validity = {str(header[f"MPVALK{n:02d}"]): header[f"MPVALV{n:02d}"] for n in range(1, int(header["MPNVAL"]) + 1)}
    print(f"{header['MPNFRM']} frames stacked, shape {hdul[0].data.shape}")
    print("validity: " + ", ".join(f"{k} {v}" for k, v in validity.items()))


def show_command(args: argparse.Namespace) -> int:
    """Print the kind, statuses, provenance and per-fiber summary of a product."""
    path = Path(args.product)
    with fits.open(path) as hdul:
        kind = _product_kind(hdul)
        print(f"{path.name}: {kind}")
        prov = Provenance.from_header(hdul[0].header)
        if kind == "spectra":
            _show_spectra(path)
        elif kind == "reference":
            _show_reference(hdul)
        elif kind == "traces":
            _show_traces(hdul)
        elif kind == "combined":
            _show_combined(hdul)
        elif kind.startswith("master"):
            _show_master(hdul)
        _print_provenance(prov)
    return 0
