"""The ``quicklook``, ``solve`` and ``show`` sub-commands.

Each registrar adds one sub-parser to :data:`mastpipe.cli.COMMANDS`; the command functions
turn the parsed arguments into calls of :func:`~mastpipe.pipeline.runner.run_quicklook`,
:func:`~mastpipe.pipeline.identify.run_solve` (seeded by the calibration library) and a
reader of any product file, and print what happened.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits

from mastpipe.config.loader import build_effective, parse_override
from mastpipe.config.models import InstrumentConfig
from mastpipe.instruments.registry import load_instrument
from mastpipe.io import products
from mastpipe.pipeline.callib import CalibLibrary, registry_geometry
from mastpipe.pipeline.identify import CALIB_FOLDER, run_solve
from mastpipe.pipeline.provenance import Provenance
from mastpipe.pipeline.runner import TIME_POLICIES, RunOptions, run_quicklook

RUN_IDENTIFY_FIRST = "run identify first"


def _resolve_camera(instrument: InstrumentConfig, camera: str | None) -> str:
    if camera is not None:
        return camera
    if len(instrument.cameras) != 1:
        raise SystemExit(
            f"{instrument.name} has cameras {sorted(instrument.cameras)}; pass --camera"
        )
    return next(iter(instrument.cameras))


def _add_set(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--set",
        action="append",
        default=[],
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
    parser.add_argument(
        "--configuration", help="configuration of frames whose header has none, e.g. Mg"
    )
    parser.add_argument("--camera", help="reduce this camera only")
    parser.add_argument(
        "--layout",
        help="fiber layout preset (default: a manifest's fiber_preset, else the instrument's)",
    )
    parser.add_argument(
        "--manifest", type=Path, help="CSV of metadata overrides by path, pattern or file name"
    )
    parser.add_argument(
        "--on-error",
        choices=["continue", "raise"],
        default="continue",
        help="record a frame's failure and go on (default), or stop at the first",
    )
    parser.add_argument(
        "--force", action="store_true", help="reduce products that exist and rebuild masters"
    )
    parser.add_argument(
        "--reference", metavar="ID", help="use this reference solution instead of the one valid"
    )
    parser.add_argument(
        "--time-policy", choices=list(TIME_POLICIES), help="how drift arcs are paired in time"
    )
    parser.add_argument("--library", type=Path, help="calibration library (default: OUT/calib)")
    _add_set(parser)
    parser.set_defaults(func=quicklook_command)


def quicklook_command(args: argparse.Namespace) -> int:
    """Run ``quicklook`` from parsed arguments; exit code 2 when anything failed or the
    night folder does not exist."""
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
    )
    try:
        summary = run_quicklook(args.night_dir, args.out, options)
    except FileNotFoundError as error:
        print(f"FAILED: {error}")
        return 2
    print(
        f"{len(summary.products)} products, {len(summary.qc)} QC files, "
        f"{len(summary.skipped)} skipped, {len(summary.reduced_again)} reduced again, "
        f"{len(summary.failed)} failed"
    )
    for item in summary.skipped:
        print(f"  skipped {item}")
    for source in summary.reduced_again:
        print(f"  reduced again {source}")
    for source, why in summary.failed:
        print(f"  FAILED {source}: {why}")
    return 0 if not summary.failed else 2


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
    ref = run_solve(args.full_arc, effective, layout, args.out, previous)
    path = Path(args.out) / CALIB_FOLDER / products.reference_name(ref)
    if root.resolve() != path.parent.resolve():
        _, prov = products.read_reference(path, registry_geometry(camera))
        library.add_reference(ref, prov, camera)
    print(
        f"reference {ref.identifier}: {len(ref.solutions)} fibers solved, "
        f"seeded by {previous.identifier}"
    )
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
    print(
        f"wavelength {prov.wavelength_status}, flux {prov.flux_status}, "
        f"telluric {prov.telluric_status}"
    )
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
    print(
        f"frame {meta.source} ({meta.exposure_type}), mid-time {meta.mid_time}, "
        f"targets {_targets_text(meta.targets)}"
    )
    if drift is not None:
        print(
            f"drift arc {drift.arc_mid_time}: scatter {drift.scatter_px:.3f} px, "
            f"usable {drift.usable}"
        )
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
        + ", ".join(
            f"fiber {int(row['FIBER'])} ({row['TELESCOPE']}, w {float(row['WEIGHT']):.2f})"
            for row in contributors
        )
    )
    for row in hdul["EXCLUDED"].data:
        print(f"excluded fiber {int(row['FIBER'])} ({row['TELESCOPE']}): {row['REASON']}")


def _show_master(hdul: fits.HDUList) -> None:
    header = hdul[0].header
    validity = {
        str(header[f"MPVALK{n:02d}"]): header[f"MPVALV{n:02d}"]
        for n in range(1, int(header["MPNVAL"]) + 1)
    }
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
