"""The quick-look recipe: a night folder to masters, traces, references, drift models,
wavelength-calibrated and sky-subtracted spectra and per-target combinations, with QC.

:func:`run_quicklook` scans the folder and its subfolders (manifest applied; frames of
another instrument and manifest rows matching no frame go into ``skipped``), drops the
frames whose type is ``unknown`` into ``skipped`` with their classification rule, gives
``options.configuration`` to every frame whose header and manifest name no grating (so its
products, masters and library entries carry that name), and reduces every
``(camera, configuration)`` group in turn:

1. **masters** — bias, dark and background frames are stacked (dark and background one
   exposure time at a time, the largest batch serving the night), the flats calibrated
   and stacked, and every master is added to the calibration library (``out/calib`` or
   ``options.library``) — unless the library already holds one stacked from the same
   frames under the same stacking parameters and detector (``calib_hash``; and, for the
   flat, calibrated with the same masters), which is reused rather than rebuilt and
   overwritten; ``options.force`` rebuilds. The frames are loaded one at a time: one that
   cannot be read, is not binned as the group's science frames (the
   most common binning among them) or has another shape than the detector's (after crop
   and orientation) is recorded in ``failed`` and left out; one read out otherwise than the
   group (another readout mode or amplifier) is recorded in ``skipped`` and left out,
   while a science, ``light`` or arc frame of another binning or readout fails. A kind the night lacks, or
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

import logging
import re
from collections import defaultdict
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from mastpipe.__version__ import __version__
from mastpipe.config.loader import build_effective, calib_hash, config_hash
from mastpipe.config.models import InstrumentConfig
from mastpipe.core.model import FrameRef
from mastpipe.instruments.registry import layout_hash, load_header_map, load_instrument, load_layout
from mastpipe.io.ingest import ASSUMED_UTC, apply_manifest, read_manifest, scan_folder
from mastpipe.pipeline.arc_stage import _arcs
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.group import (
    LOG_FOLDER,
    NO_FRAMES,
    SPEC_FOLDER,
    TIME_POLICIES,
    RunOptions,
    RunSummary,
    _describe,
    _Group,
    _group_binning,
    _group_readout,
    _per_frame,
    _time_key,
    _timezone,
    _unlisted_amplifiers,
    log,
)
from mastpipe.pipeline.identify import CALIB_FOLDER
from mastpipe.pipeline.masters import _masters
from mastpipe.pipeline.science import _science_frame
from mastpipe.pipeline.traces import _resolve_light, _traces

__all__ = [
    "EARLIER_PRODUCT",
    "NO_FRAMES",
    "TIME_POLICIES",
    "RunOptions",
    "RunSummary",
    "run_quicklook",
]


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
    temperatures = [m.detector_temperature for m in (r.meta for r in refs) if m.detector_temperature is not None]
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


EARLIER_PRODUCT = re.compile(r"_(?:\d{8}T\d{6}|notime)_(?:spec|comb_[A-Za-z0-9_-]*)\.fits$")
"""A science or combined product name without the checksum segment after the stamp (an
earlier version's naming); :func:`~mastpipe.io.products.product_name` always puts a six
character checksum or ``nosum`` there."""


def _warn_earlier_products(spec_dir: Path) -> None:
    """Warn once when ``spec_dir`` holds products named by an earlier version: they are
    never found at their new names, so the run writes new files beside them."""
    if not spec_dir.is_dir():
        return
    earlier = sorted(path.name for path in spec_dir.iterdir() if EARLIER_PRODUCT.search(path.name))
    if earlier:
        log.warning(
            "%s holds %d product file(s) written by an earlier mastpipe version (no checksum in the "
            "name, e.g. %s); they are not replaced and the new products are written beside them: "
            "use a fresh --out, or clear its spec and qa folders",
            spec_dir,
            len(earlier),
            earlier[0],
        )


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
        _warn_earlier_products(out_dir / SPEC_FOLDER)
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
            if meta.configuration is None and default_conf is not None:
                # the run's --configuration names the grating of a frame whose header does
                # not, so its products, masters and library entries carry that name too
                ref = replace(ref, meta=replace(meta, configuration=default_conf))
            groups[(meta.camera, ref.meta.configuration)].append(ref)
        for (camera, configuration), group_refs in sorted(groups.items(), key=str):
            label = f"group {camera}/{configuration}"
            if options.camera is not None and camera != options.camera:
                summary.skipped.append(f"{label}: not camera {options.camera}")
                log.info("skipped %s: not camera %s", label, options.camera)
                continue
            if configuration is None:
                choices = sorted(instrument.configurations)
                summary.failed.append((label, f"configuration unknown; pass --configuration one of {choices}"))
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
