"""Stage 1 of the quick-look recipe: the night's bias, dark, background and flat masters,
stacked into the calibration library or reused from it, and the group's calibration set.
The whole recipe is described in :mod:`mastpipe.pipeline.runner`.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

from astropy.io import fits

from mastpipe.core.hashing import HASH_LENGTH
from mastpipe.core.model import Frame, FrameRef, MasterCalib
from mastpipe.io.ingest import expected_shape, file_checksum, load_frame
from mastpipe.pipeline.group import (
    MASTER_KINDS,
    _by_exposure,
    _check_binning,
    _Group,
    _master_provenance,
    _per_frame,
    _same_readout,
    log,
)
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame, calibration_steps


def _usable_frames(group: _Group, refs: list[FrameRef], calibrate: bool) -> list[Frame]:
    """The frames of ``refs`` that can go into one master, loaded one at a time (and, when
    ``calibrate``, calibrated as a flat is: background subtracted when
    ``params.calib.flat_subtract_background`` and the set has a background with an
    exposure time). A frame that fails to
    load or calibrate, whose binning is not ``group.binning`` or whose loaded shape is not
    the detector's at that binning (:func:`~mastpipe.io.ingest.expected_shape`) is dropped
    under the error policy of :func:`~mastpipe.pipeline.group._per_frame`: recorded in
    ``failed`` with the reason, or raised under ``on_error='raise'``. A frame whose readout key is not
    ``group.readout_key`` is left out first and recorded in ``skipped``
    (:func:`~mastpipe.pipeline.group._same_readout`)."""
    detector = group.effective.detector
    shape = expected_shape(detector, group.binning)
    frames: list[Frame] = []

    def use(group: _Group, ref: FrameRef) -> None:
        _check_binning(group, ref)
        frame = load_frame(ref, detector)
        if frame.data.shape != shape:
            raise ValueError(f"shape {frame.data.shape} is not the detector's {shape}")
        if calibrate:
            calib = group.calib
            assert calib is not None
            subtract = group.params.calib.flat_subtract_background and calib.has_background
            frame = calibrate_frame(frame, calib, subtract_background=subtract, saturation=detector.saturation)
        frames.append(frame)

    for ref in _same_readout(group, refs):
        _per_frame(group, ref, use)
    return frames


def _masters_hash(masters: dict[str, tuple[MasterCalib, Path]]) -> str | None:
    """SHA-256 over the sorted file checksums of ``masters``' files (first
    :data:`~mastpipe.core.hashing.HASH_LENGTH` hex characters), part of a science
    product's idempotence key; ``None`` without masters."""
    if not masters:
        return None
    checksums = sorted(file_checksum(path) for _, path in masters.values())
    return hashlib.sha256("\n".join(checksums).encode()).hexdigest()[:HASH_LENGTH]


def _batch_score(frames: list[Frame]) -> tuple[int, float]:
    """How well a batch of surviving frames (one exposure time) serves the night: the
    number of frames, then the exposure time."""
    return len(frames), float(frames[0].meta.exposure_time or 0.0)


def _library_master(
    group: _Group, kind: str, frames: list[Frame], masters_hash: str | None
) -> tuple[MasterCalib, Path] | None:
    """The library's master of ``kind`` stacked from exactly ``frames`` (by checksum) under
    the group's ``calib_hash`` and calibrated with masters of exactly the content
    ``masters_hash`` says (:func:`_masters_hash` of the masters applied, the file's
    ``MPMKH``; ``None`` for a master stacked from raw frames), to reuse instead of
    building it again; ``None`` under ``options.force`` or without one. Comparing content
    rather than file names catches a master rewritten under the same name; a flat written
    before ``MPMKH`` was recorded reads as ``None`` and is built again once."""
    if group.options.force:
        return None
    checksums = [frame.meta.checksum for frame in frames]
    found = group.library.find_by_inputs(kind, checksums, group.effective.camera, group.calib_hash)
    if found is None:
        return None
    stored = Provenance.from_header(fits.getheader(found[1]))
    if stored.masters_hash != masters_hash:
        log.info(
            "master %s in the library was calibrated with other masters (%s, hash %s, now %s); building it again",
            kind,
            ", ".join(stored.calibrations) or "none",
            stored.masters_hash,
            masters_hash,
        )
        return None
    log.info("master %s reused from the library: %s", kind, found[1].name)
    return found


def _masters(group: _Group, by_type: dict[str, list[FrameRef]]) -> None:
    """Stage 1: stack the night's masters into the library and build the calibration set.
    Every frame goes through :func:`_usable_frames`; a master the library already holds,
    stacked from the same surviving frames under the same ``calib_hash`` (and, for the
    flat, calibrated with masters of the same content, :func:`_library_master`, none of
    them built in this run), is reused rather than built and overwritten unless
    ``options.force``; a kind none of whose frames survive is sought in the library.
    Every kind is batched by exposure time (the bias in one batch) and the batch with the
    most surviving frames serves the night (the longer exposure on a tie); the dark and
    background batches are each kept in the library, the flat's losing batches are not
    stacked. Sets ``group.masters`` and ``group.masters_hash``."""
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
            found = _library_master(group, kind, frames, None)
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
            score = _batch_score(frames)
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
    applied_hash = _masters_hash(chosen)

    flats = by_type.get("flat", [])
    frames: list[Frame] = []
    for batch in _by_exposure(flats):
        survivors = _usable_frames(group, batch, calibrate=True)
        if survivors and (not frames or _batch_score(survivors) > _batch_score(frames)):
            frames = survivors
    found = None if not frames or built else _library_master(group, "flat", frames, applied_hash)
    if found is not None:
        chosen["flat"] = found
    elif frames:
        subtract = params.calib.flat_subtract_background and group.calib.has_background
        master = build_master(frames, "flat", params.calib)
        steps = calibration_steps(group.calib, float(frames[0].meta.exposure_time or 0.0), subtract)
        steps.append(f"stack {params.calib.stack_method} of {len(frames)} flat frames")
        prov = _master_provenance(group, frames, applied, steps, masters_hash=applied_hash)
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
