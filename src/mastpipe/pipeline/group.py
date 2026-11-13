"""The state shared by the quick-look stages: the run options and summary, the state of
one ``(camera, configuration)`` group (:class:`_Group`), the output folder and policy
constants, the ``mastpipe`` logger, and the helpers shared by the stages and the runner
(error policy, binning and readout checks, calibration of one frame, QC writing).

The stage modules (:mod:`~mastpipe.pipeline.masters`, :mod:`~mastpipe.pipeline.traces`,
:mod:`~mastpipe.pipeline.arc_stage`, :mod:`~mastpipe.pipeline.science`) import this module
and never each other. The whole recipe is described in :mod:`mastpipe.pipeline.runner`.
"""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Literal

import numpy as np
import pandas as pd

from mastpipe.__version__ import __version__
from mastpipe.config.models import DriftParams, EffectiveConfig, InstrumentConfig, StepParams
from mastpipe.core.errors import RunError
from mastpipe.core.model import (
    CalibrationSet,
    DriftModel,
    Frame,
    FrameMeta,
    FrameRef,
    MasterCalib,
    QCReport,
    ReferenceSolution,
    Spectrum,
    TraceSet,
)
from mastpipe.io import products
from mastpipe.io.ingest import load_frame
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.calib import calibrate_frame, calibration_steps
from mastpipe.steps.classify import ARC_LIKE, CONTINUUM
from mastpipe.steps.cosmics import flag_cosmic_rays
from mastpipe.steps.geometry import GratingGeometry, grating_geometry

log = logging.getLogger("mastpipe")

SPEC_FOLDER = "spec"
QC_FOLDER = "qa"
LOG_FOLDER = "logs"
MASTER_KINDS = ("bias", "dark", "background")
"""The master kinds stacked from raw frames, in the order they are built."""
TIME_POLICIES = ("nearest_before", "nearest", "interpolate")
ALREADY_REDUCED = "already reduced"
"""Start of the ``skipped`` reason of a science product whose idempotence key is unchanged."""
CONFIGURATION_FROM_RUN = "run"
"""QC ``configuration_origin`` of a frame given ``options.configuration`` by the run."""
CONFIGURATION_FROM_FRAME = "frame"
"""QC ``configuration_origin`` of a frame whose header, file name or manifest named its grating."""


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
    ``out_dir/calib``); ``quiet`` keeps the console to warnings and errors (the run's log
    file still has every ``INFO`` line).
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
    quiet: bool = False


Scope = Literal["frame", "group", "manifest", "run"]
"""What a summary item is about: one frame (``source`` its path relative to the night
folder), one ``(camera, configuration)`` group (``source`` is ``<camera>/<configuration>``),
one manifest row (``source`` its source cell) or the whole run (``source`` is ``None``)."""


def _label(scope: Scope, source: str | None) -> str:
    if scope == "run":
        return "run"
    if scope == "group":
        return f"group {source}"
    return str(source)


@dataclass(frozen=True)
class Skipped:
    """Something a run left out on purpose, and why."""

    scope: Scope
    source: str | None
    reason: str

    @property
    def label(self) -> str:
        """``source`` for a frame or manifest row, ``group <source>`` for a group, ``run``
        for the run."""
        return _label(self.scope, self.source)


@dataclass(frozen=True)
class Failure:
    """Something a run could not do: ``message`` says what went wrong, ``kind`` sorts the
    failure (``"error"`` unless a site names it), ``error_type`` names the exception's
    class when an exception caused it."""

    scope: Scope
    source: str | None
    message: str
    kind: str = "error"
    error_type: str | None = None

    @property
    def label(self) -> str:
        """As :attr:`Skipped.label`."""
        return _label(self.scope, self.source)


@dataclass
class RunSummary:
    """What a run produced: product and QC paths, the things it skipped
    (:class:`Skipped`) or failed on (:class:`Failure`), the sources of the science products
    it reduced again, and the run's facts: its start and end (UTC,
    ``YYYY-MM-DDTHH:MM:SSZ``), the number of frames the scan found (after the manifest),
    and its log file and calibration library."""

    products: list[Path] = field(default_factory=list)
    skipped: list[Skipped] = field(default_factory=list)
    failed: list[Failure] = field(default_factory=list)
    qc: list[Path] = field(default_factory=list)
    reduced_again: list[str] = field(default_factory=list)
    started: str | None = None
    finished: str | None = None
    frames_scanned: int | None = None
    log: Path | None = None
    library: Path | None = None


def stopped_by_run_error(summary: RunSummary, error: BaseException | None) -> bool:
    """Whether ``error`` is a run error of the result (``failed``, exit 4, with its own
    ``kind`` and ``retryable``): a :class:`~mastpipe.core.errors.RunError` raised before
    the calibration library was open (``summary.library`` still ``None``), when nothing
    has been reduced. Once the library is open, any exception leaving the run -- a
    ``RunError`` re-raised from a group under ``--on-error raise`` included -- is an
    unexpected ``error``. One predicate for the exit code and the result document
    (:mod:`mastpipe.pipeline.summary`), stderr (:mod:`mastpipe.pipeline.commands`) and the
    log (:mod:`mastpipe.pipeline.runner`)."""
    return isinstance(error, RunError) and summary.library is None


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
    calib_hash: str
    """:func:`~mastpipe.config.loader.calib_hash` of ``effective``: masters are reused under it."""
    timezone: str | None = None
    """How the group's frame times were read (:func:`_timezone` of their metadata)."""
    layout_hash: str | None = None
    """:func:`~mastpipe.instruments.registry.layout_hash` of the group's preset."""
    masters: dict[str, tuple[MasterCalib, Path]] = field(default_factory=dict)
    masters_hash: str | None = None
    """:func:`~mastpipe.pipeline.masters._masters_hash` of ``masters``, set once they are chosen."""
    calib: CalibrationSet | None = None
    traceset: TraceSet | None = None
    traces_path: Path | None = None
    throughput: dict[int, float] = field(default_factory=dict)
    arcs: list[tuple[FrameMeta, DriftModel, str]] = field(default_factory=list)
    """The night's drift models with the identifier of the reference each was fitted against."""
    run_configured: frozenset[tuple[str, int]] = frozenset()
    """``(source, extension)`` of the frames whose header and manifest named no grating and
    that the run gave ``options.configuration`` (:func:`_configuration_steps`)."""

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


def _time_key(ref: FrameRef) -> tuple[bool, str, str]:
    return (ref.meta.mid_time is None, ref.meta.mid_time or "", ref.meta.source)


def _record_skip(summary: RunSummary, scope: Scope, source: str | None, reason: str, level: int) -> None:
    item = Skipped(scope, source, reason)
    summary.skipped.append(item)
    log.log(level, "skipped %s: %s", item.label, reason)


def skip(summary: RunSummary, scope: Scope, source: str | None, reason: str) -> None:
    """Record in ``summary`` that ``source`` was skipped for ``reason``, with a warning."""
    _record_skip(summary, scope, source, reason, logging.WARNING)


def fail(
    summary: RunSummary,
    scope: Scope,
    source: str | None,
    error: BaseException | str,
    *,
    prefix: str = "",
    kind: str | None = None,
) -> None:
    """Record in ``summary`` that ``source`` failed, and log it (with the traceback when
    ``error`` is an exception). An exception gives the message ``prefix + str(error)``
    (its class name when ``str(error)`` is empty) and its class name as ``error_type``; a
    string is the message itself, without ``prefix``. An explicit ``kind`` wins; else a
    :class:`~mastpipe.core.errors.RunError` gives its own ``kind``; else ``"error"``."""
    if kind is None:
        kind = error.kind if isinstance(error, RunError) else "error"
    if isinstance(error, BaseException):
        message = prefix + (str(error) or type(error).__name__)
        item = Failure(scope, source, message, kind, type(error).__name__)
    else:
        item = Failure(scope, source, error, kind)
    summary.failed.append(item)
    log.error(
        "failed %s: %s",
        item.label,
        item.message,
        exc_info=error if isinstance(error, BaseException) else None,
    )


def _skip(group: _Group, source: str, reason: str) -> None:
    """Record that the frame ``source`` was skipped for ``reason`` (logged at INFO)."""
    _record_skip(group.summary, "frame", source, reason, logging.INFO)


def _per_frame(group: _Group, ref: FrameRef, action: Callable[[_Group, FrameRef], None]) -> None:
    """Run ``action`` on one frame under the error policy."""
    try:
        action(group, ref)
    except Exception as error:
        if group.options.on_error == "raise":
            raise
        fail(group.summary, "frame", ref.meta.source, error)


def _by_exposure(refs: list[FrameRef]) -> list[list[FrameRef]]:
    """``refs`` batched by exposure time, shortest first (frames without one last)."""
    batches: dict[float | None, list[FrameRef]] = {}
    for ref in refs:
        batches.setdefault(ref.meta.exposure_time, []).append(ref)
    return [batches[t] for t in sorted(batches, key=lambda t: (t is None, t or 0.0))]


def _median_snr(spectrum: Spectrum) -> float | None:
    good = ~spectrum.mask & np.isfinite(spectrum.flux) & np.isfinite(spectrum.error) & (spectrum.error > 0)
    if not good.any():
        return None
    return float(np.median(spectrum.flux[good] / spectrum.error[good]))


def _write_qc(group: _Group, name: str, stage: str, values: dict[str, Any]) -> Path:
    path = products.write_qc(group.qc_dir / name, QCReport(stage, values))
    group.summary.qc.append(path)
    return path


def _master_meta(group: _Group, master: MasterCalib) -> FrameMeta:
    """A :class:`~mastpipe.core.model.FrameMeta` standing for a master (from the library,
    whose frames are not at hand), for naming and for the frames built from it. Its
    ``checksum`` is the master's first input frame's, the checksum the master's own name
    carries, so a product named after it (the trace set of a flat) carries it too; empty
    (``nosum`` in a name) for a master that lists no inputs."""
    validity = master.validity
    inputs = master.provenance.get("inputs") or []
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
        checksum=str(inputs[0]) if inputs else "",
    )


def _calibrated(group: _Group, ref: FrameRef, cosmics: bool = False) -> tuple[Frame, list[str]]:
    """``ref`` loaded and calibrated with the group's masters (background subtracted when
    ``params.calib.subtract_background`` and the set has a background master with an
    exposure time to scale it by, :attr:`~mastpipe.core.model.CalibrationSet.has_background`;
    the bias/dark path otherwise), cosmic rays flagged into the mask when ``cosmics``; with
    the steps applied."""
    detector, params = group.effective.detector, group.params
    calib = group.calib
    assert calib is not None
    frame = load_frame(ref, detector)
    subtract = params.calib.subtract_background and calib.has_background
    out = calibrate_frame(frame, calib, subtract_background=subtract, saturation=detector.saturation)
    steps = calibration_steps(calib, float(ref.meta.exposure_time or 0.0), subtract)
    if cosmics:
        assert out.variance is not None
        before = 0 if out.mask is None else int(out.mask.sum())
        mask = flag_cosmic_rays(out.data, out.variance, out.mask, params.calib.cosmic_sigma, params.calib.cosmic_size)
        out = Frame(data=out.data, meta=out.meta, variance=out.variance, mask=mask)
        steps.append(
            f"cosmic rays: {int(mask.sum()) - before} px flagged "
            f"(sigma {params.calib.cosmic_sigma:g}, size {params.calib.cosmic_size})"
        )
    return out, steps


def _timezone(metas: Iterable[FrameMeta]) -> str | None:
    """The ``FrameMeta.timezone`` labels of ``metas`` for ``Provenance.timezone``: the one
    label they share, the distinct labels joined by commas when they differ, ``None`` when
    no frame has a time."""
    zones = sorted({meta.timezone for meta in metas if meta.timezone is not None})
    return ",".join(zones) if zones else None


def _master_provenance(
    group: _Group,
    frames: list[Frame],
    calibrations: list[str],
    steps: list[str],
    masters_hash: str | None = None,
) -> Provenance:
    """The provenance of a master stacked from ``frames``; ``masters_hash`` is the
    :func:`~mastpipe.pipeline.masters._masters_hash` of the masters applied to them (a
    flat's), ``None`` for one stacked from raw frames."""
    return Provenance(
        pipeline_version=__version__,
        config_hash=group.cfg_hash,
        inputs=[frame.meta.checksum for frame in frames],
        calibrations=calibrations,
        steps=steps,
        timezone=_timezone(frame.meta for frame in frames),
        masters_hash=masters_hash,
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


def _readout_mismatch(group: _Group, ref: FrameRef) -> str | None:
    """Why ``ref`` was not read out as the group (another readout mode or output amplifier:
    :attr:`~mastpipe.core.model.FrameMeta.readout_key`); ``None`` when it was."""
    if ref.meta.readout_key != group.readout_key:
        return f"readout {ref.meta.readout_key} is not the group's {group.readout_key}"
    return None


def _check_readout(group: _Group, ref: FrameRef) -> None:
    """``ValueError`` when ``ref`` was not read out as the group (:func:`_readout_mismatch`):
    for science, ``light`` and arc frames, which fail."""
    reason = _readout_mismatch(group, ref)
    if reason is not None:
        raise ValueError(reason)


def _same_readout(group: _Group, refs: list[FrameRef]) -> list[FrameRef]:
    """The calibration frames of ``refs`` read out as the group; each other one is recorded
    in ``skipped`` with :func:`_readout_mismatch`'s reason, not as a failure (a calibration
    frame of another readout belongs to another set, not to a broken one)."""
    kept = []
    for ref in refs:
        reason = _readout_mismatch(group, ref)
        if reason is None:
            kept.append(ref)
        else:
            _skip(group, ref.meta.source, reason)
    return kept


def _content_steps(meta: FrameMeta) -> list[str]:
    """The provenance step recording a type resolved from the frame's content, if it was."""
    if meta.classification_rule in (ARC_LIKE, CONTINUUM):
        return [f"exposure type {meta.exposure_type} from the frame's content ({meta.classification_rule})"]
    return []


def _configuration_origin(group: _Group, meta: FrameMeta) -> str:
    """Where the frame's configuration came from, for its QC: :data:`CONFIGURATION_FROM_RUN`
    (``options.configuration``) or :data:`CONFIGURATION_FROM_FRAME` (its header, file name or
    manifest)."""
    configured = (meta.source, meta.extension) in group.run_configured
    return CONFIGURATION_FROM_RUN if configured else CONFIGURATION_FROM_FRAME


def _configuration_steps(group: _Group, meta: FrameMeta) -> list[str]:
    """The provenance step recording a configuration the run gave the frame, if it did."""
    if _configuration_origin(group, meta) == CONFIGURATION_FROM_RUN:
        origin = "from the run's --configuration (none in the header, file name or manifest)"
        return [f"configuration {meta.configuration} {origin}"]
    return []


def _reference_for(group: _Group, at_time: str | None) -> ReferenceSolution | None:
    """The reference a frame at ``at_time`` uses: the pinned ``options.reference_id``
    (``KeyError`` when the library lacks it), else the one valid at that time (the latest
    overall for a frame without a time)."""
    effective, library = group.effective, group.library
    keys = (effective.instrument, effective.configuration, group.preset)
    if group.options.reference_id is not None:
        ref = library.find_reference(*keys, identifier=group.options.reference_id)
        if ref is None:
            raise KeyError(f"reference {group.options.reference_id!r} is not in the library {library.root}")
        return ref
    return library.find_reference(*keys, at_time=at_time)
