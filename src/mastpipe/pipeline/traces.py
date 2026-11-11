"""Stages 1b and 2 of the quick-look recipe: ``light`` frames resolved from their content
into arcs and science frames, and the traces found on the flat master or a full arc and
matched to the fiber layout. The whole recipe is described in
:mod:`mastpipe.pipeline.runner`.
"""

from __future__ import annotations

from dataclasses import replace
from functools import partial

from mastpipe.__version__ import __version__
from mastpipe.core.model import Frame, FrameRef
from mastpipe.io import products
from mastpipe.pipeline.group import (
    _calibrated,
    _check_binning,
    _check_readout,
    _Group,
    _master_meta,
    _per_frame,
    _same_readout,
    _skip,
    _time_key,
    _write_qc,
    log,
)
from mastpipe.pipeline.identify import CALIB_FOLDER
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.classify import classify_light
from mastpipe.steps.trace import find_traces, layout_indexed, match_layout, relative_throughput, trace_source


def _resolve_light(group: _Group, by_type: dict[str, list[FrameRef]]) -> None:
    """Stage 1b: every ``light`` frame, in time order and under the error policy, calibrated
    with the group's masters (no cosmic-ray flagging) and classified from its content: moved
    into ``by_type["arc"]`` or ``by_type["science"]`` with its ``classification_rule``, or
    skipped with its reason; a frame not binned or read out as the group fails
    (:func:`~mastpipe.pipeline.group._check_binning`,
    :func:`~mastpipe.pipeline.group._check_readout`). ``by_type`` loses its ``light``
    entry."""

    def resolve(group: _Group, ref: FrameRef) -> None:
        _check_binning(group, ref)
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


def _full_arc_candidates(group: _Group, arc_refs: list[FrameRef]) -> list[Frame]:
    """Calibrated arcs that may serve as the trace source, the one lighting the most traces
    first: the arcs the header calls full, else the arcs of unknown kind. A candidate not
    read out as the group is a calibration frame of another readout: recorded in
    ``skipped`` (:func:`~mastpipe.pipeline.group._same_readout`) and removed from
    ``arc_refs``. Each other candidate runs under the error policy of
    :func:`~mastpipe.pipeline.group._per_frame`: one not binned as the group, that cannot
    be calibrated or on which :func:`~mastpipe.steps.trace.find_traces` raises is recorded
    in ``failed`` and removed from ``arc_refs``, so it is withheld from the arc stage too;
    either way the arc stage does not record it again, and the others still compete."""
    known = [ref for ref in arc_refs if ref.meta.arc_kind == "full"]
    candidates = known or [ref for ref in arc_refs if ref.meta.arc_kind is None]
    trace_params = group.params.trace.model_copy(update={"collapse": "percentile"})
    scored: list[tuple[int, int, Frame]] = []
    same = _same_readout(group, candidates)
    for ref in candidates:
        if ref not in same:
            arc_refs.remove(ref)
    candidates = same

    def score(group: _Group, ref: FrameRef, n: int) -> None:
        _check_binning(group, ref)
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
    """Stage 2: find and match the traces, write the trace product and its QC. Without a
    flat master the full arcs compete as the trace source
    (:func:`_full_arc_candidates`), which removes from ``by_type["arc"]`` every candidate
    skipped for its readout or failed, so the arc stage never sees them."""
    params = group.params
    flat_frame: Frame | None = None
    flat_name: list[str] = []
    inputs: list[str] = []
    if "flat" in group.masters:
        master, path = group.masters["flat"]
        flat_frame = Frame(data=master.data, meta=_master_meta(group, master), variance=master.variance)
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
    matched = match_layout(find_traces(frame.data, trace_params), group.layout, params.trace.layout_tolerance)
    extras = [trace for trace in matched if trace.role is None]
    traceset = layout_indexed(matched)
    where = source.replace("_", " ")
    if len(traceset) == 0:
        raise RuntimeError(f"no trace of layout {group.preset} found on the {where} ({len(matched)} detected)")
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
    group.throughput = relative_throughput(flat_frame, traceset, params.extract.halfwidth) if flat_frame else {}
    log.info(
        "traces from the %s: %d of %d detections matched to layout %s, shift %s px, missing %s, %d extra dropped -> %s",
        where,
        len(traceset),
        len(matched),
        group.preset,
        "none" if traceset.layout_shift is None else f"{traceset.layout_shift:.2f}",
        list(traceset.missing),
        len(extras),
        path.name,
    )
