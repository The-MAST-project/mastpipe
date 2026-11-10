"""Stage 3 of the quick-look recipe: the night's arcs, a full arc re-solved into a new
reference solution and a calibration-fiber arc fitted for its drift (the time pairing of
exposures with drift arcs is :mod:`mastpipe.pipeline.arcs`). The whole recipe is described
in :mod:`mastpipe.pipeline.runner`.
"""

from __future__ import annotations

import numpy as np

from mastpipe.core.model import FrameRef, TraceSet
from mastpipe.io import products
from mastpipe.pipeline.group import (
    _calibrated,
    _check_binning,
    _check_readout,
    _configuration_origin,
    _configuration_steps,
    _content_steps,
    _Group,
    _per_frame,
    _reference_for,
    _skip,
    _time_key,
    _write_qc,
    log,
)
from mastpipe.pipeline.identify import _arc_lines, _provenance, _solve_and_package
from mastpipe.steps.drift import drift_unusable_reason, fit_drift
from mastpipe.steps.trace import classify_arc


def _arcs(group: _Group, arc_refs: list[FrameRef]) -> None:
    """Stage 3: every arc in time order, under the error policy."""
    for ref in sorted(arc_refs, key=_time_key):
        _per_frame(group, ref, _arc)


def _arc(group: _Group, ref: FrameRef) -> None:
    _check_binning(group, ref)
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
        [*_content_steps(meta), *_configuration_steps(group, meta), f"wavecal seeded by reference {previous.identifier}"],
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
            "configuration_origin": _configuration_origin(group, meta),
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
            "configuration_origin": _configuration_origin(group, meta),
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
