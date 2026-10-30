"""Pairing a science exposure with the calibration-fiber arcs that bracket it in time.

A night's drift arcs each carry a :class:`~mastpipe.core.model.DriftModel` against the
reference solution; which of them a science exposure uses, and with what weights, is the
time policy of :class:`~mastpipe.config.models.DriftParams`. :func:`pair_arcs` applies it
to the frame metadata alone (no file access): candidates are the arcs of the exposure's
configuration whose drift fit is usable, taken no earlier than the reference's
``valid_from`` and within ``max_age_hours`` of the exposure. ``nearest_before`` takes the
latest arc before the exposure (the nearest after when none precedes it), ``nearest`` the
closest in either direction and ``interpolate`` the arcs before and after weighted by
their distance in time (the nearest alone when only one side exists). Times are ISO-8601
strings; a naive time is UTC. :func:`status_for` names the wavelength status of the
outcome: ``"transferred"`` with at least one arc, ``"reference"`` without.
"""

from __future__ import annotations

import pandas as pd

from mastpipe.config.models import DriftParams
from mastpipe.core.model import DriftModel, FrameMeta, ReferenceSolution

SECONDS_PER_HOUR = 3600.0


def _timestamp(value: str) -> pd.Timestamp:
    """``value`` as a UTC timestamp; a naive time is taken to be UTC."""
    stamp = pd.Timestamp(value)
    if stamp.tzinfo is None:
        return stamp.tz_localize("UTC")
    return stamp.tz_convert("UTC")


def pair_arcs(
    exposure: FrameMeta,
    arcs: list[tuple[FrameMeta, DriftModel]],
    reference: ReferenceSolution,
    params: DriftParams,
) -> list[tuple[DriftModel, float]]:
    """The drift models an exposure uses, with their weights, under ``params.time_policy``.

    An arc is a candidate when its ``configuration`` is the exposure's, its model is
    ``usable``, its mid-time is not before ``reference.valid_from`` and it lies within
    ``params.max_age_hours`` of the exposure's mid-time. Weights sum to one; the list is
    empty when nothing qualifies (or the exposure has no mid-time).
    """
    if exposure.mid_time is None:
        return []
    t_exposure = _timestamp(exposure.mid_time)
    valid_from = _timestamp(reference.valid_from)
    candidates: list[tuple[float, DriftModel]] = []
    for meta, model in arcs:
        if meta.configuration != exposure.configuration or not model.usable:
            continue
        if meta.mid_time is None:
            continue
        t_arc = _timestamp(meta.mid_time)
        if t_arc < valid_from:
            continue
        dt_hours = (t_arc - t_exposure).total_seconds() / SECONDS_PER_HOUR
        if abs(dt_hours) > params.max_age_hours:
            continue
        candidates.append((dt_hours, model))
    if not candidates:
        return []

    nearest = min(candidates, key=lambda c: abs(c[0]))
    if params.time_policy == "nearest":
        return [(nearest[1], 1.0)]
    before = [c for c in candidates if c[0] <= 0.0]
    after = [c for c in candidates if c[0] > 0.0]
    latest_before = max(before, key=lambda c: c[0]) if before else None
    first_after = min(after, key=lambda c: c[0]) if after else None
    if params.time_policy == "nearest_before":
        chosen = latest_before if latest_before is not None else first_after
        assert chosen is not None
        return [(chosen[1], 1.0)]
    # interpolate
    if latest_before is None or first_after is None:
        return [(nearest[1], 1.0)]
    dt_before, dt_after = latest_before[0], first_after[0]
    fraction = -dt_before / (dt_after - dt_before)
    if fraction == 0.0:
        return [(latest_before[1], 1.0)]
    return [(latest_before[1], 1.0 - fraction), (first_after[1], fraction)]


def status_for(pairs: list[tuple[DriftModel, float]]) -> str:
    """The wavelength status an exposure gets from its pairing."""
    return "transferred" if pairs else "reference"
