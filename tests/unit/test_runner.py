"""Unit tests of the quick-look runner's pure helpers: the drift models a frame may be
paired with (those fitted against its own reference) and the blend of two drift models
under the ``interpolate`` policy."""

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.models import DriftParams
from mastpipe.core.model import PER_FIBER_COLUMNS, DriftModel, FrameMeta, ReferenceSolution
from mastpipe.pipeline.arcs import pair_arcs
from mastpipe.pipeline.runner import _arc_weights, _arcs_for, _blend_drifts

N_PIX = 1598


def make_meta(mid_time: str, checksum: str, exposure_type: str = "arc") -> FrameMeta:
    return FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type=exposure_type,
        configuration="mgb518",
        mid_time=mid_time,
        checksum=checksum,
        source=f"{checksum}.fits",
    )


def make_reference(identifier: str, valid_from: str) -> ReferenceSolution:
    return ReferenceSolution(
        instrument="highspec",
        configuration="mgb518",
        fiber_preset="survey45",
        solutions={},
        created=valid_from,
        valid_from=valid_from,
        source_checksums=(identifier,),
        identifier=identifier,
        hash="h" + identifier,
    )


def make_drift(
    checksum: str,
    mid_time: str,
    shift_coeffs: list[float],
    stretch_coeffs: list[float],
    fibers: dict[int, tuple[float, float, int, float, bool]],
    scatter_px: float = 0.02,
    usable: bool = True,
) -> DriftModel:
    rows = [
        {
            "fiber": fiber,
            "y_center": 10.0 + 8.5 * fiber,
            "shift_px": shift,
            "stretch": stretch,
            "n_lines": n_lines,
            "rms_px": rms,
            "usable": ok,
        }
        for fiber, (shift, stretch, n_lines, rms, ok) in fibers.items()
    ]
    return DriftModel(
        arc_checksums=(checksum,),
        arc_mid_time=mid_time,
        per_fiber=pd.DataFrame(rows, columns=list(PER_FIBER_COLUMNS)),
        shift_coeffs=np.array(shift_coeffs),
        stretch_coeffs=np.array(stretch_coeffs),
        scatter_px=scatter_px,
        usable=usable,
        n_pix=N_PIX,
    )


# ------------------------------------------------------------- arcs tied to their reference


@pytest.mark.parametrize("policy", ["nearest", "interpolate", "nearest_before"])
def test_frame_never_receives_drift_fitted_against_another_reference(policy):
    """A full arc solved mid-night makes reference B; the arcs after it are fitted against
    B. A science frame before that time uses A and must pair with A's arcs only, however
    close B's arcs are in time."""
    ref_a = make_reference("A", "2026-10-01T18:00:00")
    ref_b = make_reference("B", "2026-10-01T21:00:00")
    fibers = {0: (0.3, 0.0, 20, 0.02, True), 44: (0.3, 0.0, 20, 0.02, True)}
    arc_a1 = make_drift("a1", "2026-10-01T19:00:00", [0.3], [0.0], fibers)
    arc_b1 = make_drift("b1", "2026-10-01T21:10:00", [-2.0], [0.0], fibers)
    arc_b2 = make_drift("b2", "2026-10-01T21:30:00", [-2.0], [0.0], fibers)
    arcs = [
        (make_meta("2026-10-01T19:00:00", "a1"), arc_a1, "A"),
        (make_meta("2026-10-01T21:10:00", "b1"), arc_b1, "B"),
        (make_meta("2026-10-01T21:30:00", "b2"), arc_b2, "B"),
    ]
    # the frame sits just before B's validity, between an A arc and the B arcs
    frame = make_meta("2026-10-01T20:55:00", "sci", exposure_type="science")
    params = DriftParams(time_policy=policy)

    mine = _arcs_for(arcs, ref_a)
    assert [meta.checksum for meta, _ in mine] == ["a1"]
    pairs = pair_arcs(frame, mine, ref_a, params)
    assert [model.arc_checksums for model, _ in pairs] == [("a1",)]
    assert pairs[0][1] == pytest.approx(1.0)
    # without the filter the nearest / interpolate policies would have taken the B arc
    # fifteen minutes away (nearest_before alone is safe: no B arc precedes the frame)
    unfiltered = pair_arcs(frame, [(m, d) for m, d, _ in arcs], ref_a, params)
    taken = sorted(c for model, _ in unfiltered for c in model.arc_checksums)
    assert taken == (
        ["a1"]
        if policy == "nearest_before"
        else ["a1", "b1"]
        if policy == "interpolate"
        else ["b1"]
    )
    # a frame after B's time uses B and only B's arcs
    later = make_meta("2026-10-01T21:20:00", "sci2", exposure_type="science")
    theirs = _arcs_for(arcs, ref_b)
    assert sorted(meta.checksum for meta, _ in theirs) == ["b1", "b2"]
    pairs = pair_arcs(later, theirs, ref_b, params)
    assert pairs and all(c.startswith("b") for model, _ in pairs for c in model.arc_checksums)


# --------------------------------------------------------------------------- the blend


def test_blend_drifts_single_model_is_itself():
    drift = make_drift("x", "2026-10-01T20:00:00", [0.1], [0.0], {0: (0.1, 0.0, 10, 0.01, True)})
    blended, weights = _blend_drifts([(drift, 1.0)])
    assert blended is drift and weights is None


def test_blend_drifts_weights_two_models():
    fibers_a = {0: (0.2, 1e-5, 20, 0.02, True), 20: (0.4, 1e-5, 18, 0.03, True)}
    fibers_b = {0: (1.2, 3e-5, 15, 0.04, True), 44: (1.6, 3e-5, 16, 0.01, False)}
    first = make_drift("aaa", "2026-10-01T20:00:00", [0.2], [1e-5], fibers_a, scatter_px=0.02)
    second = make_drift(
        "bbb", "2026-10-01T21:00:00", [2e-3, 1.2], [0.0, 3e-5], fibers_b, scatter_px=0.05
    )
    blended, weights = _blend_drifts([(first, 0.3), (second, 0.7)])
    for y in (10.0, 300.0):
        assert blended.shift_at(y) == pytest.approx(
            0.3 * first.shift_at(y) + 0.7 * second.shift_at(y)
        )
        assert blended.stretch_at(y) == pytest.approx(
            0.3 * first.stretch_at(y) + 0.7 * second.stretch_at(y)
        )
    # right-aligned coefficients: the constant model contributes to the last slot only
    assert blended.shift_coeffs == pytest.approx([0.7 * 2e-3, 0.3 * 0.2 + 0.7 * 1.2])
    assert blended.stretch_coeffs == pytest.approx([0.0, 0.3 * 1e-5 + 0.7 * 3e-5])
    assert blended.arc_checksums == ("aaa", "bbb")
    assert blended.arc_mid_time == "2026-10-01T20:42:00"  # 0.7 of the hour after 20:00
    assert blended.n_pix == N_PIX
    assert blended.scatter_px == pytest.approx(0.05)
    assert blended.usable is True
    table = blended.per_fiber.set_index("fiber")
    assert sorted(table.index) == [0, 20, 44]
    assert table.loc[0, "shift_px"] == pytest.approx(0.3 * 0.2 + 0.7 * 1.2)
    assert table.loc[0, "stretch"] == pytest.approx(0.3 * 1e-5 + 0.7 * 3e-5)
    assert table.loc[0, "n_lines"] == 15 and table.loc[0, "usable"]
    assert table.loc[0, "rms_px"] == pytest.approx(0.3 * 0.02 + 0.7 * 0.04)
    # a fiber only one model holds takes that model's values and its weight
    assert table.loc[20, "shift_px"] == pytest.approx(0.4)
    assert table.loc[44, "shift_px"] == pytest.approx(1.6) and not table.loc[44, "usable"]
    assert weights == {0: pytest.approx(1.0), 20: pytest.approx(0.3), 44: pytest.approx(0.7)}
    # unnormalised weights are normalised; an unusable model makes the blend unusable
    same, _ = _blend_drifts([(first, 3.0), (second, 7.0)])
    assert same.shift_at(100.0) == pytest.approx(blended.shift_at(100.0))
    bad = make_drift("ccc", "2026-10-01T22:00:00", [0.5], [0.0], fibers_a, usable=False)
    assert not _blend_drifts([(first, 0.5), (bad, 0.5)])[0].usable
    assert _arc_weights([(first, 0.3), (second, 0.7)], ["aaa", "bbb"]) == pytest.approx([0.3, 0.7])
