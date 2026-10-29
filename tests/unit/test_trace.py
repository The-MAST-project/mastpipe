"""Tests for the trace step: trace finding, comb matching, layout assignment, arc
classification and relative throughput, against simulator truth and hand-built combs."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

from mastpipe.config.models import CalibParams, TraceParams
from mastpipe.core.model import Frame, FrameMeta, Trace, TraceSet
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.sim.simulator import FWHM_PX, GAUSSIAN_FWHM, render_traces
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.trace import (
    classify_arc,
    comb_shift,
    find_traces,
    match_layout,
    relative_throughput,
    robust_sigma,
    trace_source,
)

SURVEY_PARAMS = TraceParams(trace_halfwidth=4)
"""The 45-trace layout has an 8.5 px pitch: the HighSpec default (``trace_halfwidth = 8``,
``min_separation = 16``) cannot separate neighbours, so the survey tests use a 4 px halfwidth."""
APERTURE_HALFWIDTH = 3.0
"""Rows summed around a trace on the survey layout (a 3 px FWHM image, 8.5 px pitch)."""


# ------------------------------------------------------------------------------- fixtures


def calibrated_night(sim, highspec_effective) -> SimpleNamespace:
    """Calibrated simulated frames: master flat, full arc, calibration-fiber arc, science."""
    out_dir, truth = sim
    detector = highspec_effective.detector
    refs = scan_folder(out_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    params = CalibParams()
    bias = build_master(by_type["bias"], "bias", params)
    dark = build_master(by_type["dark"], "dark", params)
    background = build_master(by_type["background"], "background", params)
    calib = build_calibration_set(bias, dark, background, detector)

    def clean(frame: Frame, subtract_background: bool = True) -> Frame:
        return calibrate_frame(
            frame, calib, subtract_background=subtract_background, saturation=detector.saturation
        )

    flats = [clean(frame, subtract_background=False) for frame in by_type["flat"]]
    master = build_master(flats, "flat", params)
    flat = Frame(data=master.data, meta=flats[0].meta, variance=master.variance)
    kind_of = dict(zip(truth.frames["file"], truth.frames["arc_kind"], strict=True))
    arcs: dict[str, Frame] = {}
    for frame in by_type["arc"]:
        arcs.setdefault(kind_of[frame.meta.source], frame)
    return SimpleNamespace(
        truth=truth,
        flat=flat,
        full_arc=clean(arcs["full"]),
        cal_arc=clean(arcs["calfiber"]),
        science=clean(by_type["science"][0]),
    )


@pytest.fixture(scope="module")
def night(sim_night, highspec_effective) -> SimpleNamespace:
    return calibrated_night(sim_night, highspec_effective)


@pytest.fixture(scope="module")
def night_lab3(sim_night_lab3, highspec_effective) -> SimpleNamespace:
    return calibrated_night(sim_night_lab3, highspec_effective)


@pytest.fixture(scope="module")
def survey_traces(night, survey_layout) -> TraceSet:
    return match_layout(find_traces(night.flat.data, SURVEY_PARAMS), survey_layout, 3.0)


def make_trace(index: int, y: float, x_center: float, slope: float = 0.0) -> Trace:
    coeffs = np.array([0.0, slope, y - slope * x_center])
    return Trace(
        index=index,
        coeffs=coeffs,
        y_center=y,
        slope=slope,
        rms=0.0,
        n_points=10,
        snr=100.0,
        x_samples=np.arange(10.0),
        y_samples=np.full(10, y),
        kept=np.ones(10, dtype=bool),
    )


def toy_frame(data: np.ndarray, exposure_type: str = "arc") -> Frame:
    meta = FrameMeta(instrument="toy", camera="main", exposure_type=exposure_type)
    return Frame(data=data.astype(np.float32), meta=meta)


# ----------------------------------------------------------------------------- robust_sigma


def test_robust_sigma_matches_gaussian_width():
    rng = np.random.default_rng(1)
    values = rng.normal(3.0, 2.0, size=20000)
    assert abs(robust_sigma(values) - 2.0) < 0.1
    values[:200] = 1e6  # outliers barely move it
    assert abs(robust_sigma(values) - 2.0) < 0.1
    assert robust_sigma(np.full(10, 7.0)) == 0.0
    assert robust_sigma(np.array([1.0, np.nan, 3.0])) == pytest.approx(1.4826)


# ------------------------------------------------------------------------------ find_traces


def test_find_traces_on_simulated_flat_recovers_45(night):
    truth = night.truth
    traces = find_traces(night.flat.data, SURVEY_PARAMS)
    assert traces.shape == night.flat.data.shape
    assert len(traces) == len(truth.trace_coeffs) == 45
    assert [trace.index for trace in traces] == list(range(45))
    centers = np.array([trace.y_center for trace in traces])
    assert np.all(np.diff(centers) > 0)
    nx = night.flat.data.shape[1]
    worst = 0.0
    for x in (100, nx // 2, nx - 100):
        found = np.array([trace.y(x) for trace in traces])
        expected = np.polyval(truth.trace_coeffs.T, x)
        worst = max(worst, float(np.max(np.abs(found - expected))))
    assert worst < 0.3
    for trace in traces:
        assert trace.n_points >= SURVEY_PARAMS.min_points
        assert trace.rms < 0.1
        assert trace.snr > 20
        assert trace.kept.sum() == trace.n_points
        assert len(trace.x_samples) == len(trace.y_samples) == len(trace.kept)
        assert trace.role is None and trace.layout_index is None
    # slope at the centre column follows the simulator's 0.002
    slopes = np.array([trace.slope for trace in traces])
    assert np.allclose(slopes, 0.002, atol=2e-4)


def test_find_traces_lab3_like_curvature():
    shape = (398, 1598)
    nx = shape[1]
    xc = nx / 2.0
    rng = np.random.default_rng(3)
    y0 = np.array([17.4, 189.0, 360.7])
    slope = np.full(3, 0.002)
    curvature = np.array([2e-6, -1.5e-6, 2e-6])
    coeffs = np.column_stack(
        [curvature, slope - 2 * curvature * xc, y0 - slope * xc + curvature * xc**2]
    )
    column_factor = FWHM_PX / GAUSSIAN_FWHM * np.sqrt(2 * np.pi)
    fluxes = np.full((3, nx), 15000.0 * column_factor)
    image = render_traces(shape, coeffs, fluxes, FWHM_PX)
    image += rng.normal(0.0, 20.0, size=shape)

    traces = find_traces(image, TraceParams())  # the HighSpec defaults: halfwidth 8, degree 2
    assert len(traces) == 3
    x = np.array([0, 400, 799, 1200, 1597], dtype=float)
    for trace, coeff in zip(traces, coeffs, strict=True):
        assert trace.rms < 0.02
        assert np.max(np.abs(trace.y(x) - np.polyval(coeff, x))) < 0.05
        assert trace.slope == pytest.approx(np.polyval(np.polyder(coeff), nx // 2), abs=5e-5)
    # the sag over the detector is ~1.3 px: a straight line cannot follow it
    linear = find_traces(image, TraceParams(degree=1))
    assert len(linear) == 3
    assert all(trace.rms > 0.1 for trace in linear)


def test_find_traces_on_simulated_lab3_flat_gives_exactly_3(night_lab3, highspec_effective):
    truth = night_lab3.truth
    params = highspec_effective.params.trace  # the HighSpec configuration's values
    assert params.trace_halfwidth == 8 and params.collapse == "median"
    traces = find_traces(night_lab3.flat.data, params)
    assert len(traces) == len(truth.trace_coeffs) == 3  # no noise peaks
    nx = night_lab3.flat.data.shape[1]
    for x in (100, nx // 2, nx - 100):
        found = np.array([trace.y(x) for trace in traces])
        assert np.max(np.abs(found - np.polyval(truth.trace_coeffs.T, x))) < 0.3
    assert all(trace.rms < 0.05 for trace in traces)


def test_find_traces_on_simulated_full_arc_with_percentile_collapse(night):
    truth = night.truth
    arc = night.full_arc.data
    nx = arc.shape[1]
    expected = np.polyval(truth.trace_coeffs.T, nx // 2)
    # the median along x of a line-only spectrum is the level between lines: no traces
    assert len(find_traces(arc, SURVEY_PARAMS)) < 45
    params = TraceParams(trace_halfwidth=4, collapse="percentile")
    assert params.collapse_percentile == 90.0
    traces = find_traces(arc, params)
    assert len(traces) == 45
    found = np.array([trace.y_center for trace in traces])
    assert np.max(np.abs(found - expected)) < 0.3
    assert all(trace.n_points >= params.min_points for trace in traces)


def test_trace_params_collapse_is_validated():
    assert TraceParams().collapse == "median"
    assert TraceParams(collapse="mean").collapse == "mean"
    with pytest.raises(ValidationError):
        TraceParams(collapse="maximum")
    with pytest.raises(ValueError, match="collapse_percentile"):
        find_traces(np.zeros((20, 50)), TraceParams(collapse="percentile", collapse_percentile=100))


def test_find_traces_empty_image_gives_no_traces():
    rng = np.random.default_rng(5)
    image = rng.normal(0.0, 5.0, size=(200, 600))
    traces = find_traces(image, TraceParams())
    assert len(traces) == 0
    assert traces.shape == (200, 600)


# ------------------------------------------------------------------------------- comb_shift


def test_comb_shift_subpixel():
    a = np.array([10.0, 25.5, 40.0, 70.2, 95.0, 130.3])
    for true_shift in (2.37, -1.84, 0.0):
        b = a + true_shift
        assert abs(comb_shift(a, b, max_shift=5.0, step=0.1) - true_shift) < 0.02
    # an unmatched member on each side does not pull the estimate
    b = np.concatenate([a[:-1] + 2.37, [57.0]])
    assert abs(comb_shift(a, b, max_shift=5.0) - 2.37) < 0.02
    # beyond the search window the best shift inside it is returned, not the true one
    assert abs(comb_shift(a, a + 2.37, max_shift=1.0)) <= 1.0
    assert comb_shift(np.array([]), a, max_shift=5.0) == 0.0


# ----------------------------------------------------------------------------- match_layout


def test_match_layout_with_global_shift_and_missing(night, survey_layout):
    truth = night.truth
    image = np.roll(night.flat.data, 2, axis=0)  # every trace 2 rows further down
    dropped = 17
    nx = image.shape[1]
    x = np.arange(nx)
    rows = np.rint(np.polyval(truth.trace_coeffs[dropped], x)).astype(int) + 2
    for offset in range(-4, 5):
        image[rows + offset, x] = 0.0
    traces = find_traces(image, SURVEY_PARAMS)
    assert len(traces) == 44

    matched = match_layout(traces, survey_layout, SURVEY_PARAMS.layout_tolerance)
    assert matched.layout_shift == pytest.approx(2.0, abs=0.15)
    assert matched.missing == (dropped,)
    assert matched.extra == ()
    assert matched.shape == traces.shape
    assert len(matched) == 44
    for trace in matched:
        row = survey_layout.iloc[trace.layout_index]
        assert trace.role == row["role"]
        assert trace.telescope == row["telescope"]
        assert trace.aperture is None  # the survey layout leaves the aperture to the parameters
        assert abs(trace.y_center - (row["y_nominal"] + matched.layout_shift)) < 1.0
    assert [trace.layout_index for trace in matched] == [
        i for i in range(len(survey_layout)) if i != dropped
    ]
    # the input is left untouched
    assert all(trace.role is None for trace in traces)


def test_match_layout_reports_extra_and_copies_aperture(lab3_layout):
    x_center = 799.0
    nominal = lab3_layout["y_nominal"].to_numpy()
    ys = np.concatenate([nominal - 1.3, [250.0]])  # one spurious detection
    traces = [make_trace(i, y, x_center) for i, y in enumerate(np.sort(ys))]
    traceset = TraceSet(traces=traces, shape=(398, 1598))
    matched = match_layout(traceset, lab3_layout, 3.0)
    assert matched.layout_shift == pytest.approx(-1.3, abs=0.1)
    assert matched.missing == ()
    assert matched.extra == (2,)
    roles = [trace.role for trace in matched]
    assert roles == ["arc", "arc", None, "arc"]
    assert [trace.layout_index for trace in matched] == [0, 1, None, 2]
    for trace in matched:
        if trace.layout_index is not None:
            expected = lab3_layout.iloc[trace.layout_index]["aperture"]
            assert trace.aperture is None if np.isnan(expected) else trace.aperture == expected


def test_match_layout_one_to_one_by_distance(lab3_layout):
    x_center = 799.0
    # two detections near layout row 1: the closer one wins, the other is extra
    nominal = lab3_layout["y_nominal"].to_numpy()
    ys = np.sort([nominal[0], nominal[1] - 1.0, nominal[1] + 0.4, nominal[2]])
    traces = [make_trace(i, y, x_center) for i, y in enumerate(ys)]
    matched = match_layout(TraceSet(traces=traces, shape=(398, 1598)), lab3_layout, 2.0)
    assert matched.layout_shift == pytest.approx(0.0, abs=0.1)
    assert matched.missing == ()
    assert matched.extra == (1,)
    assert [trace.layout_index for trace in matched] == [0, None, 1, 2]


# ----------------------------------------------------------------------------- classify_arc


def test_classify_arc_full_vs_calfiber_vs_unknown(night, survey_traces):
    assert survey_traces.missing == () and survey_traces.extra == ()
    arc_traces = {trace.index for trace in survey_traces.by_role("arc")}
    assert len(arc_traces) == 5

    kind, flags = classify_arc(night.full_arc, survey_traces, APERTURE_HALFWIDTH)
    assert kind == "full"
    assert set(flags) == set(range(45)) and all(flags.values())

    kind, flags = classify_arc(night.cal_arc, survey_traces, APERTURE_HALFWIDTH)
    assert kind == "calfiber"
    assert {i for i, lit in flags.items() if lit} == arc_traces

    kind, flags = classify_arc(night.science, survey_traces, APERTURE_HALFWIDTH)
    assert kind == "unknown"
    assert not any(flags[i] for i in arc_traces)
    assert all(flags[trace.index] for trace in survey_traces.by_role("science"))


def test_classify_arc_toy_frames():
    shape = (60, 200)
    x_center = 100.0
    traces = [make_trace(i, y, x_center) for i, y in enumerate((10.0, 30.0, 50.0))]
    roles = ("arc", "science", "sky")
    traces = [
        replace(trace, role=role, layout_index=i)
        for i, (trace, role) in enumerate(zip(traces, roles, strict=True))
    ]
    traceset = TraceSet(traces=traces, shape=shape)
    rng = np.random.default_rng(7)
    noise = rng.normal(0.0, 2.0, size=shape)

    lit_all = noise.copy()
    for trace in traces:
        lit_all[int(trace.y_center), :] += 100.0
    assert classify_arc(toy_frame(lit_all), traceset, 2.0)[0] == "full"

    arc_only = noise.copy()
    for start in range(5, 200, 10):  # 3-column emission lines only, nothing between them
        arc_only[10, start : start + 3] += 500.0
    kind, flags = classify_arc(toy_frame(arc_only), traceset, 2.0)
    assert kind == "calfiber"
    assert flags == {0: True, 1: False, 2: False}

    arc_hit = arc_only.copy()
    arc_hit[30, 77] += 5000.0  # one cosmic ray on the dark science trace does not light it
    kind, flags = classify_arc(toy_frame(arc_hit), traceset, 2.0)
    assert kind == "calfiber"
    assert flags == {0: True, 1: False, 2: False}

    kind, flags = classify_arc(toy_frame(noise), traceset, 2.0)
    assert kind == "unknown"
    assert not any(flags.values())

    # an unmatched detection (no role) is reported but does not decide the kind
    with_extra = traceset.traces + [make_trace(3, 40.0, x_center)]
    kind, flags = classify_arc(toy_frame(lit_all), TraceSet(with_extra, shape), 2.0)
    assert kind == "full"
    assert flags[3] is False


# ---------------------------------------------------------------------- relative_throughput


def test_relative_throughput_matches_truth(night, survey_traces):
    truth = night.truth
    relative = relative_throughput(night.flat, survey_traces, APERTURE_HALFWIDTH)
    assert sorted(relative) == list(range(45))
    expected = truth.throughput / np.median(truth.throughput)
    for trace in survey_traces:
        assert relative[trace.index] == pytest.approx(expected[trace.layout_index], rel=0.03)
    assert np.median(list(relative.values())) == pytest.approx(1.0)
    assert relative_throughput(None, survey_traces, APERTURE_HALFWIDTH) == {}


# ---------------------------------------------------------------------------- trace_source


def test_trace_source_fallback_order():
    flat = toy_frame(np.zeros((4, 4)), "flat")
    arc_a = toy_frame(np.zeros((4, 4)), "arc")
    arc_b = toy_frame(np.zeros((4, 4)), "arc")
    assert trace_source({"flat": [flat], "arc_full": [arc_a]}) == ("flat", flat)
    assert trace_source({"flat": [], "arc_full": [arc_a, arc_b]}) == ("arc_full", arc_a)
    assert trace_source({"arc_calfiber": [arc_a], "science": [flat]}) == ("none", None)
    assert trace_source({}) == ("none", None)
