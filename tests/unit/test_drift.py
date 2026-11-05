"""Tests for the drift step: the per-fiber shift and stretch fit on synthetic line lists,
the fit on a simulated calibration-fiber arc against the simulator's drift law, the
transfer of that drift to the science fibers against the science frame's truth, and the
degenerate cases of the across-detector model."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.loader import build_effective
from mastpipe.config.models import CalibParams, DriftParams, WavecalParams
from mastpipe.core.model import (
    LINE_COLUMNS,
    PER_FIBER_COLUMNS,
    DriftModel,
    Frame,
    FrameMeta,
    ReferenceSolution,
    TraceSet,
    WavelengthSolution,
)
from mastpipe.instruments.registry import load_header_map, load_instrument
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.pipeline.identify import _arc_lines
from mastpipe.sim.simulator import drift_x, simulate_night, undrift_x
from mastpipe.steps.atlas import catalogue_arrays
from mastpipe.steps.calib import build_master
from mastpipe.steps.drift import apply_drift, drift_unusable_reason, fit_drift, fit_fiber_drift
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.wavecal import build_reference, solve_full_arc
from mastpipe.steps.wavemodel import anchor_model, model_from_solution
from tests.conftest import truth_model, truth_traceset

N_PIX = 1000
ANCHOR_LOG_INTENSITY = 2.0
MAX_SHIFT_ERROR_PX = 0.05
"""Largest tolerated error of the fitted shift against the simulator's drift law."""
MAX_TRANSFER_ERROR_DISP = 0.3
"""Largest tolerated wavelength error of a transferred science fiber, in dispersions."""


# ------------------------------------------------------------------------------ synthetic


def make_solution(x_lines: np.ndarray, n_pix: int = N_PIX, fiber: int = 0) -> WavelengthSolution:
    lines = pd.DataFrame(
        {
            "x": np.asarray(x_lines, dtype=float),
            "wavelength": 500.0 + 0.03 * np.asarray(x_lines, dtype=float),
            "species": "Th",
            "dx_px": 0.0,
        }
    )[list(LINE_COLUMNS)]
    return WavelengthSolution(
        fiber=fiber,
        kind="grating",
        params=np.array([500.0, 3.0e4, 0.5]),
        n_pix=n_pix,
        lines=lines,
        rms_px=0.01,
        n_lines=len(lines),
        anchor_px=(float(lines["x"].min()), float(lines["x"].max())),
        anchor_nm=(float(lines["wavelength"].min()), float(lines["wavelength"].max())),
        trace_y_center=50.0,
        geometry={},
    )


def synthetic_lines(rng: np.random.Generator, n: int = 40) -> np.ndarray:
    return np.sort(rng.uniform(20.0, N_PIX - 20.0, n))


def test_fit_fiber_drift_recovers_shift_and_stretch():
    rng = np.random.default_rng(1)
    x_ref = synthetic_lines(rng)
    shift, stretch = 0.6, 2e-4
    x_obs = drift_x(x_ref, shift, stretch, N_PIX) + rng.normal(0.0, 0.005, x_ref.size)
    # a few reference lines missing from the observation, a few spurious detections
    x_obs = np.concatenate([x_obs[3:], [123.4, 456.7, 789.0]])
    result = fit_fiber_drift(x_obs, make_solution(x_ref), DriftParams())
    assert set(result) == {"shift_px", "stretch", "n_lines", "rms_px", "usable"}
    assert result["shift_px"] == pytest.approx(shift, abs=0.02)
    assert result["stretch"] == pytest.approx(stretch, abs=2e-5)
    assert result["n_lines"] >= x_ref.size - 3 - 2
    assert result["rms_px"] < 0.02
    assert result["usable"]


def test_fit_fiber_drift_with_large_shift_and_outliers():
    rng = np.random.default_rng(2)
    x_ref = synthetic_lines(rng)
    x_obs = x_ref - 37.3 + rng.normal(0.0, 0.01, x_ref.size)
    x_obs[5] += 0.9  # an outlier inside the matching tolerance: clipped, not fitted
    result = fit_fiber_drift(x_obs, make_solution(x_ref), DriftParams())
    assert result["shift_px"] == pytest.approx(-37.3, abs=0.02)
    assert result["stretch"] == 0.0
    assert result["n_lines"] == x_ref.size - 1
    assert result["usable"]


def test_insignificant_stretch_is_zero():
    rng = np.random.default_rng(3)
    x_ref = synthetic_lines(rng)
    x_obs = x_ref + 0.3 + rng.normal(0.0, 0.02, x_ref.size)
    result = fit_fiber_drift(x_obs, make_solution(x_ref), DriftParams())
    assert result["stretch"] == 0.0
    assert result["shift_px"] == pytest.approx(0.3, abs=0.02)
    assert result["usable"]


def test_fit_fiber_drift_too_few_lines_is_unusable():
    x_ref = np.array([100.0, 300.0, 500.0, 700.0])
    result = fit_fiber_drift(x_ref + 0.2, make_solution(x_ref), DriftParams(min_lines=5))
    assert result["n_lines"] == 4 and not result["usable"]
    empty = fit_fiber_drift(np.array([]), make_solution(x_ref), DriftParams())
    assert empty["n_lines"] == 0 and not empty["usable"] and empty["shift_px"] == 0.0


# ------------------------------------------------------------------- across the detector


def make_meta(mid_time: str = "2026-10-01T20:00:00") -> FrameMeta:
    return FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type="arc",
        configuration="Mg",
        mid_time=mid_time,
        arc_kind="calfiber",
        checksum="abc123",
    )


def make_reference(solutions: dict[int, WavelengthSolution]) -> ReferenceSolution:
    return build_reference("highspec", "Mg", "survey45", solutions, "2026-10-01T18:00:00", ())


def synthetic_tables(x_ref: np.ndarray, shifts: dict[int, float], rng: np.random.Generator) -> dict[int, pd.DataFrame]:
    tables = {}
    for fiber, shift in shifts.items():
        x_obs = x_ref + shift + rng.normal(0.0, 0.01, x_ref.size)
        tables[fiber] = pd.DataFrame({"x_center": x_obs, "use_x": True})
    return tables


def test_fit_drift_single_usable_fiber_is_a_constant():
    rng = np.random.default_rng(4)
    x_ref = synthetic_lines(rng)
    solutions = {0: make_solution(x_ref, fiber=0), 1: make_solution(x_ref, fiber=1)}
    tables = synthetic_tables(x_ref, {0: 0.4}, rng)
    tables[1] = pd.DataFrame({"x_center": x_ref[:3] + 0.1, "use_x": True})  # below min_lines
    traces = SimpleNamespace(
        traces=[SimpleNamespace(index=0, y_center=10.0), SimpleNamespace(index=1, y_center=90.0)],
        shape=(100, N_PIX),
    )
    model = fit_drift(tables, make_reference(solutions), traces, make_meta(), DriftParams())
    assert isinstance(model, DriftModel)
    assert list(model.per_fiber.columns) == list(PER_FIBER_COLUMNS)
    assert model.per_fiber["usable"].tolist() == [True, False]
    assert model.shift_coeffs.shape == (1,) and model.stretch_coeffs.shape == (1,)
    assert model.shift_at(10.0) == pytest.approx(0.4, abs=0.02)
    assert model.shift_at(90.0) == pytest.approx(0.4, abs=0.02)
    # one fiber gives a constant with no scatter to judge it by: recorded, not usable
    assert model.scatter_px == 0.0 and model.n_usable == 1 and not model.usable
    one = fit_drift(tables, make_reference(solutions), traces, make_meta(), DriftParams(min_fibers=1))
    assert one.usable
    assert model.arc_checksums == ("abc123",) and model.arc_mid_time == "2026-10-01T20:00:00"
    assert model.n_pix == N_PIX


def test_fit_drift_linear_across_fibers_and_unusable_without_fibers():
    rng = np.random.default_rng(5)
    x_ref = synthetic_lines(rng)
    y_of = {0: 10.0, 1: 50.0, 2: 90.0}
    shifts = {fiber: 0.2 + 0.005 * y for fiber, y in y_of.items()}
    solutions = {fiber: make_solution(x_ref, fiber=fiber) for fiber in y_of}
    traces = SimpleNamespace(traces=[SimpleNamespace(index=f, y_center=y) for f, y in y_of.items()], shape=(100, N_PIX))
    model = fit_drift(
        synthetic_tables(x_ref, shifts, rng),
        make_reference(solutions),
        traces,
        make_meta(),
        DriftParams(degree=1),
    )
    assert model.shift_coeffs.shape == (2,)
    assert model.shift_at(30.0) == pytest.approx(0.2 + 0.005 * 30.0, abs=0.02)
    assert model.scatter_px < 0.02 and model.usable

    none = fit_drift({}, make_reference(solutions), traces, make_meta(), DriftParams())
    assert not none.usable and len(none.per_fiber) == 0
    assert none.shift_at(30.0) == 0.0 and none.stretch_at(30.0) == 0.0


def test_drift_with_two_fibers_is_not_usable():
    """Two usable fibers fix a degree-1 shift model exactly, so the scatter gate cannot fail:
    the model needs ``min_fibers`` (three by default) and records how many it had."""
    rng = np.random.default_rng(6)
    x_ref = synthetic_lines(rng)
    y_of = {0: 10.0, 1: 90.0, 2: 50.0}
    solutions = {fiber: make_solution(x_ref, fiber=fiber) for fiber in y_of}
    traces = SimpleNamespace(traces=[SimpleNamespace(index=f, y_center=y) for f, y in y_of.items()], shape=(100, N_PIX))
    two = synthetic_tables(x_ref, {0: 0.3, 1: 0.9}, rng)
    model = fit_drift(two, make_reference(solutions), traces, make_meta(), DriftParams())
    assert model.per_fiber["usable"].all()
    assert model.scatter_px < 1e-9  # the line runs through both points
    assert model.n_usable == 2 and not model.usable
    assert drift_unusable_reason(model, DriftParams()) == "2 usable arc fibers, min_fibers is 3"

    # a configuration that accepts two fibers makes the same model usable
    relaxed = DriftParams(min_fibers=2)
    model = fit_drift(two, make_reference(solutions), traces, make_meta(), relaxed)
    assert model.n_usable == 2 and model.usable
    assert drift_unusable_reason(model, relaxed) is None

    three = synthetic_tables(x_ref, {0: 0.3, 1: 0.9, 2: 0.6}, rng)
    model = fit_drift(three, make_reference(solutions), traces, make_meta(), DriftParams())
    assert model.n_usable == 3 and model.usable

    # the scatter gate names itself
    bad = synthetic_tables(x_ref, {0: 0.3, 1: 0.9, 2: 2.0}, rng)
    model = fit_drift(bad, make_reference(solutions), traces, make_meta(), DriftParams())
    assert model.n_usable == 3 and not model.usable
    reason = drift_unusable_reason(model, DriftParams())
    assert reason is not None and reason.startswith("scatter") and "max_scatter_px 0.2" in reason


# ------------------------------------------------------------------------ simulated night


@pytest.fixture(scope="module")
def reference(sim_full_arc, sim_arc_lines, sim_arc_fluxes, mg_lines) -> ReferenceSolution:
    """The reference solution of the simulated full arc, solved from a truth anchor."""
    truth = sim_full_arc.truth
    geom = sim_full_arc.geom
    cat_wl, cat_sp = catalogue_arrays(mg_lines, usable_only=True)
    usable = np.isin(truth.lines_nm, cat_wl)
    bright = usable & truth.lines_on_chip & (truth.line_logI > ANCHOR_LOG_INTENSITY)
    first, last = np.flatnonzero(bright)[[0, -1]]
    pairs = [(float(truth.lines_nm[i]), round(float(truth.line_pixels[0][i]), 1)) for i in (first, last)]
    anchor = anchor_model(pairs, geom, sim_full_arc.n_pix)
    solutions = solve_full_arc(
        sim_arc_lines[0],
        sim_arc_fluxes,
        sim_full_arc.traces,
        anchor,
        cat_wl,
        cat_sp,
        WavecalParams(),
    )
    full = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    return build_reference("highspec", "Mg", "survey45", solutions, str(full["mid_time"]), ("fullarc",))


@pytest.fixture(scope="module")
def cal_arc(sim_night, sim_full_arc, highspec_effective, survey_layout) -> SimpleNamespace:
    """The last simulated calibration-fiber arc taken to its line tables (arc fibers only)."""
    night_dir, truth = sim_night
    detector = highspec_effective.detector
    refs = scan_folder(night_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        if ref.meta.exposure_type in ("bias", "background"):
            by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    masters = {kind: build_master(frames, kind, CalibParams()) for kind, frames in by_type.items()}
    frames = truth.frames
    row = frames.loc[frames["arc_kind"] == "calfiber"].iloc[-1]
    lines = _arc_lines(
        night_dir / str(row["file"]),
        highspec_effective,
        survey_layout,
        night_dir,
        masters=masters,
        traceset=sim_full_arc.traces,
    )
    arc_fibers = [trace.index for trace in sim_full_arc.traces.by_role("arc")]
    tables = {fiber: lines.tables[fiber] for fiber in arc_fibers}
    return SimpleNamespace(meta=lines.frame.meta, tables=tables, fibers=arc_fibers, row=row, truth=truth)


@pytest.fixture(scope="module")
def drift(cal_arc, reference, sim_full_arc) -> DriftModel:
    return fit_drift(cal_arc.tables, reference, sim_full_arc.traces, cal_arc.meta, DriftParams())


def test_fit_drift_on_simulated_cal_arc(cal_arc, drift, sim_full_arc):
    truth_shift = cal_arc.truth.trace_shifts(str(cal_arc.row["file"]))  # per trace
    truth_stretch = float(cal_arc.row["stretch"])
    assert float(cal_arc.row["shift_px"]) > 0.1  # the law has moved the lines since the full arc
    per_fiber = drift.per_fiber.set_index("fiber")
    assert sorted(per_fiber.index) == sorted(cal_arc.fibers)
    assert per_fiber["usable"].all()
    assert (per_fiber["n_lines"] >= 15).all()
    expected = truth_shift[[cal_arc.truth.row_of(fiber) for fiber in per_fiber.index]]
    assert np.allclose(per_fiber["shift_px"], expected, atol=MAX_SHIFT_ERROR_PX)
    assert np.allclose(per_fiber["stretch"], truth_stretch, atol=5e-5)
    for trace in sim_full_arc.traces.by_role("arc"):
        assert per_fiber.loc[trace.index, "y_center"] == pytest.approx(trace.y_center)
        assert drift.shift_at(trace.y_center) == pytest.approx(
            truth_shift[cal_arc.truth.row_of(trace.index)], abs=MAX_SHIFT_ERROR_PX
        )
    assert drift.scatter_px <= MAX_SHIFT_ERROR_PX
    assert drift.usable
    assert drift.arc_checksums == (cal_arc.meta.checksum,)
    assert drift.arc_mid_time == cal_arc.meta.mid_time == str(cal_arc.row["mid_time"])
    assert drift.n_pix == sim_full_arc.n_pix


GRADIENT_PX = 2.0
"""Bottom-to-top shift difference of the gradient night at its last frame, px."""


@pytest.fixture(scope="module")
def gradient_night(tmp_path_factory, survey_layout, mg_lines) -> SimpleNamespace:
    """A night with a 2 px drift gradient across the detector (default jitter): its full arc
    taken to a reference from truth anchors, its calibration-fiber arc to line tables, both
    through the pipeline's own ``_arc_lines`` on the truth traces."""
    effective = build_effective(load_instrument("highspec"), "Mg", "main", layout="survey45")
    night_dir = tmp_path_factory.mktemp("gradient_night")
    truth = simulate_night(
        effective,
        survey_layout,
        mg_lines,
        night_dir,
        seed=5,
        n_bias=3,
        n_dark=1,
        n_flat=1,
        n_cal_arc=1,
        n_science=1,
        drift_gradient_px=GRADIENT_PX,
    )
    detector = effective.detector
    refs = scan_folder(night_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        if ref.meta.exposure_type in ("bias", "background"):
            by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    masters = {kind: build_master(frames, kind, CalibParams()) for kind, frames in by_type.items()}
    frames = truth.frames
    full = frames.loc[frames["arc_kind"] == "full"].iloc[0]
    cal = frames.loc[frames["arc_kind"] == "calfiber"].iloc[0]
    traces = truth_traceset(truth, tuple(truth.extras["shape"]))

    def lines_of(row):
        path = night_dir / str(row["file"])
        return _arc_lines(path, effective, survey_layout, night_dir, masters, traces)

    full_lines = lines_of(full)
    geom = grating_geometry(effective.geometry, detector.pixel_size_um)
    cat_wl, cat_sp = catalogue_arrays(mg_lines, usable_only=True)
    bright = np.isin(truth.lines_nm, cat_wl) & truth.lines_on_chip
    bright &= truth.line_logI > ANCHOR_LOG_INTENSITY
    first, last = np.flatnonzero(bright)[[0, -1]]
    row0 = truth.row_of(0)
    pairs = [(float(truth.lines_nm[i]), float(truth.line_pixels[row0][i])) for i in (first, last)]
    solutions = solve_full_arc(
        full_lines.tables,
        full_lines.fluxes,
        traces,
        anchor_model(pairs, geom, int(truth.extras["n_pix"])),
        cat_wl,
        cat_sp,
        effective.params.wavecal,
    )
    reference = build_reference("highspec", "Mg", "survey45", solutions, str(full["mid_time"]), ("fullarc",))
    cal_lines = lines_of(cal)
    arc_fibers = [trace.index for trace in traces.by_role("arc")]
    return SimpleNamespace(
        truth=truth,
        traces=traces,
        reference=reference,
        tables={fiber: cal_lines.tables[fiber] for fiber in arc_fibers},
        meta=cal_lines.frame.meta,
        row=cal,
        params=effective.params.drift,
    )


def test_fit_drift_follows_a_y_gradient(gradient_night):
    """The pipeline's drift model follows a shift that varies across the detector: fitted
    against the night's own reference, ``shift_at(y_i)`` matches the simulator's shift at
    every arc fiber's row and its slope in ``y`` has the planted sign and size. A model
    flat in ``y`` (``degree = 0``) misses the outer fibers, which is what makes this test
    sensitive to a drift step that ignores ``y`` or applies the drift at the wrong row."""
    night = gradient_night
    truth = night.truth
    gradient = float(night.row["gradient_px"])
    assert gradient > 0.5  # the cal arc is late enough for the gradient to have grown
    shifts = truth.trace_shifts(str(night.row["file"]))
    arc = night.traces.by_role("arc")
    expected = np.array([shifts[truth.row_of(t.index)] for t in arc])
    ys = np.array([t.y_center for t in arc])
    assert np.ptp(expected) > 10 * MAX_SHIFT_ERROR_PX

    drift = fit_drift(night.tables, night.reference, night.traces, night.meta, night.params)
    assert drift.usable and night.params.degree == 1
    fitted = np.array([drift.shift_at(y) for y in ys])
    np.testing.assert_allclose(fitted, expected, atol=MAX_SHIFT_ERROR_PX)
    ny = truth.extras["shape"][0]
    slope = (fitted[-1] - fitted[0]) / (ys[-1] - ys[0])
    assert slope == pytest.approx(gradient / ny, rel=0.1)  # positive: top drifts further

    # the same fit with the configuration forcing a drift model flat in y
    overrides = [{"drift": {"degree": 0}}]
    flat_params = build_effective(
        load_instrument("highspec"), "Mg", "main", overrides=overrides, layout="survey45"
    ).params.drift
    assert flat_params.degree == 0
    flat = fit_drift(night.tables, night.reference, night.traces, night.meta, flat_params)
    flat_fit = np.array([flat.shift_at(y) for y in ys])
    assert np.ptp(flat_fit) == 0.0
    assert np.max(np.abs(flat_fit - expected)) > 5 * MAX_SHIFT_ERROR_PX


def test_apply_drift_matches_science_truth(cal_arc, drift, reference, sim_full_arc):
    truth = cal_arc.truth
    frames = truth.frames
    # the science frame nearest the calibration arc in time
    science = frames.loc[frames["exposure_type"] == "science"]
    t_arc = pd.Timestamp(cal_arc.row["mid_time"])
    nearest = science.loc[(pd.to_datetime(science["mid_time"]) - t_arc).abs().idxmin()]
    shifts, stretch = truth.trace_shifts(str(nearest["file"])), float(nearest["stretch"])
    geom = sim_full_arc.geom
    n_pix = sim_full_arc.n_pix

    wavelengths, skipped = apply_drift(reference, [(drift, 1.0)], sim_full_arc.traces, geom, tolerance=3.0, n_pix=n_pix)
    assert set(wavelengths) == set(reference.solutions) and skipped == {}
    x = np.arange(n_pix, dtype=float)
    dips = np.array(truth.extras["absorption_nm"], dtype=float)
    for trace in sim_full_arc.traces.by_role("science"):
        model = truth_model(truth, trace.index, geom, n_pix)
        shift = shifts[truth.row_of(trace.index)]
        lam_true = model.wavelength(undrift_x(x, shift, stretch, n_pix))
        dispersion = model.dispersion(x)
        x_dips = drift_x(model.pixel(dips), shift, stretch, n_pix)
        lam = wavelengths[trace.index]
        assert lam.shape == (n_pix,)
        assert np.max(np.abs(lam - lam_true) / dispersion) < MAX_TRANSFER_ERROR_DISP
        # the drifted absorption dips land at their wavelengths
        lam_at_dips = np.interp(x_dips, x, lam)
        assert np.all(np.abs(lam_at_dips - dips) / np.interp(x_dips, x, dispersion) < 0.3)

    # no drift at all is the reference map itself
    plain, _ = apply_drift(reference, [], sim_full_arc.traces, geom)
    for fiber, solution in reference.solutions.items():
        assert np.allclose(plain[fiber], model_from_solution(solution, geom).wavelength(x))
    # the transfer undoes the drift: the reference map evaluated at the undrifted column
    fiber = sim_full_arc.traces.by_role("science")[0].index
    s = float(drift.shift_at(sim_full_arc.traces.traces[fiber].y_center))
    e = float(drift.stretch_at(sim_full_arc.traces.traces[fiber].y_center))
    expected = model_from_solution(reference.solutions[fiber], geom).wavelength(undrift_x(x, s, e, n_pix))
    assert np.allclose(wavelengths[fiber], expected)
    # weights combine linearly
    halves, _ = apply_drift(reference, [(drift, 0.5), (drift, 0.5)], sim_full_arc.traces, geom)
    assert np.allclose(halves[fiber], wavelengths[fiber])


def test_apply_drift_skips_a_fiber_whose_trace_moved_and_checks_the_grid(reference, sim_full_arc):
    """A trace centre far from the reference's ``trace_y_center`` (a fiber keyed to the
    wrong solution) is not transferred but reported; a reference on another pixel grid
    is refused outright."""
    geom = sim_full_arc.geom
    traces = sim_full_arc.traces
    moved = next(iter(reference.solutions))
    shifted = TraceSet(
        traces=[replace(t, y_center=t.y_center + 5.0) if t.index == moved else t for t in traces],
        shape=traces.shape,
    )
    wavelengths, skipped = apply_drift(reference, [], shifted, geom, tolerance=3.0)
    assert moved not in wavelengths and set(wavelengths) == set(reference.solutions) - {moved}
    assert list(skipped) == [moved]
    assert "5.0" in skipped[moved] and "3.0" in skipped[moved]
    # within the tolerance, or without one, the fiber is transferred
    assert moved in apply_drift(reference, [], shifted, geom, tolerance=6.0)[0]
    assert moved in apply_drift(reference, [], shifted, geom)[0]
    # a fiber the reference solves but the night did not trace is still mapped (unused)
    fewer = TraceSet(traces=[t for t in traces if t.index != moved], shape=traces.shape)
    assert moved in apply_drift(reference, [], fewer, geom, tolerance=3.0)[0]
    with pytest.raises(ValueError, match=f"{reference.identifier}.*{sim_full_arc.n_pix}.*1000"):
        apply_drift(reference, [], traces, geom, n_pix=1000)
