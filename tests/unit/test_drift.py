"""Tests for the drift step: the per-fiber shift and stretch fit on synthetic line lists,
the fit on a simulated calibration-fiber arc against the simulator's drift law, the
transfer of that drift to the science fibers against the science frame's truth, and the
degenerate cases of the across-detector model."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.models import CalibParams, DriftParams, WavecalParams
from mastpipe.core.model import (
    LINE_COLUMNS,
    PER_FIBER_COLUMNS,
    DriftModel,
    Frame,
    FrameMeta,
    ReferenceSolution,
    WavelengthSolution,
)
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.pipeline.identify import _arc_lines
from mastpipe.sim.simulator import drift_x, undrift_x
from mastpipe.steps.atlas import catalogue_arrays
from mastpipe.steps.calib import build_master
from mastpipe.steps.drift import apply_drift, fit_drift, fit_fiber_drift
from mastpipe.steps.wavecal import build_reference, solve_full_arc
from mastpipe.steps.wavemodel import GratingModel, anchor_model, model_from_solution

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
        configuration="mgb518",
        mid_time=mid_time,
        arc_kind="calfiber",
        checksum="abc123",
    )


def make_reference(solutions: dict[int, WavelengthSolution]) -> ReferenceSolution:
    return build_reference("highspec", "mgb518", "survey45", solutions, "2026-10-01T18:00:00", ())


def synthetic_tables(
    x_ref: np.ndarray, shifts: dict[int, float], rng: np.random.Generator
) -> dict[int, pd.DataFrame]:
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
    assert model.scatter_px == 0.0 and model.usable
    assert model.arc_checksums == ("abc123",) and model.arc_mid_time == "2026-10-01T20:00:00"
    assert model.n_pix == N_PIX


def test_fit_drift_linear_across_fibers_and_unusable_without_fibers():
    rng = np.random.default_rng(5)
    x_ref = synthetic_lines(rng)
    y_of = {0: 10.0, 1: 50.0, 2: 90.0}
    shifts = {fiber: 0.2 + 0.005 * y for fiber, y in y_of.items()}
    solutions = {fiber: make_solution(x_ref, fiber=fiber) for fiber in y_of}
    traces = SimpleNamespace(
        traces=[SimpleNamespace(index=f, y_center=y) for f, y in y_of.items()], shape=(100, N_PIX)
    )
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


# ------------------------------------------------------------------------ simulated night


@pytest.fixture(scope="module")
def reference(sim_full_arc, sim_arc_lines, sim_arc_fluxes, mgb518_lines) -> ReferenceSolution:
    """The reference solution of the simulated full arc, solved from a truth anchor."""
    truth = sim_full_arc.truth
    geom = sim_full_arc.geom
    cat_wl, cat_sp = catalogue_arrays(mgb518_lines, usable_only=True)
    usable = np.isin(truth.lines_nm, cat_wl)
    bright = usable & truth.lines_on_chip & (truth.line_logI > ANCHOR_LOG_INTENSITY)
    first, last = np.flatnonzero(bright)[[0, -1]]
    pairs = [
        (float(truth.lines_nm[i]), round(float(truth.line_pixels[0][i]), 1)) for i in (first, last)
    ]
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
    return build_reference(
        "highspec", "mgb518", "survey45", solutions, str(full["mid_time"]), ("fullarc",)
    )


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
    return SimpleNamespace(
        meta=lines.frame.meta, tables=tables, fibers=arc_fibers, row=row, truth=truth
    )


@pytest.fixture(scope="module")
def drift(cal_arc, reference, sim_full_arc) -> DriftModel:
    return fit_drift(cal_arc.tables, reference, sim_full_arc.traces, cal_arc.meta, DriftParams())


def test_fit_drift_on_simulated_cal_arc(cal_arc, drift, sim_full_arc):
    truth_shift = float(cal_arc.row["shift_px"])
    truth_stretch = float(cal_arc.row["stretch"])
    assert truth_shift > 0.1  # the law has moved the lines since the full arc
    per_fiber = drift.per_fiber.set_index("fiber")
    assert sorted(per_fiber.index) == sorted(cal_arc.fibers)
    assert per_fiber["usable"].all()
    assert (per_fiber["n_lines"] >= 15).all()
    assert np.allclose(per_fiber["shift_px"], truth_shift, atol=MAX_SHIFT_ERROR_PX)
    assert np.allclose(per_fiber["stretch"], truth_stretch, atol=5e-5)
    for trace in sim_full_arc.traces.by_role("arc"):
        assert per_fiber.loc[trace.index, "y_center"] == pytest.approx(trace.y_center)
        assert drift.shift_at(trace.y_center) == pytest.approx(truth_shift, abs=MAX_SHIFT_ERROR_PX)
    assert drift.scatter_px <= MAX_SHIFT_ERROR_PX
    assert drift.usable
    assert drift.arc_checksums == (cal_arc.meta.checksum,)
    assert drift.arc_mid_time == cal_arc.meta.mid_time == str(cal_arc.row["mid_time"])
    assert drift.n_pix == sim_full_arc.n_pix


def test_apply_drift_matches_science_truth(cal_arc, drift, reference, sim_full_arc):
    truth = cal_arc.truth
    frames = truth.frames
    # the science frame nearest the calibration arc in time
    science = frames.loc[frames["exposure_type"] == "science"]
    t_arc = pd.Timestamp(cal_arc.row["mid_time"])
    nearest = science.loc[(pd.to_datetime(science["mid_time"]) - t_arc).abs().idxmin()]
    shift, stretch = float(nearest["shift_px"]), float(nearest["stretch"])
    geom = sim_full_arc.geom
    n_pix = sim_full_arc.n_pix
    truth_model = GratingModel(truth.x_c, truth.B, geom.alpha_rad, None, n_pix, geom)

    wavelengths = apply_drift(reference, [(drift, 1.0)], sim_full_arc.traces, geom)
    assert set(wavelengths) == set(reference.solutions)
    x = np.arange(n_pix, dtype=float)
    lam_true = truth_model.wavelength(undrift_x(x, shift, stretch, n_pix))
    dispersion = truth_model.dispersion(x)
    dips = np.array(truth.extras["absorption_nm"], dtype=float)
    x_dips = drift_x(truth_model.pixel(dips), shift, stretch, n_pix)
    for trace in sim_full_arc.traces.by_role("science"):
        lam = wavelengths[trace.index]
        assert lam.shape == (n_pix,)
        assert np.max(np.abs(lam - lam_true) / dispersion) < MAX_TRANSFER_ERROR_DISP
        # the drifted absorption dips land at their wavelengths
        lam_at_dips = np.interp(x_dips, x, lam)
        assert np.all(np.abs(lam_at_dips - dips) / np.interp(x_dips, x, dispersion) < 0.3)

    # no drift at all is the reference map itself
    plain = apply_drift(reference, [], sim_full_arc.traces, geom)
    for fiber, solution in reference.solutions.items():
        assert np.allclose(plain[fiber], model_from_solution(solution, geom).wavelength(x))
    # the transfer undoes the drift: the reference map evaluated at the undrifted column
    fiber = sim_full_arc.traces.by_role("science")[0].index
    s = float(drift.shift_at(sim_full_arc.traces.traces[fiber].y_center))
    e = float(drift.stretch_at(sim_full_arc.traces.traces[fiber].y_center))
    expected = model_from_solution(reference.solutions[fiber], geom).wavelength(
        undrift_x(x, s, e, n_pix)
    )
    assert np.allclose(wavelengths[fiber], expected)
    # weights combine linearly
    halves = apply_drift(reference, [(drift, 0.5), (drift, 0.5)], sim_full_arc.traces, geom)
    assert np.allclose(halves[fiber], wavelengths[fiber])
