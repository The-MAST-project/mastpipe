"""Tests for the full-arc wavelength solve: line identification and the cross-correlation
shift on synthetic inputs, ``solve_fiber`` and ``solve_full_arc`` on the simulated full arc
against the simulator's grating truth, seeding from a previous reference, missing fibers,
the fiber-to-fiber agreement check and the reference hash."""

import hashlib
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.models import WavecalParams
from mastpipe.core.model import LINE_COLUMNS, ReferenceSolution, WavelengthSolution
from mastpipe.steps import wavecal
from mastpipe.steps.atlas import catalogue_arrays
from mastpipe.steps.wavecal import (
    build_reference,
    fiber_agreement,
    identify_lines,
    shift_to_reference,
    solve_fiber,
    solve_full_arc,
    wavelength_map,
)
from mastpipe.steps.wavemodel import GratingModel, anchor_model, model_from_solution

MIN_SOLVED = 43
MAX_RMS_PX = 0.05
MIN_LINES = 20
MAX_MAP_ERROR_DISP = 0.3
"""Largest tolerated ``|lambda_fit - lambda_true|`` anywhere on the detector, in dispersions."""
ANCHOR_LOG_INTENSITY = 2.0
CREATED = "2026-10-29T20:15:00"


# ------------------------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def arc(sim_full_arc, sim_arc_lines, sim_arc_fluxes, mgb518_lines) -> SimpleNamespace:
    """The simulated full arc's inputs to the solve plus the grating truth."""
    truth = sim_full_arc.truth
    geom = sim_full_arc.geom
    cat_wl, cat_sp = catalogue_arrays(mgb518_lines, usable_only=True)
    truth_model = GratingModel(truth.x_c, truth.B, geom.alpha_rad, None, sim_full_arc.n_pix, geom)
    return SimpleNamespace(
        truth=truth,
        geom=geom,
        n_pix=sim_full_arc.n_pix,
        traces=sim_full_arc.traces,
        tables=sim_arc_lines[0],
        fluxes=sim_arc_fluxes,
        cat_wl=cat_wl,
        cat_sp=cat_sp,
        truth_model=truth_model,
        params=WavecalParams(),
    )


@pytest.fixture(scope="module")
def anchor(arc) -> GratingModel:
    """The anchor model of fiber 0: two bright on-chip catalogue lines at their truth pixels,
    rounded to a tenth of a pixel as if read off the anchor panels."""
    truth = arc.truth
    usable = np.isin(truth.lines_nm, arc.cat_wl)
    bright = usable & truth.lines_on_chip & (truth.line_logI > ANCHOR_LOG_INTENSITY)
    assert bright.sum() >= 2
    first, last = np.flatnonzero(bright)[[0, -1]]
    pairs = [
        (float(truth.lines_nm[i]), round(float(truth.line_pixels[0][i]), 1)) for i in (first, last)
    ]
    return anchor_model(pairs, arc.geom, arc.n_pix)


@pytest.fixture(scope="module")
def solved(arc, anchor) -> dict[int, WavelengthSolution]:
    return solve_full_arc(
        arc.tables, arc.fluxes, arc.traces, anchor, arc.cat_wl, arc.cat_sp, arc.params
    )


def x_obs_of(table: pd.DataFrame) -> np.ndarray:
    return table.loc[table["use_x"].to_numpy(dtype=bool), "x_center"].to_numpy(dtype=float)


def map_error_in_dispersions(solution: WavelengthSolution, arc) -> float:
    lam = wavelength_map(solution, arc.geom)
    x = np.arange(arc.n_pix)
    return float(
        np.max(np.abs(lam - arc.truth_model.wavelength(x)) / arc.truth_model.dispersion(x))
    )


# ------------------------------------------------------------------- identify / shift


def test_identify_lines_nearest_within_tolerance_and_unique(arc):
    model = arc.truth_model
    disp = model.dispersion(arc.n_pix / 2)
    cat_wl = np.array([515.0, 518.0, 518.0 + 2.5 * disp, 522.0])
    cat_sp = np.array(["Th I", "Ar II", "Ar I", "Th II"], dtype=object)
    # a line 0.5 px off 515.0, one exactly on 522.0, one 3 px off (rejected by tol), one at
    # 518.0 whose runner-up 2.5 px away is too close (rejected as ambiguous)
    x_obs = np.array(
        [model.pixel(515.0) + 0.5, model.pixel(522.0), model.pixel(515.0) + 3.0, model.pixel(518.0)]
    )
    ident = identify_lines(x_obs, model, cat_wl, cat_sp, tol_px=2.0, unique_tol_px=4.0)
    assert list(ident.columns) == list(LINE_COLUMNS)
    assert ident["wavelength"].tolist() == [515.0, 522.0]
    assert ident["species"].tolist() == ["Th I", "Th II"]
    assert ident["x"].to_numpy() == pytest.approx(x_obs[:2])
    # dx_px is the offset of the model wavelength from the catalogue, in local dispersions,
    # signed by the dispersion direction
    sign = np.sign(model.B)
    assert ident["dx_px"].iloc[0] == pytest.approx(0.5 * sign, abs=1e-3)
    assert ident["dx_px"].iloc[1] == pytest.approx(0.0, abs=1e-9)
    # widening the uniqueness tolerance to nothing lets the 518.0 match through
    loose = identify_lines(x_obs, model, cat_wl, cat_sp, tol_px=2.0, unique_tol_px=0.0)
    assert loose["wavelength"].tolist() == [515.0, 522.0, 518.0]
    # no measured lines, or a catalogue of one line, must not break
    empty = identify_lines(np.array([]), model, cat_wl, cat_sp, 2.0, 4.0)
    assert len(empty) == 0 and list(empty.columns) == list(LINE_COLUMNS)
    single = identify_lines(x_obs[:1], model, cat_wl[:1], cat_sp[:1], 2.0, 4.0)
    assert single["wavelength"].tolist() == [515.0]


def test_shift_to_reference_recovers_subpixel_shift():
    rng = np.random.default_rng(5)
    x = np.arange(1600, dtype=float)
    centres = rng.uniform(50.0, 1550.0, 40)
    heights = rng.uniform(200.0, 5000.0, 40)

    def spectrum(offset: float) -> np.ndarray:
        flux = rng.normal(0.0, 3.0, x.size)
        for c, h in zip(centres, heights, strict=True):
            flux += h * np.exp(-0.5 * ((x - c - offset) / 1.3) ** 2)
        return flux

    ref = spectrum(0.0)
    assert shift_to_reference(ref, ref) == pytest.approx(0.0, abs=1e-9)
    assert shift_to_reference(spectrum(7.3), ref) == pytest.approx(7.3, abs=0.1)
    assert shift_to_reference(spectrum(-12.6), ref) == pytest.approx(-12.6, abs=0.1)
    # NaN columns (masked pixels) are tolerated
    shifted = spectrum(3.0)
    shifted[100:110] = np.nan
    assert shift_to_reference(shifted, ref) == pytest.approx(3.0, abs=0.1)
    # the search range bounds the answer
    assert abs(shift_to_reference(spectrum(50.0), ref, max_shift=20)) <= 20.0


# ------------------------------------------------------------------- the simulated arc


def test_solve_fiber_from_anchor_on_simulated_full_arc(arc, anchor):
    x_obs = x_obs_of(arc.tables[0])
    model, ident = solve_fiber(x_obs, anchor, arc.cat_wl, arc.cat_sp, arc.params)
    assert model is not None
    assert isinstance(model, GratingModel)
    assert len(ident) >= MIN_LINES
    assert list(ident.columns) == list(LINE_COLUMNS)
    resid = model.pixel(ident["wavelength"].to_numpy()) - ident["x"].to_numpy()
    rms = float(np.sqrt(np.mean(resid**2)))
    assert rms < MAX_RMS_PX
    # dx_px of the returned table is the final fit residual (fit_model evaluates the
    # distortion term at the measured column, pixel() at the model column: ~1e-6 px apart)
    np.testing.assert_allclose(ident["dx_px"].to_numpy(), resid, atol=1e-5)
    # every identified wavelength is a catalogue line and sits at its truth pixel
    assert set(ident["wavelength"]) <= set(arc.cat_wl)
    truth_px = arc.truth_model.pixel(ident["wavelength"].to_numpy())
    assert np.max(np.abs(truth_px - ident["x"].to_numpy())) < 0.2
    # the fitted map agrees with the truth across the detector
    x = np.arange(arc.n_pix)
    err = np.abs(model.wavelength(x) - arc.truth_model.wavelength(x)) / arc.truth_model.dispersion(
        x
    )
    assert err.max() < MAX_MAP_ERROR_DISP
    assert model.x_c == pytest.approx(arc.truth.x_c, abs=0.5)
    assert model.B == pytest.approx(arc.truth.B, rel=2e-3)


def test_solve_fiber_gives_up_below_min_lines(arc, anchor):
    x_obs = x_obs_of(arc.tables[0])[:4]
    model, ident = solve_fiber(x_obs, anchor, arc.cat_wl, arc.cat_sp, arc.params)
    assert model is None
    assert len(ident) <= 4
    model, ident = solve_fiber(np.array([]), anchor, arc.cat_wl, arc.cat_sp, arc.params)
    assert model is None and len(ident) == 0
    with pytest.raises(NotImplementedError):
        solve_fiber(x_obs, anchor, arc.cat_wl, arc.cat_sp, arc.params, final_model="poly")
    with pytest.raises(ValueError):
        solve_fiber(x_obs, anchor, arc.cat_wl, arc.cat_sp, arc.params, final_model="spline")


def test_solve_full_arc_all_45_from_single_model(arc, solved):
    assert len(solved) >= MIN_SOLVED
    assert set(solved) <= set(range(len(arc.traces)))
    for fiber, solution in solved.items():
        assert solution.fiber == fiber
        assert solution.kind == "grating"
        assert solution.n_pix == arc.n_pix
        assert solution.n_lines == len(solution.lines) >= arc.params.min_lines
        assert list(solution.lines.columns) == list(LINE_COLUMNS)
        assert solution.rms_px < MAX_RMS_PX
        assert solution.anchor_px == (solution.lines["x"].min(), solution.lines["x"].max())
        assert solution.anchor_nm == (
            solution.lines["wavelength"].min(),
            solution.lines["wavelength"].max(),
        )
        assert solution.trace_y_center == arc.traces.traces[fiber].y_center
        assert set(solution.geometry) == {"d_nm", "alpha_rad", "order", "B_nom"}
        assert solution.geometry["B_nom"] == arc.geom.B_nom
        assert map_error_in_dispersions(solution, arc) < MAX_MAP_ERROR_DISP
    # the wavelength map is the stored model evaluated on every pixel
    solution = solved[0]
    lam = wavelength_map(solution, arc.geom)
    assert lam.shape == (arc.n_pix,)
    rebuilt = model_from_solution(solution, arc.geom)
    np.testing.assert_allclose(lam, rebuilt.wavelength(np.arange(arc.n_pix)), atol=1e-12)
    # fibers restricts the solve
    subset = solve_full_arc(
        arc.tables,
        arc.fluxes,
        arc.traces,
        model_from_solution(solved[0], arc.geom),
        arc.cat_wl,
        arc.cat_sp,
        arc.params,
        fibers=[3, 5],
    )
    assert set(subset) <= {3, 5}


def test_solve_full_arc_from_previous_reference_is_faster_and_equal(
    arc, anchor, solved, monkeypatch
):
    reference = build_reference("highspec", "mgb518", "survey45", solved, CREATED, ("abc123",))
    calls: list[int] = []

    def counting_fit(*args, **kwargs):
        calls.append(1)
        return wavecal_fit(*args, **kwargs)

    wavecal_fit = wavecal.fit_model
    monkeypatch.setattr(wavecal, "fit_model", counting_fit)
    solve_full_arc(arc.tables, arc.fluxes, arc.traces, anchor, arc.cat_wl, arc.cat_sp, arc.params)
    n_from_anchor = len(calls)
    calls.clear()
    again = solve_full_arc(
        arc.tables,
        arc.fluxes,
        arc.traces,
        reference,
        arc.cat_wl,
        arc.cat_sp,
        arc.params,
        geom=arc.geom,
    )
    n_from_reference = len(calls)
    assert n_from_reference <= n_from_anchor
    assert set(again) == set(solved)
    for fiber, solution in solved.items():
        other = again[fiber]
        assert other.lines["wavelength"].tolist() == solution.lines["wavelength"].tolist()
        np.testing.assert_allclose(other.lines["x"], solution.lines["x"], atol=1e-9)
        # the two starts reach the same least-squares optimum to the solver's tolerance:
        # x_c within 1e-6 px, B within 1e-6 of itself, the wavelength map within 1e-6 nm
        assert other.params[0] == pytest.approx(solution.params[0], abs=1e-6)
        assert other.params[1] == pytest.approx(solution.params[1], rel=1e-6)
        assert other.params[2] == solution.params[2]
        np.testing.assert_allclose(
            wavelength_map(other, arc.geom), wavelength_map(solution, arc.geom), atol=1e-6
        )
    # a reference start needs the geometry to rebuild its models
    with pytest.raises(ValueError):
        solve_full_arc(
            arc.tables, arc.fluxes, arc.traces, reference, arc.cat_wl, arc.cat_sp, arc.params
        )


def test_missing_fiber_absent_not_crash(arc, anchor, solved):
    tables = list(arc.tables)
    tables[7] = tables[7].iloc[0:0]  # no detected lines at all
    tables[9] = tables[9].assign(use_x=False)  # lines, none usable along the dispersion
    result = solve_full_arc(
        tables, arc.fluxes, arc.traces, anchor, arc.cat_wl, arc.cat_sp, arc.params
    )
    assert 7 not in result and 9 not in result
    assert len(result) >= MIN_SOLVED - 2
    assert set(result) == set(solved) - {7, 9}
    # a fiber absent from the previous reference seeds from its nearest solved neighbour
    partial = {fiber: solution for fiber, solution in solved.items() if fiber != 7}
    reference = build_reference("highspec", "mgb518", "survey45", partial, CREATED, ())
    result = solve_full_arc(
        arc.tables,
        arc.fluxes,
        arc.traces,
        reference,
        arc.cat_wl,
        arc.cat_sp,
        arc.params,
        geom=arc.geom,
    )
    assert 7 in result
    assert map_error_in_dispersions(result[7], arc) < MAX_MAP_ERROR_DISP
    # the inputs must line up with the trace set
    with pytest.raises(ValueError):
        solve_full_arc(
            tables[:-1], arc.fluxes, arc.traces, anchor, arc.cat_wl, arc.cat_sp, arc.params
        )


# --------------------------------------------------------------- agreement and reference


def make_solution(
    fiber: int, lines: pd.DataFrame, params=(800.0, 20000.0, 0.75)
) -> WavelengthSolution:
    table = lines.copy()
    table["species"] = "Th I"
    return WavelengthSolution(
        fiber=fiber,
        kind="grating",
        params=np.asarray(params, dtype=float),
        n_pix=1598,
        lines=table[list(LINE_COLUMNS)],
        rms_px=0.01,
        n_lines=len(table),
        anchor_px=(float(table["x"].min()), float(table["x"].max())),
        anchor_nm=(float(table["wavelength"].min()), float(table["wavelength"].max())),
        trace_y_center=10.0 * fiber,
        geometry={"d_nm": 1000.0, "alpha_rad": 0.75, "order": 1, "B_nom": 20000.0},
    )


def test_fiber_agreement_small_on_truth(arc, solved):
    agreement = fiber_agreement(solved)
    assert list(agreement.columns) == ["wavelength", "n_fibers", "spread_px"]
    assert (agreement["n_fibers"] >= 2).all()
    assert agreement["n_fibers"].max() == len(solved)
    assert set(agreement["wavelength"]) <= set(arc.cat_wl)
    assert agreement["wavelength"].is_monotonic_increasing
    assert agreement["spread_px"].median() < 0.1
    assert agreement["spread_px"].max() < MAX_MAP_ERROR_DISP
    # hand-built: the spread is max - min of dx_px over the fibers that hold the line, and a
    # line seen in one fiber only is left out
    a = pd.DataFrame({"x": [100.0, 200.0, 300.0], "wavelength": [514.0, 516.0, 518.0]})
    a["dx_px"] = [0.01, -0.02, 0.05]
    b = pd.DataFrame({"x": [101.0, 201.0], "wavelength": [514.0, 516.0]})
    b["dx_px"] = [-0.03, 0.02]
    c = pd.DataFrame({"x": [102.0], "wavelength": [514.0]})
    c["dx_px"] = [0.0]
    table = fiber_agreement(
        {0: make_solution(0, a), 1: make_solution(1, b), 2: make_solution(2, c)}
    )
    assert table["wavelength"].tolist() == [514.0, 516.0]
    assert table["n_fibers"].tolist() == [3, 2]
    assert table["spread_px"].to_numpy() == pytest.approx([0.04, 0.04])
    assert len(fiber_agreement({0: make_solution(0, a)})) == 0


def test_reference_hash_changes_with_params():
    lines = pd.DataFrame({"x": [100.0, 900.0], "wavelength": [514.0, 521.0], "dx_px": [0.0, 0.0]})
    solutions = {0: make_solution(0, lines), 1: make_solution(1, lines)}
    reference = build_reference("highspec", "mgb518", "survey45", solutions, CREATED, ("c1", "c2"))
    assert isinstance(reference, ReferenceSolution)
    assert reference.identifier == "highspec_mgb518_survey45_20261029T201500"
    assert reference.valid_from == reference.created == CREATED
    assert reference.source_checksums == ("c1", "c2")
    assert reference.solutions is solutions
    assert len(reference.hash) == 12 and int(reference.hash, 16) >= 0
    # deterministic, independent of the created time and checksums
    same = build_reference("x", "y", "z", solutions, "2027-01-01T00:00:00", ())
    assert same.hash == reference.hash
    # a change of parameters, of a line position or of a line's wavelength changes it
    moved = {0: make_solution(0, lines, params=(800.001, 20000.0, 0.75)), 1: solutions[1]}
    assert build_reference("x", "y", "z", moved, CREATED, ()).hash != reference.hash
    shifted = {0: make_solution(0, lines.assign(x=[100.5, 900.0])), 1: solutions[1]}
    assert build_reference("x", "y", "z", shifted, CREATED, ()).hash != reference.hash
    relabelled = {0: make_solution(0, lines.assign(wavelength=[514.1, 521.0])), 1: solutions[1]}
    assert build_reference("x", "y", "z", relabelled, CREATED, ()).hash != reference.hash
    # the fiber number is part of it; a change below the 1e-9 rounding is not
    renumbered = {0: solutions[0], 2: make_solution(2, lines)}
    assert build_reference("x", "y", "z", renumbered, CREATED, ()).hash != reference.hash
    tiny = {0: make_solution(0, lines, params=(800.0 + 1e-12, 20000.0, 0.75)), 1: solutions[1]}
    assert build_reference("x", "y", "z", tiny, CREATED, ()).hash == reference.hash
    assert hashlib.sha256  # the digest family the hash is drawn from
