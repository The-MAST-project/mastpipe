"""Tests for the wavelength model: the grating model's pixel/wavelength round trip, its
agreement with the simulator's truth, parameter recovery by ``fit_model`` from a perturbed
start, the anchor solution, the distortion powers and the serialisation round trips."""

import math

import numpy as np
import pandas as pd
import pytest

from mastpipe.core.model import WavelengthSolution
from mastpipe.steps.geometry import GratingGeometry, grating_geometry, grating_x, nominal_solution
from mastpipe.steps.wavemodel import (
    GratingModel,
    anchor_model,
    fit_model,
    model_from_config,
    model_from_solution,
    poly_powers,
)

N_PIX = 1598
BAND = (513.5, 527.0)


@pytest.fixture
def geom(highspec_effective) -> GratingGeometry:
    return grating_geometry(highspec_effective.geometry, highspec_effective.detector.pixel_size_um)


@pytest.fixture
def nominal(geom: GratingGeometry) -> GratingModel:
    return model_from_config("grating", geom, N_PIX)


# ---------------------------------------------------------------------------- the model


def test_model_from_config_is_nominal(geom: GratingGeometry, nominal: GratingModel):
    x_c, B = nominal_solution(geom, N_PIX)
    assert nominal.kind == "grating"
    assert nominal.n_pix == N_PIX
    assert nominal.x_c == pytest.approx(x_c)
    assert pytest.approx(B) == nominal.B
    assert nominal.alpha == pytest.approx(geom.alpha_rad)
    assert nominal.poly is None
    np.testing.assert_allclose(nominal.params, [x_c, B, geom.alpha_rad])
    assert nominal.wavelength(N_PIX / 2) == pytest.approx(geom.lambda_center, abs=1e-9)


def test_grating_model_roundtrip_pixel_wavelength(geom: GratingGeometry, nominal: GratingModel):
    x = np.linspace(5.0, N_PIX - 5.0, 41)
    lam = nominal.wavelength(x)
    assert lam.shape == x.shape
    np.testing.assert_allclose(nominal.pixel(lam), x, atol=1e-6)
    # with a distortion term (cubic and quadratic, no constant or linear) the inverse must
    # still land on the same pixel
    distorted = GratingModel(nominal.x_c, nominal.B, nominal.alpha, np.array([0.4, -0.3, 0.0, 0.0]), N_PIX, geom)
    lam_d = distorted.wavelength(x)
    assert np.max(np.abs(lam_d - lam)) > 1e-4  # the distortion is not a no-op
    np.testing.assert_allclose(distorted.pixel(lam_d), x, atol=1e-6)
    # scalars in, floats out
    assert isinstance(distorted.wavelength(700.0), float)
    assert isinstance(distorted.pixel(518.0), float)
    assert distorted.pixel(distorted.wavelength(700.0)) == pytest.approx(700.0, abs=1e-6)


def test_wavelength_matches_notebook_formula(geom: GratingGeometry, nominal: GratingModel):
    poly = np.array([0.5, 0.2, 0.0, 0.0])
    model = nominal.with_params(np.concatenate([nominal.params, poly]))
    x = np.array([100.0, 800.0, 1500.0])
    t = (x - N_PIX / 2) / (N_PIX / 2)
    corr = np.polyval(poly, t)
    theta = np.arctan((x - corr - model.x_c) / model.B)
    expected = geom.d_nm * (math.sin(model.alpha) + np.sin(model.alpha + theta)) / geom.order
    np.testing.assert_allclose(model.wavelength(x), expected, rtol=0, atol=1e-12)


def test_dispersion_is_central_difference(geom: GratingGeometry, nominal: GratingModel):
    x = np.array([200.0, N_PIX / 2, 1400.0])
    expected = np.abs(nominal.wavelength(x + 0.5) - nominal.wavelength(x - 0.5))
    np.testing.assert_allclose(nominal.dispersion(x), expected, atol=1e-12)
    # close to the plate scale of the geometry, positive whatever the dispersion sign
    assert nominal.dispersion(N_PIX / 2) == pytest.approx(geom.dispersion_nm_px, rel=0.05)
    flipped = GratingModel(nominal.x_c, -nominal.B, nominal.alpha, None, N_PIX, geom)
    assert np.all(flipped.dispersion(x) > 0)


def test_grating_model_matches_simulator_truth(sim_night, highspec_effective):
    _, truth = sim_night
    geom = grating_geometry(highspec_effective.geometry, highspec_effective.detector.pixel_size_um)
    n_pix = truth.science_flux.shape[-1]
    on = truth.lines_on_chip
    for fiber in (0, 17, 44):  # every trace has its own (x_c, B)
        model = GratingModel(truth.x_c_per_trace[fiber], truth.B_per_trace[fiber], geom.alpha_rad, None, n_pix, geom)
        pixels = truth.line_pixels[fiber]
        np.testing.assert_allclose(model.pixel(truth.lines_nm), pixels, atol=1e-6)
        np.testing.assert_allclose(model.wavelength(pixels[on]), truth.lines_nm[on], atol=1e-9)


def test_with_params_and_dict_roundtrip(geom: GratingGeometry, nominal: GratingModel):
    poly = np.array([0.1, -0.2, 0.0, 0.0])
    params = np.array([nominal.x_c + 2.0, nominal.B * 1.001, nominal.alpha + 1e-4, *poly])
    model = nominal.with_params(params)
    assert model is not nominal
    assert nominal.poly is None  # the original is untouched
    np.testing.assert_allclose(model.params, params)
    np.testing.assert_allclose(model.poly, poly)
    assert model.n_pix == N_PIX
    d = model.to_dict()
    assert d["kind"] == "grating"
    assert d["n_pix"] == N_PIX
    assert isinstance(d["x_c"], float)
    assert isinstance(d["poly"], list)
    back = GratingModel.from_dict(d, geom)
    np.testing.assert_allclose(back.params, model.params)
    x = np.linspace(10.0, 1500.0, 7)
    np.testing.assert_allclose(back.wavelength(x), model.wavelength(x), atol=1e-12)
    # three parameters mean no distortion term
    plain = nominal.with_params(nominal.params[:3])
    assert plain.poly is None
    assert GratingModel.from_dict(plain.to_dict(), geom).poly is None


def test_model_from_solution(geom: GratingGeometry, nominal: GratingModel):
    poly = np.array([0.25, 0.0, 0.0, 0.0])
    model = nominal.with_params(np.concatenate([nominal.params, poly]))
    lines = pd.DataFrame({"x": [100.0, 900.0], "wavelength": [514.0, 521.0], "species": ["Th I", "Ar II"]})
    lines["dx_px"] = 0.0
    solution = WavelengthSolution(
        fiber=0,
        kind=model.kind,
        params=model.params,
        n_pix=model.n_pix,
        lines=lines,
        rms_px=0.01,
        n_lines=2,
        anchor_px=(100.0, 900.0),
        anchor_nm=(514.0, 521.0),
        trace_y_center=17.4,
        geometry=geom.__dict__,
    )
    rebuilt = model_from_solution(solution, geom)
    assert isinstance(rebuilt, GratingModel)
    np.testing.assert_allclose(rebuilt.params, model.params)
    x = np.linspace(0.0, N_PIX, 9)
    np.testing.assert_allclose(rebuilt.wavelength(x), model.wavelength(x), atol=1e-12)
    bad = WavelengthSolution(
        fiber=0,
        kind="dispersion",
        params=np.zeros(3),
        n_pix=N_PIX,
        lines=lines,
        rms_px=0.0,
        n_lines=2,
        anchor_px=(0.0, 1.0),
        anchor_nm=(0.0, 1.0),
        trace_y_center=0.0,
        geometry={},
    )
    with pytest.raises(NotImplementedError):
        model_from_solution(bad, geom)


def test_dispersion_model_not_implemented(geom: GratingGeometry):
    with pytest.raises(NotImplementedError):
        model_from_config("dispersion", geom, N_PIX)
    with pytest.raises(ValueError):
        model_from_config("polynomial", geom, N_PIX)


# ----------------------------------------------------------------------------- fitting


def test_poly_powers_excludes_low_orders():
    assert poly_powers(False, 3) == [3, 2]
    assert poly_powers(True, 3) == [3]
    assert poly_powers(False, 4) == [4, 3, 2]
    assert poly_powers(True, 4) == [4, 3]
    assert poly_powers(False, 2) == [2]
    assert poly_powers(True, 2) == []
    assert poly_powers(False, 1) == []
    assert poly_powers(False, 0) == []


def _lines(model: GratingModel, n: int, noise_px: float, seed: int = 1):
    rng = np.random.default_rng(seed)
    lam = np.linspace(BAND[0], BAND[1], n)
    x = model.pixel(lam) + rng.normal(0.0, noise_px, n)
    return x, lam


def test_fit_model_recovers_params_from_lines(geom: GratingGeometry, nominal: GratingModel):
    truth = GratingModel(nominal.x_c + 4.0, nominal.B * 0.995, nominal.alpha, None, N_PIX, geom)
    x, lam = _lines(truth, 25, 0.02)
    start = GratingModel(truth.x_c - 3.0, truth.B * 1.01, truth.alpha, None, N_PIX, geom)
    fitted, resid = fit_model(start, x, lam)
    assert isinstance(fitted, GratingModel)
    assert resid.shape == (25,)
    assert fitted.x_c == pytest.approx(truth.x_c, abs=0.05)
    assert pytest.approx(truth.B, rel=1e-4) == fitted.B
    assert fitted.alpha == truth.alpha  # alpha held
    assert fitted.poly is None
    assert np.sqrt(np.mean(resid**2)) < 0.03
    np.testing.assert_allclose(resid, fitted.pixel(lam) - x, atol=1e-9)
    # the start is unchanged
    assert start.x_c == pytest.approx(truth.x_c - 3.0)


def test_fit_model_with_distortion_and_weights(geom: GratingGeometry, nominal: GratingModel):
    poly = np.array([0.6, -0.4, 0.0, 0.0])
    truth = GratingModel(nominal.x_c + 1.0, nominal.B * 1.002, nominal.alpha, poly, N_PIX, geom)
    x, lam = _lines(truth, 30, 0.0)
    fitted, resid = fit_model(nominal, x, lam, poly_deg=3)
    assert fitted.poly is not None
    assert fitted.poly.shape == (4,)
    np.testing.assert_allclose(fitted.poly[2:], 0.0)  # constant and linear never fitted
    np.testing.assert_allclose(fitted.poly[:2], poly[:2], atol=1e-3)
    assert np.max(np.abs(resid)) < 1e-3
    # weights of zero drop a line from the fit entirely
    x_bad = x.copy()
    x_bad[10] += 5.0
    weights = np.ones_like(x)
    weights[10] = 0.0
    fitted_w, resid_w = fit_model(nominal, x_bad, lam, weights=weights, poly_deg=3)
    assert np.max(np.abs(np.delete(resid_w, 10))) < 1e-3
    assert resid_w[10] == pytest.approx(-5.0, abs=1e-3)
    np.testing.assert_allclose(fitted_w.params, fitted.params, rtol=1e-6, atol=1e-6)


def test_fit_model_free_alpha_and_bounds(geom: GratingGeometry, nominal: GratingModel):
    alpha_true = nominal.alpha + math.radians(0.3)
    truth = GratingModel(nominal.x_c - 20.0, nominal.B * 0.99, alpha_true, None, N_PIX, geom)
    x, lam = _lines(truth, 30, 0.0)
    fitted, resid = fit_model(nominal, x, lam, fit_alpha=True, poly_deg=3)
    # x_c, B, alpha and a cubic P whose quadratic, linear and constant terms stay zero
    assert fitted.params.shape == (7,)
    assert fitted.poly is not None and fitted.poly.shape == (4,)
    np.testing.assert_allclose(fitted.poly[1:], 0.0)
    assert fitted.alpha == pytest.approx(alpha_true, abs=1e-6)
    assert np.max(np.abs(resid)) < 1e-3
    # bounds: alpha pinned to within 0.1 deg of the design angle, B to within 0.2 percent
    bounded, _ = fit_model(nominal, x, lam, fit_alpha=True, alpha_bound_deg=0.1, B_bound_frac=0.002)
    assert abs(bounded.alpha - geom.alpha_rad) <= math.radians(0.1) + 1e-12
    assert abs(bounded.B) <= geom.B_nom * 1.002 + 1e-9
    assert abs(bounded.B) >= geom.B_nom * 0.998 - 1e-9


# ----------------------------------------------------------------------------- anchors


def test_anchor_model_two_pairs(geom: GratingGeometry, nominal: GratingModel):
    truth = GratingModel(nominal.x_c + 7.0, nominal.B * 0.997, nominal.alpha, None, N_PIX, geom)
    waves = (515.8605, 525.8360)
    pairs = [(w, truth.pixel(w)) for w in waves]
    model = anchor_model(pairs, geom, N_PIX)
    assert isinstance(model, GratingModel)
    assert model.x_c == pytest.approx(truth.x_c, abs=1e-6)
    assert pytest.approx(truth.B, rel=1e-9) == model.B
    assert model.alpha == geom.alpha_rad
    assert model.poly is None
    assert model.n_pix == N_PIX
    # one pair fixes x_c at the nominal B
    single = anchor_model(pairs[:1], geom, N_PIX)
    assert pytest.approx(nominal.B) == single.B
    assert single.pixel(waves[0]) == pytest.approx(pairs[0][1], abs=1e-9)
    # three pairs: least squares through all of them
    three = pairs + [(520.0, truth.pixel(520.0))]
    model3 = anchor_model(three, geom, N_PIX)
    assert model3.x_c == pytest.approx(truth.x_c, abs=1e-6)
    assert grating_x(520.0, model3.x_c, model3.B, model3.alpha, geom) == pytest.approx(three[2][1], abs=1e-6)
    with pytest.raises(ValueError):
        anchor_model([], geom, N_PIX)
