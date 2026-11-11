"""Tests for the grating geometry step module."""

import math

import numpy as np
import pytest

from mastpipe.config.models import GeometryConfig
from mastpipe.steps.geometry import (
    GratingGeometry,
    grating_geometry,
    grating_x,
    lambda_of_x,
    nominal_solution,
    theta_of_lambda,
)

N_PIX = 1598


@pytest.fixture
def mg(highspec_effective) -> GratingGeometry:
    return grating_geometry(highspec_effective.geometry, highspec_effective.detector.pixel_size_um)


def test_nominal_dispersion_matches_notebook(mg: GratingGeometry):
    assert mg.dispersion_nm_px == pytest.approx(0.00866, abs=2e-4)
    assert mg.lambda_littrow == pytest.approx(519.0, abs=0.5)
    assert mg.B_nom == pytest.approx(479.41 * 1000.0 / 16.0)
    assert mg.d_nm == pytest.approx(367.0)
    assert mg.alpha_rad == pytest.approx(math.radians(45.0))
    assert mg.order == 1
    assert mg.lambda_center == pytest.approx(518.0)
    assert mg.dispersion_sign == 1


def test_grating_x_inverts_theta(mg: GratingGeometry):
    g = mg
    x_c, B = nominal_solution(g, N_PIX)
    for lam in (513.5, 518.0, 524.0):
        x = grating_x(lam, x_c, B, g.alpha_rad, g)
        assert x == pytest.approx(x_c + B * math.tan(theta_of_lambda(lam, g)), abs=1e-9)
        lam_back = lambda_of_x(x, x_c, B, g.alpha_rad, g)
        assert lam_back == pytest.approx(lam, abs=1e-9)
        assert grating_x(lam_back, x_c, B, g.alpha_rad, g) == pytest.approx(x, abs=1e-9)
    lams = np.array([513.5, 518.0, 524.0])
    xs = grating_x(lams, x_c, B, g.alpha_rad, g)
    assert xs.shape == (3,)
    np.testing.assert_allclose(lambda_of_x(xs, x_c, B, g.alpha_rad, g), lams, atol=1e-9)
    # at the Littrow wavelength the diffracted beam returns along the incidence direction
    assert theta_of_lambda(g.lambda_littrow, g) == pytest.approx(0.0, abs=1e-12)
    assert theta_of_lambda(g.lambda_littrow, g, alpha=math.radians(46.0)) != pytest.approx(0.0)


def test_theta_outside_grating_range_raises(mg: GratingGeometry):
    with pytest.raises(ValueError):
        theta_of_lambda(700.0, mg)
    with pytest.raises(ValueError):
        theta_of_lambda(np.array([518.0, 700.0]), mg)
    with pytest.raises(ValueError):
        grating_x(-200.0, 800.0, mg.B_nom, mg.alpha_rad, mg)


def test_nominal_solution_centres_lambda_center(highspec_effective, mg: GratingGeometry):
    g = mg
    x_c, B = nominal_solution(g, N_PIX)
    assert B == pytest.approx(g.B_nom)
    assert grating_x(g.lambda_center, x_c, B, g.alpha_rad, g) == pytest.approx(N_PIX / 2)
    # wavelength increases with x for these constants
    assert grating_x(520.0, x_c, B, g.alpha_rad, g) > grating_x(516.0, x_c, B, g.alpha_rad, g)
    # the dispersion sign flips B and the direction, the centre stays put
    flipped_config: GeometryConfig = highspec_effective.geometry.model_copy(update={"dispersion_sign": -1})
    flipped = grating_geometry(flipped_config, highspec_effective.detector.pixel_size_um)
    x_c2, B2 = nominal_solution(flipped, N_PIX)
    assert B2 == pytest.approx(-g.B_nom)
    assert grating_x(g.lambda_center, x_c2, B2, g.alpha_rad, g) == pytest.approx(N_PIX / 2)
    assert grating_x(520.0, x_c2, B2, g.alpha_rad, g) < grating_x(516.0, x_c2, B2, g.alpha_rad, g)
