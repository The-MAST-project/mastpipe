"""Grating geometry: the dispersion model of a HighSpec configuration.

A port of the grating equation used in the HighSpec laboratory reduction notebook. The
detector sits in the focal plane of a camera of focal length ``f_cam`` behind a grating of
groove spacing ``d`` illuminated at incidence angle ``alpha`` in order ``m``; a wavelength
``lam`` leaves the grating at ``theta_of_lambda(lam)`` from the incidence direction and
lands at ``grating_x(lam) = x_c + B tan(theta)`` on the detector, where ``x_c`` is the pixel
of the returning Littrow beam and ``B`` the camera focal length in pixels. The nominal
geometry has ``B = B_nom = f_cam / pixel`` and the configuration's central wavelength at the
middle of the frame; the wavelength solution of a fiber later fits ``(x_c, B)``.

Everything here is pure: plain numbers in, plain numbers or arrays out.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from mastpipe.config.models import GeometryConfig


@dataclass(frozen=True)
class GratingGeometry:
    """The constants of one grating configuration and its nominal plate scale.

    ``dispersion_sign`` is +1 when wavelength increases with ``x`` on the oriented frame and
    -1 otherwise; :func:`nominal_solution` folds it into the sign of ``B``.
    """

    d_nm: float
    alpha_rad: float
    order: int
    B_nom: float
    lambda_littrow: float
    dispersion_nm_px: float
    lambda_center: float
    dispersion_sign: int = 1


def grating_geometry(geom: GeometryConfig, pixel_size_um: float) -> GratingGeometry:
    """Nominal plate scale, Littrow wavelength and dispersion of a grating configuration.

    ``B_nom = f_cam / pixel`` (pixels), ``lambda_littrow = 2 d sin(alpha) / m`` and
    ``dispersion = d cos(alpha) / (m B_nom)`` (nm per pixel at Littrow). ``ValueError`` for
    a geometry whose ``model`` is not ``"grating"``.
    """
    if geom.model != "grating":
        raise ValueError(f"grating geometry needs model 'grating', got {geom.model!r}")
    d = geom.groove_spacing_nm
    alpha = math.radians(geom.incidence_angle_deg)
    m = geom.order
    B_nom = geom.camera_focal_mm * 1000.0 / pixel_size_um
    return GratingGeometry(
        d_nm=d,
        alpha_rad=alpha,
        order=m,
        B_nom=B_nom,
        lambda_littrow=2.0 * d * math.sin(alpha) / m,
        dispersion_nm_px=d * math.cos(alpha) / (m * B_nom),
        lambda_center=geom.central_wavelength_nm,
        dispersion_sign=geom.dispersion_sign,
    )


def theta_of_lambda(lam, g: GratingGeometry, alpha: float | None = None):
    """Diffraction angle from the incidence direction, ``arcsin(m lam / d - sin alpha) - alpha``.

    ``lam`` is a scalar or array of air wavelengths in nm; ``alpha`` overrides the
    incidence angle of ``g``. ``ValueError`` when any wavelength is outside the range the
    grating can diffract into this order.
    """
    if alpha is None:
        alpha = g.alpha_rad
    s = np.asarray(lam, dtype=float) * g.order / g.d_nm - math.sin(alpha)
    if np.any(np.abs(s) > 1.0):
        raise ValueError(f"wavelength outside the grating range of order {g.order}: {lam!r}")
    theta = np.arcsin(s) - alpha
    return float(theta) if theta.ndim == 0 else theta


def grating_x(lam, x_c: float, B: float, alpha: float, g: GratingGeometry):
    """Detector column of wavelength ``lam``: ``x_c + B tan(theta_of_lambda(lam))``."""
    x = x_c + B * np.tan(theta_of_lambda(lam, g, alpha))
    return float(x) if np.ndim(x) == 0 else x


def lambda_of_x(x, x_c: float, B: float, alpha: float, g: GratingGeometry):
    """The inverse of :func:`grating_x`: wavelength landing at column ``x``."""
    theta = np.arctan((np.asarray(x, dtype=float) - x_c) / B)
    lam = g.d_nm * (np.sin(theta + alpha) + math.sin(alpha)) / g.order
    return float(lam) if lam.ndim == 0 else lam


def nominal_solution(g: GratingGeometry, n_pix: int) -> tuple[float, float]:
    """``(x_c, B)`` of the nominal geometry: ``B = B_nom`` (signed by ``dispersion_sign``)
    and ``x_c`` such that ``lambda_center`` lands at ``n_pix / 2``."""
    B = g.B_nom * g.dispersion_sign
    x_c = n_pix / 2.0 - B * math.tan(theta_of_lambda(g.lambda_center, g))
    return x_c, B
