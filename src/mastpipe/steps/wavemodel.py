"""Wavelength models: the map between detector column and wavelength of one fiber.

A :class:`WavelengthModel` answers ``wavelength(x)``, ``pixel(lam)`` and ``dispersion(x)``
and is a plain parameter vector for storage. :class:`GratingModel` is the HighSpec model, a
port of the wavelength-solution functions of the laboratory reduction notebook: the grating
equation of :mod:`mastpipe.steps.geometry` plus an optional distortion polynomial ``P(t)``
in the normalised column ``t = (x - n/2) / (n/2)``,

    lam(x) = d (sin alpha + sin(alpha + theta)) / m,   theta = atan((x - P(t) - x_c) / B).

``P`` carries no constant or linear term (those are ``x_c`` and ``B``; giving them to ``P``
makes the fit degenerate) and no quadratic term when ``alpha`` is free (a change of
incidence angle is a quadratic in ``x`` to first order); :func:`poly_powers` enforces this.
``params`` is ``[x_c, B, alpha, *poly]`` with ``poly`` in ``numpy.polyval`` order.

:func:`fit_model` is the notebook's ``fit_solution`` (least squares on ``x_c, B[, alpha][,
P]`` with residuals in pixels) and :func:`anchor_model` its ``fit_anchor_solution``. The
DeepSpec dispersion model belongs to a later plan; ``model_from_config("dispersion", ...)``
raises ``NotImplementedError``.

Everything here is pure: arrays and plain values in, models and arrays out.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from typing import Any, ClassVar, Protocol, runtime_checkable

import numpy as np
from scipy.optimize import least_squares

from mastpipe.core.model import WavelengthSolution
from mastpipe.steps.geometry import GratingGeometry, grating_x, nominal_solution, theta_of_lambda

NEWTON_STEPS = 2
"""Newton iterations that invert the distortion term in :meth:`GratingModel.pixel`."""


@runtime_checkable
class WavelengthModel(Protocol):
    """What every wavelength model offers; ``kind`` names the family for storage."""

    kind: ClassVar[str]
    n_pix: int

    @property
    def params(self) -> np.ndarray: ...

    def wavelength(self, x): ...

    def pixel(self, lam): ...

    def dispersion(self, x): ...

    def with_params(self, params: np.ndarray) -> WavelengthModel: ...

    def to_dict(self) -> dict[str, Any]: ...


def tnorm(x, n_pix: int):
    """The distortion variable: column normalised to ``[-1, 1]`` across the detector."""
    return (np.asarray(x, dtype=float) - n_pix / 2.0) / (n_pix / 2.0)


def _scalar_or_array(value: np.ndarray):
    return float(value) if value.ndim == 0 else value


@dataclass(frozen=True, eq=False)
class GratingModel:
    """The grating equation with an optional distortion polynomial.

    ``x_c`` is the column of the returning Littrow beam, ``B`` the camera focal length in
    pixels (signed by the dispersion direction), ``alpha`` the incidence angle in radians
    and ``poly`` the distortion coefficients in ``numpy.polyval`` order (``None`` for no
    distortion). ``geom`` supplies the groove spacing and order.
    """

    kind: ClassVar[str] = "grating"

    x_c: float
    B: float
    alpha: float
    poly: np.ndarray | None
    n_pix: int
    geom: GratingGeometry

    def __post_init__(self) -> None:
        poly = None if self.poly is None else np.asarray(self.poly, dtype=float)
        if poly is not None and poly.ndim != 1:
            raise ValueError(f"poly must be a 1-D coefficient array, got shape {poly.shape}")
        object.__setattr__(self, "poly", poly)
        object.__setattr__(self, "x_c", float(self.x_c))
        object.__setattr__(self, "B", float(self.B))
        object.__setattr__(self, "alpha", float(self.alpha))
        object.__setattr__(self, "n_pix", int(self.n_pix))

    @property
    def params(self) -> np.ndarray:
        """``[x_c, B, alpha, *poly]``."""
        poly = [] if self.poly is None else list(self.poly)
        return np.array([self.x_c, self.B, self.alpha, *poly], dtype=float)

    @classmethod
    def from_params(cls, params, n_pix: int, geom: GratingGeometry) -> GratingModel:
        """The model with ``params`` ``[x_c, B, alpha, *poly]``; three values mean no ``poly``."""
        p = np.asarray(params, dtype=float)
        if p.ndim != 1 or p.size < 3:
            raise ValueError(f"grating params need at least [x_c, B, alpha], got {p!r}")
        poly = p[3:].copy() if p.size > 3 else None
        return cls(p[0], p[1], p[2], poly, n_pix, geom)

    def with_params(self, params) -> GratingModel:
        """A new model with the same ``n_pix`` and geometry and these parameters."""
        return self.from_params(params, self.n_pix, self.geom)

    def distortion(self, x):
        """``P(t(x))`` in pixels (zero without a distortion term)."""
        x = np.asarray(x, dtype=float)
        if self.poly is None:
            return np.zeros_like(x)
        return np.polyval(self.poly, tnorm(x, self.n_pix))

    def wavelength(self, x):
        """Wavelength (air nm) landing at column(s) ``x``."""
        x = np.asarray(x, dtype=float)
        theta = np.arctan((x - self.distortion(x) - self.x_c) / self.B)
        g = self.geom
        lam = g.d_nm * (math.sin(self.alpha) + np.sin(self.alpha + theta)) / g.order
        return _scalar_or_array(lam)

    def pixel(self, lam):
        """Column(s) of wavelength(s) ``lam``: the grating column plus the distortion term,
        the latter solved by :data:`NEWTON_STEPS` Newton steps of ``x - P(t(x)) = x_grating``.
        """
        x_grating = np.asarray(grating_x(lam, self.x_c, self.B, self.alpha, self.geom), float)
        if self.poly is None:
            return _scalar_or_array(x_grating)
        x = x_grating + self.distortion(x_grating)
        derivative = np.polyder(self.poly) if self.poly.size > 1 else np.zeros(1)
        for _ in range(NEWTON_STEPS):
            t = tnorm(x, self.n_pix)
            residual = x - np.polyval(self.poly, t) - x_grating
            slope = 1.0 - np.polyval(derivative, t) * 2.0 / self.n_pix
            x = x - residual / slope
        return _scalar_or_array(x)

    def dispersion(self, x):
        """Local dispersion in nm per pixel, ``|lam(x + 1/2) - lam(x - 1/2)|``."""
        x = np.asarray(x, dtype=float)
        return _scalar_or_array(np.abs(self.wavelength(x + 0.5) - self.wavelength(x - 0.5)))

    def to_dict(self) -> dict[str, Any]:
        """JSON-native form; ``from_dict`` with the geometry inverts it."""
        return {
            "kind": self.kind,
            "x_c": self.x_c,
            "B": self.B,
            "alpha": self.alpha,
            "poly": None if self.poly is None else [float(c) for c in self.poly],
            "n_pix": self.n_pix,
        }

    @classmethod
    def from_dict(cls, d: dict[str, Any], geometry: GratingGeometry) -> GratingModel:
        if d.get("kind", cls.kind) != cls.kind:
            raise ValueError(f"expected a {cls.kind!r} model, got {d.get('kind')!r}")
        poly = d.get("poly")
        return cls(
            d["x_c"],
            d["B"],
            d["alpha"],
            None if poly is None else np.asarray(poly, dtype=float),
            d["n_pix"],
            geometry,
        )


# ----------------------------------------------------------------------------- builders


def model_from_config(kind: str, geom: GratingGeometry, n_pix: int) -> WavelengthModel:
    """The nominal model of a configuration: ``"grating"`` at :func:`nominal_solution`.

    ``"dispersion"`` (DeepSpec) is a later plan and raises ``NotImplementedError``; any
    other kind is a ``ValueError``.
    """
    if kind == "grating":
        x_c, B = nominal_solution(geom, n_pix)
        return GratingModel(x_c, B, geom.alpha_rad, None, n_pix, geom)
    if kind == "dispersion":
        raise NotImplementedError("the dispersion wavelength model belongs to Plan 2")
    raise ValueError(f"unknown wavelength model kind {kind!r}")


def model_from_solution(solution: WavelengthSolution, geom: GratingGeometry) -> WavelengthModel:
    """Rebuild the model a :class:`~mastpipe.core.model.WavelengthSolution` stores.

    The inverse of ``WavelengthSolution(kind=model.kind, params=model.params,
    n_pix=model.n_pix, ...)``.
    """
    if solution.kind == "grating":
        return GratingModel.from_params(solution.params, solution.n_pix, geom)
    if solution.kind == "dispersion":
        raise NotImplementedError("the dispersion wavelength model belongs to Plan 2")
    raise ValueError(f"unknown wavelength model kind {solution.kind!r}")


# ------------------------------------------------------------------------------ fitting


def poly_powers(fit_alpha: bool, poly_deg: int) -> list[int]:
    """Powers of ``t`` the distortion term may use, highest first.

    The constant and linear terms are always excluded (they are ``x_c`` and ``B``), the
    quadratic as well when ``alpha`` is free, so with ``alpha`` free ``P`` starts at the
    cubic. Empty when ``poly_deg`` is below the lowest allowed power.
    """
    lowest = 3 if fit_alpha else 2
    return list(range(poly_deg, lowest - 1, -1)) if poly_deg >= lowest else []


def _poly_from_powers(powers: list[int], values: np.ndarray) -> np.ndarray | None:
    """Coefficients in ``polyval`` order from the fitted values of ``powers``."""
    if not powers:
        return None
    coeffs = np.zeros(max(powers) + 1)
    for power, value in zip(powers, values, strict=True):
        coeffs[max(powers) - power] = value
    return coeffs


def fit_model(
    model: GratingModel,
    x,
    lam,
    weights=None,
    fit_alpha: bool = False,
    poly_deg: int = 0,
    alpha_bound_deg: float | None = None,
    B_bound_frac: float | None = None,
) -> tuple[GratingModel, np.ndarray]:
    """Least squares of the grating model on identified lines; residuals in pixels.

    ``x`` and ``lam`` are the measured columns and catalogue wavelengths of the lines,
    ``weights`` multiplies their residuals (a weight of zero drops a line). ``model``
    supplies the starting ``x_c``, ``B`` and ``alpha``; the distortion term starts at zero
    with the powers of :func:`poly_powers`. ``alpha_bound_deg`` keeps the incidence angle
    within that many degrees of the design value and ``B_bound_frac`` keeps ``B`` within
    that fraction of ``B_nom``: ``alpha`` and ``B`` are nearly the same parameter here
    (correlated at 0.9993 in the laboratory data), so an unconstrained fit trades one
    against the other. Holding ``alpha`` (the default) is usually enough.

    Returns the fitted model and ``x_model(lam) - x`` (unweighted) for every line.
    """
    x = np.asarray(x, dtype=float)
    lam = np.asarray(lam, dtype=float)
    if x.shape != lam.shape:
        raise ValueError(f"x and lam differ in shape: {x.shape} vs {lam.shape}")
    w = np.ones_like(x) if weights is None else np.asarray(weights, dtype=float)
    if w.shape != x.shape:
        raise ValueError(f"weights shape {w.shape} differs from x {x.shape}")
    geom = model.geom
    n_pix = model.n_pix
    powers = poly_powers(fit_alpha, poly_deg)
    sign = np.sign(model.B) or 1.0

    p0 = [model.x_c, model.B]
    lo, hi = [-np.inf, -np.inf], [np.inf, np.inf]
    scale = [1e2, 1e4]
    if B_bound_frac is not None:
        lo[1], hi[1] = sorted(
            [sign * geom.B_nom * (1.0 - B_bound_frac), sign * geom.B_nom * (1.0 + B_bound_frac)]
        )
    if fit_alpha:
        p0.append(model.alpha)
        scale.append(1e-3)
        if alpha_bound_deg is None:
            lo.append(-np.inf)
            hi.append(np.inf)
        else:
            lo.append(geom.alpha_rad - math.radians(alpha_bound_deg))
            hi.append(geom.alpha_rad + math.radians(alpha_bound_deg))
    p0 += [0.0] * len(powers)
    scale += [1.0] * len(powers)
    lo += [-np.inf] * len(powers)
    hi += [np.inf] * len(powers)
    # least_squares needs the start strictly inside the bounds
    p0 = [
        min(max(v, lower + 1e-9), upper - 1e-9) if np.isfinite(lower) and np.isfinite(upper) else v
        for v, lower, upper in zip(p0, lo, hi, strict=True)
    ]
    t = tnorm(x, n_pix)

    def unpack(p) -> tuple[float, float, float, np.ndarray | None]:
        alpha = p[2] if fit_alpha else model.alpha
        return p[0], p[1], alpha, _poly_from_powers(powers, p[2 + int(fit_alpha) :])

    def model_x(p) -> np.ndarray:
        x_c, B, alpha, poly = unpack(p)
        predicted = np.asarray(grating_x(lam, x_c, B, alpha, geom), dtype=float)
        if poly is not None:
            predicted = predicted + np.polyval(poly, t)
        return predicted

    def residuals(p) -> np.ndarray:
        return (model_x(p) - x) * w

    result = least_squares(residuals, p0, bounds=(lo, hi), x_scale=scale)
    x_c, B, alpha, poly = unpack(result.x)
    fitted = replace(model, x_c=float(x_c), B=float(B), alpha=float(alpha), poly=poly)
    return fitted, model_x(result.x) - x


def anchor_model(
    pairs: list[tuple[float, float]], geom: GratingGeometry, n_pix: int
) -> GratingModel:
    """The grating model through ``(wavelength nm, pixel)`` anchor pairs.

    Least squares of ``x = x_c + B tan(theta(lam))`` at the design ``alpha``: two pairs fix
    ``(x_c, B)`` exactly, more are fitted, and one pair fixes ``x_c`` at the nominal
    (signed) ``B``. ``ValueError`` without any pair.
    """
    if len(pairs) == 0:
        raise ValueError("anchor_model needs at least one (wavelength, pixel) pair")
    waves = np.array([float(pair[0]) for pair in pairs])
    pixels = np.array([float(pair[1]) for pair in pairs])
    tan_theta = np.tan(np.asarray(theta_of_lambda(waves, geom), dtype=float).reshape(-1))
    if waves.size == 1:
        B = geom.B_nom * geom.dispersion_sign
        x_c = pixels[0] - B * tan_theta[0]
    else:
        design = np.vstack([np.ones_like(tan_theta), tan_theta]).T
        (x_c, B), *_ = np.linalg.lstsq(design, pixels, rcond=None)
    return GratingModel(float(x_c), float(B), geom.alpha_rad, None, n_pix, geom)
