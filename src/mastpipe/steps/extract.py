"""Aperture extraction and continuum fitting: a port of the notebook's ``extract_lsf`` and
``fit_continuum``.

:func:`aperture_weights` turns a (possibly fractional) trace centre and a halfwidth into
the detector rows of the aperture and their weights, in one of three conventions:

- ``"round"``: whole rows ``round(y - hw) .. round(y + hw)``, weight 1 each;
- ``"floor"``: whole rows ``floor(y - hw) .. floor(y + hw)``, weight 1 each;
- ``"subpixel"``: the interval of width ``2 hw + 1`` centred on ``y``, with row ``r``
  spanning ``[r - 0.5, r + 0.5)``; the two edge rows take the fraction of their height
  inside the interval, so the weights always sum to ``2 hw + 1``. A centre of 10.3 with
  halfwidth 2 covers ``[7.8, 12.8)``: rows 8..13 weighted ``0.7, 1, 1, 1, 1, 0.3``.

:func:`extract_aperture` follows every trace's ``y(x)`` across the columns of a calibrated
frame and sums (or averages) the weighted pixels. :func:`fit_continuum` fits a Chebyshev
polynomial to a spectrum with asymmetric sigma clipping and stores it as
``Spectrum.continuum``; the raw ``flux`` is never replaced, and callers subtract the two
when they need the continuum-subtracted flux. Everything here is pure: arrays, dataclasses
and a parameter model in, dataclasses out.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from mastpipe.config.models import ExtractParams
from mastpipe.core.model import Frame, Spectrum, TraceSet
from mastpipe.steps.trace import robust_sigma

APERTURES = ("round", "subpixel", "floor")
"""How :func:`aperture_weights` turns a fractional centre into rows and weights."""
METHODS = ("sum", "mean")
"""Whether :func:`extract_aperture` sums or averages the aperture."""


# ------------------------------------------------------------------------ aperture weights


def _aperture_grid(centres: np.ndarray, halfwidth: float, aperture: str) -> tuple[np.ndarray, np.ndarray]:
    """Rows ``(k, n)`` and weights ``(k, n)`` of the aperture around every centre in
    ``centres`` (``(n,)``), padded with zero-weight rows so every column has ``k`` entries.

    The conventions are those of the module docstring; ``ValueError`` for an unknown
    ``aperture`` or a negative ``halfwidth``.
    """
    if aperture not in APERTURES:
        raise ValueError(f"aperture must be one of {APERTURES}, got {aperture!r}")
    if not halfwidth >= 0:
        raise ValueError(f"halfwidth must be non-negative, got {halfwidth}")
    centres = np.asarray(centres, dtype=float)
    if aperture == "subpixel":
        lo = centres - halfwidth - 0.5
        hi = centres + halfwidth + 0.5
        first = np.floor(lo + 0.5).astype(int)  # the first row whose span reaches past lo
        k = int(np.ceil(2.0 * halfwidth + 1.0)) + 1
        rows = first[None, :] + np.arange(k)[:, None]
        overlap = np.minimum(hi[None, :], rows + 0.5) - np.maximum(lo[None, :], rows - 0.5)
        weights = np.clip(overlap, 0.0, 1.0)
        return rows, weights
    edge = np.rint if aperture == "round" else np.floor
    first = edge(centres - halfwidth).astype(int)
    last = edge(centres + halfwidth).astype(int)
    k = int(np.max(last - first)) + 1 if centres.size else 0
    rows = first[None, :] + np.arange(k)[:, None]
    weights = (rows <= last[None, :]).astype(float)
    return rows, weights


def aperture_weights(y_center: float, halfwidth: float, aperture: str = "round") -> tuple[np.ndarray, np.ndarray]:
    """Rows and per-row weights of the aperture of ``halfwidth`` around ``y_center``.

    ``aperture`` is one of :data:`APERTURES` (see the module docstring for the
    conventions). Rows of zero weight are left out, so the returned arrays are as long as
    the aperture is wide; rows may fall outside a detector, which is the caller's concern.
    """
    rows, weights = _aperture_grid(np.array([float(y_center)]), halfwidth, aperture)
    keep = weights[:, 0] > 0
    return rows[keep, 0], weights[keep, 0]


# ---------------------------------------------------------------------------- extraction


def extract_aperture(frame: Frame, traceset: TraceSet, params: ExtractParams) -> list[Spectrum]:
    """Extract one :class:`Spectrum` per trace of ``traceset`` from the calibrated ``frame``.

    Per trace the halfwidth is ``trace.aperture`` when set, else ``params.halfwidth``;
    the rows and weights of each column come from :func:`aperture_weights` at the trace's
    ``y(x)`` with ``params.aperture``. Pixels that are on the detector, finite (data and
    variance) and unmasked are *good*; the column's coverage is the weight of its good
    pixels over the aperture's total weight (rows past the detector edge count as
    uncovered, the aperture never wraps). With ``params.method == "sum"`` the flux is
    ``sum(w · data)`` over the good pixels and the error ``sqrt(sum(w² · variance))``;
    with ``"mean"`` both are divided by the weight of the good pixels (the error of the
    weighted mean). A column whose coverage is below ``params.min_coverage`` has NaN flux
    and error; the mask marks those columns and every column whose aperture holds a pixel
    of non-zero weight that ``frame.mask`` marks (a column with acceptable coverage keeps
    its partial sum, which is not rescaled). ``frame.variance`` is required
    (``ValueError`` when ``None``).

    The spectrum carries ``fiber = trace.index``, the trace's ``role`` and ``telescope``,
    ``target = frame.meta.target_for(trace.telescope)`` for traces with role ``"science"``
    (``None`` otherwise)
    and ``pixel = arange(nx)``.
    """
    if frame.variance is None:
        raise ValueError("extract_aperture needs frame.variance; calibrate the frame first")
    if params.method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {params.method!r}")
    data = np.asarray(frame.data, dtype=float)
    variance = np.asarray(frame.variance, dtype=float)
    ny, nx = data.shape
    x = np.arange(nx)
    spectra: list[Spectrum] = []
    for trace in traceset:
        halfwidth = params.halfwidth if trace.aperture is None else float(trace.aperture)
        centres = np.asarray(trace.y(x), dtype=float)
        rows, weights = _aperture_grid(centres, halfwidth, params.aperture)
        inside = (rows >= 0) & (rows < ny) & (weights > 0)
        columns = np.broadcast_to(x, rows.shape)
        safe_rows = np.where(inside, rows, 0)
        picked = data[safe_rows, columns]
        picked_var = variance[safe_rows, columns]
        good = inside & np.isfinite(picked) & np.isfinite(picked_var)
        flagged = np.zeros(nx, dtype=bool)
        if frame.mask is not None:
            masked = inside & frame.mask[safe_rows, columns]
            good &= ~masked
            flagged = masked.any(axis=0)
        w_good = np.where(good, weights, 0.0)
        total = weights.sum(axis=0)
        covered = w_good.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            coverage = np.where(total > 0, covered / total, 0.0)
        flux = (w_good * np.where(good, picked, 0.0)).sum(axis=0)
        error = np.sqrt((w_good**2 * np.where(good, picked_var, 0.0)).sum(axis=0))
        if params.method == "mean":
            with np.errstate(invalid="ignore", divide="ignore"):
                flux = np.where(covered > 0, flux / covered, np.nan)
                error = np.where(covered > 0, error / covered, np.nan)
        low = coverage < params.min_coverage
        flux[low] = np.nan
        error[low] = np.nan
        spectra.append(
            Spectrum(
                fiber=trace.index,
                pixel=x.copy(),
                flux=flux,
                error=error,
                mask=low | flagged,
                role=trace.role,
                telescope=trace.telescope,
                target=frame.meta.target_for(trace.telescope) if trace.role == "science" else None,
            )
        )
    return spectra


# ------------------------------------------------------------------------------ continuum


def _chebyshev(pixel: np.ndarray, flux: np.ndarray, kept: np.ndarray, degree: int) -> np.ndarray:
    """Chebyshev polynomial of ``degree`` through the kept points, on ``pixel`` mapped to
    ``[-1, 1]`` over the whole grid, evaluated at every pixel."""
    domain = [float(pixel[0]), float(pixel[-1])]
    fit = np.polynomial.Chebyshev.fit(pixel[kept], flux[kept], degree, domain=domain)
    return np.asarray(fit(pixel), dtype=float)


def fit_continuum(spectrum: Spectrum, degree: int, clip_above: float, clip_below: float, n_iter: int = 10) -> Spectrum:
    """``spectrum`` with its ``continuum`` set to a clipped Chebyshev fit of its flux.

    The fit is a Chebyshev polynomial of ``degree`` on the pixel grid mapped to
    ``[-1, 1]``, through the usable pixels (unmasked, finite flux). It is iterated at most
    ``n_iter`` times: after each fit the residuals ``flux - fit`` of the kept pixels give
    their median and robust sigma, and the next kept set is every usable pixel whose
    residual is at most ``clip_above · sigma`` above that median (emission lines go) and
    at most ``clip_below · sigma`` below it (cold pixels go); the iteration stops when the
    kept set no longer changes, when the sigma is zero, or when fewer than ``degree + 1``
    pixels would remain. Clipping about the median rather than about zero keeps the
    continuum pixels of a line-dominated spectrum, which the first, line-pulled fit leaves
    far below itself. The continuum is evaluated at every pixel, masked ones included.

    The returned spectrum is a new object with ``continuum`` set; ``flux``, ``error`` and
    ``mask`` are the raw extraction, left exactly as they were. The continuum-subtracted
    flux is ``flux - continuum``, computed by the caller. When fewer than ``degree + 1``
    pixels are usable the continuum is NaN everywhere.
    """
    pixel = np.asarray(spectrum.pixel, dtype=float)
    flux = np.asarray(spectrum.flux, dtype=float)
    usable = np.isfinite(flux) & ~spectrum.mask
    if usable.sum() < degree + 1 or pixel.size < 2:
        return replace(spectrum, continuum=np.full(pixel.shape, np.nan))
    kept = usable
    for _ in range(max(1, int(n_iter))):
        model = _chebyshev(pixel, flux, kept, degree)
        fitted_on = kept
        residual = flux - model
        sigma = robust_sigma(residual[kept])
        if not sigma > 0:
            break
        centre = float(np.median(residual[kept]))
        with np.errstate(invalid="ignore"):
            tighter = usable & (residual - centre <= clip_above * sigma) & (residual - centre >= -clip_below * sigma)
        if tighter.sum() < degree + 1 or np.array_equal(tighter, kept):
            break
        kept = tighter
    if not np.array_equal(fitted_on, kept):
        model = _chebyshev(pixel, flux, kept, degree)
    return replace(spectrum, continuum=model)


def continuum_subtract_all(spectra: list[Spectrum], params: ExtractParams) -> list[Spectrum]:
    """:func:`fit_continuum` applied to every spectrum with ``params.continuum_degree``,
    ``clip_above``, ``clip_below`` and ``continuum_iter``.

    Despite the name this sets each spectrum's ``continuum`` and leaves ``flux`` as the raw
    extraction; the subtraction ``flux - continuum`` is the caller's.
    """
    return [
        fit_continuum(
            spectrum,
            params.continuum_degree,
            params.clip_above,
            params.clip_below,
            params.continuum_iter,
        )
        for spectrum in spectra
    ]
