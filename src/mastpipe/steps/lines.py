"""Line finding on the 2-D frame: a port of the notebook's ``find_lines``.

Detection runs on the calibrated image following each trace, not on an extracted
spectrum. For every trace a strip of ``2 halfwidth + 1`` rows is cut and straightened
onto the fitted trace (:func:`_detilt`), a wide running median removes the continuum
from each row (:func:`_remove_row_continuum`), and the emission lines are detected on
the y-collapsed profile, where nothing blends across the dispersion. Peaks closer than
``blend_distance`` form one group that is fitted together as a sum of Gaussians
(:func:`_fit_group`); the residual-driven deblender (:func:`_deblend_group`) adds a
component wherever what is left after the fit still peaks ``add_sigma`` noise units
high. Each line is then measured in both axes inside a window clipped at the midpoint to
its neighbours, and the flags decide the two verdicts: ``use_x`` needs the line resolved
along the dispersion, ``use_y`` only needs it complete and on the trace, so a blend keeps
a valid spatial width.

Beyond the notebook, a ``saturated`` flag (raw strip peak above ``LineParams.saturation``)
bars a line from ``use_x``, and a variance image, when given, provides the per-line
``snr`` in place of the profile's robust noise. The notebook's figures and printing are
not part of the step; the detilted strips come back for the caller's QA.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.optimize import OptimizeWarning, curve_fit
from scipy.signal import find_peaks

from mastpipe.config.models import LineParams
from mastpipe.core.model import Trace, TraceSet
from mastpipe.steps.trace import robust_sigma

GAUSSIAN_FWHM = 2.3548
"""FWHM of a Gaussian in units of its sigma, as the notebook rounds it."""

LINE_COLUMNS = (
    "id",
    "group",
    "n_in_group",
    "x_peak",
    "x_centroid",
    "x_fit",
    "x_center",
    "y_center",
    "y_offset",
    "height",
    "flux",
    "baseline",
    "snr",
    "fwhm_x",
    "fwhm_x_fit",
    "sigma_x_fit",
    "amp_x_fit",
    "fwhm_y",
    "sigma_y_fit",
    "x_window",
    "neighbour",
    "isolated",
    "near_edge",
    "complete",
    "on_trace",
    "fit_ok",
    "saturated",
    "use_x",
    "use_y",
)
"""Columns of every line table, in order."""

_FLAG_COLUMNS = (
    "isolated",
    "near_edge",
    "complete",
    "on_trace",
    "fit_ok",
    "saturated",
    "use_x",
    "use_y",
)
_INT_COLUMNS = ("id", "group", "n_in_group", "x_peak", "x_window")


# ------------------------------------------------------------------------------- helpers


def _detilt(data: np.ndarray, trace: Trace, halfwidth: int) -> tuple[np.ndarray, np.ndarray]:
    """Cut a strip along ``trace`` and straighten it.

    Row ``j`` of the result holds detector row ``round(y_trace(x)) + j - halfwidth``, so
    the trace sits on the middle row in every column. Rows off the detector come back
    NaN. Returns ``(strip, row_centre)`` with ``row_centre`` the rounded detector row of
    the trace in each column.
    """
    ny, nx = data.shape
    x = np.arange(nx)
    centre = np.round(trace.y(x)).astype(int)
    offsets = np.arange(-halfwidth, halfwidth + 1)
    rows = centre[None, :] + offsets[:, None]
    ok = (rows >= 0) & (rows < ny)
    strip = np.where(ok, data[np.clip(rows, 0, ny - 1), x[None, :]], np.nan)
    return strip.astype(float), centre


def _remove_row_continuum(strip: np.ndarray, window: int) -> np.ndarray:
    """Subtract a wide running median from every row of ``strip``.

    The window is far wider than a line, so the median rides the smooth background and
    steps over the lines. NaNs are filled with the row median before filtering and
    restored afterwards.
    """
    window = int(window) | 1
    out = np.empty_like(strip)
    for j, row in enumerate(strip):
        bad = ~np.isfinite(row)
        if bad.all():
            out[j] = np.nan
            continue
        filled = np.where(bad, np.nanmedian(row), row)
        out[j] = row - median_filter(filled, size=window, mode="nearest")
    return out


def _gauss_sum(x: np.ndarray, *p: float) -> np.ndarray:
    """Constant baseline ``p[0]`` plus triplets of (amplitude, centre, sigma)."""
    y = np.full_like(x, p[0], dtype=float)
    for k in range(1, len(p), 3):
        y = y + p[k] * np.exp(-0.5 * ((x - p[k + 1]) / p[k + 2]) ** 2)
    return y


def _fit_group(
    x: np.ndarray,
    y: np.ndarray,
    centres: list[float],
    sigma_guess: float,
    max_shift: float = 3.0,
    sigma_bounds: tuple[float, float] = (0.4, 12.0),
) -> tuple[float, list[tuple[float, float, float]]] | None:
    """Fit ``len(centres)`` Gaussians on a shared baseline.

    Returns ``(baseline, [(amplitude, centre, sigma), ...])`` or ``None`` if the fit
    fails (too few finite samples, no convergence, or a covariance the optimiser cannot
    estimate). Components come back sorted by centre.
    """
    p0, lo, hi = [0.0], [-np.inf], [np.inf]
    for c in centres:
        p0 += [max(float(np.interp(c, x, y)), 1.0), float(c), float(sigma_guess)]
        lo += [0.0, c - max_shift, sigma_bounds[0]]
        hi += [np.inf, c + max_shift, sigma_bounds[1]]
    good = np.isfinite(y)
    if good.sum() < len(p0) + 2:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            popt, _ = curve_fit(_gauss_sum, x[good], y[good], p0=p0, bounds=(lo, hi), maxfev=20000)
    except (RuntimeError, ValueError, OptimizeWarning):
        return None
    comps = sorted(
        (tuple(float(v) for v in popt[k : k + 3]) for k in range(1, len(popt), 3)),
        key=lambda c: c[1],
    )
    return float(popt[0]), comps


def _deblend_group(
    x: np.ndarray,
    y: np.ndarray,
    centres: list[float],
    sigma_guess: float,
    noise: float,
    max_components: int = 4,
    add_sigma: float = 5.0,
    min_split: float | None = None,
    gain: float = 0.15,
) -> tuple[float, list[tuple[float, float, float]]] | None:
    """Fit a group, then keep adding components while the residual says one is missing.

    ``find_peaks`` cannot separate two lines closer than about one FWHM: they show up as
    a single maximum with a shoulder. What recovers them is the residual: fit one
    Gaussian, and if what is left over still has a peak worth ``add_sigma`` noise, put a
    component there and refit. A new component is kept only when it is that significant,
    sits at least ``min_split`` from every component already in the fit, stays inside the
    window and cuts the residual scatter by at least ``gain``. Returns
    ``(baseline, components)`` or ``None``.
    """
    if min_split is None:
        min_split = 1.0 * sigma_guess
    best = _fit_group(x, y, centres, sigma_guess)
    if best is None:
        return None
    good = np.isfinite(y)

    def rms(fit: tuple[float, list[tuple[float, float, float]]]) -> float:
        model = _gauss_sum(x, fit[0], *[v for c in fit[1] for v in c])
        return float(np.sqrt(np.mean((y[good] - model[good]) ** 2)))

    while len(best[1]) < max_components:
        model = _gauss_sum(x, best[0], *[v for c in best[1] for v in c])
        resid = np.where(good, y - model, -np.inf)
        i = int(np.argmax(resid))
        if resid[i] < add_sigma * noise:
            break
        trial = _fit_group(x, y, [c[1] for c in best[1]] + [x[i]], sigma_guess)
        if trial is None:
            break
        centres_new = [c[1] for c in trial[1]]
        too_close = any(
            b - a < min_split for a, b in zip(centres_new, centres_new[1:], strict=False)
        )
        outside = any(c < x[0] or c > x[-1] for c in centres_new)
        weakest = min(c[0] for c in trial[1])
        if (
            too_close
            or outside
            or weakest < add_sigma * noise
            or rms(trial) > (1 - gain) * rms(best)
        ):
            break
        best = trial
    return best


def _numeric_fwhm(y: np.ndarray, i: int, floor: float = 0.0) -> float:
    """FWHM by linear interpolation of the half-power crossings either side of the peak
    at index ``i``. NaN if either side does not come back down in the array, or when the
    peak does not rise above ``floor`` (nothing to halve)."""
    n = y.size
    if not 0 <= i < n or not y[i] > floor:
        return np.nan
    half = floor + 0.5 * (y[i] - floor)

    left = i
    while left > 0 and y[left] > half:
        left -= 1
    if y[left] > half:
        return np.nan
    x_left = np.interp(half, [y[left], y[left + 1]], [left, left + 1])

    right = i
    while right < n - 1 and y[right] > half:
        right += 1
    if y[right] > half:
        return np.nan
    x_right = np.interp(half, [y[right], y[right - 1]], [right, right - 1])
    return float(x_right - x_left)


def _profile_measures(prof: np.ndarray, i_peak: int) -> tuple[float, float]:
    """``(centroid, numeric fwhm)`` of a 1-D profile around index ``i_peak``."""
    w = np.clip(np.where(np.isfinite(prof), prof, 0.0), 0.0, None)
    total = w.sum()
    cen = float((np.arange(w.size) * w).sum() / total) if total > 0 else np.nan
    return cen, _numeric_fwhm(np.nan_to_num(prof), i_peak)


def _group_peaks(peaks: np.ndarray, blend_distance: float) -> list[list[int]]:
    """Runs of peaks whose consecutive spacing is at most ``blend_distance``."""
    groups: list[list[int]] = []
    current: list[int] = []
    for i in peaks:
        if current and i - current[-1] > blend_distance:
            groups.append(current)
            current = []
        current.append(int(i))
    if current:
        groups.append(current)
    return groups


def _empty_table() -> pd.DataFrame:
    table = pd.DataFrame({column: pd.Series(dtype=float) for column in LINE_COLUMNS})
    for column in _INT_COLUMNS:
        table[column] = table[column].astype(int)
    for column in _FLAG_COLUMNS:
        table[column] = table[column].astype(bool)
    return table


# ------------------------------------------------------------------------------ find_lines


def find_lines(
    image: np.ndarray,
    variance: np.ndarray | None,
    traceset: TraceSet,
    params: LineParams,
) -> tuple[list[pd.DataFrame], list[np.ndarray]]:
    """Detect the emission lines on every trace of ``image`` and measure each in both axes.

    ``image`` is a calibrated frame ``(ny, nx)`` with the dispersion along x; ``variance``
    (same shape) is optional. Every trace of ``traceset`` is measured, whatever its role.
    Returns ``(tables, strips)``: one :data:`LINE_COLUMNS` table and one straightened,
    continuum-subtracted strip of ``2 halfwidth + 1`` rows per trace, in trace order.

    Detection: the strip is collapsed along y, the noise is the robust sigma of that
    profile, and ``find_peaks`` takes everything ``threshold_sigma`` noise units high at
    least ``min_separation`` px apart (``max_lines`` keeps the strongest). Peaks closer
    than ``blend_distance`` (default ``5 sigma_guess``) are one group, fitted together
    inside a window bounded by the midpoints to the neighbouring groups and deblended by
    the residual up to ``max_components`` Gaussians.

    Measurement: each line's windows are clipped at the midpoint to the nearest line on
    either side; ``x_center`` is the fitted centre where the fit converged, else the
    centroid. The flags: ``isolated`` (nearest neighbour at least ``isolation`` px away,
    default ``2 x_halfwidth``), ``near_edge`` (within ``edge_margin`` px of either end),
    ``complete`` (the window holds no off-detector pixel and is not cut by the detector
    edge), ``on_trace`` (spatial centroid within ``max_y_offset`` of the trace; ``None``
    disables), ``fit_ok``, ``saturated`` (peak of the raw strip above ``saturation``;
    ``None`` disables). ``use_x`` needs ``isolated or fit_ok``, not ``near_edge``,
    ``complete``, ``on_trace`` and not ``saturated``; ``use_y`` needs not ``near_edge``,
    ``complete``, a finite ``fwhm_y`` and ``on_trace``, and stays valid for a blend.

    ``snr`` is ``height`` over the square root of the variance summed down the strip's
    column when ``variance`` is given, else over the profile's robust noise.
    """
    image = np.asarray(image, dtype=float)
    if variance is not None:
        variance = np.asarray(variance, dtype=float)
        if variance.shape != image.shape:
            raise ValueError(f"variance shape {variance.shape} differs from image {image.shape}")
    halfwidth = int(params.halfwidth)
    x_halfwidth = int(params.x_halfwidth)
    y_halfwidth = halfwidth if params.y_halfwidth is None else int(params.y_halfwidth)
    sigma_guess = float(params.sigma_guess)
    blend_distance = (
        5.0 * sigma_guess if params.blend_distance is None else float(params.blend_distance)
    )
    isolation = 2 * x_halfwidth if params.isolation is None else float(params.isolation)

    tables: list[pd.DataFrame] = []
    strips: list[np.ndarray] = []
    for trace in traceset:
        strip_raw, centre = _detilt(image, trace, halfwidth)
        strip = _remove_row_continuum(strip_raw, params.baseline_window)
        var_prof = None
        if variance is not None:
            var_prof = np.nansum(_detilt(variance, trace, halfwidth)[0], axis=0)
        nrow, nx = strip.shape
        mid = halfwidth

        # detect on the y-collapsed profile: y never blends, so this is a pure gain in
        # signal-to-noise
        prof = np.nansum(strip, axis=0)
        prof[np.all(~np.isfinite(strip), axis=0)] = np.nan
        finite = np.isfinite(prof)
        noise = robust_sigma(prof)
        prof0 = np.where(finite, prof, 0.0)

        if np.isfinite(noise) and noise > 0:
            peaks, _ = find_peaks(
                prof0, height=params.threshold_sigma * noise, distance=params.min_separation
            )
            peaks = peaks[finite[peaks]]
        else:
            peaks = np.zeros(0, dtype=int)
        if params.max_lines is not None and peaks.size > params.max_lines:
            peaks = np.sort(peaks[np.argsort(prof0[peaks])[::-1][: params.max_lines]])

        groups = _group_peaks(peaks, blend_distance)

        # pass 1: deblend each group. The fitting window is bounded by the midpoint to the
        # neighbouring group: padding blindly by x_halfwidth reaches over the next line,
        # the deblender then "discovers" it and the same line ends up in two groups.
        fitted = []
        for g_id, group in enumerate(groups):
            prev_p = groups[g_id - 1][-1] if g_id > 0 else None
            next_p = groups[g_id + 1][0] if g_id + 1 < len(groups) else None
            lo = max(0, group[0] - x_halfwidth - 2)
            hi = min(nx, group[-1] + x_halfwidth + 3)
            if prev_p is not None:
                lo = max(lo, (prev_p + group[0]) // 2 + 1)
            if next_p is not None:
                hi = min(hi, (group[-1] + next_p) // 2 + 1)
            xs = np.arange(lo, hi, dtype=float)

            fit = _deblend_group(
                xs,
                prof0[lo:hi],
                group,
                sigma_guess,
                noise,
                max_components=params.max_components,
                add_sigma=params.add_sigma,
            )
            comps = fit[1] if fit is not None else []
            if len(comps) > len(group):  # the deblender split a shoulder
                group = [int(round(c[1])) for c in comps]
            for k, i in enumerate(group):
                fitted.append(
                    (
                        g_id,
                        len(group),
                        int(i),
                        comps[k] if k < len(comps) else None,
                        fit[0] if fit is not None else np.nan,
                    )
                )

        # pass 2: measure, with every window clipped at the midpoint to the nearest line
        # on each side. A plain +/- x_halfwidth window reaches into a close neighbour, and
        # one bright pixel of it drags the centroid by pixels, which is why the centroid
        # and the fit must be bounded alike.
        centres = [(f[3][1] if f[3] is not None else float(f[2])) for f in fitted]
        rows = []
        for m, (g_id, n_grp, i, comp, base) in enumerate(fitted):
            here = centres[m]
            xlo, xhi = max(0, i - x_halfwidth), min(nx, i + x_halfwidth + 1)
            if m > 0:
                xlo = max(xlo, int(np.ceil((centres[m - 1] + here) / 2.0)))
            if m + 1 < len(fitted):
                xhi = min(xhi, int(np.floor((here + centres[m + 1]) / 2.0)) + 1)
            if xhi - xlo < 3 or not xlo <= i < xhi:  # never leave less than 3 samples
                xlo, xhi = max(0, i - 1), min(nx, i + 2)
            ylo, yhi = max(0, mid - y_halfwidth), min(nrow, mid + y_halfwidth + 1)

            y_prof = np.nansum(strip[:, xlo:xhi], axis=1)
            x_prof = np.nansum(strip[ylo:yhi, :], axis=0)
            y_cen, y_fwhm = _profile_measures(y_prof, int(np.argmax(y_prof)))
            x_cen, x_fwhm = _profile_measures(x_prof[xlo:xhi], i - xlo)
            x_cen += xlo

            yfit = _fit_group(
                np.arange(nrow, dtype=float), y_prof, [y_cen], sigma_guess, max_shift=halfwidth / 2
            )
            sigma_y = yfit[1][0][2] if yfit else np.nan

            if var_prof is not None and var_prof[i] > 0:
                snr = float(prof[i] / np.sqrt(var_prof[i]))
            else:
                snr = float(prof[i] / noise)
            raw_window = np.nan_to_num(strip_raw[ylo:yhi, xlo:xhi], nan=-np.inf)
            raw_peak = float(raw_window.max()) if raw_window.size else -np.inf
            saturated = params.saturation is not None and raw_peak > params.saturation

            # truncation by the detector edge is a defect; clipping against a neighbour is
            # deliberate, so only the former marks a line incomplete
            edge_cut = (i - x_halfwidth < 0) or (i + x_halfwidth + 1 > nx)
            rows.append(
                {
                    "group": g_id,
                    "n_in_group": n_grp,
                    "x_peak": int(i),
                    "x_centroid": x_cen,
                    "x_fit": comp[1] if comp else np.nan,
                    "amp_x_fit": comp[0] if comp else np.nan,
                    "y_center": float(centre[i] + y_cen - mid),
                    "y_offset": float(centre[i] + y_cen - mid - trace.y(i)),
                    "height": float(prof[i]),
                    "snr": snr,
                    "flux": float(np.nansum(strip[ylo:yhi, xlo:xhi])),
                    "x_window": int(xhi - xlo),
                    "fwhm_x": x_fwhm if n_grp == 1 else np.nan,
                    "baseline": base,
                    "sigma_x_fit": comp[2] if comp else np.nan,
                    "fwhm_x_fit": GAUSSIAN_FWHM * comp[2] if comp else np.nan,
                    "fwhm_y": y_fwhm,
                    "sigma_y_fit": sigma_y,
                    "complete": bool(np.all(np.isfinite(strip[:, xlo:xhi]))) and not edge_cut,
                    "fit_ok": comp is not None,
                    "saturated": bool(saturated),
                }
            )

        if not rows:
            tables.append(_empty_table())
            strips.append(strip)
            continue

        tab = pd.DataFrame(rows)
        # the fit models the line and is bounded to its own cell, so it is the better
        # centre wherever it converged; the centroid is the fallback
        tab["x_center"] = np.where(tab["fit_ok"], tab["x_fit"], tab["x_centroid"])
        xc = tab["x_center"].to_numpy(dtype=float)
        if xc.size > 1:
            d = np.abs(xc[:, None] - xc[None, :])
            np.fill_diagonal(d, np.inf)
            nearest = d.min(axis=1)
        else:
            nearest = np.full(xc.size, np.inf)
        tab["neighbour"] = nearest
        tab["isolated"] = nearest >= isolation
        tab["near_edge"] = (xc < params.edge_margin) | (xc > nx - 1 - params.edge_margin)
        # Detection runs on the y-collapsed profile, which is blind to WHERE in the strip
        # the flux sat. A cosmic ray or a hot pixel off the trace collapses into the
        # profile exactly like a line does, so the spatial centroid has to be checked
        # against the trace afterwards.
        if params.max_y_offset is not None:
            tab["on_trace"] = tab["y_offset"].abs() <= params.max_y_offset
        else:
            tab["on_trace"] = np.ones(len(tab), dtype=bool)
        # x needs the line resolved; y only needs it complete and inboard
        tab["use_x"] = (
            (tab["isolated"] | tab["fit_ok"])
            & ~tab["near_edge"]
            & tab["complete"]
            & tab["on_trace"]
            & ~tab["saturated"]
        )
        tab["use_y"] = ~tab["near_edge"] & tab["complete"] & tab["fwhm_y"].notna() & tab["on_trace"]
        tab.insert(0, "id", np.arange(len(tab)))
        tables.append(tab[list(LINE_COLUMNS)].reset_index(drop=True))
        strips.append(strip)

    return tables, strips


def lines_summary(tables: list[pd.DataFrame]) -> pd.DataFrame:
    """One row per trace: counts, and the median width in each axis."""

    def median_of(table: pd.DataFrame, column: str, flag: str) -> float:
        if len(table) and table[flag].any():
            return round(float(table.loc[table[flag], column].median()), 3)
        return np.nan

    out = []
    for i, t in enumerate(tables):
        out.append(
            {
                "trace": i,
                "found": len(t),
                "blended": int((t["n_in_group"] > 1).sum()) if len(t) else 0,
                "use_x": int(t["use_x"].sum()) if len(t) else 0,
                "use_y": int(t["use_y"].sum()) if len(t) else 0,
                "median FWHM_x [px]": median_of(t, "fwhm_x", "use_x"),
                "median FWHM_x fit [px]": median_of(t, "fwhm_x_fit", "use_x"),
                "median FWHM_y [px]": median_of(t, "fwhm_y", "use_y"),
                "median |y_offset| [px]": (
                    round(float(t.loc[t["use_y"], "y_offset"].abs().median()), 3)
                    if len(t) and t["use_y"].any()
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(out)
