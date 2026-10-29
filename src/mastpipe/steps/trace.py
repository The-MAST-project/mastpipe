"""Trace finding, layout matching and arc classification: a port of the notebook's
``find_traces`` plus the layout step the design adds on top of it.

:func:`find_traces` follows the notebook. A band of ``±center_width`` columns around the
detector centre is median-collapsed along x into a spatial profile; a running median of
``baseline_window`` rows is subtracted; peaks of the residual above ``threshold_sigma``
times the profile noise, at least ``min_separation`` rows apart, are the candidate traces
(the brightest ``max_traces`` when a limit is set). Each candidate is then followed along
x in blocks of ``x_bin`` median-combined columns: inside ``±trace_halfwidth`` rows of the
predicted position, the pixels above ``centroid_sigma`` times the block noise give a
flux-weighted centroid, one sample per block. A polynomial of ``degree`` is fitted with
``clip_sigma`` clipping in two passes, the second pass re-centroiding about the first
fit; traces keeping fewer than ``min_points`` samples are dropped.

One estimate departs from the notebook. There, the noise of the collapsed profile is its
own robust sigma, which is a noise estimate only while most rows are empty: on a layout
of 45 traces at an 8.5 px pitch the median absolute deviation measures the comb itself
(about a fifth of a peak, so a 5σ threshold exceeds every peak), and on a sparse layout
the median of integer-valued pixels is quantised, which makes the deviation collapse
towards zero. Here the noise of a median-collapsed profile is what it is statistically:
``sqrt(pi/2)`` times the typical robust scatter of a row across the collapsed columns,
divided by the square root of their number. The two agree where the notebook's estimate
is a noise estimate.

:func:`match_layout` finds one global offset between the detected traces and the layout's
nominal rows with :func:`comb_shift`, then assigns layout rows to detections one-to-one by
distance within a tolerance. :func:`classify_arc` decides from the lit traces whether an
arc frame is a full arc or a calibration-fiber arc. Everything here is pure: arrays,
dataclasses, a parameter model and the layout table in, dataclasses and plain values out.
"""

from __future__ import annotations

import warnings
from dataclasses import replace

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.stats import norm

from mastpipe.config.models import TraceParams
from mastpipe.core.model import Frame, Trace, TraceSet

COLLAPSE_STATISTICS = ("median", "mean", "percentile")
"""How a band of columns is collapsed along x into a spatial profile."""
MAD_TO_SIGMA = 1.4826
"""The median absolute deviation of Gaussian noise scaled to its standard deviation."""
MEDIAN_EFFICIENCY = 1.2533
"""Standard error of a median relative to that of a mean, ``sqrt(pi / 2)`` for Gaussian noise."""
CLIP_ITERATIONS = 5
"""Upper bound on the fit-and-clip rounds of one polynomial pass."""
RMS_FLOOR_PX = 1e-3
"""Smallest residual scatter used for clipping, so a perfect fit does not clip everything."""
COMB_KERNEL_STEPS = 2.0
"""Width of the kernel :func:`comb_shift` lays on every position difference, in grid steps."""
LAYOUT_SEARCH_PX = 20.0
"""Largest search range of the global layout shift, in rows."""
LAYOUT_SEARCH_PITCH = 0.45
"""The search range as a fraction of the smallest layout pitch, kept below one half so a
periodic layout cannot be matched one fiber off."""
LIT_BASELINE_HALFWIDTHS = 4
"""Rows of the per-column running median removed before a lit test, in aperture halfwidths."""
SPIKE_FILTER = 3
"""Columns of the running median that removes single-column spikes before a lit test."""


# ---------------------------------------------------------------------------- statistics


def robust_sigma(v: np.ndarray) -> float:
    """``1.4826`` times the median absolute deviation of the finite values of ``v``.

    ``nan`` when ``v`` has no finite value.
    """
    values = np.asarray(v, dtype=float).ravel()
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    return float(MAD_TO_SIGMA * np.median(np.abs(values - np.median(values))))


def _statistic_factor(statistic: str, percentile: float) -> float:
    """Standard error of the collapse statistic relative to that of a mean, for Gaussian
    noise: 1 for the mean, ``sqrt(pi/2)`` for the median, ``sqrt(p (1 - p)) / phi(z_p)``
    for the ``p`` quantile."""
    if statistic == "mean":
        return 1.0
    if statistic == "median":
        return MEDIAN_EFFICIENCY
    if statistic == "percentile":
        if not 0.0 < percentile < 100.0:
            raise ValueError(f"collapse_percentile must be within (0, 100), got {percentile}")
        p = percentile / 100.0
        return float(np.sqrt(p * (1.0 - p)) / norm.pdf(norm.ppf(p)))
    raise ValueError(f"collapse must be one of {COLLAPSE_STATISTICS}, got {statistic!r}")


def _collapse(
    block: np.ndarray, statistic: str = "median", percentile: float = 90.0
) -> tuple[np.ndarray, float]:
    """``block`` collapsed along x by ``statistic`` (the median, the mean or the
    ``percentile`` of each row) and the noise of that statistic.

    The noise is the statistic's standard error for Gaussian noise: the median over rows
    of each row's robust scatter across the block's columns, divided by the square root
    of the column count and scaled by :func:`_statistic_factor`. Rows without a finite
    value collapse to ``0.0``.
    """
    factor = _statistic_factor(statistic, percentile)
    block = np.asarray(block, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN rows are handled below
        centre = np.nanmedian(block, axis=1)
        if statistic == "mean":
            profile = np.nanmean(block, axis=1)
        elif statistic == "percentile":
            profile = np.nanpercentile(block, percentile, axis=1)
        else:
            profile = centre
        deviation = np.nanmedian(np.abs(block - centre[:, None]), axis=1)
    row_sigma = MAD_TO_SIGMA * deviation
    row_sigma = row_sigma[np.isfinite(row_sigma)]
    n_columns = block.shape[1]
    if row_sigma.size == 0 or n_columns == 0:
        noise = 0.0
    else:
        noise = float(factor * np.median(row_sigma) / np.sqrt(n_columns))
    return np.nan_to_num(profile, nan=0.0), noise


def _odd(window: int) -> int:
    window = max(1, int(window))
    return window if window % 2 else window + 1


# --------------------------------------------------------------------------- trace finding


def _centroid(
    profile: np.ndarray, noise: float, y_predicted: float, halfwidth: int, threshold_sigma: float
) -> float | None:
    """Flux-weighted centroid of the pixels above ``threshold_sigma · noise`` within
    ``±halfwidth`` rows of ``y_predicted``; ``None`` when fewer than two pixels qualify
    or the brightest one sits on the window's edge (the trace is not in the window)."""
    ny = profile.shape[0]
    if not np.isfinite(y_predicted):
        return None
    centre = int(np.rint(y_predicted))
    lo, hi = max(0, centre - halfwidth), min(ny, centre + halfwidth + 1)
    if hi - lo < 3:
        return None
    values = profile[lo:hi]
    brightest = int(np.argmax(values))
    if brightest in (0, hi - lo - 1):
        return None
    good = values > threshold_sigma * noise
    if good.sum() < 2:
        return None
    rows = np.arange(lo, hi, dtype=float)
    weights = values[good]
    return float(np.sum(rows[good] * weights) / np.sum(weights))


def _clipped_fit(
    x: np.ndarray, y: np.ndarray, degree: int, clip_sigma: float
) -> tuple[np.ndarray, np.ndarray, float] | None:
    """Polynomial of ``degree`` through ``(x, y)`` with iterated ``clip_sigma`` clipping.

    Returns ``(coeffs, kept, rms)`` with ``coeffs`` in ``polyval`` order, or ``None`` when
    fewer than ``degree + 1`` points remain.
    """
    kept = np.ones(x.shape[0], dtype=bool)
    coeffs = np.zeros(degree + 1)
    residual = np.zeros(x.shape[0])
    for _ in range(CLIP_ITERATIONS):
        if kept.sum() < degree + 1:
            return None
        fitted = np.polynomial.Polynomial.fit(x[kept], y[kept], degree).convert().coef
        coeffs = np.zeros(degree + 1)
        coeffs[: fitted.shape[0]] = fitted
        coeffs = coeffs[::-1]
        residual = y - np.polyval(coeffs, x)
        scale = max(float(np.std(residual[kept])), RMS_FLOOR_PX)
        tighter = kept & (np.abs(residual) <= clip_sigma * scale)
        if np.array_equal(tighter, kept):
            break
        kept = tighter
    if kept.sum() < degree + 1:
        return None
    rms = float(np.sqrt(np.mean(residual[kept] ** 2)))
    return coeffs, kept, rms


def _follow(
    blocks: list[tuple[float, np.ndarray, float]],
    predict,
    halfwidth: int,
    centroid_sigma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """One centroid per block where the trace is measurable, about ``predict(x_block)``."""
    xs, ys = [], []
    for x_block, profile, noise in blocks:
        y = _centroid(profile, noise, float(predict(x_block)), halfwidth, centroid_sigma)
        if y is not None:
            xs.append(x_block)
            ys.append(y)
    return np.array(xs, dtype=float), np.array(ys, dtype=float)


def find_traces(image: np.ndarray, params: TraceParams) -> TraceSet:
    """Find the fiber traces on ``image`` (``(ny, nx)``, dispersion along x).

    The algorithm is described in the module docstring; ``params.collapse`` picks the
    statistic (median, mean or ``params.collapse_percentile`` of each row) that collapses
    the detection band and the ``x_bin`` blocks along x, the percentile being the one
    that sees a trace carrying only emission lines. Each trace carries ``coeffs`` in
    ``polyval`` order, ``y_center`` and ``slope`` at the centre column ``nx // 2``, the
    ``rms`` of its kept centroids about the fit, ``n_points`` (kept centroids), ``snr``
    (detection peak height over the collapsed-profile noise) and the second-pass samples
    ``x_samples``, ``y_samples`` with their ``kept`` flags. Traces are sorted by
    ``y_center`` and indexed ``0 .. n - 1``; the layout fields stay ``None`` until
    :func:`match_layout`. Non-finite pixels take no part.
    """
    image = np.asarray(image, dtype=float)
    if image.ndim != 2:
        raise ValueError(f"find_traces expects a 2-D image, got {image.ndim}-D")
    ny, nx = image.shape
    x_center = nx // 2
    halfwidth = int(params.trace_halfwidth)
    window = _odd(params.baseline_window or 8 * halfwidth + 1)
    min_separation = params.min_separation or 2.0 * halfwidth
    centroid_sigma = params.centroid_sigma
    if centroid_sigma is None:
        centroid_sigma = min(3.0, params.threshold_sigma)

    # detection on the central band
    lo, hi = max(0, x_center - params.center_width), min(nx, x_center + params.center_width + 1)
    collapse = (params.collapse, params.collapse_percentile)
    profile, noise = _collapse(image[:, lo:hi], *collapse)
    residual = profile - median_filter(profile, size=window, mode="nearest")
    peaks, properties = find_peaks(
        residual,
        height=params.threshold_sigma * noise,
        distance=max(1, int(np.rint(min_separation))),
    )
    heights = properties["peak_heights"]
    if params.max_traces is not None and peaks.shape[0] > params.max_traces:
        order = np.argsort(-heights, kind="stable")[: params.max_traces]
        peaks, heights = peaks[order], heights[order]
    noise_floor = max(noise, np.finfo(float).tiny)

    # per-block profiles, shared by every trace
    x_bin = max(1, int(params.x_bin))
    blocks = []
    for start in range(0, nx - x_bin + 1, x_bin):
        block_profile, block_noise = _collapse(image[:, start : start + x_bin], *collapse)
        block_profile = block_profile - median_filter(block_profile, size=window, mode="nearest")
        blocks.append((start + (x_bin - 1) / 2.0, block_profile, block_noise))

    traces: list[Trace] = []
    for peak, height in zip(peaks, heights, strict=True):
        xs, ys = _follow(blocks, lambda _x, y0=float(peak): y0, halfwidth, centroid_sigma)
        first = _clipped_fit(xs, ys, params.degree, params.clip_sigma)
        if first is None:
            continue
        coeffs_1 = first[0]
        xs, ys = _follow(blocks, lambda x, c=coeffs_1: np.polyval(c, x), halfwidth, centroid_sigma)
        second = _clipped_fit(xs, ys, params.degree, params.clip_sigma)
        if second is None:
            continue
        coeffs, kept, rms = second
        n_points = int(kept.sum())
        if n_points < params.min_points:
            continue
        traces.append(
            Trace(
                index=0,
                coeffs=coeffs,
                y_center=float(np.polyval(coeffs, x_center)),
                slope=float(np.polyval(np.polyder(coeffs), x_center)),
                rms=rms,
                n_points=n_points,
                snr=float(height / noise_floor),
                x_samples=xs,
                y_samples=ys,
                kept=kept,
            )
        )
    traces.sort(key=lambda trace: trace.y_center)
    traces = [replace(trace, index=i) for i, trace in enumerate(traces)]
    return TraceSet(traces=traces, shape=(ny, nx))


# -------------------------------------------------------------------------- comb matching


def comb_shift(
    a: np.ndarray, b: np.ndarray, max_shift: float, step: float = 0.1, width: float | None = None
) -> float:
    """The shift ``s`` in ``[-max_shift, max_shift]`` for which ``a + s`` best overlaps ``b``.

    Every pairwise difference ``b_j - a_i`` within range is laid down as a Gaussian of
    ``width`` (default ``COMB_KERNEL_STEPS · step``; callers set it to the scatter they
    expect between matching members, such as their matching tolerance); the sum is
    sampled on a grid of ``step`` and the maximum refined by a parabola through its
    neighbours. Members without a partner add no peak and do not pull the estimate.
    ``0.0`` when either comb is empty or no difference falls within the range.
    """
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if a.size == 0 or b.size == 0 or max_shift <= 0 or step <= 0:
        return 0.0
    if width is None:
        width = COMB_KERNEL_STEPS * step
    if width <= 0:
        raise ValueError(f"comb_shift width must be positive, got {width}")
    differences = (b[None, :] - a[:, None]).ravel()
    differences = differences[np.abs(differences) <= max_shift + 4.0 * width]
    if differences.size == 0:
        return 0.0
    n_grid = int(np.rint(2.0 * max_shift / step)) + 1
    grid = -max_shift + step * np.arange(n_grid)
    score = np.exp(-0.5 * ((differences[:, None] - grid[None, :]) / width) ** 2).sum(axis=0)
    k = int(np.argmax(score))
    shift = float(grid[k])
    if 0 < k < n_grid - 1:
        y0, y1, y2 = score[k - 1 : k + 2]
        curvature = y0 - 2.0 * y1 + y2
        if curvature < 0:
            shift += 0.5 * step * float(y0 - y2) / float(curvature)
    return float(np.clip(shift, -max_shift, max_shift))


def _layout_value(value, cast):
    if value is None or (isinstance(value, float) and np.isnan(value)) or pd.isna(value):
        return None
    return cast(value)


def match_layout(traceset: TraceSet, layout: pd.DataFrame, tolerance: float) -> TraceSet:
    """Assign the layout's rows to the detected traces and report what does not match.

    One global offset ``layout_shift`` (detected minus nominal, in rows) comes from
    :func:`comb_shift` with a kernel of ``tolerance`` over
    ``±min(LAYOUT_SEARCH_PX, LAYOUT_SEARCH_PITCH · smallest layout pitch)``: a periodic
    layout matched one pitch off would still pair almost every row, so the search never
    reaches the alias. Each layout row then takes the detected trace closest to
    ``y_nominal + layout_shift`` within ``tolerance``, one-to-one, the smallest distances
    first. Matched traces receive ``role``, ``telescope``, ``layout_index`` (the row's
    position in ``layout``) and ``aperture`` (``None`` when the layout leaves it unset);
    unmatched ones keep ``None`` in all four. ``missing`` lists the layout indices without
    a detection, ``extra`` the trace indices without a layout row. The input traces are
    not modified.
    """
    detected = np.array([trace.y_center for trace in traceset.traces], dtype=float)
    nominal = layout["y_nominal"].to_numpy(dtype=float)
    shift = 0.0
    if detected.size and nominal.size:
        max_shift = LAYOUT_SEARCH_PX
        if nominal.size > 1:
            pitch = float(np.min(np.diff(np.sort(nominal))))
            max_shift = min(max_shift, LAYOUT_SEARCH_PITCH * pitch)
        if max_shift > 0:
            shift = -comb_shift(detected, nominal, max_shift, width=float(tolerance))
    predicted = nominal + shift

    pairs = sorted(
        (abs(y - predicted[i]), i, j)
        for i in range(nominal.shape[0])
        for j, y in enumerate(detected)
        if abs(y - predicted[i]) <= tolerance
    )
    layout_of: dict[int, int] = {}
    rows_taken: set[int] = set()
    for _, i, j in pairs:
        if i in rows_taken or j in layout_of:
            continue
        layout_of[j] = i
        rows_taken.add(i)

    traces = []
    for j, trace in enumerate(traceset.traces):
        i = layout_of.get(j)
        if i is None:
            unmatched = dict(role=None, telescope=None, layout_index=None, aperture=None)
            traces.append(replace(trace, **unmatched))
            continue
        row = layout.iloc[i]
        traces.append(
            replace(
                trace,
                role=_layout_value(row["role"], str),
                telescope=_layout_value(row["telescope"], str),
                layout_index=int(i),
                aperture=_layout_value(row["aperture"], float),
            )
        )
    missing = tuple(i for i in range(nominal.shape[0]) if i not in rows_taken)
    extra = tuple(j for j in range(detected.shape[0]) if j not in layout_of)
    return TraceSet(
        traces=traces,
        shape=traceset.shape,
        layout_shift=float(shift),
        missing=missing,
        extra=extra,
    )


# ------------------------------------------------------------------------ frames by trace


def trace_source(frames_by_type: dict[str, list[Frame]]) -> tuple[str, Frame | None]:
    """Where to find traces: ``("flat", frame)`` when a flat is listed, else
    ``("arc_full", frame)`` for the first full arc, else ``("none", None)``.

    A full arc is emission lines on a dark continuum, invisible to a median collapse
    along x; when the source is ``"arc_full"`` the runner passes ``find_traces``
    parameters with ``collapse="percentile"``.
    """
    for key in ("flat", "arc_full"):
        frames = frames_by_type.get(key) or []
        if frames:
            return key, frames[0]
    return "none", None


def _aperture_sum(
    data: np.ndarray, trace: Trace, halfwidth: float, mask: np.ndarray | None = None
) -> tuple[np.ndarray, int]:
    """Column sums of ``data`` over the rows within ``ceil(halfwidth)`` of the rounded trace
    position, masked and non-finite pixels counting as zero; with the row count."""
    ny, nx = data.shape
    x = np.arange(nx)
    centre = np.rint(np.asarray(trace.y(x), dtype=float)).astype(int)
    half = int(np.ceil(halfwidth))
    rows = centre[None, :] + np.arange(-half, half + 1)[:, None]
    inside = (rows >= 0) & (rows < ny)
    columns = np.broadcast_to(x, rows.shape)
    values = np.zeros(rows.shape)
    picked = np.asarray(data, dtype=float)[rows[inside], columns[inside]]
    if mask is not None:
        picked = np.where(mask[rows[inside], columns[inside]], 0.0, picked)
    values[inside] = np.nan_to_num(picked, nan=0.0, posinf=0.0, neginf=0.0)
    return values.sum(axis=0), 2 * half + 1


def classify_arc(
    frame: Frame,
    traceset: TraceSet,
    halfwidth: float,
    lit_sigma: float = 5.0,
    lit_fraction: float = 0.1,
) -> tuple[str, dict[int, bool]]:
    """Which traces carry light on ``frame``, and what kind of arc that makes it.

    The frame first loses a per-column running median along y over
    ``LIT_BASELINE_HALFWIDTHS · ceil(halfwidth)`` rows (made odd), which removes a residual
    background and the smooth wings of neighbouring fibers. A trace's level is then the
    mean along x of its aperture-summed spectrum (rows within ``halfwidth`` of the trace,
    single-column spikes removed by a running median of :data:`SPIKE_FILTER` columns); it
    is lit when that level exceeds both ``lit_sigma`` times the noise of the mean (the
    baseline-subtracted frame's robust pixel sigma scaled by the aperture's row count and
    the column count) and ``lit_fraction`` times the brightest trace's level, so that a
    frame's own lit fibers set the scale and faint leakage never counts. The mean rather
    than the median is compared because an arc spectrum is emission lines on a dark
    continuum: its median along x is the level between lines. Returns ``("full", flags)``
    when every matched trace (``role`` set) is lit, ``("calfiber", flags)`` when every
    ``arc`` trace is lit and no ``science`` or ``sky`` trace is, else ``("unknown",
    flags)``; ``flags`` maps trace index to its lit flag, unmatched traces included.
    """
    data = np.asarray(frame.data, dtype=float)
    ny, nx = data.shape
    finite = np.isfinite(data)
    if frame.mask is not None:
        finite &= ~frame.mask
    filled = np.where(finite, data, 0.0)
    window = _odd(LIT_BASELINE_HALFWIDTHS * int(np.ceil(halfwidth)) + 1)
    data = filled - median_filter(filled, size=(window, 1), mode="nearest")
    pixel_sigma = robust_sigma(data[finite])
    if not np.isfinite(pixel_sigma):
        pixel_sigma = 0.0
    levels: dict[int, float] = {}
    noises: dict[int, float] = {}
    for trace in traceset:
        spectrum, n_rows = _aperture_sum(data, trace, halfwidth, frame.mask)
        levels[trace.index] = float(
            np.mean(median_filter(spectrum, size=SPIKE_FILTER, mode="nearest"))
        )
        noises[trace.index] = pixel_sigma * np.sqrt(n_rows) / np.sqrt(nx)
    brightest = max(levels.values(), default=0.0)
    flags = {
        index: bool(level > max(lit_sigma * noises[index], lit_fraction * brightest))
        for index, level in levels.items()
    }

    matched = [trace for trace in traceset if trace.role is not None]
    arcs = [trace for trace in matched if trace.role == "arc"]
    targets = [trace for trace in matched if trace.role in ("science", "sky")]
    if matched and all(flags[trace.index] for trace in matched):
        kind = "full"
    elif (
        arcs
        and all(flags[trace.index] for trace in arcs)
        and not any(flags[trace.index] for trace in targets)
    ):
        kind = "calfiber"
    else:
        kind = "unknown"
    return kind, flags


def relative_throughput(
    flat: Frame | None, traceset: TraceSet, halfwidth: float
) -> dict[int, float]:
    """Integrated flat flux of every trace (aperture sum over all columns) divided by the
    median over traces, keyed by trace index; ``{}`` without a flat or without traces.
    ``ValueError`` when the median flux is not positive."""
    if flat is None or len(traceset) == 0:
        return {}
    totals = {
        trace.index: float(_aperture_sum(flat.data, trace, halfwidth, flat.mask)[0].sum())
        for trace in traceset
    }
    scale = float(np.median(list(totals.values())))
    if not scale > 0:
        raise ValueError("the flat carries no flux on the traces; cannot normalise throughput")
    return {index: total / scale for index, total in totals.items()}
