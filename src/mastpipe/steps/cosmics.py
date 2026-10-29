"""Cosmic-ray mask: single-pixel hits standing far above a median-filtered image.

A cosmic ray deposits its charge in one or a few pixels and has no counterpart in a
median-filtered copy of the frame; the residual ``data - median_filter(data, footprint)``
is compared with the expected noise ``sqrt(variance)`` and pixels above ``sigma`` noise
units are flagged. Only positive residuals count. The footprint runs **along the
dispersion axis** (one row, ``size`` columns): a spectral line is extended across rows by
the fiber image but is as sharp as a cosmic ray along a column, so a square footprint
flags the cores of bright emission lines (the simulated sky lines, for instance) while
the along-row one does not. The runner applies this to science frames only.

Known gap: a cosmic-ray *track* that spans several adjacent pixels along a row survives
the along-row median and is not flagged. Profile-based rejection (Horne extraction, next
plan) closes it.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import median_filter


def flag_cosmic_rays(
    data: np.ndarray,
    variance: np.ndarray,
    mask: np.ndarray | None,
    sigma: float = 8.0,
    size: int | tuple[int, int] = 5,
) -> np.ndarray:
    """Mark pixels more than ``sigma sqrt(variance)`` above the median-filtered image.

    An integer ``size`` means the along-dispersion footprint ``(1, size)``: one row,
    ``size`` columns. A ``(rows, columns)`` tuple is used as given (``(5, 5)`` for a
    square). Pixels that are non-finite or already masked never enter the median: they
    are replaced by a local estimate (the row median, refined by the footprint median of
    the filled image) before filtering, so a cosmic ray next to a bad pixel is still seen
    against its good neighbours; they are never flagged themselves. Returns a new bool
    mask equal to ``mask`` (``None`` for no prior mask) with the flagged pixels added;
    neither input is modified. Tracks of several pixels along a row are not caught (see
    the module docstring).
    """
    data = np.asarray(data, dtype=np.float64)
    noise = np.sqrt(np.clip(np.asarray(variance, dtype=np.float64), 0.0, None))
    bad = ~np.isfinite(data)
    if mask is not None:
        bad |= np.asarray(mask, dtype=bool)
    shape = footprint(size)
    filled = _fill_bad(data, bad)
    filled = np.where(bad, median_filter(filled, size=shape), filled)
    residual = data - median_filter(filled, size=shape)
    flagged = ~bad & (residual > sigma * noise)
    if mask is None:
        return flagged
    return np.asarray(mask, dtype=bool) | flagged


def _fill_bad(data: np.ndarray, bad: np.ndarray) -> np.ndarray:
    """``data`` with the ``bad`` pixels replaced by the median of their row's good pixels
    (the frame's good-pixel median for a row without any)."""
    good = np.where(bad, np.nan, data)
    row_has_good = ~bad.all(axis=1)
    frame_level = np.nanmedian(good) if row_has_good.any() else 0.0
    row_level = np.full(data.shape[0], frame_level)
    if row_has_good.any():
        row_level[row_has_good] = np.nanmedian(good[row_has_good], axis=1)
    return np.where(bad, row_level[:, None], data)


def footprint(size: int | tuple[int, int]) -> tuple[int, int]:
    """The ``(rows, columns)`` median-filter footprint of ``size``: ``(1, size)`` for an int."""
    if isinstance(size, int):
        return (1, size)
    rows, columns = size
    return (int(rows), int(columns))
