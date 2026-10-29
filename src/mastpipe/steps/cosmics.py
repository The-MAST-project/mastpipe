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
    square). Returns a new bool mask equal to ``mask`` (``None`` for no prior mask) with
    the flagged pixels added; neither input is modified. Tracks of several pixels along a
    row are not caught (see the module docstring).
    """
    data = np.asarray(data, dtype=np.float64)
    noise = np.sqrt(np.clip(np.asarray(variance, dtype=np.float64), 0.0, None))
    residual = data - median_filter(data, size=footprint(size))
    flagged = residual > sigma * noise
    if mask is None:
        return flagged
    return np.asarray(mask, dtype=bool) | flagged


def footprint(size: int | tuple[int, int]) -> tuple[int, int]:
    """The ``(rows, columns)`` median-filter footprint of ``size``: ``(1, size)`` for an int."""
    if isinstance(size, int):
        return (1, size)
    rows, columns = size
    return (int(rows), int(columns))
