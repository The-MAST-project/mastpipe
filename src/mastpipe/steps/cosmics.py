"""Cosmic-ray mask: single-pixel hits standing far above a median-filtered image.

A cosmic ray deposits its charge in one or a few pixels and has no counterpart in a
median-filtered copy of the frame; the residual ``data - median_filter(data, size)`` is
compared with the expected noise ``sqrt(variance)`` and pixels above ``sigma`` noise
units are flagged. Only positive residuals count. The runner applies this to science
frames only: the arc lines of an arc frame are sharp enough in both axes to be flagged.
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
    """Mark pixels more than ``sigma sqrt(variance)`` above the ``size x size`` median.

    ``size`` may also be a ``(rows, columns)`` footprint, e.g. ``(1, 5)`` to compare each
    pixel with its neighbours along the dispersion axis only, which spares the cores of
    bright emission lines on a trace. Returns a new bool mask equal to ``mask`` (``None``
    for no prior mask) with the flagged pixels added; neither input is modified.
    """
    data = np.asarray(data, dtype=np.float64)
    noise = np.sqrt(np.clip(np.asarray(variance, dtype=np.float64), 0.0, None))
    residual = data - median_filter(data, size=size)
    flagged = residual > sigma * noise
    if mask is None:
        return flagged
    return np.asarray(mask, dtype=bool) | flagged
