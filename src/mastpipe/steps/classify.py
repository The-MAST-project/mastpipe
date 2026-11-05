"""Content classification of a ``light`` frame: an arc, a science frame or no signal.

The control software writes ``IMAGETYP = 'light'`` for every shutter-open exposure, lamp
on or off, and no lamp card, so the header cannot tell a Th-Ar arc from a science frame.
:func:`classify_light` tells them apart from the calibrated frame itself. It runs before
any trace is known (on a night without a flat the full arc is the trace source), so the
statistic is frame-level: the frame collapsed along the spatial axis into a 1-D spectrum.

**The collapsed spectrum.** Masked and non-finite pixels count as zero. Every column first
loses its running median along the spatial axis over :data:`BASELINE_ROWS` rows: the
diffuse light that fills the detector between fibers (ambient light, scattered light, a
residual bias or background offset) and is no fiber's signal. The laboratory frames need
this: their lamp-off ambient light is about 30 times the Th-Ar line flux of a focus frame
and has broad spectral structure of its own. A fiber image is a few rows wide, so the
median over many more rows is the level between fibers and the fiber light survives
(partly, on a densely packed layout, but in proportion for continuum and lines alike).
The columns are then summed, and a running median of :data:`SPIKE_COLUMNS` columns removes
single-column spikes: the cosmic rays of an uncleaned frame, which would otherwise be
narrow peaks.

**The statistic** (:func:`arc_likeness`). The continuum is a running median of
``continuum_window`` columns of the spectrum and ``sigma`` the robust scatter
(:func:`~mastpipe.steps.trace.robust_sigma`) of the residual about it, the noise of the
spectrum between lines.

``has_signal``
    the :data:`SIGNAL_PERCENTILE` percentile of the spectrum exceeds ``signal_sigma ·
    sigma``: fiber light, continuum or lines, well above the noise in more than a few
    columns.
``n_peaks``
    the runs of columns whose residual exceeds ``signal_sigma · sigma``, each grown by
    ``peak_halfwidth_px`` columns on either side (overlapping runs merge).
``peak_fraction``
    the positive residual inside those peaks over the positive flux of the whole
    spectrum: the share of the fiber light that sits in narrow emission lines. A Th-Ar
    spectrum on a dark continuum puts nearly all of it there; a science frame with sky
    emission lines on a continuum, or a flat, a few per cent at most.

A frame without signal is ``unknown`` (``content:no_signal``); one whose peak fraction
reaches ``arc_peak_fraction`` is an ``arc`` (``content:arc_like``); any other is
``science`` (``content:continuum``). A flat is a continuum and comes out ``science``: only
its header type or a manifest tells it apart. A frame lit by diffuse light alone, as the
laboratory lamp-off background, has no signal.

Everything here is pure: arrays, a frame and a parameter model in, plain values out.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import binary_dilation, label, median_filter

from mastpipe.config.models import ClassifyParams
from mastpipe.core.model import Frame
from mastpipe.steps.trace import robust_sigma

BASELINE_ROWS = 51
"""Rows of the per-column running median removed as diffuse light before collapsing: many
times a fiber image's width (a few rows), well below the detector height."""
SPIKE_COLUMNS = 3
"""Columns of the running median that removes single-column spikes from the spectrum."""
SIGNAL_PERCENTILE = 99.0
"""Percentile of the collapsed spectrum compared with ``signal_sigma · sigma``."""

ARC_LIKE = "content:arc_like"
CONTINUUM = "content:continuum"
NO_SIGNAL = "content:no_signal"


@dataclass(frozen=True)
class ArcLikeness:
    """The content statistic of one frame (see the module docstring)."""

    has_signal: bool
    peak_fraction: float
    n_peaks: int


def collapsed_spectrum(data: np.ndarray, mask: np.ndarray | None) -> np.ndarray:
    """The spatially collapsed spectrum of a 2-D frame whose second axis is the dispersion
    axis: diffuse light removed, columns summed, single-column spikes removed."""
    values = np.asarray(data, dtype=float)
    good = np.isfinite(values)
    if mask is not None:
        good &= ~np.asarray(mask, dtype=bool)
    filled = np.where(good, values, 0.0)
    rows = min(BASELINE_ROWS, filled.shape[0])
    baseline = median_filter(filled, size=(rows, 1), mode="nearest")
    spectrum = (filled - baseline).sum(axis=0)
    return median_filter(spectrum, size=SPIKE_COLUMNS, mode="nearest")


def arc_likeness(data: np.ndarray, mask: np.ndarray | None, params: ClassifyParams) -> ArcLikeness:
    """``has_signal``, ``peak_fraction`` and ``n_peaks`` of a calibrated frame (``mask``:
    ``True`` = bad pixel), see the module docstring."""
    spectrum = collapsed_spectrum(data, mask)
    window = min(params.continuum_window, spectrum.size)
    residual = spectrum - median_filter(spectrum, size=window, mode="nearest")
    sigma = robust_sigma(residual)
    if not np.isfinite(sigma) or sigma <= 0.0:
        return ArcLikeness(has_signal=False, peak_fraction=0.0, n_peaks=0)
    threshold = params.signal_sigma * sigma
    has_signal = bool(np.percentile(spectrum, SIGNAL_PERCENTILE) > threshold)
    core = residual > threshold
    peaks = core
    if core.any() and params.peak_halfwidth_px > 0:
        grow = np.ones(2 * params.peak_halfwidth_px + 1, dtype=bool)
        peaks = binary_dilation(core, structure=grow)
    n_peaks = int(label(peaks)[1])
    total = float(spectrum.clip(0.0).sum())
    in_peaks = float(residual[peaks].clip(0.0).sum())
    fraction = in_peaks / total if total > 0.0 else 0.0
    return ArcLikeness(has_signal=has_signal, peak_fraction=fraction, n_peaks=n_peaks)


def classify_light(frame: Frame, params: ClassifyParams) -> tuple[str, str]:
    """``(exposure_type, rule)`` of a ``light`` frame from its content: ``("arc",
    "content:arc_like")``, ``("science", "content:continuum")`` or ``("unknown",
    "content:no_signal")``."""
    stat = arc_likeness(frame.data, frame.mask, params)
    if not stat.has_signal:
        return "unknown", NO_SIGNAL
    if stat.peak_fraction >= params.arc_peak_fraction:
        return "arc", ARC_LIKE
    return "science", CONTINUUM
