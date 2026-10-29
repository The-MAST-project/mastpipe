"""Master calibrations and detector calibration: a port of the notebook's calibration cell.

Frame model, in ADU: ``raw = bias + dark_rate t + stray t / t_bkg + signal``. The
notebook subtracts a *background* frame (an exposure of the science length with the
lamps off, which carries the dark current and any stray light) scaled by exposure time,

    ``out = raw - bias - (background - bias) t / t_bkg``,

and falls back to the dark current alone when there is no background,

    ``out = raw - bias - dark_rate t``,  ``dark_rate = (dark - bias) / t_dark``.

The variance of a calibrated pixel, in ADU², is ``max(signal_e, 0) / gain² + RN² / gain²``
with ``signal_e = gain (raw - bias)`` the bias-subtracted signal in electrons and ``RN``
the read noise in electrons. The mask marks saturated pixels (in the raw frame or in any
master) and non-finite values.

Everything here is pure: arrays, dataclasses and plain parameters in, dataclasses out.
"""

from __future__ import annotations

from collections.abc import Iterable
from datetime import datetime, timedelta
from typing import Any

import numpy as np

from mastpipe.config.models import CalibParams, DetectorConfig
from mastpipe.core.model import CalibrationSet, Frame, FrameMeta, MasterCalib

STACK_METHODS = ("median", "mean")
MEDIAN_EFFICIENCY = 1.253
"""Standard error of a median relative to that of a mean, ``sqrt(pi / 2)`` for Gaussian noise."""
MAD_TO_SIGMA = 1.4826
"""The median absolute deviation of Gaussian noise scaled to its standard deviation."""


# ------------------------------------------------------------------------------- stacking


def stack_frames(
    frames: list[np.ndarray], method: str = "median", clip_sigma: float | None = None
) -> tuple[np.ndarray, np.ndarray]:
    """Combine same-shaped images pixel by pixel into a master and its variance.

    ``method`` is ``"median"`` or ``"mean"``. Non-finite values never take part; with
    ``clip_sigma`` set, values further than ``clip_sigma`` robust standard deviations
    (``1.4826 MAD``) from the per-pixel median are dropped as well before combining.

    The variance of the master is estimated from the per-pixel sample standard deviation
    ``s`` (``ddof = 1``) of the values that took part: ``s² / n`` for a mean and
    ``1.253² s² / n`` for a median (the median of ``n`` Gaussian samples scatters
    ``sqrt(pi / 2)`` more than their mean). The sample standard deviation is used rather
    than a MAD-based scale because the clipping, when requested, has already removed the
    outliers, and because for the small stacks of a night (three to five frames) the MAD
    is a poorer estimate of the scatter than ``s``. Pixels with fewer than two valid
    values get a variance of ``0``; pixels with none become ``NaN`` in the master.

    Both outputs are ``float32``. ``ValueError`` for an empty list or an unknown method.
    """
    if method not in STACK_METHODS:
        raise ValueError(f"unknown stack method {method!r}, expected one of {STACK_METHODS}")
    if not frames:
        raise ValueError("stack_frames needs at least one frame")
    stack = np.stack([np.asarray(frame, dtype=np.float64) for frame in frames])
    good = np.isfinite(stack)
    if clip_sigma is not None:
        centre = _nan_median(stack, good)
        deviation = np.abs(stack - centre)
        scale = MAD_TO_SIGMA * _nan_median(deviation, good)
        good &= deviation <= clip_sigma * np.maximum(scale, np.finfo(np.float64).tiny)

    count = good.sum(axis=0)
    valid = np.where(good, stack, 0.0)
    mean = valid.sum(axis=0) / np.maximum(count, 1)
    squares = np.where(good, (stack - mean) ** 2, 0.0).sum(axis=0)
    sample_variance = np.where(count > 1, squares / np.maximum(count - 1, 1), 0.0)
    variance = sample_variance / np.maximum(count, 1)
    if method == "median":
        master = _nan_median(stack, good)
        variance = variance * MEDIAN_EFFICIENCY**2
    else:
        master = np.where(count > 0, mean, np.nan)
    return master.astype(np.float32), variance.astype(np.float32)


def _nan_median(stack: np.ndarray, good: np.ndarray) -> np.ndarray:
    """Per-pixel median over the values marked ``good``; ``NaN`` where none is."""
    any_good = good.any(axis=0)
    filled = np.where(good, stack, np.nan)
    filled[:, ~any_good] = 0.0  # keeps nanmedian silent on all-NaN slices
    median = np.nanmedian(filled, axis=0)
    median[~any_good] = np.nan
    return median


# -------------------------------------------------------------------------------- masters


def build_master(frames: list[Frame], kind: str, params: CalibParams) -> MasterCalib:
    """Stack ``frames`` with ``params.stack_method`` / ``params.clip_sigma`` into a master.

    The validity record holds the keys that decide where the master applies: ``camera``,
    ``binning`` and ``readout_mode`` (which every input must share, ``ValueError``
    otherwise), the mean ``detector_temperature``, the mean ``mid_time`` (ISO-8601) and
    the ``exposure_time`` (the common value, or the mean when they differ); a value no
    input provides is ``None``. The provenance lists the input checksums in order, the
    number of frames and the stacking method. ``ValueError`` for no frames or for frames
    of different shapes.
    """
    if not frames:
        raise ValueError(f"master {kind}: no frames to stack")
    metas = [frame.meta for frame in frames]
    for key in ("camera", "binning", "readout_mode"):
        values = {getattr(meta, key) for meta in metas}
        if len(values) != 1:
            raise ValueError(f"master {kind}: inputs mix {key} values {sorted(map(str, values))}")
    shapes = {frame.data.shape for frame in frames}
    if len(shapes) != 1:
        raise ValueError(f"master {kind}: inputs have different shapes {sorted(shapes)}")
    data, variance = stack_frames(
        [frame.data for frame in frames], params.stack_method, params.clip_sigma
    )
    validity: dict[str, Any] = {
        "camera": metas[0].camera,
        "binning": metas[0].binning,
        "readout_mode": metas[0].readout_mode,
        "detector_temperature": _mean_or_none(meta.detector_temperature for meta in metas),
        "mid_time": _mean_time_or_none(meta.mid_time for meta in metas),
        "exposure_time": _mean_or_none(meta.exposure_time for meta in metas),
    }
    provenance = {
        "inputs": [meta.checksum for meta in metas],
        "n_frames": len(frames),
        "method": params.stack_method,
    }
    return MasterCalib(
        kind=kind,
        data=data,
        validity=validity,
        n_frames=len(frames),
        provenance=provenance,
        variance=variance,
    )


def _mean_or_none(values: Iterable[float | None]) -> float | None:
    present = [float(value) for value in values if value is not None]
    return float(np.mean(present)) if present else None


def _mean_time_or_none(values: Iterable[str | None]) -> str | None:
    """Mean of ISO-8601 times as ISO-8601, ``None`` when no time is given."""
    times = [datetime.fromisoformat(value) for value in values if value is not None]
    if not times:
        return None
    origin = times[0]
    offset = sum((time - origin for time in times), timedelta(0)) / len(times)
    return (origin + offset).isoformat()


# ---------------------------------------------------------------------- calibration set


def build_calibration_set(
    bias: MasterCalib | None,
    dark: MasterCalib | None,
    background: MasterCalib | None,
    detector: DetectorConfig,
) -> CalibrationSet:
    """The calibration inputs of one detector from its masters.

    ``dark_rate = (dark - bias) / t_dark`` (ADU per second, per pixel) when both a dark and
    a bias master exist and the dark's ``exposure_time`` is positive, ``0.0`` otherwise;
    ``t_background`` is the background master's ``exposure_time``. The gain and read noise
    come from ``detector``. ``ValueError`` when no master is given at all.
    """
    masters = [master for master in (bias, dark, background) if master is not None]
    if not masters:
        raise ValueError("a calibration set needs at least one master frame")
    shape = tuple(int(n) for n in masters[0].data.shape)
    dark_rate: np.ndarray | float = 0.0
    if dark is not None and bias is not None:
        t_dark = dark.validity.get("exposure_time") or 0.0
        if t_dark > 0:
            dark_rate = ((dark.data - bias.data) / float(t_dark)).astype(np.float32)
    t_background = None if background is None else background.validity.get("exposure_time")
    return CalibrationSet(
        bias=None if bias is None else bias.data,
        dark_rate=dark_rate,
        background=None if background is None else background.data,
        t_background=None if t_background is None else float(t_background),
        gain=detector.gain,
        read_noise=detector.read_noise,
        shape=shape,  # type: ignore[arg-type]
    )


# ------------------------------------------------------------------------- calibration


def calibration_steps(
    calib: CalibrationSet, t_frame: float, subtract_background: bool
) -> list[str]:
    """The operations :func:`calibrate_frame` applies, as short strings for provenance."""
    steps: list[str] = []
    if calib.bias is not None:
        steps.append("bias")
    if subtract_background:
        _check_background(calib)
        scale = float(t_frame) / float(calib.t_background)
        steps.append(
            f"background x {scale:.6g} (t {t_frame:g} s / t_background {calib.t_background:g} s)"
        )
    elif np.any(np.asarray(calib.dark_rate) != 0.0):
        steps.append(f"dark_rate x {float(t_frame):g} s")
    return steps


def calibrate_frame(
    frame: Frame,
    calib: CalibrationSet,
    subtract_background: bool = True,
    saturation: float | None = None,
) -> Frame:
    """Bias/background- or bias/dark-subtract ``frame`` and attach its variance and mask.

    With ``subtract_background`` the output is ``raw - bias - (background - bias) t / t_bkg``
    (``ValueError`` when the set has no background with an exposure time); without it,
    ``raw - bias - dark_rate t``. A missing bias master subtracts nothing for that term.
    The variance is ``max(gain (raw - bias), 0) / gain² + RN² / gain²`` in ADU²; the mask
    marks pixels at or above ``saturation`` in the raw frame, the bias or the background
    master (when a level is given), non-finite output values and whatever ``frame.mask``
    already marks. The dark master is held as a rate, so a pixel saturated in it is not
    recoverable here and is not masked. The data are ``float32``; the input frame is left
    untouched and the output shares its metadata. ``ValueError`` when the frame has no
    exposure time or does not match the calibration shape.
    """
    t = frame.meta.exposure_time
    if t is None:
        raise ValueError(f"{_label(frame.meta)}: exposure time unknown, cannot calibrate")
    if subtract_background:
        _check_background(calib)
    raw = np.asarray(frame.data, dtype=np.float32)
    if raw.shape != tuple(calib.shape):
        raise ValueError(
            f"{_label(frame.meta)}: shape {raw.shape} differs from calibration {calib.shape}"
        )
    bias = np.zeros_like(raw) if calib.bias is None else np.asarray(calib.bias, dtype=np.float32)
    signal = raw - bias
    if subtract_background:
        background = np.asarray(calib.background, dtype=np.float32)
        out = signal - (background - bias) * (float(t) / float(calib.t_background))
    else:
        out = signal - np.asarray(calib.dark_rate, dtype=np.float32) * float(t)
    gain, read_noise = float(calib.gain), float(calib.read_noise)
    signal_e = np.clip(signal, 0.0, None) * gain
    variance = (signal_e / gain**2 + (read_noise / gain) ** 2).astype(np.float32)
    mask = ~np.isfinite(out)
    if saturation is not None:
        mask |= raw >= saturation
        for master in (calib.bias, calib.background):
            if master is not None:
                mask |= np.asarray(master) >= saturation
    if frame.mask is not None:
        mask |= frame.mask
    return Frame(data=out.astype(np.float32), meta=frame.meta, variance=variance, mask=mask)


def _check_background(calib: CalibrationSet) -> None:
    if calib.background is None or not calib.t_background:
        raise ValueError(
            "background subtraction requested but the calibration set has no background "
            "master with a positive exposure time"
        )


def _label(meta: FrameMeta) -> str:
    return meta.source or meta.exposure_type
