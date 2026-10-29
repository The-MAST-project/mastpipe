"""Tests for the calibration step: stacking, masters, the calibration set and frame
calibration, checked against hand-computed values."""

import numpy as np
import pytest

from mastpipe.config.models import CalibParams, DetectorConfig
from mastpipe.core.model import CalibrationSet, Frame, FrameMeta, MasterCalib
from mastpipe.steps.calib import (
    build_calibration_set,
    build_master,
    calibrate_frame,
    calibration_steps,
    stack_frames,
)

SHAPE = (4, 6)
DET = DetectorConfig(
    name="main",
    shape=SHAPE,
    gain=2.0,
    read_noise=4.0,
    pixel_size_um=16.0,
    saturation=65000.0,
    crop=None,
)


def frame(
    value: float,
    kind: str,
    t: float,
    temp: float = -50.0,
    cs: str = "c",
    mid_time: str = "2026-10-01T18:00:00",
    **meta_overrides,
) -> Frame:
    meta = FrameMeta(
        instrument="toy",
        camera="main",
        exposure_type=kind,
        exposure_time=t,
        detector_temperature=temp,
        checksum=cs,
        mid_time=mid_time,
        readout_mode="Image",
        **meta_overrides,
    )
    return Frame(data=np.full(SHAPE, value, dtype=np.float32), meta=meta)


# ------------------------------------------------------------------------- stack_frames


def test_stack_median_and_variance():
    values = [10.0, 12.0, 11.0, 9.0, 13.0]
    stack = [np.full((2, 2), v) for v in values]
    stack[2][0, 0] = 1000.0  # one outlier at a single pixel
    master, var = stack_frames(stack, method="median")
    assert master.shape == (2, 2) and var.shape == (2, 2)
    assert master.dtype == np.float32 and var.dtype == np.float32
    # the median ignores the outlier ...
    assert master[0, 0] == pytest.approx(np.median([10.0, 12.0, 1000.0, 9.0, 13.0]))
    assert master[1, 1] == pytest.approx(11.0)
    # ... and the variance is 1.253^2 s^2 / n with s the sample standard deviation
    s2 = np.var(values, ddof=1)  # 2.5
    assert var[1, 1] == pytest.approx(1.253**2 * s2 / 5, rel=1e-4)
    # with clipping, the outlier is dropped before the scatter is measured
    clipped, var_clipped = stack_frames(stack, method="median", clip_sigma=3.0)
    assert clipped[0, 0] == pytest.approx(np.median([10.0, 12.0, 9.0, 13.0]))
    assert var_clipped[0, 0] == pytest.approx(1.253**2 * np.var([10, 12, 9, 13], ddof=1) / 4)
    # the inputs are untouched
    assert stack[2][0, 0] == 1000.0


def test_stack_mean_and_errors():
    stack = [np.full((2, 2), v) for v in (1.0, 2.0, 3.0, 4.0)]
    master, var = stack_frames(stack, method="mean")
    assert master[0, 0] == pytest.approx(2.5)
    assert var[0, 0] == pytest.approx(np.var([1, 2, 3, 4], ddof=1) / 4)
    single, var_single = stack_frames([np.full((2, 2), 7.0)])
    assert single[0, 0] == 7.0 and var_single[0, 0] == 0.0
    with pytest.raises(ValueError):
        stack_frames(stack, method="mode")
    with pytest.raises(ValueError):
        stack_frames([])


def test_stack_ignores_non_finite_values():
    stack = [np.full((2, 2), v) for v in (1.0, 2.0, 3.0)]
    stack[0][0, 0] = np.nan
    master, var = stack_frames(stack)
    assert master[0, 0] == pytest.approx(2.5)
    assert master[1, 1] == pytest.approx(2.0)
    all_bad = [np.full((1, 1), np.nan) for _ in range(3)]
    master, var = stack_frames(all_bad)
    assert np.isnan(master[0, 0]) and var[0, 0] == 0.0


# ------------------------------------------------------------------------- build_master


def test_build_master_validity_and_provenance():
    frames = [
        frame(500.0, "bias", 0.0, temp=-50.0, cs="a", mid_time="2026-10-01T18:00:00"),
        frame(502.0, "bias", 0.0, temp=-52.0, cs="b", mid_time="2026-10-01T18:10:00"),
        frame(504.0, "bias", 0.0, temp=-54.0, cs="c", mid_time="2026-10-01T18:20:00"),
    ]
    master = build_master(frames, "bias", CalibParams())
    assert isinstance(master, MasterCalib)
    assert master.kind == "bias" and master.n_frames == 3
    assert master.data.dtype == np.float32 and master.data.shape == SHAPE
    assert master.data[0, 0] == pytest.approx(502.0)
    assert master.variance is not None and master.variance.shape == SHAPE
    assert master.validity == {
        "camera": "main",
        "binning": (1, 1),
        "readout_mode": "Image",
        "detector_temperature": pytest.approx(-52.0),
        "mid_time": "2026-10-01T18:10:00",
        "exposure_time": 0.0,
    }
    assert master.provenance == {"inputs": ["a", "b", "c"], "n_frames": 3, "method": "median"}

    darks = [frame(530.0, "dark", 300.0, cs="d"), frame(532.0, "dark", 300.0, cs="e")]
    dark = build_master(darks, "dark", CalibParams(stack_method="mean"))
    assert dark.validity["exposure_time"] == pytest.approx(300.0)
    assert dark.data[0, 0] == pytest.approx(531.0)
    assert dark.provenance["method"] == "mean"

    # a dark or background of mixed exposure times is an error, a bias averages them;
    # missing values become None
    mixed = [frame(1.0, "dark", 100.0), frame(1.0, "dark", 300.0)]
    with pytest.raises(ValueError):
        build_master(mixed, "dark", CalibParams())
    with pytest.raises(ValueError):
        build_master(mixed, "background", CalibParams())
    nearly = [frame(1.0, "dark", 300.0), frame(1.0, "dark", 300.0 * (1 + 1e-7))]
    assert build_master(nearly, "dark", CalibParams()).validity["exposure_time"] == (
        pytest.approx(300.0)
    )
    biases = [frame(1.0, "bias", 0.0), frame(1.0, "bias", 0.1)]
    assert build_master(biases, "bias", CalibParams()).validity["exposure_time"] == (
        pytest.approx(0.05)
    )
    bare = Frame(
        data=np.zeros(SHAPE, dtype=np.float32),
        meta=FrameMeta(instrument="toy", camera="main", exposure_type="bias"),
    )
    validity = build_master([bare], "bias", CalibParams()).validity
    assert validity["mid_time"] is None and validity["detector_temperature"] is None
    assert validity["exposure_time"] is None and validity["readout_mode"] is None


def test_build_master_rejects_inconsistent_inputs():
    with pytest.raises(ValueError):
        build_master([], "bias", CalibParams())
    with pytest.raises(ValueError):
        build_master(
            [frame(1.0, "bias", 0.0), frame(1.0, "bias", 0.0, binning=(2, 2))],
            "bias",
            CalibParams(),
        )
    meta = FrameMeta(instrument="toy", camera="red", exposure_type="bias")
    red = Frame(data=np.zeros(SHAPE, dtype=np.float32), meta=meta)
    with pytest.raises(ValueError):
        build_master([frame(1.0, "bias", 0.0), red], "bias", CalibParams())
    with pytest.raises(ValueError):
        build_master(
            [frame(1.0, "bias", 0.0), Frame(np.zeros((2, 2), np.float32), meta)],
            "bias",
            CalibParams(),
        )


# ----------------------------------------------------------------- build_calibration_set


def test_calibration_set_dark_rate_and_background():
    bias = build_master([frame(500.0, "bias", 0.0)], "bias", CalibParams())
    dark = build_master([frame(530.0, "dark", 300.0)], "dark", CalibParams())
    bkg = build_master([frame(560.0, "background", 6.0)], "background", CalibParams())

    calib = build_calibration_set(bias, dark, None, DET)
    assert isinstance(calib, CalibrationSet)
    assert np.allclose(calib.dark_rate, 0.1)
    assert calib.background is None and calib.t_background is None
    assert calib.gain == 2.0 and calib.read_noise == 4.0 and calib.shape == SHAPE
    assert np.array_equal(calib.bias, bias.data)

    assert build_calibration_set(bias, None, None, DET).dark_rate == 0.0
    assert build_calibration_set(None, dark, None, DET).dark_rate == 0.0
    zero_t = build_master([frame(530.0, "dark", 0.0)], "dark", CalibParams())
    assert build_calibration_set(bias, zero_t, None, DET).dark_rate == 0.0

    with_bkg = build_calibration_set(bias, dark, bkg, DET)
    assert np.array_equal(with_bkg.background, bkg.data)
    assert with_bkg.t_background == pytest.approx(6.0)
    with pytest.raises(ValueError):
        build_calibration_set(None, None, None, DET)


# -------------------------------------------------------------------- calibrate_frame


def test_calibrate_background_arithmetic_matches_notebook():
    # raw - bias - (background - bias) * t / t_bkg: 1000 - 500 - 60 * 12 / 6 = 380
    bias = build_master([frame(500.0, "bias", 0.0)], "bias", CalibParams())
    bkg = build_master([frame(560.0, "background", 6.0)], "background", CalibParams())
    dark = build_master([frame(530.0, "dark", 300.0)], "dark", CalibParams())
    calib = build_calibration_set(bias, dark, bkg, DET)
    sci = frame(1000.0, "science", 12.0)
    raw_before = sci.data.copy()

    out = calibrate_frame(sci, calib, subtract_background=True)
    assert out.data.dtype == np.float32 and out.data.shape == SHAPE
    assert np.allclose(out.data, 380.0)
    assert out.meta is sci.meta
    assert np.array_equal(sci.data, raw_before)  # the input is untouched
    # variance from the bias-subtracted signal: gain (raw - bias) = 1000 e-
    # -> 1000 / 4 + 16 / 4 = 254 ADU^2
    assert np.allclose(out.variance, 254.0)
    assert out.mask.dtype == np.bool_ and not out.mask.any()

    steps = calibration_steps(calib, 12.0, True)
    assert steps[0] == "bias" and "background" in steps[1] and "2" in steps[1]
    assert not any("dark" in step for step in steps)

    # the dark path is a different number on the same inputs
    assert np.allclose(calibrate_frame(sci, calib, subtract_background=False).data, 498.8)


def test_calibrate_dark_path_when_no_background():
    bias = build_master([frame(500.0, "bias", 0.0)], "bias", CalibParams())
    dark = build_master([frame(530.0, "dark", 300.0)], "dark", CalibParams())
    calib = build_calibration_set(bias, dark, None, DET)
    sci = frame(1500.0, "science", 100.0)

    out = calibrate_frame(sci, calib, subtract_background=False)
    assert np.allclose(out.data, 1500.0 - 500.0 - 0.1 * 100.0)  # 990 ADU
    assert np.allclose(out.variance, 2.0 * 1000.0 / 4.0 + 16.0 / 4.0)  # 504 ADU^2
    steps = calibration_steps(calib, 100.0, False)
    assert steps[0] == "bias" and "dark" in steps[1] and "100" in steps[1]

    # no dark master: bias only, dark_rate 0
    bias_only = build_calibration_set(bias, None, None, DET)
    assert np.allclose(calibrate_frame(sci, bias_only, subtract_background=False).data, 1000.0)
    assert calibration_steps(bias_only, 100.0, False) == ["bias"]

    # asking for a background that does not exist is an error, not a silent fallback
    with pytest.raises(ValueError):
        calibrate_frame(sci, bias_only, subtract_background=True)
    with pytest.raises(ValueError):
        calibration_steps(bias_only, 100.0, True)
    # so is a frame without an exposure time, or of the wrong shape
    no_t = Frame(
        data=np.zeros(SHAPE, dtype=np.float32),
        meta=FrameMeta(instrument="toy", camera="main", exposure_type="science"),
    )
    with pytest.raises(ValueError):
        calibrate_frame(no_t, calib, subtract_background=False)
    wrong = Frame(data=np.zeros((2, 2), dtype=np.float32), meta=sci.meta)
    with pytest.raises(ValueError):
        calibrate_frame(wrong, calib, subtract_background=False)


def test_variance_and_saturation_mask():
    bias = build_master([frame(500.0, "bias", 0.0)], "bias", CalibParams())
    bkg = build_master([frame(560.0, "background", 6.0)], "background", CalibParams())
    bkg.data[2, 2] = 66000.0  # a pixel saturated in the master
    calib = build_calibration_set(bias, None, bkg, DET)
    sci = frame(1000.0, "science", 12.0)
    sci.data[1, 1] = 70000.0  # saturated in the raw frame
    sci.data[0, 3] = 490.0  # below the bias: variance floors at the read-noise term
    sci.data[3, 5] = np.nan
    sci.mask = np.zeros(SHAPE, dtype=bool)
    sci.mask[3, 0] = True  # already bad on input

    out = calibrate_frame(sci, calib, subtract_background=True, saturation=65000.0)
    assert out.variance[0, 3] == pytest.approx((4.0 / 2.0) ** 2)
    assert out.variance[1, 1] == pytest.approx(2.0 * 69500.0 / 4.0 + 4.0)
    assert out.mask[1, 1] and out.mask[2, 2] and out.mask[3, 5] and out.mask[3, 0]
    assert out.mask.sum() == 4
    assert not out.mask[0, 0]
    assert np.isfinite(out.data[~out.mask]).all()

    # without a saturation level only non-finite and input-masked pixels are flagged
    out_nosat = calibrate_frame(sci, calib, subtract_background=True)
    assert out_nosat.mask.sum() == 2 and out_nosat.mask[3, 5] and out_nosat.mask[3, 0]
