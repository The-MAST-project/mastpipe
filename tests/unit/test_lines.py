"""Tests for the line-finding step: recovery of the simulated arc's line positions, deblending
of close pairs, rejection of off-trace hits, saturation and edge flags, the variance-based
signal-to-noise and the per-trace summary."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.models import LineParams
from mastpipe.core.model import CalibrationSet, Frame, Trace, TraceSet
from mastpipe.sim.simulator import FWHM_PX
from mastpipe.steps.calib import calibrate_frame
from mastpipe.steps.lines import LINE_COLUMNS, find_lines, lines_summary
from tests.conftest import ARC_SIGMA_PX, ARC_STRIP_HALFWIDTH

STRIP_HALFWIDTH = ARC_STRIP_HALFWIDTH
SIGMA_PX = ARC_SIGMA_PX
MIN_LOG_INTENSITY = 1.5
MATCH_TOLERANCE_PX = 0.15
MIN_RECOVERED_FRACTION = 0.6


# ------------------------------------------------------------------------------- fixtures


def make_trace(
    index: int,
    y: float,
    slope: float = 0.0,
    role: str | None = None,
    coeffs: np.ndarray | None = None,
) -> Trace:
    if coeffs is None:
        coeffs = np.array([0.0, slope, y])
    return Trace(
        index=index,
        coeffs=coeffs,
        y_center=y,
        slope=slope,
        rms=0.0,
        n_points=10,
        snr=100.0,
        x_samples=np.arange(10.0),
        y_samples=np.full(10, y),
        kept=np.ones(10, dtype=bool),
        role=role,
    )


def synthetic_arc(
    shape: tuple[int, int],
    trace: Trace,
    lines: list[tuple[float, float]],
    noise: float = 2.0,
    seed: int = 1,
    sigma_x: float = SIGMA_PX,
    sigma_y: float = SIGMA_PX,
) -> np.ndarray:
    """A noisy frame with Gaussian lines of ``(x, peak)`` on ``trace``."""
    ny, nx = shape
    rng = np.random.default_rng(seed)
    image = rng.normal(0.0, noise, size=shape)
    yy, xx = np.mgrid[0:ny, 0:nx]
    y_trace = trace.y(xx)
    for x0, peak in lines:
        image += peak * np.exp(-0.5 * (((xx - x0) / sigma_x) ** 2 + ((yy - y_trace) / sigma_y) ** 2))
    return image


def one_trace(shape: tuple[int, int], y: float, slope: float = 0.0) -> TraceSet:
    return TraceSet(traces=[make_trace(0, y, slope=slope)], shape=shape)


@pytest.fixture(scope="module")
def night(sim_full_arc) -> SimpleNamespace:
    """The simulated full arc, calibrated, with traces built from the truth (session-wide)."""
    return sim_full_arc


@pytest.fixture(scope="module")
def arc_lines(sim_arc_lines) -> tuple[list[pd.DataFrame], list[np.ndarray]]:
    """Line tables and strips of every trace of the simulated full arc (session-wide)."""
    return sim_arc_lines


# ------------------------------------------------------------------------- simulated arc


def test_find_lines_recovers_simulated_arc_positions(night, arc_lines):
    truth = night.truth
    tables, strips = arc_lines
    assert len(tables) == len(night.traces) == len(strips)
    assert strips[0].shape == (2 * STRIP_HALFWIDTH + 1, night.arc.shape[1])
    table = tables[0]
    assert list(table.columns) == list(LINE_COLUMNS)
    assert table["id"].tolist() == list(range(len(table)))

    # every catalogue line that is bright and on the chip should sit within the match
    # tolerance of the nearest *usable* (use_x) line centre
    bright = (truth.line_logI > MIN_LOG_INTENSITY) & truth.lines_on_chip
    expected = truth.line_pixels[0][bright]
    usable = table.loc[table["use_x"], "x_center"].to_numpy(dtype=float)
    assert usable.size > 10
    offsets = np.abs(expected[:, None] - usable[None, :]).min(axis=1)
    found = offsets <= MATCH_TOLERANCE_PX
    assert found.mean() >= MIN_RECOVERED_FRACTION, (
        f"{found.sum()} of {found.size} bright lines recovered, median |dx| = {np.median(offsets[found]):.3f} px"
    )
    # the fitted widths follow the simulator's line image
    assert table.loc[table["use_x"], "fwhm_x_fit"].median() == pytest.approx(FWHM_PX, abs=0.3)
    assert table.loc[table["use_y"], "fwhm_y"].median() == pytest.approx(FWHM_PX, abs=0.5)
    # at most a couple of spurious detections (noise or a cosmic ray) fail the trace check
    on_trace = table["on_trace"]
    assert (~on_trace).sum() <= 2
    assert (table.loc[on_trace, "y_offset"].abs() <= night.params.max_y_offset).all()
    assert not table["saturated"].any()
    assert table["snr"].max() > 100


def test_find_lines_measures_every_trace_in_order(night, arc_lines):
    tables, _ = arc_lines
    for trace, table in zip(night.traces, tables, strict=True):
        assert table["use_x"].sum() > 10
        # the spatial centre of the lines follows this trace, not a neighbour's
        on_trace = table["on_trace"]
        y_expected = trace.y(table.loc[on_trace, "x_peak"].to_numpy(dtype=float))
        assert np.abs(table.loc[on_trace, "y_center"].to_numpy() - y_expected).max() < 1.0


# ------------------------------------------------------------------------- synthetic cases


def test_deblend_two_close_lines():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0, slope=0.01)
    image = synthetic_arc(shape, traceset.traces[0], [(160.0, 3000.0), (164.0, 2000.0)])
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX)
    tables, strips = find_lines(image, None, traceset, params)
    table = tables[0]
    assert len(table) == 2
    assert table["group"].nunique() == 1
    assert (table["n_in_group"] == 2).all()
    assert table["fit_ok"].all()
    assert table["x_fit"].to_numpy() == pytest.approx([160.0, 164.0], abs=0.1)
    assert table["x_center"].to_numpy() == pytest.approx([160.0, 164.0], abs=0.1)
    # a blend has a meaningless numeric fwhm_x but a valid fitted width and y width
    assert table["fwhm_x"].isna().all()
    assert table["fwhm_x_fit"].to_numpy() == pytest.approx(FWHM_PX, abs=0.2)
    assert not table["isolated"].any()
    assert table["use_x"].all() and table["use_y"].all()
    # the strip is straightened: the trace (y = 21.6 at x = 160, rounded to row 22) sits on
    # the middle row, and the spatial centre is reported on the detector
    strip = strips[0]
    assert strip.shape == (2 * STRIP_HALFWIDTH + 1, shape[1])
    assert int(np.argmax(strip[:, 160])) == STRIP_HALFWIDTH
    assert table["y_center"].to_numpy() == pytest.approx([21.6, 21.64], abs=0.1)
    assert table["y_offset"].abs().max() < 0.1


def test_deblend_splits_a_shoulder():
    # 2.4 px apart at 3 px FWHM: find_peaks sees one maximum with a shoulder and the
    # residual-driven deblender must add the second component
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    image = synthetic_arc(shape, traceset.traces[0], [(150.0, 4000.0), (152.4, 2000.0)])
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX)
    tables, _ = find_lines(image, None, traceset, params)
    table = tables[0]
    assert len(table) == 2
    assert (table["n_in_group"] == 2).all()
    assert table["x_center"].to_numpy() == pytest.approx([150.0, 152.4], abs=0.15)


def test_cosmic_ray_off_trace_rejected():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    image = synthetic_arc(shape, traceset.traces[0], [(100.0, 3000.0)])
    image[23, 200] += 20000.0  # a hot pixel 3 rows above the trace, inside the strip
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX, max_y_offset=1.5)
    tables, _ = find_lines(image, None, traceset, params)
    table = tables[0].set_index("x_peak")
    assert set(table.index) == {100, 200}
    assert table.loc[100, "on_trace"] and table.loc[100, "use_x"] and table.loc[100, "use_y"]
    assert table.loc[200, "y_offset"] == pytest.approx(3.0, abs=0.3)
    assert not table.loc[200, "on_trace"]
    assert not table.loc[200, "use_x"] and not table.loc[200, "use_y"]
    # None disables the check
    tables, _ = find_lines(image, None, traceset, params.model_copy(update={"max_y_offset": None}))
    assert tables[0]["on_trace"].all()


def test_saturated_line_not_usable():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    image = synthetic_arc(shape, traceset.traces[0], [(100.0, 50000.0), (200.0, 5000.0)])
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX, saturation=30000.0)
    tables, _ = find_lines(image, None, traceset, params)
    table = tables[0].set_index("x_peak")
    assert table.loc[100, "saturated"]
    assert not table.loc[100, "use_x"]
    assert not table.loc[200, "saturated"]
    assert table.loc[200, "use_x"]
    # without a saturation level nothing is flagged
    tables, _ = find_lines(image, None, traceset, params.model_copy(update={"saturation": None}))
    assert not tables[0]["saturated"].any()
    assert tables[0]["use_x"].all()


def test_saturated_lines_are_flagged_through_the_mask(sim_full_arc, highspec_effective):
    """A raw pixel at the detector's saturation inside one line: ``calibrate_frame`` masks it
    and ``find_lines`` with that mask marks the line ``saturated`` (so not ``use_x``), even
    though ``LineParams.saturation`` is unset; without the mask nothing is flagged."""
    detector = highspec_effective.detector
    arc, traces, params = sim_full_arc.arc, sim_full_arc.traces, sim_full_arc.params
    assert params.saturation is None
    (table,), _ = find_lines(arc.data, arc.variance, TraceSet([traces.traces[0]], arc.shape), params)
    clean = table.loc[table["use_x"] & table["isolated"]]
    x0 = int(clean["x_peak"].iloc[len(clean) // 2])
    trace = traces.traces[0]
    raw = np.array(arc.data, dtype=np.float32)
    raw[int(np.round(trace.y(x0))), x0] = 65535.0
    calib = CalibrationSet(
        bias=None,
        dark_rate=0.0,
        background=None,
        t_background=None,
        gain=detector.gain,
        read_noise=detector.read_noise,
        shape=arc.shape,
    )
    frame = calibrate_frame(Frame(data=raw, meta=arc.meta), calib, False, saturation=detector.saturation)
    assert frame.mask is not None and int(frame.mask.sum()) == 1

    one = TraceSet([trace], arc.shape)
    (with_mask,), _ = find_lines(frame.data, frame.variance, one, params, mask=frame.mask)
    (without,), _ = find_lines(frame.data, frame.variance, one, params)
    distance = (with_mask["x_center"] - x0).abs()
    nearest = with_mask.loc[distance.idxmin()]
    assert nearest["saturated"] and not nearest["use_x"]
    # only lines whose window can reach the pixel are touched
    assert not with_mask.loc[distance > 2 * params.x_halfwidth + 1, "saturated"].any()
    assert not without["saturated"].any()
    # a mask of the wrong shape is an error
    with pytest.raises(ValueError, match="mask"):
        find_lines(frame.data, None, one, params, mask=frame.mask[:-1])


def test_edge_lines_flagged():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    lines = [(4.0, 3000.0), (150.0, 3000.0), (295.0, 3000.0)]
    image = synthetic_arc(shape, traceset.traces[0], lines)
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX, edge_margin=10)
    tables, _ = find_lines(image, None, traceset, params)
    table = tables[0].set_index("x_peak")
    assert set(table.index) == {4, 150, 295}
    assert table.loc[4, "near_edge"] and table.loc[295, "near_edge"]
    assert not table.loc[150, "near_edge"]
    assert not table.loc[4, "use_x"] and not table.loc[295, "use_x"]
    assert not table.loc[4, "use_y"] and not table.loc[295, "use_y"]
    assert table.loc[150, "use_x"] and table.loc[150, "use_y"]
    # the window of a line at the very edge is cut by the detector: incomplete
    assert not table.loc[4, "complete"]
    assert table.loc[150, "complete"]


def test_snr_uses_variance_when_given():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    image = synthetic_arc(shape, traceset.traces[0], [(150.0, 3000.0)])
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX)
    variance = np.full(shape, 25.0)
    (without,), _ = find_lines(image, None, traceset, params)
    (with_var,), _ = find_lines(image, variance, traceset, params)
    assert len(without) == len(with_var) == 1
    # the strip sums 2 hw + 1 rows of variance 25 under the peak column
    expected = with_var["height"].iloc[0] / np.sqrt(25.0 * (2 * STRIP_HALFWIDTH + 1))
    assert with_var["snr"].iloc[0] == pytest.approx(expected)
    assert with_var["snr"].iloc[0] != pytest.approx(without["snr"].iloc[0])
    assert without["snr"].iloc[0] > params.threshold_sigma


def test_empty_trace_gives_empty_table():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    image = np.random.default_rng(3).normal(0.0, 2.0, size=shape)
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX)
    tables, strips = find_lines(image, None, traceset, params)
    assert len(tables[0]) == 0
    assert list(tables[0].columns) == list(LINE_COLUMNS)
    summary = lines_summary(tables)
    assert summary["found"].tolist() == [0]
    assert np.isnan(summary["median FWHM_x [px]"].iloc[0])


def test_lines_summary_counts():
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    lines = [(60.0, 3000.0), (150.0, 3000.0), (154.0, 2000.0), (240.0, 3000.0)]
    image = synthetic_arc(shape, traceset.traces[0], lines)
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX)
    tables, _ = find_lines(image, None, traceset, params)
    summary = lines_summary(tables)
    assert isinstance(summary, pd.DataFrame)
    row = summary.iloc[0]
    assert row["trace"] == 0
    assert row["found"] == 4
    assert row["blended"] == 2
    assert row["use_x"] == 4
    assert row["use_y"] == 4
    assert row["median FWHM_x [px]"] == pytest.approx(FWHM_PX, abs=0.3)
    assert row["median FWHM_y [px]"] == pytest.approx(FWHM_PX, abs=0.5)


def test_trace_near_detector_edge_still_usable():
    # the default strip (halfwidth 9) around a trace at y = 8 reaches below row 0: the
    # off-detector rows must not bar the lines from use_x
    shape = (41, 300)
    traceset = one_trace(shape, 8.0)
    lines = [(60.0, 3000.0), (150.0, 3000.0), (240.0, 3000.0)]
    image = synthetic_arc(shape, traceset.traces[0], lines)
    params = LineParams(sigma_guess=SIGMA_PX)
    tables, _ = find_lines(image, None, traceset, params)
    table = tables[0]
    assert len(table) == 3
    assert table["fit_ok"].all() and table["on_trace"].all()
    assert table["use_x"].all()
    # with the whole strip as the measuring window the off-detector rows make it incomplete
    assert not table["complete"].any()
    assert not table["use_y"].any()
    # measuring rows that stay on the detector make the lines complete and usable in y
    tables, _ = find_lines(image, None, traceset, params.model_copy(update={"y_halfwidth": 3}))
    assert tables[0]["complete"].all()
    assert tables[0]["use_y"].all()


def test_spike_next_to_line_not_deblended():
    # a single-pixel spike 3 px from a line is a residual peak the deblender would fit as a
    # component pinned at the sigma lower bound; it must be rejected on width
    shape = (41, 300)
    traceset = one_trace(shape, 20.0)
    image = synthetic_arc(shape, traceset.traces[0], [(150.0, 3000.0)])
    image[20, 153] += 400.0
    params = LineParams(halfwidth=STRIP_HALFWIDTH, sigma_guess=SIGMA_PX)
    tables, _ = find_lines(image, None, traceset, params)
    table = tables[0]
    assert len(table) == 1
    assert table["n_in_group"].iloc[0] == 1
    assert table["x_center"].iloc[0] == pytest.approx(150.0, abs=0.1)
    assert table["use_x"].iloc[0]
