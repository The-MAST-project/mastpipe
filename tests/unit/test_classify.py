"""Tests for the content classification of ``light`` frames (``mastpipe.steps.classify``)."""

from types import SimpleNamespace

import numpy as np
import pytest

from mastpipe.config.models import CalibParams, ClassifyParams
from mastpipe.core.model import Frame, FrameMeta
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.classify import (
    ARC_LIKE,
    CONTINUUM,
    NO_SIGNAL,
    arc_likeness,
    classify_light,
)

P = ClassifyParams()
SHAPE = (200, 1000)
READ_NOISE = 5.0


def synthetic(data: np.ndarray) -> Frame:
    meta = FrameMeta(instrument="highspec", camera="main", exposure_type="light")
    return Frame(data=data.astype(np.float32), meta=meta)


def noise(seed: int = 0) -> np.ndarray:
    return np.random.default_rng(seed).normal(0.0, READ_NOISE, SHAPE)


def fibers(spectrum: np.ndarray, rows=(40, 100, 160), sigma: float = 1.3) -> np.ndarray:
    """``spectrum`` along x on Gaussian fiber images at ``rows``."""
    y = np.arange(SHAPE[0], dtype=float)[:, None]
    profile = sum(np.exp(-0.5 * ((y - row) / sigma) ** 2) for row in rows)
    return profile * spectrum[None, :]


def lines(centres, peak: float, sigma: float = 1.3) -> np.ndarray:
    x = np.arange(SHAPE[1], dtype=float)
    return sum(peak * np.exp(-0.5 * ((x - c) / sigma) ** 2) for c in centres)


LINE_CENTRES = np.linspace(50, 950, 25)


@pytest.fixture(scope="module")
def sim_frames(sim_night, highspec_effective) -> SimpleNamespace:
    """The simulated full arc, a calibration-fiber arc, a science frame and a flat, bias,
    dark and background subtracted (no cosmic-ray flagging), as the runner calibrates a
    ``light`` frame."""
    out_dir, truth = sim_night
    detector = highspec_effective.detector
    refs = {
        ref.meta.source: ref
        for ref in scan_folder(out_dir, "highspec", load_header_map("highspec"))
    }
    rows = truth.frames
    files = {
        kind: rows.loc[rows["exposure_type"] == kind, "file"].tolist()
        for kind in ("bias", "dark", "background", "flat", "science")
    }
    masters = [
        build_master(
            [load_frame(refs[name], detector) for name in files[kind]], kind, CalibParams()
        )
        for kind in ("bias", "dark", "background")
    ]
    calib = build_calibration_set(*masters, detector)

    def calibrated(name: str) -> Frame:
        return calibrate_frame(
            load_frame(refs[name], detector), calib, saturation=detector.saturation
        )

    return SimpleNamespace(
        full_arc=calibrated(rows.loc[rows["arc_kind"] == "full", "file"].iloc[0]),
        cal_arcs=[calibrated(name) for name in rows.loc[rows["arc_kind"] == "calfiber", "file"]],
        science=[calibrated(name) for name in files["science"]],
        flat=calibrated(files["flat"][0]),
    )


# ---------------------------------------------------------------------- simulated frames


def test_simulated_arcs_are_arc_like(sim_frames):
    assert classify_light(sim_frames.full_arc, P) == ("arc", ARC_LIKE)
    assert sim_frames.cal_arcs
    for frame in sim_frames.cal_arcs:  # 5 of 45 traces lit: faint, but lines
        assert classify_light(frame, P) == ("arc", ARC_LIKE)


def test_simulated_science_is_continuum(sim_frames):
    # continuum with absorption dips, sky emission lines on the sky fibers, cosmic rays
    for frame in sim_frames.science:
        assert classify_light(frame, P) == ("science", CONTINUUM)


def test_flat_is_continuum(sim_frames):
    assert classify_light(sim_frames.flat, P)[0] == "science"


def test_no_signal_is_unknown():
    noise_frame = synthetic(noise())  # a bias-subtracted bias: read noise alone
    assert classify_light(noise_frame, P) == ("unknown", NO_SIGNAL)
    assert classify_light(synthetic(np.zeros(SHAPE)), P) == ("unknown", NO_SIGNAL)


# --------------------------------------------------------------------- synthetic arrays


def test_lines_on_diffuse_light_are_arc_like():
    """Emission lines in three fibers over a diffuse light 30 times their flux, with broad
    structure along x (the laboratory focus frames): the diffuse light is not signal."""
    x = np.arange(SHAPE[1], dtype=float)
    diffuse = 300.0 + 100.0 * np.sin(x / 150.0)
    data = noise() + diffuse[None, :] + fibers(lines(LINE_CENTRES, 400.0))
    stat = arc_likeness(data, None, P)
    assert stat.has_signal and stat.peak_fraction > 0.9 and stat.n_peaks == len(LINE_CENTRES)
    assert classify_light(synthetic(data), P) == ("arc", ARC_LIKE)


def test_diffuse_light_alone_is_no_signal():
    x = np.arange(SHAPE[1], dtype=float)
    diffuse = 300.0 + 100.0 * np.sin(x / 150.0)
    data = noise() + diffuse[None, :]
    assert not arc_likeness(data, None, P).has_signal
    assert classify_light(synthetic(data), P) == ("unknown", NO_SIGNAL)


def test_continuum_with_a_few_emission_lines_is_continuum():
    spectrum = 200.0 + lines(LINE_CENTRES[::8], 300.0)
    stat = arc_likeness(noise() + fibers(spectrum), None, P)
    assert stat.has_signal and stat.n_peaks == len(LINE_CENTRES[::8])
    assert stat.peak_fraction < 0.1


def test_cosmic_rays_alone_are_no_signal():
    data = noise()
    rng = np.random.default_rng(1)
    for y, x in zip(rng.integers(0, SHAPE[0], 20), rng.integers(0, SHAPE[1], 20), strict=True):
        data[y, x] += 5000.0
    assert classify_light(synthetic(data), P) == ("unknown", NO_SIGNAL)


def test_masked_pixels_do_not_count():
    data = noise() + fibers(lines(LINE_CENTRES, 400.0))
    block = (slice(98, 103), slice(515, 520))  # a saturated blob between two lines
    data[block] += 1e5
    assert arc_likeness(data, None, P).n_peaks == len(LINE_CENTRES) + 1
    mask = np.zeros(SHAPE, dtype=bool)
    mask[block] = True
    assert arc_likeness(data, mask, P).n_peaks == len(LINE_CENTRES)


def test_thresholds_come_from_the_params():
    data = noise() + fibers(lines(LINE_CENTRES, 400.0))
    stat = arc_likeness(data, None, P)
    strict = ClassifyParams(arc_peak_fraction=min(1.0, stat.peak_fraction + 0.01))
    assert classify_light(synthetic(data), strict) == ("science", CONTINUUM)
    deaf = ClassifyParams(signal_sigma=1e6)
    assert classify_light(synthetic(data), deaf) == ("unknown", NO_SIGNAL)
