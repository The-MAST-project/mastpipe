"""Tests for the cosmic-ray mask, on hand arrays and on a simulated science frame."""

import numpy as np

from mastpipe.config.models import CalibParams
from mastpipe.core.model import Frame
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.cosmics import flag_cosmic_rays


def test_flag_cosmic_rays_on_hand_array():
    data = np.full((9, 9), 10.0)
    variance = np.full((9, 9), 4.0)  # sigma = 2
    data[4, 4] = 100.0  # 45 sigma above the background
    data[1, 1] = 20.0  # 5 sigma: below the default threshold
    mask = np.zeros((9, 9), dtype=bool)
    mask[7, 7] = True
    data_before, mask_before = data.copy(), mask.copy()

    out = flag_cosmic_rays(data, variance, mask, sigma=8.0, size=5)
    assert out.dtype == np.bool_ and out.shape == data.shape
    assert out[4, 4] and out[7, 7] and not out[1, 1]
    assert out.sum() == 2
    assert np.array_equal(data, data_before) and np.array_equal(mask, mask_before)
    assert out is not mask

    # a lower threshold catches the 5 sigma pixel too; a negative spike is never a cosmic
    assert flag_cosmic_rays(data, variance, mask, sigma=4.0)[1, 1]
    data[2, 6] = -100.0
    assert not flag_cosmic_rays(data, variance, mask)[2, 6]
    # no input mask: only the flagged pixels come back
    assert flag_cosmic_rays(data, variance, None).sum() == 1
    # a (rows, columns) footprint is accepted as well
    assert flag_cosmic_rays(data, variance, None, size=(1, 5)).sum() == 1


def test_cosmic_rays_flagged_on_simulated_science(sim_night, highspec_effective):
    out_dir, truth = sim_night
    detector = highspec_effective.detector
    refs = scan_folder(out_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list] = {}
    for ref in refs:
        by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    params = CalibParams()
    bias = build_master(by_type["bias"], "bias", params)
    dark = build_master(by_type["dark"], "dark", params)
    background = build_master(by_type["background"], "background", params)
    calib = build_calibration_set(bias, dark, background, detector)

    science = by_type["science"][0]
    hits = truth.extras["cosmic_rays"][science.meta.source]
    planted = {(int(y), int(x)) for y, x in hits}
    assert planted  # the simulator planted something to find

    clean = calibrate_frame(
        science, calib, subtract_background=True, saturation=detector.saturation
    )
    flagged = flag_cosmic_rays(
        clean.data, clean.variance, clean.mask, params.cosmic_sigma, params.cosmic_size
    )
    new = flagged & ~clean.mask
    found = sum(new[y, x] for y, x in planted)
    assert found >= 0.9 * len(planted)
    n_pixels = clean.data.size
    false_on_science = new.sum() - found
    assert false_on_science < 1e-3 * n_pixels

    # the same frame with the planted hits removed is clean: fewer than 0.1 % flagged
    spotless = science.data.copy()
    for y, x in hits:
        spotless[y, x] -= truth.extras["cosmic_ray_adu"]
    spotless_frame = calibrate_frame(
        Frame(data=spotless, meta=science.meta), calib, True, detector.saturation
    )
    spotless_flags = (
        flag_cosmic_rays(spotless_frame.data, spotless_frame.variance, spotless_frame.mask)
        & ~spotless_frame.mask
    )
    assert spotless_flags.sum() < 1e-3 * n_pixels
