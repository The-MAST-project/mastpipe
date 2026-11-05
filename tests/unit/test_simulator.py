"""Tests for the synthetic-night simulator, checked against the truth it writes."""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from scipy.ndimage import median_filter

from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.sim.simulator import SimTruth, load_truth, simulate_night
from mastpipe.steps.geometry import grating_geometry, grating_x, nominal_solution

FILE_PATTERN = re.compile(r"^sim_(bias|dark|flat|arc|science|background)_\d{2}\.fits$")


def cropped(path: Path, truth: SimTruth) -> np.ndarray:
    """A simulated frame as float, cropped the way ``load_frame`` crops it."""
    data = np.asarray(fits.getdata(path), dtype=float)
    if truth.crop is None:
        return data
    x_min, x_max, y_min, y_max = truth.crop
    return data[y_min:y_max, x_min:x_max]


def aperture_spectrum(data: np.ndarray, coeffs: np.ndarray, halfwidth: int = 4) -> np.ndarray:
    """Sum of the rows within ``halfwidth`` of the trace centre, column by column."""
    x = np.arange(data.shape[1])
    centre = np.rint(np.polyval(coeffs, x)).astype(int)
    rows = centre[None, :] + np.arange(-halfwidth, halfwidth + 1)[:, None]
    return data[rows, x[None, :]].sum(axis=0)


def centroid(spectrum: np.ndarray, x0: float, halfwidth: int = 4) -> float:
    lo, hi = int(round(x0)) - halfwidth, int(round(x0)) + halfwidth + 1
    x = np.arange(lo, hi)
    weights = spectrum[lo:hi] - np.median(spectrum)
    return float((weights * x).sum() / weights.sum())


def frame_path(out_dir: Path, truth: SimTruth, exposure_type: str, index: int = 0) -> Path:
    rows = truth.frames[truth.frames.exposure_type == exposure_type]
    return out_dir / rows.file.iloc[index]


def background_level(truth: SimTruth, path: Path) -> float:
    return truth.bias_level + truth.dark_rate * float(fits.getheader(path)["EXPOSURE"])


def sim_geometry(effective):
    return grating_geometry(effective.geometry, effective.detector.pixel_size_um)


def brightest_line(truth: SimTruth) -> int:
    """Index of the brightest catalogue line on the chip."""
    return int(np.argmax(np.where(truth.lines_on_chip, truth.line_logI, -np.inf)))


def on_chip(truth: SimTruth, trace_row: int) -> np.ndarray:
    """Undrifted pixels of the catalogue lines the truth marks as on the chip."""
    return truth.line_pixels[trace_row][truth.lines_on_chip]


def test_simulate_night_writes_expected_frames(sim_night, highspec_effective):
    out_dir, truth = sim_night
    frames = truth.frames
    assert list(frames.columns) == [
        "file",
        "exposure_type",
        "arc_kind",
        "mid_time",
        "shift_px",
        "stretch",
        "gradient_px",
    ]
    counts = frames.exposure_type.value_counts().to_dict()
    assert counts == {"bias": 5, "dark": 3, "flat": 3, "arc": 3, "science": 4, "background": 1}
    kinds = frames.arc_kind.tolist()
    assert kinds.count("full") == 1 and kinds.count("calfiber") == 2
    not_arc = frames.exposure_type != "arc"
    assert frames.arc_kind[not_arc].isna().all()
    assert all(FILE_PATTERN.match(name) for name in frames.file)
    assert all((out_dir / name).is_file() for name in frames.file)
    assert frames.file.is_unique

    # time order: bias, dark, flat, full arc, then cal arcs before their science frames
    times = pd.to_datetime(frames.mid_time)
    assert times.is_monotonic_increasing and times.is_unique
    order = frames.exposure_type.tolist()
    expected_order = ["bias"] * 5 + ["dark"] * 3 + ["flat"] * 3 + ["arc", "background"]
    expected_order += ["arc", "science", "science"] * 2
    assert order == expected_order
    assert kinds[11] == "full" and kinds[13] == "calfiber" and kinds[16] == "calfiber"
    for i in frames.index[frames.exposure_type == "science"]:
        assert "calfiber" in kinds[:i]
    assert frames.mid_time.iloc[0] == "2026-10-01T18:00:00"

    # the raw files are uncropped uint16 and carry the HighSpec cards
    science = frame_path(out_dir, truth, "science")
    with fits.open(science) as hdul:
        header = hdul[0].header
        assert hdul[0].data.shape == highspec_effective.detector.shape
        assert hdul[0].data.dtype == np.uint16
    for card in ("IMAGETYP", "EXPOSURE", "DATE-OBS", "TEMP", "HBIN", "VBIN", "READMODE"):
        assert card in header
    for card in ("GRATING", "LAMP", "OBJECT", "SHUTTER"):
        assert card in header
    assert "FIBPRE" not in header  # a MAST fiber preset is a stage position, not a layout
    assert header["IMAGETYP"] == "light" and header["LAMP"] == "off"
    assert header["EXPOSURE"] == pytest.approx(300.0) and header["SHUTTER"] == "open"
    # the header stamps the exposure start; the truth table records the mid-time
    science_row = frames.loc[frames.file == science.name].iloc[0]
    assert pd.Timestamp(header["DATE-OBS"]) == pd.Timestamp(science_row.mid_time) - pd.Timedelta(
        seconds=150
    )
    assert header["GRATING"] == "Mg"
    assert fits.getheader(frame_path(out_dir, truth, "arc"))["LAMP"] == "thar"
    assert fits.getheader(frame_path(out_dir, truth, "flat"))["LAMP"] == "flat"
    assert fits.getheader(frame_path(out_dir, truth, "bias"))["SHUTTER"] == "closed"
    assert fits.getheader(frame_path(out_dir, truth, "dark"))["SHUTTER"] == "closed"
    background = fits.getheader(frame_path(out_dir, truth, "background"))
    assert background["IMAGETYP"] == "background" and background["LAMP"] == "off"
    assert background["SHUTTER"] == "open" and background["EXPOSURE"] == pytest.approx(300.0)
    assert "ARCKIND" not in fits.getheader(frame_path(out_dir, truth, "arc"))

    # ingest classifies every frame from its header alone; a science frame is 'light', the
    # provisional type the runner resolves from its content
    refs = scan_folder(out_dir, "highspec", load_header_map("highspec"))
    by_file = {ref.meta.source: ref for ref in refs}
    assert set(by_file) == set(frames.file)
    for row in frames.itertuples():
        meta = by_file[row.file].meta
        assert meta.exposure_type == (
            "light" if row.exposure_type == "science" else row.exposure_type
        )
        assert meta.classification_rule == "header"
        assert meta.mid_time == row.mid_time
        assert meta.configuration == "Mg" and meta.fiber_preset is None
        assert meta.arc_kind is None
    loaded = load_frame(by_file[science.name], highspec_effective.detector)
    assert loaded.shape == (398, 1598)
    assert tuple(truth.extras["shape"]) == loaded.shape

    # targets: four targets over the twenty science fibers, five telescopes each
    layout = truth.layout
    science_rows = layout[layout.role == "science"]
    assert science_rows.target.nunique() == 4
    assert science_rows.groupby("target").telescope.nunique().tolist() == [5, 5, 5, 5]
    assert layout.loc[layout.role == "arc", "target"].isna().all()
    for target in science_rows.target.unique():
        assert target in header["OBJECT"]
    # every telescope's target is written per unit as OBJnn; calibration frames carry none
    telescopes = science_rows.drop_duplicates("telescope")
    for row in telescopes.itertuples():
        assert header[f"OBJ{row.telescope[-2:]}"] == row.target
    assert len(telescopes) == 20
    science_meta = by_file[science.name].meta
    expected = dict(zip(telescopes.telescope, telescopes.target, strict=True))
    assert dict(science_meta.targets) == expected
    arc_header = fits.getheader(frame_path(out_dir, truth, "arc"))
    assert not any(f"OBJ{n:02d}" in arc_header for n in range(1, 21))


def test_truth_roundtrip(sim_night, highspec_effective):
    out_dir, truth = sim_night
    assert (out_dir / "truth.npz").is_file() and (out_dir / "truth.json").is_file()
    loaded = load_truth(out_dir)
    assert isinstance(loaded, SimTruth)
    pd.testing.assert_frame_equal(loaded.layout, truth.layout)
    pd.testing.assert_frame_equal(loaded.frames, truth.frames)
    for name in (
        "trace_coeffs",
        "lines_nm",
        "line_logI",
        "line_pixels",
        "lines_on_chip",
        "throughput",
        "science_flux",
        "sky_flux",
        "sky_lines_nm",
        "x_c_per_trace",
        "B_per_trace",
    ):
        np.testing.assert_array_equal(getattr(loaded, name), getattr(truth, name), err_msg=name)
    for name in (
        "x_c",
        "B",
        "drift_gradient_px",
        "bias_level",
        "dark_rate",
        "read_noise_e",
        "gain",
        "fwhm_px",
    ):
        assert getattr(loaded, name) == getattr(truth, name), name
    assert loaded.crop == truth.crop
    assert loaded.extras == truth.extras
    assert loaded.trace_coeffs.shape == (45, 3)
    assert loaded.line_pixels.shape == (45, loaded.lines_nm.size)
    assert loaded.lines_on_chip.shape == (loaded.lines_nm.size,)
    assert loaded.lines_on_chip.dtype == np.bool_
    assert loaded.x_c_per_trace.shape == loaded.B_per_trace.shape == (45,)
    assert loaded.science_flux.shape == loaded.sky_flux.shape == (45, loaded.extras["n_pix"])
    n_pix = loaded.extras["n_pix"]
    # on the chip is judged at the nominal solution, which the per-trace pixels jitter about
    geometry = sim_geometry(highspec_effective)
    pixels = grating_x(loaded.lines_nm, loaded.x_c, loaded.B, geometry.alpha_rad, geometry)
    np.testing.assert_array_equal(loaded.lines_on_chip, (pixels >= 10) & (pixels <= n_pix - 10))
    assert 40 < loaded.lines_on_chip.sum() < loaded.lines_nm.size


def test_calfiber_arc_lights_only_arc_traces(sim_night, highspec_effective):
    out_dir, truth = sim_night
    frames = truth.frames
    path = out_dir / frames.file[frames.arc_kind == "calfiber"].iloc[0]
    data = cropped(path, truth)
    level = background_level(truth, path)
    x = np.arange(data.shape[1])
    layout = truth.layout
    science_rows = np.where(layout.role == "science")[0]
    arc_rows = np.where(layout.role == "arc")[0]
    science_pixels = np.concatenate(
        [data[np.rint(np.polyval(truth.trace_coeffs[i], x)).astype(int), x] for i in science_rows]
    )
    assert science_pixels.mean() == pytest.approx(level, abs=1.0)
    assert science_pixels.max() < level + 100.0
    for i in arc_rows:
        spectrum = aperture_spectrum(data, truth.trace_coeffs[i])
        assert spectrum.max() > level * 9 + 5000.0

    # the full arc lights the science traces too, and nothing saturates
    full_path = out_dir / frames.file[frames.arc_kind == "full"].iloc[0]
    full = cropped(full_path, truth)
    for i in science_rows[:3]:
        assert aperture_spectrum(full, truth.trace_coeffs[i]).max() > level * 9 + 5000.0
    assert fits.getdata(full_path).max() < highspec_effective.detector.saturation


def test_science_frames_have_no_arc_lines(sim_night):
    out_dir, truth = sim_night
    frames = truth.frames
    science_path = frame_path(out_dir, truth, "science")
    data = cropped(science_path, truth)
    for y, x in truth.extras["cosmic_rays"][science_path.name]:
        data[y, x] = background_level(truth, science_path)
    fiber = int(np.where(truth.layout.role == "science")[0][0])
    spectrum = aperture_spectrum(data, truth.trace_coeffs[fiber])
    residual = spectrum - median_filter(spectrum, size=31, mode="nearest")
    sigma = 1.4826 * np.median(np.abs(residual - np.median(residual)))
    pixels = np.rint(on_chip(truth, fiber)).astype(int)
    assert len(pixels) > 40
    assert (residual[pixels] < 5.0 * sigma).all()
    # the continuum carries the absorption dips at the expected depth
    assert np.median(spectrum) == pytest.approx(
        9 * background_level(truth, science_path) + 200.0 * truth.throughput[fiber], rel=0.05
    )

    # positive control: the same fiber on the full arc shows the lines
    full = cropped(out_dir / frames.file[frames.arc_kind == "full"].iloc[0], truth)
    arc = aperture_spectrum(full, truth.trace_coeffs[fiber])
    arc_residual = arc - median_filter(arc, size=31, mode="nearest")
    assert (arc_residual[pixels] > 5.0 * sigma).sum() >= 30


def test_drift_interpolates_in_time(sim_night):
    out_dir, truth = sim_night
    frames = truth.frames
    cal = frames.index[frames.arc_kind == "calfiber"].tolist()
    assert len(cal) == 2
    shifts = frames.shift_px
    assert shifts.iloc[cal[0]] < shifts.iloc[cal[1]]
    between = shifts.iloc[cal[0] : cal[1] + 1]
    assert between.is_monotonic_increasing and between.is_unique
    assert shifts.iloc[-1] == pytest.approx(0.6)
    full = frames.index[frames.arc_kind == "full"][0]
    assert shifts.iloc[full] == pytest.approx(0.0)
    assert (shifts.iloc[:full] == 0.0).all()  # before the full arc: the start value
    assert (frames.stretch == 0.0).all()
    # the law is linear in time from the full arc's mid-time to the last frame's
    seconds = (
        pd.to_datetime(frames.mid_time) - pd.Timestamp(frames.mid_time.iloc[full])
    ).dt.total_seconds()
    expected = 0.6 * np.clip(seconds, 0.0, None) / seconds.iloc[-1]
    np.testing.assert_allclose(shifts, expected, atol=1e-9)

    # and the lines on the cal arcs really sit at the truth pixel plus the recorded shift
    fiber = int(np.where(truth.layout.role == "arc")[0][0])
    brightest = int(
        np.argmax(
            np.where(
                np.isin(truth.line_pixels[fiber], on_chip(truth, fiber)), truth.line_logI, -np.inf
            )
        )
    )
    x_true = truth.line_pixels[fiber, brightest]
    for i in cal:
        spectrum = aperture_spectrum(
            cropped(out_dir / frames.file.iloc[i], truth), truth.trace_coeffs[fiber]
        )
        shift = truth.trace_shifts(frames.file.iloc[i])[fiber]  # the frame's shift at this row
        assert centroid(spectrum, x_true) == pytest.approx(x_true + shift, abs=0.1)
    spectrum = aperture_spectrum(
        cropped(out_dir / frames.file.iloc[full], truth), truth.trace_coeffs[fiber]
    )
    assert centroid(spectrum, x_true) == pytest.approx(x_true, abs=0.1)


def test_flat_traces_follow_truth(sim_night):
    out_dir, truth = sim_night
    data = cropped(frame_path(out_dir, truth, "flat"), truth)
    level = background_level(truth, frame_path(out_dir, truth, "flat"))
    x = np.array([200, 800, 1400])
    for i in range(len(truth.layout)):
        y_true = np.polyval(truth.trace_coeffs[i], x)
        for xi, yi in zip(x, y_true, strict=True):
            rows = np.arange(int(round(yi)) - 3, int(round(yi)) + 4)  # clear of the neighbours
            weights = data[rows, xi] - level
            assert (weights * rows).sum() / weights.sum() == pytest.approx(yi, abs=0.1)
    # the per-fiber throughput shows in the integrated flat flux
    integrated = np.array(
        [aperture_spectrum(data - level, truth.trace_coeffs[i]).sum() for i in range(45)]
    )
    ratio = integrated / np.median(integrated)
    expected = truth.throughput / np.median(truth.throughput)
    np.testing.assert_allclose(ratio, expected, rtol=0.03)


def test_no_cal_arcs_and_determinism(tmp_path, highspec_effective, survey_layout, mg_lines):
    kwargs = dict(seed=3, n_bias=1, n_dark=1, n_flat=1, n_full_arc=1, n_cal_arc=2, n_science=3)
    first = simulate_night(
        highspec_effective, survey_layout, mg_lines, tmp_path / "a", no_cal_arcs=True, **kwargs
    )
    assert first.frames.arc_kind.tolist().count("calfiber") == 0
    assert first.frames.arc_kind.tolist().count("full") == 1
    assert (first.frames.exposure_type == "science").sum() == 3
    assert (first.frames.exposure_type == "background").sum() == 1
    second = simulate_night(
        highspec_effective, survey_layout, mg_lines, tmp_path / "b", no_cal_arcs=True, **kwargs
    )
    pd.testing.assert_frame_equal(first.frames, second.frames)
    np.testing.assert_array_equal(first.trace_coeffs, second.trace_coeffs)
    name = first.frames.file[first.frames.exposure_type == "science"].iloc[-1]
    np.testing.assert_array_equal(
        fits.getdata(tmp_path / "a" / name), fits.getdata(tmp_path / "b" / name)
    )
    # a different seed changes the noise
    third = simulate_night(
        highspec_effective,
        survey_layout,
        mg_lines,
        tmp_path / "c",
        no_cal_arcs=True,
        **{**kwargs, "seed": 4},
    )
    assert not np.array_equal(
        fits.getdata(tmp_path / "a" / name), fits.getdata(tmp_path / "c" / name)
    )
    assert not np.array_equal(first.throughput, third.throughput)


def test_lab3_layout_simulates_arc_fibers_only(sim_night_lab3):
    out_dir, truth = sim_night_lab3
    assert len(truth.layout) == 3 and (truth.layout.role == "arc").all()
    assert truth.layout.target.isna().all()
    assert truth.trace_coeffs.shape == (3, 3)
    assert "FIBPRE" not in fits.getheader(frame_path(out_dir, truth, "science"))
    full = cropped(out_dir / truth.frames.file[truth.frames.arc_kind == "full"].iloc[0], truth)
    level = background_level(
        truth, out_dir / truth.frames.file[truth.frames.arc_kind == "full"].iloc[0]
    )
    for i in range(3):
        assert aperture_spectrum(full, truth.trace_coeffs[i]).max() > level * 9 + 5000.0
    # the quadratic curvature term is within the brief's range
    assert (np.abs(truth.trace_coeffs[:, 0]) <= 2e-6).all()


# ------------------------------------------------------------- fibers that differ


def test_simulated_traces_have_distinct_dispersions(sim_night, highspec_effective):
    """Every trace has its own ``(x_c, B)``: the central wavelength lands within
    ``dispersion_jitter_px`` of the nominal column and the relative jitter of ``B`` moves the
    chip edges by at most as much again; ``line_pixels[i]`` is trace ``i``'s own map."""
    _, truth = sim_night
    geometry = sim_geometry(highspec_effective)
    n_pix = truth.extras["n_pix"]
    jitter = truth.extras["dispersion_jitter_px"]
    assert jitter == 0.3
    assert len(set(np.round(truth.x_c_per_trace, 3))) > 1
    assert len(set(np.round(truth.B_per_trace, 3))) > 1
    assert (truth.x_c, truth.B) == nominal_solution(geometry, n_pix)
    centre = np.array(
        [
            grating_x(geometry.lambda_center, x_c, B, geometry.alpha_rad, geometry)
            for x_c, B in zip(truth.x_c_per_trace, truth.B_per_trace, strict=True)
        ]
    )
    assert np.abs(centre - n_pix / 2).max() <= jitter
    assert np.ptp(centre) > jitter  # the 45 draws spread over most of the range
    assert np.abs(truth.B_per_trace / truth.B - 1.0).max() <= 2 * jitter / n_pix
    for i in (0, 17, 44):
        expected = grating_x(
            truth.lines_nm,
            truth.x_c_per_trace[i],
            truth.B_per_trace[i],
            geometry.alpha_rad,
            geometry,
        )
        np.testing.assert_allclose(truth.line_pixels[i], expected, atol=1e-9)
    # the per-trace maps differ by up to twice the jitter on the chip
    on = truth.lines_on_chip
    spread = np.ptp(truth.line_pixels[:, on], axis=0)
    assert spread.max() <= 4 * jitter and spread.min() > 0.1


def test_per_trace_dispersion_shows_on_the_full_arc(
    tmp_path, highspec_effective, survey_layout, mg_lines
):
    """With a large jitter the brightest line sits on each arc fiber's own truth pixel, and
    the arc fibers' pixels differ by much more than the centroid's error."""
    out_dir = tmp_path / "night"
    truth = simulate_night(
        highspec_effective,
        survey_layout,
        mg_lines,
        out_dir,
        seed=1,
        n_bias=1,
        n_dark=1,
        n_flat=1,
        n_cal_arc=0,
        n_science=0,
        dispersion_jitter_px=3.0,
    )
    full = cropped(out_dir / truth.frames.file[truth.frames.arc_kind == "full"].iloc[0], truth)
    line = brightest_line(truth)
    arc_rows = np.flatnonzero(truth.layout.role == "arc")
    x_true = truth.line_pixels[arc_rows, line]
    assert np.ptp(x_true) > 1.0
    for i, x in zip(arc_rows, x_true, strict=True):
        assert centroid(aperture_spectrum(full, truth.trace_coeffs[i]), x) == pytest.approx(
            x, abs=0.1
        )


def test_drift_varies_linearly_with_trace_row(
    tmp_path, highspec_effective, survey_layout, mg_lines
):
    """The drift at trace ``i`` is the frame's shift plus ``gradient_px (y_i - ny/2) / ny``;
    the gradient ramps in time like the shift, from zero at the full arc to
    ``drift_gradient_px`` at the last frame."""
    out_dir = tmp_path / "night"
    truth = simulate_night(
        highspec_effective,
        survey_layout,
        mg_lines,
        out_dir,
        seed=2,
        n_bias=1,
        n_dark=1,
        n_flat=1,
        n_cal_arc=1,
        n_science=1,
        drift_gradient_px=3.0,
    )
    assert truth.drift_gradient_px == 3.0
    frames = truth.frames
    full = frames.index[frames.arc_kind == "full"][0]
    np.testing.assert_allclose(
        frames.gradient_px, 3.0 * frames.shift_px / frames.shift_px.iloc[-1], atol=1e-12
    )
    assert frames.gradient_px.iloc[full] == 0.0 and frames.gradient_px.iloc[-1] == 3.0
    cal = frames.loc[frames.arc_kind == "calfiber"].iloc[0]
    ny = truth.extras["shape"][0]
    y = np.array([np.polyval(c, truth.extras["n_pix"] / 2) for c in truth.trace_coeffs])
    shifts = truth.trace_shifts(cal.file)
    np.testing.assert_allclose(shifts, cal.shift_px + cal.gradient_px * (y - ny / 2) / ny)
    data = cropped(out_dir / cal.file, truth)
    line = brightest_line(truth)
    arc_rows = np.flatnonzero(truth.layout.role == "arc")
    measured = [
        centroid(aperture_spectrum(data, truth.trace_coeffs[i]), truth.line_pixels[i, line])
        - truth.line_pixels[i, line]
        for i in arc_rows
    ]
    np.testing.assert_allclose(measured, shifts[arc_rows], atol=0.1)
    assert measured[-1] - measured[0] > 1.0  # the bottom and top arc fibers drift apart


def test_swapped_fibers_swap_their_spectra_not_the_truth(
    tmp_path, highspec_effective, survey_layout, mg_lines
):
    """``swap_fibers=(a, b)`` (layout trace ids) lights fiber ``a`` with fiber ``b``'s
    spectrum and vice versa; the truth keeps every trace's own map, so the truth comparison
    flags exactly those two. On a layout whose ids are not its row positions the ids, not
    the rows, are swapped, and ``row_of`` finds a fiber's truth row."""
    kwargs = dict(seed=1, n_bias=1, n_dark=1, n_flat=1, n_cal_arc=0, n_science=0)
    kwargs["dispersion_jitter_px"] = 3.0
    plain = simulate_night(highspec_effective, survey_layout, mg_lines, tmp_path / "a", **kwargs)
    out_dir = tmp_path / "b"
    swapped = simulate_night(
        highspec_effective, survey_layout, mg_lines, out_dir, swap_fibers=(1, 3), **kwargs
    )
    np.testing.assert_array_equal(swapped.line_pixels, plain.line_pixels)
    np.testing.assert_array_equal(swapped.x_c_per_trace, plain.x_c_per_trace)
    assert swapped.extras["swap_fibers"] == [1, 3] and plain.extras["swap_fibers"] is None
    full = cropped(
        out_dir / swapped.frames.file[swapped.frames.arc_kind == "full"].iloc[0], swapped
    )
    line = brightest_line(swapped)
    lit_as = {1: 3, 3: 1, 2: 2, 4: 4}
    for trace, source in lit_as.items():
        spectrum = aperture_spectrum(full, swapped.trace_coeffs[trace])
        x = swapped.line_pixels[source, line]
        assert centroid(spectrum, x) == pytest.approx(x, abs=0.1), trace
    assert abs(swapped.line_pixels[1, line] - swapped.line_pixels[3, line]) > 0.5
    with pytest.raises(ValueError, match="swap_fibers"):
        simulate_night(
            highspec_effective,
            survey_layout,
            mg_lines,
            tmp_path / "c",
            swap_fibers=(1, 1),
            **kwargs,
        )
    # ids that are not row positions: the same two fibers swapped by id give the same frames
    renumbered = survey_layout.assign(trace=survey_layout["trace"] + 100)
    by_id = simulate_night(
        highspec_effective, renumbered, mg_lines, tmp_path / "d", swap_fibers=(101, 103), **kwargs
    )
    assert by_id.row_of(103) == 3 and swapped.row_of(3) == 3
    with pytest.raises(KeyError):
        by_id.row_of(3)
    for name in by_id.frames.file:
        np.testing.assert_array_equal(
            fits.getdata(tmp_path / "d" / name), fits.getdata(out_dir / name), err_msg=name
        )
    with pytest.raises(ValueError, match="swap_fibers"):
        simulate_night(
            highspec_effective, renumbered, mg_lines, tmp_path / "e", swap_fibers=(1, 3), **kwargs
        )
