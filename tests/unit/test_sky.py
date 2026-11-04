"""Tests for the sky-subtraction step: resampling on toy spectra, the kept raw extraction,
the pairing of science and sky fibers per telescope, and the simulated night's sky lines."""

from types import SimpleNamespace

import numpy as np
import pytest

from mastpipe.config.models import CalibParams, ExtractParams, SkyParams
from mastpipe.core.model import Frame, Spectrum, Trace, TraceSet
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.sim.simulator import undrift_x
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.extract import extract_aperture
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.sky import subtract_sky, subtract_sky_all
from mastpipe.steps.wavemodel import GratingModel
from tests.conftest import truth_traceset

N_PIX = 200
LINE_PX = 80.0
LINE_SIGMA_PX = 1.3
LINE_PEAK = 300.0
CONTINUUM = 50.0
APERTURE_HALFWIDTH = 4.0
"""Rows on either side of a simulated trace (3 px FWHM image, 8.5 px pitch)."""


# ------------------------------------------------------------------------------- fixtures


def make_trace(index: int, role: str, telescope: str | None) -> Trace:
    return Trace(
        index=index,
        coeffs=np.array([0.0, 0.0, 10.0 * index]),
        y_center=10.0 * index,
        slope=0.0,
        rms=0.0,
        n_points=10,
        snr=100.0,
        x_samples=np.arange(10.0),
        y_samples=np.full(10, 10.0 * index),
        kept=np.ones(10, dtype=bool),
        role=role,
        telescope=telescope,
    )


def make_spectrum(
    fiber: int,
    flux: np.ndarray,
    role: str | None = None,
    telescope: str | None = None,
    wavelength: np.ndarray | None = None,
    mask: np.ndarray | None = None,
    error: float = 2.0,
) -> Spectrum:
    n = len(flux)
    return Spectrum(
        fiber=fiber,
        pixel=np.arange(n),
        flux=np.asarray(flux, dtype=float),
        error=np.full(n, error),
        mask=np.zeros(n, dtype=bool) if mask is None else mask,
        role=role,
        telescope=telescope,
        wavelength=wavelength,
    )


def gaussian_line(x: np.ndarray, centre: float, peak: float = LINE_PEAK) -> np.ndarray:
    return CONTINUUM + peak * np.exp(-0.5 * ((x - centre) / LINE_SIGMA_PX) ** 2)


def grid(offset_px: float = 0.0, dispersion: float = 0.01, origin: float = 515.0) -> np.ndarray:
    """A linear wavelength grid; ``offset_px`` shifts it by that many pixels."""
    return origin + dispersion * (np.arange(N_PIX) + offset_px)


@pytest.fixture(scope="module")
def night(sim_night, highspec_effective) -> SimpleNamespace:
    """The first simulated science frame extracted, every fiber given its truth wavelengths.

    ``spectra`` (one per trace, with the layout's roles and telescopes), ``traces``,
    ``truth``, ``throughput`` (the truth throughput keyed by fiber), ``x_of`` (pixel of a
    wavelength on this frame's drifted grid) and ``line_px`` (the sky lines' pixels).
    """
    out_dir, truth = sim_night
    detector = highspec_effective.detector
    refs = scan_folder(out_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    params = CalibParams()
    masters = [build_master(by_type[kind], kind, params) for kind in ("bias", "dark", "background")]
    calib = build_calibration_set(*masters, detector)
    frame = calibrate_frame(by_type["science"][0], calib, saturation=detector.saturation)
    traces = truth_traceset(truth, frame.shape)
    spectra = extract_aperture(frame, traces, ExtractParams(halfwidth=APERTURE_HALFWIDTH))

    row = truth.frames.set_index("file").loc[frame.meta.source]
    shift, stretch = float(row["shift_px"]), float(row["stretch"])
    geom = grating_geometry(highspec_effective.geometry, detector.pixel_size_um)
    n_pix = frame.shape[1]
    model = GratingModel(truth.x_c, truth.B, geom.alpha_rad, None, n_pix, geom)
    x = np.arange(n_pix, dtype=float)
    wavelength = np.asarray(model.wavelength(undrift_x(x, shift, stretch, n_pix)), dtype=float)
    for spectrum in spectra:
        spectrum.wavelength = wavelength.copy()
    increasing = wavelength[0] < wavelength[-1]
    xp = wavelength if increasing else wavelength[::-1]
    fp = x if increasing else x[::-1]
    return SimpleNamespace(
        truth=truth,
        traces=traces,
        spectra=spectra,
        throughput={i: float(t) for i, t in enumerate(truth.throughput)},
        x_of=lambda nm: np.interp(nm, xp, fp),
        line_px=np.interp(truth.sky_lines_nm, xp, fp),
    )


# ----------------------------------------------------------------------------- subtract_sky


def test_subtract_sky_resamples_scaled_sky_onto_science_grid():
    # the same line seen by two fibers whose grids differ by 0.4 px and whose throughputs
    # differ by a factor 1.5: in wavelength the line sits at the same place
    x = np.arange(N_PIX, dtype=float)
    science_wl = grid(0.0)
    sky_wl = grid(0.4)
    science = make_spectrum(
        1, 1.5 * gaussian_line(x, LINE_PX), "science", "mast01", wavelength=science_wl
    )
    sky = make_spectrum(2, gaussian_line(x, LINE_PX - 0.4), "sky", "mast01", wavelength=sky_wl)
    out = subtract_sky(science, sky, 1.5, SkyParams())
    assert out is not science
    skysub = out.versions["skysub"]
    inside = slice(1, None)  # the first science pixel lies before the shifted sky grid
    # linear interpolation of a 3 px FWHM line between two samples 0.4 px apart leaves at
    # most ~7 % of the peak at the line centre (the second derivative there), nothing in
    # the continuum
    assert np.all(np.abs(skysub.flux[inside]) < 0.1 * 1.5 * LINE_PEAK)
    assert np.all(np.abs(skysub.flux[inside][:70]) < 1e-9)
    assert np.allclose(skysub.flux[inside], (science.flux - out.sky)[inside])
    # the resampled sky peaks where the science line is, not where the sky fiber saw it
    assert int(np.nanargmax(out.sky)) == int(LINE_PX)
    # naive pixel-by-pixel subtraction would leave the shifted line's derivative behind
    naive = science.flux - 1.5 * sky.flux
    assert np.max(np.abs(naive)) > 0.15 * 1.5 * LINE_PEAK
    assert np.max(np.abs(skysub.flux[inside])) < 0.5 * np.max(np.abs(naive))
    # errors add in quadrature: science pixel i lies 0.6 of the way from sky pixel i - 1 to
    # sky pixel i, so the sky term is 1.5 * 2 * sqrt(0.4^2 + 0.6^2); the science pixel
    # outside the sky coverage is masked
    sky_term = 1.5 * 2.0 * np.sqrt(0.4**2 + 0.6**2)
    assert np.allclose(skysub.error[inside], np.sqrt(science.error[inside] ** 2 + sky_term**2))
    assert skysub.mask[0] and not skysub.mask[inside].any()
    assert np.isnan(out.sky[0]) and np.isnan(skysub.flux[0])
    # a decreasing sky grid is interpolated the same way
    sky_flipped = make_spectrum(
        2, sky.flux[::-1], "sky", "mast01", wavelength=sky_wl[::-1], mask=sky.mask[::-1]
    )
    flipped = subtract_sky(science, sky_flipped, 1.5, SkyParams())
    assert np.allclose(flipped.sky, out.sky, equal_nan=True)


def test_sky_mask_marks_any_masked_neighbour():
    science = make_spectrum(1, np.full(N_PIX, CONTINUUM), "science", "mast01", wavelength=grid(0.0))
    sky_mask = np.zeros(N_PIX, dtype=bool)
    sky_mask[100] = True
    science_mask = np.zeros(N_PIX, dtype=bool)
    science_mask[5] = True
    science.mask = science_mask
    sky = make_spectrum(
        2, np.full(N_PIX, CONTINUUM), "sky", "mast01", wavelength=grid(0.5), mask=sky_mask
    )
    out = subtract_sky(science, sky, 1.0, SkyParams())
    skysub = out.versions["skysub"]
    # science pixel 100 falls between sky pixels 99 and 100, pixel 101 between 100 and 101
    assert skysub.mask[[100, 101]].all()
    assert not skysub.mask[[98, 99, 102, 103]].any()
    assert skysub.mask[5]  # the science mask is kept
    assert skysub.mask[0]  # before the first sky wavelength: outside the coverage
    assert np.array_equal(out.mask, science_mask)  # raw mask untouched


def test_sky_error_propagates_exactly():
    # a half-pixel offset: every science pixel i is the mean of sky pixels i - 1 and i, so
    # its sky error is sqrt(0.25 s[i-1]^2 + 0.25 s[i]^2), not the mean of the two sigmas
    x = np.arange(N_PIX, dtype=float)
    sigma = 1.0 + 0.5 * np.sin(x / 7.0) + (x % 3)  # varies pixel to pixel
    science = make_spectrum(
        1, np.zeros(N_PIX), "science", "mast01", wavelength=grid(0.0), error=0.0
    )
    sky = make_spectrum(2, np.zeros(N_PIX), "sky", "mast01", wavelength=grid(0.5))
    sky.error = sigma.copy()
    out = subtract_sky(science, sky, 2.0, SkyParams())
    error = out.versions["skysub"].error
    expected = 2.0 * np.sqrt(0.25 * sigma[:-1] ** 2 + 0.25 * sigma[1:] ** 2)
    assert np.allclose(error[1:], expected)
    assert not np.allclose(error[1:], 2.0 * 0.5 * (sigma[:-1] + sigma[1:]))  # not interpolated
    assert np.isnan(error[0])
    # a sky grid 0.25 px later puts science pixel i at t = 0.75 between sky pixels i - 1
    # and i: weights 0.25 and 0.75
    sky.wavelength = grid(0.25)
    error = subtract_sky(science, sky, 1.0, SkyParams()).versions["skysub"].error
    assert np.allclose(error[1:], np.sqrt(0.25**2 * sigma[:-1] ** 2 + 0.75**2 * sigma[1:] ** 2))
    # on the same grid every target sits on a sky pixel and takes its sigma
    sky.wavelength = grid(0.0)
    error = subtract_sky(science, sky, 1.0, SkyParams()).versions["skysub"].error
    assert np.array_equal(error, sigma)
    # the science error adds in quadrature
    science.error = np.full(N_PIX, 3.0)
    error = subtract_sky(science, sky, 1.0, SkyParams()).versions["skysub"].error
    assert np.allclose(error, np.sqrt(9.0 + sigma**2))


def test_skysub_version_and_raw_kept():
    x = np.arange(N_PIX, dtype=float)
    flux = gaussian_line(x, LINE_PX)
    science = make_spectrum(1, flux, "science", "mast01", wavelength=grid(0.0))
    sky = make_spectrum(2, gaussian_line(x, LINE_PX, 0.5 * LINE_PEAK), "sky", "mast01", grid(0.0))
    raw_flux, raw_error, raw_mask = science.flux.copy(), science.error.copy(), science.mask.copy()
    out = subtract_sky(science, sky, 2.0, SkyParams())
    assert "skysub" in out.versions and out.latest() is out.versions["skysub"]
    assert np.array_equal(out.flux, raw_flux)
    assert np.array_equal(out.error, raw_error)
    assert np.array_equal(out.mask, raw_mask)
    assert np.array_equal(science.flux, raw_flux) and science.versions == {}
    assert out.sky is not None
    assert np.allclose(out.sky, 2.0 * sky.flux)
    assert np.allclose(out.versions["skysub"].flux, raw_flux - 2.0 * sky.flux)
    assert np.allclose(out.versions["skysub"].error, np.sqrt(raw_error**2 + (2.0 * sky.error) ** 2))
    assert (out.fiber, out.role, out.telescope) == (1, "science", "mast01")
    assert np.array_equal(out.wavelength, science.wavelength)
    # NaN flux (a low-coverage column) does not warn and stays masked
    science.flux[7] = np.nan
    science.mask[7] = True
    again = subtract_sky(science, sky, 2.0, SkyParams())
    assert np.isnan(again.versions["skysub"].flux[7]) and again.versions["skysub"].mask[7]


def test_requires_wavelength():
    x = np.arange(N_PIX, dtype=float)
    science = make_spectrum(1, gaussian_line(x, LINE_PX), "science", "mast01", wavelength=grid())
    sky = make_spectrum(2, gaussian_line(x, LINE_PX), "sky", "mast01", wavelength=grid())
    with pytest.raises(ValueError, match="wavelength"):
        subtract_sky(make_spectrum(1, science.flux, "science", "mast01"), sky, 1.0, SkyParams())
    with pytest.raises(ValueError, match="wavelength"):
        subtract_sky(science, make_spectrum(2, sky.flux, "sky", "mast01"), 1.0, SkyParams())
    with pytest.raises(ValueError, match="monotonic"):
        wobbly = grid()
        wobbly[50], wobbly[51] = wobbly[51], wobbly[50]
        subtract_sky(
            science,
            make_spectrum(2, sky.flux, "sky", "mast01", wavelength=wobbly),
            1.0,
            SkyParams(),
        )


# ------------------------------------------------------------------------- subtract_sky_all


def toy_layout() -> tuple[TraceSet, list[Spectrum]]:
    """Two telescopes with science + sky, one arc fiber, one science fiber without a sky."""
    x = np.arange(N_PIX, dtype=float)
    traces = TraceSet(
        traces=[
            make_trace(0, "arc", None),
            make_trace(1, "science", "mast01"),
            make_trace(2, "sky", "mast01"),
            make_trace(3, "science", "mast02"),
            make_trace(4, "sky", "mast02"),
            make_trace(5, "science", "mast03"),
        ],
        shape=(100, N_PIX),
    )
    spectra = [
        make_spectrum(0, np.zeros(N_PIX), "arc", None, wavelength=grid()),
        make_spectrum(1, gaussian_line(x, LINE_PX), "science", "mast01", wavelength=grid()),
        make_spectrum(2, gaussian_line(x, LINE_PX), "sky", "mast01", wavelength=grid()),
        make_spectrum(3, gaussian_line(x, LINE_PX), "science", "mast02", wavelength=grid()),
        make_spectrum(4, gaussian_line(x, LINE_PX), "sky", "mast02", wavelength=grid()),
        make_spectrum(5, gaussian_line(x, LINE_PX), "science", "mast03", wavelength=grid()),
    ]
    return traces, spectra


def test_subtract_sky_all_pairs_by_telescope_and_scales():
    traces, spectra = toy_layout()
    throughput = {1: 1.2, 2: 0.8, 3: 1.0, 4: 1.0}
    out, skipped = subtract_sky_all(spectra, traces, throughput, SkyParams())
    assert [s.fiber for s in out] == [0, 1, 2, 3, 4, 5]
    assert out[0] is spectra[0] and out[2] is spectra[2] and out[4] is spectra[4]
    assert "skysub" in out[1].versions and "skysub" in out[3].versions
    assert np.allclose(out[1].sky, 1.5 * spectra[2].flux)
    assert np.allclose(out[3].sky, spectra[4].flux)
    assert skipped == [{"fiber": 5, "telescope": "mast03", "reason": "no sky fiber"}]
    assert out[5] is spectra[5]
    # no throughput for one side, or scaling switched off: ratio 1
    out, _ = subtract_sky_all(spectra, traces, {1: 1.2}, SkyParams())
    assert np.allclose(out[1].sky, spectra[2].flux)
    out, _ = subtract_sky_all(spectra, traces, throughput, SkyParams(scale_by_throughput=False))
    assert np.allclose(out[1].sky, spectra[2].flux)


def test_missing_sky_fiber_leaves_raw_only():
    traces, spectra = toy_layout()
    # the sky trace of mast01 is absent from the spectra (extraction skipped it)
    without = [s for s in spectra if s.fiber != 2]
    out, skipped = subtract_sky_all(without, traces, {}, SkyParams())
    assert [s.fiber for s in out] == [0, 1, 3, 4, 5]
    assert out[1] is without[1] and out[1].versions == {} and out[1].sky is None
    assert skipped == [
        {"fiber": 1, "telescope": "mast01", "reason": "no sky fiber"},
        {"fiber": 5, "telescope": "mast03", "reason": "no sky fiber"},
    ]
    assert "skysub" in out[2].versions
    # a science fiber without a telescope has no sky fiber either
    loose = [make_spectrum(9, np.zeros(N_PIX), "science", None, wavelength=grid())]
    out, skipped = subtract_sky_all(loose, traces, {}, SkyParams())
    assert out == loose
    assert skipped == [{"fiber": 9, "telescope": None, "reason": "no sky fiber"}]


def test_masked_sky_fiber_leaves_raw_only():
    traces, spectra = toy_layout()
    spectra[2].mask[:] = True
    out, skipped = subtract_sky_all(spectra, traces, {}, SkyParams())
    assert out[1] is spectra[1] and out[1].versions == {}
    assert {"fiber": 1, "telescope": "mast01", "reason": "sky fiber fully masked"} in skipped
    assert "skysub" in out[3].versions
    # a sky fiber without wavelength, and a science fiber without one
    traces, spectra = toy_layout()
    spectra[4].wavelength = None
    spectra[1].wavelength = None
    out, skipped = subtract_sky_all(spectra, traces, {}, SkyParams())
    assert out[1] is spectra[1] and out[3] is spectra[3]
    assert skipped == [
        {"fiber": 1, "telescope": "mast01", "reason": "no wavelength"},
        {"fiber": 3, "telescope": "mast02", "reason": "no wavelength"},
        {"fiber": 5, "telescope": "mast03", "reason": "no sky fiber"},
    ]


def test_bad_throughput_leaves_raw_only():
    traces, spectra = toy_layout()
    throughput = {1: np.nan, 2: 0.8, 3: 1.0, 4: 0.0}
    out, skipped = subtract_sky_all(spectra, traces, throughput, SkyParams())
    assert out[1] is spectra[1] and out[3] is spectra[3]
    assert out[1].versions == {} and out[3].versions == {}
    assert skipped == [
        {"fiber": 1, "telescope": "mast01", "reason": "bad throughput"},
        {"fiber": 3, "telescope": "mast02", "reason": "bad throughput"},
        {"fiber": 5, "telescope": "mast03", "reason": "no sky fiber"},
    ]
    out, skipped = subtract_sky_all(spectra, traces, {1: 1.0, 2: -0.5}, SkyParams())
    assert skipped[0] == {"fiber": 1, "telescope": "mast01", "reason": "bad throughput"}
    assert "skysub" in out[3].versions  # mast02 has no throughput at all: ratio 1
    # with scaling off the throughputs are not consulted
    out, skipped = subtract_sky_all(
        spectra, traces, throughput, SkyParams(scale_by_throughput=False)
    )
    assert [s["reason"] for s in skipped] == ["no sky fiber"]
    assert "skysub" in out[1].versions and "skysub" in out[3].versions
    assert np.all(np.isfinite(out[1].sky))


# -------------------------------------------------------------------------- simulated night


def test_subtract_sky_removes_simulated_sky_lines(night):
    """The simulator's science fibers carry no sky, so the sky fiber of one telescope is
    taken as the "science" and the sky fiber of the next telescope as its sky: both saw the
    same sky through different throughputs and with independent noise. Scaled by the truth
    throughput ratio, the subtraction leaves no line: the residual at each of the
    ``sky_lines_nm`` is below 3 sigma of the propagated error, while the unsubtracted
    "science" shows the lines at more than 10 sigma and an unscaled subtraction leaves a
    residual above 3 sigma for at least one pair."""
    skies = [s for s in night.spectra if s.role == "sky"]
    assert len(skies) >= 2
    at_lines = np.rint(night.line_px).astype(int)
    assert np.all((at_lines > 2) & (at_lines < len(skies[0].pixel) - 2))
    unscaled_fails = 0
    for science, sky in zip(skies, skies[1:], strict=False):
        ratio = night.throughput[science.fiber] / night.throughput[sky.fiber]
        assert not (science.mask[at_lines].any() or sky.mask[at_lines].any())
        continuum = np.nanmedian(science.flux[~science.mask])
        assert np.all((science.flux[at_lines] - continuum) > 10.0 * science.error[at_lines])
        out = subtract_sky(science, sky, ratio, SkyParams())
        skysub = out.versions["skysub"]
        assert not skysub.mask[at_lines].any()
        assert np.all(np.abs(skysub.flux[at_lines]) < 3.0 * skysub.error[at_lines])
        # the subtracted sky carries the lines
        assert np.all(out.sky[at_lines] > 10.0 * science.error[at_lines])
        plain = subtract_sky(science, sky, 1.0, SkyParams()).versions["skysub"]
        if np.any(np.abs(plain.flux[at_lines]) > 3.0 * plain.error[at_lines]):
            unscaled_fails += 1
    assert unscaled_fails > 0


def test_subtract_sky_all_on_simulated_frame(night):
    out, skipped = subtract_sky_all(night.spectra, night.traces, night.throughput, SkyParams())
    assert skipped == []
    assert [s.fiber for s in out] == [s.fiber for s in night.spectra]
    for before, after in zip(night.spectra, out, strict=True):
        if before.role != "science":
            assert after is before
            continue
        sky = night.traces.by_telescope(before.telescope)["sky"]
        sky_spectrum = next(s for s in night.spectra if s.fiber == sky.index)
        ratio = night.throughput[before.fiber] / night.throughput[sky.index]
        assert np.array_equal(after.flux, before.flux)
        assert np.allclose(after.sky, ratio * sky_spectrum.flux, equal_nan=True)
        skysub = after.versions["skysub"]
        good = ~skysub.mask
        assert np.allclose(skysub.flux[good], before.flux[good] - after.sky[good])
        # the science fiber had no sky: the subtraction digs the sky lines in
        at_lines = np.rint(night.line_px).astype(int)
        assert np.all(skysub.flux[at_lines] < -10.0 * skysub.error[at_lines])
