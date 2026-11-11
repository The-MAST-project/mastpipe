"""Tests for the extraction step: aperture weights against hand values, aperture extraction
against toy frames and simulator truth, continuum fitting against a synthetic spectrum."""

from types import SimpleNamespace

import numpy as np
import pytest

from mastpipe.config.models import CalibParams, ExtractParams
from mastpipe.core.model import Frame, FrameMeta, Spectrum, Trace, TraceSet
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.extract import (
    aperture_weights,
    continuum_subtract_all,
    extract_aperture,
    fit_continuum,
)
from mastpipe.steps.trace import robust_sigma

APERTURE_HALFWIDTH = 4.0
"""Rows on either side of a simulated trace (3 px FWHM image, 8.5 px pitch): the aperture
holds 99.9 % of the profile and stays 4.5 rows clear of the neighbouring fiber's centre."""


# ------------------------------------------------------------------------------- fixtures


def make_trace(
    index: int,
    y: float,
    slope: float = 0.0,
    x_center: float = 0.0,
    role: str | None = None,
    telescope: str | None = None,
    aperture: float | None = None,
    coeffs: np.ndarray | None = None,
) -> Trace:
    if coeffs is None:
        coeffs = np.array([0.0, slope, y - slope * x_center])
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
        telescope=telescope,
        aperture=aperture,
    )


def toy_frame(
    data: np.ndarray,
    variance: np.ndarray | None = None,
    mask: np.ndarray | None = None,
    targets: tuple[tuple[str, str], ...] = (),
) -> Frame:
    meta = FrameMeta(instrument="toy", camera="main", exposure_type="science", targets=targets)
    return Frame(data=data.astype(np.float32), meta=meta, variance=variance, mask=mask)


@pytest.fixture(scope="module")
def night(sim_night, highspec_effective) -> SimpleNamespace:
    """The first simulated science frame, calibrated, with traces built from the truth."""
    out_dir, truth = sim_night
    detector = highspec_effective.detector
    refs = scan_folder(out_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    params = CalibParams()
    masters = [build_master(by_type[kind], kind, params) for kind in ("bias", "dark", "background")]
    calib = build_calibration_set(*masters, detector)
    # the simulator writes its science frames as IMAGETYP 'light'
    science = calibrate_frame(by_type["light"][0], calib, saturation=detector.saturation)
    traces = [
        make_trace(
            i,
            float(np.polyval(coeffs, science.shape[1] // 2)),
            role=row["role"],
            telescope=None if isinstance(row["telescope"], float) else row["telescope"],
            coeffs=coeffs,
        )
        for i, (coeffs, (_, row)) in enumerate(zip(truth.trace_coeffs, truth.layout.iterrows(), strict=True))
    ]
    return SimpleNamespace(truth=truth, science=science, traces=TraceSet(traces=traces, shape=science.shape))


# ------------------------------------------------------------------------- aperture_weights


def test_aperture_weights_round_subpixel_floor():
    # centre 10.3, halfwidth 2: whole rows round(8.3) .. round(12.3)
    rows, weights = aperture_weights(10.3, 2.0, "round")
    assert rows.tolist() == [8, 9, 10, 11, 12]
    assert weights.tolist() == [1.0] * 5
    # floor: floor(8.3) .. floor(12.3)
    rows, weights = aperture_weights(10.3, 2.0, "floor")
    assert rows.tolist() == [8, 9, 10, 11, 12]
    assert weights.tolist() == [1.0] * 5
    # subpixel: the interval of width 2 hw + 1 = 5 centred on 10.3 is [7.8, 12.8); row r
    # spans [r - 0.5, r + 0.5), so row 8 keeps 0.7 of its height and row 13 gets 0.3
    rows, weights = aperture_weights(10.3, 2.0, "subpixel")
    assert rows.tolist() == [8, 9, 10, 11, 12, 13]
    assert weights == pytest.approx([0.7, 1.0, 1.0, 1.0, 1.0, 0.3])
    assert weights.sum() == pytest.approx(5.0)
    # the mirror case rounds the other way and the fractional weights swap
    rows, weights = aperture_weights(10.7, 2.0, "round")
    assert rows.tolist() == [9, 10, 11, 12, 13]
    rows, weights = aperture_weights(10.7, 2.0, "floor")
    assert rows.tolist() == [8, 9, 10, 11, 12]
    rows, weights = aperture_weights(10.7, 2.0, "subpixel")
    assert rows.tolist() == [8, 9, 10, 11, 12, 13]
    assert weights == pytest.approx([0.3, 1.0, 1.0, 1.0, 1.0, 0.7])
    # an integer centre makes the subpixel aperture a plain block of 2 hw + 1 rows
    rows, weights = aperture_weights(10.0, 2.0, "subpixel")
    assert rows.tolist() == [8, 9, 10, 11, 12]
    assert weights == pytest.approx([1.0] * 5)
    # a half-integer halfwidth spreads 2 hw + 1 = 6 rows of weight
    rows, weights = aperture_weights(10.0, 2.5, "subpixel")
    assert rows.tolist() == [7, 8, 9, 10, 11, 12, 13]
    assert weights == pytest.approx([0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5])
    with pytest.raises(ValueError, match="aperture"):
        aperture_weights(10.0, 2.0, "gaussian")


# ------------------------------------------------------------------------- extract_aperture


def test_extract_sum_mean_and_metadata():
    ny, nx = 12, 6
    data = np.arange(ny * nx, dtype=float).reshape(ny, nx)
    variance = np.full((ny, nx), 2.0)
    frame = toy_frame(data, variance, targets=(("*", "star"),))
    traces = TraceSet(
        traces=[
            make_trace(0, 5.0, role="science", telescope="mast01"),
            make_trace(1, 8.0, role="sky", telescope="mast01", aperture=1.0),
        ],
        shape=(ny, nx),
    )
    spectra = extract_aperture(frame, traces, ExtractParams(halfwidth=2.0, aperture="round"))
    assert [s.fiber for s in spectra] == [0, 1]
    science, sky = spectra
    assert (science.role, science.telescope, science.target) == ("science", "mast01", "star")
    assert (sky.role, sky.telescope, sky.target) == ("sky", "mast01", None)
    assert science.pixel.tolist() == list(range(nx))
    # rows 3..7 summed per column, error sqrt(5 * 2)
    assert science.flux == pytest.approx(data[3:8].sum(axis=0))
    assert science.error == pytest.approx(np.full(nx, np.sqrt(10.0)))
    assert not science.mask.any()
    # the trace's own aperture (1.0) wins over the parameter: rows 7..9
    assert sky.flux == pytest.approx(data[7:10].sum(axis=0))
    assert sky.error == pytest.approx(np.full(nx, np.sqrt(6.0)))
    # the mean divides the sum and its error by the weight
    mean = extract_aperture(frame, traces, ExtractParams(halfwidth=2.0, method="mean"))[0]
    assert mean.flux == pytest.approx(data[3:8].mean(axis=0))
    assert mean.error == pytest.approx(np.full(nx, np.sqrt(10.0) / 5.0))
    # a subpixel aperture weights the edge rows
    sub = extract_aperture(
        frame,
        TraceSet(traces=[make_trace(0, 5.3)], shape=(ny, nx)),
        ExtractParams(halfwidth=2.0, aperture="subpixel"),
    )[0]
    expected = 0.7 * data[3] + data[4:8].sum(axis=0) + 0.3 * data[8]
    assert sub.flux == pytest.approx(expected)
    assert sub.error == pytest.approx(np.full(nx, np.sqrt(2.0 * (0.7**2 + 4.0 + 0.3**2))))
    with pytest.raises(ValueError, match="variance"):
        extract_aperture(toy_frame(data), traces, ExtractParams(halfwidth=2.0))


def test_extract_assigns_target_per_telescope():
    ny, nx = 16, 4
    data = np.ones((ny, nx))
    targets = (("*", "target_z"), ("mast01", "target_a"), ("mast02", "target_b"))
    frame = toy_frame(data, np.ones((ny, nx)), targets=targets)
    traces = TraceSet(
        traces=[
            make_trace(0, 3.0, role="science", telescope="mast01"),
            make_trace(1, 6.0, role="sky", telescope="mast01"),
            make_trace(2, 9.0, role="science", telescope="mast02"),
            make_trace(3, 12.0, role="science", telescope="mast03"),
        ],
        shape=(ny, nx),
    )
    spectra = extract_aperture(frame, traces, ExtractParams(halfwidth=1.0))
    assert [s.target for s in spectra] == ["target_a", None, "target_b", "target_z"]
    only_one = toy_frame(data, np.ones((ny, nx)), targets=(("mast01", "target_a"),))
    spectra = extract_aperture(only_one, traces, ExtractParams(halfwidth=1.0))
    assert [s.target for s in spectra] == ["target_a", None, None, None]


def test_extract_recovers_simulated_flux(night):
    truth = night.truth
    params = ExtractParams(halfwidth=APERTURE_HALFWIDTH, aperture="round", min_coverage=0.99)
    spectra = extract_aperture(night.science, night.traces, params)
    assert len(spectra) == len(night.traces) == 45
    ratios = []
    pulls = []
    for spectrum, trace in zip(spectra, night.traces, strict=True):
        assert spectrum.fiber == trace.index
        assert spectrum.role == trace.role
        if trace.role != "science":
            assert spectrum.target is None
            continue
        assert spectrum.target == night.science.meta.target_for(trace.telescope)
        row = truth.layout.loc[truth.layout.trace == trace.index]
        assert spectrum.target == row.target.item()
        row_i = truth.row_of(trace.index)
        expected = truth.throughput[row_i] * truth.science_flux[row_i]
        good = ~spectrum.mask
        assert good.mean() > 0.95
        ratios.append(float(np.median(spectrum.flux[good] / expected[good])))
        pull = (spectrum.flux[good] - expected[good]) / spectrum.error[good]
        pulls.append(robust_sigma(pull))
    assert len(ratios) == 20
    assert np.max(np.abs(np.array(ratios) - 1.0)) < 0.02
    # the propagated error tracks the scatter about the truth. The calibrated variance
    # carries the raw frame's Poisson and read noise but not that of the single background
    # frame subtracted from it, which on a 200 ADU continuum over 9 rows (25 ADU² read
    # noise and 15 ADU of dark per pixel) widens the true scatter by about 1.28; the
    # propagation itself is pinned exactly on the toy frames above
    assert 0.8 < np.median(pulls) < 1.5
    # sky fibers carry the sky spectrum, not the science one
    sky = [s for s in spectra if s.role == "sky"][0]
    assert np.median(sky.flux) < 0.5 * np.median(truth.science_flux[truth.row_of(sky.fiber)])


def test_masked_pixel_masks_column():
    ny, nx = 10, 8
    data = np.ones((ny, nx))
    variance = np.ones((ny, nx))
    mask = np.zeros((ny, nx), dtype=bool)
    mask[4, 3] = True  # inside the aperture of a trace at row 5 (rows 3..7)
    mask[0, 5] = True  # outside it
    frame = toy_frame(data, variance, mask)
    traces = TraceSet(traces=[make_trace(0, 5.0)], shape=(ny, nx))
    spectrum = extract_aperture(frame, traces, ExtractParams(halfwidth=2.0))[0]
    assert spectrum.mask.tolist() == [i == 3 for i in range(nx)]
    # coverage 4 / 5 is below the default 0.99: the column is NaN
    assert np.isnan(spectrum.flux[3]) and np.isnan(spectrum.error[3])
    assert spectrum.flux[5] == pytest.approx(5.0)
    # with a lower coverage floor the column keeps its partial sum but stays masked
    partial = extract_aperture(frame, traces, ExtractParams(halfwidth=2.0, min_coverage=0.5))[0]
    assert partial.mask.tolist() == [i == 3 for i in range(nx)]
    assert partial.flux[3] == pytest.approx(4.0)
    assert partial.error[3] == pytest.approx(2.0)


def test_coverage_below_min_is_masked():
    ny, nx = 10, 6
    data = np.ones((ny, nx))
    variance = np.ones((ny, nx))
    data[6, 2] = np.nan  # a non-finite pixel counts as uncovered
    frame = toy_frame(data, variance)
    # a trace one row from the top edge: rows -1 .. 3, four of five on the detector
    traces = TraceSet(traces=[make_trace(0, 1.0)], shape=(ny, nx))
    strict = extract_aperture(frame, traces, ExtractParams(halfwidth=2.0, min_coverage=0.99))[0]
    assert strict.mask.all()
    assert np.isnan(strict.flux).all()
    loose = extract_aperture(frame, traces, ExtractParams(halfwidth=2.0, min_coverage=0.8))[0]
    assert not loose.mask.any()
    assert loose.flux == pytest.approx(np.full(nx, 4.0))
    assert loose.error == pytest.approx(np.full(nx, 2.0))
    # the NaN pixel on a trace at row 6 drops that column's coverage to 4 / 5
    inner = extract_aperture(
        frame,
        TraceSet(traces=[make_trace(0, 6.0)], shape=(ny, nx)),
        ExtractParams(halfwidth=2.0, min_coverage=0.99),
    )[0]
    assert inner.mask.tolist() == [i == 2 for i in range(nx)]
    assert np.isnan(inner.flux[2])
    assert inner.flux[0] == pytest.approx(5.0)
    # at the bottom edge (rows 6..10, row 10 off the detector) the aperture does not wrap
    # around; column 2 also loses the NaN pixel at row 6 and keeps three of five rows
    bottom = extract_aperture(
        frame,
        TraceSet(traces=[make_trace(0, 8.0)], shape=(ny, nx)),
        ExtractParams(halfwidth=2.0, min_coverage=0.5),
    )[0]
    assert bottom.flux == pytest.approx([4.0, 4.0, 3.0, 4.0, 4.0, 4.0])
    assert not bottom.mask.any()


# ---------------------------------------------------------------------------- fit_continuum


def synthetic_spectrum(seed: int = 3) -> tuple[Spectrum, np.ndarray, np.ndarray]:
    """A smooth continuum with noise, strong emission lines, a cold pixel, a masked pixel
    and a NaN pixel; returns the spectrum, the true continuum and the line-free selection."""
    rng = np.random.default_rng(seed)
    n = 2000
    pixel = np.arange(n)
    x = pixel / (n - 1)
    continuum = 100.0 + 30.0 * np.cos(2.0 * np.pi * 1.3 * x) + 20.0 * x
    centres = np.linspace(60.0, n - 60.0, 15)
    peaks = rng.uniform(500.0, 5000.0, size=centres.size)
    lines = (peaks[:, None] * np.exp(-0.5 * ((pixel[None, :] - centres[:, None]) / 1.3) ** 2)).sum(axis=0)
    flux = continuum + lines + rng.normal(0.0, 1.0, size=n)
    flux[700] = -2000.0  # a cold pixel
    flux[1300] = np.nan
    mask = np.zeros(n, dtype=bool)
    mask[900] = True
    flux[900] = 1e5  # masked, so the fit never sees it
    error = np.ones(n)
    away = np.min(np.abs(pixel[:, None] - centres[None, :]), axis=1) > 8.0
    away &= ~mask & np.isfinite(flux)
    away[700] = False
    return Spectrum(fiber=0, pixel=pixel, flux=flux, error=error, mask=mask), continuum, away


def test_continuum_fit_ignores_lines():
    spectrum, continuum, away = synthetic_spectrum()
    flux_before = spectrum.flux.copy()
    fitted = fit_continuum(spectrum, degree=10, clip_above=2.0, clip_below=3.0, n_iter=10)
    assert fitted is not spectrum
    assert fitted.continuum is not None and fitted.continuum.shape == spectrum.flux.shape
    # the raw flux is untouched, on both the input and the output
    np.testing.assert_array_equal(spectrum.flux, flux_before)
    np.testing.assert_array_equal(fitted.flux, flux_before)
    assert fitted.mask is not spectrum.mask or fitted.mask.tolist() == spectrum.mask.tolist()
    # the continuum is evaluated everywhere, within 1 % of the truth away from the lines
    assert np.isfinite(fitted.continuum).all()
    relative = np.abs(fitted.continuum[away] / continuum[away] - 1.0)
    assert relative.max() < 0.01
    # under a line the fit interpolates the continuum rather than following the line
    peak = int(np.argmax(spectrum.flux[:600]))
    assert abs(fitted.continuum[peak] / continuum[peak] - 1.0) < 0.02
    # the continuum-subtracted flux is the caller's arithmetic
    subtracted = fitted.flux - fitted.continuum
    assert abs(np.median(subtracted[away])) < 0.5


def test_continuum_subtract_all_uses_params():
    spectrum, continuum, away = synthetic_spectrum(seed=5)
    params = ExtractParams(continuum_degree=12, clip_above=2.0, clip_below=3.0, continuum_iter=8)
    out = continuum_subtract_all([spectrum, spectrum], params)
    assert len(out) == 2
    for fitted in out:
        assert fitted.continuum is not None
        assert np.abs(fitted.continuum[away] / continuum[away] - 1.0).max() < 0.01
    assert spectrum.continuum is None
    # a spectrum with too few usable pixels gets a NaN continuum rather than an error
    dead = Spectrum(
        fiber=1,
        pixel=np.arange(50),
        flux=np.full(50, np.nan),
        error=np.ones(50),
        mask=np.zeros(50, dtype=bool),
    )
    (fitted,) = continuum_subtract_all([dead], params)
    assert fitted.continuum is not None and np.isnan(fitted.continuum).all()
