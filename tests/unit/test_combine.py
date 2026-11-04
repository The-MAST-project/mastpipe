"""Tests for the within-exposure combination step: the common grid, the weighted mean on
toy spectra, the include/exclude and rejection bookkeeping, versions combined separately
and the grouping of science fibers by target."""

import numpy as np
import pytest

from mastpipe.config.models import CombineParams
from mastpipe.core.model import CombinedSpectrum, Spectrum, SpectrumVersion
from mastpipe.steps.combine import combine_spectra, common_grid, group_by_target
from mastpipe.vocab import unit_name

N_PIX = 200
CONTINUUM = 100.0
SIGMA = 5.0


def grid(dispersion: float = 0.01, origin: float = 515.0, n: int = N_PIX) -> np.ndarray:
    return origin + dispersion * np.arange(n)


def make_spectrum(
    fiber: int,
    flux: np.ndarray,
    wavelength: np.ndarray | None,
    role: str | None = "science",
    telescope: str | None = None,
    target: str | None = "star",
    error: float = SIGMA,
    mask: np.ndarray | None = None,
    versions: dict[str, SpectrumVersion] | None = None,
) -> Spectrum:
    n = len(flux)
    return Spectrum(
        fiber=fiber,
        pixel=np.arange(n),
        flux=np.asarray(flux, dtype=float),
        error=np.full(n, float(error)),
        mask=np.zeros(n, dtype=bool) if mask is None else mask,
        role=role,
        telescope=telescope,
        target=target,
        wavelength=wavelength,
        versions={} if versions is None else versions,
    )


def noisy_copies(rng: np.random.Generator, n_copies: int = 4) -> list[Spectrum]:
    wave = grid()
    return [
        make_spectrum(
            fiber=i,
            flux=CONTINUUM + rng.normal(0.0, SIGMA, N_PIX),
            wavelength=wave,
            telescope=unit_name(i + 1),
        )
        for i in range(n_copies)
    ]


# ---------------------------------------------------------------------------- common grid


def test_common_grid_log_step():
    coarse = make_spectrum(0, np.ones(N_PIX), grid(dispersion=0.02, origin=515.0))
    fine = make_spectrum(1, np.ones(N_PIX), grid(dispersion=0.01, origin=516.0))
    out = common_grid([coarse, fine], log=True)
    assert out[0] == pytest.approx(515.0)
    assert out[-1] == pytest.approx(515.0 + 0.02 * (N_PIX - 1))
    log_steps = np.diff(np.log(out))
    assert np.allclose(log_steps, log_steps[0])
    finest = np.min(np.diff(fine.wavelength) / fine.wavelength[:-1])
    assert log_steps[0] <= finest
    assert log_steps[0] > 0.9 * finest

    linear = common_grid([coarse, fine], log=False)
    steps = np.diff(linear)
    assert np.allclose(steps, steps[0])
    assert steps[0] <= 0.01 + 1e-12
    assert linear[0] == pytest.approx(515.0)


# -------------------------------------------------------------------------- weighted mean


def test_combine_weighted_mean_reduces_noise():
    rng = np.random.default_rng(1)
    spectra = noisy_copies(rng)
    cs = combine_spectra(spectra, "star", CombineParams(log_grid=False))
    assert isinstance(cs, CombinedSpectrum)
    assert cs.target == "star"
    assert cs.grouping == "exposure"
    assert set(cs.versions) == {"raw"}
    raw = cs.versions["raw"]
    assert not raw.mask.any()
    assert np.std(raw.flux) == pytest.approx(SIGMA / 2, rel=0.15)
    assert np.median(raw.error) == pytest.approx(SIGMA / 2, rel=0.05)
    assert list(cs.contributors["fiber"]) == [0, 1, 2, 3]
    assert np.allclose(cs.contributors["weight"], 0.25)
    assert list(cs.contributors["source"]) == [""] * 4
    assert len(cs.excluded) == 0
    assert cs.barycentric == {}


def test_outside_coverage_masked():
    rng = np.random.default_rng(2)
    short = make_spectrum(0, CONTINUUM + rng.normal(0, SIGMA, 100), grid(n=100), telescope="mast01")
    long = make_spectrum(1, CONTINUUM + rng.normal(0, SIGMA, N_PIX), grid(), telescope="mast02")
    cs = combine_spectra([short, long], "star", CombineParams(log_grid=False))
    raw = cs.versions["raw"]
    assert not raw.mask.any()
    covered = cs.wavelength <= short.wavelength[-1]
    assert np.median(raw.error[covered]) < np.median(raw.error[~covered])
    assert np.median(raw.error[~covered]) == pytest.approx(SIGMA, rel=0.05)


def test_all_masked_pixel_is_masked_in_combination():
    rng = np.random.default_rng(3)
    spectra = noisy_copies(rng, n_copies=2)
    for spectrum in spectra:
        spectrum.mask[50] = True
    cs = combine_spectra(spectra, "star", CombineParams(log_grid=False))
    raw = cs.versions["raw"]
    assert raw.mask[50]
    assert raw.mask.sum() == 1


# --------------------------------------------------------------------- include / exclude


def test_exclude_list_and_reason():
    rng = np.random.default_rng(4)
    spectra = noisy_copies(rng)
    cs = combine_spectra(spectra, "star", CombineParams(exclude=["mast03"], log_grid=False))
    assert list(cs.contributors["telescope"]) == ["mast01", "mast02", "mast04"]
    assert list(cs.excluded["telescope"]) == ["mast03"]
    assert list(cs.excluded["fiber"]) == [2]
    assert list(cs.excluded["reason"]) == ["excluded by configuration"]

    cs = combine_spectra(
        spectra, "star", CombineParams(include=["mast01", "mast02"], log_grid=False)
    )
    assert list(cs.contributors["telescope"]) == ["mast01", "mast02"]
    assert sorted(cs.excluded["telescope"]) == ["mast03", "mast04"]
    assert set(cs.excluded["reason"]) == {"excluded by configuration"}


def test_rejection_low_snr_listed_not_dropped_silently():
    rng = np.random.default_rng(5)
    spectra = noisy_copies(rng)
    spectra[1] = make_spectrum(
        1, CONTINUUM + rng.normal(0, 60.0, N_PIX), grid(), telescope="mast02", error=60.0
    )
    params = CombineParams(reject=True, reject_snr_min=3.0, log_grid=False)
    cs = combine_spectra(spectra, "star", params)
    assert list(cs.contributors["telescope"]) == ["mast01", "mast03", "mast04"]
    assert list(cs.excluded["fiber"]) == [1]
    reason = cs.excluded["reason"].iloc[0]
    assert reason.startswith("low S/N (")
    assert "< 3.0" in reason
    # Without rejection the noisy spectrum contributes, with a small weight.
    kept = combine_spectra(spectra, "star", CombineParams(log_grid=False))
    assert list(kept.contributors["fiber"]) == [0, 1, 2, 3]
    assert kept.contributors["weight"].iloc[1] < 0.05


def test_rejection_flux_ratio_and_drift():
    rng = np.random.default_rng(6)
    spectra = noisy_copies(rng)
    spectra[2] = make_spectrum(
        2, 0.1 * CONTINUUM + rng.normal(0, 0.5, N_PIX), grid(), telescope="mast03", error=0.5
    )
    qc = {0: {"drift_px": 0.1}, 1: {"drift_px": 0.1}, 2: {"drift_px": 0.1}, 3: {"drift_px": 2.5}}
    params = CombineParams(reject=True, reject_drift_px=1.0, log_grid=False)
    cs = combine_spectra(spectra, "star", params, qc=qc)
    assert list(cs.contributors["fiber"]) == [0, 1]
    reasons = dict(zip(cs.excluded["fiber"], cs.excluded["reason"], strict=True))
    assert reasons[2].startswith("flux ratio ")
    assert "below" in reasons[2]
    assert reasons[3].startswith("drift ")
    assert "px above limit" in reasons[3]


def test_clip_sigma_rejects_outlier_pixel():
    rng = np.random.default_rng(7)
    spectra = noisy_copies(rng)
    spectra[0].flux[100] += 10 * SIGMA
    unclipped = combine_spectra(spectra, "star", CombineParams(clip_sigma=0.0, log_grid=False))
    assert unclipped.versions["raw"].flux[100] > CONTINUUM + 2 * SIGMA
    cs = combine_spectra(spectra, "star", CombineParams(clip_sigma=4.0, log_grid=False))
    raw = cs.versions["raw"]
    assert abs(raw.flux[100] - CONTINUUM) < 2 * SIGMA
    assert raw.error[100] == pytest.approx(SIGMA / np.sqrt(3))
    assert not raw.mask[100]
    # A gross outlier that drags the mean past every contributor does not empty the pixel.
    spectra[0].flux[120] += 1000.0
    cs = combine_spectra(spectra, "star", CombineParams(clip_sigma=4.0, log_grid=False))
    assert np.isfinite(cs.versions["raw"].flux[120])
    assert not cs.versions["raw"].mask[120]


def test_zero_error_pixel_does_not_warn_in_clipping():
    rng = np.random.default_rng(11)
    params = CombineParams(clip_sigma=4.0, log_grid=False)
    unmasked = noisy_copies(rng)
    unmasked[0].error[10] = 0.0
    cs = combine_spectra(unmasked, "star", params)
    raw = cs.versions["raw"]
    assert not raw.mask[10]
    assert raw.error[10] == pytest.approx(SIGMA / np.sqrt(3))
    masked = noisy_copies(rng)
    masked[0].error[10] = 0.0
    masked[0].mask[10] = True
    cs = combine_spectra(masked, "star", params)
    assert not cs.versions["raw"].mask[10]
    assert cs.versions["raw"].error[10] == pytest.approx(SIGMA / np.sqrt(3))


def test_default_params_log_grid_fractional_offsets():
    rng = np.random.default_rng(12)
    spectra = [
        make_spectrum(
            fiber=i,
            flux=CONTINUUM + rng.normal(0.0, SIGMA, N_PIX),
            wavelength=grid() + 0.01 * (0.37 * i),
            telescope=unit_name(i + 1),
        )
        for i in range(4)
    ]
    cs = combine_spectra(spectra, "star", CombineParams())
    assert np.allclose(np.diff(np.log(cs.wavelength)), np.diff(np.log(cs.wavelength))[0])
    raw = cs.versions["raw"]
    overlap = (cs.wavelength >= spectra[3].wavelength[0]) & (
        cs.wavelength <= spectra[0].wavelength[-1]
    )
    assert overlap.sum() > 150
    assert not raw.mask.any()  # every grid pixel is covered by at least one contributor
    assert np.std(raw.flux[overlap]) == pytest.approx(SIGMA / 2, rel=0.25)
    assert np.std(raw.flux[overlap]) < 0.6 * SIGMA
    assert np.median(raw.error[overlap]) < 0.5 * SIGMA
    assert np.allclose(cs.contributors["weight"], 0.25, atol=0.05)


# ---------------------------------------------------------------------------- versions


def test_versions_combined_separately():
    rng = np.random.default_rng(8)
    spectra = noisy_copies(rng)
    for spectrum in spectra:
        spectrum.versions["skysub"] = SpectrumVersion(
            flux=spectrum.flux - 10.0, error=spectrum.error.copy(), mask=spectrum.mask.copy()
        )
    cs = combine_spectra(spectra, "star", CombineParams(log_grid=False))
    assert set(cs.versions) == {"raw", "skysub"}
    assert np.median(cs.versions["raw"].flux) == pytest.approx(CONTINUUM, abs=1.0)
    assert np.median(cs.versions["skysub"].flux) == pytest.approx(CONTINUUM - 10.0, abs=1.0)

    del spectra[2].versions["skysub"]
    cs = combine_spectra(spectra, "star", CombineParams(log_grid=False))
    assert set(cs.versions) == {"raw"}


def test_version_cap_limits_combined_versions():
    rng = np.random.default_rng(9)
    spectra = noisy_copies(rng)
    for spectrum in spectra:
        version = SpectrumVersion(spectrum.flux.copy(), spectrum.error.copy(), spectrum.mask.copy())
        spectrum.versions["skysub"] = version
        spectrum.versions["fluxcal"] = version
    cs = combine_spectra(spectra, "star", CombineParams(log_grid=False), version="skysub")
    assert set(cs.versions) == {"raw", "skysub"}
    cs = combine_spectra(spectra, "star", CombineParams(log_grid=False))
    assert set(cs.versions) == {"raw", "skysub", "fluxcal"}


def test_version_argument_validated():
    rng = np.random.default_rng(13)
    spectra = noisy_copies(rng)
    for bad in ("raw", "latest", ""):
        with pytest.raises(ValueError, match="skysub"):
            combine_spectra(spectra, "star", CombineParams(), version=bad)


# -------------------------------------------------------------------------- bookkeeping


def test_group_by_target_science_only():
    wave = grid()
    spectra = [
        make_spectrum(0, np.ones(N_PIX), wave, role="science", target="A"),
        make_spectrum(1, np.ones(N_PIX), wave, role="sky", target="A"),
        make_spectrum(2, np.ones(N_PIX), wave, role="science", target="B"),
        make_spectrum(3, np.ones(N_PIX), wave, role="science", target=None),
        make_spectrum(4, np.ones(N_PIX), wave, role="arc", target=None),
        make_spectrum(5, np.ones(N_PIX), wave, role="science", target="A"),
    ]
    groups = group_by_target(spectra)
    assert list(groups) == ["A", "B", "unknown"]
    assert [s.fiber for s in groups["A"]] == [0, 5]
    assert [s.fiber for s in groups["B"]] == [2]
    assert [s.fiber for s in groups["unknown"]] == [3]


def test_other_targets_and_no_wavelength_are_not_combined():
    rng = np.random.default_rng(10)
    spectra = noisy_copies(rng, n_copies=2)
    other = make_spectrum(7, np.ones(N_PIX), grid(), telescope="mast08", target="other")
    sky = make_spectrum(8, np.ones(N_PIX), grid(), role="sky", telescope="mast01", target="star")
    bare = make_spectrum(9, np.ones(N_PIX), None, telescope="mast10", target="star")
    cs = combine_spectra(
        [*spectra, other, sky, bare], "star", CombineParams(log_grid=False), barycentric={"rv": 1.5}
    )
    assert list(cs.contributors["fiber"]) == [0, 1]
    assert list(cs.excluded["fiber"]) == [9]
    assert cs.excluded["reason"].iloc[0] == "no wavelength"
    assert cs.barycentric == {"rv": 1.5}
