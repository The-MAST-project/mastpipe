"""Tests for the flux-calibration and telluric-correction interfaces: the pass-through
defaults leave the spectrum and the statuses untouched; a response or transmission adds a
version."""

import numpy as np
import pytest

from mastpipe.core.model import Spectrum, SpectrumVersion
from mastpipe.steps.fluxcal import apply_response, flux_status
from mastpipe.steps.telluric import correct_telluric, telluric_status

N_PIX = 50


def make_spectrum(fiber: int = 0, role: str = "science") -> Spectrum:
    return Spectrum(
        fiber=fiber,
        pixel=np.arange(N_PIX),
        flux=np.full(N_PIX, 100.0),
        error=np.full(N_PIX, 4.0),
        mask=np.zeros(N_PIX, dtype=bool),
        role=role,
        wavelength=515.0 + 0.01 * np.arange(N_PIX),
    )


def test_passthrough_without_calibration_leaves_spectrum_and_status():
    spectra = [make_spectrum(0), make_spectrum(1, role="sky")]
    out = [apply_response(spectrum) for spectrum in spectra]
    assert out[0] is spectra[0]
    assert out[1] is spectra[1]
    assert out[0].versions == {}
    assert flux_status(out) == "uncalibrated"
    out = [correct_telluric(spectrum) for spectrum in out]
    assert out[0] is spectra[0]
    assert out[0].versions == {}
    assert telluric_status(out) == "uncorrected"


def test_passthrough_with_response_adds_version():
    spectrum = make_spectrum()
    spectrum.versions["skysub"] = SpectrumVersion(
        flux=np.full(N_PIX, 90.0), error=np.full(N_PIX, 3.0), mask=np.zeros(N_PIX, dtype=bool)
    )
    response = np.full(N_PIX, 2.0)
    response[3] = 0.0
    response[4] = np.nan
    out = apply_response(spectrum, response)
    assert out is not spectrum
    assert "fluxcal" not in spectrum.versions
    assert set(out.versions) == {"skysub", "fluxcal"}
    calibrated = out.versions["fluxcal"]
    assert calibrated.flux[0] == pytest.approx(45.0)
    assert calibrated.error[0] == pytest.approx(1.5)
    assert calibrated.mask[3] and calibrated.mask[4]
    assert calibrated.mask.sum() == 2
    assert np.array_equal(out.flux, spectrum.flux)
    assert flux_status([out]) == "calibrated"
    assert flux_status([out, make_spectrum(1)]) == "uncalibrated"
    assert flux_status([out, make_spectrum(2, role="sky")]) == "calibrated"

    transmission = np.full(N_PIX, 0.5)
    corrected = correct_telluric(out, transmission)
    assert set(corrected.versions) == {"skysub", "fluxcal", "telluric"}
    assert corrected.versions["telluric"].flux[0] == pytest.approx(90.0)
    assert corrected.versions["telluric"].error[0] == pytest.approx(3.0)
    assert telluric_status([corrected]) == "corrected"
    assert telluric_status([corrected, out]) == "uncorrected"


def test_negative_response_keeps_error_positive():
    spectrum = make_spectrum()
    out = apply_response(spectrum, np.full(N_PIX, -2.0))
    assert out.versions["fluxcal"].flux[0] == pytest.approx(-50.0)
    assert out.versions["fluxcal"].error[0] == pytest.approx(2.0)
