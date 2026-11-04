"""Tests for the atlas tools: pedestal removal, convolution to a resolving power, the
usable-line flag on a synthetic atlas, the catalogue arrays of a shipped line table, the
reader on a tiny synthetic atlas, and (when ``MASTPIPE_ATLAS_DIR`` is set) the real atlas."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.steps.atlas import (
    catalogue_arrays,
    convolve_to_resolving_power,
    flag_usable_lines,
    load_murphy_atlas,
    remove_pedestal,
    synthetic_atlas,
)

R = 23081.0
BAND = (513.0, 527.5)
GAUSSIAN_FWHM = 2.0 * np.sqrt(2.0 * np.log(2.0))


def make_catalogue(rows: list[tuple[float, float, str]]) -> pd.DataFrame:
    """A catalogue table from ``(wavelength nm, logI, species)`` rows."""
    table = pd.DataFrame(rows, columns=["wavelength", "logI", "species"])
    table["wavenumber"] = 1e7 / table["wavelength"]
    table["element"] = table["species"].str.split().str[0]
    table["ion"] = table["species"].str.split().str[1]
    table["source"] = "P"
    return table


# ------------------------------------------------------------------------------ tools


def test_remove_pedestal():
    wave = np.arange(0.0, 100.0, 0.1)
    pedestal = 50.0 + 0.2 * wave
    flux = pedestal.copy()
    flux[300] += 1000.0
    flux[700] += 500.0
    cleaned = remove_pedestal(flux, window=201, q=5)
    assert cleaned.shape == flux.shape
    # the ramp is gone to within the window's slope, the lines keep their height
    off = np.ones_like(flux, dtype=bool)
    off[290:311] = off[690:711] = False
    assert np.max(np.abs(cleaned[off])) < 0.2 * 0.1 * 101
    assert cleaned[300] == pytest.approx(1000.0, abs=5.0)
    assert cleaned[700] == pytest.approx(500.0, abs=5.0)
    assert remove_pedestal(np.array([1.0, 2.0, 3.0]), window=5, q=0).tolist() == [0.0, 1.0, 2.0]


def test_convolve_fwhm():
    wave = np.arange(517.0, 519.0, 0.0005)
    flux = np.zeros_like(wave)
    flux[len(wave) // 2] = 100.0
    conv, fwhm_nm = convolve_to_resolving_power(wave, flux, R)
    assert fwhm_nm == pytest.approx(np.mean(wave) / R)
    assert conv.shape == wave.shape
    assert conv.sum() == pytest.approx(flux.sum(), rel=1e-6)
    peak = conv.argmax()
    assert peak == len(wave) // 2
    above = wave[conv >= conv[peak] / 2.0]
    measured = above[-1] - above[0] + 0.0005
    assert measured == pytest.approx(fwhm_nm, abs=0.001)
    sigma_nm = fwhm_nm / GAUSSIAN_FWHM
    expected = 100.0 * np.exp(-0.5 * ((wave - wave[peak]) / sigma_nm) ** 2)
    expected *= flux.sum() / expected.sum()
    np.testing.assert_allclose(conv, expected, atol=0.05)


def test_flag_usable_lines_on_synthetic_atlas():
    cat = make_catalogue(
        [
            (515.0, 3.5, "Th I"),  # isolated and strong
            (518.0, 4.0, "Ar II"),  # isolated and strong
            (520.0, 3.5, "Th I"),  # a blended pair, 0.4 FWHM apart
            (520.009, 3.4, "Th II"),
            (523.0, 1.0, "XX"),  # far too faint for the atlas noise
        ]
    )
    # sampled like the real atlas (about 0.003 nm, seven samples per FWHM) so that the noise
    # estimate from adjacent differences means the same thing
    wave, flux = synthetic_atlas(cat, BAND, R, step_nm=0.003, noise=50.0)
    assert wave[0] >= BAND[0] and wave[-1] <= BAND[1]
    assert np.all(np.diff(wave) > 0)
    conv, fwhm_nm = convolve_to_resolving_power(wave, remove_pedestal(flux), R)
    flagged = flag_usable_lines(wave, conv, cat, fwhm_nm, max_shift=0.02, min_snr=20)
    assert list(flagged.columns[-4:]) == ["shift_fwhm", "shift_pm", "snr", "usable"]
    assert len(flagged) == len(cat)
    assert flagged["usable"].dtype == np.bool_
    usable = dict(zip(flagged["wavelength"], flagged["usable"], strict=True))
    assert usable[515.0]
    assert usable[518.0]
    assert not usable[520.0]
    assert not usable[520.009]
    assert not usable[523.0]
    by = flagged.set_index("wavelength")
    assert abs(by.loc[518.0, "shift_fwhm"]) < 0.005
    assert by.loc[518.0, "snr"] > 20
    assert abs(by.loc[520.0, "shift_fwhm"]) > 0.02
    assert by.loc[523.0, "snr"] < 20
    np.testing.assert_allclose(by["shift_pm"], by["shift_fwhm"] * fwhm_nm * 1e3)
    # the catalogue is returned with its rows in order and untouched
    assert flagged["wavelength"].tolist() == cat["wavelength"].tolist()
    assert "usable" not in cat.columns


def test_catalogue_arrays():
    table = pd.DataFrame(
        {
            "wavelength_nm": [520.0, 514.0, 518.0],
            "species": ["Th I", "Ar II", "Th II"],
            "logI": [3.0, 2.0, 1.0],
            "usable": [True, False, True],
        }
    )
    wl, sp = catalogue_arrays(table, usable_only=False)
    np.testing.assert_allclose(wl, [514.0, 518.0, 520.0])
    assert sp.tolist() == ["Ar II", "Th II", "Th I"]
    assert wl.dtype == float
    wl_u, sp_u = catalogue_arrays(table, usable_only=True)
    np.testing.assert_allclose(wl_u, [518.0, 520.0])
    assert sp_u.tolist() == ["Th II", "Th I"]


def test_catalogue_arrays_of_shipped_table(mg_lines):
    wl, sp = catalogue_arrays(mg_lines, usable_only=True)
    assert len(wl) == int(mg_lines["usable"].sum())
    assert np.all(np.diff(wl) > 0)
    assert len(sp) == len(wl)


# ----------------------------------------------------------------------------- reader


def write_tiny_atlas(folder: Path) -> None:
    """A log-linear two-row spectrum and a catalogue file in the atlas layout."""
    crval = np.log10(5100.0)
    cdelt = 2.5e-6
    n = 4000
    wave_a = 10 ** (crval + np.arange(n) * cdelt)
    flux = 10.0 + 1000.0 * np.exp(-0.5 * ((wave_a - 5180.0) / 0.1) ** 2)
    data = np.vstack([flux, np.ones(n)])
    header = fits.Header()
    header["CRVAL1"] = crval
    header["CDELT1"] = cdelt
    header["CRPIX1"] = 1
    fits.PrimaryHDU(data=data, header=header).writeto(folder / "thar_spec_MM201006.fits")
    lines = [
        "19512.1951  5125.000000  2.500 Th I  P",
        "19305.0193  5180.000000  3.900 Th II  P",
        "19230.7692  5200.000000  1.100 XX 0   P",
        "19047.6190  5250.000000  2.000 Ar II  W",
    ]
    (folder / "thar_MM201006.dat").write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_load_murphy_atlas_parses_header_and_catalogue(tmp_path: Path):
    write_tiny_atlas(tmp_path)
    spectrum, catalogue = load_murphy_atlas(tmp_path, (513.0, 522.0), pad=0.5)
    assert list(spectrum.columns) == ["wavelength", "flux"]
    assert spectrum["wavelength"].min() >= 512.5
    assert spectrum["wavelength"].max() <= 522.5
    assert spectrum["wavelength"].min() < 513.0  # the pad is applied to the spectrum
    assert spectrum["flux"].max() == pytest.approx(1010.0, abs=3.0)  # grid misses the peak
    assert spectrum["wavelength"].iloc[spectrum["flux"].argmax()] == pytest.approx(518.0, abs=1e-3)
    assert list(catalogue.columns) == [
        "wavenumber",
        "wavelength",
        "logI",
        "element",
        "ion",
        "source",
        "species",
    ]
    assert catalogue["wavelength"].tolist() == [518.0, 520.0]  # band is not padded
    assert catalogue["species"].tolist() == ["Th II", "XX"]
    assert catalogue["ion"].tolist() == ["II", "0"]
    assert catalogue.index.tolist() == [0, 1]
    assert catalogue["logI"].tolist() == [3.9, 1.1]


@pytest.mark.skipif("MASTPIPE_ATLAS_DIR" not in os.environ, reason="needs the Murphy atlas")
def test_load_murphy_atlas_real():
    spectrum, catalogue = load_murphy_atlas(Path(os.environ["MASTPIPE_ATLAS_DIR"]), BAND, pad=0.5)
    assert len(catalogue) == 73
    assert catalogue["wavelength"].is_monotonic_increasing
    assert spectrum["wavelength"].iloc[0] == pytest.approx(BAND[0] - 0.5, abs=0.01)
    assert spectrum["wavelength"].iloc[-1] == pytest.approx(BAND[1] + 0.5, abs=0.01)
    wave = spectrum["wavelength"].to_numpy()
    flux = remove_pedestal(spectrum["flux"].to_numpy(), window=201, q=5)
    conv, fwhm_nm = convolve_to_resolving_power(wave, flux, R)
    assert fwhm_nm == pytest.approx(0.02254, abs=1e-4)
    flagged = flag_usable_lines(wave, conv, catalogue, fwhm_nm, max_shift=0.02, min_snr=20)
    n_usable = int(flagged["usable"].sum())
    assert 28 <= n_usable <= 34
    wl, sp = catalogue_arrays(flagged.rename(columns={"wavelength": "wavelength_nm"}), True)
    assert len(wl) == n_usable
