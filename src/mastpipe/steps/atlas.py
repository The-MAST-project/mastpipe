"""Th-Ar atlas tools: which catalogue lines a spectrograph can use at its resolution.

A port of the atlas functions of the HighSpec laboratory reduction notebook, working on the
UVES Th-Ar atlas of Murphy et al. 2007, MNRAS 378, 221 (line list) and Murphy et al. 2010,
ApJ 723, 89 (spectrum). Air wavelengths in nm throughout, as the catalogue gives them.

The chain run by ``scripts/make_line_table.py`` once per configuration is
:func:`load_murphy_atlas` -> :func:`remove_pedestal` -> :func:`convolve_to_resolving_power`
-> :func:`flag_usable_lines`; the result is the shipped line table that the pipeline reads
at run time (``mastpipe.instruments.registry.load_line_table``) and reduces to arrays with
:func:`catalogue_arrays`. :func:`load_murphy_atlas` is the one function in the step modules
that reads files: it takes the atlas folder as an explicit ``Path`` and is only ever called
by that script and by the ``identify`` command's panels, never by the runner.
:func:`synthetic_atlas` exists for the tests.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

SPECTRUM_FILE = "thar_spec_MM201006.fits"
LINE_LIST_FILE = "thar_MM201006.dat"
LINE_LIST_COLUMNS = ("wavenumber", "wavelength", "logI", "element", "ion", "source")
GAUSSIAN_FWHM = 2.3548200
"""FWHM of a unit Gaussian, as the notebook rounds it."""
NOISE_FROM_MAD = 1.4826
"""Standard deviation per median absolute difference of adjacent samples."""
CENTROID_HALFWIDTH_FWHM = 0.8
"""Half-width, in FWHM, of the window whose centroid tests a catalogue line."""
BASELINE_WINDOW = 2001
BASELINE_PERCENTILE = 10


def load_murphy_atlas(
    folder: Path, band: tuple[float, float], pad: float = 0.5
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Atlas spectrum and catalogue over ``band`` (air nm) from the atlas ``folder``.

    The spectrum file is a 2-D image whose row 0 is the flux on a log-linear wavelength
    axis in Angstrom (``wave = 10**(CRVAL1 + i CDELT1) / 10``); it is returned as a
    ``wavelength, flux`` table over ``band`` widened by ``pad`` on both sides, so that a
    line at the band edge keeps its wings for the convolution. The line list is returned as
    ``wavenumber, wavelength, logI, element, ion, source, species`` over ``band`` itself
    (inclusive), ``species`` being ``element ion`` with a trailing ``0``/``1`` ion dropped.

    This is the only file reader among the step modules; it runs offline in
    ``scripts/make_line_table.py`` and in the ``identify`` command's atlas panels, never in
    the pipeline runner, which only reads the shipped line table.
    """
    folder = Path(folder)
    with fits.open(folder / SPECTRUM_FILE) as hdul:
        header = hdul[0].header
        n = int(header["NAXIS1"])
        wave = 10.0 ** (header["CRVAL1"] + np.arange(n) * header["CDELT1"]) / 10.0
        flux = np.asarray(hdul[0].data[0][:n], dtype=float)
    inside = (wave >= band[0] - pad) & (wave <= band[1] + pad)
    spectrum = pd.DataFrame({"wavelength": wave[inside], "flux": flux[inside]})

    catalogue = pd.read_csv(
        folder / LINE_LIST_FILE,
        sep=r"\s+",
        header=None,
        names=list(LINE_LIST_COLUMNS),
        dtype={"ion": str},
        encoding="utf-8",
    )
    catalogue["wavelength"] = catalogue["wavelength"].astype(float) / 10.0
    catalogue["species"] = (catalogue["element"] + " " + catalogue["ion"]).str.replace(
        r" [01]$", "", regex=True
    )
    in_band = (catalogue["wavelength"] >= band[0]) & (catalogue["wavelength"] <= band[1])
    return spectrum, catalogue.loc[in_band].reset_index(drop=True)


def remove_pedestal(flux, window: int = 201, q: float = 5) -> np.ndarray:
    """``flux`` minus its centred rolling ``q``-th percentile over ``window`` samples."""
    flux = np.asarray(flux, dtype=float)
    pedestal = pd.Series(flux).rolling(window, center=True, min_periods=1).quantile(q / 100.0)
    return flux - pedestal.to_numpy()


def convolve_to_resolving_power(wave, flux, R: float) -> tuple[np.ndarray, float]:
    """``flux`` convolved to resolving power ``R`` and the Gaussian's FWHM in nm.

    The FWHM is the mean wavelength over ``R``; the kernel is a normalised Gaussian on the
    (uniform) sample spacing, truncated at three sigma, applied with ``numpy.convolve`` in
    ``same`` mode.
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    fwhm_nm = float(np.mean(wave)) / R
    step = float(np.median(np.diff(wave)))
    sigma_px = fwhm_nm / GAUSSIAN_FWHM / step
    n = int(np.ceil(6.0 * sigma_px)) | 1
    offsets = np.arange(n) - n // 2
    kernel = np.exp(-0.5 * (offsets / sigma_px) ** 2)
    return np.convolve(flux, kernel / kernel.sum(), mode="same"), fwhm_nm


def flag_usable_lines(
    wave,
    flux_conv,
    cat: pd.DataFrame,
    fwhm_nm: float,
    max_shift: float = 0.02,
    min_snr: float = 20.0,
) -> pd.DataFrame:
    """Which catalogue lines still mark their wavelength in the convolved spectrum.

    ``flux_conv`` is the pedestal-free atlas convolved to the instrument's resolution with
    :func:`convolve_to_resolving_power` (``fwhm_nm`` its FWHM). For every catalogue line the
    centroid of the positive flux above a baseline (rolling :data:`BASELINE_PERCENTILE`-th
    percentile over :data:`BASELINE_WINDOW` samples) within
    ``+/-`` :data:`CENTROID_HALFWIDTH_FWHM` FWHM of the catalogue wavelength gives
    ``shift_fwhm`` (in FWHM) and ``shift_pm`` (in pm); ``snr`` is the peak above the baseline
    over the noise (:data:`NOISE_FROM_MAD` times the median absolute sample-to-sample
    difference). ``usable`` is ``|shift_fwhm| < max_shift`` and ``snr > min_snr``: the
    feature the spectrograph records still sits on the catalogue wavelength, no neighbour
    has dragged it away, and it stands out of the atlas noise. Returns a copy of ``cat``
    with the four columns appended.
    """
    wave = np.asarray(wave, dtype=float)
    y = remove_pedestal(flux_conv, window=BASELINE_WINDOW, q=BASELINE_PERCENTILE)
    noise = float(np.median(np.abs(np.diff(y))) * NOISE_FROM_MAD)
    halfwidth = CENTROID_HALFWIDTH_FWHM * fwhm_nm

    out = cat.copy()
    shift, snr = [], []
    for wl in out["wavelength"].to_numpy(dtype=float):
        window = (wave > wl - halfwidth) & (wave < wl + halfwidth)
        yy = np.clip(y[window], 0.0, None)
        centroid = float(np.sum(wave[window] * yy) / yy.sum()) if yy.sum() > 0 else wl
        shift.append((centroid - wl) / fwhm_nm)
        snr.append(float(y[window].max() / noise) if window.any() else 0.0)
    out["shift_fwhm"] = np.asarray(shift, dtype=float)
    out["shift_pm"] = out["shift_fwhm"] * fwhm_nm * 1e3
    out["snr"] = np.asarray(snr, dtype=float)
    out["usable"] = (out["shift_fwhm"].abs() < max_shift) & (out["snr"] > min_snr)
    return out


def synthetic_atlas(
    cat: pd.DataFrame,
    band: tuple[float, float],
    R: float,
    step_nm: float = 0.0005,
    noise: float = 1.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """A test atlas: Gaussian lines of height ``10**logI`` at the catalogue wavelengths.

    The lines are unresolved at ``R`` (their FWHM is a quarter of the band centre over
    ``R``) on a pedestal of 10 with Gaussian ``noise``; ``wave`` runs over ``band`` in steps
    of ``step_nm``. For the tests only.
    """
    rng = np.random.default_rng(seed)
    wave = np.arange(band[0], band[1], step_nm)
    sigma = 0.5 * (band[0] + band[1]) / R / GAUSSIAN_FWHM / 4.0
    flux = np.full_like(wave, 10.0)
    for wl, logI in zip(
        cat["wavelength"].to_numpy(dtype=float), cat["logI"].to_numpy(dtype=float), strict=True
    ):
        flux += 10.0**logI * np.exp(-0.5 * ((wave - wl) / sigma) ** 2)
    flux += rng.normal(0.0, noise, wave.size)
    return wave, flux


def catalogue_arrays(table: pd.DataFrame, usable_only: bool) -> tuple[np.ndarray, np.ndarray]:
    """``(cat_wl, cat_sp)`` of a shipped line table (``wavelength_nm, species, usable``).

    Sorted by wavelength; ``usable_only`` keeps the rows flagged usable.
    """
    rows = table.loc[table["usable"].astype(bool)] if usable_only else table
    rows = rows.sort_values("wavelength_nm", kind="stable")
    return rows["wavelength_nm"].to_numpy(dtype=float), rows["species"].to_numpy(dtype=object)
