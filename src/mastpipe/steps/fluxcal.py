"""Flux calibration: the interface of stage 5, whose method is still open.

The default is a pass-through that adds nothing and leaves the exposure's ``flux_status`` at
``"uncalibrated"``. When a response curve is given, the latest version up to ``"skysub"``
is divided by it into ``versions["fluxcal"]`` (see :data:`mastpipe.core.model.VERSIONS`);
the raw extraction is never replaced. Pure functions: dataclasses in, dataclasses out.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from mastpipe.core.model import Spectrum, SpectrumVersion

INPUT_VERSION = "skysub"
"""The most advanced version a response is applied to."""
OUTPUT_VERSION = "fluxcal"


def divide_version(latest: SpectrumVersion, divisor: np.ndarray) -> SpectrumVersion:
    """``latest`` divided by ``divisor`` pixel by pixel: flux ``f / d``, error ``σ / |d|``,
    with every pixel where the divisor is non-finite or zero masked (its flux and error are
    NaN). ``ValueError`` when the lengths differ."""
    divisor = np.asarray(divisor, dtype=float)
    if len(divisor) != len(latest.flux):
        raise ValueError(f"divisor has length {len(divisor)}, expected {len(latest.flux)}")
    bad = ~np.isfinite(divisor) | (divisor == 0)
    safe = np.where(bad, 1.0, divisor)
    flux = np.where(bad, np.nan, np.asarray(latest.flux, dtype=float) / safe)
    error = np.where(bad, np.nan, np.asarray(latest.error, dtype=float) / np.abs(safe))
    return SpectrumVersion(flux=flux, error=error, mask=latest.mask | bad)


def apply_response(spectrum: Spectrum, response: np.ndarray | None = None) -> Spectrum:
    """``spectrum`` itself when ``response`` is ``None``; otherwise a copy with
    ``versions["fluxcal"]`` = ``latest(upto="skysub")`` divided by ``response`` as in
    :func:`divide_version`."""
    if response is None:
        return spectrum
    version = divide_version(spectrum.latest(upto=INPUT_VERSION), response)
    return replace(spectrum, versions={**spectrum.versions, OUTPUT_VERSION: version})


def flux_status(spectra: list[Spectrum]) -> str:
    """``"calibrated"`` when every ``"science"`` spectrum carries a ``"fluxcal"`` version and
    there is at least one, else ``"uncalibrated"``."""
    science = [spectrum for spectrum in spectra if spectrum.role == "science"]
    if science and all(OUTPUT_VERSION in spectrum.versions for spectrum in science):
        return "calibrated"
    return "uncalibrated"
