"""Telluric correction: the interface of stage 6, whose method is still open.

The default is a pass-through that adds nothing and leaves the exposure's
``telluric_status`` at ``"uncorrected"``. When a transmission curve is given, the latest
version up to ``"fluxcal"`` is divided by it into ``versions["telluric"]`` (see
:data:`mastpipe.core.model.VERSIONS`); the raw extraction is never replaced. Pure
functions: dataclasses in, dataclasses out.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from mastpipe.core.model import Spectrum
from mastpipe.steps.fluxcal import divide_version

INPUT_VERSION = "fluxcal"
"""The most advanced version a transmission is applied to."""
OUTPUT_VERSION = "telluric"


def correct_telluric(spectrum: Spectrum, transmission: np.ndarray | None = None) -> Spectrum:
    """``spectrum`` itself when ``transmission`` is ``None``; otherwise a copy with
    ``versions["telluric"]`` = ``latest(upto="fluxcal")`` divided by ``transmission`` as in
    :func:`mastpipe.steps.fluxcal.divide_version`."""
    if transmission is None:
        return spectrum
    version = divide_version(spectrum.latest(upto=INPUT_VERSION), transmission)
    return replace(spectrum, versions={**spectrum.versions, OUTPUT_VERSION: version})


def telluric_status(spectra: list[Spectrum]) -> str:
    """``"corrected"`` when every ``"science"`` spectrum carries a ``"telluric"`` version and
    there is at least one, else ``"uncorrected"``."""
    science = [spectrum for spectrum in spectra if spectrum.role == "science"]
    if science and all(OUTPUT_VERSION in spectrum.versions for spectrum in science):
        return "corrected"
    return "uncorrected"
