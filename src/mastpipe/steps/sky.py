"""Sky subtraction: each science fiber minus the sky fiber of its own telescope.

The two fibers of a telescope see the same patch of sky through different throughputs and
with different wavelength solutions, so the sky fiber's extraction is scaled by the
throughput ratio (science over sky) and resampled in wavelength onto the science fiber's
grid before it is subtracted. The raw extraction (``flux``/``error``/``mask``) is never
replaced: the resampled, scaled sky is kept under ``Spectrum.sky`` and the result goes into
``versions["skysub"]`` (see :data:`mastpipe.core.model.VERSIONS`), so later stages and the
products can show both. A science fiber whose sky fiber is missing, fully masked or without
a wavelength grid keeps its raw extraction only and is reported with the reason. Everything
here is pure: dataclasses and a parameter model in, dataclasses out.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from mastpipe.config.models import SkyParams
from mastpipe.core.model import Spectrum, SpectrumVersion, TraceSet

RESAMPLERS = ("linear",)
"""How :func:`subtract_sky` moves the sky onto the science wavelength grid."""
REASONS = ("no sky fiber", "sky fiber fully masked", "no wavelength", "bad throughput")
"""Why :func:`subtract_sky_all` leaves a science fiber at its raw extraction."""


def _increasing(wavelength: np.ndarray, name: str) -> tuple[np.ndarray, bool]:
    """``wavelength`` as an increasing array and whether it was flipped to get there;
    ``ValueError`` when it is not strictly monotonic."""
    steps = np.diff(wavelength)
    if np.all(steps > 0):
        return wavelength, False
    if np.all(steps < 0):
        return wavelength[::-1], True
    raise ValueError(f"the {name} wavelength grid must be strictly monotonic")


def resample_sky(
    sky: Spectrum, wavelength: np.ndarray, scale: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``sky``'s flux and error scaled by ``scale`` and linearly interpolated onto
    ``wavelength``, with the mask they carry.

    A target wavelength between sky pixels ``j`` and ``j + 1`` at fraction ``t`` takes the
    flux ``(1 - t) f_j + t f_{j+1}`` and the error of that combination,
    ``sqrt((1 - t)² σ_j² + t² σ_{j+1}²)`` (a target exactly on a sky pixel takes that
    pixel's flux and error). The resampled mask marks every target pixel whose
    interpolation gives non-zero weight to a masked sky pixel, and every target pixel
    outside the sky's wavelength coverage, where the flux and error are NaN. A decreasing
    sky grid is flipped for the interpolation. ``ValueError`` when the sky has no
    wavelength grid or it is not monotonic.
    """
    if sky.wavelength is None:
        raise ValueError(f"sky fiber {sky.fiber} has no wavelength grid; cannot resample")
    xp, flipped = _increasing(np.asarray(sky.wavelength, dtype=float), "sky")
    flux = np.asarray(sky.flux, dtype=float) * scale
    error = np.asarray(sky.error, dtype=float) * abs(scale)
    bad = sky.mask
    if flipped:
        flux, error, bad = flux[::-1], error[::-1], bad[::-1]
    target = np.asarray(wavelength, dtype=float)
    inside = (target >= xp[0]) & (target <= xp[-1])
    j = np.clip(np.searchsorted(xp, target, side="right") - 1, 0, len(xp) - 2)
    k = j + 1
    t = np.clip((target - xp[j]) / (xp[k] - xp[j]), 0.0, 1.0)
    on_left = t == 0.0
    on_right = t == 1.0
    between = ~on_left & ~on_right
    flux_out = np.where(on_left, flux[j], np.where(on_right, flux[k], np.nan))
    error_out = np.where(on_left, error[j], np.where(on_right, error[k], np.nan))
    flux_out[between] = (1.0 - t[between]) * flux[j][between] + t[between] * flux[k][between]
    error_out[between] = np.sqrt(
        (1.0 - t[between]) ** 2 * error[j][between] ** 2 + t[between] ** 2 * error[k][between] ** 2
    )
    mask_out = (bad[j] & ~on_right) | (bad[k] & ~on_left)
    flux_out[~inside] = np.nan
    error_out[~inside] = np.nan
    mask_out[~inside] = True
    return flux_out, error_out, mask_out


def subtract_sky(
    science: Spectrum, sky: Spectrum, throughput_ratio: float, params: SkyParams
) -> Spectrum:
    """``science`` with the scaled, resampled ``sky`` subtracted into ``versions["skysub"]``.

    The sky flux is scaled by ``throughput_ratio`` (science throughput over sky throughput)
    and resampled onto the science wavelength grid with :func:`resample_sky`
    (``params.resample`` must be one of :data:`RESAMPLERS`); the returned spectrum is a
    copy of ``science`` with ``sky`` set to that resampled sky and a ``"skysub"`` version
    of ``flux - sky``, error ``sqrt(error² + error_sky²)`` and mask ``mask | mask_sky``.
    The raw ``flux``, ``error`` and ``mask`` are left as they were. ``ValueError`` when
    either spectrum has no wavelength grid, when the sky grid is not monotonic, or when the
    resampler is unknown.
    """
    if params.resample not in RESAMPLERS:
        raise ValueError(f"resample must be one of {RESAMPLERS}, got {params.resample!r}")
    if science.wavelength is None:
        raise ValueError(f"science fiber {science.fiber} has no wavelength grid")
    sky_flux, sky_error, sky_mask = resample_sky(sky, science.wavelength, float(throughput_ratio))
    flux = np.asarray(science.flux, dtype=float) - sky_flux
    error = np.sqrt(np.asarray(science.error, dtype=float) ** 2 + sky_error**2)
    version = SpectrumVersion(flux=flux, error=error, mask=science.mask | sky_mask)
    return replace(science, sky=sky_flux, versions={**science.versions, "skysub": version})


def subtract_sky_all(
    spectra: list[Spectrum], traceset: TraceSet, throughput: dict[int, float], params: SkyParams
) -> tuple[list[Spectrum], list[dict]]:
    """:func:`subtract_sky` applied to every science spectrum with its telescope's sky fiber.

    The sky fiber of a science spectrum is the ``"sky"`` trace of
    ``traceset.by_telescope(spectrum.telescope)``, and its spectrum the one of ``spectra``
    with that ``fiber`` index. The throughput ratio is ``throughput[science] /
    throughput[sky]`` when ``params.scale_by_throughput`` and both fibers are in
    ``throughput``, else 1.0. Spectra whose role is not
    ``"science"`` pass through unchanged, as does a science spectrum that cannot be
    subtracted; those are listed in the second element as ``{"fiber", "telescope",
    "reason"}`` with a reason from :data:`REASONS`: ``"no sky fiber"`` (the spectrum has no
    telescope, the layout has no sky trace for it, or no spectrum carries that fiber),
    ``"sky fiber fully masked"``, ``"no wavelength"`` (either side), or ``"bad throughput"``
    (scaling is on, both fibers are in ``throughput``, and either value is non-finite or not
    positive). The order of ``spectra`` is kept.
    """
    by_fiber = {spectrum.fiber: spectrum for spectrum in spectra}
    out: list[Spectrum] = []
    skipped: list[dict] = []
    for spectrum in spectra:
        if spectrum.role != "science":
            out.append(spectrum)
            continue
        telescope = spectrum.telescope
        sky_trace = None if telescope is None else traceset.by_telescope(telescope).get("sky")
        sky = None if sky_trace is None else by_fiber.get(sky_trace.index)
        if sky is None:
            reason = "no sky fiber"
        elif bool(np.all(sky.mask)):
            reason = "sky fiber fully masked"
        elif spectrum.wavelength is None or sky.wavelength is None:
            reason = "no wavelength"
        else:
            ratio = 1.0
            if (
                params.scale_by_throughput
                and spectrum.fiber in throughput
                and sky.fiber in throughput
            ):
                own = float(throughput[spectrum.fiber])
                other = float(throughput[sky.fiber])
                if not (np.isfinite(own) and np.isfinite(other) and own > 0 and other > 0):
                    reason = "bad throughput"
                    skipped.append(
                        {"fiber": spectrum.fiber, "telescope": telescope, "reason": reason}
                    )
                    out.append(spectrum)
                    continue
                ratio = own / other
            out.append(subtract_sky(spectrum, sky, ratio, params))
            continue
        skipped.append({"fiber": spectrum.fiber, "telescope": telescope, "reason": reason})
        out.append(spectrum)
    return out, skipped
