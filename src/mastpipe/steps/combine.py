"""Within-exposure combination: the science spectra of one target, taken through several
telescopes in the same exposure, combined on a common wavelength grid.

Every contributor is resampled onto the grid of :func:`common_grid` (the union of the
contributors' coverage at the finest dispersion among them, log-uniform by default) and the
contributors are averaged pixel by pixel with inverse-variance weights, with one pass of
sigma clipping about the weighted mean where at least three contributors are unmasked.
Every version present in all contributors is combined separately and stored under the same
key, the raw extraction under ``"raw"`` (see
:data:`mastpipe.core.model.COMBINED_VERSIONS`), so the combined product carries the same
stages as its inputs. Spectra that do not take part are never dropped silently: the
``excluded`` table lists each one with its reason, whether it was excluded by configuration
(``include`` / ``exclude`` telescope lists), rejected on quality (``reject``), or had no
wavelength grid. ``Spectrum`` carries no source file name, so the ``source`` column of
both tables is empty in this plan; the runner knows the exposure the combination came
from. Everything here is pure: dataclasses and a parameter model in, a dataclass out.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from mastpipe.config.models import CombineParams
from mastpipe.core.model import (
    CONTRIBUTOR_COLUMNS,
    EXCLUDED_COLUMNS,
    VERSIONS,
    CombinedSpectrum,
    Spectrum,
    SpectrumVersion,
)
from mastpipe.steps.sky import resample_sky

GROUPING = "exposure"
"""The ``grouping`` of every :class:`CombinedSpectrum` this module builds."""
UNKNOWN_TARGET = "unknown"
"""Key of :func:`group_by_target` for science spectra without a target."""
CONFIG_REASON = "excluded by configuration"
NO_WAVELENGTH_REASON = "no wavelength"
MIN_CLIP_CONTRIBUTORS = 3
"""Unmasked contributors a grid pixel needs before sigma clipping is applied to it."""
STEP_TOLERANCE = 1e-6
"""Fraction of a step below which the span's remainder does not add a grid pixel, so a grid
built from identical contributors reproduces their pixel count."""
GRID_TOLERANCE = 1e-9
"""Relative tolerance within which a contributor's grid counts as the common grid and is
used as is."""


def group_by_target(spectra: list[Spectrum]) -> dict[str, list[Spectrum]]:
    """The ``"science"`` spectra of ``spectra`` keyed by ``target``, in first-seen order;
    a ``None`` target is keyed :data:`UNKNOWN_TARGET`."""
    groups: dict[str, list[Spectrum]] = {}
    for spectrum in spectra:
        if spectrum.role != "science":
            continue
        key = UNKNOWN_TARGET if spectrum.target is None else spectrum.target
        groups.setdefault(key, []).append(spectrum)
    return groups


def _coverage(spectrum: Spectrum) -> tuple[float, float, float, float] | None:
    """``(min, max, step, log_step)`` of a spectrum's unmasked wavelengths, ``None`` when
    fewer than two pixels are unmasked."""
    wavelength = np.asarray(spectrum.wavelength, dtype=float)
    good = ~spectrum.mask & np.isfinite(wavelength)
    if good.sum() < 2:
        return None
    kept = wavelength[good]
    steps = np.abs(np.diff(wavelength))
    log_steps = steps / np.minimum(wavelength[:-1], wavelength[1:])
    step, log_step = float(np.median(steps)), float(np.median(log_steps))
    return float(kept.min()), float(kept.max()), step, log_step


def common_grid(spectra: list[Spectrum], log: bool = True) -> np.ndarray:
    """A wavelength grid covering the union of the spectra's unmasked wavelength ranges at the
    finest dispersion among them.

    The range runs from the smallest to the largest unmasked wavelength over all spectra. The
    step is the minimum over spectra of the median pixel step, measured as ``Δλ`` for a
    linear grid and as ``Δλ/λ`` for a log-uniform one (``log``); the number of pixels is
    rounded up so the grid step never exceeds that finest step (a remainder below
    :data:`STEP_TOLERANCE` of a step is ignored). ``ValueError`` when no
    spectrum has a wavelength grid with at least two unmasked pixels.
    """
    coverages = [
        _coverage(spectrum)
        for spectrum in spectra
        if spectrum.wavelength is not None and len(spectrum.wavelength) >= 2
    ]
    coverages = [coverage for coverage in coverages if coverage is not None]
    if not coverages:
        raise ValueError("common_grid needs at least one spectrum with a wavelength grid")
    low = min(coverage[0] for coverage in coverages)
    high = max(coverage[1] for coverage in coverages)
    if log:
        step = min(coverage[3] for coverage in coverages)
        span = np.log(high / low)
    else:
        step = min(coverage[2] for coverage in coverages)
        span = high - low
    if not (step > 0 and np.isfinite(step)):
        raise ValueError("common_grid needs a positive, finite dispersion")
    n = int(np.ceil(span / step - STEP_TOLERANCE)) + 1
    if log:
        return np.geomspace(low, high, n)
    return np.linspace(low, high, n)


def _median_snr(version: SpectrumVersion) -> float:
    good = ~version.mask & np.isfinite(version.flux) & np.isfinite(version.error)
    good &= version.error > 0
    if not good.any():
        return 0.0
    return float(np.median(np.abs(version.flux[good]) / version.error[good]))


def _median_flux(version: SpectrumVersion) -> float:
    good = ~version.mask & np.isfinite(version.flux)
    if not good.any():
        return np.nan
    return float(np.median(version.flux[good]))


def _reject(
    candidates: list[Spectrum],
    params: CombineParams,
    version: str | None,
    qc: dict[int, dict] | None,
) -> tuple[list[Spectrum], list[tuple[Spectrum, str]]]:
    """``candidates`` split into the kept spectra and ``(spectrum, reason)`` pairs rejected on
    S/N, flux ratio to the group median, or drift."""
    latest = {spectrum.fiber: spectrum.latest(upto=version) for spectrum in candidates}
    medians = np.array([_median_flux(latest[spectrum.fiber]) for spectrum in candidates])
    finite = np.isfinite(medians)
    group_median = float(np.median(medians[finite])) if finite.any() else np.nan
    kept: list[Spectrum] = []
    rejected: list[tuple[Spectrum, str]] = []
    for spectrum, median in zip(candidates, medians, strict=True):
        snr = _median_snr(latest[spectrum.fiber])
        if snr < params.reject_snr_min:
            rejected.append((spectrum, f"low S/N ({snr:.3g} < {params.reject_snr_min})"))
            continue
        if np.isfinite(group_median) and group_median > 0:
            ratio = median / group_median
            if not (ratio >= params.reject_flux_ratio):
                rejected.append(
                    (spectrum, f"flux ratio {ratio:.3g} below min {params.reject_flux_ratio}")
                )
                continue
        if params.reject_drift_px is not None and qc is not None:
            drift = (qc.get(spectrum.fiber) or {}).get("drift_px")
            if drift is not None and abs(float(drift)) > params.reject_drift_px:
                rejected.append(
                    (
                        spectrum,
                        f"drift {abs(float(drift)):.3g} px above limit {params.reject_drift_px}",
                    )
                )
                continue
        kept.append(spectrum)
    return kept, rejected


def _resampled(spectrum: Spectrum, version: SpectrumVersion, wavelength: np.ndarray) -> tuple:
    """The version's flux, error and mask on ``wavelength`` (see
    :func:`mastpipe.steps.sky.resample_sky`); the arrays themselves when the spectrum's grid
    is ``wavelength`` within :data:`GRID_TOLERANCE`."""
    own = np.asarray(spectrum.wavelength, dtype=float)
    if len(own) == len(wavelength) and np.allclose(own, wavelength, rtol=GRID_TOLERANCE, atol=0):
        return (
            np.asarray(version.flux, dtype=float),
            np.asarray(version.error, dtype=float),
            np.asarray(version.mask, dtype=bool),
        )
    carrier = Spectrum(
        fiber=spectrum.fiber,
        pixel=spectrum.pixel,
        flux=version.flux,
        error=version.error,
        mask=version.mask,
        wavelength=spectrum.wavelength,
    )
    return resample_sky(carrier, wavelength, 1.0)


def _stack(
    contributors: list[Spectrum], key: str | None, wavelength: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(flux, error, mask)`` stacks of shape ``(n_contributors, n_grid)`` of the
    contributors' version ``key`` (``None`` for the raw extraction) on ``wavelength``."""
    rows = []
    for spectrum in contributors:
        version = (
            SpectrumVersion(spectrum.flux, spectrum.error, spectrum.mask)
            if key is None
            else spectrum.versions[key]
        )
        rows.append(_resampled(spectrum, version, wavelength))
    flux = np.array([row[0] for row in rows], dtype=float)
    error = np.array([row[1] for row in rows], dtype=float)
    mask = np.array([row[2] for row in rows], dtype=bool)
    return flux, error, mask


def _weighted_mean(
    flux: np.ndarray, error: np.ndarray, mask: np.ndarray, clip_sigma: float | None
) -> tuple[SpectrumVersion, np.ndarray]:
    """Inverse-variance weighted mean of the stack rows per column, with one pass of sigma
    clipping about the weighted mean where at least :data:`MIN_CLIP_CONTRIBUTORS` rows are
    unmasked (a pixel where clipping would reject every row keeps its unclipped mean); the
    version and the weight share of every row (``(n_rows, n_grid)``, zero where the row does
    not contribute)."""
    usable = ~mask & np.isfinite(flux) & np.isfinite(error) & (error > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        weight = np.where(usable, 1.0 / error**2, 0.0)
    mean, combined_error, share = _average(flux, weight)
    if clip_sigma is not None and clip_sigma > 0:
        n_usable = usable.sum(axis=0)
        safe_error = np.where(usable, error, 1.0)
        deviation = np.where(usable, np.abs(flux - mean[None, :]) / safe_error, 0.0)
        outlier = usable & (deviation > clip_sigma) & (n_usable[None, :] >= MIN_CLIP_CONTRIBUTORS)
        outlier &= (outlier.sum(axis=0) < n_usable)[None, :]
        if outlier.any():
            weight = np.where(outlier, 0.0, weight)
            mean, combined_error, share = _average(flux, weight)
    out_mask = ~(weight > 0).any(axis=0)
    return SpectrumVersion(flux=mean, error=combined_error, mask=out_mask), share


def _average(flux: np.ndarray, weight: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    total = weight.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = np.where(total > 0, (weight * np.nan_to_num(flux)).sum(axis=0) / total, np.nan)
        error = np.where(total > 0, 1.0 / np.sqrt(total), np.nan)
        share = np.where(total[None, :] > 0, weight / total[None, :], 0.0)
    return mean, error, share


def _rows(pairs: list[tuple[Spectrum, str]], column: str, columns: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "fiber": spectrum.fiber,
                "telescope": spectrum.telescope,
                "source": "",
                column: value,
            }
            for spectrum, value in pairs
        ],
        columns=list(columns),
    )


def combine_spectra(
    spectra: list[Spectrum],
    target: str,
    params: CombineParams,
    version: str | None = None,
    barycentric: dict[str, float] | None = None,
    qc: dict[int, dict] | None = None,
) -> CombinedSpectrum:
    """The science spectra of ``target`` among ``spectra`` combined on their common grid.

    Candidates are the spectra with ``role == "science"`` and that ``target``. When
    ``params.include`` is non-empty only the telescopes it lists take part; the telescopes
    of ``params.exclude`` never do; either way the spectra left out go to ``excluded`` with
    the reason ``"excluded by configuration"``, and a candidate without a wavelength grid
    goes there with ``"no wavelength"``. When ``params.reject`` a candidate is excluded with
    the reason when the median S/N of its ``latest(upto=version)`` is below
    ``params.reject_snr_min`` (``"low S/N (x < min)"``), when its median flux is below
    ``params.reject_flux_ratio`` times the group's median of those medians (``"flux ratio x
    below min"``), or, when ``params.reject_drift_px`` and ``qc`` are given, when
    ``|qc[fiber]["drift_px"]|`` exceeds the limit (``"drift x px above limit"``).

    The contributors are resampled onto :func:`common_grid` (log-uniform when
    ``params.log_grid``) exactly as :func:`mastpipe.steps.sky.resample_sky` does (linear
    interpolation with the error of the interpolated pair, a pixel masked when a masked
    neighbour or the outside of the coverage enters it), and averaged with inverse-variance
    weights: error ``1 / sqrt(Σ w)``, masked where no contributor is unmasked, with one pass
    of ``params.clip_sigma`` clipping about the weighted mean at pixels with at least three
    unmasked contributors (never emptying a pixel; a non-positive ``clip_sigma`` disables
    the clipping). The raw extractions are combined under ``"raw"`` and every
    version of :data:`mastpipe.core.model.VERSIONS` (up to ``version`` when given) that
    all contributors carry is combined separately under its own key. ``contributors`` holds
    ``fiber, telescope, source, weight`` with the weight the contributor's mean share of
    the raw combination's weights (``1 / n`` when no pixel was combined); ``source`` is
    empty, since a :class:`Spectrum` carries no file name. ``grouping`` is
    :data:`GROUPING` and ``barycentric`` a copy of the given terms. With no contributor
    the grid is empty and ``versions`` holds an empty ``"raw"``. ``ValueError`` when
    ``version`` is neither ``None`` nor one of :data:`mastpipe.core.model.VERSIONS`.
    """
    if version is not None and version not in VERSIONS:
        raise ValueError(f"version must be None or one of {VERSIONS}, got {version!r}")
    candidates = [s for s in spectra if s.role == "science" and s.target == target]
    excluded: list[tuple[Spectrum, str]] = []
    kept: list[Spectrum] = []
    for spectrum in candidates:
        telescope = spectrum.telescope
        if (params.include and telescope not in params.include) or telescope in params.exclude:
            excluded.append((spectrum, CONFIG_REASON))
        elif spectrum.wavelength is None:
            excluded.append((spectrum, NO_WAVELENGTH_REASON))
        else:
            kept.append(spectrum)
    if params.reject and kept:
        kept, rejected = _reject(kept, params, version, qc)
        excluded.extend(rejected)

    keys: list[str] = []
    if kept:
        wavelength = common_grid(kept, log=params.log_grid)
        allowed = VERSIONS if version is None else VERSIONS[: VERSIONS.index(version) + 1]
        keys = [key for key in allowed if all(key in s.versions for s in kept)]
    else:
        wavelength = np.zeros(0, dtype=float)

    versions: dict[str, SpectrumVersion] = {}
    shares = np.zeros((len(kept), len(wavelength)))
    if kept:
        for key in (None, *keys):
            combined, share = _weighted_mean(
                *_stack(kept, key, wavelength), clip_sigma=params.clip_sigma
            )
            versions["raw" if key is None else key] = combined
            if key is None:
                shares = share
    else:
        empty = np.zeros(0, dtype=float)
        versions["raw"] = SpectrumVersion(empty, empty.copy(), np.zeros(0, dtype=bool))

    combined_pixels = shares.sum(axis=0) > 0
    if combined_pixels.any():
        weights = shares[:, combined_pixels].mean(axis=1)
    else:
        weights = np.full(len(kept), 1.0 / len(kept) if kept else 0.0)
    contributors = _rows(
        list(zip(kept, (float(w) for w in weights), strict=True)), "weight", CONTRIBUTOR_COLUMNS
    )
    return CombinedSpectrum(
        target=target,
        wavelength=wavelength,
        versions=versions,
        contributors=contributors,
        excluded=_rows(excluded, "reason", EXCLUDED_COLUMNS),
        grouping=GROUPING,
        barycentric=dict(barycentric or {}),
    )
