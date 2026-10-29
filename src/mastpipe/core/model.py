"""Core data model: the dataclasses that every pipeline layer exchanges.

Plain dataclasses over numpy arrays and pandas tables. Nothing here reads a file
or a configuration. Conventions:

- Masks are boolean arrays with ``True`` marking a bad pixel.
- A :class:`Spectrum` keeps its raw extraction in ``flux``/``error``/``mask``;
  later stages never replace them but add a :class:`SpectrumVersion` under one
  of the keys in :data:`VERSIONS`, in that processing order.
- Times are ISO-8601 UTC strings; wavelengths are air wavelengths in nm.
- Polynomial coefficients follow the ``numpy.polyval`` order (highest degree first).

Dataclasses that hold arrays are declared with ``eq=False`` so that comparing two
instances falls back to identity instead of an ambiguous element-wise test.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

EXPOSURE_TYPES = ("bias", "dark", "flat", "arc", "science", "background", "unknown")
ARC_KINDS = ("full", "calfiber", "unknown")
ROLES = ("arc", "science", "sky", "unknown")
WAVELENGTH_STATUSES = ("none", "nominal", "reference", "transferred", "solved")
FLUX_STATUSES = ("uncalibrated", "calibrated")
TELLURIC_STATUSES = ("uncorrected", "corrected")
VERSIONS = ("skysub", "fluxcal", "telluric")
COMBINED_VERSIONS = ("raw", *VERSIONS)

LINE_COLUMNS = ("x", "wavelength", "species", "dx_px")
PER_FIBER_COLUMNS = ("fiber", "y_center", "shift_px", "stretch", "n_lines", "rms_px", "usable")
CONTRIBUTOR_COLUMNS = ("fiber", "telescope", "source", "weight")
EXCLUDED_COLUMNS = ("fiber", "telescope", "source", "reason")


def _check_choice(name: str, value: Any, choices: tuple[str, ...]) -> None:
    if value not in choices:
        raise ValueError(f"{name} must be one of {choices}, got {value!r}")


def _check_length(name: str, array: np.ndarray, n: int) -> None:
    if len(array) != n:
        raise ValueError(f"{name} has length {len(array)}, expected {n}")


def _check_mask(name: str, mask: np.ndarray) -> None:
    if mask.dtype != np.bool_:
        raise ValueError(f"{name} must be a bool array (True = bad), got dtype {mask.dtype}")


def _check_columns(name: str, table: pd.DataFrame, columns: tuple[str, ...]) -> None:
    missing = [column for column in columns if column not in table.columns]
    if missing:
        raise ValueError(f"{name} is missing columns {missing}; expected {list(columns)}")


@dataclass(frozen=True)
class FrameMeta:
    """Pipeline vocabulary for one raw frame, filled by ingest from the header map."""

    instrument: str
    camera: str
    exposure_type: str
    configuration: str | None = None
    exposure_time: float | None = None
    mid_time: str | None = None
    binning: tuple[int, int] = (1, 1)
    readout_mode: str | None = None
    detector_temperature: float | None = None
    lamp_state: str | None = None
    fiber_preset: str | None = None
    target: str | None = None
    arc_kind: str | None = None
    checksum: str = ""
    source: str = ""
    extension: int = 0
    unresolved: tuple[str, ...] = ()
    classification_rule: str | None = None

    def __post_init__(self) -> None:
        _check_choice("exposure_type", self.exposure_type, EXPOSURE_TYPES)
        if self.arc_kind is not None:
            _check_choice("arc_kind", self.arc_kind, ARC_KINDS)


@dataclass(frozen=True)
class FrameRef:
    """A frame on disk: file, extension and its classified metadata."""

    path: Path
    extension: int
    meta: FrameMeta


@dataclass(eq=False)
class Frame:
    """A 2-D image with optional variance and bad-pixel mask of the same shape."""

    data: np.ndarray
    meta: FrameMeta
    variance: np.ndarray | None = None
    mask: np.ndarray | None = None

    def __post_init__(self) -> None:
        for name in ("variance", "mask"):
            array = getattr(self, name)
            if array is not None and array.shape != self.data.shape:
                raise ValueError(f"{name} shape {array.shape} differs from data {self.data.shape}")
        if self.mask is not None:
            _check_mask("mask", self.mask)

    @property
    def shape(self) -> tuple[int, ...]:
        return self.data.shape


@dataclass(eq=False)
class MasterCalib:
    """A stacked calibration image with the keys that decide where it applies."""

    kind: str
    data: np.ndarray
    validity: dict[str, Any]
    n_frames: int
    provenance: dict[str, Any] = field(default_factory=dict)
    variance: np.ndarray | None = None


@dataclass(eq=False)
class CalibrationSet:
    """The calibration inputs for one detector: masters, scaling times and noise model."""

    bias: np.ndarray | None
    dark_rate: np.ndarray | float
    background: np.ndarray | None
    t_background: float | None
    gain: float
    read_noise: float
    shape: tuple[int, int]


@dataclass(eq=False)
class Trace:
    """One fiber trace: polynomial ``y(x)`` plus the centroid samples it was fitted to."""

    index: int
    coeffs: np.ndarray
    y_center: float
    slope: float
    rms: float
    n_points: int
    snr: float
    x_samples: np.ndarray
    y_samples: np.ndarray
    kept: np.ndarray
    role: str | None = None
    telescope: str | None = None
    layout_index: int | None = None
    aperture: float | None = None

    def y(self, x: float | np.ndarray) -> float | np.ndarray:
        """Trace centre at column(s) ``x``."""
        return np.polyval(self.coeffs, x)


@dataclass(eq=False)
class TraceSet:
    """All traces found on one frame, with the layout-matching outcome."""

    traces: list[Trace]
    shape: tuple[int, int]
    layout_shift: float | None = None
    missing: tuple[int, ...] = ()
    extra: tuple[int, ...] = ()

    def __len__(self) -> int:
        return len(self.traces)

    def __iter__(self) -> Iterator[Trace]:
        return iter(self.traces)

    def by_role(self, role: str) -> list[Trace]:
        """Traces with the given role, in trace order."""
        return [trace for trace in self.traces if trace.role == role]

    def by_telescope(self, telescope: str) -> dict[str, Trace]:
        """Traces of one telescope keyed by role (one trace per role in a layout)."""
        return {trace.role: trace for trace in self.traces if trace.telescope == telescope}


@dataclass(eq=False)
class SpectrumVersion:
    """Flux, error and mask of one processing stage on the spectrum's pixel grid."""

    flux: np.ndarray
    error: np.ndarray
    mask: np.ndarray

    def __post_init__(self) -> None:
        n = len(self.flux)
        _check_length("error", self.error, n)
        _check_length("mask", self.mask, n)
        _check_mask("mask", self.mask)


@dataclass(eq=False)
class Spectrum:
    """One fiber's extraction: raw arrays, optional grids, and later versions."""

    fiber: int
    pixel: np.ndarray
    flux: np.ndarray
    error: np.ndarray
    mask: np.ndarray
    role: str | None = None
    telescope: str | None = None
    target: str | None = None
    wavelength: np.ndarray | None = None
    sky: np.ndarray | None = None
    continuum: np.ndarray | None = None
    versions: dict[str, SpectrumVersion] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n = len(self.pixel)
        for name in ("flux", "error", "mask", "wavelength", "sky", "continuum"):
            array = getattr(self, name)
            if array is not None:
                _check_length(name, array, n)
        _check_mask("mask", self.mask)
        for key, version in self.versions.items():
            _check_choice("versions key", key, VERSIONS)
            _check_length(f"versions[{key!r}]", version.flux, n)

    def latest(self, upto: str | None = None) -> SpectrumVersion:
        """The most advanced version present, in :data:`VERSIONS` order.

        ``upto`` caps the versions considered (``"skysub"`` ignores ``fluxcal`` and
        ``telluric``); ``None`` means no cap. When no version qualifies the raw
        arrays are returned wrapped in a :class:`SpectrumVersion`.
        """
        if upto is None:
            candidates = VERSIONS
        else:
            _check_choice("upto", upto, VERSIONS)
            candidates = VERSIONS[: VERSIONS.index(upto) + 1]
        for key in reversed(candidates):
            if key in self.versions:
                return self.versions[key]
        return SpectrumVersion(self.flux, self.error, self.mask)


@dataclass(eq=False)
class ExposureSpectra:
    """Every spectrum of one exposure, with the exposure's metadata and statuses."""

    meta: FrameMeta
    spectra: list[Spectrum]
    wavelength_status: str = "none"
    flux_status: str = "uncalibrated"
    telluric_status: str = "uncorrected"
    reference_id: str | None = None
    reference_hash: str | None = None
    drift_arcs: tuple[str, ...] = ()
    provenance: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        _check_choice("wavelength_status", self.wavelength_status, WAVELENGTH_STATUSES)
        _check_choice("flux_status", self.flux_status, FLUX_STATUSES)
        _check_choice("telluric_status", self.telluric_status, TELLURIC_STATUSES)


@dataclass(eq=False)
class WavelengthSolution:
    """One fiber's wavelength model and the identified lines that constrain it."""

    fiber: int
    kind: str
    params: np.ndarray
    n_pix: int
    lines: pd.DataFrame
    rms_px: float
    n_lines: int
    anchor_px: tuple[float, float]
    anchor_nm: tuple[float, float]
    trace_y_center: float
    geometry: dict[str, Any]

    def __post_init__(self) -> None:
        _check_columns("lines", self.lines, LINE_COLUMNS)


@dataclass(eq=False)
class ReferenceSolution:
    """Per-fiber wavelength solutions from one full arc, valid from its mid-time on."""

    instrument: str
    configuration: str
    fiber_preset: str
    solutions: dict[int, WavelengthSolution]
    created: str
    valid_from: str
    source_checksums: tuple[str, ...]
    identifier: str
    hash: str = ""


def make_identifier(instrument: str, configuration: str, fiber_preset: str, created: str) -> str:
    """Reference identifier ``{instrument}_{configuration}_{fiber_preset}_{YYYYmmddTHHMMSS}``."""
    stamp = datetime.fromisoformat(created).strftime("%Y%m%dT%H%M%S")
    return f"{instrument}_{configuration}_{fiber_preset}_{stamp}"


@dataclass(eq=False)
class DriftModel:
    """Shift and stretch of the calibration fibers against a reference, smoothed across y."""

    arc_checksums: tuple[str, ...]
    arc_mid_time: str
    per_fiber: pd.DataFrame
    shift_coeffs: np.ndarray
    stretch_coeffs: np.ndarray
    scatter_px: float
    usable: bool
    n_pix: int

    def __post_init__(self) -> None:
        _check_columns("per_fiber", self.per_fiber, PER_FIBER_COLUMNS)

    def shift_at(self, y: float | np.ndarray) -> float | np.ndarray:
        """Pixel shift at trace centre ``y``."""
        return np.polyval(self.shift_coeffs, y)

    def stretch_at(self, y: float | np.ndarray) -> float | np.ndarray:
        """Dispersion stretch at trace centre ``y``."""
        return np.polyval(self.stretch_coeffs, y)


@dataclass(eq=False)
class CombinedSpectrum:
    """Spectra of one target combined on a common grid, version by version."""

    target: str
    wavelength: np.ndarray
    versions: dict[str, SpectrumVersion]
    contributors: pd.DataFrame
    excluded: pd.DataFrame
    grouping: str
    barycentric: dict[str, float]

    def __post_init__(self) -> None:
        n = len(self.wavelength)
        for key, version in self.versions.items():
            _check_choice("versions key", key, COMBINED_VERSIONS)
            _check_length(f"versions[{key!r}]", version.flux, n)
        _check_columns("contributors", self.contributors, CONTRIBUTOR_COLUMNS)
        _check_columns("excluded", self.excluded, EXCLUDED_COLUMNS)


@dataclass(eq=False)
class QCReport:
    """Scalar quality numbers of one stage, serialisable to JSON."""

    stage: str
    values: dict[str, Any]
