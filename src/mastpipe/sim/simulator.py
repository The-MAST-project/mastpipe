"""A synthetic HighSpec night with the truth it was drawn from.

:func:`simulate_night` writes raw-shaped ``uint16`` FITS frames (bias, dark, flat,
one full arc, a background, calibration-fiber arcs and science frames) whose headers the HighSpec
header map classifies, together with ``truth.npz`` + ``truth.json`` that
:func:`load_truth` reads back as a :class:`SimTruth`. The simulation model:

- the detector reads ``bias_level`` plus Poisson photon/dark noise and Gaussian read
  noise; dark current is ``dark_rate`` ADU per second inside the cropped area, the padded
  border carries bias and read noise only;
- trace ``i`` follows ``y = y_nominal + slope (x - n/2) + c (x - n/2)^2`` with a random
  curvature ``|c| <= 2e-6``, imaged with a Gaussian spatial profile of ``fwhm_px``;
- flats are a quadratic blaze times a per-fiber throughput (0.8-1.2) on every trace;
- the **full arc** lights every trace with the catalogue lines at the exact
  grating-equation pixel of the nominal solution ``(x_c, B)``, drifted like every other
  frame by the law at its own time;
- a **calibration-fiber arc** lights only the ``arc`` traces, its lines moved by
  ``shift + stretch (x - n/2)`` where shift and stretch interpolate linearly in time from
  ``drift_px[0]`` / ``drift_stretch[0]`` at the full arc's mid-time to ``drift_px[1]`` /
  ``drift_stretch[1]`` at the last frame's mid-time (frames before the full arc carry the
  start value, so with the default ``drift_px = (0.0, 0.6)`` the full arc is undrifted);
- science frames carry no arc light: science fibers show a flat continuum with Gaussian
  absorption dips at fixed wavelengths, sky fibers a faint continuum plus emission lines,
  both drifted by the same law at the frame's own time and scaled by the throughput, plus
  ``cosmic_rays`` single-pixel hits of +5000 ADU.

All truth coordinates (traces, line pixels, cosmic rays) are in the cropped and oriented
frame that :func:`mastpipe.io.ingest.load_frame` returns, ``x`` the dispersion column.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.config.models import DetectorConfig, EffectiveConfig
from mastpipe.steps.geometry import (
    GratingGeometry,
    grating_geometry,
    grating_x,
    lambda_of_x,
    nominal_solution,
)

BIAS_LEVEL_ADU = 500.0
DARK_RATE_ADU_S = 0.05
FWHM_PX = 3.0
TRACE_CURVATURE = 2e-6
THROUGHPUT_RANGE = (0.8, 1.2)
FLAT_PEAK_ADU = 20000.0
BLAZE_EDGE_DROP = 0.3
ARC_PEAK_ADU = 40000.0
SCIENCE_CONTINUUM_ADU = 200.0
ABSORPTION_NM = (514.0, 516.0, 518.0, 520.0, 522.0, 524.0)
ABSORPTION_DEPTH = 0.3
ABSORPTION_FWHM_NM = 0.05
SKY_CONTINUUM_ADU = 20.0
SKY_LINES_NM = (515.0, 519.0, 523.0)
SKY_LINE_PEAK_ADU = 300.0
COSMIC_RAY_ADU = 5000.0
EXPOSURE_S = {"bias": 0.0, "dark": 300.0, "flat": 5.0, "arc": 30.0}
TEMPERATURE_C = -60.0
READ_MODE = "Image"
TARGETS = ("target_a", "target_b", "target_c", "target_d")
TELESCOPES_PER_TARGET = 5
PROFILE_HALFWIDTH = 6
ADU_MAX = 65535
GAUSSIAN_FWHM = 2.0 * math.sqrt(2.0 * math.log(2.0))

FRAME_COLUMNS = ("file", "exposure_type", "arc_kind", "mid_time", "shift_px", "stretch")
IMAGE_TYPES = {
    "bias": "bias",
    "dark": "dark",
    "flat": "flat",
    "arc": "arc",
    "science": "light",
    "background": "background",
}
LAMPS = {"arc": "thar", "flat": "flat"}
SHUTTER_CLOSED = ("bias", "dark")
ON_CHIP_MARGIN_PX = 10
TRUTH_ARRAYS = (
    "trace_coeffs",
    "lines_nm",
    "line_logI",
    "line_pixels",
    "lines_on_chip",
    "throughput",
    "science_flux",
    "sky_flux",
    "sky_lines_nm",
)
TRUTH_SCALARS = ("x_c", "B", "bias_level", "dark_rate", "read_noise_e", "gain", "fwhm_px")


@dataclass(eq=False)
class SimTruth:
    """What the simulated frames were drawn from.

    ``layout`` is the fiber layout with a ``target`` column; ``trace_coeffs`` is
    ``(n_traces, 3)`` in ``numpy.polyval`` order, ``y(x)`` in cropped coordinates;
    ``line_pixels`` is ``(n_traces, n_lines)``: the undrifted column of every catalogue line,
    and ``lines_on_chip`` (``(n_lines,)``, bool) marks the lines whose nominal pixel lies
    within ``[ON_CHIP_MARGIN_PX, n_pix - ON_CHIP_MARGIN_PX]``; ``science_flux`` and
    ``sky_flux`` are the undrifted column spectra (ADU per column, throughput 1) of a
    science and a sky fiber; ``frames`` has :data:`FRAME_COLUMNS`, one row per file in
    time order, with the shift and stretch applied to that frame's content (the drift law
    at its mid-time); ``extras`` holds JSON-native scalars and the cosmic-ray positions.
    """

    layout: pd.DataFrame
    trace_coeffs: np.ndarray
    x_c: float
    B: float
    lines_nm: np.ndarray
    line_logI: np.ndarray
    line_pixels: np.ndarray
    lines_on_chip: np.ndarray
    throughput: np.ndarray
    science_flux: np.ndarray
    sky_flux: np.ndarray
    sky_lines_nm: np.ndarray
    bias_level: float
    dark_rate: float
    read_noise_e: float
    gain: float
    fwhm_px: float
    crop: tuple[int, int, int, int] | None
    frames: pd.DataFrame
    extras: dict[str, Any]


@dataclass(frozen=True)
class _Scene:
    """Everything fixed for the night: geometry, traces, throughputs and line positions."""

    layout: pd.DataFrame
    shape: tuple[int, int]
    coeffs: np.ndarray
    throughput: np.ndarray
    geometry: GratingGeometry
    x_c: float
    B: float
    lines_nm: np.ndarray
    line_logI: np.ndarray
    line_x: np.ndarray
    sky_x: np.ndarray

    @property
    def n_pix(self) -> int:
        return self.shape[1]

    @property
    def sigma_px(self) -> float:
        return FWHM_PX / GAUSSIAN_FWHM

    @property
    def column_factor(self) -> float:
        """Column flux that gives a peak pixel of 1 ADU on the trace centre."""
        return self.sigma_px * math.sqrt(2.0 * math.pi)

    def wavelength_of(self, x: np.ndarray) -> np.ndarray:
        return lambda_of_x(x, self.x_c, self.B, self.geometry.alpha_rad, self.geometry)


# --------------------------------------------------------------------------- frame geometry


def oriented_shape(detector: DetectorConfig) -> tuple[int, int]:
    """Shape of a frame after ``load_frame``'s crop and orientation."""
    ny, nx = detector.shape
    if detector.crop is not None:
        x_min, x_max, y_min, y_max = detector.crop
        ny, nx = y_max - y_min, x_max - x_min
    return (nx, ny) if detector.transpose else (ny, nx)


def embed(oriented: np.ndarray, detector: DetectorConfig) -> np.ndarray:
    """The inverse of ``load_frame``: undo the orientation and pad to the raw shape."""
    image = oriented
    if detector.flip_y:
        image = image[::-1, :]
    if detector.flip_x:
        image = image[:, ::-1]
    if detector.transpose:
        image = image.T
    if detector.crop is None:
        if image.shape != tuple(detector.shape):
            raise ValueError(f"image shape {image.shape} differs from detector {detector.shape}")
        return np.array(image, dtype=float)
    raw = np.zeros(detector.shape, dtype=float)
    x_min, x_max, y_min, y_max = detector.crop
    raw[y_min:y_max, x_min:x_max] = image
    return raw


# --------------------------------------------------------------------------------- the scene


def _trace_coeffs(layout: pd.DataFrame, n_pix: int, rng: np.random.Generator) -> np.ndarray:
    """``y = y0 + s (x - xc) + c (x - xc)^2`` expanded to ``polyval`` order, ``xc = n_pix/2``."""
    xc = n_pix / 2.0
    y0 = layout["y_nominal"].to_numpy(dtype=float)
    s = layout["slope_nominal"].to_numpy(dtype=float)
    c = rng.uniform(-TRACE_CURVATURE, TRACE_CURVATURE, size=len(layout))
    return np.column_stack([c, s - 2.0 * c * xc, y0 - s * xc + c * xc**2])


def _assign_targets(layout: pd.DataFrame) -> pd.DataFrame:
    """The layout with a ``target`` column: telescopes in layout order, five per target."""
    telescopes = [t for t in layout["telescope"].tolist() if t is not None]
    unique = list(dict.fromkeys(telescopes))
    target_of = {
        tel: TARGETS[(i // TELESCOPES_PER_TARGET) % len(TARGETS)] for i, tel in enumerate(unique)
    }
    table = layout.copy()
    targets = [None if t is None else target_of[t] for t in layout["telescope"]]
    table["target"] = pd.Series(targets, index=table.index, dtype=object)
    return table


def _build_scene(
    effective: EffectiveConfig,
    layout: pd.DataFrame,
    line_table: pd.DataFrame,
    rng: np.random.Generator,
) -> _Scene:
    shape = oriented_shape(effective.detector)
    n_pix = shape[1]
    geometry = grating_geometry(effective.geometry, effective.detector.pixel_size_um)
    x_c, B = nominal_solution(geometry, n_pix)
    lines_nm = line_table["wavelength_nm"].to_numpy(dtype=float)
    return _Scene(
        layout=_assign_targets(layout),
        shape=shape,
        coeffs=_trace_coeffs(layout, n_pix, rng),
        throughput=rng.uniform(*THROUGHPUT_RANGE, size=len(layout)),
        geometry=geometry,
        x_c=x_c,
        B=B,
        lines_nm=lines_nm,
        line_logI=line_table["logI"].to_numpy(dtype=float),
        line_x=np.asarray(grating_x(lines_nm, x_c, B, geometry.alpha_rad, geometry)),
        sky_x=np.asarray(grating_x(np.array(SKY_LINES_NM), x_c, B, geometry.alpha_rad, geometry)),
    )


# ------------------------------------------------------------------------- column spectra


def drift_x(x_true: np.ndarray, shift: float, stretch: float, n_pix: int) -> np.ndarray:
    """Where a feature truly at ``x_true`` appears: ``x + shift + stretch (x - n_pix/2)``."""
    return x_true + shift + stretch * (x_true - n_pix / 2.0)


def undrift_x(x: np.ndarray, shift: float, stretch: float, n_pix: int) -> np.ndarray:
    """The inverse of :func:`drift_x`."""
    return (x - shift + stretch * n_pix / 2.0) / (1.0 + stretch)


def _gaussian_comb(x: np.ndarray, centres: np.ndarray, amplitudes: np.ndarray, sigma: float):
    """Sum of Gaussians of peak ``amplitudes`` at ``centres``, sampled at ``x``."""
    profile = np.exp(-0.5 * ((x[None, :] - centres[:, None]) / sigma) ** 2)
    return amplitudes @ profile


def _flat_spectrum(scene: _Scene) -> np.ndarray:
    x = np.arange(scene.n_pix, dtype=float)
    half = scene.n_pix / 2.0
    blaze = 1.0 - BLAZE_EDGE_DROP * ((x - half) / half) ** 2
    return FLAT_PEAK_ADU * scene.column_factor * blaze


def _arc_spectrum(scene: _Scene, shift: float, stretch: float) -> np.ndarray:
    x = np.arange(scene.n_pix, dtype=float)
    centres = drift_x(scene.line_x, shift, stretch, scene.n_pix)
    peaks = ARC_PEAK_ADU * 10.0 ** (scene.line_logI - scene.line_logI.max())
    return _gaussian_comb(x, centres, peaks * scene.column_factor, scene.sigma_px)


def _science_spectrum(scene: _Scene, shift: float, stretch: float) -> np.ndarray:
    x = np.arange(scene.n_pix, dtype=float)
    lam = scene.wavelength_of(undrift_x(x, shift, stretch, scene.n_pix))
    sigma_nm = ABSORPTION_FWHM_NM / GAUSSIAN_FWHM
    dips = np.exp(-0.5 * ((lam[None, :] - np.array(ABSORPTION_NM)[:, None]) / sigma_nm) ** 2)
    return SCIENCE_CONTINUUM_ADU * np.prod(1.0 - ABSORPTION_DEPTH * dips, axis=0)


def _sky_spectrum(scene: _Scene, shift: float, stretch: float) -> np.ndarray:
    x = np.arange(scene.n_pix, dtype=float)
    centres = drift_x(scene.sky_x, shift, stretch, scene.n_pix)
    peaks = np.full(centres.size, SKY_LINE_PEAK_ADU * scene.column_factor)
    return SKY_CONTINUUM_ADU + _gaussian_comb(x, centres, peaks, scene.sigma_px)


def _fiber_fluxes(
    scene: _Scene, exposure_type: str, arc_kind: str | None, shift: float, stretch: float
) -> np.ndarray | None:
    """Column spectrum of every trace for one frame, ``None`` when no light reaches the chip."""
    roles = scene.layout["role"].to_numpy()
    fluxes = np.zeros((len(roles), scene.n_pix))
    if exposure_type == "flat":
        fluxes[:] = _flat_spectrum(scene)
    elif exposure_type == "arc":
        lit = np.ones(len(roles), dtype=bool) if arc_kind == "full" else roles == "arc"
        fluxes[lit] = _arc_spectrum(scene, shift, stretch)
    elif exposure_type == "science":
        fluxes[roles == "science"] = _science_spectrum(scene, shift, stretch)
        fluxes[roles == "sky"] = _sky_spectrum(scene, shift, stretch)
    else:
        return None
    return fluxes * scene.throughput[:, None]


# ------------------------------------------------------------------------- frame builders


def render_traces(
    shape: tuple[int, int], coeffs: np.ndarray, fluxes: np.ndarray, fwhm_px: float
) -> np.ndarray:
    """Image of column spectra ``fluxes`` (one row per trace) spread along the traces with a
    Gaussian profile of ``fwhm_px``, normalised so a column's sum equals its flux."""
    ny, nx = shape
    sigma = fwhm_px / GAUSSIAN_FWHM
    norm = 1.0 / (sigma * math.sqrt(2.0 * math.pi))
    image = np.zeros(shape, dtype=float)
    x = np.arange(nx)
    for coeff, flux in zip(coeffs, fluxes, strict=True):
        if not flux.any():
            continue
        y_centre = np.polyval(coeff, x)
        y_round = np.rint(y_centre).astype(int)
        for k in range(-PROFILE_HALFWIDTH, PROFILE_HALFWIDTH + 1):
            rows = y_round + k
            ok = (rows >= 0) & (rows < ny)
            weight = norm * np.exp(-0.5 * ((rows - y_centre) / sigma) ** 2)
            image[rows[ok], x[ok]] += weight[ok] * flux[ok]
    return image


def _read_out(signal: np.ndarray, detector: DetectorConfig, rng: np.random.Generator):
    """Bias plus Poisson noise on the signal (in electrons) plus Gaussian read noise, in ADU."""
    electrons = rng.poisson(np.clip(signal, 0.0, None) * detector.gain)
    read = rng.normal(0.0, detector.read_noise / detector.gain, size=signal.shape)
    return BIAS_LEVEL_ADU + electrons / detector.gain + read


def _cosmic_rays(shape: tuple[int, int], count: int, rng: np.random.Generator) -> list[list[int]]:
    ny, nx = shape
    rows = rng.integers(0, ny, size=count)
    cols = rng.integers(0, nx, size=count)
    return [[int(y), int(x)] for y, x in zip(rows, cols, strict=True)]


def _build_frame(
    scene: _Scene,
    detector: DetectorConfig,
    exposure_type: str,
    arc_kind: str | None,
    exposure: float,
    shift: float,
    stretch: float,
    cosmic_rays: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, list[list[int]]]:
    """A raw ``uint16`` frame and the cropped positions of its cosmic rays."""
    signal = np.full(scene.shape, DARK_RATE_ADU_S * exposure)
    fluxes = _fiber_fluxes(scene, exposure_type, arc_kind, shift, stretch)
    if fluxes is not None:
        signal += render_traces(scene.shape, scene.coeffs, fluxes, FWHM_PX)
    data = _read_out(embed(signal, detector), detector, rng)
    hits: list[list[int]] = []
    if exposure_type == "science" and cosmic_rays > 0:
        hits = _cosmic_rays(scene.shape, cosmic_rays, rng)
        delta = np.zeros(scene.shape)
        for y, x in hits:
            delta[y, x] += COSMIC_RAY_ADU
        data += embed(delta, detector)
    return np.clip(np.rint(data), 0, ADU_MAX).astype(np.uint16), hits


# --------------------------------------------------------------------------- the schedule


def _schedule(
    n_bias: int,
    n_dark: int,
    n_background: int,
    n_flat: int,
    n_full_arc: int,
    n_cal_arc: int,
    n_science: int,
) -> list[tuple[str, str | None]]:
    """``(exposure_type, arc_kind)`` in time order: bias, dark, flat, the full arc, the
    background, then groups of one calibration-fiber arc followed by its share of the
    science frames."""
    order: list[tuple[str, str | None]] = []
    order += [("bias", None)] * n_bias
    order += [("dark", None)] * n_dark
    order += [("flat", None)] * n_flat
    order += [("arc", "full")] * n_full_arc
    order += [("background", None)] * n_background
    if n_cal_arc <= 0:
        order += [("science", None)] * n_science
        return order
    base, extra = divmod(n_science, n_cal_arc)
    for group in range(n_cal_arc):
        order.append(("arc", "calfiber"))
        order += [("science", None)] * (base + (1 if group < extra else 0))
    return order


def _drift_at(seconds: float, origin: float, end: float, endpoints: tuple[float, float]) -> float:
    """The drift law: ``endpoints[0]`` up to ``origin`` seconds, linear to ``endpoints[1]`` at
    ``end`` seconds."""
    if seconds <= origin or end <= origin:
        return float(endpoints[0])
    fraction = (seconds - origin) / (end - origin)
    return float(endpoints[0] + (endpoints[1] - endpoints[0]) * fraction)


def _cards(
    exposure_type: str,
    exposure: float,
    mid_time: str,
    configuration: str,
    fiber_preset: str,
    target: str,
) -> dict[str, Any]:
    return {
        "IMAGETYP": IMAGE_TYPES[exposure_type],
        "EXPOSURE": exposure,
        "DATE-OBS": mid_time,
        "TEMP": TEMPERATURE_C,
        "HBIN": 1,
        "VBIN": 1,
        "READMODE": READ_MODE,
        "GRATING": configuration,
        "LAMP": LAMPS.get(exposure_type, "off"),
        "OBJECT": target,
        "FIBPRE": fiber_preset,
        "SHUTTER": "closed" if exposure_type in SHUTTER_CLOSED else "open",
    }


def _write_fits(path: Path, data: np.ndarray, cards: dict[str, Any]) -> None:
    hdu = fits.PrimaryHDU(data)
    hdu.header.update(cards)
    hdu.writeto(path, overwrite=True)


def _frame_table(rows: list[dict[str, Any]]) -> pd.DataFrame:
    """The frames table with ``arc_kind`` as an object column holding ``None`` when unset."""
    table = pd.DataFrame(rows, columns=list(FRAME_COLUMNS))
    kinds = table["arc_kind"]
    table["arc_kind"] = kinds.astype(object).where(kinds.notna(), None)
    return table


def simulate_night(
    effective: EffectiveConfig,
    layout: pd.DataFrame,
    line_table: pd.DataFrame,
    out_dir: Path | str,
    seed: int = 0,
    n_bias: int = 5,
    n_dark: int = 3,
    n_flat: int = 3,
    n_full_arc: int = 1,
    n_cal_arc: int = 2,
    n_science: int = 4,
    exposure_science: float = 300.0,
    drift_px: tuple[float, float] = (0.0, 0.6),
    drift_stretch: tuple[float, float] = (0.0, 0.0),
    night_start: str = "2026-10-01T18:00:00",
    cadence_min: float = 12.0,
    cosmic_rays: int = 20,
    fiber_preset: str = "survey45",
    no_cal_arcs: bool = False,
    n_background: int = 1,
) -> SimTruth:
    """Write a synthetic night into ``out_dir`` and return (and save) its truth.

    Frames are ``sim_<type>_<nn>.fits`` at ``detector.shape``, one every ``cadence_min``
    minutes from ``night_start`` in the order bias, dark, flat, full arc, background, then
    one calibration-fiber arc followed by its share of the science frames, repeated
    ``n_cal_arc`` times (``no_cal_arcs`` writes the science frames with no such arc). The
    drift law runs from the (first) full arc's mid-time to the last frame's; without a
    full arc it runs from ``night_start``. Every random number comes from
    ``numpy.random.default_rng(seed)``.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    detector = effective.detector
    scene = _build_scene(effective, layout, line_table, rng)
    schedule = _schedule(
        n_bias, n_dark, n_background, n_flat, n_full_arc, 0 if no_cal_arcs else n_cal_arc, n_science
    )
    start = datetime.fromisoformat(night_start)
    kinds = [kind for _, kind in schedule]
    origin_seconds = (kinds.index("full") if "full" in kinds else 0) * cadence_min * 60.0
    end_seconds = (len(schedule) - 1) * cadence_min * 60.0
    exposures = {**EXPOSURE_S, "science": exposure_science, "background": exposure_science}
    targets = sorted({t for t in scene.layout["target"] if t is not None})
    science_object = "+".join(targets) if targets else "sim"

    rows: list[dict[str, Any]] = []
    hits: dict[str, list[list[int]]] = {}
    counters: dict[str, int] = {}
    for k, (exposure_type, arc_kind) in enumerate(schedule):
        seconds = k * cadence_min * 60.0
        shift = _drift_at(seconds, origin_seconds, end_seconds, drift_px)
        stretch = _drift_at(seconds, origin_seconds, end_seconds, drift_stretch)
        exposure = exposures[exposure_type]
        mid_time = (start + timedelta(seconds=seconds)).isoformat(timespec="seconds")
        index = counters.get(exposure_type, 0)
        counters[exposure_type] = index + 1
        name = f"sim_{exposure_type}_{index:02d}.fits"
        data, cosmic = _build_frame(
            scene, detector, exposure_type, arc_kind, exposure, shift, stretch, cosmic_rays, rng
        )
        target = science_object if exposure_type == "science" else exposure_type
        cards = _cards(
            exposure_type, exposure, mid_time, effective.configuration, fiber_preset, target
        )
        _write_fits(out_dir / name, data, cards)
        if exposure_type == "science":
            hits[name] = cosmic
        rows.append(
            {
                "file": name,
                "exposure_type": exposure_type,
                "arc_kind": arc_kind,
                "mid_time": mid_time,
                "shift_px": shift,
                "stretch": stretch,
            }
        )

    truth = SimTruth(
        layout=scene.layout,
        trace_coeffs=scene.coeffs,
        x_c=scene.x_c,
        B=scene.B,
        lines_nm=scene.lines_nm,
        line_logI=scene.line_logI,
        line_pixels=np.tile(scene.line_x, (len(scene.layout), 1)),
        lines_on_chip=(scene.line_x >= ON_CHIP_MARGIN_PX)
        & (scene.line_x <= scene.n_pix - ON_CHIP_MARGIN_PX),
        throughput=scene.throughput,
        science_flux=_science_spectrum(scene, 0.0, 0.0),
        sky_flux=_sky_spectrum(scene, 0.0, 0.0),
        sky_lines_nm=np.array(SKY_LINES_NM, dtype=float),
        bias_level=BIAS_LEVEL_ADU,
        dark_rate=DARK_RATE_ADU_S,
        read_noise_e=detector.read_noise,
        gain=detector.gain,
        fwhm_px=FWHM_PX,
        crop=detector.crop,
        frames=_frame_table(rows),
        extras={
            "seed": seed,
            "configuration": effective.configuration,
            "fiber_preset": fiber_preset,
            "shape": list(scene.shape),
            "n_pix": scene.n_pix,
            "night_start": night_start,
            "cadence_min": cadence_min,
            "drift_px": list(drift_px),
            "drift_stretch": list(drift_stretch),
            "exposures": exposures,
            "absorption_nm": list(ABSORPTION_NM),
            "absorption_depth": ABSORPTION_DEPTH,
            "absorption_fwhm_nm": ABSORPTION_FWHM_NM,
            "science_continuum_adu": SCIENCE_CONTINUUM_ADU,
            "sky_continuum_adu": SKY_CONTINUUM_ADU,
            "sky_line_peak_adu": SKY_LINE_PEAK_ADU,
            "arc_peak_adu": ARC_PEAK_ADU,
            "flat_peak_adu": FLAT_PEAK_ADU,
            "cosmic_ray_adu": COSMIC_RAY_ADU,
            "cosmic_rays": hits,
            "targets": targets,
        },
    )
    save_truth(out_dir, truth)
    return truth


# ------------------------------------------------------------------------------ truth I/O


def _records(table: pd.DataFrame) -> list[dict[str, Any]]:
    """Rows as JSON-native dicts: NaN becomes ``None``, numpy scalars become Python ones."""
    clean = table.astype(object).where(table.notna(), None)
    return [
        {
            key: (value.item() if isinstance(value, np.generic) else value)
            for key, value in row.items()
        }
        for row in clean.to_dict("records")
    ]


def _layout_table(records: list[dict[str, Any]]) -> pd.DataFrame:
    table = pd.DataFrame(records)
    table["trace"] = table["trace"].astype(int)
    for column in ("telescope", "target"):
        table[column] = table[column].astype(object).where(table[column].notna(), None)
    for column in ("y_nominal", "slope_nominal", "aperture"):
        table[column] = table[column].astype(float)
    return table


def save_truth(out_dir: Path | str, truth: SimTruth) -> None:
    """Write ``truth.npz`` (the arrays) and ``truth.json`` (scalars, tables, extras)."""
    out_dir = Path(out_dir)
    np.savez(out_dir / "truth.npz", **{name: getattr(truth, name) for name in TRUTH_ARRAYS})
    payload = {name: getattr(truth, name) for name in TRUTH_SCALARS}
    payload["crop"] = None if truth.crop is None else list(truth.crop)
    payload["layout"] = _records(truth.layout)
    payload["frames"] = _records(truth.frames)
    payload["extras"] = truth.extras
    with (out_dir / "truth.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=1)


def load_truth(out_dir: Path | str) -> SimTruth:
    """Read back what :func:`save_truth` wrote."""
    out_dir = Path(out_dir)
    with (out_dir / "truth.json").open(encoding="utf-8") as handle:
        payload = json.load(handle)
    with np.load(out_dir / "truth.npz") as arrays:
        loaded = {name: arrays[name] for name in TRUTH_ARRAYS}
    crop = payload["crop"]
    return SimTruth(
        layout=_layout_table(payload["layout"]),
        frames=_frame_table(payload["frames"]),
        crop=None if crop is None else tuple(int(v) for v in crop),
        extras=payload["extras"],
        **{name: float(payload[name]) for name in TRUTH_SCALARS},
        **loaded,
    )
