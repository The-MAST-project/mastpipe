"""Configuration models: instruments are data, validated here, never code.

Every model forbids unknown keys, so a misspelled field in an instrument YAML file or
in a command-line override is an error rather than a silently ignored setting, and
every model is frozen so that nothing downstream can edit a configuration in place.

Step functions never receive these objects; they receive the ``*Params`` models or
plain values. The pipeline layer builds an :class:`EffectiveConfig` per (camera,
configuration) pair with :func:`mastpipe.config.loader.build_effective`.
"""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class DetectorConfig(_Strict):
    """One camera: raw frame geometry, noise model and orientation on read."""

    name: str
    shape: tuple[int, int]
    gain: float
    read_noise: float
    pixel_size_um: float
    saturation: float
    crop: tuple[int, int, int, int] | None
    dispersion_axis: int = 1
    transpose: bool = False
    flip_x: bool = False
    flip_y: bool = False


class GeometryConfig(_Strict):
    """Dispersion geometry of one configuration (the grating model in this plan)."""

    model: Literal["grating", "dispersion"]
    groove_spacing_nm: float
    incidence_angle_deg: float
    order: int
    camera_focal_mm: float
    central_wavelength_nm: float
    dispersion_sign: int = 1


class CalibParams(_Strict):
    subtract_background: bool = True
    flat_subtract_background: bool = False
    stack_method: str = "median"
    clip_sigma: float | None = None
    cosmic_sigma: float = 8.0
    cosmic_size: int | tuple[int, int] = 5  # int: (1, size) along dispersion


class TraceParams(_Strict):
    center_width: int = 100
    trace_halfwidth: int = 8
    threshold_sigma: float = 5.0
    centroid_sigma: float | None = None
    min_separation: float | None = None
    baseline_window: int | None = None
    x_bin: int = 25
    degree: int = 2
    clip_sigma: float = 3.0
    min_points: int = 4
    layout_tolerance: float = 3.0
    max_traces: int | None = None
    collapse: Literal["median", "mean", "percentile"] = "median"
    collapse_percentile: float = 90.0


class ExtractParams(_Strict):
    halfwidth: float = 9.0
    aperture: Literal["round", "subpixel", "floor"] = "round"
    method: Literal["sum", "mean"] = "sum"
    min_coverage: float = 0.99
    continuum_degree: int = 15
    clip_above: float = 2.0
    clip_below: float = 3.0
    continuum_iter: int = 10


class LineParams(_Strict):
    halfwidth: int = 9
    threshold_sigma: float = 5.0
    min_separation: int = 2
    x_halfwidth: int = 6
    y_halfwidth: int | None = None
    blend_distance: float | None = None
    isolation: float | None = None
    edge_margin: int = 10
    baseline_window: int = 51
    sigma_guess: float = 1.3
    max_lines: int | None = None
    max_components: int = 4
    add_sigma: float = 5.0
    max_y_offset: float = 1.5
    saturation: float | None = None


class WavecalParams(_Strict):
    resolving_power: float = 20000.0
    max_shift_fwhm: float = 0.02
    min_snr_atlas: float = 20.0
    usable_only: bool = True
    tol_first_px: float = 2.0
    tol_later_px: float = 1.0
    unique_tol_px: float = 4.0
    n_passes: int = 3
    converge_px: float = 1e-3
    clip_sigma: float = 3.0
    fit_alpha: bool = False
    poly_deg: int = 3
    min_lines: int = 6


class DriftParams(_Strict):
    tol_px: float = 1.5
    max_shift_px: int = 200
    clip_sigma: float = 3.0
    min_lines: int = 5
    stretch_significance: float = 3.0
    degree: int = 1
    max_scatter_px: float = 0.2
    min_fibers: int = 3
    max_age_hours: float = 12.0
    time_policy: Literal["nearest_before", "nearest", "interpolate"] = "nearest_before"
    warn_shift_px: float = 2.0


class SkyParams(_Strict):
    scale_by_throughput: bool = True
    resample: Literal["linear"] = "linear"


class CombineParams(_Strict):
    across_telescopes: bool = True
    include: list[str] = Field(default_factory=list)
    exclude: list[str] = Field(default_factory=list)
    reject: bool = False
    reject_snr_min: float = 3.0
    reject_flux_ratio: float = 0.3
    reject_drift_px: float | None = None
    clip_sigma: float = 4.0
    log_grid: bool = True


class ClassifyParams(_Strict):
    """How a ``light`` frame is resolved into ``arc``, ``science`` or ``unknown`` from its
    content (:mod:`mastpipe.steps.classify`).

    The defaults are the brief's starting values, kept because they separate every frame
    class with a wide margin (statistic of :func:`~mastpipe.steps.classify.arc_likeness`
    measured with them, calibrated as the runner does it, no cosmic-ray flagging):

    - ``signal_sigma = 5``: the 99th percentile of the collapsed spectrum is 12-18 sigma on
      the simulated science frames and 70-4000 on the arcs, against 2-3 on a bias-level
      noise frame or a lamp-off dark and 3.6 on the laboratory lamp-off background (its
      ambient light is diffuse and removed as such);
    - ``arc_peak_fraction = 0.5``: the share of the fiber light in narrow peaks is 0.99 on
      the simulated full and calibration-fiber arcs and 0.69-0.74 on the eleven
      laboratory Th-Ar focus frames, against at most 0.03 on the simulated science frames
      (sky lines included) and 0.00 on the simulated and laboratory flats;
    - ``continuum_window = 101`` columns: wide against a line (a few columns, broader on
      the defocused focus frames) and narrow against the blaze;
    - ``peak_halfwidth_px = 4`` columns grown on either side of a line's core, so that its
      wings count with it.
    """

    continuum_window: int = 101
    signal_sigma: float = 5.0
    arc_peak_fraction: float = 0.5
    peak_halfwidth_px: int = 4


class StepParams(_Strict):
    """The parameters of every step, each group at its defaults unless set."""

    calib: CalibParams = Field(default_factory=CalibParams)
    trace: TraceParams = Field(default_factory=TraceParams)
    extract: ExtractParams = Field(default_factory=ExtractParams)
    lines: LineParams = Field(default_factory=LineParams)
    wavecal: WavecalParams = Field(default_factory=WavecalParams)
    drift: DriftParams = Field(default_factory=DriftParams)
    sky: SkyParams = Field(default_factory=SkyParams)
    combine: CombineParams = Field(default_factory=CombineParams)
    classify: ClassifyParams = Field(default_factory=ClassifyParams)


class ConfigurationBlock(_Strict):
    """One grating (HighSpec) or band (DeepSpec).

    ``params`` holds partial :class:`StepParams` overrides as a plain mapping; they are
    validated when layered onto the instrument parameters by ``build_effective``.
    """

    name: str
    band_nm: tuple[float, float]
    resolving_power: float
    geometry: GeometryConfig
    line_table: str | None = None
    params: dict[str, Any] = Field(default_factory=dict)


class InstrumentConfig(_Strict):
    """An instrument: its cameras, configurations, step defaults and fiber layouts.

    ``layout_params`` holds, per layout preset, partial :class:`StepParams` overrides as a
    plain mapping (the halfwidths a dense fiber head needs, say); ``build_effective`` layers
    the chosen preset's block between the configuration block and the command-line
    overrides. Every key must name a preset listed in ``layouts``.
    """

    name: str
    cameras: dict[str, DetectorConfig]
    configurations: dict[str, ConfigurationBlock]
    params: StepParams
    default_layout: str
    layouts: list[str]
    layout_params: dict[str, dict[str, Any]] = Field(default_factory=dict)
    arc_classification: Literal["header", "lit_traces"] = "lit_traces"

    @model_validator(mode="after")
    def _layout_params_name_known_layouts(self) -> InstrumentConfig:
        unknown = sorted(set(self.layout_params) - set(self.layouts))
        if unknown:
            raise ValueError(f"layout_params for presets not in layouts: {unknown}")
        return self


class EffectiveConfig(_Strict):
    """Everything a run needs for one instrument, camera and configuration.

    ``layers`` names the layers that were merged, in order; it is provenance for the
    log and is not part of the configuration hash.
    """

    instrument: str
    camera: str
    configuration: str
    detector: DetectorConfig
    geometry: GeometryConfig
    band_nm: tuple[float, float]
    resolving_power: float
    line_table: str | None
    params: StepParams
    layers: list[str]
