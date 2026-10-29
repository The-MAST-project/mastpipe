# Changes

## 0.0.1 (unreleased)

- Package skeleton: `pyproject.toml` with the `src` layout, `mastpipe` console
  script (`--version`, help), CI matrix (Linux, Windows, macOS × Python 3.11,
  3.12) and the privacy scan for tracked files.
- Core data model (`mastpipe.core.model`): the dataclasses exchanged between
  layers (`FrameMeta`, `Frame`, `MasterCalib`, `CalibrationSet`, `Trace`,
  `TraceSet`, `Spectrum` with kept `versions`, `ExposureSpectra`,
  `WavelengthSolution`, `ReferenceSolution`, `DriftModel`, `CombinedSpectrum`,
  `QCReport`) and the vocabulary constants they validate against.
- Configuration layer (`mastpipe.config`): frozen Pydantic models that forbid unknown
  keys (`DetectorConfig`, `GeometryConfig`, the per-step `*Params` under `StepParams`,
  `ConfigurationBlock`, `InstrumentConfig`, `EffectiveConfig`) and the loader
  (`load_yaml`, `deep_merge`, `build_effective` layering instrument → configuration
  block → overrides, `parse_override` for `key.path=value`, `config_hash`).
- Instrument registry (`mastpipe.instruments.registry`) and the HighSpec definition as
  package data: `config.yaml` (placeholder gain and read noise), `header_map.yaml`, the
  `mgb518` grating block, `lab3` and generated `survey45` layouts, the Mg b line table
  (73 Murphy et al. catalogue lines, 30 usable at R = 23081) and the `mgb518`/`lab3`
  anchor pairs; `scripts/make_survey45_layout.py` and `scripts/make_line_table.py`.
- Ingest (`mastpipe.io`): `HeaderMap` (`from_dict`, `resolve` with `keys`/`type`/`lower`/`map`/
  `default`/`filename_regex`, `classify` with `header`/`rule:<i>`/`unclassified` provenance),
  `scan_folder` (one `FrameRef` per image HDU, primary header merged under extensions,
  unreadable files collected on request), `load_frame` (float32, crop, then transpose /
  flip_x / flip_y), `read_manifest`/`apply_manifest` (CSV overrides by `source` [+ `extension`])
  and `file_checksum` (SHA-256, 16 hex).
- Grating geometry step (`mastpipe.steps.geometry`: `GratingGeometry`, `grating_geometry`,
  `theta_of_lambda`, `grating_x`, `lambda_of_x`, `nominal_solution`) and the synthetic-night
  simulator (`mastpipe.sim.simulator`: `simulate_night` writing raw-shaped HighSpec frames with
  a full arc, calibration-fiber arcs drifting linearly in time, drifted science and sky spectra,
  cosmic rays and a background frame; `SimTruth` saved as `truth.npz` + `truth.json`, read back
  by `load_truth`); session-scoped `sim_night` / `sim_night_lab3` test fixtures.
- Clean step (`mastpipe.steps.calib`, `mastpipe.steps.cosmics`): `stack_frames` (median or
  mean with optional sigma clipping, per-pixel variance of the master), `build_master`
  (validity keys and input-checksum provenance), `build_calibration_set` (dark rate from the
  dark and bias masters), `calibrate_frame` (the notebook's bias/background or bias/dark
  arithmetic, variance in ADU² from gain and read noise, saturation and non-finite mask),
  `calibration_steps` and `flag_cosmic_rays` (median-filter residual above `sigma` noise units).
- Trace step (`mastpipe.steps.trace`): `find_traces` (the notebook's band collapse, running-median
  baseline, `find_peaks`, per-block centroids and two-pass clipped polynomial; the profile noise
  is the statistical noise of the median across the collapsed columns), `comb_shift` (sub-pixel
  shift between two position combs), `match_layout` (global shift, one-to-one assignment by
  distance, `missing`/`extra`), `trace_source` (flat, else full arc), `classify_arc` (`full` /
  `calfiber` / `unknown` from the traces lit above a noise floor and a fraction of the
  brightest, after a spatial running-median baseline) and `relative_throughput`; `TraceParams.collapse`
  (`median` / `mean` / `percentile`) picks the band statistic, the percentile seeing a
  line-only arc when there is no flat.
- Extraction step (`mastpipe.steps.extract`): `aperture_weights` (`round`, `floor` and
  `subpixel` apertures carrying `2·hw + 1` rows of weight), `extract_aperture` (weighted sum
  or mean per column along each trace with its propagated error; a column is masked below
  `min_coverage` or when a masked pixel falls in the aperture; per-trace `aperture` override),
  `fit_continuum` (Chebyshev fit with asymmetric clipping about the median residual, stored
  as `Spectrum.continuum` beside the untouched raw flux) and `continuum_subtract_all`.
