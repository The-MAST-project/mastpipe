# Changes

## 0.0.1 (unreleased)

- Final fix wave of Plan 1: fibers are keyed by the layout's trace id
  (`steps.trace.layout_indexed`, used by the runner and `identify`; detections
  outside the layout are dropped, logged and listed in the trace QC), with a
  trace-centre guard in `apply_drift` (`wavelength_skipped` in the science QC)
  and a pixel-grid check against the reference; header-map options
  `missing_values` (Andor `TEMP = -999`) and `reference` (`start | mid | end`)
  with the mid-time derived by ingest from the exposure start, the simulator
  stamping `DATE-OBS` with the start; per-preset parameters
  (`InstrumentConfig.layout_params`, `build_effective(layout=...)`,
  `survey45` halfwidths and aperture of 4 px); a non-existent night folder
  raises and an empty one records a failure; the layout preset joins the
  science idempotence key (`MPLAYOUT`); the `lt` classification operator
  ignores non-numeric values; README, crop comment and `slow` marker corrected.
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
- Line step (`mastpipe.steps.lines`): `find_lines` (the notebook's 2-D line finder: strip
  detilted onto each trace, running-median continuum per row, detection on the y-collapsed
  profile, groups fitted as sums of Gaussians with the residual-driven deblender, windows
  clipped at the midpoint to the neighbours, the `isolated`/`near_edge`/`complete`/`on_trace`/
  `fit_ok` flags and the `use_x`/`use_y` verdicts; adds a `saturated` flag that bars `use_x`
  and a variance-based `snr` when a variance image is given) and `lines_summary`.
- Wavelength model and atlas tools (`mastpipe.steps.wavemodel`, `mastpipe.steps.atlas`):
  the `WavelengthModel` protocol and the HighSpec `GratingModel` (the notebook's grating
  equation plus a distortion polynomial without constant or linear term; `pixel` inverts the
  distortion by Newton steps), `model_from_config` / `model_from_solution`, `fit_model` (the
  notebook's `fit_solution` with `poly_powers` and the `alpha` / `B` bounds) and `anchor_model`;
  `load_murphy_atlas`, `remove_pedestal`, `convolve_to_resolving_power`, `flag_usable_lines`,
  `synthetic_atlas` and `catalogue_arrays`. The `mgb518` line table is regenerated by
  `scripts/make_line_table.py` from these: 31 of 73 lines usable (514.6056 nm Th I now passes
  the SNR 20 cut at 20.4 with the notebook's truncated Gaussian kernel; the earlier 30 came
  from a `gaussian_filter1d` kernel that put it at 19.996).
- Full-arc wavelength solve (`mastpipe.steps.wavecal`): `identify_lines` (nearest catalogue
  line within a tolerance, runner-up rejection), `shift_to_reference` (cross-correlation of
  continuum-subtracted fluxes with a parabolic peak), `solve_fiber` (the notebook's identify ->
  fit -> clip passes with the convergence rule and the `min_lines` give-up), `solve_full_arc`
  (every fiber seeded from one anchored model slid by the cross-correlation shift, or from its
  own solution in a previous reference), `build_reference` (identifier, `valid_from`, 12-hex
  SHA-256 over the solutions), `fiber_agreement` (the same line, fiber to fiber) and
  `wavelength_map`; session-scoped calibrated-full-arc, line-table and flux test fixtures.
- Provenance and FITS products (`mastpipe.pipeline.provenance`, `mastpipe.io.products`):
  `Provenance` as `MP*` header cards (`MPVER`, `MPCFGH`, the three statuses, `MPWREF`/`MPWREFH`,
  numbered `MPIN0001`/`MPCAL001`/`MPDARC01`/`MPSTEP01` lists with count cards) and back;
  `product_name`/`reference_name`/`combined_name`; writer/reader pairs for masters (image +
  `VARIANCE`, validity as `MPVALKnn`/`MPVALVnn`), trace sets (`TRACES`, `SAMPLES`), exposure
  spectra (`SPECTRA` with kept versions as `FLUX_<V>`/`ERROR_<V>`/`MASK_<V>` columns and
  presence cards, `FrameMeta` as `MPMKnn`/`MPMVnn`, optional `DRIFT` table), reference
  solutions (`SOLUTIONS`, `LINES`; hash and geometry checked on read) and combined spectra
  (`COMBINED`, `CONTRIB`, `EXCLUDED`, barycentric as `MPBCKnn`/`MPBCVnn`); `write_qc`/`read_qc`
  JSON. No raw instrument keyword reaches a product header.
- `identify` and `solve` (`mastpipe.pipeline.identify`, `mastpipe.io.qa`, the `mastpipe identify`
  command): one shared frame-to-lines chain (masters from `out/calib` when present, traces from
  the flat master or the arc itself), the notebook's two-panel anchor figure with numbered
  candidates (`anchor_panels.png`, `anchor_candidates.csv`; atlas template with `--atlas-dir`,
  synthetic otherwise), `run_identify` solving every fiber from anchor pairs (given, from a YAML
  file, or stored with the instrument, whose `trace:` is read by `load_anchor_trace`) and
  writing the reference product with `fiber_agreement.csv`, and `run_solve` seeded by a previous
  reference. matplotlib runs on the `Agg` backend.
- Drift (`mastpipe.steps.drift`, `mastpipe.pipeline.arcs`): `fit_fiber_drift` (comb-shift
  guess, one-to-one match within `tol_px`, clipped linear shift + stretch fit against the
  reference's identified lines, stretch zeroed below `stretch_significance` sigma),
  `fit_drift` (polynomial in `y_center` across the usable arc fibers into a `DriftModel`,
  a constant for one fiber), `apply_drift` (weighted shift and stretch undone per fiber
  before the reference model gives the wavelength), `pair_arcs` (`nearest_before`,
  `nearest`, `interpolate` within `max_age_hours`, same configuration, usable, not before
  the reference's `valid_from`) and `status_for`.
- Sky subtraction (`mastpipe.steps.sky`): `subtract_sky` (the telescope's sky fiber scaled by
  the throughput ratio, linearly resampled in wavelength onto the science grid and subtracted
  into `versions["skysub"]`, the raw extraction kept and the resampled sky under `Spectrum.sky`)
  and `subtract_sky_all` (science/sky pairing per telescope, with the fibers left raw listed by
  reason: no sky fiber, sky fiber fully masked, no wavelength).
- Within-exposure combination (`mastpipe.steps.combine`: `common_grid`, `combine_spectra`
  with include/exclude lists, optional S/N, flux-ratio and drift rejection listed in
  `excluded`, inverse-variance weighted mean with one-pass clipping, every shared version
  combined separately with the raw extraction under `"raw"`, `group_by_target`) and the
  flux-calibration and telluric-correction interfaces (`mastpipe.steps.fluxcal.apply_response`
  / `flux_status`, `mastpipe.steps.telluric.correct_telluric` / `telluric_status`) with
  pass-through defaults that add nothing and leave the statuses `"uncalibrated"` /
  `"uncorrected"`.
- Calibration library, quick-look runner and commands (`mastpipe.pipeline.callib`,
  `mastpipe.pipeline.runner`, `mastpipe.pipeline.commands`): `CalibLibrary` keeps masters and
  reference solutions under one folder with an `index.json` (`add`/`find` by camera, binning,
  readout mode, configuration, temperature tolerance and nearest time; `add_reference` with the
  camera, `find_reference` by the latest `valid_from` at or before a time or by identifier,
  `reference_for_source`); `identify`/`solve` register their reference in the output folder's
  library. `run_quicklook` reduces a night per `(camera, configuration)` group: masters (dark and
  background stacked one exposure time at a time, library fallback), traces from the flat or a
  full arc, arcs classified by lit traces, full arcs re-solved from the library's reference,
  drift arcs fitted against the reference valid at their time, science frames calibrated with
  cosmic rays flagged, extracted, given `transferred`/`reference`/`nominal` wavelengths by the
  time policy, sky-subtracted, passed through the flux and telluric stages, written with
  provenance, drift table and QC JSON, and combined per target; idempotence key (input checksum,
  configuration hash, reference identifier, sorted drift-arc checksums) with `--force` and the
  re-reduce of products whose key changed (`reduced_again`); `unknown` frames skipped with their
  classification rule; `--on-error continue|raise`; `out/logs/quicklook_<stamp>.log`. Commands
  `mastpipe quicklook`, `mastpipe solve` (exit 2 with "run identify first") and `mastpipe show`.
  `_arc_lines` also accepts a scanned `FrameRef`, so manifest overrides reach the arc chain.
  pytest now runs with `filterwarnings = error` and a registered `slow` marker.
- Quick-look: drift models are tied to the reference they were fitted against and only
  those are paired with a frame; the `interpolate` blend is tested and the per-arc blend
  weights are written to the science QC JSON (`drift_weights`).
- Regression on the laboratory data (`tests/regression/test_lab_highspec.py`, `slow`,
  skipped unless `MASTPIPE_TEST_DATA` points at the frames): the lab3 calibration frames
  and the best-focus Th-Ar frame staged with a manifest; trace finding on the calibrated
  flat reproduces the notebook's three trace centres within 0.3 px with quadratic rms
  below 0.01 px; `run_identify` with the stored anchors gives three solutions at rms
  0.016 / 0.033 / 0.012 px with 23 / 25 / 20 lines (the notebook: 0.016 / 0.038 / 0.010,
  23 / 25 / 19), coverage 513.50–527.32 nm within 0.03 nm and a median fiber-to-fiber
  spread of 0.023 px; the arc's own lines fitted for drift against its reference give
  zero shift; `run_quicklook` before and after `identify` fails on nothing. README:
  commands with examples, the three-cadence wavelength scheme, the product layout and
  the spectra table columns, the environment variables. DECISIONS: the rulings of Plan 1
  not yet recorded.

