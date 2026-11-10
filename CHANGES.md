# Changes

## 0.0.1 (unreleased)

- Provenance, QA and messages (2026-11-11):
  - a frame given the run's `--configuration` (its header and manifest name no grating)
    says so: its science or reference product gains the provenance step `configuration
    <name> from the run's --configuration (none in the header or manifest)`, and the QC
    of every science frame, full arc and drift arc gains `configuration_origin` (`run`
    or `frame`). Neither is part of an idempotence key or a reference hash, so no
    existing product is reduced again; a product keeps its old steps until it is
    rewritten for another reason;
  - a group whose frames name no configuration is labelled `group main/(no
    configuration)` rather than `group main/None`;
  - the anchor-panel PNG of `identify` is written atomically, like the FITS and JSON
    products: a temporary next to it moved into place, so a failed drawing leaves no
    partial file;
- Ingest, time and metadata fixes (2026-11-11):
  - a header-map `timezone` naming UTC by one of its IANA names (`UTC`, `Etc/UTC`, `UCT`,
    `Universal`, `Zulu` and their `Etc/` forms; `GMT` is left alone) is read as `utc`, and its frames are labelled `utc` rather than `UTC` /
    `Etc/UTC`. The label is part of the science idempotence key, so a science product
    reduced under such a custom map is reduced again once (its `MPTZ` becomes `utc`); the
    shipped HighSpec map uses `assume-utc`, whose labels (`utc` for a stamp with an
    offset, `assumed-utc` otherwise) do not change, so no HighSpec product is reduced
    again. The label logic lives in `core.times` (`read_time`, `zone_label`, and
    `ASSUMED_UTC`, still importable from `io.ingest`); a frame's stamp is parsed once;
  - the simulator reads `night_start` with `core.times.to_utc`: a start with an offset is
    converted to UTC instead of being written into the truth and headers with its offset;
  - a manifest `source` with `\` folder separators (written on Windows) now matches;
  - `scan_folder(manifest=..., unmatched=[...])` reports the manifest rows matching no
    frame, as `apply_manifest` does (the runner's own path is unchanged);
  - a manifest row naming an existing file that is not a FITS file is skipped as
    `manifest row N (name): not a FITS file (.fits, .fit, .fits.fz)` rather than as
    matching no frame (`io.ingest.names_non_fits_file`; `apply_manifest` reads no file);
    a row matching no frame is skipped as `manifest row N (name): matches no frame`
    (was `manifest row N: name: matches no frame`);
  - a manifest row that makes a calibration frame a science (or any non-calibration)
    frame without a target lists `targets` among its unresolved fields, as ingest does
    for such a header;
  - a multi-extension file of another instrument is skipped once, not once per HDU;
  - `OBJECT = 0` is the target `0`: only a missing card, a blank string or a logical card
    (`OBJECT = F`, likewise an `OBJnn`) names no target;
  - a calibration frame (bias, dark, flat, arc, background, from its header or a manifest
    row) no longer lists `targets` among its unresolved fields; a header map whose
    `targets` field is not `per_unit` is rejected when it is loaded;
  - a `FrameMeta` naming one telescope twice in `targets` is rejected (ingest and the
    manifest never produce one; a hand-built one made `target_for` take the last entry);
  - `calibrate_frame` raises on the dark path for a hand-built `CalibrationSet` with a
    `dark_rate_variance` but no `t_dark`, instead of dropping the dark variance;
  - `vocab.is_unit_name` and `unit_number` accept ASCII digits only and take the same
    optional `project` argument (default `vocab.PROJECT`, `mast`) as `unit_name`.
- Internal restructuring, no behaviour change: the quick-look runner is split into stage
  modules (`pipeline.group`, `masters`, `traces`, `arc_stage`, `science`); the public names
  stay importable from `mastpipe.pipeline.runner`.
- Masters and library fixes (2026-11-10):
  - a library flat is reused only when the masters it was calibrated with have the same
    content (the hash of their file checksums, now recorded in the flat as `MPMKH`), not
    merely the same file names; a flat written before this change carries no `MPMKH` and
    is built again once, on the first run that would have reused it. The rebuilt flat's
    file checksum changes, and with it the masters hash (`MPMKH`) of the science key, so
    every science product of an existing output folder is reduced again once;
  - the flat batch (one exposure time) is chosen among the frames that survive loading
    and calibration, as for the other kinds, not by its raw frame count; every flat batch
    is now loaded and calibrated, so an unusable frame of a batch that is not chosen is
    recorded in `failed` (as for darks and backgrounds);
  - a trace set built on a flat is named with the flat's first input checksum instead of
    `nosum`, so two trace sets of one stamp no longer share a name; an output folder's
    earlier `calib/*_nosum_traces.fits` and `qa/*_nosum_traces.json` files stay beside the
    new ones (the rename alone reduces no science product again: the trace set's name is
    not part of the science key).
- Review fixes (2026-11-08):
  - a background master whose inputs carried no exposure time no longer fails every
    frame of the group with "background subtraction requested ..."; the frames take the
    bias/dark path and `identify` / `solve` record why the background was not subtracted
    (`CalibrationSet.has_background`);
  - QC JSON files are strict JSON: NaN and infinities (a drift fiber without lines, say)
    are written as `null`, not the bare `NaN` token;
  - a manifest row that sets a field removes it from the frame's `unresolved` list, so
    product headers no longer report a manifest-supplied `mid_time`, `targets` or
    `arc_kind` as unknown;
  - the relative throughput and the arc lit-test sum each trace over its layout
    `aperture` when the layout sets one, the rows extraction sums, instead of the global
    `extract.halfwidth`.
- **Plan 1.5, hardening at the boundaries — upgrade notes.** These changes break
  compatibility with Plan 1 outputs and invocations:
  - the Mg b grating is `Mg` (the `MAST_common` name), no longer `mgb518`; `Halpha` and
    `Ca` likewise. `--configuration`, `--set`, grating files (`gratings/Mg.yaml`,
    `lines/Mg.csv`, `anchors/Mg_lab3.yaml`) and product names use the new names; the old
    spellings are read only from laboratory file names;
  - telescopes are `mast01`…`mast20` in layouts, products and QC;
  - `FrameMeta.target` is replaced by `FrameMeta.targets`, a mapping of telescope to
    target (`"*"` for every fiber);
  - product names gain a six-character input-checksum segment after the time stamp
    (`highspec_Mg_main_20261001T221500_3f9a1c_spec.fits`) and the new configuration
    name, so a Plan 1 product is never found at its new path: it is not reduced again
    in place, a new file is written beside it and the old `spec/` and `qa/` files remain
    (a glob over the folder then sees both). Use a fresh `--out` for a Plan 1 output
    folder, or clear its `spec/` and `qa/` folders first; the calibration library
    (`calib/`) may be kept, since its old entries (keyed `…/unknown` for the readout)
    stay in the index but never match, and the masters are rebuilt beside them.
    `quicklook` warns once per run when `spec/` holds products named the old way;
  - products, masters and library entries of a night whose frames name no grating
    (no `GRATING` card, no laboratory file-name alias), run with `--configuration Mg`,
    are now named and indexed `Mg`, not `noconf`, so a later night, `identify` and
    `solve` find that night's masters;
  - `identify` and `solve` take their masters from the calibration library (`--library`,
    `OUT/calib` by default) by camera, binning, readout, temperature and time, the flat
    by configuration too, as `quicklook` does, rather than the latest file of the
    folder; a master no longer in the index, or of another readout or temperature, is
    not used;
  - `SimTruth.science_flux` and `sky_flux` are per trace (`(n_traces, n_pix)`), as are
    `line_pixels`;
  - the configuration hash and the masters' calibration hash change, and the idempotence
    key gains parts: after upgrading, existing masters are rebuilt, and a product written
    under the new names is reduced again once when its key changes;
  - `tzdata` is a new runtime dependency (the IANA zones on Windows);
  - `quicklook` needs `--configuration` when no header card or file name names the
    grating (no silent default; the run fails with the list of the instrument's
    configurations);
  - the `FIBPRE` card is no longer read (in MAST a fiber preset is a stage position, not
    a fiber layout): a survey night that relied on `FIBPRE` now needs `--layout
    survey45` (or a manifest's `fiber_preset`).
- Fleet vocabulary: `mastpipe.vocab` mirrors `MAST_common`'s grating, disperser,
  instrument and unit names and is what the code uses; a `fleet` CI job checks
  `MAST_common` out as a sibling (paired branch, else `master`) and a test fails when the
  mirror diverges (`MASTPIPE_FLEET_CHECK=1`). Grating ids are case-sensitive; the
  registry refuses a grating block whose name is not a fleet grating name.
- Time: one helper, `mastpipe.core.times`, parses every ISO stamp. Header-map time fields
  carry a `timezone` option (`utc`, an IANA zone, `assume-utc`); a stamp with an offset is
  honoured. HighSpec reads `DATE-OBS`, then `FRAME`, then `DATE`, under `assume-utc`,
  which logs one warning per run and writes `MPTZ = 'assumed-utc'` on every product.
  Mid-times stay naive ISO UTC strings.
- Ingest: night folders are scanned recursively (`<night>/highspec/acquisition-NNNN/
  exposure-NNN.fits`), a frame's `source` is its path relative to the night folder, and
  files without a FITS suffix are not scanned; a FITS file that cannot be opened is a
  failure. Frames whose `INSTRUME` names another instrument are skipped with the reason.
  Manifest rows name frames by relative path, glob pattern (`**` crosses folders) or
  unique file name (an ambiguous name is an error), and rows matching nothing are
  reported among the skipped.
- Frame type: `IMAGETYP = 'light'` maps to the provisional type `light`, refined at ingest
  by any conditional header-map rule (zero exposure → bias, lamp cards) and otherwise
  resolved by the runner from the frame's content (`steps.classify.classify_light`,
  `ClassifyParams`) into `arc` (`content:arc_like`) or `science` (`content:continuum`), or
  skipped (`content:no_signal`), before the trace stage; a manifest's `exposure_type`
  wins. `FIBPRE` is no longer read: the layout comes from `--layout`, a manifest's
  `fiber_preset` or the instrument's default.
- Targets per fiber: a header-map field kind `per_unit` (`OBJ{nn}`, one card per
  telescope unit), a manifest column `targets` (`mast01=HD 1234;mast02=…`; the old
  `target` column is still read), and a single `OBJECT` as the fallback for every science
  fiber (an empty or blank `OBJECT` names none); combination groups by target across
  fibers.
- Runner robustness: master building loads calibration frames one at a time, and a frame
  that cannot be read, is binned otherwise than the group's science frames, has another
  shape than the detector's (crop and orientation, scaled by binning:
  `io.ingest.expected_shape`) is dropped and recorded as a failure (exit status 2); the
  night is still reduced. A calibration frame of another readout (a full-arc trace
  candidate included) is left out and recorded in `skipped`, not as a failure; a
  science, `light` or arc frame of another binning or readout fails. A failing full-arc
  candidate is recorded and skipped. Products, masters, references, indexes and QC files are written
  atomically (`io.atomic.atomic_write`, a `<name>.tmp-<pid>` temporary moved into place),
  and the calibration library's index is read and updated under a lock on `index.lock`
  (`io.atomic.locked`).
- Idempotence: a science product's key covers the frame's checksum, the configuration
  hash (time policy and every override included), the reference solution's hash
  (`MPWREFH`), the drift arcs' checksums, the layout preset and its file's hash
  (`MPLAYH`, computed with CRLF read as LF, and `.gitattributes` keeps the CSV files LF),
  the masters' checksums (`MPMKH`) and the time zone (`MPTZ`). A master the
  library already holds from the same frames, stacking parameters and detector
  (`calib_hash`) is reused instead of rebuilt; `--force` rebuilds.
- Fibers are keyed by the layout's `trace` column, not by the row position, so a layout
  with shuffled rows keys the same physical fibers.
- Simulator: every fiber has its own dispersion (`dispersion_jitter_px`, 0.3 px by
  default) and the drift varies with the row (`drift_gradient_px`, 0.05 px); `swap_fibers`
  exchanges two traces' light while the truth stays unswapped; `SimTruth` gains
  `x_c_per_trace`, `B_per_trace`, `row_of(trace_id)` and `trace_shifts(file)`. New tests
  prove the comparisons are sensitive: an end-to-end swap test and a drift fit against a
  2 px y-gradient. `simulate_night(site_layout=True)` writes a night in the site's shape
  (`highspec/acquisition-NNNN/exposure-NNN.fits`, `seq.txt`, the cards the control
  software writes, `IMAGETYP = 'light'` for arcs, science frames and the background, no
  grating or fiber-preset card), which quick-look reduces with `--configuration Mg
  --layout survey45` and no manifest (without `FIBPRE` nothing else selects the survey
  layout; the background, which content cannot type, is skipped); the laboratory
  regression's manifest types only what content cannot (the flat and the background).
- Numerics: the master variances enter the calibrated frame's variance with exact
  propagation (the bias correlated through both subtractions); lines flagged saturated
  are excluded from the wavelength fit; a drift model is usable only with at least
  `drift.min_fibers` (3) usable arc fibers, recorded in the QC and `MPDUSE`; a combined
  spectrum keeps every version some contributor carries, listing the fibers without it
  in `MPVMISS`.
- Readout: `FrameMeta.readout_key` combines `READMODE` and the output amplifier
  (`OUTPTAMP`, `image/em`), the calibration library matches masters on it, and the
  detector configuration gains an optional `readouts` table of gain and read noise per
  amplifier (one warning per run for an amplifier missing from it).
- Tooling: the fleet's ruff rules (`F,E,W,I,N,B,UP,C90,SIM,TID` plus `NPY`), line length
  125, `ruff==0.16.0` pinned in the `dev` extra, `ruff format --check` before
  `ruff check`; CI gains a separate `lint` job and a concurrency group per ref.
- Documents: the design documents under `docs/design/` are published, with roles in
  place of names; the night-folder guide describes the site's folder shape, the
  provisional `light` type, `--configuration`, the time-zone option, manifest patterns
  and per-fiber targets; README and install pages say the PypeIt adapter is planned and
  does not exist yet.
- The documentation is hosted on Read the Docs at <https://mastpipe.readthedocs.io/>
  (project `mastpipe`, built from `.readthedocs.yaml`); the README links to it.
- The repository is public under `The-MAST-project/mastpipe`. Documentation, comments
  and test identifiers refer to the HighSpec laboratory reduction notebook rather than
  to individuals; the repository `CLAUDE.md` records the public visibility.
- Install documentation (`docs/running/install.md`, README) separates a plain install
  (`pip install "git+<repository>"` into a venv or conda environment) from the developer
  install (clone, `pip install -e ".[dev]"`), and says the tests run from the root of the
  clone or with the `tests/` path, since the suite is not installed with the package.
- Tutorial notebooks under `docs/tutorials/` (committed without outputs, executed by the
  documentation build through `myst-nb`): *The pipeline step by step* walks the Python API on
  a simulated night with a hand-made nine-fiber layout, from raw frames to a written product;
  *Inspecting a night's products* reduces a simulated survey night with `identify` and
  `quicklook` and reads back the spectra, drift, combined spectra, reference, QC and log.
- `io.qa` draws the anchor panels on an `Agg` canvas through matplotlib's object-oriented
  interface instead of switching the global backend, so importing the pipeline no longer
  disables inline figures in a notebook.
- Documentation site (`docs/`, Sphinx with MyST and Furo; `pip install -e ".[docs]"`,
  `sphinx-build -W -b html docs docs/_build/html`): an observer guide (install, night
  folder and manifest, first reference, quick-look, products, troubleshooting), a
  description of the method (stages, three-cadence wavelengths, traces and extraction,
  sky and combination, products policy, limits), a developer guide (architecture, data
  model, adding a grating or an instrument, simulator and tests, conventions) and an
  autodoc API reference; `.readthedocs.yaml`; a CI job building with warnings as errors;
  `scripts/make_doc_figures.py` producing the two committed figures from the simulator.
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
  `Mg` grating block, `lab3` and generated `survey45` layouts, the Mg b line table
  (73 Murphy et al. catalogue lines, 30 usable at R = 23081) and the `Mg`/`lab3`
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
  `synthetic_atlas` and `catalogue_arrays`. The `Mg` line table is regenerated by
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

