# mastpipe Plan 1 — HighSpec end to end: raw frames to sky-subtracted, wavelength-calibrated spectra

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A pip-installable `mastpipe` package that takes a night of HighSpec frames to per-fiber spectra that are cleaned, extracted, wavelength-calibrated (reference solution plus calibration-fiber drift), sky-subtracted and combined across telescopes, with every earlier version kept, FITS products with provenance, a manual `identify` command for the first reference solution, and `mastpipe quicklook` running the whole chain on a simulated night and on the notebook's laboratory frames.

**Architecture:** Three strictly separated layers (September spec §3): declarative instrument configuration; pure step functions on arrays and dataclasses; pipeline and I/O (ingest, calibration library, runner, products, CLI). The HighSpec laboratory reduction notebook functions are ported one-to-one into the step layer (calibration, traces, extraction, continuum, line finding, grating-equation wavelength solution). Wavelengths follow the three-cadence scheme of revision 1 §2: a stored reference solution per fiber from a full arc, a linear drift per exposure from calibration-fiber arcs, paired by time. A simulator produces synthetic nights with known truth so every step and the end-to-end path are tested before survey data exist.

**Tech Stack:** Python ≥ 3.11 · numpy · scipy · astropy · pandas · matplotlib · pydantic v2 · pyyaml · pytest · ruff · GitHub Actions (Linux/Windows/macOS matrix).

**Spec:** `docs/design/2026-09-24-mastpipe-design.md` as amended by `docs/design/2026-10-29-mastpipe-design-revision-1.md` (same folder as this plan; the revision wins where they differ). Draft repository `CLAUDE.md`: `docs/design/2026-09-24-mastpipe-claude-md-draft.md`. This plan supersedes `docs/design/2026-09-24-mastpipe-plan-1-foundation.md`, whose data-model, configuration, ingest, simulator, calibration, trace, extraction and product designs it carries forward with the changes listed in revision 1 §7.

## Roadmap (three plans)

| Plan | Delivers |
|---|---|
| **1 — HighSpec end to end (this plan)** | skeleton + CI + privacy scan; data model; configuration; HighSpec definition with the Mg b line table; header-map ingest with arc sub-classification; simulator with full and calibration-fiber arcs and drift; cleaning (masters, variance, cosmic rays); traces + layout matching; aperture extraction + continuum; line finding; grating wavelength model; full-arc solve from a previous reference; `identify`; drift fit, transfer and time policy; sky subtraction; within-exposure combination; flux/telluric pass-through stages; products with kept versions; calibration library with time-valid references; `quicklook` with re-reduce; regression on the notebook's lab data |
| 2 — DeepSpec and better extraction | DeepSpec definition (orientation, per-camera files, bands with Kr/Th-Ar catalogues), `dispersion` model, Kr-then-ThAr `identify`, regression against the DeepSpec instrument team's solutions, DeepSpec science path, Horne extraction, across-exposure combination, DeepSpec drift once answered |
| 3 — Operations | QA figures and pages, `watch` mode, focus recipe, PypeIt adapter, flux-calibration and telluric recipes once the observing policy exists |

## Global Constraints

- Python `>=3.11`; hard dependencies only: numpy, scipy, astropy, pandas, matplotlib, pydantic (v2), pyyaml. No new top-level dependency without asking the user.
- Pure Python, cross-platform: `pathlib` only, no shell-outs, no POSIX-only calls, explicit `encoding="utf-8"` on every text open, `newline=""` on CSV writes.
- **No local absolute paths, hostnames, IP addresses, usernames, e-mail addresses, tokens** in any tracked file (code, tests, docs, notebooks, CI). Real test data are located via the environment variable `MASTPIPE_TEST_DATA`; the Murphy UVES Th-Ar atlas via `MASTPIPE_ATLAS_DIR`; neither enters git.
- Tests use synthetic data only (the simulator or tiny hand-made arrays), fixed seeds, no network. Tests needing real data skip when the variable is unset.
- Step functions (`src/mastpipe/steps/`) never read files, never touch configuration objects, never log to files; they take arrays/dataclasses/plain parameters and return dataclasses. Figures are produced only by `io/qa.py` (this plan: the `identify` anchor panels) and never by step functions.
- No raw FITS keyword name appears outside `src/mastpipe/io/ingest.py` and the instrument `header_map.yaml` files. Raw headers are never copied wholesale into products or logs.
- The fiber layout is an input table; the number of traces is never inferred. The wavelength model kind is set by configuration, never inferred.
- Air wavelengths in nm internally; barycentric corrections stored, never applied.
- Science exposures carry no arc light. The wavelength of a science fiber is the reference solution corrected by drift transferred from the calibration fibers (`transferred`), else the reference alone (`reference`), else geometry (`nominal`).
- A reference solution is never edited in place; it is valid from its full arc's mid-time until the next full arc of the same configuration and preset.
- Every per-exposure product keeps the raw extraction and every later version as separate columns on one pixel grid; every product header carries pipeline version, configuration hash, input checksums, calibrations used, reference-solution identifier and hash, drift-arc checksums, `wavelength_status`, `flux_status`, `telluric_status`.
- `src/` layout, package `mastpipe`, console script `mastpipe`.
- Git: local repository only. **No remote, no push** — the GitHub repository is created only when the user says so. Commit after every task; `docs/design/` and `spike/` are git-ignored.
- Ruff for lint and format (line length 100). CI must be green on Linux, Windows and macOS before a task is considered done.
- Writing style in code comments and docs: no first-person AI voice. `DECISIONS.md` records design decisions as they are made.

## Review Focus

1. **A night with no calibration-fiber arc at all** — quick-look must still produce every science product with `wavelength_status = "reference"` and say so in the QC file, never crash or skip the frames (Task 18 test `test_quicklook_without_drift_arcs_gives_reference_status`).
2. **A science frame whose sky fiber is missing or fully masked** — the science fiber keeps its raw extraction, gets no `skysub` version, and the QC lists the telescope and the reason (Task 15 test `test_missing_sky_fiber_leaves_raw_only`).
3. **A drift arc from another configuration or older than `drift.max_age`** — it must never be paired with the exposure (Task 14 test `test_pairing_ignores_other_configuration_and_stale_arcs`).
4. **Saturated or blended arc lines** — excluded from the solve and from the drift fit through `use_x`, so one hot line cannot pull the solution (Task 10 test `test_saturated_line_not_usable`).
5. **A frame with no exposure time in the header or manifest** — classified `unknown`, skipped with a recorded reason; never scaled by a guessed time (Task 5 test `test_missing_exposure_time_is_unknown`, Task 18 test `test_unknown_frames_are_skipped_with_reason`).

## File structure (this plan)

```
mastpipe/
  pyproject.toml  CLAUDE.md  README.md  CHANGES.md  DECISIONS.md  LICENSE  .gitignore
  .github/workflows/ci.yml
  scripts/privacy_scan.py  scripts/make_line_table.py  scripts/make_survey45_layout.py
  src/mastpipe/
    __init__.py  __version__.py  cli.py
    core/model.py           # dataclasses (Task 2)
    config/models.py  config/loader.py
    instruments/__init__.py  registry.py
    instruments/highspec/config.yaml  header_map.yaml
    instruments/highspec/layouts/lab3.csv  survey45.csv
    instruments/highspec/gratings/mgb518.yaml
    instruments/highspec/lines/mgb518.csv          # Murphy catalogue lines in band + usable flag
    instruments/highspec/anchors/mgb518_lab3.yaml  # the notebook's two anchor pairs
    io/headers.py  io/ingest.py  io/products.py  io/qa.py
    steps/geometry.py  steps/calib.py  steps/cosmics.py  steps/trace.py  steps/extract.py
    steps/lines.py  steps/atlas.py  steps/wavemodel.py  steps/wavecal.py  steps/drift.py
    steps/sky.py  steps/combine.py  steps/fluxcal.py  steps/telluric.py
    sim/simulator.py
    pipeline/provenance.py  pipeline/callib.py  pipeline/arcs.py  pipeline/runner.py
  tests/conftest.py  tests/unit/test_*.py  tests/e2e/test_quicklook.py  tests/regression/test_lab_highspec.py
```

---

### Task 1: Repository skeleton, packaging, CI, privacy scan

**Files:**
- Create: `pyproject.toml`, `src/mastpipe/__init__.py`, `src/mastpipe/__version__.py`, `src/mastpipe/cli.py`, `CLAUDE.md`, `README.md`, `CHANGES.md`, `DECISIONS.md`, `LICENSE`, `.github/workflows/ci.yml`, `scripts/privacy_scan.py`, `tests/conftest.py`, `tests/unit/test_cli.py`
- Existing: `.gitignore` (keep; it ignores `docs/design/` and `spike/`)

**Interfaces:**
- Produces: console script `mastpipe --version`; package importable; `scripts/privacy_scan.py` exits 1 on a hit.

- [ ] **Step 1: Initialise git and write `pyproject.toml`** — `git init -b main` in the repository root. Metadata: name `mastpipe`, BSD-3-Clause, `requires-python >= 3.11`, dynamic version from `mastpipe.__version__.__version__`, dependencies `numpy>=1.26 scipy>=1.11 astropy>=6.0 pandas>=2.0 matplotlib>=3.8 pydantic>=2.5 pyyaml>=6.0`; extras `dev = [pytest>=8, pytest-cov, ruff>=0.6, nbstripout]`, `pypeit = [pypeit>=2.0]`; script `mastpipe = "mastpipe.cli:main"`; setuptools `src` layout with package data `instruments/**/*.yaml`, `instruments/**/*.csv`; ruff line length 100, target py311, select `E F I UP B NPY`; pytest `testpaths = ["tests"]`.
- [ ] **Step 2: Package files** — `__version__ = "0.0.1"`; `cli.py` with `build_parser() -> argparse.ArgumentParser` (prog `mastpipe`, `--version`, sub-parsers under `dest="command"`) and `main(argv: Sequence[str] | None = None) -> int` printing help and returning 0 when no command. `CLAUDE.md` = the draft without its italic banner. `README.md`: purpose, `pip install -e ".[dev]"`, `pytest`, the two environment variables and the expected layout under `MASTPIPE_TEST_DATA` (`highspec/lab3/calib/`, `highspec/lab3/scan/`) and `MASTPIPE_ATLAS_DIR` (`thar_spec_MM201006.fits`, `thar_MM201006.dat`). `CHANGES.md` with `## 0.0.1 (unreleased)`. `DECISIONS.md` with two dated entries: 2026-09-24 (independent package; soxspipe rejected; PypeIt deferred) and 2026-10-29 (the notebook's engine confirmed; three-cadence wavelength scheme; HighSpec-first with the whole path in Plan 1; six stages; kept versions; the eight calls of revision 1 §8). `LICENSE` BSD-3-Clause, holder "The MAST project".
- [ ] **Step 3: Privacy scan** — `scripts/privacy_scan.py`: patterns `local path` (the Linux mount and home roots, the macOS and Windows user folders), `e-mail`, `ipv4` (allow `0.0.0.0`, `127.0.0.1`), `token` (`ghp_…`, `AKIA…`, `sk-…`); scans `git ls-files` with text suffixes; skips itself; prints hits; exit 1 on any.
- [ ] **Step 4: CI** — `.github/workflows/ci.yml`: matrix `ubuntu-latest, windows-latest, macos-latest` × Python `3.11, 3.12`; steps: checkout, setup-python, `pip install -e ".[dev]"`, `ruff check src tests scripts`, `ruff format --check src tests scripts`, `python scripts/privacy_scan.py`, `pytest`.
- [ ] **Step 5: Failing tests** — `tests/unit/test_cli.py::test_version_flag` (SystemExit code 0, version in stdout), `::test_no_command_prints_help` (returns 0, `usage: mastpipe` in stdout).
- [ ] **Step 6: Run** — `python -m pip install -e ".[dev]" && pytest tests/unit/test_cli.py -v` → 2 passed.
- [ ] **Step 7: Lint, scan, commit** — `ruff check src tests scripts && ruff format src tests scripts && python scripts/privacy_scan.py` → exit 0. `git add -A && git commit -m "chore: package skeleton, CLI, CI matrix, privacy scan"`.

---

### Task 2: Core data model

**Files:**
- Create: `src/mastpipe/core/__init__.py`, `src/mastpipe/core/model.py`
- Test: `tests/unit/test_model.py`

**Interfaces (used by every later task):**
- Constants: `EXPOSURE_TYPES = ("bias","dark","flat","arc","science","background","unknown")`; `ARC_KINDS = ("full","calfiber","unknown")`; `ROLES = ("arc","science","sky","unknown")`; `WAVELENGTH_STATUSES = ("none","nominal","reference","transferred","solved")`; `FLUX_STATUSES = ("uncalibrated","calibrated")`; `TELLURIC_STATUSES = ("uncorrected","corrected")`; `VERSIONS = ("skysub","fluxcal","telluric")`.
- `FrameMeta(instrument, camera, exposure_type, configuration=None, exposure_time=None, mid_time=None, binning=(1,1), readout_mode=None, detector_temperature=None, lamp_state=None, fiber_preset=None, target=None, arc_kind=None, checksum="", source="", extension=0, unresolved=(), classification_rule=None)` frozen; validates `exposure_type` and `arc_kind`; `mid_time` ISO-8601 UTC string.
- `FrameRef(path: Path, extension: int, meta: FrameMeta)`.
- `Frame(data, meta, variance=None, mask=None)` with `.shape`; `mask` bool, True = bad; shape mismatch raises `ValueError`.
- `MasterCalib(kind, data, validity: dict, n_frames, provenance: dict = {}, variance=None)`.
- `CalibrationSet(bias, dark_rate, background, t_background, gain, read_noise, shape)`.
- `Trace(index, coeffs, y_center, slope, rms, n_points, snr, x_samples, y_samples, kept, role=None, telescope=None, layout_index=None, aperture=None)` with `.y(x)` (numpy polyval order).
- `TraceSet(traces, shape, layout_shift=None, missing=(), extra=())` with `__len__`, `__iter__`, `by_role(role)`, `by_telescope(telescope) -> dict[str, Trace]` (role → trace).
- `SpectrumVersion(flux, error, mask)` — same length as the pixel grid.
- `Spectrum(fiber, pixel, flux, error, mask, role=None, telescope=None, target=None, wavelength=None, sky=None, continuum=None, versions: dict[str, SpectrumVersion] = {})`; `flux/error/mask` are the raw extraction and are never replaced; `versions` keys must be in `VERSIONS`; `.latest(upto: str | None = None) -> SpectrumVersion` returns the most advanced version present in `VERSIONS` order (or the raw arrays as a `SpectrumVersion` when none).
- `ExposureSpectra(meta, spectra, wavelength_status="none", flux_status="uncalibrated", telluric_status="uncorrected", reference_id=None, reference_hash=None, drift_arcs: tuple[str, ...] = (), provenance=None)`; validates the three statuses.
- `WavelengthSolution(fiber, kind, params: np.ndarray, n_pix, lines: pd.DataFrame, rms_px, n_lines, anchor_px: tuple[float,float], anchor_nm: tuple[float,float], trace_y_center, geometry: dict)`; `lines` columns `x, wavelength, species, dx_px`.
- `ReferenceSolution(instrument, configuration, fiber_preset, solutions: dict[int, WavelengthSolution], created: str, valid_from: str, source_checksums: tuple[str, ...], identifier: str, hash: str = "")`; `make_identifier(instrument, configuration, fiber_preset, created) -> str` = `"{instrument}_{configuration}_{fiber_preset}_{YYYYmmddTHHMMSS}"`.
- `DriftModel(arc_checksums: tuple[str,...], arc_mid_time: str, per_fiber: pd.DataFrame, shift_coeffs: np.ndarray, stretch_coeffs: np.ndarray, scatter_px: float, usable: bool, n_pix: int)`; `per_fiber` columns `fiber, y_center, shift_px, stretch, n_lines, rms_px, usable`; `.shift_at(y)`, `.stretch_at(y)` (polyval over trace `y_center`).
- `CombinedSpectrum(target, wavelength, versions: dict[str, SpectrumVersion], contributors: pd.DataFrame, excluded: pd.DataFrame, grouping: str, barycentric: dict[str, float])`; `contributors` columns `fiber, telescope, source, weight`; `excluded` columns `fiber, telescope, source, reason`.
- `QCReport(stage, values: dict)`.

- [ ] **Step 1: Failing tests** — `tests/unit/test_model.py`: `test_frame_meta_rejects_unknown_exposure_type`, `test_frame_meta_rejects_unknown_arc_kind`, `test_frame_rejects_shape_mismatch`, `test_trace_evaluates_polynomial` (coeffs `[0.001, 0.5, 10.0]` → `y(2.0) == 11.004`), `test_traceset_by_role_and_telescope`, `test_spectrum_versions_validate_and_latest` (a `Spectrum` with `versions={"skysub": v}` returns `v` from `.latest()`, raw from `.latest(upto=None)` when empty, raises on key `"banana"` and on a version of the wrong length), `test_exposure_spectra_statuses` (bad `flux_status` raises), `test_reference_identifier_format` (`make_identifier("highspec","mgb518","lab3","2026-09-09T12:34:56")` == `"highspec_mgb518_lab3_20260909T123456"`), `test_drift_model_evaluates_at_y`.
- [ ] **Step 2: Run** — `pytest tests/unit/test_model.py -v` → ImportError.
- [ ] **Step 3: Implement** the dataclasses above in `core/model.py` (frozen `FrameMeta`, `FrameRef`; `__post_init__` validation elsewhere).
- [ ] **Step 4: Run** → all pass. **Step 5: Commit** `feat(core): data model`.

---

### Task 3: Configuration models and layered loader

**Files:**
- Create: `src/mastpipe/config/__init__.py`, `src/mastpipe/config/models.py`, `src/mastpipe/config/loader.py`
- Test: `tests/unit/test_config.py`

**Interfaces:**
- Pydantic v2 models (all `extra="forbid"`):
  - `DetectorConfig(name, shape: tuple[int,int], gain: float, read_noise: float, pixel_size_um: float, saturation: float, crop: tuple[int,int,int,int] | None, dispersion_axis: int = 1, transpose: bool = False, flip_x: bool = False, flip_y: bool = False)`.
  - `GeometryConfig(model: Literal["grating","dispersion"], groove_spacing_nm, incidence_angle_deg, order, camera_focal_mm, central_wavelength_nm, dispersion_sign: int = 1)` (dispersion-model fields are added in Plan 2).
  - `CalibParams(subtract_background: bool = True, flat_subtract_background: bool = False, stack_method: str = "median", clip_sigma: float | None = None, cosmic_sigma: float = 8.0, cosmic_size: int = 5)`.
  - `TraceParams(center_width=100, trace_halfwidth=8, threshold_sigma=5.0, centroid_sigma=None, min_separation=None, baseline_window=None, x_bin=25, degree=2, clip_sigma=3.0, min_points=4, layout_tolerance=3.0, max_traces=None)`.
  - `ExtractParams(halfwidth=9.0, aperture: Literal["round","subpixel","floor"]="round", method: Literal["sum","mean"]="sum", min_coverage=0.99, continuum_degree=15, clip_above=2.0, clip_below=3.0, continuum_iter=10)`.
  - `LineParams(halfwidth=9, threshold_sigma=5.0, min_separation=2, x_halfwidth=6, y_halfwidth=None, blend_distance=None, isolation=None, edge_margin=10, baseline_window=51, sigma_guess=1.3, max_lines=None, max_components=4, add_sigma=5.0, max_y_offset=1.5, saturation=None)`.
  - `WavecalParams(resolving_power=20000.0, max_shift_fwhm=0.02, min_snr_atlas=20.0, usable_only=True, tol_first_px=2.0, tol_later_px=1.0, unique_tol_px=4.0, n_passes=3, converge_px=1e-3, clip_sigma=3.0, fit_alpha=False, poly_deg=3, min_lines=6, min_arc_fibers_lit=1.0)`.
  - `DriftParams(tol_px=1.5, max_shift_px=200, clip_sigma=3.0, min_lines=5, stretch_significance=3.0, degree=1, max_scatter_px=0.2, max_age_hours=12.0, time_policy: Literal["nearest_before","nearest","interpolate"]="nearest_before", warn_shift_px=2.0)`.
  - `SkyParams(scale_by_throughput=True, resample: Literal["linear"]="linear")`.
  - `CombineParams(across_telescopes=True, include: list[str] = [], exclude: list[str] = [], reject=False, reject_snr_min=3.0, reject_flux_ratio=0.3, reject_drift_px=None, clip_sigma=4.0, log_grid=True)`.
  - `StepParams(calib, trace, extract, lines, wavecal, drift, sky, combine)` with defaults.
  - `ConfigurationBlock(name, band_nm: tuple[float,float], resolving_power, geometry: GeometryConfig, line_table: str | None = None, params: dict = {})`.
  - `InstrumentConfig(name, cameras: dict[str, DetectorConfig], configurations: dict[str, ConfigurationBlock], params: StepParams, default_layout: str, layouts: list[str], arc_classification: Literal["header","lit_traces"] = "lit_traces")`.
  - `EffectiveConfig(instrument, camera, configuration, detector, geometry, band_nm, resolving_power, line_table, params: StepParams, layers: list[str])`.
- `deep_merge(base, override) -> dict` (override wins; dicts recurse; lists replace); `load_yaml(path) -> dict`; `build_effective(instrument, configuration, camera, overrides: list[dict] = ()) -> EffectiveConfig` (merge order: instrument params → block params → overrides); `config_hash(effective) -> str` (12 hex of sha256 over canonical JSON); `parse_override("a.b.c=1.5") -> dict`.

- [ ] **Step 1: Failing tests** — `test_deep_merge_override_wins_and_lists_replace`, `test_build_effective_layers_block_over_instrument`, `test_config_hash_stable_and_sensitive`, `test_unknown_key_rejected` (`extra="forbid"`), `test_parse_override_types` (`"x=1"`→int, `"x=1.5"`→float, `"x=true"`→bool, `"x=abc"`→str, `"x=[1,2]"`→list via JSON).
- [ ] **Step 2: Run** → fail. **Step 3: Implement.** **Step 4: Run** → pass. **Step 5: Commit** `feat(config): models and layered loader`.

---

### Task 4: Instrument registry, HighSpec definition, Mg b line table

**Files:**
- Create: `src/mastpipe/instruments/__init__.py`, `registry.py`, `highspec/config.yaml`, `highspec/header_map.yaml`, `highspec/gratings/mgb518.yaml`, `highspec/layouts/lab3.csv`, `highspec/layouts/survey45.csv`, `highspec/lines/mgb518.csv`, `highspec/anchors/mgb518_lab3.yaml`, `scripts/make_survey45_layout.py`, `scripts/make_line_table.py`
- Test: `tests/unit/test_registry.py`

**Interfaces:**
- `instrument_dir(name) -> Path`; `list_instruments() -> list[str]`; `load_instrument(name) -> InstrumentConfig` (merges every `gratings/*.yaml` into `configurations`; a new grating is one file, no code change); `load_header_map(name) -> dict`; `load_layout(name, preset) -> pd.DataFrame` with columns `trace:int, role, telescope, y_nominal, slope_nominal, aperture` sorted by `y_nominal`, roles validated; `load_line_table(name, configuration) -> pd.DataFrame` with columns `wavelength_nm, species, logI, usable` sorted by wavelength; `load_anchors(name, configuration, preset) -> list[tuple[float, float]]` (wavelength nm, pixel) or `[]`.
- HighSpec numbers (the HighSpec laboratory reduction notebook): detector `main` shape `[400, 1600]`, crop `[1, 1599, 1, 399]`, 16 µm pixels, gain/read-noise **placeholders** flagged in comments; grating `mgb518`: groove spacing 367.0 nm (2724.8 l/mm), α 45°, order 1, f_cam 479.41 mm, λc 518.0, band `[513.0, 527.5]`, R 20000, `model: grating`, `line_table: lines/mgb518.csv`; trace halfwidth 8, degree 2; extract halfwidth 9, continuum degree 15; line finding as in `LineParams` defaults with `sigma_guess = 3.0/2.3548`; `default_layout: lab3`.
- `lab3.csv`: three traces, roles `arc`, `y_nominal` 17.4 / 189.0 / 360.7 (cropped coordinates), slope 0.0017 / 0.0020 / 0.0023, telescope empty. `survey45.csv` generated by `scripts/make_survey45_layout.py`: 45 traces at pitch 8.5 px from y 12.0, order arc/science/sky interleaved per telescope `T01..T20` with the 5 arc fibers at indices 0, 11, 22, 33, 44.
- `lines/mgb518.csv`: Murphy et al. (2007) catalogue lines in `[513.0, 527.5]` nm (73 lines) with `usable` from `scripts/make_line_table.py --instrument highspec --configuration mgb518 --resolving-power 23081` run against `MASTPIPE_ATLAS_DIR` (criteria: convolved centroid within 0.02 FWHM, S/N > 20 — Task 11's `flag_usable_lines`). The script is committed; the CSV it produced is committed; the atlas is not. Header comment cites Murphy et al. 2007, MNRAS 378, 221 and 2010, ApJ 723, 89.
- `anchors/mgb518_lab3.yaml`: `trace: 0`, `pairs: [[515.8605, 267.22], [525.8360, 1423.22]]` (the notebook's anchors, air nm and cropped pixel).
- `header_map.yaml` (Andor Solis keywords): `exposure_time: {keys: [EXPOSURE, EXPTIME], type: float}`, `mid_time: {keys: [DATE-OBS, DATE], type: str}`, `detector_temperature: {keys: [TEMP, CCD-TEMP]}`, `binning_x: {keys: [HBIN]}`, `binning_y: {keys: [VBIN]}`, `readout_mode: {keys: [READMODE, OUTPTAMP], lower: true}`, `exposure_type: {keys: [IMAGETYP, OBSTYPE], lower: true, map: {bias: bias, zero: bias, dark: dark, flat: flat, arc: arc, comp: arc, object: science, light: science, science: science, background: background}}`, `arc_kind: {keys: [ARCKIND], lower: true}`, `lamp_state: {keys: [LAMP]}`, `configuration: {keys: [GRATING], lower: true, filename_regex: "(?P<value>mgb518|halpha|caiihk)"}`, `fiber_preset: {keys: [FIBPRE]}`, `target: {keys: [OBJECT]}`; classification rules in order: `exposure_time < 0.001 → bias`, `lamp_state == "off" and shutter closed → dark`, `lamp_state == "thar" → arc`, `lamp_state == "flat" → flat`, `exposure_time missing → unknown`, else `science`.

- [ ] **Step 1: Failing tests** — `test_highspec_is_registered`, `test_highspec_config_loads_with_grating_blocks` (shape `(400,1600)`, `mgb518` λc 518.0, `geometry.model == "grating"`), `test_layout_tables_parse_and_validate` (lab3: 3 rows; survey45: 45 rows, 5 arc, 20 science, 20 sky, monotonic `y_nominal`, telescopes set for non-arc), `test_line_table_loads` (73 rows, ≥ 25 usable, sorted, species strings like `"Th I"`), `test_anchors_load` (two pairs, first wavelength 515.8605), `test_header_map_yaml_has_required_sections`, `test_unknown_instrument` (KeyError), `test_extra_grating_file_is_picked_up` (temporary `halpha.yaml` copy is discovered).
- [ ] **Step 2: Run** → fail. **Step 3: Write data files and `registry.py`.** `make_line_table.py` imports Task 11's atlas functions; until Task 11 exists, commit the CSV produced by a throwaway run of the notebook criteria, and re-generate it in Task 11 (the diff must be empty). **Step 4: Run** → pass. **Step 5: Commit** `feat(instruments): registry, HighSpec definition, Mg b line table`.

---

### Task 5: Header map resolution and ingest

**Files:**
- Create: `src/mastpipe/io/__init__.py`, `src/mastpipe/io/headers.py`, `src/mastpipe/io/ingest.py`
- Test: `tests/unit/test_headers.py`, `tests/unit/test_ingest.py`

**Interfaces:**
- `HeaderMap.from_dict(d)`; `.resolve(header: Mapping, filename: str) -> tuple[dict, list[str]]` (`keys` priority, `type`, `lower`, `map`, `default`, `filename_regex` with group `value`; unresolved names listed); `.classify(fields) -> tuple[str, str]` (`exposure_type`, rule name `"header"` / `"rule:<i>"` / `"unclassified"`).
- `scan_folder(folder, instrument, header_map, pattern="*.fit*", manifest=None, unreadable=None) -> list[FrameRef]`; `load_frame(ref, detector) -> Frame` (float32, crop, then orientation: `transpose`, `flip_x`, `flip_y`, in that order; no variance); `read_manifest(path) -> pd.DataFrame`; `apply_manifest(refs, manifest) -> list[FrameRef]` (override any `FrameMeta` field by `source` [+ `extension`]); `file_checksum(path) -> str` (sha256, 16 hex).
- `arc_kind` is set from the header map when the keyword resolves; otherwise `None`, to be set by the runner from lit traces (Task 8).

- [ ] **Step 1: Failing tests** — `test_resolve_priority_and_transforms`, `test_filename_regex_fallback`, `test_classify_header_wins_then_rules`, `test_missing_exposure_time_is_unknown` (no `EXPOSURE`/`EXPTIME`, no manifest → `exposure_type == "unknown"`, rule name recorded), `test_scan_folder_builds_refs_with_checksums` (two tiny FITS in `tmp_path`), `test_unreadable_file_collected_not_raised`, `test_load_frame_crops_and_orients` (a 6×4 array with `transpose=True, flip_x=True` → expected array), `test_manifest_overrides_fields` (sets `exposure_type`, `target`, `arc_kind`), `test_mef_extensions_become_separate_refs`.
- [ ] **Step 2: Run** → fail. **Step 3: Implement.** **Step 4: Run** → pass. **Step 5: Commit** `feat(io): header map and ingest`.

---

### Task 6: Grating geometry and the synthetic-night simulator

**Files:**
- Create: `src/mastpipe/steps/__init__.py`, `src/mastpipe/steps/geometry.py`, `src/mastpipe/sim/__init__.py`, `src/mastpipe/sim/simulator.py`
- Test: `tests/unit/test_geometry.py`, `tests/unit/test_simulator.py`; fixtures in `tests/conftest.py`

**Interfaces:**
- `geometry.py`: `GratingGeometry(d_nm, alpha_rad, order, B_nom, lambda_littrow, dispersion_nm_px, lambda_center)`; `grating_geometry(geom: GeometryConfig, pixel_size_um) -> GratingGeometry` (port of the notebook's: `B_nom = f_cam/pixel`, `lambda_littrow = 2 d sin α / m`, `dispersion = d cos α / (m B_nom)`); `theta_of_lambda(lam, g, alpha=None)`; `grating_x(lam, x_c, B, alpha, g)`; `nominal_solution(g, n_pix) -> tuple[float, float]` = `(x_c, B)` placing `lambda_center` at `n_pix/2` with `B = B_nom`.
- `simulator.py`: `SimTruth` (`layout`, `trace_coeffs`, `x_c`, `B`, `lines_nm`, `line_logI`, `throughput`, `science_flux`, `sky_flux`, `sky_lines_nm`, `bias_level`, `dark_rate`, `read_noise_e`, `gain`, `fwhm_px`, `crop`, `frames: pd.DataFrame` with columns `file, exposure_type, arc_kind, mid_time, shift_px, stretch`, `extras`).
- `simulate_night(effective, layout, line_table, out_dir, seed=0, n_bias=5, n_dark=3, n_flat=3, n_full_arc=1, n_cal_arc=2, n_science=4, exposure_science=300.0, drift_px=(0.0, 0.6), drift_stretch=(0.0, 0.0), night_start="2026-10-01T18:00:00", cadence_min=12.0, cosmic_rays=20, fiber_preset="survey45", no_cal_arcs=False) -> SimTruth` — writes `sim_<type>_<nn>.fits` with header cards `IMAGETYP, EXPOSURE, DATE-OBS, TEMP, HBIN, VBIN, READMODE, GRATING, LAMP, OBJECT, FIBPRE` matching the HighSpec header map (no `ARCKIND` card: arc kind is detected from lit traces), uncropped `detector.shape`; `truth.npz` + `truth.json`; `load_truth(out_dir) -> SimTruth`.
- Simulation model: bias 500 ADU; read noise and gain from detector; dark 0.05 ADU/s; trace `i`: `y(x) = y_nominal_i + slope_i (x − xc) + c_i (x − xc)²`, `c_i ∈ ±2e-6`; Gaussian spatial profile with `fwhm_px = 3.0`; flat: smooth quadratic blaze × per-fiber throughput (0.8–1.2) × 20000 ADU peak on every trace; **full arc** (`arc_kind="full"`): every trace lit with the catalogue lines of `line_table` at `10**logI` relative intensities placed by the exact grating equation at `(x_c, B)` from `nominal_solution`; **calibration-fiber arc** (`"calfiber"`): only `role == "arc"` traces lit, the same lines shifted by `shift_px + stretch (x − n/2)` where shift and stretch are interpolated linearly in time between `drift_px[0]` at `night_start` and `drift_px[1]` at the end of the night (likewise `drift_stretch`); **science**: no arc light; science fibers: continuum 200 ADU/px with 6 Gaussian absorption dips (depth 0.3, FWHM 0.05 nm) at fixed wavelengths, drifted with the same time law; sky fibers: continuum 20 ADU/px plus 3 emission lines at `sky_lines_nm`, drifted likewise; science and sky are scaled by `throughput`; cosmic rays: `cosmic_rays` random single pixels + 5000 ADU in science frames; Poisson + Gaussian read noise; `uint16`, clipped. Frame order in time: bias, dark, flat, full arc, then alternating cal-arc / science groups at `cadence_min`, so every science frame has a cal arc before it except the first science frame when `n_cal_arc == 0` or `no_cal_arcs`.
- `tests/conftest.py` fixtures: `highspec_effective` (mgb518, main), `survey_layout`, `lab3_layout`, `mgb518_lines`, `sim_night(tmp_path_factory)` (session-scoped, seed 0, `n_science=4`), `sim_night_lab3`.

- [ ] **Step 1: Failing tests** — `test_nominal_dispersion_matches_notebook` (`dispersion_nm_px` ≈ 0.00866 ± 0.0002 for mgb518, and `lambda_littrow` ≈ 519.0 ± 0.5), `test_grating_x_inverts_theta` (round trip at three wavelengths to 1e-9 px), `test_simulate_night_writes_expected_frames` (counts per type, `arc_kind` column has one `full` and two `calfiber`, files exist, headers carry `EXPOSURE` and `DATE-OBS`), `test_truth_roundtrip`, `test_calfiber_arc_lights_only_arc_traces` (sum over science-trace rows of a cal arc ≈ background level), `test_science_frames_have_no_arc_lines` (no catalogue line above 5σ in a science fiber), `test_drift_interpolates_in_time` (`frames.shift_px` monotonic between the two cal arcs).
- [ ] **Step 2: Run** → fail. **Step 3: Implement.** **Step 4: Run** → pass (`pytest tests/unit/test_geometry.py tests/unit/test_simulator.py -v`). **Step 5: Commit** `feat(sim): grating geometry and synthetic night with full and calibration-fiber arcs`.

---

### Task 7: Clean — master calibrations, variance, cosmic rays (port of the notebook's calibration cell)

**Files:**
- Create: `src/mastpipe/steps/calib.py`, `src/mastpipe/steps/cosmics.py`
- Test: `tests/unit/test_calib.py`, `tests/unit/test_cosmics.py`

**Interfaces:**
- `stack_frames(frames: list[np.ndarray], method="median", clip_sigma=None) -> tuple[np.ndarray, np.ndarray]`; `build_master(frames: list[Frame], kind, params: CalibParams) -> MasterCalib` (validity: camera, binning, readout_mode, mean temperature, mean mid_time, exposure_time; provenance: input checksums); `build_calibration_set(bias, dark, background, detector) -> CalibrationSet` (`dark_rate = (dark − bias)/t_dark` when both exist else `0.0`); `calibrate_frame(frame, calib, subtract_background=True, saturation=None) -> Frame` (the notebook's arithmetic: `raw − bias − (background − bias)·t/t_bkg` or `raw − bias − dark_rate·t`; variance `max(signal_e, 0)/gain² + RN²/gain²` in ADU²; mask saturated + non-finite); `calibration_steps(calib, t_frame, subtract_background) -> list[str]`.
- `flag_cosmic_rays(data, variance, mask, sigma=8.0, size=5) -> np.ndarray` — pixels exceeding `sigma·sqrt(variance)` above the `size × size` median-filtered image are marked; returns the updated bool mask. Applied by the runner to science frames only.

- [ ] **Step 1: Failing tests** — `test_stack_median_and_variance`, `test_build_master_validity_and_provenance`, `test_calibrate_background_arithmetic_matches_notebook` (hand arrays), `test_calibrate_dark_path_when_no_background`, `test_variance_and_saturation_mask`, `test_cosmic_rays_flagged_on_simulated_science` (≥ 90 % of the simulator's cosmic-ray pixels masked, < 0.1 % false positives on the clean frame).
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): calibration masters, variance, cosmic-ray mask`.

---

### Task 8: Traces, layout matching, arc classification (port of `find_traces`)

**Files:**
- Create: `src/mastpipe/steps/trace.py`
- Test: `tests/unit/test_trace.py`

**Interfaces:**
- `robust_sigma(v) -> float`; `find_traces(image, params: TraceParams) -> TraceSet` (notebook algorithm: detection band `center_width` → baseline median filter → `find_peaks` → per-block centroids in `x_bin` columns → clipped polynomial fit of `degree`, two passes; sorted by `y_center`).
- `comb_shift(a: np.ndarray, b: np.ndarray, max_shift: float, step: float = 0.1) -> float` — shift maximising the overlap of two delta combs (positions), sub-pixel by parabola; used here for the layout and in Task 14 for lines.
- `match_layout(traceset, layout, tolerance) -> TraceSet` (global shift via `comb_shift`, assign `role, telescope, layout_index, aperture`, fill `missing`/`extra`).
- `trace_source(frames_by_type: dict[str, list[Frame]]) -> tuple[str, Frame | None]` — returns `("flat", master_flat)` when a flat exists, else `("arc_full", first full arc)`, else `("none", None)`.
- `classify_arc(frame: Frame, traceset: TraceSet, halfwidth: float, lit_sigma: float = 5.0) -> tuple[str, dict]` — per trace, median of the aperture-summed spectrum above `lit_sigma`·robust noise counts as lit; returns `"full"` when every matched trace is lit, `"calfiber"` when all `arc` traces are lit and no `science`/`sky` trace is, else `"unknown"`, plus the per-trace lit flags.
- `relative_throughput(flat: Frame, traceset, halfwidth) -> dict[int, float]` (integrated flat flux per trace, normalised to the median) — `{}` when no flat.

- [ ] **Step 1: Failing tests** — `test_find_traces_on_simulated_flat_recovers_45` (all 45 within 0.3 px of truth at three x positions), `test_find_traces_lab3_like_curvature` (quadratic rms < 0.02 px), `test_match_layout_with_global_shift_and_missing` (drop one trace from the image → `missing == (k,)`), `test_comb_shift_subpixel`, `test_classify_arc_full_vs_calfiber_vs_unknown` (simulator full arc, cal arc, science frame), `test_relative_throughput_matches_truth` (within 3 %), `test_trace_source_fallback_order`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): trace finding, layout matching, arc classification`.

---

### Task 9: Aperture extraction and continuum (port of `extract_lsf`, `fit_continuum`)

**Files:**
- Create: `src/mastpipe/steps/extract.py`
- Test: `tests/unit/test_extract.py`

**Interfaces:**
- `aperture_weights(y_center, halfwidth, aperture) -> tuple[np.ndarray, np.ndarray]`; `extract_aperture(frame, traceset, params: ExtractParams) -> list[Spectrum]` (flux = weighted sum or mean; error = sqrt of weighted variance; mask where coverage < `min_coverage` or any masked pixel in the aperture; per-trace `aperture` override honoured; `role/telescope/target` copied, `target` from `frame.meta.target` for science traces); `fit_continuum(spectrum, degree, clip_above, clip_below, n_iter=10) -> Spectrum` (new Spectrum with `continuum` set; raw `flux` untouched — the continuum-subtracted flux is `flux − continuum`, computed by callers); `continuum_subtract_all(spectra, params) -> list[Spectrum]`.

- [ ] **Step 1: Failing tests** — `test_aperture_weights_round_subpixel_floor`, `test_extract_recovers_simulated_flux` (science fibers: median `flux/truth` within 2 %), `test_masked_pixel_masks_column`, `test_continuum_fit_ignores_lines` (synthetic continuum + emission lines: continuum within 1 %), `test_coverage_below_min_is_masked`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): aperture extraction and continuum`.

---

### Task 10: Line finding (port of `find_lines`)

**Files:**
- Create: `src/mastpipe/steps/lines.py`
- Test: `tests/unit/test_lines.py`

**Interfaces:**
- `find_lines(image: np.ndarray, variance: np.ndarray | None, traceset: TraceSet, params: LineParams) -> tuple[list[pd.DataFrame], list[np.ndarray]]` — one table and one detilted strip per trace, notebook columns: `id, group, n_in_group, x_peak, x_centroid, x_fit, x_center, y_center, y_offset, height, flux, baseline, snr, fwhm_x, fwhm_x_fit, sigma_x_fit, amp_x_fit, fwhm_y, sigma_y_fit, x_window, neighbour, isolated, near_edge, complete, on_trace, fit_ok, saturated, use_x, use_y`; `saturated` when the raw peak exceeds `params.saturation`; `use_x` requires `fit_ok or isolated`, not `near_edge`, `on_trace`, not `saturated`; private helpers `_detilt`, `_remove_row_continuum`, `_fit_group`, `_deblend_group` ported as in the notebook.
- `lines_summary(tables) -> pd.DataFrame` (per trace: n lines, n use_x, n blended, median fwhm_x).

- [ ] **Step 1: Failing tests** — `test_find_lines_recovers_simulated_arc_positions` (full-arc frame, trace 0: ≥ 60 % of catalogue lines with `logI > 1.5` found within 0.15 px of truth), `test_deblend_two_close_lines` (synthetic pair 4 px apart → two components, `n_in_group == 2`), `test_cosmic_ray_off_trace_rejected` (`max_y_offset`), `test_saturated_line_not_usable` (peak above `saturation` → `use_x False`), `test_edge_lines_flagged`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): line finding with deblending`.

---

### Task 11: Wavelength model, atlas and usable lines

**Files:**
- Create: `src/mastpipe/steps/wavemodel.py`, `src/mastpipe/steps/atlas.py`; finalize `scripts/make_line_table.py`
- Test: `tests/unit/test_wavemodel.py`, `tests/unit/test_atlas.py`

**Interfaces:**
- `WavelengthModel` (Protocol): `kind: str`; `n_pix: int`; `params: np.ndarray`; `wavelength(x) -> np.ndarray`; `pixel(lam) -> np.ndarray`; `dispersion(x) -> np.ndarray` (nm/px); `with_params(p) -> WavelengthModel`; `to_dict() / from_dict(d, geometry)`.
- `GratingModel(x_c, B, alpha, poly: np.ndarray | None, n_pix, geom: GratingGeometry)` implementing the protocol with the notebook's `wave_from_solution` (`λ = d (sin α + sin(α + atan((x − P(t) − x_c)/B)))/m`, `t = (x − n/2)/(n/2)`), `pixel` by `grating_x` plus the distortion term solved iteratively (two Newton steps), `dispersion` as the central difference; `params` = `[x_c, B, alpha, *poly]`.
- `model_from_config(kind, geom, n_pix) -> WavelengthModel` — `"grating"` → `GratingModel` at `nominal_solution`; `"dispersion"` → `NotImplementedError("Plan 2")`.
- `model_from_solution(solution: WavelengthSolution, geom: GratingGeometry) -> WavelengthModel` — rebuilds the model from `solution.kind`, `solution.params`, `solution.n_pix` (the inverse of `WavelengthSolution(..., params=model.params, ...)`); used by Tasks 12, 14, 17 and 18 wherever a stored solution is evaluated.
- `fit_model(model, x, lam, weights=None, fit_alpha=False, poly_deg=0, alpha_bound_deg=None, B_bound_frac=None) -> tuple[WavelengthModel, np.ndarray]` — the notebook's `fit_solution` (least squares on `x_c, B[, alpha][, P]` with `poly_powers`), returns fitted model and residuals in px.
- `anchor_model(pairs: list[tuple[float, float]], geom, n_pix) -> GratingModel` — the notebook's `fit_anchor_solution`.
- `atlas.py`: `load_murphy_atlas(folder: Path, band, pad=0.5) -> tuple[pd.DataFrame, pd.DataFrame]` (spectrum `wavelength, flux` in air nm from the log-linear header; catalogue `wavenumber, wavelength, logI, element, ion, source, species`); `remove_pedestal(flux, window=201, q=5)`; `convolve_to_resolving_power(wave, flux, R) -> tuple[np.ndarray, float]`; `flag_usable_lines(wave, flux_conv, cat, fwhm_nm, max_shift=0.02, min_snr=20) -> pd.DataFrame` (adds `shift_fwhm, shift_pm, snr, usable`); `synthetic_atlas(cat, band, R, step_nm=0.0005) -> tuple[np.ndarray, np.ndarray]` (Gaussian lines at `10**logI`, for tests only); `catalogue_arrays(table, usable_only) -> tuple[np.ndarray, np.ndarray]` (`cat_wl` sorted, `cat_sp`).

- [ ] **Step 1: Failing tests** — `test_grating_model_roundtrip_pixel_wavelength` (1e-6 px), `test_grating_model_matches_simulator_truth` (model at truth `(x_c, B)` reproduces `line_pixels` within 1e-6), `test_fit_model_recovers_params_from_lines` (perturbed start, 25 lines with 0.02 px noise → `x_c` within 0.05 px, rms < 0.03 px), `test_anchor_model_two_pairs`, `test_poly_powers_excludes_low_orders`, `test_dispersion_model_not_implemented`, `test_flag_usable_lines_on_synthetic_atlas` (a blended pair becomes unusable, an isolated line usable), `test_convolve_fwhm`, `test_load_murphy_atlas_real` (skipped unless `MASTPIPE_ATLAS_DIR`; 73 lines in the Mg b band, 31 ± 3 usable at R 23081 — the notebook's numbers).
- [ ] **Step 2–4:** run → fail; implement; run → pass. **Step 5: Regenerate the line table** with the script against the real atlas when `MASTPIPE_ATLAS_DIR` is set and confirm `git diff --stat src/mastpipe/instruments/highspec/lines/mgb518.csv` is empty (otherwise commit the regenerated table and note it in CHANGES). **Step 6: Commit** `feat(steps): grating wavelength model, atlas tools, usable-line table`.

---

### Task 12: Full-arc solve and reference solutions (port of `solve_trace`)

**Files:**
- Create: `src/mastpipe/steps/wavecal.py`
- Test: `tests/unit/test_wavecal.py`

**Interfaces:**
- `identify_lines(x_obs, model, cat_wl, cat_sp, tol_px, unique_tol_px) -> pd.DataFrame` (`x, wavelength, species, dx_px`).
- `shift_to_reference(flux, flux_ref, max_shift=200) -> float` (the notebook's cross-correlation with parabolic peak).
- `solve_fiber(x_obs: np.ndarray, model0: WavelengthModel, cat_wl, cat_sp, params: WavecalParams) -> tuple[WavelengthModel | None, pd.DataFrame]` — identify → fit → clip, `n_passes`, convergence `converge_px`, `min_lines`; returns the fitted model and the identified lines.
- `solve_full_arc(line_tables: list[pd.DataFrame], flux_list: list[np.ndarray], traceset: TraceSet, start: ReferenceSolution | WavelengthModel, cat_wl, cat_sp, params, fibers: list[int] | None = None) -> dict[int, WavelengthSolution]` — per fiber: start from the previous reference's model for that fiber when `start` is a `ReferenceSolution` (missing fibers seed from the nearest solved fiber plus `shift_to_reference`), else from the single model slid by `shift_to_reference` against fiber 0; fibers with no `use_x` lines or `< min_lines` matches are absent from the result; `WavelengthSolution.anchor_px/anchor_nm` = range of identified lines; `trace_y_center` from the trace set.
- `build_reference(instrument, configuration, fiber_preset, solutions, created, source_checksums) -> ReferenceSolution` (identifier via `make_identifier`, hash = sha256 over the solutions' params and lines, 12 hex).
- `fiber_agreement(solutions) -> pd.DataFrame` — for every catalogue line identified in ≥ 2 fibers, the spread of implied wavelengths (the notebook's external accuracy check): columns `wavelength, n_fibers, spread_px`.
- `wavelength_map(solution) -> np.ndarray` (nm per pixel, `n_pix`).

- [ ] **Step 1: Failing tests** — `test_solve_fiber_from_anchor_on_simulated_full_arc` (fiber 0 of the simulated full arc with the two anchors transformed to the simulated geometry → rms < 0.05 px, ≥ 20 lines), `test_solve_full_arc_all_45_from_single_model` (≥ 43 solved; for every solution `|λ_fit − λ_true| < 0.3·dispersion` at all pixels), `test_solve_full_arc_from_previous_reference_is_faster_and_equal` (second solve seeded by the first: same line sets, equal params within 1e-6), `test_missing_fiber_absent_not_crash` (drop one fiber's lines → key absent), `test_fiber_agreement_small_on_truth`, `test_reference_hash_changes_with_params`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): full-arc wavelength solve and reference solutions`.

---

### Task 13: `identify` and `solve` commands with anchor panels

**Files:**
- Create: `src/mastpipe/io/qa.py`, `src/mastpipe/pipeline/__init__.py`, `src/mastpipe/pipeline/identify.py`; modify `src/mastpipe/cli.py`
- Test: `tests/unit/test_identify.py`

**Interfaces:**
- `qa.plot_anchor_panels(atlas_wave, atlas_flux_conv, usable: pd.DataFrame, obs_flux, obs_lines: pd.DataFrame, band, n_labels=8, min_sep=0.04, flip_x=False, out: Path) -> tuple[pd.DataFrame, pd.DataFrame]` — the notebook's two-panel figure written as PNG, returning the numbered template and measured peaks (`id, wavelength` / `id, x_center`).
- `identify.run_identify(full_arc: Path, effective, layout, out_dir, anchors: list[tuple[float, float]] | None, anchor_trace=0, atlas_dir: Path | None, write_reference=True) -> tuple[ReferenceSolution | None, Path]` — loads and calibrates the frame (bias/background from `out_dir/calib` if present, else none), finds traces (flat if present in the night folder, else the arc itself), extracts, continuum-subtracts, finds lines; **without anchors**: writes `anchor_panels.png` and `anchor_candidates.csv` and returns `(None, png_path)`; **with anchors** (from the argument, else the instrument's stored anchors for the configuration/preset): builds the anchor model, `solve_full_arc` from it, writes the reference product (Task 17) to `out_dir/calib/` and `fiber_agreement.csv`; the atlas spectrum is used for the panels only when `atlas_dir` is given, else `synthetic_atlas` from the line table.
- `identify.run_solve(full_arc: Path, effective, layout, out_dir, previous: ReferenceSolution) -> ReferenceSolution` — same chain seeded by `previous` (`solve_full_arc` with a `ReferenceSolution` start); the caller finds `previous` (Task 18 wires the library lookup and the `solve` CLI command, which raises "run identify first" when the library has none).
- CLI (this task): `mastpipe identify <full_arc.fits> --out <dir> --instrument highspec --configuration mgb518 [--layout lab3] [--anchors 515.8605:267.22,525.8360:1423.22 | --anchors-file f.yaml | --stored-anchors] [--trace 0] [--atlas-dir DIR]`.

- [ ] **Step 1: Failing tests** — `test_identify_without_anchors_writes_panels_and_candidates`, `test_identify_with_anchors_writes_reference` (simulated full arc, anchors from truth → reference file exists, ≥ 43 fibers, `valid_from == frame mid_time`), `test_run_solve_seeded_by_previous_reference` (same line sets, params equal within 1e-6), `test_cli_identify_parses_anchor_pairs`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(pipeline): identify and solve commands`.

---

### Task 14: Drift — fit on calibration fibers, transfer across the detector, time pairing

**Files:**
- Create: `src/mastpipe/steps/drift.py`, `src/mastpipe/pipeline/arcs.py`
- Test: `tests/unit/test_drift.py`, `tests/unit/test_arcs.py`

**Interfaces:**
- `drift.fit_fiber_drift(x_obs: np.ndarray, reference: WavelengthSolution, params: DriftParams) -> dict` — predicted pixels are `reference.lines.x`; global shift guess by `comb_shift(x_obs, reference.lines.x, params.max_shift_px)`; match within `tol_px`; clipped linear fit `x_obs − x_ref = shift + stretch·(x_ref − n_pix/2)`; `stretch` set to 0 when `|stretch| < stretch_significance·σ_stretch`; returns `shift_px, stretch, n_lines, rms_px, usable` (`usable` = `n_lines ≥ min_lines and rms_px ≤ max_scatter_px`).
- `drift.fit_drift(line_tables_by_fiber: dict[int, pd.DataFrame], reference: ReferenceSolution, traceset, arc_meta: FrameMeta, params) -> DriftModel` — per arc fiber `fit_fiber_drift` on `use_x` lines; polynomial of `min(params.degree, n_usable − 1)` in `y_center` for shift and for stretch; `scatter_px` = rms of usable fibers about the shift model; `usable` = ≥ 1 usable fiber and `scatter_px ≤ max_scatter_px`.
- `drift.apply_drift(reference: ReferenceSolution, drifts: list[tuple[DriftModel, float]], traceset, geom: GratingGeometry) -> dict[int, np.ndarray]` — for every fiber with a reference solution: weighted shift `s = Σ w_k·shift_k(y)`, stretch `e` likewise; `x_ref = (x − s + e·n/2)/(1 + e)`; wavelength = `model_from_solution(reference.solutions[f], geom).wavelength(x_ref)` (Task 11); the geometry is passed in as `geom: GratingGeometry`.
- `arcs.pair_arcs(exposure: FrameMeta, arcs: list[tuple[FrameMeta, DriftModel]], reference: ReferenceSolution, params) -> list[tuple[DriftModel, float]]` — candidates: same `configuration`, `usable`, `mid_time ≥ reference.valid_from`, `|Δt| ≤ max_age_hours`; `nearest_before` → the latest before (else nearest after); `nearest` → smallest `|Δt|`; `interpolate` → before and after weighted `(1 − f, f)` by time when both exist, else nearest; `[]` when none.
- `arcs.status_for(pairs) -> str` → `"transferred"` when pairs non-empty else `"reference"`.

- [ ] **Step 1: Failing tests** — `test_fit_fiber_drift_recovers_shift_and_stretch` (synthetic lines shifted 0.6 px, stretched 2e-4 → within 0.02 px / 2e-5), `test_insignificant_stretch_is_zero`, `test_fit_drift_on_simulated_cal_arc` (truth shift at that time within 0.05 px, `usable`), `test_apply_drift_matches_science_truth` (drifted absorption dips land at their wavelengths within 0.3·dispersion), `test_pairing_policies` (three arcs at t−2h, t−10min, t+15min: `nearest_before` → t−10min; `nearest` → t−10min; `interpolate` → weights ≈ (0.6, 0.4)), `test_pairing_ignores_other_configuration_and_stale_arcs`, `test_no_candidates_gives_reference_status`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): drift fit, transfer and time pairing`.

---

### Task 15: Sky subtraction

**Files:**
- Create: `src/mastpipe/steps/sky.py`
- Test: `tests/unit/test_sky.py`

**Interfaces:**
- `subtract_sky(science: Spectrum, sky: Spectrum, throughput_ratio: float, params: SkyParams) -> Spectrum` — requires both `wavelength` arrays; sky flux scaled by `throughput_ratio` (science/sky), resampled linearly onto the science wavelength grid, subtracted; returns a copy of `science` with `sky` set to the resampled scaled sky and `versions["skysub"] = SpectrumVersion(flux − sky, sqrt(err² + err_sky²), mask | mask_sky_resampled)`; raises `ValueError` when either wavelength is missing.
- `subtract_sky_all(spectra: list[Spectrum], traceset: TraceSet, throughput: dict[int, float], params) -> tuple[list[Spectrum], list[dict]]` — pairs each science fiber with its telescope's sky fiber via `traceset.by_telescope`; a science fiber whose sky fiber is missing, absent from `spectra`, fully masked or without wavelength is returned unchanged and listed in the second element as `{"fiber", "telescope", "reason"}`; throughput ratio 1.0 when `throughput` lacks either fiber.

- [ ] **Step 1: Failing tests** — `test_subtract_sky_removes_simulated_sky_lines` (residual at `sky_lines_nm` < 3σ), `test_skysub_version_and_raw_kept` (raw `flux` unchanged, `versions["skysub"]` present), `test_requires_wavelength`, `test_missing_sky_fiber_leaves_raw_only` (reason `"no sky fiber"`), `test_masked_sky_fiber_leaves_raw_only`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): sky subtraction from the telescope's sky fiber`.

---

### Task 16: Within-exposure combination, flux and telluric pass-through

**Files:**
- Create: `src/mastpipe/steps/combine.py`, `src/mastpipe/steps/fluxcal.py`, `src/mastpipe/steps/telluric.py`
- Test: `tests/unit/test_combine.py`, `tests/unit/test_passthrough.py`

**Interfaces:**
- `combine.common_grid(spectra, log=True) -> np.ndarray` — from the union wavelength range, step = the minimum local dispersion (log-uniform when `log`).
- `combine.combine_spectra(spectra: list[Spectrum], target: str, params: CombineParams, version: str | None = None, barycentric: dict[str, float] | None = None, qc: dict[int, dict] | None = None) -> CombinedSpectrum` — uses `.latest(upto=version)` of each spectrum; include/exclude telescope lists; optional rejection (`reject`): a spectrum is excluded with a reason when its median S/N < `reject_snr_min`, or its median flux ratio to the group median < `reject_flux_ratio`, or (when `qc` gives `drift_px`) `|drift| > reject_drift_px`; inverse-variance weighted mean on the common grid with `clip_sigma` outlier rejection when ≥ 3 contributors; **every version present in all contributors is combined separately** and stored under the same key, raw under `"raw"`; `contributors`/`excluded` tables filled; `grouping = "exposure"`.
- `combine.group_by_target(spectra) -> dict[str, list[Spectrum]]` (science fibers only, by `target`; `None` target → key `"unknown"`).
- `fluxcal.apply_response(spectrum, response: np.ndarray | None = None) -> Spectrum` — returns the input unchanged when `response is None`; otherwise adds `versions["fluxcal"]` = latest ÷ response. `fluxcal.flux_status(spectra) -> str`.
- `telluric.correct_telluric(spectrum, transmission: np.ndarray | None = None) -> Spectrum` — same pattern with `versions["telluric"]`; `telluric.telluric_status(spectra) -> str`.

- [ ] **Step 1: Failing tests** — `test_common_grid_log_step`, `test_combine_weighted_mean_reduces_noise` (4 noisy copies → noise ÷ ≈ 2), `test_exclude_list_and_reason`, `test_rejection_low_snr_listed_not_dropped_silently`, `test_versions_combined_separately` (raw and skysub both present), `test_group_by_target_science_only`, `test_passthrough_without_calibration_leaves_spectrum_and_status` (`flux_status == "uncalibrated"`, `telluric_status == "uncorrected"`, no new version), `test_passthrough_with_response_adds_version`.
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(steps): within-exposure combination; flux and telluric pass-through stages`.

---

### Task 17: Provenance and FITS products

**Files:**
- Create: `src/mastpipe/pipeline/provenance.py`, `src/mastpipe/io/products.py`
- Test: `tests/unit/test_products.py`

**Interfaces:**
- `Provenance(pipeline_version, config_hash, inputs: list[str], calibrations: list[str], wavelength_status="none", flux_status="uncalibrated", telluric_status="uncorrected", reference_id=None, reference_hash=None, drift_arcs: list[str] = [], steps: list[str] = [])` with `to_header_cards() -> list[tuple[str,str,str]]` (`MPVER, MPCFGH, MPWSTAT, MPFSTAT, MPTSTAT, MPWREF, MPWREFH, MPNIN, MPIN0001…, MPNCAL, MPCAL001…, MPNDARC, MPDARC01…, MPSTEPnn`) and `from_header(header)`.
- `product_name(meta, kind, ext=".fits") -> str` = `"{instrument}_{configuration or 'noconf'}_{camera}_{YYYYmmddTHHMMSS or 'notime'}_{kind}{ext}"`; `reference_name(ref) -> str` = `f"{ref.identifier}_wref.fits"`; `combined_name(meta, target) -> str`.
- `write_master / read_master`; `write_traceset / read_traceset` (tables `TRACES`, `SAMPLES`, cards `MPSHIFT, MPMISS, MPEXTRA, MPNY, MPNX`); `write_spectra(path, es, prov) / read_spectra(path)` — table `SPECTRA`, one row per fiber, vector columns `PIXEL, FLUX, ERROR, MASK, WAVELENGTH, SKY, CONTINUUM` and per version `FLUX_<V>, ERROR_<V>, MASK_<V>` for `V in SKYSUB, FLUXCAL, TELLURIC` (NaN/True-filled with presence cards `MPHAS<V>` per row-set), scalar columns `FIBER, ROLE, TELESCOPE, TARGET`; `FrameMeta` as numbered cards `MPMKnn/MPMVnn`; optional table `DRIFT` (`DriftModel.per_fiber` plus cards `MPDSCAT`, `MPDUSE`) when the exposure used drift arcs; `write_reference(path, ref, prov) / read_reference(path, geometry_lookup)` — tables `SOLUTIONS` (`FIBER, KIND, NPARAMS, PARAMS(12), NPIX, RMSPX, NLINES, ANCHPX0, ANCHPX1, ANCHNM0, ANCHNM1, YCENTER`) and `LINES` (`FIBER, X, WAVELENGTH, SPECIES, DXPX`), cards `MPREFID, MPREFH, MPCREAT, MPVALID, MPINST, MPCONF, MPPRESET`; `write_combined(path, cs, meta, prov) / read_combined(path)` — table `COMBINED` with `WAVELENGTH` and per version `FLUX_<V>, ERROR_<V>, MASK_<V>` (`RAW` included), tables `CONTRIB`, `EXCLUDED`, cards `MPTARGET, MPGROUP`, barycentric as `MPBCnn` pairs; `write_qc(path, report)` JSON.

- [ ] **Step 1: Failing tests** — `test_provenance_cards_roundtrip` (including drift arcs and reference id), `test_product_name_format`, `test_master_roundtrip`, `test_traceset_roundtrip_degrees`, `test_spectra_roundtrip_with_versions_and_drift` (raw and `skysub` arrays equal, statuses equal, `DRIFT` table present), `test_spectra_roundtrip_without_versions` (presence cards false, `versions == {}`), `test_reference_roundtrip` (params, lines, `valid_from`, hash equal; `GratingModel` reconstructed and evaluates equally), `test_combined_roundtrip`, `test_no_raw_keywords_in_products` (no card named `EXPOSURE`, `IMAGETYP`, `OBJECT` in any product header).
- [ ] **Step 2–5:** run → fail; implement; run → pass; commit `feat(io): FITS products with provenance, kept versions, reference and combined spectra`.

---

### Task 18: Calibration library, quick-look runner, CLI

**Files:**
- Create: `src/mastpipe/pipeline/callib.py`, `src/mastpipe/pipeline/runner.py`; modify `src/mastpipe/cli.py`
- Test: `tests/unit/test_callib.py`, `tests/e2e/test_quicklook.py`

**Interfaces:**
- `CalibLibrary(root)`: `add(master, prov, meta_for_name) -> Path`; `find(kind, camera, binning, readout_mode, configuration=None, temperature=None, temperature_tolerance=2.0, near_time=None)`; `add_reference(ref, prov) -> Path`; `find_reference(instrument, configuration, fiber_preset, at_time: str | None = None, identifier: str | None = None) -> ReferenceSolution | None` — with `identifier` that one; else the reference with the latest `valid_from ≤ at_time` (latest overall when `at_time is None`); `list(kind=None)`; index in `index.json`.
- `RunOptions(instrument, configuration=None, camera=None, layout=None, overrides: list[dict] = (), on_error: Literal["continue","raise"] = "continue", force=False, manifest=None, reference_id=None, time_policy=None, library: Path | None = None)`.
- `run_quicklook(night_dir, out_dir, options) -> RunSummary(products, skipped, failed: list[tuple[str,str]], qc, reduced_again: list[str])`.
- Recipe per `(camera, configuration)` group:
  1. Scan + manifest. Frames of type `unknown` are recorded in `skipped` with their classification rule.
  2. Masters (bias, dark, background, flat) → library at `out_dir/calib` (or `options.library`).
  3. Traces from `trace_source`; `match_layout`; throughput from the flat or `{}`; write trace product + QC.
  4. Arcs: calibrate, `classify_arc` when `arc_kind is None`; **full arcs** → if a reference exists for the configuration/preset: `run_solve`-equivalent (seeded solve) and `add_reference`; else logged "no reference; run identify" and skipped; **cal arcs** → extract, continuum, `find_lines`, `fit_drift` against the reference valid at the arc's time → drift models kept in memory and written as `…_drift.fits` QC JSON.
  5. Science frames: calibrate + `flag_cosmic_rays`, extract, reference by `find_reference(at_time=mid_time)` (or `options.reference_id`), `pair_arcs` with the policy (option overrides the mode default `nearest_before`), `apply_drift` or reference wavelengths, status via `status_for`; `nominal` wavelengths from `model_from_config` when no reference; sky subtraction; `apply_response` / `correct_telluric` with `None`; `write_spectra` + QC (per fiber median S/N, masked fraction, drift shift/stretch used, sky reasons); `combine_spectra` per target with `version=None` → `write_combined`.
  6. Idempotence key per science product: input checksum + config hash + reference id + sorted drift-arc checksums; existing product with the same key and not `force` → `skipped`. A product whose stored key has status `reference` while the current pairing yields drift arcs is reduced again and listed in `reduced_again`.
  7. Error policy and logging as in the September plan (per-frame `continue` catches; group-level failures recorded; `out_dir/logs/quicklook_<stamp>.log`).
- CLI: `mastpipe quicklook <night_dir> --out <dir> --instrument highspec [--configuration mgb518] [--camera main] [--layout survey45] [--manifest f.csv] [--on-error continue|raise] [--force] [--reference ID] [--time-policy …] [--library DIR] [--set key.path=value …]`; `mastpipe solve <full_arc.fits> --out <dir> --instrument … --configuration … [--library DIR]` (finds the previous reference in the library, calls `run_solve`, adds the result; exits 2 with "run identify first" when none); `mastpipe show <product.fits>` prints kind, statuses, provenance, and per-fiber summary.

- [ ] **Step 1: Failing tests** — `test_callib_find_reference_by_time_and_identifier` (two references, `at_time` between → the earlier; after → the later; identifier override), `test_quicklook_end_to_end_on_simulated_night` (session simulated night with a reference made by `run_identify` from truth anchors: every science product has `wavelength_status == "transferred"`, science-fiber wavelengths within `0.3·dispersion` of truth, `skysub` present, combined product per target exists, QC JSON lists drift shift), `test_quicklook_without_drift_arcs_gives_reference_status` (`no_cal_arcs=True` night), `test_quicklook_without_reference_gives_nominal`, `test_quicklook_idempotent_then_reduced_again` (run; add a later cal arc file; run again → `reduced_again` non-empty for the `reference`-status products only), `test_unknown_frames_are_skipped_with_reason`, `test_quicklook_continues_after_bad_frame`, `test_cli_quicklook_and_show`.
- [ ] **Step 2–5:** run → fail; implement; run → pass (`pytest tests -v`); commit `feat(pipeline): calibration library with time-valid references, quick-look runner, CLI`.

---

### Task 19: Regression on the notebook's laboratory data, documentation, full verification

**Files:**
- Create: `tests/regression/test_lab_highspec.py`; modify `README.md`, `CHANGES.md`, `DECISIONS.md`
- Expected layout under `MASTPIPE_TEST_DATA`: `highspec/lab3/calib/{image_bias.fits, background.fits, flat_518.fits}` and `highspec/lab3/scan/FOCUS_05_13.30004_mm.fits` (plus the other scan frames, unused here); the Murphy atlas under `MASTPIPE_ATLAS_DIR`.

- [ ] **Step 1: Regression tests** (all `pytest.mark.skipif` without the variables) — `test_lab_traces_match_notebook` (3 traces at `y_center` 17.40 / 188.96 / 360.72 ± 0.3 px, quadratic rms < 0.01 px), `test_lab_reference_from_stored_anchors` (`run_identify` with `--stored-anchors` on `FOCUS_05`: three solutions, rms ≤ 0.05 px each, coverage 513.50–527.32 nm ± 0.03, ≥ 19 lines each, `fiber_agreement.spread_px` median < 0.03), `test_lab_drift_against_itself_is_zero` (`fit_fiber_drift` of the same frame's lines against its own reference → shift < 0.01 px), `test_lab_quicklook_runs` (`run_quicklook` on a folder with the calib frames and `FOCUS_05` as a manual-classified full arc via manifest → reference written, no failures).
- [ ] **Step 2: Run with data** — `MASTPIPE_TEST_DATA=… MASTPIPE_ATLAS_DIR=… pytest tests/regression -v` → pass; without the variables → skipped.
- [ ] **Step 3: Docs** — README: commands (`identify`, `solve`, `quicklook`, `show`), the three-cadence wavelength scheme in one paragraph, product layout and columns, environment variables; CHANGES 0.0.1 entries per task; DECISIONS: anchors reproduced instead of imported; `comb_shift` shared by layout and drift; versions as columns; drift key in idempotence.
- [ ] **Step 4: Full verification** — `ruff check src tests scripts && ruff format --check src tests scripts && python scripts/privacy_scan.py && pytest -q` → all green; push nothing; confirm `git status` clean and `git log --oneline | wc -l ≥ 19`.
- [ ] **Step 5: Commit** `docs: README, changelog, decisions; regression suite on laboratory data`.

---

## Self-review notes (checked against the spec and revision 1)

- **Coverage.** Revision §2.1 cadences → Tasks 6 (simulator), 8 (`classify_arc`), 12–13 (full arc → reference), 14 (drift, pairing). §2.2 time policy, statuses, re-reduce → Tasks 14, 18. §2.3 validity, never edited in place, anchors reproduced → Tasks 12, 17, 18, 19. §2.4 model interface with the `dispersion` slot → Task 11. §2.5 `identify` → Task 13. §3 stages 1–6 → Tasks 7, 9, 14, 15, 16 (pass-through), combination Task 16; quick-look scope → Task 18. §4 DeepSpec → Plan 2 (orientation and per-camera fields already in Tasks 3, 5). §5 kept versions and status cards → Tasks 2, 17. §8 calls → DECISIONS in Task 1, trace fallback Task 8, pairing bounds Task 14, time-keyed references Task 18, re-reduce Task 18, versions as columns Task 17.
- **Types.** `SpectrumVersion`, `versions` keys, `WavelengthSolution.lines` columns, `DriftModel.per_fiber` columns, `ReferenceSolution.identifier` format and `comb_shift` are defined once (Tasks 2, 8) and used by name afterwards.
- **Review Focus** items 1–5 each have a named test in Tasks 18, 15, 14, 10, 5/18.
- **Proportion.** Signatures, test names and assertions only; the only algorithm text is where the notebook port leaves a choice (drift formula, pairing weights, grid step).
