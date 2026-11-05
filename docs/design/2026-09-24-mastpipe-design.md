# mastpipe — design specification

*Data-reduction pipeline for the MAST spectrographs (HighSpec and DeepSpec).*
Written 2026-09-24 (Claude Fable 5.1, from a design session with the project lead).
This document is written to be copied verbatim into the `mastpipe` repository:
it contains no local paths, hostnames or credentials.

## 1. Purpose and decision

`mastpipe` turns raw frames from the MAST fiber-fed spectrographs into
wavelength-calibrated, sky-subtracted spectra per fiber, combined spectra per
target, and quality-assurance products, in two modes: an unattended quick-look
at the observatory and a careful manual reduction by a researcher.

### 1.1 Options considered

| Option | Verdict | Reason |
|---|---|---|
| Build on **soxspipe** (SOXS / X-shooter pipeline) | Rejected | Echelle-order data model with no fiber concept; wavelength solution requires ESO single- and multi-pinhole frames plus predicted-line tables; 34 hard-coded instrument branches; frame classification bound to ESO `DPR` keywords; effectively single-maintainer, GPL-3, no contribution policy. |
| Build on **PypeIt** | Deferred | Strong community project (BSD-3, ~5 core developers). Its new `Fiber` pipeline (per-fiber optimal extraction, joint sky fit over sky fibers) is the right shape but unreleased at design time, still carries instrument-specific logic, and has no concept of arc fibers exposed alongside science fibers. A framework tracking a development branch is a fragile base for an unattended service. |
| **Independent package** seeded from the HighSpec laboratory reduction notebook (`HighSpec_data_reduction_V3`) | **Chosen** | The notebook already contains the physics that matters most: per-fiber grating-equation wavelength solutions, line-width and focus analysis. What it lacks (packaging, science extraction, sky, flats, on-disk products) is generic. |

The chosen design keeps a **future PypeIt integration open**: reduction steps
are pure functions on arrays, so a PypeIt spectrograph class could wrap them
once the `Fiber` pipeline is mature. PypeIt array-level routines may be used
today through one optional adapter module.

### 1.2 Scope

In scope: both instruments; bias/dark/background/flat calibration; tracing;
aperture and optimal extraction; per-fiber wavelength solutions from full arcs;
transfer of the solution to science fibers via the simultaneous arc fibers; sky
subtraction from sky fibers; combination across telescopes and exposures;
FITS products with provenance; quality-assurance figures and pages; a focus-scan
diagnostic recipe; a synthetic-frame simulator; command-line and Python
interfaces.

Out of scope for the first release: absolute flux calibration and telluric
correction; pixel-to-pixel flat-fielding (designed in, switched off until a
properly illuminated flat exists); integration with the observatory's process
supervisor and dashboard (the products are designed so that this is easy
later); radial-velocity measurement (belongs to downstream analysis tools).

## 2. Constraints

- **Instruments.** HighSpec: one detector, straight nearly parallel traces,
  narrow-band gratings each covering ~14 nm at R ≈ 20 000, Th-Ar arcs. Three
  gratings exist today (Mg b 518 nm, Hα, Ca II H&K) and **more will be added**,
  so a grating must be a configuration entry (band, geometry, atlas, line
  table), never code. DeepSpec: four cameras (bands U, G, R, I), R ≈ 500.
  Survey layout planned: 45 traces per detector, 5 arc, 20 science, 20 sky (one
  science and one sky fiber per telescope). All telescopes expose
  simultaneously with the same exposure time, not necessarily on the same
  target. Arcs are planned daily and at every instrument switch.
- **Headers will change.** The control software's FITS keywords are provisional.
  Nothing downstream of ingest may depend on a raw keyword name.
- **Platforms.** Must run on Windows, Linux and macOS. Pure Python, `pathlib`
  everywhere, no POSIX-only calls, no shell-outs.
- **Standalone.** A self-contained Python package with no dependency on other observatory software stacks.
- **Privacy.** The repository contains no local paths, hostnames, IP addresses,
  usernames, credentials or unpublished data. Real test data live outside the
  repository, located through an environment variable.
- **Python** 3.11 or newer. Hard dependencies: numpy, scipy, astropy, pandas,
  matplotlib, pydantic (v2), pyyaml. Optional extras: `pypeit`, `dev`.

## 3. Architecture

Three layers, kept strictly separate.

1. **Instruments** — declarative configuration only: detector parameters,
   header map, fiber-layout tables, one block per grating or band with the
   dispersion geometry and atlas. No code beyond validation models. Adding a
   HighSpec grating means adding one YAML file under `gratings/` (plus its
   atlas band and, once identified, its stored line table); the registry picks
   it up without any code change.
2. **Steps** — pure functions: arrays and small dataclasses in, dataclasses
   out. No file access, no global state, no configuration lookups (the runner
   passes parameters explicitly). Each step is unit-testable on synthetic
   arrays and could later be wrapped by another framework.
3. **Pipeline and I/O** — ingest (classification through the header map), the
   runner (chains steps per recipe, caches master calibrations), product
   writers and loaders, quality-assurance figures and pages, the command-line
   interface.

The **PypeIt adapter** is a single module with lazy imports; the package
installs and runs without PypeIt. Any use of PypeIt goes through it, behind
the same function signatures as the native step, so it can be swapped.

A **reference notebook** under `notebooks/` reproduces the notebook's interactive
workflow by importing the package. The logic lives in the package, never in
notebooks.

### 3.1 Repository layout

```
mastpipe/
  pyproject.toml
  CLAUDE.md  README.md  LICENSE  CHANGES.md  DECISIONS.md
  .github/workflows/ci.yml          # Linux, Windows, macOS matrix
  src/mastpipe/
    __init__.py  __version__.py
    instruments/
      __init__.py  registry.py      # name -> InstrumentConfig
      highspec/  config.yaml  header_map.yaml  layouts/*.csv  gratings/*.yaml
      deepspec/  config.yaml  header_map.yaml  layouts/*.csv  bands/*.yaml
    config/       models.py  loader.py          # Pydantic models, layered merge
    data/         atlases/  linelists/           # public reference data only
    steps/
      calib.py    trace.py   extract.py  lines.py  widths.py  focus.py
      wavecal.py  drift.py   sky.py      combine.py  qc.py
    io/           ingest.py  manifest.py  products.py  loaders.py  qa.py
    pipeline/     runner.py  recipes.py  callib.py   provenance.py
    ext/          pypeit_adapter.py
    sim/          simulator.py  spectra.py
    cli.py
  tests/          unit/  e2e/  regression/  conftest.py
  notebooks/      highspec_reference.ipynb
  docs/           design/  (this document, the implementation plan)
```

### 3.2 Core data model (dataclasses, all serialisable to FITS)

- `Frame` — 2-D data, variance, mask, plus a `FrameMeta` record (instrument,
  camera, configuration, exposure type, exposure time, mid-time, binning,
  readout mode, temperature, lamp state, fiber-stage preset, checksum,
  source filename).
- `MasterCalib` — a calibrated 2-D product with validity keys and provenance.
- `Trace` / `TraceSet` — per trace: index, role, telescope id, polynomial
  coefficients, aperture half-width, centroid samples with clipping mask,
  fit rms, relative throughput.
- `Spectrum` — per fiber: pixel, flux, error, mask, optional wavelength,
  optional subtracted sky, plus the fiber's role, telescope id, target id.
- `ExposureSpectra` — all `Spectrum` objects of one exposure with the
  exposure's `FrameMeta` and wavelength status.
- `WavelengthSolution` — per fiber: kind (`grating` or `poly`), geometry
  parameters, distortion coefficients, identified lines with residuals, rms,
  line count, wavelength map.
- `DriftModel` — per science exposure: per-arc-fiber shift and stretch, the
  smooth model across the detector, and its residuals.
- `CombinedSpectrum` — wavelength grid, flux, error, mask, contributing
  exposures and fibers, per-telescope barycentric corrections, exclusions with
  reasons, the grouping rule.
- `QCReport` — the scalar quality numbers of a stage, serialisable to JSON.

## 4. Ingest and configuration

### 4.1 Header map

One YAML file per instrument maps the pipeline vocabulary to raw keywords.
Vocabulary: `instrument`, `camera`, `configuration` (grating or band),
`exposure_type`, `exposure_time`, `mid_time`, `binning`, `readout_mode`,
`detector_temperature`, `lamp_state`, `fiber_preset`, `target`, `telescope
pointing` (per telescope, when available). Each entry lists candidate keywords
in priority order, optional value transforms (unit conversion, enumerations),
and a filename regular expression as last resort. Unresolved fields are
recorded as `unknown`; nothing is guessed silently.

### 4.2 Frame classification

`exposure_type` comes from the header when present. Otherwise explicit,
configurable rules apply (zero exposure time → bias; shutter closed → dark;
lamp on → arc; continuum lamp → flat; else science). Every decision is logged
with the rule that made it. A **manifest** table (CSV) can override any field
frame by frame; it is the way laboratory data without instrument keywords
enter the pipeline.

### 4.3 Multi-camera frames

DeepSpec frames arrive as multi-extension FITS. Ingest splits them into one
`Frame` per camera immediately. Downstream code never sees extensions; each
camera is a detector with its own band configuration.

### 4.4 Layered configuration

Merged in order, later wins: instrument defaults → configuration block
(grating or band) → optional site/run file → command-line overrides. Validated
by Pydantic models. The configuration block carries the geometry constants
(groove density, incidence angle, diffraction order, camera focal length,
pixel size, central wavelength, dispersion sign), the atlas file, band limits,
and per-step parameters (aperture half-width, detection thresholds, polynomial
degrees, clipping limits). The effective configuration is hashed for
provenance.

### 4.5 Fiber layout table

One CSV per instrument and fiber-stage preset: `trace`, `role`
(`arc`/`science`/`sky`), `telescope`, `y_nominal` (position at detector
centre), `slope_nominal`, `aperture` (optional override). Trace finding
matches detected traces to this table allowing a global shift and small
per-trace offsets, and reports missing or extra traces. The number of traces
is never inferred from data. The three-fiber laboratory setup and the
45-trace survey setup are two files.

### 4.6 Calibration library

Master frames are stored on disk with validity keys (date, temperature,
binning, readout mode, configuration) and matched to science frames on those
keys; a mismatch is an error in manual mode and a logged warning in
quick-look. Arcs are matched by nearest mid-time within the same
configuration.

## 5. Calibration and wavelength flow

1. **Detector calibration.** Master bias (sigma-clipped stack), master dark
   rate, background model (the notebook's method), applied per camera. Variance is
   built from gain and read noise. Bad-pixel mask from the detector
   parameters plus outliers in the master frames.
2. **Master flat.** Used for three purposes: trace finding, per-fiber relative
   throughput (integrated flat flux per trace, normalised), and the spatial
   profile for optimal extraction. Pixel-to-pixel flat-fielding within the
   apertures is implemented behind a switch, off by default.
3. **Traces.** Collapse the flat along the dispersion axis, detect peaks, match
   to the layout table, fit each trace with a quadratic (captures the
   conical-diffraction curvature of the outer fibers). Products: `TraceSet`
   with per-trace rms, the diagnostic to watch.
4. **Extraction.** Two extractors with one signature: aperture sum with
   fractional-pixel weights (the notebook's), and Horne optimal extraction using the
   flat-derived spatial profile, with cosmic-ray rejection from profile
   residuals. Both propagate variance and mask.
5. **Full-arc wavelength solution** (every fiber sees the lamp). Per fiber:
   continuum subtraction (asymmetrically clipped Chebyshev), line finding with
   deblending, matching to the atlas convolved to the instrument resolution,
   initial solution from the grating equation and geometry constants,
   automatic anchoring by cross-correlation against the convolved atlas,
   fallback to a stored per-configuration table of previously identified lines,
   then a low-order distortion refinement with sigma clipping. The incidence
   angle is held fixed (it is degenerate with the plate scale). Air wavelengths
   throughout. Product: `WavelengthSolution` per fiber, with the
   fiber-to-fiber comparison as the external accuracy check.
6. **Transfer to science fibers.** In a science exposure only the arc fibers
   see the lamp. Their lines are measured against the reference full-arc
   solution, giving a shift and, if significant, a stretch per arc fiber. A
   smooth low-order model of the drift across the detector is fitted through
   the arc-fiber positions and applied to every science and sky fiber's
   reference solution. The wavelength status of the exposure becomes
   `transferred`. If no arc fibers are usable the status is `reference` (last
   full arc) or `nominal` (geometry only), and the products say so.
7. **DeepSpec.** Same interfaces. Low resolution and wide band mean few lines
   per camera; template cross-correlation is the primary identification
   method, with the grating-equation bootstrap from geometry when available.
   The PypeIt adapter's automated line identification is the first candidate
   for use here.
8. **Sky subtraction.** Each science fiber is paired with the sky fiber of the
   same telescope (from the layout table). The sky spectrum is scaled by the
   relative throughput, resampled onto the science fiber's wavelength grid,
   and subtracted; the subtracted sky is kept in the product. Subtraction is
   in wavelength space because the two fibers have different solutions.
9. **Combination.** See §7.
10. **Focus recipe.** the notebook's line-width analysis (four width models: Gaussian,
    top-hat, disk-convolved-Gaussian, enclosed energy) and focus-scan fitting
    (global parabola, core refit, per-fiber, per-line, focal-plane tilt) run as
    a separate diagnostic recipe over a scan folder, independent of the
    control software's autofocus.

Conventions: air wavelengths, with vacuum conversion at output as an option;
barycentric corrections computed per telescope and stored, never applied to
the wavelength grid.

## 6. Run modes

Both modes use the same runner and steps; they differ in recipe, configuration
layer and failure policy.

**Quick-look.** `mastpipe quicklook <night-folder>` or `mastpipe watch
<folder>`. Non-interactive. Uses the most recent valid calibration-library
entries; degrades gracefully (wavelength status `reference`/`nominal`) rather
than stopping. One bad frame is logged and skipped. Idempotent, keyed on frame
checksums. Produces per-exposure spectra and quality page, and a per-night
summary page. Designed to be launched by the observatory's supervisor; that
wiring is out of scope.

**Manual.** `mastpipe reduce --config run.yaml [--manifest frames.csv]`, plus
the Python API and notebooks. `mastpipe identify` plots a spectrum against the
template and reads anchor pairs from a small text file, replacing the
notebook's by-eye step. Re-runs resume from cached products; `--force <stage>`
redoes a stage. Errors raise immediately.

**Provenance.** Every product header carries the pipeline version, the
configuration hash, input-frame checksums, the calibration files used, and the
wavelength status.

**Logging.** Structured per-run log file next to the products, plus a console
summary.

## 7. Products

**Layout.** One tree per night, by instrument and configuration:
`calib/`, `spec/`, `qa/`, `logs/`. Filenames are built from ingest metadata
(instrument, configuration, UTC mid-time, product kind), never from raw
filenames.

**Files (FITS, astropy-written).**

1. Master calibrations: images with validity keys and provenance.
2. Trace set: binary table, one row per trace.
3. Wavelength solution: per-fiber table, identified-lines table, wavelength
   map image (fibers × pixels); science exposures add a drift table.
4. Extracted spectra: one file per raw exposure, one table with columns
   `fiber, role, telescope, target, wavelength, flux, error, mask, sky`;
   header states the wavelength status.
5. Combined spectra: one file per target and group: `wavelength, flux, error,
   mask`, per-telescope barycentric corrections, contributing exposures and
   fibers, exclusions with reasons, the grouping rule.

**Combination options.**

- Across telescopes within an exposure: on/off.
- Grouping of exposures in time: `exposure` (none), `count: N` (every N
  consecutive exposures), `window: <minutes>`, `night`, or an explicit list
  from the manifest.
- Telescope subset: include or exclude lists; default is all telescopes with a
  science fiber on the target.
- Automatic rejection, **off by default**: each exposure and telescope
  spectrum is scored on signal-to-noise per fiber, arc-fiber drift, trace
  matching failures, and flux relative to the median of the group; anything
  beyond configurable thresholds is excluded and listed with the triggering
  metric. Excluded spectra remain in the per-exposure files.
- Resampling onto a common logarithmic wavelength grid; weighted mean with
  outlier rejection.

**Quality assurance.** PNG figures per stage (the notebook's figure set, mapped to
stages), an HTML page per exposure and per night with trace residuals, lines
used and wavelength rms per fiber, drift, resolving power and signal-to-noise
per fiber, and a JSON file per exposure with the same numbers for later
dashboard use.

**Access.** Loaders return the same dataclasses the steps produce; `mastpipe
show <product>` prints a summary.

## 8. Testing

**Simulator.** `mastpipe.sim` builds a synthetic night from an instrument
configuration, a layout table and an atlas: bias, dark, flat, arc and science
frames, with a model stellar spectrum in science fibers, sky in sky fibers,
detector noise, cosmic rays, trace curvature and a controllable wavelength
drift. Deterministic seeds. It is the first deliverable after the skeleton and
the only way to exercise the 45-trace layout before survey data exist.

**Tiers.**
1. Unit tests per step on small synthetic arrays (the bulk of the suite).
2. End-to-end: quick-look on a synthetic night, every product validated
   against the simulator's truth.
3. Regression against the notebook's published numbers on the laboratory data
   (per-fiber wavelength rms, resolving power, best focus of the eleven-frame
   scan). Data located via `MASTPIPE_TEST_DATA`; tests skip when absent.

**CI.** GitHub Actions matrix on Linux, Windows and macOS running tiers 1–2,
plus lint (ruff) and formatting checks.

## 9. Development workflow and build order

- Design and implementation both by Claude Fable 5.1 (user preference for this
  project; no downgrade to Sonnet unless asked). Test-first, small tasks, each
  reviewed before the next.
- Build order: package skeleton and configuration model → simulator → detector
  calibration and traces → extraction → full-arc wavelength solution →
  arc-fiber transfer → sky and combination → products and quality pages →
  quick-look command → DeepSpec configuration and regression on available
  frames → focus recipe port → PypeIt adapter (only where it earns its place).
- The repository is created only when the user says so, as a **private**
  repository under the MAST GitHub organisation, replacing the stale fork
  that currently holds the name. No push without per-push approval.
- `DECISIONS.md` in the repository records design decisions as they are made.

## 10. Repository hygiene

- No local absolute paths, hostnames, IP addresses, usernames, e-mail
  addresses, tokens or credentials anywhere in tracked files. A pre-commit
  scan enforces this (patterns for common path roots, e-mail and IP formats).
- Real data never enter git; `tests/` use synthetic data only.
- The repository `CLAUDE.md` (drafted alongside this document) sets the
  conventions, privacy and security rules for agents working in it.

## 11. Open questions and risks

- **DeepSpec specifics** (arc lamp species, per-band geometry, fiber layout)
  are not yet known; the DeepSpec configuration will be filled when the
  instrument team provides them, and the laboratory frames from 2024 are not a
  reference.
- **Trace density.** 45 traces on a 400-row detector leaves about 9 rows per
  trace; the aperture half-width and deblending windows from the three-fiber
  work will need retuning. The simulator is the tool to find the limits.
- **Arc-fiber transfer accuracy.** Whether a smooth drift model across five
  arc fibers reaches the few-km/s goal is an empirical question; the
  fiber-to-fiber comparison figure is the acceptance test.
- **Windows CI** may surface path and line-ending issues early; that is the
  point of running it from the first commit.
- **PypeIt adapter validity.** PypeIt core routines are called as plain
  functions; their signatures may change between releases, so the adapter
  pins a version and is covered by its own tests.
