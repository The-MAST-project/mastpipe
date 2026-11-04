# mastpipe

[![CI](https://github.com/The-MAST-project/mastpipe/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/The-MAST-project/mastpipe/actions/workflows/ci.yml)
[![Documentation](https://readthedocs.org/projects/mastpipe/badge/?version=latest)](https://mastpipe.readthedocs.io/en/latest/)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://github.com/The-MAST-project/mastpipe/blob/main/pyproject.toml)
[![License: BSD-3-Clause](https://img.shields.io/github/license/The-MAST-project/mastpipe.svg)](https://github.com/The-MAST-project/mastpipe/blob/main/LICENSE)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Data-reduction pipeline for the MAST fiber-fed spectrographs, **HighSpec**
(R ≈ 20 000, narrow-band gratings, one detector) and **DeepSpec** (R ≈ 500,
four cameras). It takes raw frames to wavelength-calibrated, sky-subtracted
spectra per fiber, combined spectra per target, and quality-assurance
products. It runs as an unattended quick-look at the observatory and as a
manual reduction tool.

**Documentation:** <https://mastpipe.readthedocs.io/en/latest/> (observer guide, how it
works, developer guide, API reference and executed tutorials; built from `docs/`).

The design documents are kept locally under `docs/design/` (git-ignored);
the decision log is `DECISIONS.md` and the release history is `CHANGES.md`.
The reduction steps
are a port of the HighSpec laboratory reduction notebook,
written as pure functions on arrays behind a declarative instrument layer.

## Installation

Python 3.11 or newer is required. The package is not on PyPI yet; pip
installs it from the repository into a fresh environment:

```bash
python -m venv mastpipe-env
source mastpipe-env/bin/activate   # on Windows: mastpipe-env\Scripts\activate
pip install "git+<repository>"
mastpipe --version
```

or, with conda,

```bash
conda create -n mastpipe python=3.12
conda activate mastpipe
pip install "git+<repository>"
mastpipe --version
```

(`pip install .` from a downloaded copy does the same.) The optional `pypeit`
extra installs PypeIt for the adapter module; the package runs without it.

### Developer install

To work on the code or run the tests, clone and install in editable mode:

```bash
git clone <repository>
cd mastpipe
python -m venv .venv                # or: conda create -n mastpipe python=3.12
source .venv/bin/activate           #     conda activate mastpipe
pip install -e ".[dev]"
```

The `dev` extra installs `pytest`, `pytest-cov`, `ruff` and `nbstripout`;
the `docs` extra builds the documentation.

## Running the tests

From the root of the clone (the tests are not installed with the package; from
another folder pass the path, `pytest path/to/mastpipe/tests`):

```bash
pytest                 # the whole suite, about six minutes
pytest -m "not slow"   # the unit tests alone, about three minutes
```

The suite is offline and uses synthetic data (a simulated HighSpec night) by
default. Tests that need the real laboratory frames or the Th-Ar reference
atlas skip cleanly unless the environment variables below point at the data;
with them set, `pytest tests/regression` reproduces the notebook's results on
the laboratory frames (trace centres, per-fiber solution rms and line counts,
band coverage, fiber-to-fiber agreement) in a few seconds.

## Checks run in CI

```bash
ruff check src tests scripts
ruff format --check src tests scripts
python scripts/privacy_scan.py
pytest
```

The privacy scan fails when any tracked text file contains a local absolute
path, an e-mail address, an IPv4 address (other than `0.0.0.0` and
`127.0.0.1`) or a credential token. Real data never enter the repository: the
tests reach them only through the environment variables.

## Environment variables

Two optional variables locate the real data on each machine. The tests that
need them are skipped when they are unset.

### `MASTPIPE_TEST_DATA`

The folder holding the laboratory frames. The expected layout, relative to that
folder, is:

```text
highspec/
  lab3/
    calib/    image_bias.fits, background.fits, flat_518.fits (and image_dark.fits)
    scan/     FOCUS_00_13.05001_mm.fits … FOCUS_10_13.55007_mm.fits (the Th-Ar focus scan)
```

`tests/regression/test_lab_highspec.py` uses the three calibration frames and
`FOCUS_05_13.30004_mm.fits`, the frame at best focus, as the full arc; the
other scan frames are not read. The frames carry the Andor Solis keywords
(`EXPOSURE`, `TEMP`, `HBIN`, `VBIN`, `READMODE`, `DATE`) but no frame type,
lamp, grating or fiber preset, so the test writes a manifest for them (see
`quicklook --manifest` below).

### `MASTPIPE_ATLAS_DIR`

The folder holding the Th-Ar reference atlas of Murphy et al.:

```text
thar_spec_MM201006.fits    the reference spectrum
thar_MM201006.dat          the line list
```

It is read by `scripts/make_line_table.py`, by the real-atlas test in
`tests/unit/test_atlas.py` and by `identify` for the anchor panels (the
`--atlas-dir` default); without it the panels use a synthetic template built
from the shipped line table.

## Wavelength calibration in three cadences

A **full arc**, taken every few months and at every grating change, lights
every fiber; `identify` (the first time, from two anchor pairs) or `solve`
(later, seeded by the previous reference) fits an independent grating-equation
solution per fiber and stores the result as a **reference solution**, valid
from the arc's mid-time until the next one. A **calibration-fiber arc**, taken
one to three times a night and at every instrument switch, lights the arc
fibers only; the quick-look measures the shift and stretch of its lines in
pixel space against the reference valid at its time and fits them as a
polynomial across the detector, the **drift model**. **Science exposures**
carry no arc light: each takes its wavelengths from the reference valid at its
mid-time, corrected by the drift model of the nearest calibration-fiber arc
before it (`drift.time_policy`: `nearest_before`, `nearest` or `interpolate`,
within `drift.max_age_hours`), and its product records the status
(`transferred`; `reference` when no arc pairs; `nominal`, from the configured
geometry, when no reference exists) together with the reference identifier and
the drift arcs' checksums. HighSpec uses the grating formula and DeepSpec the
dispersion equation behind one model interface.

## Command line

```bash
mastpipe --version
mastpipe --help
mastpipe <command> --help
```

The reduction commands (`identify`, `solve`, `quicklook`) take `--instrument`
(e.g. `highspec`), `--configuration` where the frames' headers do not say which
grating was used (e.g. `Mg`), `--camera` when the instrument has several,
`--layout` to name the fiber preset (`lab3`, `survey45`; default: the frames'
preset, else the instrument's), `--out` for the output folder, `--library` for
a calibration library elsewhere than `OUT/calib`, and repeatable
`--set group.key=value` overrides of any step parameter
(`trace.trace_halfwidth=4`, `drift.time_policy=interpolate`, …). `show` takes
the product file alone. The configuration actually used is hashed into every
product.

Parameters are layered: the instrument's defaults, then the grating block,
then the **layout preset's block** (`layout_params` in `config.yaml`; for
`survey45`, whose 45 fibers sit at an 8.5 px pitch, the trace, line and
extraction halfwidths narrow to 4 px, and the layout's `aperture` column fixes
the extraction halfwidth per fiber), then the `--set` overrides. Choosing a
layout therefore selects its parameters without any override.

Fibers are identified by the layout's `trace` id throughout: traces are
matched to the layout and re-keyed by that id, detections that match no layout
row are dropped (logged, and listed as `extra` with their positions in the
trace QC), and a fiber missing on a night is simply absent from its products.
A science fiber whose trace centre lies farther than `trace.layout_tolerance`
from the reference solution's is given no wavelength and listed under
`wavelength_skipped` in the science QC.

### `identify` — bootstrap a reference solution

Without anchors, the command reduces a full arc to its lines and draws the
two-panel anchor figure (template spectrum above, measured spectrum below,
strong lines numbered in both) for the user to read a few wavelength–pixel
pairs off:

```bash
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg
#   anchor panels: out/anchor_panels.png
#   candidates: out/anchor_candidates.csv
```

With anchors — given on the command line (`--anchors`), read from a YAML file
(`--anchors-file`, with `pairs: [[nm, px], ...]` and `trace:`), or the pairs
stored with the instrument for the configuration and layout
(`--stored-anchors`) — it solves every fiber from them and writes the
reference; `--trace` names the layout trace id the pairs were measured on:

```bash
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg \
    --anchors 515.8605:267.22,525.836:1423.22 --trace 0
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg \
    --stored-anchors
#   reference highspec_Mg_lab3_20260909T113659: 3 fibers solved
#   written: out/calib/highspec_Mg_lab3_20260909T113659_wref.fits
```

Masters already in `out/calib` (`*_bias.fits`, `*_background.fits`,
`*_flat.fits`) are applied; traces come from the flat master when there is
one and from the arc itself otherwise. `fiber_agreement.csv` beside the
reference lists, per catalogue line, how far the independent per-fiber
solutions disagree.

### `solve` — re-solve a full arc from the previous reference

```bash
mastpipe solve night/arc_full.fits --out out --instrument highspec --configuration Mg
#   reference highspec_Mg_lab3_20261112T010203: 3 fibers solved, seeded by highspec_Mg_lab3_20260909T113659
```

The seed is the latest reference of the configuration and layout in the
calibration library (`--library`, default `out/calib`); the command exits with
status 2 and "run identify first" when there is none. The new reference does
not replace the old one: each is valid from its own arc's time.

### `quicklook` — reduce a night folder

```bash
mastpipe quicklook night/ --out out --instrument highspec
mastpipe quicklook night/ --out out --instrument highspec --manifest night/manifest.csv \
    --time-policy interpolate --on-error raise
```

Frames are classified from their headers (`IMAGETYP`, `LAMP`, exposure time)
and the folder is reduced per camera and configuration: masters are stacked
into the library (a kind the night lacks is taken from it), traces are found on
the flat master (or a full arc), full arcs are re-solved from the library's
reference, calibration-fiber arcs are fitted for drift, and science frames are
calibrated with cosmic rays flagged, extracted, given wavelengths by the time
policy, sky-subtracted with the telescope's sky fiber, passed through the flux
and telluric stages, written with provenance and QC, and combined per target.
A product already on disk with the same idempotence key (input checksum,
configuration hash, reference identifier, drift-arc checksums, layout preset)
is skipped; one whose key changed is reduced again; `--force` reduces
everything. `--reference ID` pins one reference solution instead of the one
valid at each frame's time. Frames the headers cannot classify are skipped
with their rule, failures are recorded and the run goes on (`--on-error
continue`, the default). The exit status is 2 when anything failed, when the
night folder holds no readable frame, or when it does not exist.

Frame times: the header map declares which instant a time keyword records
(`reference: start | mid | end`; HighSpec's `DATE-OBS` / `DATE` is the exposure
start) and ingest derives each frame's mid-time from it and the exposure time;
a frame without an exposure time keeps the header time as given and notes
`mid_time` among its unresolved fields. A header value listed under a field's
`missing_values` (the Andor `TEMP = -999` of an unread sensor) counts as
absent.

The night folder is scanned with its subfolders (the control software files
each acquisition in a folder of its own); every `*.fits`, `*.fit` or
`*.fits.fz` file (any case) is read, names starting with `.` are ignored, and a
frame is known by its path relative to the night folder. A frame whose
`INSTRUME` names another instrument is skipped with that reason.

A **manifest** (`--manifest`) is a CSV keyed by `source` (optional `extension`)
whose non-empty cells override any metadata field. A `source` is a relative
path (`highspec/acquisition-0001/exposure-001.fits`), a glob pattern on it
(`highspec/acquisition-0001/*.fits`; `*` stays within a folder, `**` crosses
folders) or a bare file name that only one frame has (two or more is an error);
a row matching no frame is reported among the skipped. The fields are: `exposure_type`
(`bias`, `dark`, `flat`, `arc`, `science`, `background`, `unknown`),
`arc_kind` (`full`, `calfiber`), `configuration`, `fiber_preset`, `target`,
`exposure_time`, `mid_time` (taken as the mid-exposure time, UTC unless it
carries an offset), `binning`,
`readout_mode`, `detector_temperature`, `lamp_state`, `camera`, `instrument`,
`classification_rule`. It is how frames without instrument keywords, such as
the laboratory data, enter the pipeline, and how a night names its fiber layout
when `--layout` does not (the headers carry none):

```csv
source,exposure_type,configuration,fiber_preset,arc_kind
image_bias.fits,bias,Mg,lab3,
background.fits,background,Mg,lab3,
flat_518.fits,flat,Mg,lab3,
FOCUS_05_13.30004_mm.fits,arc,Mg,lab3,full
```

### `show` — print what a product holds

```bash
mastpipe show out/spec/highspec_Mg_main_20261001T221500_3f9a1c_spec.fits
mastpipe show out/calib/highspec_Mg_lab3_20260909T113659_wref.fits
```

Prints the product kind, the frame and its statuses, a per-fiber summary
(role, telescope, target, median S/N, masked fraction, kept versions, drift)
or the per-fiber solutions of a reference, the traces, the contributors of a
combined spectrum or the validity of a master, followed by the provenance:
pipeline version, configuration hash, inputs, calibrations, drift arcs and the
steps applied.

## Products

Everything is written under the output folder:

```text
out/
  calib/   masters, trace sets, reference solutions, fiber_agreement.csv, index.json
  spec/    per-exposure spectra and per-target combined spectra
  qa/      one QC JSON per trace set, reference, drift arc and science frame
  logs/    quicklook_<UTC stamp>.log, one per run
```

`out/calib` is also the **calibration library**: `index.json` lists every
master (camera, binning, readout mode, configuration, temperature, time) and
reference (identifier, `valid_from`) so later runs find them. Files are named
`{instrument}_{configuration}_{camera}_{YYYYmmddTHHMMSS}_{checksum}_{kind}.fits`
after the frame's UTC mid-time and the first six hex digits of its file's
SHA-256 (so equal mid-times never collide), with `kind` one of `bias`, `dark`, `background`,
`flat`, `traces`, `spec` or `comb_<target>`; references are
`{instrument}_{configuration}_{preset}_{stamp}_wref.fits`, the stamp being the
full arc's mid-time, and that name without the suffix is the reference
identifier.

Every FITS product carries its provenance as `MP*` header cards: pipeline
version (`MPVER`), configuration hash (`MPCFGH`), the wavelength, flux and
telluric statuses, the reference identifier and hash (`MPWREF`, `MPWREFH`),
the layout preset (`MPLAYOUT`), and numbered lists of inputs, calibrations,
drift arcs and steps. No raw instrument keyword reaches a product header.

### The spectra table

A science product holds one table, `SPECTRA`, with **one row per fiber** on
one pixel grid:

| column | content |
| --- | --- |
| `FIBER`, `ROLE`, `TELESCOPE`, `TARGET` | the fiber's layout identity (`science`, `sky`, `arc` roles) |
| `PIXEL`, `FLUX`, `ERROR`, `MASK` | the raw aperture extraction, never overwritten |
| `WAVELENGTH` | the fiber's wavelengths by the time policy (`HAS_WAVELENGTH` per row) |
| `SKY`, `CONTINUUM` | the sky spectrum resampled onto the fiber, the fitted continuum |
| `FLUX_SKYSUB`, `ERROR_SKYSUB`, `MASK_SKYSUB` | the sky-subtracted version |
| `FLUX_FLUXCAL`, …, `FLUX_TELLURIC`, … | flux-calibrated and telluric-corrected versions, once those stages have calibrations |

Each **kept version** is a set of three columns on the same row, present only
when its stage ran (`HAS_<version>` per row, presence cards on the table), so
a reader finds the raw extraction and every later stage of a fiber side by
side and the table keeps one row per fiber. A `DRIFT` table (per arc fiber:
shift, stretch, lines, rms, usable flag, blend weight) is attached when a
drift model was applied. Combined products hold `COMBINED` (`WAVELENGTH` and
`FLUX_<V>`, `ERROR_<V>`, `MASK_<V>` for `RAW` and every shared version),
`CONTRIB` and `EXCLUDED`. Reference products hold `SOLUTIONS` (one row per
fiber: model kind and parameters, rms, line count, anchors, trace centre) and
`LINES` (every identified line with its residual).

## License

BSD 3-Clause; see `LICENSE`.
