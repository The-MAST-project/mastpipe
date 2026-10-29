# mastpipe

Data-reduction pipeline for the MAST fiber-fed spectrographs, **HighSpec**
(R ≈ 20 000, narrow-band gratings, one detector) and **DeepSpec** (R ≈ 500,
four cameras). It takes raw frames to wavelength-calibrated, sky-subtracted
spectra per fiber, combined spectra per target, and quality-assurance
products. It runs as an unattended quick-look at the observatory and as a
manual reduction tool.

The design documents live outside the repository; the decision log is
`DECISIONS.md` and the release history is `CHANGES.md`.

## Installation

Python 3.11 or newer is required.

```bash
python -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
pip install -e ".[dev]"
```

The `dev` extra installs `pytest`, `pytest-cov`, `ruff` and `nbstripout`.
The optional `pypeit` extra installs PypeIt for the adapter module; the
package runs without it.

## Running the tests

```bash
pytest
```

The default suite is offline and uses synthetic data only. Tests that need real
laboratory frames or a reference line atlas skip cleanly unless the environment
variables below point at the data.

## Checks run in CI

```bash
ruff check src tests scripts
ruff format --check src tests scripts
python scripts/privacy_scan.py
pytest
```

The privacy scan fails when any tracked text file contains a local absolute
path, an e-mail address, an IPv4 address (other than `0.0.0.0` and
`127.0.0.1`) or a credential token.

## Environment variables

Real observatory data are never stored in the repository. Two environment
variables locate them on each machine; both are optional, and the tests that
need them are skipped when they are unset.

### `MASTPIPE_TEST_DATA`

The folder holding real test frames. The expected layout, relative to that
folder, is:

```text
highspec/
  lab3/
    calib/    laboratory calibration frames (bias, dark, flat, arc)
    scan/     focus-scan frames
```

### `MASTPIPE_ATLAS_DIR`

The folder holding the Th-Ar reference atlas files:

```text
thar_spec_MM201006.fits    the reference spectrum
thar_MM201006.dat          the line list
```

## Command line

```bash
mastpipe --version
mastpipe --help
```

Sub-commands are added as the pipeline stages land; see `CHANGES.md`.

## License

BSD 3-Clause; see `LICENSE`.
