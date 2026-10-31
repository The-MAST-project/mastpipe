# Install

mastpipe is a pure-Python package. It runs on Linux, macOS and Windows with
Python 3.11 or newer, and depends only on numpy, scipy, astropy, pandas,
matplotlib, pydantic and PyYAML.

```bash
python -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
pip install -e ".[dev]"
mastpipe --version
```

The `dev` extra adds the test and lint tools (`pytest`, `pytest-cov`, `ruff`,
`nbstripout`); the `docs` extra builds this documentation; the `pypeit` extra
installs PypeIt for the optional adapter, which this release does not use.

## Checking the installation

```bash
pytest -m "not slow"   # the unit tests, about three minutes
pytest                 # the whole suite incl. simulated nights, about six minutes
```

The suite is offline and synthetic. Two optional environment variables point
at real data on your machine; the tests that need them skip when they are
unset:

`MASTPIPE_TEST_DATA`
: The folder holding the HighSpec laboratory frames, laid out as
  `highspec/lab3/calib/` (bias, background, flat) and `highspec/lab3/scan/`
  (the Th-Ar focus scan). With it set, `pytest tests/regression` reproduces
  the laboratory reduction in a few seconds.

`MASTPIPE_ATLAS_DIR`
: The folder holding the Th-Ar reference atlas of Murphy et al.
  (`thar_spec_MM201006.fits`, `thar_MM201006.dat`). `identify` uses it for the
  anchor panels; without it the panels show a synthetic template built from the
  shipped line table, which is enough to read anchors off.

The reduction itself never reads either folder; the tests, the line-table
script and `identify`'s anchor panels (by default from `MASTPIPE_ATLAS_DIR`, or
`--atlas-dir`) do.

## Where things go

Every command writes under the folder you pass with `--out`:

```text
out/
  calib/   masters, trace sets, reference solutions, the library index
  spec/    per-exposure spectra and per-target combined spectra
  qa/      one QC JSON per trace set, reference, drift arc and science frame
  logs/    one log per quicklook run
```

`out/calib` is also the **calibration library** that later runs read from, so
keep one output folder per instrument and reuse it from night to night. The
layout is described in [Products](products.md).
