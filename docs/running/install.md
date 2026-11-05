# Install

mastpipe is a pure-Python package. It runs on Linux, macOS and Windows with
Python 3.11 or newer, and depends only on numpy, scipy, astropy, pandas,
matplotlib, pydantic and PyYAML.

## Installing the pipeline

The package is not on PyPI or conda-forge yet; pip installs it straight from
the repository. Use a fresh environment, either a `venv` or a conda
environment:

::::{tab-set}

:::{tab-item} venv
```bash
python -m venv mastpipe-env
source mastpipe-env/bin/activate   # on Windows: mastpipe-env\Scripts\activate
pip install "git+<repository>"
mastpipe --version
```
:::

:::{tab-item} conda
```bash
conda create -n mastpipe python=3.12
conda activate mastpipe
pip install "git+<repository>"
mastpipe --version
```
Inside a conda environment pip installs the package and its dependencies the
same way; to have conda manage the scientific stack instead, create the
environment with the dependencies from conda-forge first and then install the
package alone:
```bash
conda create -n mastpipe -c conda-forge python=3.12 numpy scipy astropy pandas matplotlib pydantic pyyaml
conda activate mastpipe
pip install --no-deps "git+<repository>"
```
:::

::::

The same command with `pip install --upgrade` moves an existing installation
to the current `main`. From a downloaded or cloned copy of the repository,
`pip install .` run in its top folder does the same without touching git.

`mastpipe --version` printing a version and `mastpipe --help` listing the four
commands (`identify`, `quicklook`, `solve`, `show`) is the whole installation
check for this kind of install; the test suite is part of the repository, not
of the package, and needs the developer install below. The real check is the
first quick-look, which the next pages walk through.

## Developer install

To work on the code, run the tests or build this documentation, clone the
repository and install it in editable mode with the extras:

::::{tab-set}

:::{tab-item} venv
```bash
git clone <repository>
cd mastpipe
python -m venv .venv
source .venv/bin/activate          # on Windows: .venv\Scripts\activate
pip install -e ".[dev]"
mastpipe --version
```
:::

:::{tab-item} conda
```bash
git clone <repository>
cd mastpipe
conda create -n mastpipe python=3.12
conda activate mastpipe
pip install -e ".[dev]"
mastpipe --version
```
:::

::::

The `dev` extra adds the test and lint tools (`pytest`, `pytest-cov`, `ruff`,
`nbstripout`); the `docs` extra builds this documentation (`pip install -e
".[dev,docs]"`); the `pypeit` extra pre-installs the dependency of a planned
PypeIt adapter; no adapter module exists yet.

### Checking the installation

The tests live in the repository's `tests/` folder and are not installed with
the package, so run them from the root of the clone (the folder holding
`pyproject.toml`), with the environment activated:

```bash
cd mastpipe
pytest -m "not slow"   # the unit tests, about three minutes
pytest                 # the whole suite incl. simulated nights, about six minutes
```

Run from anywhere else, a bare `pytest` finds no tests (`collected 0 items`,
`no tests ran`); point it at the folder instead, as in
`pytest path/to/mastpipe/tests`, and it picks up the repository's configuration
from there.

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
