# mastpipe Plan 1 — Foundation: raw frames to extracted spectra on a synthetic night

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A pip-installable `mastpipe` package that ingests a night of frames, builds master calibrations, finds and identifies traces against a fiber-layout table, aperture-extracts every fiber with errors and masks, writes FITS products with provenance, and does all of that end-to-end on a simulated night through `mastpipe quicklook`.

**Architecture:** Three strictly separated layers (see the spec §3): declarative instrument configuration, pure step functions on arrays and dataclasses, and a pipeline/I-O layer (ingest, runner, products, CLI). The HighSpec laboratory reduction notebook functions for calibration, trace finding, aperture extraction and continuum fitting are ported one-to-one into the step layer. A simulator produces synthetic nights with known truth so every step and the end-to-end path are testable before real survey data exist.

**Tech Stack:** Python ≥ 3.11 · numpy · scipy · astropy · pandas · matplotlib · pydantic v2 · pyyaml · pytest · ruff · GitHub Actions (Linux/Windows/macOS matrix).

**Spec:** `docs/design/2026-09-24-mastpipe-design.md` (same folder as this plan). Draft repository `CLAUDE.md`: `docs/design/2026-09-24-mastpipe-claude-md-draft.md`.

## Roadmap (three plans)

| Plan | Delivers | Spec sections |
|---|---|---|
| **1 — Foundation (this plan)** | skeleton, config + data model, HighSpec instrument definition, header-map ingest, simulator, detector calibration, trace finding + layout matching, aperture extraction + continuum, FITS products, calibration library, runner, `quicklook` in pixel space, CI on three platforms | §2, §3, §4, §5.1–5.4, §6 (partial), §7 (products 1, 2, 4), §8 |
| **2 — Wavelength, sky, combination** | full-arc wavelength solution (the notebook's grating-equation method, automatic anchoring, stored line tables), `identify` command, arc-fiber drift transfer, Horne optimal extraction with flat profiles, sky subtraction, combination with grouping/subset/rejection, wavelength-status provenance | §5.4–5.9, §7 (products 3, 5, combination options) |
| **3 — Operations and DeepSpec** | QA figures + HTML/JSON pages, `watch` mode, idempotence polish, DeepSpec instrument definition (multi-camera), focus recipe port, PypeIt adapter, regression tests on lab data | §5.7, §5.10, §6 (quick-look polish), §7 (QA), §8 tier 3, §11 |

Each plan produces working, tested software on its own. Plan 2 and Plan 3 are written after Plan 1 is executed, against the code that then exists.

## Global Constraints

- Python `>=3.11`; hard dependencies only: numpy, scipy, astropy, pandas, matplotlib, pydantic (v2), pyyaml. No new top-level dependency without asking the user.
- Pure Python, cross-platform: `pathlib` only, no shell-outs, no POSIX-only calls, explicit `encoding="utf-8"` on every text open, `newline=""` on CSV writes.
- **No local absolute paths, hostnames, IP addresses, usernames, e-mail addresses, tokens** in any tracked file (code, tests, docs, notebooks, CI). Real test data are located via the environment variable `MASTPIPE_TEST_DATA` and never enter git.
- Tests use synthetic data only (the simulator or tiny hand-made arrays), fixed seeds, no network.
- Step functions (`src/mastpipe/steps/`) never read files, never touch configuration objects, never log to files; they take arrays/dataclasses/plain parameters and return dataclasses.
- No raw FITS keyword name appears outside `src/mastpipe/io/ingest.py` and the instrument `header_map.yaml` files.
- The fiber layout is an input table; the number of traces is never inferred.
- Air wavelengths internally (Plan 2); barycentric corrections stored, never applied.
- Every product header carries: pipeline version, configuration hash, input checksums, calibration files used, wavelength status.
- `src/` layout, package `mastpipe`, console script `mastpipe`.
- Git: local repository only. **No remote, no push** — the GitHub repository is created only when the user says so. Commit after every task; `docs/design/` is git-ignored.
- Ruff for lint and format (line length 100). CI must be green on Linux, Windows and macOS before a task is considered done.
- Writing style in code comments and docs: no first-person AI voice.

## File structure (this plan)

```
mastpipe/
  pyproject.toml            # metadata, deps, console script, ruff, pytest config
  CLAUDE.md                 # copied from docs/design draft
  README.md  CHANGES.md  DECISIONS.md  LICENSE (BSD-3-Clause)
  .gitignore                # already present (docs/design/ ignored)
  .github/workflows/ci.yml  # 3-OS matrix: ruff + pytest
  scripts/privacy_scan.py   # fails on paths / e-mails / IPs / tokens in tracked files
  src/mastpipe/
    __init__.py  __version__.py  cli.py
    core/model.py           # FrameMeta, FrameRef, Frame, MasterCalib, CalibrationSet, Trace, TraceSet, Spectrum, ExposureSpectra, QCReport
    config/models.py        # Pydantic: DetectorConfig, GeometryConfig, CalibParams, TraceParams, ExtractParams, StepParams, ConfigurationBlock, InstrumentConfig, EffectiveConfig
    config/loader.py        # YAML loading, deep merge, layering, config hash
    instruments/__init__.py registry.py
    instruments/highspec/config.yaml header_map.yaml layouts/lab3.csv layouts/survey45.csv gratings/mgb518.yaml
    io/headers.py           # HeaderMap model + resolve(header, filename) -> fields
    io/ingest.py            # scan folder -> FrameRef list; load_frame; classification; manifest
    io/products.py          # FITS writers/readers for MasterCalib, TraceSet, ExposureSpectra; naming
    steps/geometry.py       # grating equation: nominal dispersion and wavelength
    steps/calib.py          # master frames, CalibrationSet, calibrate_frame
    steps/trace.py          # find_traces, match_layout
    steps/extract.py        # aperture_weights, extract_aperture, fit_continuum
    sim/simulator.py        # simulate_night + SimTruth
    pipeline/provenance.py  # provenance header records
    pipeline/callib.py      # calibration library on disk
    pipeline/runner.py      # quicklook recipe
  tests/
    conftest.py             # shared fixtures: tmp night, highspec config
    unit/test_model.py test_config.py test_registry.py test_headers.py test_ingest.py test_geometry.py test_simulator.py test_calib.py test_trace.py test_extract.py test_products.py test_callib.py
    e2e/test_quicklook.py
```

---

### Task 1: Repository skeleton, packaging, CI, privacy scan

**Files:**
- Create: `pyproject.toml`, `src/mastpipe/__init__.py`, `src/mastpipe/__version__.py`, `src/mastpipe/cli.py`, `CLAUDE.md` (copy of `docs/design/2026-09-24-mastpipe-claude-md-draft.md` without its italic draft banner), `README.md`, `CHANGES.md`, `DECISIONS.md`, `LICENSE`, `.github/workflows/ci.yml`, `scripts/privacy_scan.py`, `tests/conftest.py`, `tests/unit/test_cli.py`
- Existing: `.gitignore` (already in the folder; keep it)

**Interfaces:**
- Produces: console script `mastpipe` with `--version`; package importable; `scripts/privacy_scan.py` exit code 1 on a hit.

- [ ] **Step 1: Initialise git and write `pyproject.toml`**

```bash
cd <repo-root>   # the mastpipe folder under code/python
git init -b main
```

`pyproject.toml`:

```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "mastpipe"
description = "Data-reduction pipeline for the MAST fiber-fed spectrographs (HighSpec, DeepSpec)"
readme = "README.md"
license = {text = "BSD-3-Clause"}
requires-python = ">=3.11"
dynamic = ["version"]
dependencies = [
  "numpy>=1.26",
  "scipy>=1.11",
  "astropy>=6.0",
  "pandas>=2.0",
  "matplotlib>=3.8",
  "pydantic>=2.5",
  "pyyaml>=6.0",
]

[project.optional-dependencies]
dev = ["pytest>=8", "pytest-cov", "ruff>=0.6", "nbstripout"]
pypeit = ["pypeit>=2.0"]

[project.scripts]
mastpipe = "mastpipe.cli:main"

[tool.setuptools.dynamic]
version = {attr = "mastpipe.__version__.__version__"}

[tool.setuptools.packages.find]
where = ["src"]

[tool.setuptools.package-data]
mastpipe = ["instruments/*/*.yaml", "instruments/*/*/*.yaml", "instruments/*/layouts/*.csv", "data/**/*"]

[tool.ruff]
line-length = 100
target-version = "py311"
src = ["src", "tests"]

[tool.ruff.lint]
select = ["E", "F", "I", "UP", "B", "NPY"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-q"
```

- [ ] **Step 2: Package files**

`src/mastpipe/__version__.py`:
```python
__version__ = "0.0.1"
```

`src/mastpipe/__init__.py`:
```python
"""mastpipe: data-reduction pipeline for the MAST spectrographs."""

from mastpipe.__version__ import __version__

__all__ = ["__version__"]
```

`src/mastpipe/cli.py` (the `quicklook` sub-command is wired in Task 11; keep the parser extensible):
```python
"""Command-line interface."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from mastpipe.__version__ import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="mastpipe", description="MAST spectrograph pipeline")
    parser.add_argument("--version", action="version", version=f"mastpipe {__version__}")
    parser.add_subparsers(dest="command")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
```

`README.md` (short: what it is, install with `pip install -e ".[dev]"`, run `pytest`, pointer to `CLAUDE.md`). `CHANGES.md` with a `## 0.0.1 (unreleased)` heading. `DECISIONS.md` with one entry dated 2026-09-24 summarising the spec decision (independent package; soxspipe rejected; PypeIt deferred; docs/design git-ignored for now). `LICENSE`: BSD-3-Clause text with copyright holder "The MAST project". `CLAUDE.md`: copy of the draft with the first italic paragraph removed.

- [ ] **Step 3: Privacy scan script**

`scripts/privacy_scan.py`:
```python
"""Fail if tracked files contain local paths, e-mails, IPs or token-like strings.

Run: python scripts/privacy_scan.py   (exit 1 on any hit)
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

PATTERNS = {
    "local path": re.compile(r"(/mnt[/]|/home[/]|/Users[/]|[A-Za-z]:\\Users\\)"),
    "e-mail": re.compile(r"[\w.+-]+@[\w-]+\.[\w.-]+"),
    "ipv4": re.compile(r"\b(?:\d{1,3}\.){3}\d{1,3}\b"),
    "token": re.compile(r"(ghp_[A-Za-z0-9]{20,}|AKIA[0-9A-Z]{16}|sk-[A-Za-z0-9]{20,})"),
}
ALLOW = {"ipv4": re.compile(r"\b(0\.0\.0\.0|127\.0\.0\.1)\b")}
TEXT_SUFFIXES = {".py", ".md", ".toml", ".yaml", ".yml", ".csv", ".txt", ".cfg", ".ini", ".json"}


def tracked_files() -> list[Path]:
    out = subprocess.run(["git", "ls-files"], capture_output=True, text=True, check=True).stdout
    return [Path(p) for p in out.splitlines() if Path(p).suffix in TEXT_SUFFIXES]


def scan(paths: list[Path]) -> list[str]:
    hits = []
    for path in paths:
        if path.name == "privacy_scan.py":
            continue
        for lineno, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            for label, pat in PATTERNS.items():
                if pat.search(line) and not (label in ALLOW and ALLOW[label].search(line)):
                    hits.append(f"{path}:{lineno}: {label}: {line.strip()[:80]}")
    return hits


def main() -> int:
    hits = scan(tracked_files())
    for h in hits:
        print(h)
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: CI workflow**

`.github/workflows/ci.yml`:
```yaml
name: ci
on: [push, pull_request]
jobs:
  test:
    strategy:
      fail-fast: false
      matrix:
        os: [ubuntu-latest, windows-latest, macos-latest]
        python: ["3.11", "3.12"]
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python }}
      - run: python -m pip install -e ".[dev]"
      - run: ruff check src tests && ruff format --check src tests
      - run: python scripts/privacy_scan.py
      - run: pytest
```

- [ ] **Step 5: Write the failing CLI test**

`tests/conftest.py` (empty for now, fixtures are added in later tasks):
```python
"""Shared fixtures."""
```

`tests/unit/test_cli.py`:
```python
from mastpipe import __version__
from mastpipe.cli import main


def test_version_flag(capsys):
    try:
        main(["--version"])
    except SystemExit as exc:  # argparse exits 0 after printing the version
        assert exc.code == 0
    assert __version__ in capsys.readouterr().out


def test_no_command_prints_help(capsys):
    assert main([]) == 0
    assert "usage: mastpipe" in capsys.readouterr().out
```

- [ ] **Step 6: Install and run**

Run: `python -m pip install -e ".[dev]" && pytest tests/unit/test_cli.py -v`
Expected: 2 passed.

- [ ] **Step 7: Lint, scan, commit**

Run: `ruff check src tests scripts && ruff format src tests scripts && python scripts/privacy_scan.py`
Expected: no output from the scan, exit 0.

```bash
git add -A
git commit -m "chore: package skeleton, CLI, CI matrix, privacy scan"
```

---

### Task 2: Core data model

**Files:**
- Create: `src/mastpipe/core/__init__.py`, `src/mastpipe/core/model.py`
- Test: `tests/unit/test_model.py`

**Interfaces:**
- Produces (used by every later task):
  - `FrameMeta(instrument, camera, exposure_type, configuration=None, exposure_time=None, mid_time=None, binning=(1,1), readout_mode=None, detector_temperature=None, lamp_state=None, fiber_preset=None, target=None, checksum="", source="", extension=0, unresolved=(), classification_rule=None)`; `EXPOSURE_TYPES = ("bias","dark","flat","arc","science","background","unknown")`.
  - `FrameRef(path: Path, extension: int, meta: FrameMeta)`.
  - `Frame(data, meta, variance=None, mask=None)` with `.shape`; `mask` is bool, True = bad.
  - `MasterCalib(kind, data, validity: dict, n_frames, provenance: dict, variance=None)`.
  - `CalibrationSet(bias, dark_rate, background, t_background, gain, read_noise, shape)`.
  - `Trace(index, coeffs, y_center, slope, rms, n_points, snr, x_samples, y_samples, kept, role=None, telescope=None, layout_index=None, aperture=None)` with `.y(x)`.
  - `TraceSet(traces, shape, layout_shift=None, missing=(), extra=())` with `__len__`, `__iter__`, `by_role(role)`.
  - `Spectrum(fiber, pixel, flux, error, mask, role=None, telescope=None, target=None, wavelength=None, sky=None, continuum=None)`.
  - `ExposureSpectra(meta, spectra, wavelength_status="none", provenance=None)`; `WAVELENGTH_STATUSES = ("none","nominal","reference","transferred","solved")`.
  - `QCReport(stage, values: dict)`.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_model.py`:
```python
import numpy as np
import pytest

from mastpipe.core.model import (
    EXPOSURE_TYPES,
    ExposureSpectra,
    Frame,
    FrameMeta,
    Spectrum,
    Trace,
    TraceSet,
)


def make_meta(**kw):
    base = dict(instrument="highspec", camera="main", exposure_type="science", checksum="abc")
    base.update(kw)
    return FrameMeta(**base)


def test_frame_meta_rejects_unknown_exposure_type():
    with pytest.raises(ValueError):
        make_meta(exposure_type="banana")
    assert "science" in EXPOSURE_TYPES


def test_frame_defaults_variance_and_mask():
    f = Frame(data=np.zeros((4, 5)), meta=make_meta())
    assert f.shape == (4, 5)
    assert f.variance is None and f.mask is None


def test_frame_rejects_shape_mismatch():
    with pytest.raises(ValueError):
        Frame(data=np.zeros((4, 5)), meta=make_meta(), mask=np.zeros((5, 4), bool))


def test_trace_evaluates_polynomial():
    t = Trace(index=0, coeffs=np.array([0.001, 0.5, 10.0]), y_center=10.0, slope=0.5,
              rms=0.01, n_points=10, snr=50.0, x_samples=np.arange(3.0),
              y_samples=np.arange(3.0), kept=np.ones(3, bool))
    assert t.y(0.0) == pytest.approx(10.0)
    assert np.allclose(t.y(np.array([0.0, 2.0])), [10.0, 11.004])


def test_traceset_by_role():
    def mk(i, role):
        return Trace(index=i, coeffs=np.array([float(i)]), y_center=float(i), slope=0.0, rms=0.0,
                     n_points=1, snr=1.0, x_samples=np.zeros(1), y_samples=np.zeros(1),
                     kept=np.ones(1, bool), role=role)
    ts = TraceSet(traces=[mk(0, "arc"), mk(1, "science"), mk(2, "sky")], shape=(400, 1600))
    assert len(ts) == 3
    assert [t.index for t in ts.by_role("science")] == [1]


def test_spectrum_and_exposure_spectra_validate():
    s = Spectrum(fiber=0, pixel=np.arange(5), flux=np.ones(5), error=np.ones(5),
                 mask=np.zeros(5, bool))
    es = ExposureSpectra(meta=make_meta(), spectra=[s])
    assert es.wavelength_status == "none"
    with pytest.raises(ValueError):
        ExposureSpectra(meta=make_meta(), spectra=[s], wavelength_status="guessed")
    with pytest.raises(ValueError):
        Spectrum(fiber=0, pixel=np.arange(5), flux=np.ones(4), error=np.ones(5),
                 mask=np.zeros(5, bool))
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/unit/test_model.py -v` — Expected: ImportError / failures.

- [ ] **Step 3: Implement**

`src/mastpipe/core/__init__.py`: empty docstring module.

`src/mastpipe/core/model.py`:
```python
"""Core dataclasses shared by all layers. All arrays are numpy; masks are bool with True = bad."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

EXPOSURE_TYPES = ("bias", "dark", "flat", "arc", "science", "background", "unknown")
WAVELENGTH_STATUSES = ("none", "nominal", "reference", "transferred", "solved")
ROLES = ("arc", "science", "sky")


@dataclass(frozen=True)
class FrameMeta:
    """Pipeline-vocabulary description of one raw frame (one camera, one extension)."""

    instrument: str
    camera: str
    exposure_type: str
    configuration: str | None = None
    exposure_time: float | None = None
    mid_time: str | None = None  # ISO-8601 UTC
    binning: tuple[int, int] = (1, 1)
    readout_mode: str | None = None
    detector_temperature: float | None = None
    lamp_state: str | None = None
    fiber_preset: str | None = None
    target: str | None = None
    checksum: str = ""
    source: str = ""  # file name only, never a directory
    extension: int = 0
    unresolved: tuple[str, ...] = ()
    classification_rule: str | None = None

    def __post_init__(self) -> None:
        if self.exposure_type not in EXPOSURE_TYPES:
            raise ValueError(f"exposure_type {self.exposure_type!r} not in {EXPOSURE_TYPES}")


@dataclass(frozen=True)
class FrameRef:
    """Where a frame lives plus its resolved metadata; pixel data are loaded on demand."""

    path: Path
    extension: int
    meta: FrameMeta


def _check_same_shape(name: str, arr: np.ndarray | None, shape: tuple[int, ...]) -> None:
    if arr is not None and arr.shape != shape:
        raise ValueError(f"{name} shape {arr.shape} != data shape {shape}")


@dataclass
class Frame:
    data: np.ndarray
    meta: FrameMeta
    variance: np.ndarray | None = None
    mask: np.ndarray | None = None

    def __post_init__(self) -> None:
        self.data = np.asarray(self.data)
        _check_same_shape("variance", self.variance, self.data.shape)
        _check_same_shape("mask", self.mask, self.data.shape)

    @property
    def shape(self) -> tuple[int, ...]:
        return self.data.shape


@dataclass
class MasterCalib:
    kind: str  # "bias" | "dark_rate" | "background" | "flat"
    data: np.ndarray
    validity: dict[str, Any]
    n_frames: int
    provenance: dict[str, Any] = field(default_factory=dict)
    variance: np.ndarray | None = None


@dataclass
class CalibrationSet:
    """Master frames applied to a science frame. dark_rate is ADU/s (array or scalar)."""

    bias: np.ndarray | None
    dark_rate: np.ndarray | float
    background: np.ndarray | None
    t_background: float | None
    gain: float  # e-/ADU
    read_noise: float  # e-
    shape: tuple[int, int]


@dataclass
class Trace:
    index: int
    coeffs: np.ndarray  # numpy polyval order (highest power first), y(x) in pixels
    y_center: float
    slope: float
    rms: float
    n_points: int
    snr: float
    x_samples: np.ndarray
    y_samples: np.ndarray
    kept: np.ndarray
    role: str | None = None
    telescope: str | None = None
    layout_index: int | None = None
    aperture: float | None = None

    def y(self, x: float | np.ndarray) -> np.ndarray:
        return np.polyval(self.coeffs, np.asarray(x, dtype=float))


@dataclass
class TraceSet:
    traces: list[Trace]
    shape: tuple[int, int]
    layout_shift: float | None = None
    missing: tuple[int, ...] = ()  # layout indices with no detected trace
    extra: tuple[int, ...] = ()  # detected trace indices with no layout row

    def __len__(self) -> int:
        return len(self.traces)

    def __iter__(self) -> Iterator[Trace]:
        return iter(self.traces)

    def by_role(self, role: str) -> list[Trace]:
        return [t for t in self.traces if t.role == role]


@dataclass
class Spectrum:
    fiber: int
    pixel: np.ndarray
    flux: np.ndarray
    error: np.ndarray
    mask: np.ndarray
    role: str | None = None
    telescope: str | None = None
    target: str | None = None
    wavelength: np.ndarray | None = None
    sky: np.ndarray | None = None
    continuum: np.ndarray | None = None

    def __post_init__(self) -> None:
        n = self.pixel.shape
        for name in ("flux", "error", "mask", "wavelength", "sky", "continuum"):
            _check_same_shape(name, getattr(self, name), n)


@dataclass
class ExposureSpectra:
    meta: FrameMeta
    spectra: list[Spectrum]
    wavelength_status: str = "none"
    provenance: dict[str, Any] | None = None

    def __post_init__(self) -> None:
        if self.wavelength_status not in WAVELENGTH_STATUSES:
            raise ValueError(
                f"wavelength_status {self.wavelength_status!r} not in {WAVELENGTH_STATUSES}"
            )


@dataclass
class QCReport:
    stage: str
    values: dict[str, Any]
```

- [ ] **Step 4: Run tests** — `pytest tests/unit/test_model.py -v` — Expected: 6 passed.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat(core): data model dataclasses"`

---

### Task 3: Configuration models and layered loader

**Files:**
- Create: `src/mastpipe/config/__init__.py`, `src/mastpipe/config/models.py`, `src/mastpipe/config/loader.py`
- Test: `tests/unit/test_config.py`

**Interfaces:**
- Produces:
  - Pydantic models `DetectorConfig`, `GeometryConfig`, `CalibParams`, `TraceParams`, `ExtractParams`, `StepParams`, `ConfigurationBlock`, `InstrumentConfig`, `EffectiveConfig`.
  - `deep_merge(base: dict, override: dict) -> dict` (override wins, dicts recurse, lists replace).
  - `load_yaml(path: Path) -> dict`.
  - `build_effective(instrument: InstrumentConfig, configuration: str, camera: str, overrides: list[dict] = ()) -> EffectiveConfig`.
  - `config_hash(effective: EffectiveConfig) -> str` (12 hex chars of sha256 over canonical JSON).

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_config.py`:
```python
import pytest
from pydantic import ValidationError

from mastpipe.config.loader import build_effective, config_hash, deep_merge
from mastpipe.config.models import (
    ConfigurationBlock,
    DetectorConfig,
    GeometryConfig,
    InstrumentConfig,
    StepParams,
)


def make_instrument():
    return InstrumentConfig(
        name="toy",
        cameras={"main": DetectorConfig(name="main", shape=(400, 1600), gain=1.0,
                                        read_noise=5.0, pixel_size_um=16.0)},
        params=StepParams(),
        configurations={
            "g518": ConfigurationBlock(
                name="g518",
                geometry=GeometryConfig(groove_spacing_nm=367.0, incidence_angle_deg=45.0,
                                        order=1, camera_focal_mm=479.41,
                                        central_wavelength_nm=518.0),
                band_nm=(511.0, 525.0),
                params={"trace": {"trace_halfwidth": 4}},
            )
        },
    )


def test_deep_merge_recurses_and_lists_replace():
    base = {"a": {"x": 1, "y": [1, 2]}, "b": 2}
    out = deep_merge(base, {"a": {"y": [9]}, "c": 3})
    assert out == {"a": {"x": 1, "y": [9]}, "b": 2, "c": 3}
    assert base["a"]["y"] == [1, 2]  # not mutated


def test_effective_layers_in_order():
    inst = make_instrument()
    eff = build_effective(inst, "g518", "main")
    assert eff.params.trace.trace_halfwidth == 4  # configuration block wins over defaults
    eff2 = build_effective(inst, "g518", "main", overrides=[{"trace": {"trace_halfwidth": 6}}])
    assert eff2.params.trace.trace_halfwidth == 6
    assert eff2.detector.shape == (400, 1600)
    assert eff2.geometry.central_wavelength_nm == 518.0


def test_effective_unknown_configuration_or_camera():
    inst = make_instrument()
    with pytest.raises(KeyError):
        build_effective(inst, "nope", "main")
    with pytest.raises(KeyError):
        build_effective(inst, "g518", "nope")


def test_config_hash_is_stable_and_sensitive():
    inst = make_instrument()
    a = config_hash(build_effective(inst, "g518", "main"))
    b = config_hash(build_effective(inst, "g518", "main"))
    c = config_hash(build_effective(inst, "g518", "main", overrides=[{"trace": {"x_bin": 10}}]))
    assert a == b and a != c and len(a) == 12


def test_validation_rejects_bad_values():
    with pytest.raises(ValidationError):
        DetectorConfig(name="m", shape=(400, 1600), gain=-1.0, read_noise=5.0, pixel_size_um=16)
    with pytest.raises(ValidationError):
        StepParams(extract={"aperture": "triangle"})
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_config.py -v`

- [ ] **Step 3: Implement**

`src/mastpipe/config/__init__.py`: docstring only.

`src/mastpipe/config/models.py`:
```python
"""Pydantic configuration models. Instruments are data, not code."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, PositiveFloat, PositiveInt


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class DetectorConfig(_Strict):
    name: str
    shape: tuple[PositiveInt, PositiveInt]  # (ny, nx) unbinned
    gain: PositiveFloat  # e-/ADU
    read_noise: float = Field(ge=0)  # e-
    pixel_size_um: PositiveFloat
    saturation: float | None = None  # ADU
    crop: tuple[int, int, int, int] | None = None  # x_min, x_max, y_min, y_max (exclusive max)
    dispersion_axis: Literal[0, 1] = 1  # 1 = along columns (x)


class GeometryConfig(_Strict):
    groove_spacing_nm: PositiveFloat
    incidence_angle_deg: float
    order: PositiveInt
    camera_focal_mm: PositiveFloat
    central_wavelength_nm: PositiveFloat
    dispersion_sign: Literal[-1, 1] = 1
    magnification: PositiveFloat = 1.6  # fiber image / fiber diameter
    fiber_diameter_um: PositiveFloat = 30.0


class CalibParams(_Strict):
    subtract_background: bool = True
    flat_subtract_background: bool = False
    stack_method: Literal["median", "mean"] = "median"
    clip_sigma: float | None = None  # sigma clipping in mean stacks
    temperature_tolerance: float = 2.0  # K, for calibration-library matching


class TraceParams(_Strict):
    center_width: PositiveInt = 100
    trace_halfwidth: PositiveInt = 8
    threshold_sigma: PositiveFloat = 5.0
    centroid_sigma: float | None = None
    min_separation: float | None = None
    baseline_window: int | None = None
    x_bin: PositiveInt = 25
    degree: int = Field(default=2, ge=0, le=4)
    clip_sigma: float | None = 3.0
    min_points: PositiveInt = 4
    max_traces: int | None = None
    layout_tolerance: PositiveFloat = 3.0  # px, matching detections to the layout table


class ExtractParams(_Strict):
    halfwidth: PositiveFloat = 9.0
    aperture: Literal["round", "subpixel", "floor"] = "round"
    method: Literal["sum", "mean"] = "sum"
    min_coverage: float = Field(default=0.99, gt=0, le=1)
    continuum_degree: int = Field(default=15, ge=0)
    clip_above: PositiveFloat = 2.0
    clip_below: PositiveFloat = 3.0


class StepParams(_Strict):
    calib: CalibParams = CalibParams()
    trace: TraceParams = TraceParams()
    extract: ExtractParams = ExtractParams()


class ConfigurationBlock(_Strict):
    """One grating (HighSpec) or band (DeepSpec)."""

    name: str
    band_nm: tuple[PositiveFloat, PositiveFloat]
    geometry: GeometryConfig | None = None
    atlas: str | None = None  # package data name, Plan 2
    linelist: str | None = None  # package data name, Plan 2
    resolving_power: float | None = None
    params: dict[str, Any] = Field(default_factory=dict)  # partial StepParams overrides


class InstrumentConfig(_Strict):
    name: str
    cameras: dict[str, DetectorConfig]
    params: StepParams = StepParams()
    configurations: dict[str, ConfigurationBlock]
    layouts: dict[str, str] = Field(default_factory=dict)  # preset -> csv file name
    default_layout: str | None = None


class EffectiveConfig(_Strict):
    """Everything a run needs for one instrument, camera and configuration."""

    instrument: str
    camera: str
    configuration: str
    detector: DetectorConfig
    geometry: GeometryConfig | None
    band_nm: tuple[float, float]
    params: StepParams
    layout: str | None = None
```

`src/mastpipe/config/loader.py`:
```python
"""YAML loading, deep merge, layering and hashing of configurations."""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import yaml

from mastpipe.config.models import EffectiveConfig, InstrumentConfig, StepParams


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Return a new dict: override wins; nested dicts merge; lists and scalars replace."""
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def load_yaml(path: Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as fh:
        data = yaml.safe_load(fh) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: top level must be a mapping")
    return data


def build_effective(
    instrument: InstrumentConfig,
    configuration: str,
    camera: str,
    overrides: Iterable[dict[str, Any]] = (),
    layout: str | None = None,
) -> EffectiveConfig:
    """Layer: instrument params -> configuration block params -> overrides (in order)."""
    if configuration not in instrument.configurations:
        raise KeyError(f"unknown configuration {configuration!r} for {instrument.name}")
    if camera not in instrument.cameras:
        raise KeyError(f"unknown camera {camera!r} for {instrument.name}")
    block = instrument.configurations[configuration]
    merged = instrument.params.model_dump()
    merged = deep_merge(merged, block.params)
    for ov in overrides:
        merged = deep_merge(merged, ov)
    return EffectiveConfig(
        instrument=instrument.name,
        camera=camera,
        configuration=configuration,
        detector=instrument.cameras[camera],
        geometry=block.geometry,
        band_nm=block.band_nm,
        params=StepParams.model_validate(merged),
        layout=layout or instrument.default_layout,
    )


def config_hash(effective: EffectiveConfig) -> str:
    payload = json.dumps(effective.model_dump(mode="json"), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:12]
```

- [ ] **Step 4: Run tests** — `pytest tests/unit/test_config.py -v` — Expected: 5 passed.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat(config): pydantic models, deep merge, layering, config hash"`

---

### Task 4: Instrument registry and the HighSpec definition

**Files:**
- Create: `src/mastpipe/instruments/__init__.py`, `src/mastpipe/instruments/registry.py`, `src/mastpipe/instruments/highspec/config.yaml`, `src/mastpipe/instruments/highspec/header_map.yaml`, `src/mastpipe/instruments/highspec/gratings/mgb518.yaml`, `src/mastpipe/instruments/highspec/layouts/lab3.csv`, `src/mastpipe/instruments/highspec/layouts/survey45.csv`, `scripts/make_survey45_layout.py`
- Test: `tests/unit/test_registry.py`

**Interfaces:**
- Produces:
  - `instrument_dir(name) -> Path` (package data folder).
  - `load_instrument(name) -> InstrumentConfig` (reads `config.yaml`, merges every `gratings/*.yaml` or `bands/*.yaml` into `configurations`). HighSpec will get more gratings over time: each is one YAML file in `gratings/`, discovered by glob, so adding a grating is a data change with no code edit. The test below proves that with a temporary extra file.
  - `load_header_map(name) -> dict` (raw YAML; parsed by Task 5's `HeaderMap`).
  - `load_layout(name, preset) -> pandas.DataFrame` with columns `trace:int, role:str, telescope:str|None, y_nominal:float, slope_nominal:float, aperture:float|None`, sorted by `y_nominal`, validated (roles in ROLES, unique trace ids).
  - `list_instruments() -> list[str]`.

Numbers below come from the HighSpec laboratory reduction notebook (grating 2724.8 l/mm → 367 nm spacing; α = 45°; f_cam = 479.41 mm; 16 µm pixels; magnification 4.8/3; 30 µm fibers). Gain and read noise are **placeholders to be measured**; say so in the YAML comments.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_registry.py`:
```python
import pytest

from mastpipe.core.model import ROLES
from mastpipe.instruments.registry import (
    list_instruments,
    load_header_map,
    load_instrument,
    load_layout,
)


def test_highspec_is_registered():
    assert "highspec" in list_instruments()


def test_highspec_config_loads_with_grating_blocks():
    inst = load_instrument("highspec")
    assert inst.name == "highspec"
    assert "main" in inst.cameras
    assert inst.cameras["main"].shape == (400, 1600)
    assert "mgb518" in inst.configurations
    g = inst.configurations["mgb518"].geometry
    assert g is not None and g.central_wavelength_nm == pytest.approx(518.0)
    assert inst.default_layout == "lab3"


def test_layout_tables_parse_and_validate():
    lab3 = load_layout("highspec", "lab3")
    assert list(lab3.columns) == ["trace", "role", "telescope", "y_nominal", "slope_nominal",
                                  "aperture"]
    assert len(lab3) == 3 and set(lab3.role) <= set(ROLES)
    survey = load_layout("highspec", "survey45")
    assert len(survey) == 45
    assert (survey.role == "arc").sum() == 5
    assert (survey.role == "science").sum() == 20 and (survey.role == "sky").sum() == 20
    assert survey.y_nominal.is_monotonic_increasing
    assert survey.loc[survey.role != "arc", "telescope"].notna().all()


def test_header_map_yaml_has_required_sections():
    hm = load_header_map("highspec")
    assert "fields" in hm and "classification" in hm
    assert "exposure_type" in hm["fields"] and "exposure_time" in hm["fields"]


def test_unknown_instrument():
    with pytest.raises(KeyError):
        load_instrument("nope")


def test_extra_grating_file_is_picked_up(tmp_path, monkeypatch):
    """A new HighSpec grating is one YAML file; no code change may be needed."""
    import shutil

    from mastpipe.instruments import registry

    copy = tmp_path / "highspec"
    shutil.copytree(registry.instrument_dir("highspec"), copy)
    (copy / "gratings" / "halpha.yaml").write_text(
        "name: halpha\nband_nm: [649.0, 664.0]\nresolving_power: 20000\n"
        "geometry: {groove_spacing_nm: 367.0, incidence_angle_deg: 45.0, order: 1, "
        "camera_focal_mm: 479.41, central_wavelength_nm: 656.3}\nparams: {}\n",
        encoding="utf-8")
    monkeypatch.setattr(registry, "instrument_dir", lambda name: copy)
    inst = registry.load_instrument("highspec")
    assert set(inst.configurations) >= {"mgb518", "halpha"}
    assert inst.configurations["halpha"].geometry.central_wavelength_nm == pytest.approx(656.3)
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_registry.py -v`

- [ ] **Step 3: Instrument data files**

`src/mastpipe/instruments/highspec/config.yaml`:
```yaml
name: highspec
cameras:
  main:
    name: main
    shape: [400, 1600]        # ny, nx (Andor Newton DU971P-BV EMCCD)
    gain: 1.0                 # e-/ADU  PLACEHOLDER: measure from flat pairs
    read_noise: 5.0           # e-      PLACEHOLDER: measure from bias pairs
    pixel_size_um: 16.0
    saturation: 65000
    crop: [1, 1599, 1, 399]   # x_min, x_max, y_min, y_max (exclusive max); drops the edge row/column
    dispersion_axis: 1
params:
  calib:
    subtract_background: true
    flat_subtract_background: false
  trace:
    center_width: 100
    trace_halfwidth: 8        # ~2.5 x the 3-px fiber image
    threshold_sigma: 5.0
    x_bin: 25
    degree: 2                 # quadratic absorbs the conical-diffraction curvature
    clip_sigma: 3.0
    min_points: 4
    layout_tolerance: 3.0
  extract:
    halfwidth: 9.0            # 3 x the fiber image
    aperture: round
    method: sum
    min_coverage: 0.99
    continuum_degree: 15
    clip_above: 2.0
    clip_below: 3.0
layouts:
  lab3: lab3.csv
  survey45: survey45.csv
default_layout: lab3
```

`src/mastpipe/instruments/highspec/gratings/mgb518.yaml`:
```yaml
name: mgb518
band_nm: [511.0, 526.0]
resolving_power: 20000
geometry:
  groove_spacing_nm: 367.0        # 2724.8 lines/mm
  incidence_angle_deg: 45.0       # held fixed: degenerate with the plate scale
  order: 1
  camera_focal_mm: 479.41
  central_wavelength_nm: 518.0
  dispersion_sign: 1
  magnification: 1.6              # 4.8 / 3
  fiber_diameter_um: 30.0
params: {}
```

`src/mastpipe/instruments/highspec/header_map.yaml` (provisional keywords: Andor SOLIS names that exist today plus the names the pipeline will look for once the control software writes them; edit here and nowhere else when they change):
```yaml
fields:
  exposure_time:        {keys: [EXPOSURE, EXPTIME, T_EXP], type: float}
  mid_time:             {keys: [DATE-OBS, T_MID, DATE], type: str}
  detector_temperature: {keys: [TEMP, CCD-TEMP, TEMP_GOAL], type: float}
  readout_mode:         {keys: [READMODE, READOUT_SPEED], type: str}
  binning_x:            {keys: [HBIN, CDELT1], type: int, default: 1}
  binning_y:            {keys: [VBIN, CDELT2], type: int, default: 1}
  exposure_type:
    keys: [IMAGETYP, TYPE]
    type: str
    lower: true
    map: {bias: bias, zero: bias, dark: dark, flat: flat, arc: arc, comp: arc, object: science,
          light: science, exp: science, science: science, background: background, bkg: background}
  configuration:        {keys: [GRATING, CONFIG], type: str, lower: true,
                         filename_regex: "(?i)(?P<value>g\\d{3})", map: {g518: mgb518}}
  lamp_state:           {keys: [LAMP], type: str, lower: true}
  fiber_preset:         {keys: [FIBPRE, FIBERPRE], type: str}
  target:               {keys: [OBJECT, TARGET], type: str}
  camera:               {keys: [BAND, CAMERA], type: str, default: main}
classification:
  - {when: {exposure_time: {eq: 0}}, type: bias}
  - {when: {lamp_state: {in: [thar, th-ar, arc]}}, type: arc}
  - {when: {lamp_state: {in: [flat, halogen, led, continuum]}}, type: flat}
  - {when: {lamp_state: {in: [off, none]}, exposure_time: {gt: 0}}, type: science}
```

`src/mastpipe/instruments/highspec/layouts/lab3.csv` (the notebook's three-fiber bench setup, uncropped detector coordinates; slopes from the measured tilts +0.100°, +0.117°, +0.133°):
```csv
trace,role,telescope,y_nominal,slope_nominal,aperture
0,arc,,18.4,0.00175,
1,arc,,190.0,0.00204,
2,arc,,361.7,0.00232,
```

`scripts/make_survey45_layout.py` writes `survey45.csv` (45 traces evenly spaced from y = 12 to y = 388, arc fibers at slots 0, 11, 22, 33, 44, the remaining 40 slots alternating science/sky for telescopes T01…T20); run it once and commit the CSV:
```python
"""Generate the provisional 45-trace survey layout (placeholder geometry until measured)."""

import csv
from pathlib import Path

import numpy as np

out = Path(__file__).resolve().parents[1] / "src/mastpipe/instruments/highspec/layouts/survey45.csv"
y = np.linspace(12.0, 388.0, 45)
arc_slots = {0, 11, 22, 33, 44}
rows, tel = [], 1
for i, yy in enumerate(y):
    if i in arc_slots:
        rows.append((i, "arc", "", f"{yy:.1f}", "0.0020", ""))
    else:
        role = "science" if (i - sum(s < i for s in arc_slots)) % 2 == 0 else "sky"
        rows.append((i, role, f"T{tel:02d}", f"{yy:.1f}", "0.0020", ""))
        if role == "sky":
            tel += 1
with out.open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["trace", "role", "telescope", "y_nominal", "slope_nominal", "aperture"])
    w.writerows(rows)
print(out)
```

- [ ] **Step 4: Registry implementation**

`src/mastpipe/instruments/__init__.py`: docstring only.

`src/mastpipe/instruments/registry.py`:
```python
"""Locate and load instrument definitions shipped as package data."""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import pandas as pd

from mastpipe.config.loader import load_yaml
from mastpipe.config.models import ConfigurationBlock, InstrumentConfig
from mastpipe.core.model import ROLES

LAYOUT_COLUMNS = ["trace", "role", "telescope", "y_nominal", "slope_nominal", "aperture"]


def _root() -> Path:
    return Path(str(resources.files("mastpipe.instruments")))


def list_instruments() -> list[str]:
    return sorted(p.name for p in _root().iterdir() if (p / "config.yaml").is_file())


def instrument_dir(name: str) -> Path:
    path = _root() / name
    if not (path / "config.yaml").is_file():
        raise KeyError(f"unknown instrument {name!r}; known: {list_instruments()}")
    return path


def load_instrument(name: str) -> InstrumentConfig:
    folder = instrument_dir(name)
    raw = load_yaml(folder / "config.yaml")
    blocks = dict(raw.get("configurations", {}))
    for sub in ("gratings", "bands"):
        for path in sorted((folder / sub).glob("*.yaml")):
            block = ConfigurationBlock.model_validate(load_yaml(path))
            blocks[block.name] = block.model_dump()
    raw["configurations"] = blocks
    return InstrumentConfig.model_validate(raw)


def load_header_map(name: str) -> dict:
    return load_yaml(instrument_dir(name) / "header_map.yaml")


def load_layout(name: str, preset: str) -> pd.DataFrame:
    inst = load_instrument(name)
    if preset not in inst.layouts:
        raise KeyError(f"unknown layout preset {preset!r} for {name}; known: {list(inst.layouts)}")
    path = instrument_dir(name) / "layouts" / inst.layouts[preset]
    df = pd.read_csv(path, dtype={"telescope": "string"})
    missing = [c for c in LAYOUT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing layout columns {missing}")
    df = df[LAYOUT_COLUMNS].copy()
    df["telescope"] = df["telescope"].astype(object).where(df["telescope"].notna(), None)
    bad_roles = set(df.role) - set(ROLES)
    if bad_roles:
        raise ValueError(f"{path.name}: unknown roles {sorted(bad_roles)}")
    if df.trace.duplicated().any():
        raise ValueError(f"{path.name}: duplicate trace ids")
    return df.sort_values("y_nominal", kind="stable").reset_index(drop=True)
```

- [ ] **Step 5: Generate the survey layout, run tests**

Run: `python scripts/make_survey45_layout.py && pytest tests/unit/test_registry.py -v` — Expected: 6 passed.

- [ ] **Step 6: Commit** — `git add -A && git commit -m "feat(instruments): registry + HighSpec definition (config, header map, layouts, Mg b grating)"`

---

### Task 5: Header map resolution and ingest

**Files:**
- Create: `src/mastpipe/io/__init__.py`, `src/mastpipe/io/headers.py`, `src/mastpipe/io/ingest.py`
- Test: `tests/unit/test_headers.py`, `tests/unit/test_ingest.py`

**Interfaces:**
- Consumes: `FrameMeta`, `FrameRef`, `Frame` (Task 2); `DetectorConfig` (Task 3); `load_header_map` (Task 4).
- Produces:
  - `HeaderMap.from_dict(d) -> HeaderMap`; `HeaderMap.resolve(header: Mapping, filename: str) -> tuple[dict, list[str]]` (resolved fields, unresolved names). Handles `keys` priority, `type` cast, `lower`, `map`, `default`, `filename_regex` (named group `value`).
  - `HeaderMap.classify(fields: dict) -> tuple[str, str]` → `(exposure_type, rule_name)`; `rule_name` is `"header"` when the header value was valid, else `"rule:<index>"`, else `"unclassified"` with type `unknown`.
  - `scan_folder(folder: Path, instrument: str, header_map: HeaderMap, pattern: str = "*.fit*", manifest: pd.DataFrame | None = None, unreadable: list[Path] | None = None) -> list[FrameRef]` — one `FrameRef` per image HDU (primary if it has data, else every `ImageHDU`); MEF extension headers are resolved after merging the primary header underneath them. A file that `astropy` cannot open raises `OSError`, unless an `unreadable` list is passed, in which case the path is appended to it and scanning continues (the quick-look runner uses this).
  - `load_frame(ref: FrameRef, detector: DetectorConfig) -> Frame` — float32 data, crop applied, no variance yet.
  - `read_manifest(path) -> pd.DataFrame` (columns: `source`, optional `extension`, plus any FrameMeta field to override).
  - `file_checksum(path: Path) -> str` (sha256 hex, first 16 chars).

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_headers.py`:
```python
import pytest

from mastpipe.instruments.registry import load_header_map
from mastpipe.io.headers import HeaderMap


@pytest.fixture
def hm():
    return HeaderMap.from_dict(load_header_map("highspec"))


def test_resolve_andor_like_header(hm):
    header = {"EXPOSURE": 6.0, "TEMP": -50.0, "HBIN": 1, "VBIN": 1, "READMODE": "Image",
              "DATE": "2026-09-09T10:00:00"}
    fields, unresolved = hm.resolve(header, "FOCUS_05_13.30004_mm.fits")
    assert fields["exposure_time"] == 6.0
    assert fields["detector_temperature"] == -50.0
    assert fields["binning_x"] == 1 and fields["camera"] == "main"
    assert "exposure_type" in unresolved and "configuration" in unresolved


def test_map_and_lower_and_filename_regex(hm):
    fields, _ = hm.resolve({"IMAGETYP": "OBJECT", "EXPOSURE": 300.0}, "g518_thar_001.fits")
    assert fields["exposure_type"] == "science"
    assert fields["configuration"] == "mgb518"


def test_classify_prefers_header_then_rules(hm):
    assert hm.classify({"exposure_type": "flat", "exposure_time": 0.0}) == ("flat", "header")
    assert hm.classify({"exposure_time": 0.0}) == ("bias", "rule:0")
    assert hm.classify({"exposure_time": 6.0, "lamp_state": "thar"}) == ("arc", "rule:1")
    assert hm.classify({"exposure_time": 6.0}) == ("unknown", "unclassified")


def test_bad_type_cast_is_unresolved(hm):
    fields, unresolved = hm.resolve({"EXPOSURE": "six"}, "x.fits")
    assert "exposure_time" in unresolved and "exposure_time" not in fields
```

`tests/unit/test_ingest.py`:
```python
import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.config.models import DetectorConfig
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.headers import HeaderMap
from mastpipe.io.ingest import load_frame, read_manifest, scan_folder


@pytest.fixture
def hm():
    return HeaderMap.from_dict(load_header_map("highspec"))


def write_simple(path, data, **cards):
    hdu = fits.PrimaryHDU(data.astype(np.uint16))
    for k, v in cards.items():
        hdu.header[k] = v
    hdu.writeto(path)


def test_scan_single_extension_frames(tmp_path, hm):
    write_simple(tmp_path / "bias1.fits", np.full((400, 1600), 500), EXPOSURE=0.0, TEMP=-50.0)
    write_simple(tmp_path / "g518_arc.fits", np.full((400, 1600), 600), EXPOSURE=6.0,
                 LAMP="ThAr", TEMP=-50.0)
    refs = scan_folder(tmp_path, "highspec", hm)
    by_name = {r.meta.source: r for r in refs}
    assert by_name["bias1.fits"].meta.exposure_type == "bias"
    arc = by_name["g518_arc.fits"].meta
    assert arc.exposure_type == "arc" and arc.configuration == "mgb518"
    assert arc.lamp_state == "thar" and arc.classification_rule == "rule:1"
    assert arc.checksum and arc.source == "g518_arc.fits"  # never a directory


def test_scan_multi_extension_splits_cameras(tmp_path, hm):
    primary = fits.PrimaryHDU()
    primary.header["EXPOSURE"] = 10.0
    primary.header["IMAGETYP"] = "science"
    hdus = [primary]
    for band in ("U", "G"):
        h = fits.ImageHDU(np.zeros((10, 12), dtype=np.float32))
        h.header["BAND"] = band
        hdus.append(h)
    fits.HDUList(hdus).writeto(tmp_path / "mef.fits")
    refs = scan_folder(tmp_path, "highspec", hm)
    assert sorted((r.extension, r.meta.camera) for r in refs) == [(1, "U"), (2, "G")]
    assert all(r.meta.exposure_time == 10.0 for r in refs)


def test_manifest_overrides_fields(tmp_path, hm):
    write_simple(tmp_path / "lab.fits", np.zeros((400, 1600)), EXPOSURE=6.0)
    manifest = pd.DataFrame([{"source": "lab.fits", "exposure_type": "arc",
                              "configuration": "mgb518", "fiber_preset": "lab3"}])
    (ref,) = scan_folder(tmp_path, "highspec", hm, manifest=manifest)
    assert ref.meta.exposure_type == "arc" and ref.meta.classification_rule == "manifest"
    assert ref.meta.fiber_preset == "lab3"


def test_load_frame_crops_and_casts(tmp_path, hm):
    write_simple(tmp_path / "f.fits", np.arange(400 * 1600).reshape(400, 1600) % 60000,
                 EXPOSURE=1.0, IMAGETYP="flat")
    (ref,) = scan_folder(tmp_path, "highspec", hm)
    det = DetectorConfig(name="main", shape=(400, 1600), gain=1.0, read_noise=5.0,
                         pixel_size_um=16.0, crop=(1, 1599, 1, 399))
    frame = load_frame(ref, det)
    assert frame.shape == (398, 1598) and frame.data.dtype == np.float32
    assert frame.meta is ref.meta


def test_read_manifest_requires_source_column(tmp_path):
    p = tmp_path / "m.csv"
    p.write_text("file,exposure_type\na.fits,arc\n", encoding="utf-8")
    with pytest.raises(ValueError):
        read_manifest(p)


def test_scan_collects_unreadable_files(tmp_path, hm):
    (tmp_path / "junk.fits").write_bytes(b"nope")
    bad = []
    assert scan_folder(tmp_path, "highspec", hm, unreadable=bad) == []
    assert bad == [tmp_path / "junk.fits"]
    with pytest.raises(OSError):
        scan_folder(tmp_path, "highspec", hm)
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_headers.py tests/unit/test_ingest.py -v`

- [ ] **Step 3: Implement `headers.py`**

```python
"""Translate raw FITS headers into the pipeline vocabulary. The only place raw keywords live."""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from mastpipe.core.model import EXPOSURE_TYPES

_CASTS = {"float": float, "int": int, "str": str, "bool": bool}


@dataclass(frozen=True)
class FieldSpec:
    name: str
    keys: tuple[str, ...]
    type: str = "str"
    lower: bool = False
    map: dict[str, Any] = field(default_factory=dict)
    default: Any = None
    filename_regex: str | None = None

    def cast(self, value: Any) -> Any:
        if isinstance(value, str):
            value = value.strip()
        value = _CASTS[self.type](value)
        if self.lower and isinstance(value, str):
            value = value.lower()
        if self.map:
            key = value.lower() if isinstance(value, str) else value
            value = self.map.get(key, value)
        return value


@dataclass(frozen=True)
class Rule:
    when: dict[str, dict[str, Any]]
    type: str

    def matches(self, fields: Mapping[str, Any]) -> bool:
        for name, cond in self.when.items():
            value = fields.get(name)
            for op, ref in cond.items():
                if op == "eq" and not (value is not None and value == ref):
                    return False
                if op == "gt" and not (value is not None and value > ref):
                    return False
                if op == "lt" and not (value is not None and value < ref):
                    return False
                if op == "in" and value not in ref:
                    return False
                if op == "isnull" and (value is None) != bool(ref):
                    return False
        return True


@dataclass(frozen=True)
class HeaderMap:
    fields: tuple[FieldSpec, ...]
    rules: tuple[Rule, ...]

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> HeaderMap:
        specs = []
        for name, spec in d["fields"].items():
            spec = dict(spec)
            specs.append(
                FieldSpec(
                    name=name,
                    keys=tuple(spec.pop("keys", ())),
                    type=spec.pop("type", "str"),
                    lower=spec.pop("lower", False),
                    map={str(k).lower(): v for k, v in spec.pop("map", {}).items()},
                    default=spec.pop("default", None),
                    filename_regex=spec.pop("filename_regex", None),
                )
            )
            if spec:
                raise ValueError(f"header map field {name}: unknown options {sorted(spec)}")
        rules = []
        for r in d.get("classification", []):
            if r["type"] not in EXPOSURE_TYPES:
                raise ValueError(f"classification rule type {r['type']!r} unknown")
            rules.append(Rule(when=dict(r["when"]), type=r["type"]))
        return cls(fields=tuple(specs), rules=tuple(rules))

    def resolve(self, header: Mapping[str, Any], filename: str) -> tuple[dict[str, Any], list[str]]:
        fields: dict[str, Any] = {}
        unresolved: list[str] = []
        for spec in self.fields:
            value = None
            for key in spec.keys:
                if key in header and header[key] not in (None, ""):
                    try:
                        value = spec.cast(header[key])
                    except (TypeError, ValueError):
                        value = None
                    break
            if value is None and spec.filename_regex:
                m = re.search(spec.filename_regex, filename)
                if m:
                    try:
                        value = spec.cast(m.group("value"))
                    except (TypeError, ValueError):
                        value = None
            if value is None and spec.default is not None:
                value = spec.cast(spec.default)
            if value is None:
                unresolved.append(spec.name)
            else:
                fields[spec.name] = value
        return fields, unresolved

    def classify(self, fields: Mapping[str, Any]) -> tuple[str, str]:
        et = fields.get("exposure_type")
        if et in EXPOSURE_TYPES and et != "unknown":
            return et, "header"
        for i, rule in enumerate(self.rules):
            if rule.matches(fields):
                return rule.type, f"rule:{i}"
        return "unknown", "unclassified"
```

- [ ] **Step 4: Implement `ingest.py`**

```python
"""Scan a folder of raw frames into FrameRefs; load pixel data on demand."""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.config.models import DetectorConfig
from mastpipe.core.model import Frame, FrameMeta, FrameRef
from mastpipe.io.headers import HeaderMap

META_FIELDS = {
    "configuration", "exposure_time", "mid_time", "readout_mode", "detector_temperature",
    "lamp_state", "fiber_preset", "target", "camera", "exposure_type",
}


def file_checksum(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def read_manifest(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    if "source" not in df.columns:
        raise ValueError(f"{path}: manifest needs a 'source' column")
    if "extension" not in df.columns:
        df["extension"] = ""
    return df


def _manifest_overrides(manifest: pd.DataFrame | None, source: str, ext: int) -> dict[str, Any]:
    if manifest is None:
        return {}
    if "extension" not in manifest.columns:
        manifest = manifest.assign(extension="")
    manifest = manifest.fillna("")
    rows = manifest[manifest["source"] == source]
    rows = rows[(rows["extension"] == "") | (rows["extension"].astype(str) == str(ext))]
    if rows.empty:
        return {}
    row = rows.iloc[0]
    out = {}
    for col in rows.columns:
        if col in ("source", "extension") or row[col] == "":
            continue
        if col not in META_FIELDS:
            raise ValueError(f"manifest column {col!r} is not a FrameMeta field")
        out[col] = float(row[col]) if col in ("exposure_time", "detector_temperature") else row[col]
    return out


def _to_meta(instrument: str, fields: Mapping[str, Any], unresolved: list[str], rule: str,
             checksum: str, source: str, ext: int) -> FrameMeta:
    return FrameMeta(
        instrument=instrument,
        camera=str(fields.get("camera", "main")),
        exposure_type=str(fields.get("exposure_type", "unknown")),
        configuration=fields.get("configuration"),
        exposure_time=fields.get("exposure_time"),
        mid_time=fields.get("mid_time"),
        binning=(int(fields.get("binning_x", 1)), int(fields.get("binning_y", 1))),
        readout_mode=fields.get("readout_mode"),
        detector_temperature=fields.get("detector_temperature"),
        lamp_state=fields.get("lamp_state"),
        fiber_preset=fields.get("fiber_preset"),
        target=fields.get("target"),
        checksum=checksum,
        source=source,
        extension=ext,
        unresolved=tuple(unresolved),
        classification_rule=rule,
    )


def _image_hdus(hdul: fits.HDUList) -> list[tuple[int, fits.Header]]:
    out = []
    for i, hdu in enumerate(hdul):
        is_image = isinstance(hdu, fits.PrimaryHDU | fits.ImageHDU)
        if is_image and hdu.header.get("NAXIS", 0) >= 2:
            header = fits.Header(hdul[0].header)
            if i > 0:
                header.update(hdu.header)  # extension cards win over primary
            out.append((i, header))
    return out


def scan_folder(folder: Path, instrument: str, header_map: HeaderMap, pattern: str = "*.fit*",
                manifest: pd.DataFrame | None = None,
                unreadable: list[Path] | None = None) -> list[FrameRef]:
    refs: list[FrameRef] = []
    for path in sorted(Path(folder).glob(pattern)):
        checksum = file_checksum(path)
        try:
            hdul = fits.open(path, memmap=False)
        except OSError:
            if unreadable is None:
                raise
            unreadable.append(path)
            continue
        with hdul:
            for ext, header in _image_hdus(hdul):
                fields, unresolved = header_map.resolve(header, path.name)
                overrides = _manifest_overrides(manifest, path.name, ext)
                fields.update(overrides)
                unresolved = [u for u in unresolved if u not in overrides]
                if "exposure_type" in overrides:
                    et, rule = overrides["exposure_type"], "manifest"
                else:
                    et, rule = header_map.classify(fields)
                fields["exposure_type"] = et
                meta = _to_meta(instrument, fields, unresolved, rule,
                                f"{checksum}:{ext}", path.name, ext)
                refs.append(FrameRef(path=path, extension=ext, meta=meta))
    return refs


def load_frame(ref: FrameRef, detector: DetectorConfig) -> Frame:
    with fits.open(ref.path, memmap=False) as hdul:
        data = np.asarray(hdul[ref.extension].data, dtype=np.float32)
    if detector.crop is not None:
        x0, x1, y0, y1 = detector.crop
        data = data[y0:y1, x0:x1]
    return Frame(data=data, meta=ref.meta)
```

- [ ] **Step 5: Run tests** — `pytest tests/unit/test_headers.py tests/unit/test_ingest.py -v` — Expected: 10 passed.

- [ ] **Step 6: Commit** — `git add -A && git commit -m "feat(io): header map resolution, classification rules, folder scan, manifest, frame loading"`

---

### Task 6: Grating geometry and the synthetic-night simulator

**Files:**
- Create: `src/mastpipe/steps/__init__.py`, `src/mastpipe/steps/geometry.py`, `src/mastpipe/sim/__init__.py`, `src/mastpipe/sim/simulator.py`
- Test: `tests/unit/test_geometry.py`, `tests/unit/test_simulator.py`; add fixtures to `tests/conftest.py`

**Interfaces:**
- Consumes: `GeometryConfig`, `DetectorConfig`, `EffectiveConfig` (Task 3); layout DataFrame (Task 4).
- Produces:
  - `diffraction_angle(geom, wavelength_nm) -> float` (radians, β from mλ = d(sin α + sin β)).
  - `nominal_dispersion_nm_per_px(geom, pixel_size_um, wavelength_nm=None) -> float`.
  - `nominal_wavelength(geom, pixel_size_um, x, x_center) -> np.ndarray` (exact grating equation, not linear).
  - `SimTruth` dataclass: `layout` (DataFrame), `trace_coeffs` (list of arrays, numpy polyval order, cropped-frame coordinates), `dispersion_nm_per_px`, `x_center`, `lines_nm` (array), `line_pixels` (list per fiber), `throughput` (per fiber), `science_flux` (per fiber, ADU per pixel row summed), `sky_flux`, `bias_level`, `dark_rate`, `read_noise_e`, `gain`, `fwhm_px`, `crop`, `extras` (dict with `drift_px`, `dips_nm`, `sky_lines_nm`, `line_amp` = arc line amplitudes, same order as `lines_nm`).
  - `simulate_night(effective: EffectiveConfig, layout: pd.DataFrame, out_dir: Path, seed: int = 0, n_bias=5, n_dark=3, n_flat=3, n_arc=2, n_science=2, exposure_science=300.0, drift_px=0.0, cosmic_rays=20, fiber_preset="lab3") -> SimTruth` — writes FITS frames named `sim_<type>_<nn>.fits` with header cards `IMAGETYP, EXPOSURE, DATE-OBS, TEMP, HBIN, VBIN, READMODE, GRATING, LAMP, OBJECT, FIBPRE` matching the HighSpec header map; saves `truth.npz` + `truth.json`; frames are **uncropped** (`detector.shape`) so the crop path is exercised.
  - `load_truth(out_dir) -> SimTruth`.
- Simulation model: bias level 500 ADU; read noise `detector.read_noise` e; gain `detector.gain`; dark 0.05 ADU/s; trace `i`: `y(x) = y_nominal_i + slope_i (x - xc) + c_i (x - xc)^2` with `c_i` drawn in ±2e-6 (conical curvature); spatial profile Gaussian with `fwhm_px = magnification * fiber_diameter_um / pixel_size_um` (3.0 px for HighSpec); flat: smooth quadratic blaze × per-fiber throughput (0.8–1.2) × 20 000 ADU peak; arc: `n_lines=40` lines at seeded wavelengths in the band with log-uniform intensities, line FWHM in x = `fwhm_px`, placed at `x = x_center + (λ - λc)/D` (linear approximation of the nominal solution is fine for truth; record `line_pixels` exactly as placed); science fibers: continuum 200 ADU/px with 6 Gaussian absorption dips of depth 0.3 and FWHM 0.05 nm; sky fibers: continuum 20 ADU/px plus 3 emission lines; `drift_px` shifts all lines and spectra in x for science frames (Plan 2 will use this); cosmic rays: random single pixels + 5000 ADU in science frames. Poisson noise on everything, Gaussian read noise, `uint16` output clipped at 65535.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_geometry.py`:
```python
import numpy as np
import pytest

from mastpipe.config.models import GeometryConfig
from mastpipe.steps.geometry import (
    diffraction_angle,
    nominal_dispersion_nm_per_px,
    nominal_wavelength,
)

GEOM = GeometryConfig(groove_spacing_nm=367.0, incidence_angle_deg=45.0, order=1,
                      camera_focal_mm=479.41, central_wavelength_nm=518.0)


def test_diffraction_angle_matches_grating_equation():
    beta = diffraction_angle(GEOM, 518.0)
    lhs = GEOM.order * 518.0
    rhs = GEOM.groove_spacing_nm * (np.sin(np.radians(45.0)) + np.sin(beta))
    assert lhs == pytest.approx(rhs)


def test_nominal_dispersion_matches_notebook():
    # The notebook's fitted solution runs 0.00887 -> 0.00834 nm/px across the detector.
    d = nominal_dispersion_nm_per_px(GEOM, 16.0)
    assert 0.0083 < d < 0.0090


def test_nominal_wavelength_is_monotonic_and_centred():
    x = np.arange(1598.0)
    lam = nominal_wavelength(GEOM, 16.0, x, x_center=799.0)
    assert lam[799] == pytest.approx(518.0, abs=1e-9)
    assert np.all(np.diff(lam) > 0)
    assert lam[-1] - lam[0] == pytest.approx(1598 * 0.0086, rel=0.05)
```

`tests/unit/test_simulator.py`:
```python
import numpy as np
from astropy.io import fits

from mastpipe.sim.simulator import load_truth, simulate_night


def test_simulate_night_writes_frames_and_truth(sim_night):
    out_dir, truth = sim_night
    names = sorted(p.name for p in out_dir.glob("*.fits"))
    assert sum(n.startswith("sim_bias") for n in names) == 5
    assert sum(n.startswith("sim_flat") for n in names) == 3
    assert sum(n.startswith("sim_arc") for n in names) == 2
    assert sum(n.startswith("sim_science") for n in names) == 2
    h = fits.getheader(out_dir / "sim_arc_00.fits")
    assert h["IMAGETYP"] == "arc" and h["LAMP"] == "thar" and h["GRATING"] == "g518"
    assert fits.getdata(out_dir / "sim_arc_00.fits").shape == (400, 1600)
    assert len(truth.trace_coeffs) == len(truth.layout) == 3
    assert load_truth(out_dir).lines_nm.shape == truth.lines_nm.shape


def test_simulator_is_deterministic(tmp_path, highspec_effective, lab3_layout):
    a = simulate_night(highspec_effective, lab3_layout, tmp_path / "a", seed=3, n_bias=1,
                       n_dark=1, n_flat=1, n_arc=1, n_science=1)
    b = simulate_night(highspec_effective, lab3_layout, tmp_path / "b", seed=3, n_bias=1,
                       n_dark=1, n_flat=1, n_arc=1, n_science=1)
    assert np.array_equal(fits.getdata(tmp_path / "a/sim_arc_00.fits"),
                          fits.getdata(tmp_path / "b/sim_arc_00.fits"))
    assert np.allclose(a.lines_nm, b.lines_nm)


def test_flat_peaks_sit_on_truth_traces(sim_night):
    out_dir, truth = sim_night
    flat = fits.getdata(out_dir / "sim_flat_00.fits").astype(float)
    x0, x1, y0, y1 = 1, 1599, 1, 399
    flat = flat[y0:y1, x0:x1]
    col = np.median(flat[:, 700:900], axis=1)
    xc = 799.0
    for coeffs in truth.trace_coeffs:
        y_true = np.polyval(coeffs, xc)
        lo, hi = int(y_true) - 6, int(y_true) + 7
        assert abs(lo + np.argmax(col[lo:hi]) - y_true) <= 1.0
```

`tests/conftest.py` additions:
```python
import pytest

from mastpipe.config.loader import build_effective
from mastpipe.instruments.registry import load_instrument, load_layout


@pytest.fixture(scope="session")
def highspec_instrument():
    return load_instrument("highspec")


@pytest.fixture(scope="session")
def highspec_effective(highspec_instrument):
    return build_effective(highspec_instrument, "mgb518", "main", layout="lab3")


@pytest.fixture(scope="session")
def lab3_layout():
    return load_layout("highspec", "lab3")


@pytest.fixture(scope="session")
def sim_night(tmp_path_factory, highspec_effective, lab3_layout):
    from mastpipe.sim.simulator import simulate_night

    out = tmp_path_factory.mktemp("night")
    truth = simulate_night(highspec_effective, lab3_layout, out, seed=42)
    return out, truth
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_geometry.py tests/unit/test_simulator.py -v`

- [ ] **Step 3: Implement `geometry.py`**

```python
"""Grating-equation geometry: the given constants, never fitted here."""

from __future__ import annotations

import numpy as np

from mastpipe.config.models import GeometryConfig


def diffraction_angle(geom: GeometryConfig, wavelength_nm: float | np.ndarray) -> np.ndarray:
    """beta such that m*lambda = d*(sin alpha + sin beta)."""
    alpha = np.radians(geom.incidence_angle_deg)
    s = geom.order * np.asarray(wavelength_nm, dtype=float) / geom.groove_spacing_nm - np.sin(alpha)
    if np.any(np.abs(s) > 1):
        raise ValueError("wavelength outside the grating's physical range")
    return np.arcsin(s)


def nominal_dispersion_nm_per_px(geom: GeometryConfig, pixel_size_um: float,
                                 wavelength_nm: float | None = None) -> float:
    lam = geom.central_wavelength_nm if wavelength_nm is None else wavelength_nm
    beta = float(diffraction_angle(geom, lam))
    pixel_mm = pixel_size_um * 1e-3
    return geom.groove_spacing_nm * np.cos(beta) * pixel_mm / (geom.order * geom.camera_focal_mm)


def nominal_wavelength(geom: GeometryConfig, pixel_size_um: float, x: np.ndarray,
                       x_center: float) -> np.ndarray:
    """Exact grating equation: beta(x) = beta_c + atan((x - xc) * pixel / f)."""
    beta_c = float(diffraction_angle(geom, geom.central_wavelength_nm))
    dx_mm = (np.asarray(x, dtype=float) - x_center) * pixel_size_um * 1e-3 * geom.dispersion_sign
    beta = beta_c + np.arctan(dx_mm / geom.camera_focal_mm)
    alpha = np.radians(geom.incidence_angle_deg)
    return geom.groove_spacing_nm * (np.sin(alpha) + np.sin(beta)) / geom.order
```

- [ ] **Step 4: Implement `simulator.py`**

```python
"""Synthetic nights with known truth. Frames are written uncropped in the detector's shape."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.config.models import EffectiveConfig
from mastpipe.steps.geometry import nominal_dispersion_nm_per_px

BIAS_LEVEL = 500.0
DARK_RATE = 0.05  # ADU/s
FLAT_PEAK = 20000.0
SCI_CONTINUUM = 200.0
SKY_CONTINUUM = 20.0
ARC_PEAK = 30000.0
FWHM_TO_SIGMA = 1.0 / 2.3548


@dataclass
class SimTruth:
    layout: pd.DataFrame
    trace_coeffs: list[np.ndarray]
    dispersion_nm_per_px: float
    x_center: float
    lines_nm: np.ndarray
    line_pixels: list[np.ndarray]
    throughput: np.ndarray
    science_flux: list[np.ndarray]
    sky_flux: list[np.ndarray]
    bias_level: float
    dark_rate: float
    read_noise_e: float
    gain: float
    fwhm_px: float
    crop: tuple[int, int, int, int] | None
    extras: dict = field(default_factory=dict)


def _gauss(y: np.ndarray, mu: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((y[:, None] - mu[None, :]) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))


def _render(shape, coeffs_list, fluxes, sigma_y):
    """Sum of traces: each trace has a per-column flux (ADU) spread over a Gaussian in y."""
    ny, nx = shape
    x = np.arange(nx, dtype=float)
    y = np.arange(ny, dtype=float)
    img = np.zeros(shape, dtype=float)
    for coeffs, flux in zip(coeffs_list, fluxes, strict=True):
        img += _gauss(y, np.polyval(coeffs, x), sigma_y) * flux[None, :]
    return img


def _write(path: Path, image: np.ndarray, rng, gain, read_noise_e, **cards) -> None:
    electrons = rng.poisson(np.clip(image, 0, None) * gain)
    adu = electrons / gain + rng.normal(0.0, read_noise_e / gain, image.shape) + BIAS_LEVEL
    data = np.clip(np.round(adu), 0, 65535).astype(np.uint16)
    hdu = fits.PrimaryHDU(data)
    for k, v in cards.items():
        hdu.header[k] = v
    hdu.writeto(path, overwrite=True)


def simulate_night(effective: EffectiveConfig, layout: pd.DataFrame, out_dir: Path, seed: int = 0,
                   n_bias: int = 5, n_dark: int = 3, n_flat: int = 3, n_arc: int = 2,
                   n_science: int = 2, exposure_science: float = 300.0, drift_px: float = 0.0,
                   cosmic_rays: int = 20, fiber_preset: str = "lab3") -> SimTruth:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    det, geom = effective.detector, effective.geometry
    if geom is None:
        raise ValueError("simulate_night needs a configuration with geometry")
    ny, nx = det.shape
    crop = det.crop
    x0 = crop[0] if crop else 0
    y0 = crop[2] if crop else 0
    nx_c = (crop[1] - crop[0]) if crop else nx
    x_center = (nx_c - 1) / 2.0
    fwhm_px = geom.magnification * geom.fiber_diameter_um / det.pixel_size_um
    sigma = fwhm_px * FWHM_TO_SIGMA
    disp = nominal_dispersion_nm_per_px(geom, det.pixel_size_um)
    x_full = np.arange(nx, dtype=float)
    x_crop = x_full - x0  # cropped-frame column coordinate of every full-frame column

    # traces in cropped coordinates (what the pipeline sees), rendered in full coordinates
    n = len(layout)
    curv = rng.uniform(-2e-6, 2e-6, n)
    coeffs_crop = [np.array([curv[i], layout.slope_nominal[i] - 2 * curv[i] * x_center,
                             layout.y_nominal[i] - y0 - layout.slope_nominal[i] * x_center
                             + curv[i] * x_center**2]) for i in range(n)]
    coeffs_full = [np.polyfit(x_full, np.polyval(c, x_crop) + y0, 2) for c in coeffs_crop]
    throughput = rng.uniform(0.8, 1.2, n)

    lo, hi = effective.band_nm
    lam_c = geom.central_wavelength_nm
    lines_nm = np.sort(rng.uniform(lo + 0.5, hi - 0.5, 40))
    line_amp = 10 ** rng.uniform(-1.5, 0, lines_nm.size) * ARC_PEAK
    line_px_crop = x_center + (lines_nm - lam_c) / disp  # cropped coordinate

    def arc_flux(shift: float) -> np.ndarray:
        prof = np.zeros(nx)
        for px, amp in zip(line_px_crop + shift + x0, line_amp, strict=True):
            prof += amp * np.exp(-0.5 * ((x_full - px) / sigma) ** 2)
        return prof

    blaze = 1.0 - 0.3 * ((x_full - nx / 2) / (nx / 2)) ** 2
    dips_nm = rng.uniform(lo + 1, hi - 1, 6)
    lam_full = lam_c + (x_crop - x_center) * disp
    sci = SCI_CONTINUUM * blaze
    for d in dips_nm:
        sci *= 1 - 0.3 * np.exp(-0.5 * ((lam_full - d) / (0.05 * FWHM_TO_SIGMA)) ** 2)
    sky_lines_nm = rng.uniform(lo + 1, hi - 1, 3)
    sky = SKY_CONTINUUM * np.ones(nx)
    for s in sky_lines_nm:
        sky += 200 * np.exp(-0.5 * ((lam_full - s) / (sigma * disp)) ** 2)

    roles = layout.role.to_numpy()
    sci_fluxes = [sci * throughput[i] if roles[i] == "science" else np.zeros(nx) for i in range(n)]
    sky_fluxes = [sky * throughput[i] if roles[i] != "arc" else np.zeros(nx) for i in range(n)]

    common = dict(TEMP=-50.0, HBIN=1, VBIN=1, READMODE="Image", GRATING="g518",
                  FIBPRE=fiber_preset)
    t0 = pd.Timestamp("2026-10-01T18:00:00")

    def stamp(i):
        return (t0 + pd.Timedelta(minutes=5 * i)).strftime("%Y-%m-%dT%H:%M:%S")

    k = 0
    for i in range(n_bias):
        _write(out_dir / f"sim_bias_{i:02d}.fits", np.zeros((ny, nx)), rng, det.gain,
               det.read_noise, IMAGETYP="bias", EXPOSURE=0.0, LAMP="off", OBJECT="bias",
               **{"DATE-OBS": stamp(k)}, **common)
        k += 1
    for i in range(n_dark):
        _write(out_dir / f"sim_dark_{i:02d}.fits", np.full((ny, nx), DARK_RATE * exposure_science),
               rng, det.gain, det.read_noise, IMAGETYP="dark", EXPOSURE=exposure_science,
               LAMP="off", OBJECT="dark", **{"DATE-OBS": stamp(k)}, **common)
        k += 1
    flat_fluxes = [FLAT_PEAK * blaze * throughput[i] for i in range(n)]
    for i in range(n_flat):
        img = _render((ny, nx), coeffs_full, flat_fluxes, sigma) + DARK_RATE * 5.0
        _write(out_dir / f"sim_flat_{i:02d}.fits", img, rng, det.gain, det.read_noise,
               IMAGETYP="flat", EXPOSURE=5.0, LAMP="flat", OBJECT="flat",
               **{"DATE-OBS": stamp(k)}, **common)
        k += 1
    arc_all = [arc_flux(0.0) * throughput[i] for i in range(n)]
    for i in range(n_arc):
        img = _render((ny, nx), coeffs_full, arc_all, sigma) + DARK_RATE * 6.0
        _write(out_dir / f"sim_arc_{i:02d}.fits", img, rng, det.gain, det.read_noise,
               IMAGETYP="arc", EXPOSURE=6.0, LAMP="thar", OBJECT="thar",
               **{"DATE-OBS": stamp(k)}, **common)
        k += 1
    arc_only = [arc_flux(drift_px) * throughput[i] if roles[i] == "arc" else np.zeros(nx)
                for i in range(n)]
    for i in range(n_science):
        fluxes = [np.roll(sci_fluxes[j] + sky_fluxes[j], int(round(drift_px))) + arc_only[j]
                  for j in range(n)]
        img = _render((ny, nx), coeffs_full, fluxes, sigma) + DARK_RATE * exposure_science
        for _ in range(cosmic_rays):
            img[rng.integers(0, ny), rng.integers(0, nx)] += 5000.0
        _write(out_dir / f"sim_science_{i:02d}.fits", img, rng, det.gain, det.read_noise,
               IMAGETYP="science", EXPOSURE=exposure_science, LAMP="off",
               OBJECT=f"SIMSTAR{i}", **{"DATE-OBS": stamp(k)}, **common)
        k += 1

    truth = SimTruth(layout=layout.copy(), trace_coeffs=coeffs_crop, dispersion_nm_per_px=disp,
                     x_center=x_center, lines_nm=lines_nm,
                     line_pixels=[line_px_crop.copy() for _ in range(n)], throughput=throughput,
                     science_flux=[f[x0:x0 + nx_c] for f in sci_fluxes],
                     sky_flux=[f[x0:x0 + nx_c] for f in sky_fluxes], bias_level=BIAS_LEVEL,
                     dark_rate=DARK_RATE, read_noise_e=det.read_noise, gain=det.gain,
                     fwhm_px=fwhm_px, crop=crop,
                     extras={"drift_px": drift_px, "dips_nm": dips_nm.tolist(),
                             "sky_lines_nm": sky_lines_nm.tolist(),
                             "line_amp": line_amp.tolist()})
    _save_truth(out_dir, truth)
    return truth


def _save_truth(out_dir: Path, t: SimTruth) -> None:
    np.savez(out_dir / "truth.npz", trace_coeffs=np.array(t.trace_coeffs), lines_nm=t.lines_nm,
             line_pixels=np.array(t.line_pixels), throughput=t.throughput,
             science_flux=np.array(t.science_flux), sky_flux=np.array(t.sky_flux))
    scalars = {k: v for k, v in asdict(t).items()
               if k not in ("layout", "trace_coeffs", "lines_nm", "line_pixels", "throughput",
                            "science_flux", "sky_flux")}
    scalars["crop"] = list(t.crop) if t.crop else None
    (out_dir / "truth.json").write_text(json.dumps(scalars, indent=1), encoding="utf-8")
    t.layout.to_csv(out_dir / "truth_layout.csv", index=False)


def load_truth(out_dir: Path) -> SimTruth:
    out_dir = Path(out_dir)
    arr = np.load(out_dir / "truth.npz")
    scalars = json.loads((out_dir / "truth.json").read_text(encoding="utf-8"))
    crop = tuple(scalars.pop("crop")) if scalars.get("crop") else None
    scalars.pop("crop", None)
    layout = pd.read_csv(out_dir / "truth_layout.csv", dtype={"telescope": "string"})
    return SimTruth(layout=layout, trace_coeffs=list(arr["trace_coeffs"]), lines_nm=arr["lines_nm"],
                    line_pixels=list(arr["line_pixels"]), throughput=arr["throughput"],
                    science_flux=list(arr["science_flux"]), sky_flux=list(arr["sky_flux"]),
                    crop=crop, **scalars)
```

- [ ] **Step 5: Run tests** — `pytest tests/unit/test_geometry.py tests/unit/test_simulator.py -v` — Expected: 6 passed. If the flat-peak test fails by more than a pixel, check the full-frame/cropped coordinate bookkeeping (`x0`, `y0`) before touching tolerances.

- [ ] **Step 6: Commit** — `git add -A && git commit -m "feat(sim): grating geometry helpers and synthetic-night simulator with truth"`

---

### Task 7: Detector calibration step (port of the notebook's calibration cell)

**Files:**
- Create: `src/mastpipe/steps/calib.py`
- Test: `tests/unit/test_calib.py`

**Interfaces:**
- Consumes: `Frame`, `MasterCalib`, `CalibrationSet` (Task 2); `CalibParams`, `DetectorConfig` (Task 3).
- Produces:
  - `stack_frames(frames: list[np.ndarray], method="median", clip_sigma=None) -> tuple[np.ndarray, np.ndarray]` → (master, per-pixel variance of the master).
  - `build_master(frames: list[Frame], kind: str, params: CalibParams) -> MasterCalib` — validity = camera, binning, readout_mode, mean detector temperature, mean mid_time, exposure_time (for dark/background); provenance lists input checksums.
  - `build_calibration_set(bias: MasterCalib | None, dark: MasterCalib | None, background: MasterCalib | None, detector: DetectorConfig) -> CalibrationSet` — `dark_rate = (dark - bias)/t_dark` when both exist and `t_dark > 0`, else `0.0`; `t_background` from the background master's validity.
  - `calibrate_frame(frame: Frame, calib: CalibrationSet, subtract_background: bool = True, saturation: float | None = None) -> Frame` — output float32 data, `variance` = `max(data_e,0)/gain² + RN²/gain²` in ADU² computed from the bias-subtracted signal, `mask` = saturated pixels (in the raw frame or in any master) and non-finite values; the same arithmetic as the notebook: `out = raw − bias − (background − bias)·t/t_bkg` when subtracting background, else `out = raw − bias − dark_rate·t`.
  - `calibration_steps(calib, t_frame, subtract_background) -> list[str]` (human-readable, for provenance).

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_calib.py`:
```python
import numpy as np
import pytest

from mastpipe.config.models import CalibParams, DetectorConfig
from mastpipe.core.model import Frame, FrameMeta
from mastpipe.steps.calib import (
    build_calibration_set,
    build_master,
    calibrate_frame,
    calibration_steps,
    stack_frames,
)

DET = DetectorConfig(name="main", shape=(4, 6), gain=2.0, read_noise=4.0, pixel_size_um=16.0)


def frame(value, kind, t, temp=-50.0, cs="c"):
    meta = FrameMeta(instrument="toy", camera="main", exposure_type=kind, exposure_time=t,
                     detector_temperature=temp, checksum=cs, mid_time="2026-10-01T18:00:00")
    return Frame(data=np.full((4, 6), value, dtype=np.float32), meta=meta)


def test_stack_median_rejects_outlier_and_reports_variance():
    stack = [np.full((2, 2), 10.0) for _ in range(5)]
    stack[2][0, 0] = 1000.0
    master, var = stack_frames(stack, method="median")
    assert master[0, 0] == 10.0 and var.shape == (2, 2)


def test_build_master_collects_validity_and_provenance():
    m = build_master([frame(500, "bias", 0.0, cs="a"), frame(502, "bias", 0.0, cs="b")], "bias",
                     CalibParams())
    assert m.kind == "bias" and m.n_frames == 2
    assert m.validity["camera"] == "main" and m.validity["exposure_time"] == 0.0
    assert m.provenance["inputs"] == ["a", "b"]
    assert m.data[0, 0] == pytest.approx(501.0)


def test_calibration_set_dark_rate():
    bias = build_master([frame(500, "bias", 0.0)], "bias", CalibParams())
    dark = build_master([frame(530, "dark", 300.0)], "dark", CalibParams())
    cs = build_calibration_set(bias, dark, None, DET)
    assert np.allclose(cs.dark_rate, 0.1)
    assert cs.background is None and cs.t_background is None
    cs0 = build_calibration_set(bias, None, None, DET)
    assert cs0.dark_rate == 0.0


def test_calibrate_frame_subtracts_and_propagates_variance():
    bias = build_master([frame(500, "bias", 0.0)], "bias", CalibParams())
    dark = build_master([frame(530, "dark", 300.0)], "dark", CalibParams())
    cs = build_calibration_set(bias, dark, None, DET)
    sci = frame(1500, "science", 100.0)
    out = calibrate_frame(sci, cs, subtract_background=False)
    assert out.data.dtype == np.float32
    assert np.allclose(out.data, 1500 - 500 - 0.1 * 100)  # 990 ADU
    expected_var = 990.0 / DET.gain + (DET.read_noise / DET.gain) ** 2  # ADU^2
    assert np.allclose(out.variance, expected_var, rtol=1e-3)
    assert not out.mask.any()
    assert out.meta is sci.meta


def test_calibrate_frame_background_scaling_and_saturation_mask():
    bias = build_master([frame(500, "bias", 0.0)], "bias", CalibParams())
    bkg = build_master([frame(560, "background", 6.0)], "background", CalibParams())
    cs = build_calibration_set(bias, None, bkg, DET)
    sci = frame(1000, "science", 12.0)
    sci.data[1, 1] = 70000.0
    out = calibrate_frame(sci, cs, subtract_background=True, saturation=65000)
    assert out.data[0, 0] == pytest.approx(1000 - 500 - 60 * 2)
    assert out.mask[1, 1] and not out.mask[0, 0]
    assert "background" in " ".join(calibration_steps(cs, 12.0, True))


def test_calibrate_requires_background_when_asked():
    bias = build_master([frame(500, "bias", 0.0)], "bias", CalibParams())
    cs = build_calibration_set(bias, None, None, DET)
    with pytest.raises(ValueError):
        calibrate_frame(frame(1000, "science", 1.0), cs, subtract_background=True)
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_calib.py -v`

- [ ] **Step 3: Implement**

```python
"""Master frames and detector calibration. Port of the notebook's calibration cell.

Frame model: raw = bias + dark_rate*t + background_ambient*t/t_bkg + signal.
"""

from __future__ import annotations

import numpy as np

from mastpipe.config.models import CalibParams, DetectorConfig
from mastpipe.core.model import CalibrationSet, Frame, MasterCalib


def stack_frames(frames: list[np.ndarray], method: str = "median",
                 clip_sigma: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    stack = np.stack([np.asarray(f, dtype=np.float32) for f in frames])
    n = stack.shape[0]
    if method == "median":
        master = np.nanmedian(stack, axis=0)
        scatter = 1.4826 * np.nanmedian(np.abs(stack - master), axis=0) if n > 2 else np.nanstd(
            stack, axis=0)
        var = (1.2533 * scatter) ** 2 / max(n, 1)  # median variance ~ pi/2 * sigma^2 / n
    elif method == "mean":
        good = np.isfinite(stack)
        if clip_sigma:
            med = np.nanmedian(stack, axis=0)
            sig = 1.4826 * np.nanmedian(np.abs(stack - med), axis=0)
            good &= np.abs(stack - med) <= clip_sigma * np.maximum(sig, 1e-6)
        cnt = good.sum(axis=0)
        master = np.where(good, stack, 0).sum(axis=0) / np.maximum(cnt, 1)
        var = np.nanvar(np.where(good, stack, np.nan), axis=0) / np.maximum(cnt, 1)
    else:
        raise ValueError(f"unknown stack method {method!r}")
    return master.astype(np.float32), np.nan_to_num(var).astype(np.float32)


def _mean_or_none(values):
    vals = [v for v in values if v is not None]
    return float(np.mean(vals)) if vals else None


def build_master(frames: list[Frame], kind: str, params: CalibParams) -> MasterCalib:
    if not frames:
        raise ValueError(f"no frames to build master {kind}")
    metas = [f.meta for f in frames]
    if len({m.camera for m in metas}) != 1 or len({m.binning for m in metas}) != 1:
        raise ValueError(f"master {kind}: frames mix cameras or binnings")
    data, var = stack_frames([f.data for f in frames], params.stack_method, params.clip_sigma)
    validity = {
        "camera": metas[0].camera,
        "binning": metas[0].binning,
        "readout_mode": metas[0].readout_mode,
        "detector_temperature": _mean_or_none(m.detector_temperature for m in metas),
        "exposure_time": _mean_or_none(m.exposure_time for m in metas),
        "mid_time": sorted(m.mid_time for m in metas if m.mid_time)[len(metas) // 2]
        if any(m.mid_time for m in metas) else None,
        "configuration": metas[0].configuration,
    }
    return MasterCalib(kind=kind, data=data, validity=validity, n_frames=len(frames),
                       provenance={"inputs": [m.checksum for m in metas],
                                   "method": params.stack_method}, variance=var)


def build_calibration_set(bias: MasterCalib | None, dark: MasterCalib | None,
                          background: MasterCalib | None,
                          detector: DetectorConfig) -> CalibrationSet:
    ref = next((m for m in (bias, dark, background) if m is not None), None)
    if ref is None:
        raise ValueError("a calibration set needs at least one master frame")
    shape = tuple(ref.data.shape)
    dark_rate: np.ndarray | float = 0.0
    if dark is not None and bias is not None:
        t_dark = dark.validity.get("exposure_time") or 0.0
        if t_dark > 0:
            dark_rate = ((dark.data - bias.data) / float(t_dark)).astype(np.float32)
    return CalibrationSet(
        bias=None if bias is None else bias.data,
        dark_rate=dark_rate,
        background=None if background is None else background.data,
        t_background=None if background is None else background.validity.get("exposure_time"),
        gain=detector.gain,
        read_noise=detector.read_noise,
        shape=shape,  # type: ignore[arg-type]
    )


def calibration_steps(calib: CalibrationSet, t_frame: float, subtract_background: bool) -> list[str]:
    steps = []
    if calib.bias is not None:
        steps.append("bias")
    if subtract_background:
        scale = float(t_frame) / float(calib.t_background)
        steps.append("background" if np.isclose(scale, 1.0) else f"background x {scale:.4g}")
    elif np.any(np.asarray(calib.dark_rate) * float(t_frame)):
        steps.append(f"dark_rate x {t_frame:g} s")
    return steps


def calibrate_frame(frame: Frame, calib: CalibrationSet, subtract_background: bool = True,
                    saturation: float | None = None) -> Frame:
    t = frame.meta.exposure_time
    if t is None:
        raise ValueError(f"{frame.meta.source}: exposure time unknown, cannot calibrate")
    if subtract_background and (calib.background is None or not calib.t_background):
        raise ValueError(f"{frame.meta.source}: background subtraction requested but no "
                         "background master with an exposure time is available")
    raw = np.asarray(frame.data, dtype=np.float32)
    if raw.shape != tuple(calib.shape):
        raise ValueError(f"{frame.meta.source}: shape {raw.shape} != calibration {calib.shape}")
    out = raw.copy()
    if calib.bias is not None:
        out -= calib.bias
    if subtract_background:
        bkg = calib.background if calib.bias is None else calib.background - calib.bias
        out -= bkg * (float(t) / float(calib.t_background))
    else:
        out -= np.asarray(calib.dark_rate, dtype=np.float32) * float(t)
    signal_e = np.clip(out, 0, None) * calib.gain
    variance = (signal_e / calib.gain**2 + (calib.read_noise / calib.gain) ** 2).astype(np.float32)
    mask = ~np.isfinite(out)
    if saturation is not None:
        mask |= raw >= saturation
        for extra in (calib.bias, calib.background):
            if extra is not None:
                mask |= np.asarray(extra) >= saturation
    if frame.mask is not None:
        mask |= frame.mask
    return Frame(data=out, meta=frame.meta, variance=variance, mask=mask)
```

- [ ] **Step 4: Run tests** — `pytest tests/unit/test_calib.py -v` — Expected: 6 passed.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat(steps): master frames, calibration set, frame calibration with variance and masks"`

---

### Task 8: Trace finding and layout matching (port of `find_traces`)

**Files:**
- Create: `src/mastpipe/steps/trace.py`
- Test: `tests/unit/test_trace.py`

**Interfaces:**
- Consumes: `Frame`, `Trace`, `TraceSet` (Task 2); `TraceParams` (Task 3); layout DataFrame (Task 4); simulator fixtures (Task 6); `calibrate_frame` (Task 7).
- Produces:
  - `robust_sigma(v) -> float`.
  - `find_traces(image: np.ndarray, params: TraceParams) -> TraceSet` — the notebook algorithm verbatim (detection band → baseline median filter → `find_peaks` → per-block centroids → clipped polynomial fit, two passes), traces sorted by `y_center`, `index` assigned 0..n-1, roles unset.
  - `match_layout(traceset: TraceSet, layout: pd.DataFrame, tolerance: float) -> TraceSet` — estimates one global shift (median of detected − nominal for the closest pairs, robust to missing traces via a coarse 1-D cross-correlation of delta-combs), assigns `role`, `telescope`, `layout_index`, `aperture` to matched traces within `tolerance` px, fills `missing` (layout indices without a match) and `extra` (detected without a match). Unmatched detections keep `role=None`.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_trace.py`:
```python
import numpy as np
import pytest
from astropy.io import fits

from mastpipe.config.models import TraceParams
from mastpipe.core.model import Trace, TraceSet
from mastpipe.instruments.registry import load_layout
from mastpipe.steps.trace import find_traces, match_layout, robust_sigma


def cropped_flat(out_dir):
    f = fits.getdata(out_dir / "sim_flat_00.fits").astype(float)
    return f[1:399, 1:1599] - 500.0  # bias off, background stays (as in the notebook)


def test_robust_sigma_ignores_outliers():
    v = np.concatenate([np.random.default_rng(0).normal(0, 1, 1000), [1e6]])
    assert 0.8 < robust_sigma(v) < 1.2


def test_find_traces_recovers_truth(sim_night):
    out_dir, truth = sim_night
    ts = find_traces(cropped_flat(out_dir), TraceParams(trace_halfwidth=8, degree=2))
    assert len(ts) == 3 and ts.shape == (398, 1598)
    x = np.linspace(0, 1597, 50)
    for t, coeffs in zip(ts.traces, truth.trace_coeffs, strict=True):
        assert np.max(np.abs(t.y(x) - np.polyval(coeffs, x))) < 0.1
        assert t.rms < 0.05 and t.n_points >= 40 and t.snr > 20


def test_find_traces_blank_image_raises():
    with pytest.raises(ValueError):
        find_traces(np.zeros((50, 200)), TraceParams(trace_halfwidth=3))


def test_find_traces_max_traces_keeps_brightest():
    img = np.zeros((100, 300))
    y = np.arange(100.0)
    for yc, amp in ((20, 100.0), (50, 1000.0), (80, 500.0)):
        img += amp * np.exp(-0.5 * ((y[:, None] - yc) / 1.3) ** 2)
    img += np.random.default_rng(1).normal(0, 1, img.shape)
    ts = find_traces(img, TraceParams(trace_halfwidth=4, center_width=50, max_traces=2,
                                      x_bin=20, degree=1))
    assert sorted(round(t.y_center) for t in ts) == [50, 80]


def mk_trace(i, y):
    return Trace(index=i, coeffs=np.array([0.0, y]), y_center=y, slope=0.0, rms=0.01,
                 n_points=10, snr=100.0, x_samples=np.zeros(1), y_samples=np.zeros(1),
                 kept=np.ones(1, bool))


def test_match_layout_assigns_roles_and_reports_missing_and_extra():
    layout = load_layout("highspec", "lab3")  # nominal 18.4, 190.0, 361.7 (uncropped)
    shift = -1.0 + 0.4  # crop offset plus a small global drift
    ts = TraceSet(traces=[mk_trace(0, 18.4 + shift), mk_trace(1, 190.0 + shift),
                          mk_trace(2, 250.0)], shape=(398, 1598))
    out = match_layout(ts, layout, tolerance=3.0)
    assert out.layout_shift == pytest.approx(shift, abs=0.05)
    roles = {t.index: (t.role, t.layout_index) for t in out}
    assert roles[0] == ("arc", 0) and roles[1] == ("arc", 1) and roles[2] == (None, None)
    assert out.missing == (2,) and out.extra == (2,)


def test_match_layout_survey45_with_two_missing():
    layout = load_layout("highspec", "survey45")
    keep = [i for i in range(45) if i not in (7, 30)]
    ts = TraceSet(traces=[mk_trace(k, float(layout.y_nominal[i]) + 2.0) for k, i in
                          enumerate(keep)], shape=(400, 1600))
    out = match_layout(ts, layout, tolerance=3.0)
    assert out.missing == (7, 30) and out.extra == ()
    assert out.layout_shift == pytest.approx(2.0, abs=0.05)
    sci = out.by_role("science")
    assert len(sci) == 19 and all(t.telescope for t in sci)
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_trace.py -v`

- [ ] **Step 3: Implement**

```python
"""Find spectral traces on a flat and identify them against the fiber layout table."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.signal import find_peaks

from mastpipe.config.models import TraceParams
from mastpipe.core.model import Trace, TraceSet


def robust_sigma(v: np.ndarray) -> float:
    v = np.asarray(v, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return 0.0
    return 1.4826 * float(np.median(np.abs(v - np.median(v))))


def _collapse(img: np.ndarray, has_nan: bool) -> np.ndarray:
    if not has_nan:
        return np.median(img, axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(img, axis=1)


def _polyfit_clipped(x, y, deg, clip_sigma=3.0, n_iter=3):
    keep = np.ones(x.size, dtype=bool)
    for _ in range(n_iter if clip_sigma else 0):
        if keep.sum() < deg + 2:
            break
        resid = y - np.polyval(np.polyfit(x[keep], y[keep], deg), x)
        s = robust_sigma(resid[keep])
        if s <= 0:
            break
        new = np.abs(resid - np.median(resid[keep])) < clip_sigma * s
        if new.sum() < deg + 2 or np.array_equal(new, keep):
            break
        keep = new
    return np.polyfit(x[keep], y[keep], deg), keep


def _centroid(profile, y_guess, halfwidth, margin, min_height, min_finite=0.75):
    lo = int(round(y_guess)) - halfwidth - margin
    hi = int(round(y_guess)) + halfwidth + margin + 1
    if lo < 0 or hi > profile.size:
        return None
    seg = profile[lo:hi]
    edges = np.concatenate([seg[:margin], seg[-margin:]])
    inner = seg[margin:seg.size - margin]
    finite = np.isfinite(inner)
    if finite.mean() < min_finite or not np.isfinite(edges).any():
        return None
    base = np.nanmedian(edges)
    w = np.clip(np.where(finite, inner - base, 0.0), 0.0, None)
    total = w.sum()
    if total <= 0 or w.max() < min_height:
        return None
    y_pix = np.arange(lo + margin, lo + margin + w.size, dtype=float)
    return float((y_pix * w).sum() / total), float(w.max())


def _fit_one_trace(peak_y, x_samples, blocks, halfwidth, margin, min_height, deg, clip_sigma,
                   min_points, n_passes=2):
    coeffs = np.array([float(peak_y)])
    result, n_last = None, 0
    for _ in range(n_passes):
        xs, ys = [], []
        for xb, prof in zip(x_samples, blocks, strict=True):
            got = _centroid(prof, np.polyval(coeffs, xb), halfwidth, margin, min_height)
            if got is not None:
                xs.append(xb)
                ys.append(got[0])
        n_last = len(xs)
        if n_last < max(min_points, deg + 2):
            return result, n_last
        xs, ys = np.asarray(xs), np.asarray(ys)
        coeffs, keep = _polyfit_clipped(xs, ys, deg, clip_sigma)
        result = (coeffs, xs, ys, keep)
    return result, n_last


def find_traces(image: np.ndarray, params: TraceParams) -> TraceSet:
    flat = np.asarray(image, dtype=float)
    ny, nx = flat.shape
    x_center = nx // 2
    has_nan = bool(np.isnan(flat).any())
    hw = params.trace_halfwidth
    min_sep = params.min_separation if params.min_separation is not None else 2 * hw
    bw = params.baseline_window if params.baseline_window is not None else 8 * hw + 1
    bw = int(bw) | 1
    margin = max(2, hw)

    lo, hi = max(0, x_center - params.center_width), min(nx, x_center + params.center_width)
    profile = _collapse(flat[:, lo:hi], has_nan)
    baseline = median_filter(profile, size=bw, mode="nearest")
    residual = profile - baseline
    noise = robust_sigma(residual)
    if noise <= 0:
        raise ValueError("the detection profile has zero noise: blank or constant image")
    height = params.threshold_sigma * noise
    peaks, _ = find_peaks(residual, height=height, distance=min_sep)
    if params.max_traces is not None and peaks.size > params.max_traces:
        brightest = np.argsort(residual[peaks])[::-1][: params.max_traces]
        peaks = np.sort(peaks[brightest])

    x_starts = np.arange(0, nx - params.x_bin + 1, params.x_bin)
    blocks = np.array([_collapse(flat[:, s:s + params.x_bin], has_nan) for s in x_starts])
    x_samples = x_starts + (params.x_bin - 1) / 2.0
    block_noise = noise * np.sqrt(max(1.0, (hi - lo) / params.x_bin))
    c_sigma = params.centroid_sigma if params.centroid_sigma is not None else min(
        3.0, params.threshold_sigma)
    min_height = c_sigma * block_noise

    traces: list[Trace] = []
    for peak_y in peaks:
        fitted, _ = _fit_one_trace(peak_y, x_samples, blocks, hw, margin, min_height,
                                   params.degree, params.clip_sigma, params.min_points)
        if fitted is None:
            continue
        coeffs, xs, ys, keep = fitted
        slope = float(np.polyval(np.polyder(coeffs), x_center)) if coeffs.size > 1 else 0.0
        traces.append(Trace(
            index=-1, coeffs=coeffs, y_center=float(np.polyval(coeffs, x_center)), slope=slope,
            rms=float(np.std(ys[keep] - np.polyval(coeffs, xs[keep]))), n_points=int(keep.sum()),
            snr=float(residual[peak_y] / noise), x_samples=xs, y_samples=ys, kept=keep))
    traces.sort(key=lambda t: t.y_center)
    for i, t in enumerate(traces):
        t.index = i
    return TraceSet(traces=traces, shape=(ny, nx))


def _global_shift(detected: np.ndarray, nominal: np.ndarray, search: float = 50.0) -> float:
    """Shift maximising coincidences between two position combs, refined by the median offset."""
    grid = np.arange(-search, search + 0.5, 0.5)
    score = [np.sum(np.min(np.abs(detected[:, None] - (nominal[None, :] + s)), axis=1) < 1.5)
             for s in grid]
    s0 = float(grid[int(np.argmax(score))])
    d = detected[:, None] - (nominal[None, :] + s0)  # (n_detected, n_layout)
    nearest = np.abs(d).argmin(axis=1)
    offsets = d[np.arange(d.shape[0]), nearest]
    close = np.abs(offsets) < 3.0
    if not close.any():
        return s0
    return s0 + float(np.median(offsets[close]))


def match_layout(traceset: TraceSet, layout: pd.DataFrame, tolerance: float) -> TraceSet:
    nominal = layout.y_nominal.to_numpy(dtype=float)
    detected = np.array([t.y_center for t in traceset.traces], dtype=float)
    if detected.size == 0:
        return TraceSet(traceset.traces, traceset.shape, None, tuple(range(len(layout))), ())
    shift = _global_shift(detected, nominal)
    target = nominal + shift
    assigned: dict[int, int] = {}
    for li in np.argsort(np.abs(nominal - np.median(nominal))):  # centre outwards
        d = np.abs(detected - target[li])
        d[list(assigned.values())] = np.inf
        j = int(np.argmin(d))
        if d[j] <= tolerance:
            assigned[int(li)] = j
    traces = []
    for t in traceset.traces:
        traces.append(Trace(**{**t.__dict__}))
    for li, j in assigned.items():
        row = layout.iloc[li]
        t = traces[j]
        t.role = str(row.role)
        t.telescope = None if row.telescope is None or pd.isna(row.telescope) else str(row.telescope)
        t.layout_index = int(row.trace)
        t.aperture = None if pd.isna(row.aperture) else float(row.aperture)
    missing = tuple(int(layout.trace[i]) for i in range(len(layout)) if i not in assigned)
    extra = tuple(t.index for t in traces if t.role is None)
    return TraceSet(traces=traces, shape=traceset.shape, layout_shift=shift, missing=missing,
                    extra=extra)
```

- [ ] **Step 4: Run tests** — `pytest tests/unit/test_trace.py -v` — Expected: 6 passed.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat(steps): trace finding (notebook port) and layout matching"`

---

### Task 9: Aperture extraction and continuum fit (port of `extract_lsf`, `fit_continuum`)

**Files:**
- Create: `src/mastpipe/steps/extract.py`
- Test: `tests/unit/test_extract.py`

**Interfaces:**
- Consumes: `Frame`, `TraceSet`, `Spectrum` (Task 2); `ExtractParams` (Task 3); `robust_sigma` (Task 8).
- Produces:
  - `aperture_weights(y_center: np.ndarray, halfwidth: float, aperture: str) -> tuple[np.ndarray, np.ndarray]` (rows, weights[rows, columns]).
  - `extract_aperture(frame: Frame, traceset: TraceSet, params: ExtractParams) -> list[Spectrum]` — per trace: `flux` = weighted sum (or mean), `error` = sqrt of weighted-sum of variance (or of mean), `mask` = columns with coverage < `min_coverage` or with any masked pixel weight; per-trace `aperture` override from the layout is honoured; `role/telescope/target` copied (target from `frame.meta.target` for science traces).
  - `fit_continuum(spectrum: Spectrum, degree: int, clip_above: float, clip_below: float, n_iter: int = 10) -> Spectrum` — returns a **new** Spectrum with `continuum` set and `flux` replaced by the continuum-subtracted flux; masked pixels are ignored in the fit.
  - `continuum_subtract_all(spectra: list[Spectrum], params: ExtractParams) -> list[Spectrum]`.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_extract.py`:
```python
import numpy as np
import pytest
from astropy.io import fits

from mastpipe.config.models import CalibParams, DetectorConfig, ExtractParams, TraceParams
from mastpipe.core.model import Frame, FrameMeta, Spectrum
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.extract import (
    aperture_weights,
    continuum_subtract_all,
    extract_aperture,
    fit_continuum,
)
from mastpipe.steps.trace import find_traces, match_layout

DET = DetectorConfig(name="main", shape=(400, 1600), gain=1.0, read_noise=5.0,
                     pixel_size_um=16.0, crop=(1, 1599, 1, 399))


def test_aperture_weights_round_and_subpixel():
    yc = np.array([10.0, 10.4, 10.6])
    rows, w = aperture_weights(yc, 2.0, "round")
    assert w.sum(axis=0).tolist() == [5.0, 5.0, 5.0]
    rows2, w2 = aperture_weights(yc, 2.0, "subpixel")
    assert np.allclose(w2.sum(axis=0), 4.0)
    with pytest.raises(ValueError):
        aperture_weights(yc, 2.0, "triangle")


def _calibrated(out_dir, name, kind, exptime, layout):
    def load(n, k, t):
        d = fits.getdata(out_dir / n).astype(np.float32)[1:399, 1:1599]
        return Frame(data=d, meta=FrameMeta(instrument="highspec", camera="main",
                                            exposure_type=k, exposure_time=t, checksum=n))
    bias = build_master([load(f"sim_bias_{i:02d}.fits", "bias", 0.0) for i in range(3)], "bias",
                        CalibParams())
    cs = build_calibration_set(bias, None, None, DET)
    return calibrate_frame(load(name, kind, exptime), cs, subtract_background=False)


def test_extract_arc_recovers_line_positions_and_errors(sim_night, lab3_layout):
    out_dir, truth = sim_night
    flat = _calibrated(out_dir, "sim_flat_00.fits", "flat", 5.0, lab3_layout)
    arc = _calibrated(out_dir, "sim_arc_00.fits", "arc", 6.0, lab3_layout)
    ts = match_layout(find_traces(flat.data, TraceParams()), lab3_layout, 3.0)
    spectra = extract_aperture(arc, ts, ExtractParams(halfwidth=9.0))
    assert len(spectra) == 3 and all(s.role == "arc" for s in spectra)
    s = spectra[1]
    assert s.flux.shape == (1598,) and np.all(s.error[~s.mask] > 0)
    # the brightest simulated line is a local maximum of the extracted arc within 1 px
    px = int(round(truth.line_pixels[1][int(np.argmax(truth.extras["line_amp"]))]))
    win = s.flux[px - 3:px + 4]
    assert np.argmax(win) in (2, 3, 4)
    # aperture sum error ~ sqrt(sum variance): check one column against direct computation
    col = 400
    rows, w = aperture_weights(np.atleast_1d(ts.traces[1].y(col)), 9.0, "round")
    direct = np.sqrt(np.sum(arc.variance[rows, col] * w[:, 0] ** 2))
    assert s.error[col] == pytest.approx(direct, rel=1e-5)


def test_extract_masks_low_coverage_and_bad_pixels():
    data = np.ones((20, 30), dtype=np.float32)
    var = np.ones_like(data)
    mask = np.zeros_like(data, bool)
    mask[10, 5] = True
    frame = Frame(data=data, meta=FrameMeta(instrument="t", camera="m", exposure_type="science",
                                            checksum="x"), variance=var, mask=mask)
    from mastpipe.core.model import Trace, TraceSet
    edge = Trace(index=0, coeffs=np.array([0.0, 1.0]), y_center=1.0, slope=0, rms=0, n_points=1,
                 snr=1, x_samples=np.zeros(1), y_samples=np.zeros(1), kept=np.ones(1, bool))
    mid = Trace(index=1, coeffs=np.array([0.0, 10.0]), y_center=10.0, slope=0, rms=0, n_points=1,
                snr=1, x_samples=np.zeros(1), y_samples=np.zeros(1), kept=np.ones(1, bool))
    spectra = extract_aperture(frame, TraceSet([edge, mid], (20, 30)), ExtractParams(halfwidth=3))
    assert spectra[0].mask.all()  # aperture runs off the top of the frame
    assert spectra[1].mask[5] and not spectra[1].mask[6]
    assert spectra[1].flux[6] == pytest.approx(7.0)


def test_fit_continuum_subtracts_lower_envelope():
    x = np.arange(1000.0)
    cont = 100 + 0.05 * x
    flux = cont.copy()
    flux[[100, 300, 500, 700]] += 5000  # emission lines
    s = Spectrum(fiber=0, pixel=x, flux=flux, error=np.ones_like(x), mask=np.zeros(1000, bool))
    out = fit_continuum(s, degree=3, clip_above=2.0, clip_below=3.0)
    assert out is not s and out.continuum is not None
    assert np.allclose(out.continuum, cont, atol=0.5)
    assert out.flux[100] == pytest.approx(5000, abs=1)
    assert continuum_subtract_all([s], ExtractParams(continuum_degree=3))[0].continuum is not None
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_extract.py -v`

- [ ] **Step 3: Implement**

```python
"""Aperture extraction along traces and continuum fitting. Port of the notebook's LSF cell."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from mastpipe.config.models import ExtractParams
from mastpipe.core.model import Frame, Spectrum, TraceSet
from mastpipe.steps.trace import robust_sigma

APERTURES = ("round", "subpixel", "floor")


def aperture_weights(y_center: np.ndarray, halfwidth: float, aperture: str = "round"):
    if aperture not in APERTURES:
        raise ValueError(f"aperture must be one of {APERTURES}")
    y_center = np.asarray(y_center, dtype=float)
    lo = int(np.floor(np.min(y_center) - halfwidth)) - 1
    hi = int(np.ceil(np.max(y_center) + halfwidth)) + 2
    rows = np.arange(lo, hi)
    r = rows[:, None].astype(float)
    if aperture == "subpixel":
        w = np.clip(np.minimum(r + 0.5, y_center + halfwidth)
                    - np.maximum(r - 0.5, y_center - halfwidth), 0.0, 1.0)
    else:
        centre = np.round(y_center) if aperture == "round" else np.floor(y_center)
        w = (np.abs(r - centre) <= round(halfwidth)).astype(float)
    return rows, w


def extract_aperture(frame: Frame, traceset: TraceSet, params: ExtractParams) -> list[Spectrum]:
    data = np.asarray(frame.data, dtype=float)
    ny, nx = data.shape
    var = np.zeros_like(data) if frame.variance is None else np.asarray(frame.variance, float)
    bad = np.zeros_like(data, bool) if frame.mask is None else frame.mask
    bad = bad | ~np.isfinite(data)
    x = np.arange(nx, dtype=float)
    out: list[Spectrum] = []
    for t in traceset.traces:
        hw = params.halfwidth if t.aperture is None else t.aperture
        rows, w = aperture_weights(t.y(x), hw, params.aperture)
        total = w.sum(axis=0)
        on = (rows >= 0) & (rows < ny)
        ww = w[on]
        strip, vstrip, bstrip = data[rows[on]], var[rows[on]], bad[rows[on]]
        good = ~bstrip
        flux = np.sum(np.where(good, strip * ww, 0.0), axis=0)
        fvar = np.sum(np.where(good, vstrip * ww**2, 0.0), axis=0)
        used = np.sum(np.where(good, ww, 0.0), axis=0)
        coverage = np.divide(used, total, out=np.zeros_like(total), where=total > 0)
        if params.method == "mean":
            flux = np.divide(flux, used, out=np.full_like(flux, np.nan), where=used > 0)
            fvar = np.divide(fvar, used**2, out=np.full_like(fvar, np.nan), where=used > 0)
        mask = (coverage < params.min_coverage) | np.any(bstrip & (ww > 0), axis=0)
        out.append(Spectrum(
            fiber=t.index, pixel=np.arange(nx), flux=flux, error=np.sqrt(fvar), mask=mask,
            role=t.role, telescope=t.telescope,
            target=frame.meta.target if t.role == "science" else None))
    return out


def fit_continuum(spectrum: Spectrum, degree: int, clip_above: float, clip_below: float,
                  n_iter: int = 10) -> Spectrum:
    y = np.asarray(spectrum.flux, dtype=float)
    x = np.asarray(spectrum.pixel, dtype=float)
    allowed = np.isfinite(y) & ~spectrum.mask
    keep = allowed.copy()
    cheb = None
    for _ in range(n_iter):
        if keep.sum() < degree + 2:
            break
        cheb = np.polynomial.Chebyshev.fit(x[keep], y[keep], degree)
        resid = np.where(allowed, y - cheb(x), np.nan)
        sigma = robust_sigma(resid[keep])
        if sigma <= 0:
            break
        new = allowed & (resid < clip_above * sigma) & (resid > -clip_below * sigma)
        if new.sum() < degree + 2 or np.array_equal(new, keep):
            break
        keep = new
    if cheb is None:
        raise ValueError(f"fiber {spectrum.fiber}: too few usable points for a degree "
                         f"{degree} continuum")
    continuum = cheb(x)
    return replace(spectrum, flux=np.where(np.isfinite(y), y - continuum, np.nan),
                   continuum=continuum)


def continuum_subtract_all(spectra: list[Spectrum], params: ExtractParams) -> list[Spectrum]:
    return [fit_continuum(s, params.continuum_degree, params.clip_above, params.clip_below)
            for s in spectra]
```

- [ ] **Step 4: Run tests** — `pytest tests/unit/test_extract.py -v` — Expected: 4 passed.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "feat(steps): aperture extraction with variance/masks and Chebyshev continuum fit"`

---

### Task 10: Provenance and FITS products

**Files:**
- Create: `src/mastpipe/pipeline/__init__.py`, `src/mastpipe/pipeline/provenance.py`, `src/mastpipe/io/products.py`
- Test: `tests/unit/test_products.py`

**Interfaces:**
- Consumes: `MasterCalib`, `TraceSet`, `Trace`, `ExposureSpectra`, `Spectrum`, `FrameMeta`, `QCReport` (Task 2); `EffectiveConfig`, `config_hash` (Task 3).
- Produces:
  - `Provenance(pipeline_version, config_hash, inputs: list[str], calibrations: list[str], wavelength_status: str = "none", steps: list[str] = ())` with `to_header_cards() -> list[tuple[str, str, str]]` (keys `MPVER`, `MPCFGH`, `MPWSTAT`, `MPNIN`, `MPIN0001…`, `MPNCAL`, `MPCAL001…`, `MPSTEPnn`) and `from_header(header) -> Provenance`.
  - `product_name(meta: FrameMeta, kind: str, ext: str = ".fits") -> str` → `"{instrument}_{configuration or 'noconf'}_{camera}_{YYYYmmddTHHMMSS or 'notime'}_{kind}{ext}"` (mid_time parsed with pandas, UTC).
  - `write_master(path, master: MasterCalib, prov: Provenance) -> Path`; `read_master(path) -> tuple[MasterCalib, Provenance]` — primary image + `VARIANCE` extension; validity as numbered key/value cards `MPVALK01`/`MPVALV01`… (name, JSON-encoded value). FITS keywords are at most eight characters, so every scheme here is numbered rather than abbreviated.
  - `write_traceset(path, ts: TraceSet, prov: Provenance) -> Path`; `read_traceset(path) -> tuple[TraceSet, Provenance]` — binary table `TRACES` with columns `index, layout_index, role, telescope, y_center, slope, rms, n_points, snr, aperture, degree, coeffs(4)` (coeffs right-aligned to degree 3, padded with zeros on the left; `degree` says how many are real) plus a second table `SAMPLES` with `index, x, y, kept`; header `MPSHIFT`, `MPMISS`, `MPEXTRA` (comma-joined), shape as `MPNY`, `MPNX`.
  - `write_spectra(path, es: ExposureSpectra, prov: Provenance) -> Path`; `read_spectra(path) -> tuple[ExposureSpectra, Provenance]` — one binary table `SPECTRA`, one row per fiber, vector columns `pixel, flux, error, mask, wavelength, sky, continuum` (NaN-filled when absent, plus boolean presence cards `MPHASWAV`, `MPHASSKY`, `MPHASCON`), scalar columns `fiber, role, telescope, target`; `FrameMeta` fields as numbered cards `MPMK01`/`MPMV01`… (field name, JSON-encoded value), so `None`, tuples and numbers round-trip exactly.
  - `write_qc(path, report: QCReport) -> Path` (JSON, `indent=1`, `sort_keys=True`).

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_products.py`:
```python
import json

import numpy as np
from astropy.io import fits

from mastpipe.core.model import (
    ExposureSpectra,
    FrameMeta,
    MasterCalib,
    QCReport,
    Spectrum,
    Trace,
    TraceSet,
)
from mastpipe.io.products import (
    product_name,
    read_master,
    read_spectra,
    read_traceset,
    write_master,
    write_qc,
    write_spectra,
    write_traceset,
)
from mastpipe.pipeline.provenance import Provenance

META = FrameMeta(instrument="highspec", camera="main", exposure_type="science",
                 configuration="mgb518", exposure_time=300.0, mid_time="2026-10-01T20:15:30",
                 checksum="abc:0", source="sim_science_00.fits", target="SIMSTAR0")
PROV = Provenance(pipeline_version="0.0.1", config_hash="deadbeef0123", inputs=["abc:0"],
                  calibrations=["highspec_mgb518_main_20261001T180000_bias.fits"],
                  steps=["bias", "dark_rate x 300 s"])


def test_product_name_from_metadata():
    assert product_name(META, "spec") == "highspec_mgb518_main_20261001T201530_spec.fits"
    meta2 = FrameMeta(instrument="x", camera="c", exposure_type="bias", checksum="k")
    assert product_name(meta2, "bias") == "x_noconf_c_notime_bias.fits"


def test_provenance_roundtrip_through_header():
    h = fits.Header()
    for k, v, c in PROV.to_header_cards():
        h[k] = (v, c)
    back = Provenance.from_header(h)
    assert back == PROV


def test_master_roundtrip(tmp_path):
    m = MasterCalib(kind="bias", data=np.full((4, 5), 500.0, np.float32),
                    validity={"camera": "main", "binning": (1, 1), "readout_mode": "Image",
                              "detector_temperature": -50.0, "exposure_time": 0.0,
                              "mid_time": "2026-10-01T18:00:00", "configuration": "mgb518"},
                    n_frames=5, provenance={"inputs": ["a", "b"]},
                    variance=np.ones((4, 5), np.float32))
    p = write_master(tmp_path / "m.fits", m, PROV)
    back, prov = read_master(p)
    assert back.kind == "bias" and back.n_frames == 5
    assert back.validity == m.validity and np.array_equal(back.data, m.data)
    assert np.array_equal(back.variance, m.variance) and prov == PROV


def test_traceset_roundtrip(tmp_path):
    t = Trace(index=0, coeffs=np.array([1e-6, 0.002, 18.0]), y_center=19.6, slope=0.002, rms=0.01,
              n_points=60, snr=80.0, x_samples=np.array([12.0, 37.0]), y_samples=np.array([18.0,
              18.1]), kept=np.array([True, False]), role="arc", telescope=None, layout_index=0,
              aperture=None)
    u = Trace(**{**t.__dict__, "index": 1, "role": "science", "telescope": "T01",
                 "layout_index": 3, "aperture": 7.0})
    ts = TraceSet([t, u], (398, 1598), layout_shift=-0.6, missing=(2,), extra=())
    back, prov = read_traceset(write_traceset(tmp_path / "t.fits", ts, PROV))
    assert len(back) == 2 and back.shape == (398, 1598)
    assert back.layout_shift == -0.6 and back.missing == (2,) and back.extra == ()
    b0, b1 = back.traces
    assert np.allclose(b0.coeffs, t.coeffs) and b0.role == "arc" and b0.telescope is None
    assert b1.telescope == "T01" and b1.aperture == 7.0 and b1.layout_index == 3
    assert np.array_equal(b0.kept, t.kept) and np.allclose(b0.x_samples, t.x_samples)
    assert prov == PROV


def test_spectra_roundtrip(tmp_path):
    n = 10
    s = Spectrum(fiber=0, pixel=np.arange(n), flux=np.arange(n, dtype=float),
                 error=np.ones(n), mask=np.arange(n) % 3 == 0, role="science", telescope="T01",
                 target="SIMSTAR0", continuum=np.ones(n))
    es = ExposureSpectra(meta=META, spectra=[s], wavelength_status="none")
    back, prov = read_spectra(write_spectra(tmp_path / "s.fits", es, PROV))
    assert back.meta == META and back.wavelength_status == "none" and prov == PROV
    b = back.spectra[0]
    assert np.array_equal(b.flux, s.flux) and np.array_equal(b.mask, s.mask)
    assert b.wavelength is None and b.sky is None and np.array_equal(b.continuum, s.continuum)
    assert b.telescope == "T01" and b.target == "SIMSTAR0"


def test_write_qc(tmp_path):
    p = write_qc(tmp_path / "q.json", QCReport(stage="trace", values={"n": 3, "rms": [0.1]}))
    assert json.loads(p.read_text(encoding="utf-8")) == {"stage": "trace",
                                                        "values": {"n": 3, "rms": [0.1]}}
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_products.py -v`

- [ ] **Step 3: Implement `provenance.py`**

```python
"""Provenance records carried in every product header."""

from __future__ import annotations

from dataclasses import dataclass, field

from astropy.io import fits


@dataclass(frozen=True)
class Provenance:
    pipeline_version: str
    config_hash: str
    inputs: list[str] = field(default_factory=list)
    calibrations: list[str] = field(default_factory=list)
    wavelength_status: str = "none"
    steps: list[str] = field(default_factory=list)

    def to_header_cards(self) -> list[tuple[str, str, str]]:
        cards = [
            ("MPVER", self.pipeline_version, "mastpipe version"),
            ("MPCFGH", self.config_hash, "effective configuration hash"),
            ("MPWSTAT", self.wavelength_status, "wavelength status"),
            ("MPNIN", str(len(self.inputs)), "number of input frames"),
        ]
        cards += [(f"MPIN{i + 1:04d}", v, "input checksum") for i, v in enumerate(self.inputs)]
        cards.append(("MPNCAL", str(len(self.calibrations)), "number of calibration files"))
        cards += [(f"MPCAL{i + 1:03d}", v, "calibration file") for i, v in
                  enumerate(self.calibrations)]
        cards += [(f"MPSTEP{i + 1:02d}", v, "processing step") for i, v in enumerate(self.steps)]
        return cards

    @classmethod
    def from_header(cls, header: fits.Header) -> Provenance:
        n_in = int(header.get("MPNIN", 0))
        n_cal = int(header.get("MPNCAL", 0))
        steps = []
        i = 1
        while f"MPSTEP{i:02d}" in header:
            steps.append(header[f"MPSTEP{i:02d}"])
            i += 1
        return cls(
            pipeline_version=header["MPVER"],
            config_hash=header["MPCFGH"],
            inputs=[header[f"MPIN{i + 1:04d}"] for i in range(n_in)],
            calibrations=[header[f"MPCAL{i + 1:03d}"] for i in range(n_cal)],
            wavelength_status=header.get("MPWSTAT", "none"),
            steps=steps,
        )
```

- [ ] **Step 4: Implement `products.py`**

```python
"""FITS product writers and readers. Every product carries a Provenance record."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.core.model import (
    ExposureSpectra,
    FrameMeta,
    MasterCalib,
    QCReport,
    Spectrum,
    Trace,
    TraceSet,
)
from mastpipe.pipeline.provenance import Provenance

MAX_DEGREE = 3
META_KEYS = ("instrument", "camera", "exposure_type", "configuration", "exposure_time", "mid_time",
             "binning", "readout_mode", "detector_temperature", "lamp_state", "fiber_preset",
             "target", "checksum", "source", "extension", "unresolved", "classification_rule")


def _write_kv(header: fits.Header, prefix: str, items: dict) -> None:
    """Numbered key/value cards: <prefix>K01 = name, <prefix>V01 = JSON value (8-char keys)."""
    for i, (key, value) in enumerate(items.items(), 1):
        header[f"{prefix}K{i:02d}"] = key
        header[f"{prefix}V{i:02d}"] = json.dumps(value)


def _read_kv(header: fits.Header, prefix: str) -> dict:
    out = {}
    i = 1
    while f"{prefix}K{i:02d}" in header:
        value = json.loads(header[f"{prefix}V{i:02d}"])
        out[header[f"{prefix}K{i:02d}"]] = tuple(value) if isinstance(value, list) else value
        i += 1
    return out


def product_name(meta: FrameMeta, kind: str, ext: str = ".fits") -> str:
    if meta.mid_time:
        stamp = pd.Timestamp(meta.mid_time).strftime("%Y%m%dT%H%M%S")
    else:
        stamp = "notime"
    return f"{meta.instrument}_{meta.configuration or 'noconf'}_{meta.camera}_{stamp}_{kind}{ext}"


def _add_prov(header: fits.Header, prov: Provenance) -> None:
    for k, v, c in prov.to_header_cards():
        header[k] = (v, c)


def _none_to_str(v):
    return "" if v is None else v


def _str_to_none(v):
    return None if v in ("", None) else v


# ---------------------------------------------------------------- masters
def write_master(path: Path, master: MasterCalib, prov: Provenance) -> Path:
    primary = fits.PrimaryHDU(master.data.astype(np.float32))
    h = primary.header
    h["MPKIND"] = master.kind
    h["MPNFRM"] = master.n_frames
    _write_kv(h, "MPVAL", master.validity)
    h["MPPROV"] = json.dumps(master.provenance)
    _add_prov(h, prov)
    hdus = [primary]
    if master.variance is not None:
        hdus.append(fits.ImageHDU(master.variance.astype(np.float32), name="VARIANCE"))
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return Path(path)


def read_master(path: Path) -> tuple[MasterCalib, Provenance]:
    with fits.open(path, memmap=False) as hdul:
        h = hdul[0].header
        validity = _read_kv(h, "MPVAL")
        var = hdul["VARIANCE"].data.copy() if "VARIANCE" in hdul else None
        master = MasterCalib(kind=h["MPKIND"], data=hdul[0].data.copy(), validity=validity,
                             n_frames=int(h["MPNFRM"]), provenance=json.loads(h["MPPROV"]),
                             variance=var)
        return master, Provenance.from_header(h)


# ---------------------------------------------------------------- trace sets
def write_traceset(path: Path, ts: TraceSet, prov: Provenance) -> Path:
    n = len(ts)
    coeffs = np.zeros((n, MAX_DEGREE + 1))
    for i, t in enumerate(ts.traces):
        coeffs[i, MAX_DEGREE + 1 - t.coeffs.size:] = t.coeffs
    cols = [
        fits.Column("index", "J", array=[t.index for t in ts]),
        fits.Column("layout_index", "J", array=[-1 if t.layout_index is None else t.layout_index
                                                for t in ts]),
        fits.Column("role", "8A", array=[_none_to_str(t.role) for t in ts]),
        fits.Column("telescope", "16A", array=[_none_to_str(t.telescope) for t in ts]),
        fits.Column("y_center", "D", array=[t.y_center for t in ts]),
        fits.Column("slope", "D", array=[t.slope for t in ts]),
        fits.Column("rms", "D", array=[t.rms for t in ts]),
        fits.Column("n_points", "J", array=[t.n_points for t in ts]),
        fits.Column("snr", "D", array=[t.snr for t in ts]),
        fits.Column("aperture", "D", array=[np.nan if t.aperture is None else t.aperture
                                            for t in ts]),
        fits.Column("degree", "J", array=[t.coeffs.size - 1 for t in ts]),
        fits.Column("coeffs", f"{MAX_DEGREE + 1}D", array=coeffs),
    ]
    traces_hdu = fits.BinTableHDU.from_columns(cols, name="TRACES")
    idx = np.concatenate([np.full(t.x_samples.size, t.index) for t in ts]) if n else np.array([])
    xs = np.concatenate([t.x_samples for t in ts]) if n else np.array([])
    ys = np.concatenate([t.y_samples for t in ts]) if n else np.array([])
    kept = np.concatenate([t.kept for t in ts]) if n else np.array([], bool)
    samples_hdu = fits.BinTableHDU.from_columns([
        fits.Column("index", "J", array=idx), fits.Column("x", "D", array=xs),
        fits.Column("y", "D", array=ys), fits.Column("kept", "L", array=kept)], name="SAMPLES")
    primary = fits.PrimaryHDU()
    h = primary.header
    h["MPNY"], h["MPNX"] = ts.shape
    h["MPSHIFT"] = -999.0 if ts.layout_shift is None else float(ts.layout_shift)
    h["MPMISS"] = ",".join(str(i) for i in ts.missing)
    h["MPEXTRA"] = ",".join(str(i) for i in ts.extra)
    _add_prov(h, prov)
    fits.HDUList([primary, traces_hdu, samples_hdu]).writeto(path, overwrite=True)
    return Path(path)


def _split_ints(s: str) -> tuple[int, ...]:
    return tuple(int(v) for v in s.split(",") if v != "")


def read_traceset(path: Path) -> tuple[TraceSet, Provenance]:
    with fits.open(path, memmap=False) as hdul:
        h = hdul[0].header
        tab = hdul["TRACES"].data
        smp = hdul["SAMPLES"].data
        traces = []
        for row in tab:
            c = np.asarray(row["coeffs"], dtype=float)[MAX_DEGREE - int(row["degree"]):]
            sel = smp["index"] == row["index"]
            traces.append(Trace(
                index=int(row["index"]), coeffs=c, y_center=float(row["y_center"]),
                slope=float(row["slope"]), rms=float(row["rms"]), n_points=int(row["n_points"]),
                snr=float(row["snr"]), x_samples=np.asarray(smp["x"][sel], float),
                y_samples=np.asarray(smp["y"][sel], float), kept=np.asarray(smp["kept"][sel], bool),
                role=_str_to_none(str(row["role"]).strip()),
                telescope=_str_to_none(str(row["telescope"]).strip()),
                layout_index=None if row["layout_index"] < 0 else int(row["layout_index"]),
                aperture=None if np.isnan(row["aperture"]) else float(row["aperture"])))
        shift = None if h["MPSHIFT"] == -999.0 else float(h["MPSHIFT"])
        ts = TraceSet(traces, (int(h["MPNY"]), int(h["MPNX"])), shift, _split_ints(h["MPMISS"]),
                      _split_ints(h["MPEXTRA"]))
        return ts, Provenance.from_header(h)


# ---------------------------------------------------------------- spectra
def _vec(spectra: list[Spectrum], name: str, n: int, dtype) -> np.ndarray:
    out = np.full((len(spectra), n), np.nan if dtype is float else 0, dtype=dtype)
    for i, s in enumerate(spectra):
        v = getattr(s, name)
        if v is not None:
            out[i] = v
    return out


def write_spectra(path: Path, es: ExposureSpectra, prov: Provenance) -> Path:
    n = es.spectra[0].pixel.size if es.spectra else 0
    has = {k: any(getattr(s, k) is not None for s in es.spectra)
           for k in ("wavelength", "sky", "continuum")}
    cols = [
        fits.Column("fiber", "J", array=[s.fiber for s in es.spectra]),
        fits.Column("role", "8A", array=[_none_to_str(s.role) for s in es.spectra]),
        fits.Column("telescope", "16A", array=[_none_to_str(s.telescope) for s in es.spectra]),
        fits.Column("target", "64A", array=[_none_to_str(s.target) for s in es.spectra]),
        fits.Column("pixel", f"{n}J", array=_vec(es.spectra, "pixel", n, int)),
        fits.Column("flux", f"{n}D", array=_vec(es.spectra, "flux", n, float)),
        fits.Column("error", f"{n}D", array=_vec(es.spectra, "error", n, float)),
        fits.Column("mask", f"{n}L", array=_vec(es.spectra, "mask", n, bool)),
        fits.Column("wavelength", f"{n}D", array=_vec(es.spectra, "wavelength", n, float)),
        fits.Column("sky", f"{n}D", array=_vec(es.spectra, "sky", n, float)),
        fits.Column("continuum", f"{n}D", array=_vec(es.spectra, "continuum", n, float)),
    ]
    table = fits.BinTableHDU.from_columns(cols, name="SPECTRA")
    primary = fits.PrimaryHDU()
    h = primary.header
    _write_kv(h, "MPM", {key: getattr(es.meta, key) for key in META_KEYS})
    h["MPWSTAT"] = es.wavelength_status
    h["MPHASWAV"], h["MPHASSKY"], h["MPHASCON"] = has["wavelength"], has["sky"], has["continuum"]
    _add_prov(h, prov)
    fits.HDUList([primary, table]).writeto(path, overwrite=True)
    return Path(path)


def read_spectra(path: Path) -> tuple[ExposureSpectra, Provenance]:
    with fits.open(path, memmap=False) as hdul:
        h = hdul[0].header
        meta = FrameMeta(**_read_kv(h, "MPM"))
        tab = hdul["SPECTRA"].data
        spectra = []
        for row in tab:
            spectra.append(Spectrum(
                fiber=int(row["fiber"]), pixel=np.asarray(row["pixel"]),
                flux=np.asarray(row["flux"], float), error=np.asarray(row["error"], float),
                mask=np.asarray(row["mask"], bool),
                role=_str_to_none(str(row["role"]).strip()),
                telescope=_str_to_none(str(row["telescope"]).strip()),
                target=_str_to_none(str(row["target"]).strip()),
                wavelength=np.asarray(row["wavelength"], float) if h["MPHASWAV"] else None,
                sky=np.asarray(row["sky"], float) if h["MPHASSKY"] else None,
                continuum=np.asarray(row["continuum"], float) if h["MPHASCON"] else None))
        es = ExposureSpectra(meta=meta, spectra=spectra, wavelength_status=h["MPWSTAT"])
        return es, Provenance.from_header(h)


def write_qc(path: Path, report: QCReport) -> Path:
    Path(path).write_text(json.dumps({"stage": report.stage, "values": report.values}, indent=1,
                                     sort_keys=True), encoding="utf-8")
    return Path(path)
```

- [ ] **Step 5: Run tests** — `pytest tests/unit/test_products.py -v` — Expected: 6 passed.

- [ ] **Step 6: Commit** — `git add -A && git commit -m "feat(io): provenance records and FITS products for masters, trace sets, spectra, QC JSON"`

---

### Task 11: Calibration library, quick-look runner and the `quicklook` command

**Files:**
- Create: `src/mastpipe/pipeline/callib.py`, `src/mastpipe/pipeline/runner.py`; Modify: `src/mastpipe/cli.py`
- Test: `tests/unit/test_callib.py`, `tests/e2e/test_quicklook.py`

**Interfaces:**
- Consumes: everything above.
- Produces:
  - `CalibLibrary(root: Path)` with `add(master: MasterCalib, prov: Provenance, meta_for_name: FrameMeta) -> Path` (writes with `write_master`, appends to `index.json`), `find(kind, camera, binning, readout_mode, configuration=None, temperature=None, temperature_tolerance=2.0, near_time=None) -> tuple[MasterCalib, Path] | None` (exact match on camera/binning/readout_mode; configuration must match when given for `flat`; temperature within tolerance when both known; nearest `mid_time` when several), `list(kind=None) -> list[dict]`.
  - `RunOptions(instrument: str, configuration: str | None = None, camera: str | None = None, layout: str | None = None, overrides: list[dict] = (), on_error: Literal["continue","raise"] = "continue", force: bool = False, manifest: Path | None = None)`.
  - `run_quicklook(night_dir: Path, out_dir: Path, options: RunOptions) -> RunSummary` where `RunSummary(products: list[Path], skipped: list[str], failed: list[tuple[str, str]], qc: list[Path])`.
  - CLI: `mastpipe quicklook <night_dir> --out <dir> --instrument highspec [--configuration mgb518] [--camera main] [--layout lab3] [--manifest file.csv] [--on-error continue|raise] [--force] [--set key.path=value ...]`.
- Recipe (per `(camera, configuration)` group; configuration falls back to the option, then to the only configuration defined by the instrument, else the group is skipped with a reason):
  1. Scan folder → `FrameRef`s; apply manifest.
  2. Masters: bias (all bias frames), dark (all dark frames), background (if any), flat (all flats). Each master: load frames, `build_master`, store in the calibration library under `out_dir/calib`. Flats are calibrated with `subtract_background=params.calib.flat_subtract_background` **before** stacking (bias + dark only by default).
  3. Trace set: `find_traces(master_flat.data, params.trace)` → `match_layout(..., layout, params.trace.layout_tolerance)`; write `..._traces.fits`; QC JSON with per-trace rms, snr, `layout_shift`, `missing`, `extra`.
  4. For each arc and science frame: `load_frame` → `calibrate_frame(subtract_background=params.calib.subtract_background if a background master exists else False)` → `extract_aperture` → arcs additionally `continuum_subtract_all` (continuum stored in the product, flux stays the **raw** extraction for science) → `ExposureSpectra(wavelength_status="none")` → `write_spectra` to `out_dir/spec/` + QC JSON (per-fiber median S/N, masked fraction).
  5. Idempotence: before writing a spectra product, if the file exists and its `MPIN0001` equals the frame checksum and `MPCFGH` equals the current hash and not `force`, skip and record in `skipped`.
  6. Error policy: with `on_error="continue"` any exception in step 4 for one frame is caught, logged (`logging.getLogger("mastpipe")`) and recorded in `failed`; steps 2–3 failures abort the group with a recorded failure. With `"raise"` exceptions propagate.
  7. Logging: `out_dir/logs/quicklook_<UTC stamp>.log` via a `FileHandler` plus console `INFO`.

- [ ] **Step 1: Write the failing tests**

`tests/unit/test_callib.py`:
```python
import numpy as np

from mastpipe.core.model import FrameMeta, MasterCalib
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.provenance import Provenance

PROV = Provenance(pipeline_version="0.0.1", config_hash="abc")


def master(kind, temp, t, mid, conf="mgb518"):
    return MasterCalib(kind=kind, data=np.zeros((2, 2), np.float32),
                       validity={"camera": "main", "binning": (1, 1), "readout_mode": "Image",
                                 "detector_temperature": temp, "exposure_time": t,
                                 "mid_time": mid, "configuration": conf}, n_frames=1)


def meta(kind, mid):
    return FrameMeta(instrument="highspec", camera="main", exposure_type=kind, checksum="c",
                     configuration="mgb518", mid_time=mid)


def test_add_find_nearest_in_time_and_temperature(tmp_path):
    lib = CalibLibrary(tmp_path)
    lib.add(master("bias", -50.0, 0.0, "2026-10-01T18:00:00"), PROV, meta("bias", "2026-10-01T18:00:00"))
    lib.add(master("bias", -50.0, 0.0, "2026-10-02T18:00:00"), PROV, meta("bias", "2026-10-02T18:00:00"))
    lib.add(master("bias", -30.0, 0.0, "2026-10-02T19:00:00"), PROV, meta("bias", "2026-10-02T19:00:00"))
    got = lib.find("bias", "main", (1, 1), "Image", temperature=-50.0,
                   near_time="2026-10-02T20:00:00")
    assert got is not None
    m, path = got
    assert m.validity["mid_time"] == "2026-10-02T18:00:00" and path.exists()
    assert lib.find("bias", "main", (2, 2), "Image") is None
    assert len(lib.list("bias")) == 3


def test_flat_requires_matching_configuration(tmp_path):
    lib = CalibLibrary(tmp_path)
    lib.add(master("flat", -50.0, 5.0, "2026-10-01T18:00:00", conf="halpha"), PROV,
            meta("flat", "2026-10-01T18:00:00"))
    assert lib.find("flat", "main", (1, 1), "Image", configuration="mgb518") is None
    assert lib.find("flat", "main", (1, 1), "Image", configuration="halpha") is not None


def test_index_survives_reopen(tmp_path):
    lib = CalibLibrary(tmp_path)
    lib.add(master("bias", -50.0, 0.0, "2026-10-01T18:00:00"), PROV, meta("bias", "2026-10-01T18:00:00"))
    assert len(CalibLibrary(tmp_path).list()) == 1
```

`tests/e2e/test_quicklook.py`:
```python
import json

import numpy as np

from mastpipe.cli import main
from mastpipe.io.products import read_spectra, read_traceset
from mastpipe.pipeline.runner import RunOptions, run_quicklook


def test_quicklook_end_to_end(sim_night, tmp_path):
    night, truth = sim_night
    out = tmp_path / "out"
    summary = run_quicklook(night, out, RunOptions(instrument="highspec", layout="lab3"))
    assert summary.failed == []
    kinds = sorted(p.name.split("_")[-1] for p in summary.products)
    assert kinds.count("spec.fits") == 4  # 2 arcs + 2 science
    assert "traces.fits" in kinds and "bias.fits" in kinds and "flat.fits" in kinds

    (traces_path,) = [p for p in summary.products if p.name.endswith("_traces.fits")]
    ts, prov = read_traceset(traces_path)
    assert len(ts) == 3 and ts.missing == () and prov.config_hash
    x = np.linspace(0, 1597, 20)
    for t, coeffs in zip(ts.traces, truth.trace_coeffs, strict=True):
        assert np.max(np.abs(t.y(x) - np.polyval(coeffs, x))) < 0.1

    spec_paths = sorted(p for p in summary.products if p.name.endswith("_spec.fits"))
    es, prov = read_spectra(spec_paths[0])
    assert es.wavelength_status == "none" and len(es.spectra) == 3
    assert prov.inputs and any("bias" in c for c in prov.calibrations)
    arc = [e for e in (read_spectra(p)[0] for p in spec_paths) if e.meta.exposure_type == "arc"][0]
    assert all(s.continuum is not None for s in arc.spectra)
    # brightest simulated line is a local maximum of the extracted arc in every fiber
    k_line = int(np.argmax(truth.extras["line_amp"]))
    for k, s in enumerate(arc.spectra):
        px = int(round(truth.line_pixels[k][k_line]))
        assert np.argmax(s.flux[px - 3:px + 4]) in (2, 3, 4)

    qc = [json.loads(p.read_text(encoding="utf-8")) for p in summary.qc]
    assert any(q["stage"] == "trace" and q["values"]["missing"] == [] for q in qc)

    # idempotent: second run skips the spectra
    again = run_quicklook(night, out, RunOptions(instrument="highspec", layout="lab3"))
    assert len(again.skipped) == 4 and again.failed == []
    logs = list((out / "logs").glob("quicklook_*.log"))
    assert logs and "skipped" in logs[-1].read_text(encoding="utf-8").lower()


def test_quicklook_continues_past_a_bad_frame(sim_night, tmp_path):
    night, _ = sim_night
    bad = tmp_path / "night2"
    bad.mkdir()
    for p in night.glob("*.fits"):
        (bad / p.name).write_bytes(p.read_bytes())
    (bad / "sim_science_00.fits").write_bytes(b"not a fits file")
    out = tmp_path / "out2"
    summary = run_quicklook(bad, out, RunOptions(instrument="highspec", layout="lab3"))
    assert len(summary.failed) == 1 and "sim_science_00.fits" in summary.failed[0][0]
    assert sum(p.name.endswith("_spec.fits") for p in summary.products) == 3


def test_cli_quicklook(sim_night, tmp_path, capsys):
    night, _ = sim_night
    rc = main(["quicklook", str(night), "--out", str(tmp_path / "o"), "--instrument", "highspec",
               "--layout", "lab3", "--set", "extract.halfwidth=8"])
    assert rc == 0
    assert "products" in capsys.readouterr().out
```

- [ ] **Step 2: Run to verify failure** — `pytest tests/unit/test_callib.py tests/e2e -v`

- [ ] **Step 3: Implement `callib.py`**

```python
"""On-disk calibration library: master frames with validity keys and a JSON index."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from mastpipe.core.model import FrameMeta, MasterCalib
from mastpipe.io.products import product_name, read_master, write_master
from mastpipe.pipeline.provenance import Provenance


class CalibLibrary:
    def __init__(self, root: Path):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.index_path = self.root / "index.json"
        self._index: list[dict] = []
        if self.index_path.exists():
            self._index = json.loads(self.index_path.read_text(encoding="utf-8"))

    def _save(self) -> None:
        self.index_path.write_text(json.dumps(self._index, indent=1), encoding="utf-8")

    def add(self, master: MasterCalib, prov: Provenance, meta_for_name: FrameMeta) -> Path:
        name = product_name(meta_for_name, master.kind)
        path = write_master(self.root / name, master, prov)
        entry = {"file": name, "kind": master.kind}
        entry.update({k: (list(v) if isinstance(v, tuple) else v)
                      for k, v in master.validity.items()})
        self._index = [e for e in self._index if e["file"] != name] + [entry]
        self._save()
        return path

    def list(self, kind: str | None = None) -> list[dict]:
        return [e for e in self._index if kind is None or e["kind"] == kind]

    def find(self, kind: str, camera: str, binning: tuple[int, int], readout_mode: str | None,
             configuration: str | None = None, temperature: float | None = None,
             temperature_tolerance: float = 2.0, near_time: str | None = None):
        cands = []
        for e in self.list(kind):
            if e.get("camera") != camera or tuple(e.get("binning", ())) != tuple(binning):
                continue
            if e.get("readout_mode") != readout_mode:
                continue
            if kind == "flat" and configuration is not None and e.get("configuration") != configuration:
                continue
            t_e = e.get("detector_temperature")
            if temperature is not None and t_e is not None and abs(t_e - temperature) > temperature_tolerance:
                continue
            cands.append(e)
        if not cands:
            return None
        if near_time is not None:
            ref = pd.Timestamp(near_time)
            cands.sort(key=lambda e: abs(pd.Timestamp(e["mid_time"]) - ref)
                       if e.get("mid_time") else pd.Timedelta.max)
        else:
            cands.sort(key=lambda e: e.get("mid_time") or "", reverse=True)
        path = self.root / cands[0]["file"]
        master, _ = read_master(path)
        return master, path
```

- [ ] **Step 4: Implement `runner.py`**

```python
"""The quick-look recipe: night folder -> masters, traces, extracted spectra, QC."""

from __future__ import annotations

import logging
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.__version__ import __version__
from mastpipe.config.loader import build_effective, config_hash, deep_merge
from mastpipe.config.models import EffectiveConfig
from mastpipe.core.model import ExposureSpectra, FrameRef, QCReport, TraceSet
from mastpipe.instruments.registry import load_header_map, load_instrument, load_layout
from mastpipe.io.headers import HeaderMap
from mastpipe.io.ingest import load_frame, read_manifest, scan_folder
from mastpipe.io.products import product_name, write_qc, write_spectra, write_traceset
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.calib import (
    build_calibration_set,
    build_master,
    calibrate_frame,
    calibration_steps,
)
from mastpipe.steps.extract import continuum_subtract_all, extract_aperture
from mastpipe.steps.trace import find_traces, match_layout

log = logging.getLogger("mastpipe")


@dataclass
class RunOptions:
    instrument: str
    configuration: str | None = None
    camera: str | None = None
    layout: str | None = None
    overrides: list[dict] = field(default_factory=list)
    on_error: Literal["continue", "raise"] = "continue"
    force: bool = False
    manifest: Path | None = None


@dataclass
class RunSummary:
    products: list[Path] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)
    failed: list[tuple[str, str]] = field(default_factory=list)
    qc: list[Path] = field(default_factory=list)


def _setup_logging(out_dir: Path) -> logging.Handler:
    (out_dir / "logs").mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S")
    handler = logging.FileHandler(out_dir / "logs" / f"quicklook_{stamp}.log", encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(handler)
    log.setLevel(logging.INFO)
    return handler


def _already_done(path: Path, checksum: str, cfg_hash: str) -> bool:
    if not path.exists():
        return False
    h = fits.getheader(path)
    return h.get("MPIN0001") == checksum and h.get("MPCFGH") == cfg_hash


def _group_key(ref: FrameRef, default_conf: str | None) -> tuple[str, str | None]:
    return ref.meta.camera, ref.meta.configuration or default_conf


def run_quicklook(night_dir: Path, out_dir: Path, options: RunOptions) -> RunSummary:
    night_dir, out_dir = Path(night_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    handler = _setup_logging(out_dir)
    summary = RunSummary()
    try:
        instrument = load_instrument(options.instrument)
        header_map = HeaderMap.from_dict(load_header_map(options.instrument))
        manifest = read_manifest(options.manifest) if options.manifest else None
        unreadable: list[Path] = []
        refs = scan_folder(night_dir, options.instrument, header_map, manifest=manifest,
                           unreadable=unreadable)
        for path in unreadable:
            summary.failed.append((path.name, "unreadable FITS"))
            log.error("unreadable FITS file skipped: %s", path.name)
        log.info("scanned %d frames in %s", len(refs), night_dir.name)
        default_conf = options.configuration or (
            next(iter(instrument.configurations)) if len(instrument.configurations) == 1 else None)
        groups: dict[tuple[str, str | None], list[FrameRef]] = defaultdict(list)
        for ref in refs:
            groups[_group_key(ref, default_conf)].append(ref)
        for (camera, conf), group in sorted(groups.items(), key=str):
            if options.camera and camera != options.camera:
                continue
            if conf is None:
                summary.failed.append((f"group {camera}", "configuration unknown; use --configuration"))
                log.error("camera %s: configuration unknown, group skipped", camera)
                continue
            try:
                eff = build_effective(instrument, conf, camera, options.overrides, options.layout)
                _run_group(group, eff, out_dir, options, summary)
            except Exception as exc:  # noqa: BLE001 - group-level failure policy
                if options.on_error == "raise":
                    raise
                summary.failed.append((f"group {camera}/{conf}", repr(exc)))
                log.exception("group %s/%s failed", camera, conf)
        log.info("done: %d products, %d skipped, %d failed", len(summary.products),
                 len(summary.skipped), len(summary.failed))
    finally:
        log.removeHandler(handler)
        handler.close()
    return summary


def _run_group(group: list[FrameRef], eff: EffectiveConfig, out_dir: Path, options: RunOptions,
               summary: RunSummary) -> None:
    cfg_hash = config_hash(eff)
    det, params = eff.detector, eff.params
    by_type: dict[str, list[FrameRef]] = defaultdict(list)
    for ref in group:
        by_type[ref.meta.exposure_type].append(ref)
    lib = CalibLibrary(out_dir / "calib")
    spec_dir, qc_dir = out_dir / "spec", out_dir / "qa"
    spec_dir.mkdir(exist_ok=True)
    qc_dir.mkdir(exist_ok=True)

    def prov(inputs, calibs, steps=()):
        return Provenance(__version__, cfg_hash, list(inputs), list(calibs), "none", list(steps))

    # --- masters -----------------------------------------------------------------
    masters: dict[str, tuple] = {}
    for kind in ("bias", "dark", "background"):
        if by_type.get(kind):
            frames = [load_frame(r, det) for r in by_type[kind]]
            m = build_master(frames, kind, params.calib)
            path = lib.add(m, prov([f.meta.checksum for f in frames], []), frames[0].meta)
            masters[kind] = (m, path)
            summary.products.append(path)
            log.info("master %s from %d frames -> %s", kind, len(frames), path.name)
    if "bias" not in masters:
        raise RuntimeError("no bias frames: cannot calibrate")
    calib = build_calibration_set(masters["bias"][0], masters.get("dark", (None,))[0],
                                  masters.get("background", (None,))[0], det)
    calib_files = [p.name for _, p in masters.values()]

    if not by_type.get("flat"):
        raise RuntimeError("no flat frames: cannot find traces")
    flats = [calibrate_frame(load_frame(r, det), calib,
                             subtract_background=params.calib.flat_subtract_background
                             and calib.background is not None, saturation=det.saturation)
             for r in by_type["flat"]]
    mflat = build_master(flats, "flat", params.calib)
    flat_path = lib.add(mflat, prov([f.meta.checksum for f in flats], calib_files), flats[0].meta)
    summary.products.append(flat_path)

    # --- traces -------------------------------------------------------------------
    layout = load_layout(eff.instrument, eff.layout) if eff.layout else None
    if layout is None:
        raise RuntimeError("no fiber layout selected; pass --layout or set default_layout")
    ts: TraceSet = match_layout(find_traces(mflat.data, params.trace), layout,
                                params.trace.layout_tolerance)
    traces_path = out_dir / "calib" / product_name(flats[0].meta, "traces")
    write_traceset(traces_path, ts, prov([f.meta.checksum for f in flats],
                                         calib_files + [flat_path.name]))
    summary.products.append(traces_path)
    qc_path = write_qc(qc_dir / traces_path.with_suffix(".json").name, QCReport("trace", {
        "n_traces": len(ts), "layout_shift": ts.layout_shift, "missing": list(ts.missing),
        "extra": list(ts.extra), "rms": [t.rms for t in ts], "snr": [t.snr for t in ts],
        "roles": [t.role for t in ts]}))
    summary.qc.append(qc_path)
    log.info("traces: %d found, shift %.2f px, missing %s, extra %s", len(ts), ts.layout_shift or 0,
             ts.missing, ts.extra)

    # --- arcs and science ----------------------------------------------------------
    for ref in by_type.get("arc", []) + by_type.get("science", []):
        target = spec_dir / product_name(ref.meta, "spec")
        if not options.force and _already_done(target, ref.meta.checksum, cfg_hash):
            summary.skipped.append(ref.meta.source)
            log.info("skipped %s (already reduced)", ref.meta.source)
            continue
        try:
            frame = load_frame(ref, det)
            sub_bkg = params.calib.subtract_background and calib.background is not None
            cal = calibrate_frame(frame, calib, subtract_background=sub_bkg,
                                  saturation=det.saturation)
            spectra = extract_aperture(cal, ts, params.extract)
            if ref.meta.exposure_type == "arc":
                spectra = continuum_subtract_all(spectra, params.extract)
            es = ExposureSpectra(meta=ref.meta, spectra=spectra, wavelength_status="none")
            steps = calibration_steps(calib, ref.meta.exposure_time or 0.0, sub_bkg)
            steps += [f"extract {params.extract.aperture} +/-{params.extract.halfwidth}"]
            write_spectra(target, es, prov([ref.meta.checksum], calib_files + [flat_path.name,
                                                                                traces_path.name],
                                           steps))
            summary.products.append(target)
            snr = [float(np.nanmedian(np.abs(s.flux[~s.mask]) / s.error[~s.mask]))
                   if (~s.mask).any() else float("nan") for s in spectra]
            qc_path = write_qc(qc_dir / target.with_suffix(".json").name, QCReport("extract", {
                "source": ref.meta.source, "exposure_type": ref.meta.exposure_type,
                "median_snr": snr, "masked_fraction": [float(s.mask.mean()) for s in spectra]}))
            summary.qc.append(qc_path)
            log.info("reduced %s -> %s", ref.meta.source, target.name)
        except Exception as exc:  # noqa: BLE001 - per-frame failure policy
            if options.on_error == "raise":
                raise
            summary.failed.append((ref.meta.source, repr(exc)))
            log.exception("failed %s", ref.meta.source)


def parse_set(values: list[str]) -> list[dict]:
    """--set extract.halfwidth=8 -> [{"extract": {"halfwidth": 8}}], typed via YAML rules."""
    import yaml

    out = []
    for item in values:
        key, _, raw = item.partition("=")
        node: dict = {}
        cur = node
        parts = key.split(".")
        for p in parts[:-1]:
            cur = cur.setdefault(p, {})
        cur[parts[-1]] = yaml.safe_load(raw)
        out.append(node)
    return out
```

- [ ] **Step 5: Wire the CLI**

In `cli.py`, add the sub-command:
```python
    ql = sub.add_parser("quicklook", help="reduce a night folder non-interactively")
    ql.add_argument("night_dir", type=Path)
    ql.add_argument("--out", type=Path, required=True)
    ql.add_argument("--instrument", required=True)
    ql.add_argument("--configuration")
    ql.add_argument("--camera")
    ql.add_argument("--layout")
    ql.add_argument("--manifest", type=Path)
    ql.add_argument("--on-error", choices=["continue", "raise"], default="continue")
    ql.add_argument("--force", action="store_true")
    ql.add_argument("--set", action="append", default=[], metavar="KEY=VALUE")
```
and in `main`:
```python
    if args.command == "quicklook":
        from mastpipe.pipeline.runner import RunOptions, parse_set, run_quicklook

        opts = RunOptions(instrument=args.instrument, configuration=args.configuration,
                          camera=args.camera, layout=args.layout, overrides=parse_set(args.set),
                          on_error=args.on_error, force=args.force, manifest=args.manifest)
        s = run_quicklook(args.night_dir, args.out, opts)
        print(f"{len(s.products)} products, {len(s.skipped)} skipped, {len(s.failed)} failed")
        for name, why in s.failed:
            print(f"  FAILED {name}: {why}")
        return 0 if not s.failed else 2
```
(`sub = parser.add_subparsers(dest="command")` replaces the bare call from Task 1; add `logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")` at the top of `main`.)

- [ ] **Step 6: Run everything** — `pytest -v` — Expected: all green. `test_cli_quicklook` returns 0 because the simulated night has no failures.

- [ ] **Step 7: Commit** — `git add -A && git commit -m "feat(pipeline): calibration library, quick-look runner, mastpipe quicklook command"`

---

### Task 12: Documentation, changelog, decisions, full verification

**Files:**
- Modify: `README.md`, `CHANGES.md`, `DECISIONS.md`
- Create: `notebooks/README.md` (one paragraph: notebooks import the package; the reference notebook is ported in Plan 2/3 as the API grows; strip outputs with `nbstripout` before committing)

- [ ] **Step 1: README** — sections: what mastpipe is (two sentences), install, `mastpipe quicklook` example on a simulated night (`python -c "..."` snippet that calls `simulate_night` into a temp folder, then the CLI), layout of the output tree (`calib/ spec/ qa/ logs/`), pointer to `CLAUDE.md` and to the design documents (say they live in `docs/design/`, which is not tracked yet).
- [ ] **Step 2: CHANGES.md** — under `0.0.1 (unreleased)`: bullets for every task above.
- [ ] **Step 3: DECISIONS.md** — add dated entries: (a) header map + classification rules as the only place raw keywords live; (b) fiber layout as input; (c) products keyed by metadata not raw filenames; (d) provenance keys `MP*`; (e) quick-look failure policy; (f) `docs/design/` git-ignored for now.
- [ ] **Step 4: Full verification**

```bash
ruff check src tests scripts && ruff format --check src tests scripts
python scripts/privacy_scan.py
pytest -v
```
Expected: lint clean, scan silent, all tests pass. Then run the test-suite once more with `-p no:cacheprovider` from a different current directory (e.g. the parent folder, `pytest <repo>/tests`) to catch any implicit dependence on the working directory.

- [ ] **Step 5: Commit** — `git add -A && git commit -m "docs: README, changelog, decisions for the foundation release"`

- [ ] **Step 6: Report to the user** with: the test count, the product tree of a simulated night, the three open items carried into Plan 2 (Horne extraction, wavelength solution, sky) and any deviations from this plan.

---

## Self-review notes (checked against the spec)

- **Spec coverage in this plan:** §2 constraints (Task 1 CI matrix, privacy scan, `pathlib`); §3 layers and layout (Tasks 2–11); §4.1–4.6 (Tasks 4, 5, 11: header map, classification, manifest, MEF split, layered config, layout table, calibration library); §5.1–5.4 (Tasks 7–9; Horne optimal extraction deferred to Plan 2 as stated in the roadmap); §6 quick-look mode (Task 11: non-interactive, degrade-and-continue, idempotent, logging); §7 products 1, 2, 4 and QC JSON (Tasks 10–11); §8 simulator and tiers 1–2 (Tasks 6, 11); §10 hygiene (Task 1). Not in this plan, by design: §5.5–5.10, §7 products 3 and 5, QA figures/HTML, `watch`, `reduce`/`identify`, DeepSpec, PypeIt adapter, tier-3 regression.
- **Type consistency:** `FrameMeta.exposure_type` values are the `EXPOSURE_TYPES` tuple everywhere; `TraceSet.missing/extra` are tuples of ints in Tasks 2, 8, 10; `Spectrum.mask` is bool with True = bad in Tasks 2, 9, 10, 11; `Provenance` is the same dataclass in Tasks 10 and 11; `EffectiveConfig.layout` (Task 3) is what `run_quicklook` uses to pick the layout (Task 11); `parse_set` values are YAML-typed so `--set extract.halfwidth=8` becomes an int, accepted by `PositiveFloat`.
- **Known rough edges to watch during execution:** the survey-45 layout is a placeholder geometry until measured; gain and read noise in the HighSpec definition are placeholders; the corrupt-file end-to-end test relies on `scan_folder`'s `unreadable` list (Task 5) and on the runner recording those paths as failures (Task 11).
