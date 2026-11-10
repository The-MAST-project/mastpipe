# mastpipe tidy-up plan (after Plan 1.5)

**Goal:** split `pipeline/runner.py` into stage modules and close the review items
deferred during Plan 1.5, before DeepSpec work (Plan 2) builds on the runner.

**Source of the items:** a triage of the Plan 1.5 deferred list against `main` after the
second merged pull request. Line numbers below are from that triage and are hints only:
**confirm each item in the current code before changing it**; if an item turns out to be
already fixed or wrong, say so in the report and skip it.

## Global Constraints

- Branch `tidy-up`. Commit there; never push, never merge, never touch `main`.
- Privacy: no local absolute paths, usernames, hostnames, e-mails, tokens or **personal
  names** in any tracked file (code, tests, docs, commit messages). The privacy scan
  (`python scripts/privacy_scan.py`)
  must pass.
- Follow the repository `CLAUDE.md` and `docs/DECISIONS.md`.
- **Docs in sync:** every behaviour change updates, in the same commit, the user docs
  under `docs/` that describe it, the README where it lists the feature, and `CHANGES.md`
  (one line under the unreleased section). A decision that changes a recorded rule gets a
  `docs/DECISIONS.md` bullet.
- Behaviour changes are test-first: write the failing test, see it fail, fix, see it pass.
- Checks before the last commit of a task: `ruff check .` and `ruff format --check .`
  clean; the focused tests for what changed; the whole suite once (it takes more than ten
  minutes — run it in the background with a long timeout, never twice in parallel);
  `sphinx-build -E -W -b html docs <tmp>` clean when anything under `docs/` or any
  docstring changed.
- `filterwarnings = error` stays; no new dependency; Python 3.11 and 3.12; Windows, Linux
  and macOS paths (use `pathlib`, POSIX-style relative sources).
- Public names stay importable where they are today:
  `from mastpipe.pipeline.runner import RunOptions, RunSummary, run_quicklook, NO_FRAMES,
  TIME_POLICIES`.
- Not in scope (do not start): splitting `io/products.py` or `sim/simulator.py`; a
  read-only `CalibLibrary`; caching calibrated `light` frames; avoiding frame loads on
  master reuse; rewording the historical design documents.

## Task 1: Split `pipeline/runner.py` into stage modules (move only)

No behaviour change, no renames of functions, no signature changes, no edits to logic.
Only moves, imports, and the docstrings needed to introduce each module.

Target modules in `src/mastpipe/pipeline/`:

| Module | Receives |
|---|---|
| `group.py` | `RunOptions`, `RunSummary`, `_Group`, the folder/kind/policy constants (`SPEC_FOLDER`, `QC_FOLDER`, `LOG_FOLDER`, `MASTER_KINDS`, `TIME_POLICIES`, `ALREADY_REDUCED`, `NO_FRAMES`), the logger, and the helpers shared by more than one stage: `_describe`, `_time_key`, `_skip`, `_per_frame`, `_by_exposure`, `_median_snr`, `_write_qc`, `_master_meta`, `_calibrated`, `_timezone`, `_master_provenance`, `_science_or_all`, `_group_binning`, `_group_readout`, `_check_binning`, `_unlisted_amplifiers`, `_readout_mismatch`, `_check_readout`, `_same_readout`, `_content_steps` |
| `masters.py` | `_usable_frames`, `_masters_hash`, `_library_master`, `_masters` |
| `traces.py` | `_resolve_light`, `_full_arc_candidates`, `_traces` |
| `arc_stage.py` | `_reference_for`, `_arcs`, `_arc`, `_full_arc`, `_cal_arc` |
| `science.py` | `_arcs_for`, `_arc_weights`, `_blend_drifts`, `_wavelengths`, `_stored_key`, `_science_frame`, `_combine_targets` |
| `runner.py` (stays) | the module docstring (the description of the whole recipe stays here), `_start_logging`, `_run_group`, `EARLIER_PRODUCT`, `_warn_earlier_products`, `run_quicklook`, and re-exports of the public names listed in Global Constraints (with `__all__`) |

Rules:

- If a helper listed under `group.py` is in fact used by one stage only, it may live in
  that stage's module; if a stage function needs a helper from another stage, move the
  helper to `group.py` rather than importing across stages. No import cycles; stage
  modules never import `runner`.
- `_reference_for` is used by the arc stage and (check) the science stage: if both, it
  goes to `group.py`.
- Each new module gets a short module docstring saying which step of the recipe it holds
  and pointing at `mastpipe.pipeline.runner` for the whole picture. Cross-references in
  existing docstrings (`:func:` roles naming moved functions) are updated so Sphinx
  resolves them.
- Tests: update imports of private helpers to their new modules
  (`tests/unit/test_runner.py` and any other test). Any `monkeypatch` target string or
  attribute that names `mastpipe.pipeline.runner.<moved name>` must point at the module
  where the name is now *looked up*. Do not otherwise edit tests.
- Docs: add an API page per new module beside `docs/api/mastpipe-pipeline-runner.md` and
  list them in `docs/api/index.md`, following the existing pages. Update
  `docs/how-it-works/overview.md` (and `pipeline/__init__.py`'s docstring, the README's
  or `CLAUDE.md`'s package map **only if** they list the modules of `pipeline/`;
  `CLAUDE.md` changes are *proposed in the report*, not made).
- `CHANGES.md`: one line (internal restructuring, no behaviour change).
  `docs/DECISIONS.md`: one bullet recording the split and the "stage modules import
  `group`, never each other or `runner`" rule.

Acceptance: whole suite passes with the same pass/skip counts as before the task (record
both); lab regression passes when the lab data variables are set; `ruff` clean;
`sphinx-build -E -W` clean; `git diff --stat` shows no file under `src/` other than
`pipeline/` modules; no module in `pipeline/` longer than about 450 lines except where a
single function forces it (report the line counts).

## Task 2: Masters, library and idempotence fixes

Each with a test first, in the stage modules produced by Task 1.

1. **Flat reuse compares names, not content** (`_library_master`): a library flat is
   reused when its recorded masters match the group's by *name*; compare the hash of the
   masters' file checksums (`_masters_hash`) instead, so a flat calibrated with a
   different bias of the same name is rebuilt. Keep reuse working for flats written
   before this change only if their stored information allows the comparison; otherwise
   they are rebuilt once (say so in `CHANGES.md`). Add the `DECISIONS.md` sentence that
   the comparison is by checksum.
2. **Flat batch chosen by raw count** (`_masters`): the flat frames to stack are chosen
   before unusable frames are dropped, while the other kinds choose among survivors.
   Choose the flat batch among the frames that survive loading and calibration, the same
   way as the other kinds. Add the missing test of a flat frame whose calibration raises.
3. **Flat-derived trace sets are always named `nosum`** (`_master_meta` / `_traces`): give
   the trace product the checksum of the flat master's first input (or of the master
   file) so two trace sets of one stamp cannot share a name; update
   `docs/how-it-works/products-policy.md` where it describes `nosum`, and fix the
   alignment of its fenced block (two columns off).
4. `_library_master` asserts `group.calib_hash is not None` in production code: make the
   field non-optional (it is always set in `_run_group`) and drop the assert.
5. The idempotence key tuple is assembled twice (`_stored_key`, `_science_frame`): one
   helper builds it for both.
6. The 12-hex-character truncation of sha256 digests is written in four places
   (`config/loader.py`, `instruments/registry.py`, the masters hash, `steps/wavecal.py`):
   one named constant in a module all four may import without a cycle (suggest
   `core/`), used by all.
7. `frames = []` in `_masters` gets its annotation; `DriftModel.n_usable` typed
   `int | None` is read with `or 0` in two runner places — make the field `int`
   (default 0) if every constructor can supply it, else leave and say why.
8. Docstrings: `_full_arc_candidates` (a candidate whose `find_traces` raises is withheld
   from the arc stage) and `_traces` (it shrinks `by_type["arc"]`).
9. Add the e2e step: master reuse under a non-calibration `--set` override (today only
   unit-tested), and a unit test that `calib_hash` changes with `readouts`.

## Task 3: Ingest, times and metadata fixes

Each with a test first.

1. **Backslash manifest sources never match** (`io/ingest.py`, source matching): normalise
   `\` to `/` in manifest sources before matching.
2. **`scan_folder(manifest=…)` drops unmatched rows silently**: accept and fill an
   `unmatched` list as `apply_manifest` does (or forward to it); the runner's own path is
   unchanged in behaviour.
3. **Zone labels**: a header-map `timezone` of `"UTC"` or `"Etc/UTC"` labels frames
   `"UTC"`, not `"utc"`, and the label is part of the idempotence key. Normalise to
   `utc` in one place (`core/times.py`).
4. **`OBJECT = 0` is treated as no target**: only a missing card or a blank string means
   no target.
5. `ASSUMED_UTC` lives in `io/ingest.py`; move it to `core/times.py` (keep it importable
   from `ingest` if tests or docs use it there). The frame stamp is parsed twice in
   ingest: parse once and return whether the zone was explicit.
6. `sim/simulator.py` has its own `fromisoformat`: use `core.times.to_utc`.
7. `FrameMeta` accepts duplicate telescope keys in its per-unit mappings: reject in
   `__post_init__` with a clear `ValueError` (confirm first how duplicates could arise;
   if the type makes them impossible, skip and report).
8. `targets` is listed in `unresolved` for frames where nothing can resolve it, including
   calibration frames: do not list it for calibration exposure types. `_build_meta`
   assumes the header map's `targets` entry is per-unit: validate that when the header
   map is loaded, with a clear error.
9. A multi-extension file of another instrument adds one `skipped` entry per HDU: one per
   file. A manifest row naming an existing file with an unsupported extension is reported
   "matches no frame": give it its own reason.
10. `steps/calib.py`: a hand-built `CalibrationSet` with `dark_rate` but no `t_dark`
    drops the dark variance silently: raise `ValueError`.
11. `vocab.py`: the unit-name regex uses `\d` (Unicode digits) — use `[0-9]`;
    `is_unit_name` / `unit_number` hard-wire the project prefix while `unit_name` takes
    it as an argument — make the three consistent (same optional argument, same
    default). Pin the three untested tuples in `tests/unit/test_vocab.py`.
12. Tests only: a DST/ambiguous local time through `to_utc` (fold 0); a FITS card holding
    astropy's `Undefined`; the manifest glob row test asserts the value it sets; the
    duplicated `highspec_map` fixture moves to `conftest.py`.
13. Docs: `io/ingest.py` module docstring says manifest `mid_time` is "in UTC" — an
    offset is honoured; `docs/running/night-folder.md` gets one sentence each on
    symlinked subfolders (not descended) and the limits of the manifest `targets` syntax
    (`;` and `=`); `docs/how-it-works/` products page documents that `MPTZ` may hold a
    comma-joined label for a group of mixed zones.

## Task 4: Provenance, atomic QA plots and messages

1. **Configuration origin** (test first): when `run_quicklook` gives the run's
   `--configuration` to a frame whose header and manifest name no grating, record that in
   the product's provenance steps (a step naming the origin, in the style of the
   content-classification step) and in the frame's QC; document in
   `docs/running/quicklook.md`. Add the missing test of a frame that names *another*
   grating under `--configuration` (it keeps its own).
2. **QA PNGs are written in place** (`io/qa.py`, `savefig`): write through
   `io/atomic.py`'s atomic writer with `format="png"` (test first: no partial file at the
   final name; no stray temporary after success).
3. Messages (assert the new text in the nearest existing test):
   - the group label `group main/None` → `group main/(no configuration)`;
   - `row N: x: matches no frame` double colon → one colon;
   - `instruments/registry.py` configuration error gains the `{path.name}:` prefix, and
     its test matches the message text rather than `GRATING_NAMES|Mg`;
   - `sim/simulator.py`: `load_truth` on an old truth file names the missing array;
     `trace_shifts` for an unknown file raises `KeyError` naming it; `row_of` names a
     duplicated id;
   - `_cal_arc`'s warning repeats the fiber count: say it once;
   - the `MPDUSE` header comment names the gate that failed when a drift model is unused
     (the reason is already in the QC JSON; keep within the FITS comment length).
4. `pipeline/identify.py`: log at INFO which masters were applied (file names), and warn
   when a library bias was rejected only because of the temperature tolerance (test with
   `caplog`).
5. `pipeline/callib.py`: `find_reference` docstring gains the unlocked-read note;
   document in the class docstring that constructing a `CalibLibrary` creates the folder
   and its lock file (write access needed).

## Task 5: Tests, lint, CI and docs leftovers

One batch of small edits; no behaviour change in `src/`.

- Tests: `tests/fleet/test_vocab_matches_common.py` — the `importorskip` reason names
  the package actually missing, and the configuration-list regex gets `re.S`;
  `tests/unit/test_callib.py` two-process tests terminate workers still alive after the
  join timeout; `tests/e2e/test_quicklook.py` — delete the dead `index =`, and the
  background-skip assertion checks the whole rule (`content:no_signal`), not just
  `content:`; `tests/regression/test_lab_highspec.py` — the unused
  `NOTEBOOK_RMS_PX` / `NOTEBOOK_N_LINES` constants are either asserted against or
  deleted (prefer asserting if the test already computes the numbers);
  `tests/unit/test_simulator.py` jitter comment says 2× while the assert is 4× — make
  the comment true; tests that index simulator truth by fiber position
  (`test_extract.py`, `test_sky.py`) use `row_of`; `tests/unit/test_runner.py` test that
  leans on a whole simulated night for one frame — leave, unless trivial.
- Lint: the per-file `N` ignores in `pyproject.toml` are narrowed to the files that need
  them; the seven "yoda" comparisons ruff's `SIM300` produced
  (`approx(x) == value` in `test_wavemodel.py`, `test_geometry.py`, `test_wavecal.py`)
  are restored to `value == approx(x)` and `SIM300` ignored for `tests/`.
- CI (`.github/workflows/ci.yml`): ruff is pinned in the workflow *and* in the `dev`
  extra — install it from the extra so there is one pin; the `fleet` job's
  `git ls-remote` branch probe falls back to the default branch on *any* failure — fall
  back only on the "no such ref" exit status (2) and fail the job on other errors.
- Docs: `docs/running/quicklook.md` lists `Halpha` and `Ca` under `--configuration`
  though the instrument ships only `Mg` — say which exist today;
  `docs/DECISIONS.md` readout bullet wording ("a full arc competing as the trace source"
  also covers arcs of unknown kind); `CHANGES.md` still says `.tmp-<pid>` for the atomic
  writer's temporary name — match the code; lines over 100 characters in `README.md`,
  `docs/running/first-reference.md` and `docs/running/night-folder.md` are rewrapped.
