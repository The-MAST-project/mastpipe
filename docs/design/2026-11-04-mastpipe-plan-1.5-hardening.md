# mastpipe Plan 1.5 — Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the boundary gaps the 2026-11-02 review found, so that a night folder shaped the way `MAST_spec` writes it reduces without a manifest, the package speaks `MAST_common`'s vocabulary, and the unattended run survives bad frames, interrupted writes and concurrent writers.

**Architecture:** The three layers stay: declarative `instruments/`, pure `steps/`, and `pipeline/` + `io/` for files. New pieces are one vocabulary module (`mastpipe.vocab`), one time helper (`mastpipe.core.times`), one atomic-file helper (`mastpipe.io.atomic`), one content classifier (`mastpipe.steps.classify`), and a `light` exposure type that the runner resolves before the trace stage. Everything else is a change inside an existing module.

**Tech Stack:** Python 3.11/3.12, numpy, scipy, astropy, pandas, Pydantic v2, pytest, ruff 0.16.0, GitHub Actions.

**Spec:** `docs/design/2026-11-04-mastpipe-plan-1.5-hardening-design.md` (read it first; section numbers below refer to it).

## Global Constraints

- No person's name, local path, hostname, e-mail or token in any tracked file (`CLAUDE.md`; `scripts/privacy_scan.py` runs in CI). The lead supplies the list of names to scrub; it is never written down in the repository.
- `CLAUDE.md` is edited only with the lead's explicit approval: a task that needs it stops and asks.
- Windows, Linux and macOS: `pathlib` only, explicit encodings, no shell-outs, no POSIX-only calls outside the one guarded lock helper.
- `filterwarnings = error` in pytest; no new warnings.
- The laboratory regression numbers (trace centres, rms 0.016 / 0.033 / 0.012 px, line counts 23 / 25 / 20) are unchanged. Regression tests run only when `MASTPIPE_TEST_DATA` and `MASTPIPE_ATLAS_DIR` are set; the lead stages the data.
- Grating identifiers are case-sensitive and spelled as `MAST_common` does: `Mg`, `Halpha`, `Ca` (`Future`, `Empty` known but unsupported). Telescopes are `mast01`…`mast20`.
- Commit after every task with a message in the repository's style (`feat:` / `fix:` / `docs:` / `test:` prefixes are not used there; plain imperative subjects are). Push nothing: pushes need per-push approval.
- Quick checks: `pytest tests/unit -q` (seconds). Full: `pytest -q` (about 6 minutes; includes the slow end-to-end tests). `ruff format --check . && ruff check .` before every commit.

## Review Focus

1. **Two frames with the same mid-time to the second** (two acquisition folders, or two frames with no time): both must get products, never one overwriting the other. Pinned by `test_product_name_differs_for_equal_times` in Task 5.
2. **A header time with an explicit UTC offset under `assume-utc`**: the offset is honoured and converted, and the product is not flagged as assumed. Pinned by `test_explicit_offset_is_honoured_under_assume_utc` in Task 3.
3. **A manifest row whose bare file name matches two files in different folders**: an error naming both paths, not a silent first-match. Pinned by `test_bare_name_ambiguous_is_an_error` in Task 4.
4. **A `light` frame with a zero exposure time** (a bias the control software took as `light`): a bias, never an `unknown` or a content-classified frame. Pinned by `test_light_with_zero_exposure_is_bias` in Task 6.
5. **A lock held past the timeout** (a hung sibling process): a clear `TimeoutError` naming the lock file, not an indefinite wait. Pinned by `test_lock_times_out_with_a_message` in Task 9.

---

### Task 1: The vocabulary mirror and the fleet job

**Files:**
- Create: `src/mastpipe/vocab.py`
- Create: `tests/fleet/__init__.py`, `tests/fleet/test_vocab_matches_common.py`
- Modify: `.github/workflows/ci.yml` (add the `fleet` job)
- Test: `tests/unit/test_vocab.py`

**Interfaces:**
- Produces: `FRAME_TYPES: tuple[str, ...] = ("light", "bias", "dark", "flat")`; `GRATING_NAMES = ("Ca", "Halpha", "Mg", "Future")`; `DISPERSERS = ("Ca", "Mg", "Halpha", "Empty")`; `SPEC_INSTRUMENTS = ("deepspec", "highspec")`; `LAMP_NAMES = ("ThAr", "qTh")`; `DEEPSPEC_BANDS = ("I", "G", "R", "U")`; `UNIT_COUNT = 20`; `def unit_name(number: int, project: str = "mast") -> str` (`unit_name(1) == "mast01"`, `ValueError` outside 1..20); `def is_unit_name(text: str) -> bool`; `def unit_number(name: str) -> int`.
- The module imports only the standard library.

- [ ] **Step 1: Write the unit tests**

```python
# tests/unit/test_vocab.py
def test_unit_name_is_two_digit():
    assert unit_name(1) == "mast01" and unit_name(20) == "mast20"

@pytest.mark.parametrize("n", [0, 21])
def test_unit_name_rejects_out_of_range(n):
    with pytest.raises(ValueError):
        unit_name(n)

def test_unit_number_round_trips():
    assert all(unit_number(unit_name(n)) == n for n in range(1, 21))
    assert not is_unit_name("T01")

def test_vocabulary_values():
    assert GRATING_NAMES == ("Ca", "Halpha", "Mg", "Future")
    assert FRAME_TYPES == ("light", "bias", "dark", "flat")
    assert LAMP_NAMES == ("ThAr", "qTh")
```

- [ ] **Step 2: Run them to see them fail** — `pytest tests/unit/test_vocab.py -q`; expected: ImportError.

- [ ] **Step 3: Write `src/mastpipe/vocab.py`** with the constants and functions in the Interfaces block, a module docstring saying each constant mirrors which name in `MAST_common` (`common/spec.py`: `FrameType`, `Disperser`, `GratingNames`, `SpecInstruments`, `WheelNames`, `DeepspecBands`; the unit rule from `common/config/site.py`), and the sentence that `mastpipe`'s own exposure types (`arc`, `background`, `unknown`, `science`) are pipeline-only and live in `mastpipe.core.model`.

- [ ] **Step 4: Run the unit tests** — expected PASS.

- [ ] **Step 5: Write the fleet test**

```python
# tests/fleet/test_vocab_matches_common.py
REQUIRED = os.environ.get("MASTPIPE_FLEET_CHECK") == "1"
common_spec = pytest.importorskip("common.spec") if not REQUIRED else importlib.import_module("common.spec")

def _literal_values(alias) -> tuple[str, ...]: ...  # typing.get_args on a Literal

def test_frame_types_match(): assert tuple(m.value for m in common_spec.FrameType) == FRAME_TYPES
def test_grating_names_match(): assert _literal_values(common_spec.GratingNames) == GRATING_NAMES
def test_dispersers_match(): assert _literal_values(common_spec.Disperser) == DISPERSERS
def test_instruments_match(): assert _literal_values(common_spec.SpecInstruments) == SPEC_INSTRUMENTS
def test_lamps_match(): assert _literal_values(common_spec.WheelNames) == LAMP_NAMES
def test_bands_match(): assert _literal_values(common_spec.DeepspecBands) == DEEPSPEC_BANDS

def test_config_lists_match_as_text():
    # common/config/specs.py keeps its lists as locals inside a validator; read the file
    text = Path(common_spec.__file__).with_name("config").joinpath("specs.py").read_text(encoding="utf-8")
    assert set(re.search(r"valid_grating_names = \[(.*?)\]", text).group(1).replace('"', "").replace(" ", "").split(",")) == set(GRATING_NAMES)
    assert set(re.search(r"valid_wheel_names = \[(.*?)\]", text).group(1).replace('"', "").replace(" ", "").split(",")) == set(LAMP_NAMES)

def test_unit_names_are_canonical():
    canonic = importlib.import_module("common.utils").canonic_unit_name
    assert all(canonic(unit_name(n)) == unit_name(n) for n in range(1, 21))
```

When `REQUIRED` and `common` is not importable the module-level import must raise, so the job fails rather than skips.

- [ ] **Step 6: Add the `fleet` job to `.github/workflows/ci.yml`**, copying `MAST_control`'s shape: check this repo out at `path: mastpipe`; resolve the paired `MAST_common` branch (`github.head_ref || github.ref_name`, `git ls-remote --heads`, else `master`); check `The-MAST-project/MAST_common` out at `path: common`; Python 3.12; `pip install -e "./mastpipe[dev]" pydantic`; run `python -m pytest tests/fleet -q -ra` with `working-directory: mastpipe`, `PYTHONPATH: ${{ github.workspace }}` and `MASTPIPE_FLEET_CHECK: "1"`.

- [ ] **Step 7: Run the fleet test locally against a clone** — with a `common/` clone next to a parent folder on `PYTHONPATH` and `MASTPIPE_FLEET_CHECK=1`: `pytest tests/fleet -q`; expected PASS. Without the variable: skipped with the sibling-layout message.

- [ ] **Step 8: Commit** — `git add src/mastpipe/vocab.py tests/unit/test_vocab.py tests/fleet .github/workflows/ci.yml` · `git commit -m "Mirror MAST_common's vocabulary in mastpipe.vocab, checked against the fleet in CI"`.

---

### Task 2: Canonical identifiers (D1)

**Files:**
- Rename: `src/mastpipe/instruments/highspec/gratings/mgb518.yaml` → `Mg.yaml` (its `name:` → `Mg`), `lines/mgb518.csv` → `lines/Mg.csv`, `anchors/mgb518_lab3.yaml` → `anchors/Mg_lab3.yaml`
- Modify: `src/mastpipe/instruments/highspec/header_map.yaml` (`configuration` field), `src/mastpipe/instruments/registry.py` (validate block names), `scripts/make_survey45_layout.py` + regenerated `layouts/survey45.csv`, `tests/conftest.py` (`mgb518_lines` → `mg_lines`, `"mgb518"` → `"Mg"`), every test and doc that names `mgb518` or `T01` (the list: `grep -rn "mgb518\|\"T0\|'T0" src tests docs README.md CLAUDE.md`), `src/mastpipe/pipeline/commands.py` and `identify.py` help strings
- Test: `tests/unit/test_registry.py`, `tests/unit/test_headers.py`

**Interfaces:**
- Produces: configuration ids `Mg` / `Halpha` / `Ca` everywhere; `load_instrument` raises `ValueError` for a grating block whose `name` is not in `vocab.GRATING_NAMES`; `survey45.csv` `telescope` values `mast01`…`mast20`.

- [ ] **Step 1: Write the failing tests**

```python
# tests/unit/test_registry.py
def test_highspec_configurations_use_fleet_names():
    assert set(load_instrument("highspec").configurations) == {"Mg"}

def test_grating_block_with_foreign_name_is_rejected(tmp_path, monkeypatch): ...
    # a copied instrument dir with gratings/mgb518.yaml → pytest.raises(ValueError, match="GRATING_NAMES|Mg")

def test_survey45_telescopes_are_unit_names():
    layout = load_layout("highspec", "survey45")
    tels = [t for t in layout["telescope"] if isinstance(t, str)]
    assert sorted(set(tels)) == [unit_name(n) for n in range(1, 21)]

# tests/unit/test_headers.py
def test_configuration_keeps_case_and_maps_lab_aliases():
    hm = HeaderMap.from_dict(load_header_map("highspec"))
    assert hm.resolve({"GRATING": "Mg"}, "x.fits")[0]["configuration"] == "Mg"
    assert hm.resolve({}, "FOCUS_00_mgb518.fits")[0]["configuration"] == "Mg"
    assert hm.resolve({}, "x_halpha.fits")[0]["configuration"] == "Halpha"
```

- [ ] **Step 2: Run** — `pytest tests/unit/test_registry.py tests/unit/test_headers.py -q`; expected FAIL.

- [ ] **Step 3: Rename the package data files and update the header map**: `configuration: {keys: [GRATING], filename_regex: "(?P<value>mgb518|halpha|caiihk)", map: {mgb518: Mg, halpha: Halpha, caiihk: Ca}}` (no `lower`). In `registry.load_instrument`, after merging blocks, raise `ValueError(f"grating block {name!r} is not a MAST_common grating name {GRATING_NAMES}")`.

- [ ] **Step 4: Regenerate `survey45.csv`** — in `scripts/make_survey45_layout.py` replace `f"T{telescope:02d}"` with `unit_name(telescope)`; run `python scripts/make_survey45_layout.py`.

- [ ] **Step 5: Remove the silent default.** In `runner.run_quicklook` delete the branch that sets `default_conf` from a single-configuration instrument; a group whose configuration is `None` fails with `f"configuration unknown; pass --configuration one of {GRATING_NAMES}"`. Add `test_group_without_configuration_fails_naming_the_gratings` to `tests/unit/test_runner.py` (a copied simulated frame with no `GRATING` card and a file name without an alias → `summary.failed` holds that message, nothing reduced).

- [ ] **Step 6: Replace every remaining `mgb518` and `T01`-style name** in `src/`, `tests/`, `docs/` (including the two notebooks' source cells), `README.md`, `CHANGES.md`; leave `DECISIONS.md` history and the design documents as they are. The simulator's `TARGETS` and `TELESCOPES_PER_TARGET` are unaffected.

- [ ] **Step 7: Run the full suite** — `pytest -q`; expected PASS (tests that relied on the default now pass `configuration="Mg"` or write a `GRATING` card). Then `grep -rn "mgb518" src tests docs README.md` must show only the header-map alias, the changelog line you add in Task 15, and the design documents.

- [ ] **Step 8: Commit** — `git commit -am "Adopt MAST_common's grating and unit names as mastpipe's identifiers"` (after `git add -A src tests docs scripts README.md`).

---

### Task 3: One time helper and the `timezone` option (§5.2.4)

**Files:**
- Create: `src/mastpipe/core/times.py`
- Modify: `src/mastpipe/io/headers.py` (`FIELD_OPTIONS`, `FieldSpec.timezone`, validation), `src/mastpipe/io/ingest.py` (`_mid_time`, `PASSTHROUGH_FIELDS`), `src/mastpipe/core/model.py` (`FrameMeta.timezone`, `make_identifier`), `src/mastpipe/pipeline/callib.py` and `pipeline/arcs.py` (drop `_timestamp`, import the helper), `src/mastpipe/steps/calib.py` (`_mean_time_or_none`), `src/mastpipe/pipeline/provenance.py` (`timezone` → `MPTZ`), `src/mastpipe/pipeline/runner.py` (one warning per run; provenance field), `src/mastpipe/instruments/highspec/header_map.yaml`
- Test: `tests/unit/test_times.py`, `tests/unit/test_headers.py`, `tests/unit/test_ingest.py`

**Interfaces:**
- Produces: `times.to_utc(value: str, timezone: str = "utc") -> pd.Timestamp` (aware, UTC; `timezone` is `"utc"`, `"assume-utc"` or an IANA zone name; a value with its own offset is converted whatever `timezone` says; `ValueError` for an unparsable value or an unknown zone); `times.utc_iso(stamp) -> str` (naive ISO-8601 UTC, seconds, the format every `mid_time` uses today); `times.is_explicit(value: str) -> bool`.
- `FieldSpec.timezone: str = "utc"`; `FrameMeta.timezone: str | None` (`"utc"`, `"assumed-utc"`, or the zone the naive stamp was read in; `None` when no time); `Provenance.timezone: str | None` → card `MPTZ`.
- HighSpec header map: `mid_time: {keys: [DATE-OBS, FRAME, DATE], type: str, reference: start, timezone: assume-utc}`.

- [ ] **Step 1: Write the tests**

```python
# tests/unit/test_times.py
def test_naive_is_read_in_the_named_zone():
    assert utc_iso(to_utc("2026-09-09T11:35:24", "Asia/Jerusalem")) == "2026-09-09T08:35:24"
def test_assume_utc_keeps_the_clock():
    assert utc_iso(to_utc("2026-09-09T11:35:24", "assume-utc")) == "2026-09-09T11:35:24"
def test_explicit_offset_is_honoured_under_assume_utc():
    assert utc_iso(to_utc("2026-11-04T20:00:00+02:00", "assume-utc")) == "2026-11-04T18:00:00"
    assert is_explicit("2026-11-04T20:00:00+02:00") and not is_explicit("2026-11-04T20:00:00")
def test_unknown_zone_is_an_error():
    with pytest.raises(ValueError): to_utc("2026-01-01T00:00:00", "Mars/Olympus")

# tests/unit/test_ingest.py
def test_frame_meta_records_the_timezone_used(...):
    # FRAME card naive under assume-utc → meta.timezone == "assumed-utc"; a DATE-OBS with +02:00 → "utc"
# tests/unit/test_headers.py
def test_timezone_option_is_validated():  # "Mars/Olympus" → ValueError at HeaderMap.from_dict
```

- [ ] **Step 2: Run** — `pytest tests/unit/test_times.py tests/unit/test_ingest.py tests/unit/test_headers.py -q`; expected FAIL.

- [ ] **Step 3: Implement `core/times.py`** (zoneinfo + pandas), add the option to `headers.py` (validate zone names with `zoneinfo.ZoneInfo` at parse time), thread `timezone` through `ingest._mid_time(value, exposure_time, reference, timezone) -> tuple[str | None, bool, str | None]` returning the zone label, add `timezone` to `FrameMeta` (after `mid_time`) and to `PASSTHROUGH_FIELDS`, replace the four parsers listed in Files with `to_utc`, add `timezone` to `Provenance` (card `MPTZ`, comment `"time zone of the frame's clock"`, read as unset when absent), and in `runner.run_quicklook` log one warning `"frame times carry no time zone; taken as UTC (header map 'assume-utc')"` when any scanned frame has `timezone == "assumed-utc"`.

- [ ] **Step 4: Run the unit suite, then the full suite** — `pytest tests/unit -q`, then `pytest -q`; expected PASS (the simulator writes naive `DATE-OBS`, which under `assume-utc` yields the same strings as before).

- [ ] **Step 5: Commit** — `git commit -m "One time helper; header maps declare the zone of a naive stamp (HighSpec: assume-utc, recorded as MPTZ)"`.

---

### Task 4: Recursive scan, instrument gate, manifest patterns (§5.2.1)

**Files:**
- Modify: `src/mastpipe/io/ingest.py` (`scan_folder`, `apply_manifest`, `_row_overrides`, new `_match_sources`), `src/mastpipe/instruments/highspec/header_map.yaml` (`instrument` field), `src/mastpipe/core/model.py` (`FrameMeta.instrument` stays; nothing new), `src/mastpipe/pipeline/runner.py` (`run_quicklook` passes `skipped`, records unmatched manifest rows), `src/mastpipe/pipeline/identify.py` (`_frame_ref` uses `recursive=False`)
- Test: `tests/unit/test_ingest.py`

**Interfaces:**
- Produces: `scan_folder(folder, instrument, header_map, pattern="*.fit*", manifest=None, unreadable=None, recursive=True, skipped=None)`; `source` is the POSIX path relative to `folder`; files are matched case-insensitively against `*.fits`, `*.fit`, `*.fits.fz` (the `pattern` argument keeps its meaning for a single-file scan); a frame whose resolved `instrument` field is set and differs from `instrument` (case-insensitive) is appended to `skipped` as `(source, f"instrument {value!r} is not {instrument!r}")` and left out; hidden entries (`.`-prefixed) are ignored.
- `apply_manifest(refs, manifest, unmatched: list[str] | None = None)`: a row's `source` matches a ref by exact relative path, by glob pattern (`fnmatch` on the relative POSIX path, `*` not crossing `/`, `**` crossing it), or by bare file name when exactly one ref has that name; two or more → `ValueError` naming the candidates; a row matching nothing is appended to `unmatched` as `"manifest row {index}: {source}"`.
- Header map gains `instrument: {keys: [INSTRUME], lower: true}` and **loses** `fiber_preset` (C4: in MAST a fiber preset is a stage position). `fiber_preset` stays a `FrameMeta` field and a manifest column; the layout comes from `--layout`, a manifest, or the instrument default. The simulator stops writing `FIBPRE`; the end-to-end `options()` helper in `tests/e2e/test_quicklook.py` passes `layout="survey45"` (today the frames' `FIBPRE` card chose it).

- [ ] **Step 1: Write the tests**

```python
def test_fibpre_card_is_ignored(tmp_path, highspec_map):  # FIBPRE='survey45' in the header → meta.fiber_preset is None

def test_scan_is_recursive_with_relative_sources(tmp_path, highspec_map):
    # night/highspec/acquisition-0001/exposure-001.fits and acquisition-0002/exposure-001.fits
    sources = {r.meta.source for r in scan_folder(night, "highspec", highspec_map)}
    assert sources == {"highspec/acquisition-0001/exposure-001.fits", "highspec/acquisition-0002/exposure-001.fits"}

def test_non_fits_and_hidden_files_are_not_scanned(...):  # seq.txt, .hidden.fits, plan file → not in refs, not in unreadable

def test_instrument_gate_skips_other_instruments(...):
    # INSTRUME='Highspec' kept; INSTRUME='Deepspec' → in skipped with reason; no INSTRUME → kept

def test_manifest_matches_relative_path_glob_and_unique_name(...):
    # rows: "highspec/acquisition-0001/*.fits" → arc; "exposure-002.fits" unique → flat

def test_bare_name_ambiguous_is_an_error(...):
    with pytest.raises(ValueError, match="acquisition-0001.*acquisition-0002"): ...

def test_unmatched_manifest_rows_are_reported(...):
    unmatched = []; apply_manifest(refs, manifest, unmatched); assert unmatched == ["manifest row 0: nothing.fits"]
```

- [ ] **Step 2: Run** — `pytest tests/unit/test_ingest.py -q`; expected FAIL.

- [ ] **Step 3: Implement** the signatures above. In `run_quicklook`, pass `skipped=` and `unmatched=` lists and append both to `summary.skipped` (with `log.warning`). In `identify._frame_ref`, call `scan_folder(..., recursive=False)` with the file's name as pattern (unchanged behaviour).

- [ ] **Step 4: Run the full suite** — `pytest -q`; expected PASS (the simulator writes flat folders; nested ones are only in the new tests).

- [ ] **Step 5: Commit** — `git commit -m "Scan night folders recursively with relative sources, gate on INSTRUME, match manifest rows by path, glob or unique name"`.

---

### Task 5: Product names that cannot collide (§5.2.1)

**Files:**
- Modify: `src/mastpipe/io/products.py` (`product_name`, `combined_name`), `src/mastpipe/pipeline/identify.py` (`_latest_master`), docs and tests that spell a product name (`grep -rn "_spec.fits\|_traces.fits\|_bias.fits" tests docs`)
- Test: `tests/unit/test_products.py`, `tests/unit/test_identify.py`

**Interfaces:**
- Produces: `product_name(meta, kind, ext=".fits")` = `{instrument}_{configuration|noconf}_{camera}_{stamp|notime}_{meta.checksum[:6] or "nosum"}_{kind}{ext}`; `_latest_master` splits the middle on `_` and compares the stamp part only.

- [ ] **Step 1: Write the tests**

```python
def test_product_name_differs_for_equal_times():
    a = replace(meta, checksum="abcdef0123456789"); b = replace(meta, checksum="123456abcdef0000")
    assert product_name(a, "spec") != product_name(b, "spec")
    assert product_name(a, "spec") == "highspec_Mg_main_20261001T180230_abcdef_spec.fits"

def test_latest_master_reads_the_stamp_before_the_checksum(tmp_path, highspec_effective): ...
```

- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement.** **Step 4:** `pytest -q` — update the literal names in `tests/e2e/test_quicklook.py` (`test_quicklook_continues_after_bad_frame`) and any doc that prints a product name (`docs/running/products.md`, `first-reference.md`); expected PASS.

- [ ] **Step 5: Commit** — `git commit -m "Product names carry a short input checksum so equal mid-times never collide"`.

---

### Task 6: The provisional `light` type and content-based arc detection (D4, §5.2.3)

**Files:**
- Create: `src/mastpipe/steps/classify.py`
- Modify: `src/mastpipe/core/model.py` (`EXPOSURE_TYPES` gains `"light"` before `"science"`), `src/mastpipe/io/headers.py` (`classify`), `src/mastpipe/instruments/highspec/header_map.yaml` (`light: light`; bias rule), `src/mastpipe/config/models.py` (`ClassifyParams`, `StepParams.classify`), `src/mastpipe/pipeline/runner.py` (`_resolve_light`, called in `_run_group` between `_masters` and `_traces`; `science` list and `first` computed after it; `run_quicklook` no longer treats `light` as reducible-as-is), `src/mastpipe/sim/simulator.py` (`IMAGE_TYPES["arc"]` stays `"arc"`; a new keyword argument `site_headers: bool = False`, see Task 16, is not needed yet)
- Test: `tests/unit/test_classify.py`, `tests/unit/test_headers.py`, `tests/unit/test_runner.py`, `tests/regression/test_lab_highspec.py` (one new data-gated test)

**Interfaces:**
- Produces: `classify.arc_likeness(data: np.ndarray, mask: np.ndarray | None, params: ClassifyParams) -> ArcLikeness` with fields `has_signal: bool`, `peak_fraction: float` (fraction of the above-continuum flux of the spatially collapsed spectrum that lies in narrow peaks), `n_peaks: int`; `classify.classify_light(frame: Frame, params: ClassifyParams) -> tuple[str, str]` returning `("arc", "content:arc_like")`, `("science", "content:continuum")` or `("unknown", "content:no_signal")`.
- `ClassifyParams(_Strict)`: `continuum_window: int = 101`, `signal_sigma: float = 5.0`, `arc_peak_fraction: float = 0.5`, `peak_halfwidth_px: int = 4`.
- `HeaderMap.classify`: a header value in `EXPOSURE_TYPES` other than `light` wins as today; `light` runs the rules and takes the first matching rule **that has conditions**; with no such match it stays `("light", "header")`. The unconditional fallback applies only when the header has no valid type.
- Runner: `_resolve_light(group, by_type) -> None` moves each `light` ref (calibrated with the group's masters, no cosmics) into `by_type["arc"]` / `["science"]` with `meta.classification_rule` set, or skips it with its reason; a manifest row that set `exposure_type` never reaches it (its type is not `light`).

- [ ] **Step 1: Write the tests**

```python
# tests/unit/test_classify.py  (uses the session sim_night fixture: a simulated full arc, cal arc, science and flat frame, loaded and bias-subtracted)
def test_simulated_arcs_are_arc_like(sim_frames): assert classify_light(sim_frames.full_arc, P)[0] == "arc" and classify_light(sim_frames.cal_arc, P)[0] == "arc"
def test_simulated_science_is_continuum(sim_frames): assert classify_light(sim_frames.science, P) == ("science", "content:continuum")
def test_flat_is_continuum(sim_frames): assert classify_light(sim_frames.flat, P)[0] == "science"
def test_no_signal_is_unknown(): # bias-level noise frame
    assert classify_light(noise_frame, P) == ("unknown", "content:no_signal")

# tests/unit/test_headers.py
def test_light_is_provisional_and_rules_refine_it():
    hm = HeaderMap.from_dict(load_header_map("highspec"))
    assert hm.classify({"exposure_type": "light", "exposure_time": 30.0}) == ("light", "header")
    assert hm.classify({"exposure_type": "light", "exposure_time": 0.0})[0] == "bias"
def test_light_with_zero_exposure_is_bias(tmp_path, highspec_map): # through scan_folder on a written frame
def test_flat_header_type_is_never_refined(): assert hm.classify({"exposure_type": "flat", "exposure_time": 0.0}) == ("flat", "header")

# tests/unit/test_runner.py
def test_resolve_light_sorts_frames_by_content(sim_night, tmp_path):
    # copy the simulated night, rewrite IMAGETYP of every arc and science frame to 'light', drop LAMP;
    # run_quicklook with layout survey45 → spec products for every science frame, QC "reference"/"drift" for the arcs,
    # every science QC has classification_rule "content:continuum"

# tests/regression/test_lab_highspec.py
def test_lab_focus_frames_are_arc_like(night, effective):  # data-gated: every FOCUS frame → "arc"; the flat → "science"; the background → "unknown"
```

- [ ] **Step 2: Run** — `pytest tests/unit/test_classify.py tests/unit/test_headers.py -q`; expected FAIL.

- [ ] **Step 3: Implement `steps/classify.py`.** Collapse along rows (`nansum` over unmasked pixels per column), continuum = `scipy.ndimage.median_filter(spectrum, continuum_window, mode="nearest")`, residual = spectrum − continuum, `sigma = robust_sigma(residual)` (reuse the helper in `steps/lines.py`), `has_signal = max(spectrum − median) > signal_sigma·sigma·sqrt(rows)`; peaks = columns with residual > `signal_sigma·sigma` grown by `peak_halfwidth_px`; `peak_fraction = residual[peaks].clip(0).sum() / (spectrum − median(spectrum)).clip(0).sum()`; `arc` when `peak_fraction >= arc_peak_fraction`. Tune the defaults until the four unit tests and the lab test pass; record the chosen defaults in the `ClassifyParams` docstring.

- [ ] **Step 4: Implement the header-map rule and the runner stage**, then `pytest -q`; expected PASS. With the lab data staged: `pytest tests/regression -q`; expected PASS.

- [ ] **Step 5: Commit** — `git commit -m "A 'light' frame is provisional: bias by exposure time, else arc or science from its content before the trace stage"`.

---

### Task 7: Targets per fiber (§5.2.5)

**Files:**
- Modify: `src/mastpipe/core/model.py` (`FrameMeta.targets`, `target_for`), `src/mastpipe/io/headers.py` (`per_unit` option), `src/mastpipe/io/ingest.py` (`targets` resolution, manifest column `targets`), `src/mastpipe/steps/extract.py` (`target=frame.meta.target_for(trace.telescope)`), `src/mastpipe/io/products.py` (`_add_meta` / `_read_meta` serialise the tuple as JSON), `src/mastpipe/pipeline/runner.py` and `commands.py` (QC and `show` print the mapping), `src/mastpipe/sim/simulator.py` (`_cards` writes `OBJ{nn}` per telescope from the scene's targets; `OBJECT` stays as the joined list for readability), `src/mastpipe/instruments/highspec/header_map.yaml` (`targets: {per_unit: "OBJ{nn}"}` and `target: {keys: [OBJECT]}`)
- Test: `tests/unit/test_model.py`, `tests/unit/test_ingest.py`, `tests/unit/test_extract.py`, `tests/e2e/test_quicklook.py`

**Interfaces:**
- Produces: `FrameMeta.targets: tuple[tuple[str, str], ...] = ()` replacing `target`; `FrameMeta.target_for(telescope: str | None) -> str | None` (exact telescope match, else the `"*"` entry, else `None`); ingest fills `targets` from the `per_unit` field (`{nn}` expanded for `n` in `1..UNIT_COUNT`, keyed by `unit_name(n)`), else from a single `target` value as `(("*", value),)`; manifest column `targets` parses `mast01=HD 1234;mast02=HD 5678` (and a bare value as `*`).
- `FieldSpec.per_unit: str | None` (a template containing `{nn}`; mutually exclusive with `keys`).

- [ ] **Step 1: Write the tests**

```python
def test_target_for_prefers_exact_then_star():
    meta = FrameMeta(..., targets=(("*", "all"), ("mast02", "two")))
    assert meta.target_for("mast02") == "two" and meta.target_for("mast07") == "all" and FrameMeta(...).target_for("mast01") is None

def test_per_unit_cards_become_targets(tmp_path, highspec_map):  # OBJ01='HD 1', OBJ03='HD 3' → (("mast01","HD 1"),("mast03","HD 3"))
def test_single_object_applies_to_all(...):  # OBJECT only → (("*", value),)
def test_manifest_targets_column(...):  # "mast01=HD 1;mast02=HD 2"
def test_extract_assigns_target_per_telescope(...):
# tests/e2e: in test_quicklook_end_to_end_on_simulated_night assert len(group_by_target(spectra)) == 4 combined products per exposure
```

- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement** (the `_add_meta` card value for `targets` is `_json(list of pairs)`; `_read_meta` restores the tuple). **Step 4:** `pytest -q`; expected PASS, including four combined products per simulated exposure.

- [ ] **Step 5: Commit** — `git commit -m "Targets per fiber: per-unit header cards, a manifest column, and a single OBJECT as the fallback"`.

---

### Task 8: Per-frame isolation in the calibration stage (§5.3)

**Files:**
- Modify: `src/mastpipe/pipeline/runner.py` (`_masters`, new `_usable_frames`, `_run_group` binning)
- Test: `tests/e2e/test_quicklook.py`, `tests/unit/test_runner.py`

**Interfaces:**
- Produces: `_usable_frames(group, refs: list[FrameRef], calibrate: bool) -> list[Frame]`: loads (and calibrates when asked) one frame at a time; a frame that raises, whose `data.shape` differs from the first good frame's, or whose `meta.binning != group.binning` is appended to `group.summary.failed` as `(source, reason)` and dropped; `on_error == "raise"` re-raises. `_masters` uses it for every kind and for the flats; a kind with no survivor falls back to the library as today.
- `_run_group`: `group.binning` is the most common binning among the science frames (after Task 6's resolution), else among all refs.

- [ ] **Step 1: Write the tests**

```python
def test_quicklook_survives_bad_calibration_frames(tmp_path, night):
    # copy the night; add a truncated bias (20x30), a flat with HBIN=2, and a corrupt arc (garbage bytes named .fits)
    summary = run_quicklook(folder, out, options())
    assert {s for s, _ in summary.failed} == {"sim_bias_bad.fits", "sim_flat_bad.fits", "sim_arc_bad.fits"}
    assert len(spec_products(out)) == N_SCIENCE          # the night still yields science products
    assert main([...quicklook argv...]) == 2               # a dropped frame is a failure to the supervisor

def test_group_binning_follows_the_science_frames(...):  # unit, on a hand-made ref list
```

- [ ] **Step 2: Run** — expected FAIL (today the group fails). **Step 3: Implement.** **Step 4:** `pytest -q`; expected PASS.

- [ ] **Step 5: Commit** — `git commit -m "One bad calibration frame is dropped and recorded; the night is still reduced"`.

---

### Task 9: Atomic writes and the index lock (§5.3)

**Files:**
- Create: `src/mastpipe/io/atomic.py`
- Modify: `src/mastpipe/io/products.py` (`_write`, `write_qc`), `src/mastpipe/pipeline/callib.py` (`_save`, `_replace`, `add`, `add_reference`), `src/mastpipe/sim/simulator.py` (`save_truth`, optional)
- Test: `tests/unit/test_atomic.py`, `tests/unit/test_callib.py`

**Interfaces:**
- Produces: `atomic.atomic_write(path: Path) -> ContextManager[Path]` (yields `path.with_name(path.name + ".tmp-<pid>")`; on normal exit `os.replace(tmp, path)`; on exception the temporary is removed and the error re-raised); `atomic.locked(path: Path, timeout_s: float = 30.0) -> ContextManager[None]` (an exclusive advisory lock on `path`, created if missing, `fcntl.flock` on POSIX and `msvcrt.locking` on Windows, polled every 50 ms; `TimeoutError(f"could not lock {path} within {timeout_s} s")`).
- `CalibLibrary`: `LOCK_FILE = "index.lock"`; `_replace` runs `with locked(self.root / LOCK_FILE): self._entries = self._load(); <drop + append>; self._save()`; `_save` writes through `atomic_write`.

- [ ] **Step 1: Write the tests**

```python
def test_atomic_write_leaves_no_partial_file(tmp_path):  # raise inside the block → target absent, no .tmp- file left
def test_lock_times_out_with_a_message(tmp_path):
    with locked(tmp_path / "index.lock"):
        with pytest.raises(TimeoutError, match="index.lock"):
            with locked(tmp_path / "index.lock", timeout_s=0.2): pass
def test_two_processes_add_to_one_library(tmp_path):
    # multiprocessing (spawn) two workers each adding 25 distinct masters (product_name differs by checksum)
    assert len(CalibLibrary(tmp_path).list("bias")) == 50 and json.loads((tmp_path / "index.json").read_text())
```

- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement.** **Step 4:** `pytest tests/unit -q` then `pytest -q`; expected PASS on the platform at hand; CI proves Windows.

- [ ] **Step 5: Commit** — `git commit -m "Atomic product and index writes; the library index is updated under a lock"`.

---

### Task 10: The complete idempotence key and master reuse (§5.4)

**Files:**
- Modify: `src/mastpipe/pipeline/provenance.py` (`layout_hash` → `MPLAYH`, `masters_hash` → `MPMKH`), `src/mastpipe/pipeline/runner.py` (`_run_group` builds `effective` with the time policy as an override; `_Group.drift_params` returns `params.drift`; `_science_frame` key; `_stored_key`; `_masters` reuse), `src/mastpipe/pipeline/callib.py` (entries gain `inputs` and `file_checksum`; `find_by_inputs`), `src/mastpipe/instruments/registry.py` (`layout_hash(instrument, preset) -> str`, sha256 of the CSV bytes, 12 hex)
- Test: `tests/unit/test_provenance.py` (create if absent), `tests/unit/test_callib.py`, `tests/e2e/test_quicklook.py`

**Interfaces:**
- Produces: key = `(checksum, config_hash, reference_hash, tuple(sorted(drift_arcs)), layout, layout_hash, masters_hash, timezone)`; `masters_hash` = sha256 over the sorted file checksums of `group.masters` values (12 hex); `CalibLibrary.find_by_inputs(kind: str, inputs: Sequence[str], camera: str) -> tuple[MasterCalib, Path] | None` matches an entry whose `inputs` equals the sorted list.
- `_masters`: before building a master from `frames`, `library.find_by_inputs(kind, [f.meta.checksum ...], camera)`; reuse it unless `options.force`.

- [ ] **Step 1: Write the tests**

```python
def test_provenance_round_trips_new_cards():  # layout_hash, masters_hash, timezone survive to_header_cards/from_header; absent cards read as None
def test_find_by_inputs_matches_sorted_checksums(tmp_path): ...
def test_quicklook_rereduces_when_layout_reference_master_or_policy_changes(tmp_path, transferred_run, ...):
    # (a) unchanged rerun: reduced_again == []; (b) --time-policy nearest: every science frame reduced again;
    # (c) a byte appended to the layout CSV copy used via a temp instrument dir, or monkeypatched layout_hash: reduced again;
    # (d) the reference file rewritten with a different hash: reduced again; (e) --force rebuilds masters, otherwise they are reused (index unchanged, mtime unchanged)
```

- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement** (note: with the time policy in the effective config its hash changes, so the config hash covers it; keep `MPTZ` from Task 3 in the key). **Step 4:** `pytest -q`; expected PASS.

- [ ] **Step 5: Commit** — `git commit -m "Idempotence key covers reference hash, layout, masters, time policy and time zone; equal masters are reused"`.

---

### Task 11: Fibers keyed by the layout's trace id; a simulator that can tell fibers apart (§5.5)

**Files:**
- Modify: `src/mastpipe/steps/trace.py` (`match_layout`: `layout_index = int(layout["trace"].iloc[i])`, `missing` lists trace ids), `src/mastpipe/sim/simulator.py` (`_Scene.x_c`, `B` become per-trace arrays; `line_x` and `sky_x` shape `(n_traces, n_lines)`; new `simulate_night` arguments `dispersion_jitter_px: float = 0.3`, `drift_gradient_px: float = 0.05` (shift difference top to bottom of the detector), `swap_fibers: tuple[int, int] | None = None`; `SimTruth` gains `x_c_per_trace`, `B_per_trace`, `drift_gradient_px`), `tests/conftest.py` (`truth_traceset` unchanged), `tests/e2e/test_quicklook.py` (truth comparison uses per-trace `line_pixels`; tolerances re-derived)
- Test: `tests/unit/test_trace.py`, `tests/unit/test_simulator.py`, `tests/e2e/test_quicklook.py`

**Interfaces:**
- Produces: `Trace.layout_index` is the layout's `trace` id; the simulator's truth `line_pixels[i]` is trace `i`'s undrifted line columns; a frame's drift at trace `i` is `shift + drift_gradient_px · (y_i − ny/2) / ny`.

- [ ] **Step 1: Write the tests**

```python
def test_match_layout_keys_by_trace_column_not_row(survey_layout):
    shuffled = survey_layout.sample(frac=1, random_state=3).reset_index(drop=True)
    a = layout_indexed(match_layout(ts, survey_layout, 2.0)); b = layout_indexed(match_layout(ts, shuffled, 2.0))
    assert [t.index for t in a] == [t.index for t in b] and [t.telescope for t in a] == [t.telescope for t in b]

def test_simulated_traces_have_distinct_dispersions(sim_night): assert len(set(np.round(truth.x_c_per_trace, 3))) > 1
def test_swapped_fibers_are_detected(tmp_path, ...):
    # simulate_night(..., swap_fibers=(1, 3)); the e2e truth comparison helper reports map errors > MAX_MAP_ERROR_DISP for exactly those fibers
```

- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement.** **Step 4:** `pytest -q`; re-derive `MAX_MAP_ERROR_DISP` and the drift tolerances from the new truth only if they fail for a reason you can state in the test's docstring; expected PASS. With the lab data staged: `pytest tests/regression -q`; expected PASS unchanged.

- [ ] **Step 5: Commit** — `git commit -m "Key fibers by the layout's trace id; the simulator gives every fiber its own dispersion and a y-dependent drift"`.

---

### Task 12: Four bounded numerical fixes (§5.6)

**Files:**
- Modify: `src/mastpipe/core/model.py` (`CalibrationSet` variances; `DriftModel.n_usable`; `CombinedSpectrum.version_missing`), `src/mastpipe/steps/calib.py` (`build_calibration_set`, `calibrate_frame`), `src/mastpipe/steps/lines.py` (`find_lines(..., mask=None)`), `src/mastpipe/pipeline/identify.py` (`_arc_lines` passes `frame.mask`), `src/mastpipe/config/models.py` (`DriftParams.min_fibers: int = 3`), `src/mastpipe/steps/drift.py` (`fit_drift`), `src/mastpipe/pipeline/runner.py` (`_blend_drifts` carries `n_usable`; QC), `src/mastpipe/io/products.py` (`MPDUSE` comment names the reason; `MPVMISS`), `src/mastpipe/steps/combine.py`
- Test: `tests/unit/test_calib.py`, `tests/unit/test_drift.py`, `tests/unit/test_combine.py`, `tests/unit/test_lines.py`

**Interfaces:**
- `CalibrationSet` gains `bias_variance`, `background_variance`, `dark_rate_variance` (`np.ndarray | None`); `calibrate_frame` adds `bias_variance + (t / t_bkg)² · (background_variance + bias_variance)` or `bias_variance + t² · dark_rate_variance` to the pixel variance when present.
- `find_lines(image, variance, traceset, params, mask: np.ndarray | None = None)`: a line is `saturated` when `params.saturation` says so **or** when `mask` (the calibrated frame's mask, which `calibrate_frame` sets on raw pixels at or above the detector's saturation) marks any pixel of the line's window on its trace. `find_lines` works on bias-subtracted images, so a threshold in raw ADU cannot be applied to them; the mask carries the raw-level fact. `_arc_lines` passes the frame's mask.
- `fit_drift`: `usable = n_usable >= params.min_fibers and scatter <= params.max_scatter_px`; `DriftModel.n_usable: int`.
- `combine_spectra`: every version carried by **any** contributor is combined from the contributors that carry it; `CombinedSpectrum.version_missing: dict[str, list[int]]` (version → fibers lacking it); `write_combined` writes it as `MPVMISS` (JSON string); `read_combined` restores it.

- [ ] **Step 1: Write the tests**

```python
def test_master_variance_enters_the_error_budget():  # same frame, calib with and without bias_variance → variance differs by exactly bias_variance
def test_saturated_lines_are_flagged_through_the_mask(sim_full_arc):
    # a raw pixel of one line set to 65535 → calibrate_frame masks it → find_lines(..., mask=frame.mask) marks that line saturated, use_x False; without the mask it is not flagged
def test_drift_with_two_fibers_is_not_usable():  # per_fiber tables for fibers 0 and 1 only → usable False, n_usable == 2
def test_combine_keeps_a_version_some_contributors_lack():
    # three science spectra, two with skysub → "skysub" in versions, version_missing == {"skysub": [fiber of the third]}
```

- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement.** **Step 4:** `pytest -q`; expected PASS (the lab3 layout has three arc fibers, the survey layout five).

- [ ] **Step 5: Commit** — `git commit -m "Master variances in the error budget, saturated lines flagged through the mask, a drift gate that needs three fibers, versions never dropped in combination"`.

---

### Task 13: The readout key includes the amplifier (§5.7)

**Files:**
- Modify: `src/mastpipe/instruments/highspec/header_map.yaml` (`amplifier` field), `src/mastpipe/instruments/highspec/config.yaml` (`readouts:`), `src/mastpipe/config/models.py` (`ReadoutNoise`, `DetectorConfig.readouts`), `src/mastpipe/core/model.py` (`FrameMeta.amplifier`, `readout_key` property), `src/mastpipe/io/ingest.py` (`PASSTHROUGH_FIELDS`, `MANIFEST_FIELDS`), `src/mastpipe/pipeline/callib.py` (`find(..., readout_key)`; entries store `readout_key`), `src/mastpipe/pipeline/runner.py` (`_Group.readout_key`; `build_calibration_set(..., amplifier=...)`), `src/mastpipe/steps/calib.py`, `src/mastpipe/sim/simulator.py` (`OUTPTAMP: "Conventional"`)
- Test: `tests/unit/test_headers.py`, `tests/unit/test_model.py`, `tests/unit/test_callib.py`, `tests/unit/test_calib.py`

**Interfaces:**
- `amplifier: {keys: [OUTPTAMP], lower: true, map: {conventional: conventional, "electron multiplying": em}}`; `FrameMeta.amplifier: str | None`; `FrameMeta.readout_key -> str` = `f"{readout_mode or 'unknown'}/{amplifier or 'unknown'}"`.
- `ReadoutNoise(_Strict)`: `gain: float`, `read_noise: float`; `DetectorConfig.readouts: dict[str, ReadoutNoise] = {}` keyed by amplifier; `build_calibration_set(bias, dark, background, detector, amplifier: str | None = None)` takes gain and read noise from `detector.readouts[amplifier]` when present, else the detector defaults.
- HighSpec `config.yaml`: `readouts: {conventional: {gain: 1.0, read_noise: 5.0}, em: {gain: 1.0, read_noise: 5.0}}` with the `PLACEHOLDER` comments kept.

- [ ] **Step 1: Write the tests** — `test_amplifier_is_read_and_keyed` (`OUTPTAMP='Electron Multiplying'` → `amplifier == "em"`, `readout_key == "image/em"`); `test_callib_find_separates_amplifiers` (a bias indexed under `image/conventional` is not found for `image/em`); `test_calibration_set_uses_the_amplifier_noise`.
- [ ] **Step 2: Run** — expected FAIL. **Step 3: Implement.** **Step 4:** `pytest -q`; expected PASS. **Step 5: Commit** — `git commit -m "Readout key includes the amplifier; per-amplifier gain and read noise"`.

---

### Task 14: Fleet lint and CI layout (§5.8)

**Files:**
- Modify: `pyproject.toml` (`[tool.ruff]`, `dev` extra pins `ruff==0.16.0`), `.github/workflows/ci.yml` (separate `lint` job, `concurrency`, format check first), source files ruff reports, `CLAUDE.md` (only after approval)

- [ ] **Step 1: Set the configuration**: `line-length = 125`, `target-version = "py311"`, `select = ["F","E","W","I","N","B","UP","C90","SIM","TID","NPY"]`, `[tool.ruff.lint.per-file-ignores] "src/mastpipe/steps/*" = ["N806", "N803"]` with a comment (physics-named variables), `[tool.ruff.lint.mccabe] max-complexity = <the worst current function, found by running ruff; ratchet comment as in MAST_common>`, `[tool.ruff.format]` `quote-style = "double"`, `docstring-code-format = true`.
- [ ] **Step 2: Run `ruff format . && ruff check --fix .`**, then fix what remains by hand (no new `noqa` without a one-line reason). `pytest tests/unit -q` must still pass.
- [ ] **Step 3: Restructure CI**: `concurrency: {group: ci-${{ github.ref }}, cancel-in-progress: true}`; a `lint` job on `ubuntu-latest` running `ruff format --check .` then `ruff check .`; remove the lint steps from the matrix job; keep `privacy_scan` in the matrix job; keep the `docs` and `fleet` jobs.
- [ ] **Step 4: Ask the lead** to approve a `CLAUDE.md` diff that references `mast-claude-config` (its `code-validation-and-ci-guidelines.md`) and states the deliberate departures (three platforms, Python 3.11 kept). Apply only on approval.
- [ ] **Step 5: Commit** — `git commit -m "Fleet ruff rules and CI layout: separate lint job, pinned ruff, concurrency group"`.

---

### Task 15: Documents (D3, §5.9)

**Files:**
- Modify: `.gitignore` (remove `docs/design/`), `README.md`, `docs/running/night-folder.md`, `docs/running/quicklook.md`, `docs/extending/adding-a-grating.md`, `docs/extending/simulator-and-tests.md`, `DECISIONS.md`, `CHANGES.md`, `CLAUDE.md` (the `ext/pypeit_adapter.py` sentence; approval needed)
- Scan: every file under `docs/design/`

- [ ] **Step 1: Scrub the design documents.** With the name list the lead supplies: `git grep -n -i -E "<names>" docs/design` must return nothing after replacing each occurrence with a role ("the HighSpec laboratory reduction notebook", "the DeepSpec instrument team", "the project lead", "a MAST team member"). Then `python scripts/privacy_scan.py` extended to `docs/design` must pass (edit its path list if needed).
- [ ] **Step 2: Remove `docs/design/` from `.gitignore`** and `git add docs/design`.
- [ ] **Step 3: Rewrite `docs/running/night-folder.md`**: the `MAST_spec` folder shape (`<night>/highspec/acquisition-NNNN/exposure-NNN.fits`), recursive scanning and relative sources, the `INSTRUME` gate, `light` and its resolution (bias rule, content rule, manifest override), the required `--configuration Mg` when no card names the grating, the `timezone` option and `MPTZ`, manifest sources as paths, globs or unique names, `targets` per unit. Update `quicklook.md` (`--configuration`), `adding-a-grating.md` (names come from `MAST_common`; the registry refuses others), `simulator-and-tests.md` (per-fiber dispersion, `swap_fibers`; Task 16 adds `site_layout`).
- [ ] **Step 4: README and CLAUDE.md**: replace the adapter sentences with "An optional PypeIt adapter is planned; the `pypeit` extra pre-installs its dependency. No adapter module exists yet." (CLAUDE.md edit → ask the lead).
- [ ] **Step 5: `DECISIONS.md`**: one dated section `## 2026-11-04 — Plan 1.5: hardening at the boundaries` with one bullet per decision D1–D4 and per designer call C1–C7 (the design document has the wording). **`CHANGES.md`**: a bullet list under the unreleased version covering Tasks 1–16.
- [ ] **Step 6: Build the docs** — `sphinx-build -W --keep-going -b html docs docs/_build/html` from the `.venv` that has `myst-nb`; expected: no warnings; the notebooks execute with the new names.
- [ ] **Step 7: Commit** — `git commit -m "Publish the design documents; night-folder guide for the site contract; decisions and changelog for Plan 1.5"`.

---

### Task 16: The site-shaped night and the laboratory manifest (acceptance 1 and 2)

**Files:**
- Modify: `src/mastpipe/sim/simulator.py` (`simulate_night(..., site_layout: bool = False)`), `tests/e2e/test_quicklook.py`, `tests/regression/test_lab_highspec.py`, `docs/extending/simulator-and-tests.md`

**Interfaces:**
- `site_layout=True` writes `highspec/acquisition-NNNN/exposure-NNN.fits` (one acquisition folder per schedule entry, numbered from 1, with a `seq.txt` in each parent as `MAST_spec` leaves) and **site cards only**: `IMAGETYP` from `{"bias": "bias", "dark": "dark", "flat": "flat", "arc": "light", "science": "light", "background": "light"}`, `EXPOSURE`, `EXPREQ`, `FRAME` (start, with `.000`), `DATE` (same instant), `CCD-TEMP`, `HBIN`, `VBIN`, `READMODE`, `OUTPTAMP`, `INSTRUME = "Highspec"`, `PROGRAM = "MAST"`; no `GRATING`, `LAMP`, `OBJECT`, `OBJnn`, `FIBPRE`, `SHUTTER`, `DATE-OBS`. The truth table's `file` column holds the relative path.

- [ ] **Step 1: Write the tests**

```python
def test_site_shaped_night_reduces_without_a_manifest(tmp_path, highspec_effective, survey_layout, mg_lines, identified_library):
    truth = simulate_night(..., out_dir=site, seed=0, n_science=4, site_layout=True)
    summary = run_quicklook(site, out, RunOptions("highspec", configuration="Mg", layout="survey45", library=identified_library))
    assert not summary.failed and len(spec_products(out)) == 4
    # spectra equal to the labelled night's within noise: same fibers, wavelength status "transferred", flux arrays allclose
    # the background (300 s, lamp off) is "unknown" by content and listed in skipped

def test_site_night_without_configuration_stops_with_a_message(tmp_path, ...):
    assert main(["quicklook", str(site), "--out", ..., "--instrument", "highspec"]) == 2  # output names the valid grating names
```

Regression (`tests/regression/test_lab_highspec.py`): the manifest keeps `exposure_type` for `flat_518.fits` (`flat`) and `background.fits` (`background`) and `image_bias.fits` (its header already says `bias`; keep the row for the configuration), sets `configuration` to `Mg` for every row, and gives the focus frames **no** `exposure_type` and no `arc_kind`; `test_lab_traces_match_notebook`'s assertion on the full arc becomes `exposure_type == "arc"` with `classification_rule == "content:arc_like"` after the runner's resolution (assert on the science QC / summary rather than on `scan_folder` alone).

- [ ] **Step 2: Run** — `pytest tests/e2e -q -k site`; expected FAIL. **Step 3: Implement `site_layout`.** **Step 4:** `pytest -q`; expected PASS. With the lab data staged: `pytest tests/regression -q`; expected PASS with the numbers unchanged.

- [ ] **Step 5: Commit** — `git commit -m "A simulated night in MAST_spec's shape reduces with no manifest; the lab manifest types only what content cannot"`.

---

## Closing

After Task 16: `ruff format --check . && ruff check . && python scripts/privacy_scan.py && pytest -q`, the docs build, and the regression suite with the data staged. Then the whole-branch review (subagent-driven development's final step). Nothing is pushed without the lead's per-push approval.
