# mastpipe Plan 1.5 — hardening at the boundaries (design)

Written 2026-11-04 (Claude Fable 5.1, with the project lead) from the first outside
review of the public repository (2026-11-02, at `main` @ `90ad1aa`) and a reading of
the MAST control-software repositories it was reviewed against (`MAST_common`
`a78cfe7`, `MAST_spec` `f0eaae4`, `MAST_control`, `mast-claude-config` `8f522d3`).
No local paths, no individuals named.

## 1. Intent

Plan 1 produced a HighSpec pipeline that is correct on the laboratory frames and on
the simulator. The review's verdict was that the algorithms are not the problem: the
work between this version and a real night is **at the boundaries** — what the
control software writes, how the night folder is shaped, which clock stamped the
frames, and what happens when one frame is bad or two processes write at once.

Plan 1.5 closes the boundary gaps that are **mastpipe's to close**, adopts the
fleet's shared vocabulary so the two sides can meet, and makes the unattended run
robust. It does not wait for the header contract with `MAST_spec` (that is a team
conversation, see §9) but prepares every seam for it.

Evidence that this is needed now, not later: the laboratory Th-Ar frames themselves
were taken through `MAST_spec`'s Newton driver and carry `IMAGETYP = 'light'`, no
`OBJECT`, `GRATING` or `LAMP` card, and a naive `DATE` with no time zone. The lab
regression works only because of its manifest.

## 2. Goals, non-goals, anti-goals

**Goals**

1. A night folder shaped the way `MAST_spec` writes it reduces **without a manifest**,
   given only the grating on the command line.
2. mastpipe speaks `MAST_common`'s vocabulary: frame types, grating names, instrument
   names, lamp names, unit names.
3. Quick-look survives one bad frame of any kind, an interrupted write, and a
   concurrent `identify`.
4. A product that is stale is reduced again; a product that is current is skipped.
   The idempotence key covers everything the product depends on.
5. The test suite can fail on a fiber-identity bug.
6. The repository follows the fleet's CI and lint conventions where they do not
   conflict with being a pip-installable, cross-platform package.
7. The repository's documents are readable by anyone who reads `DECISIONS.md`.

**Non-goals** (Plan 2 or the header contract)

- DeepSpec; Horne extraction; flux and telluric calibration beyond pass-through.
- The full-arc **campaign** (about twenty partial frames from the portable lamp unit):
  dormant until the unit exists. `solve_full_arc` keeps assuming one frame.
- The physical slit order of the survey head (which fiber is which unit): needs the
  instrument team. Only the naming of the telescope column changes here.
- An integration hook (`AssignmentNotification`, the per-plan products folder).
- Measured gain and read noise.

**Anti-goals** (what would make this a failure even if it worked)

- Any change to the numbers the lab regression reproduces (trace centres, per-fiber
  rms, line counts). The algorithms are not touched.
- A runtime dependency on MongoDB, the observatory share, or a sibling clone.
  `pip install mastpipe` on a laptop must keep working.
- A silent default anywhere a frame's metadata is unknown.
- Reintroducing a person's name, a local path or a hostname into a tracked file.

## 3. Constraints

- Same model throughout (planning, implementation, review) per the vault rule.
- Repository rules in `CLAUDE.md`: no names, no paths, no secrets; `CLAUDE.md`
  itself changes only with the lead's explicit approval.
- Windows, Linux and macOS; `pathlib` only; `filterwarnings = error`.
- Test data only through `MASTPIPE_TEST_DATA` / `MASTPIPE_ATLAS_DIR`.
- CI must stay green on every push; the docs build treats warnings as errors.

## 4. Decisions taken with the lead (2026-11-04)

| # | Decision | Chosen |
|---|---|---|
| D1 | Canonical identifiers | **Adopt `MAST_common`'s spellings**: gratings `Mg`, `Halpha`, `Ca` (and `Future`, `Empty` as known-but-unsupported); instruments `highspec`, `deepspec`; telescopes `mast01`…`mast20`. `mgb518`, `halpha`, `caiihk` survive only as file-name aliases. |
| D2 | Binding to `MAST_common` | **Mirror + fleet-style CI.** `mastpipe.vocab` restates the literals and is what the code uses. A CI job checks `MAST_common` out as a sibling (paired branch, else `master`, as `MAST_control` does) and a test fails when the mirror diverges. No runtime import. |
| D3 | Design documents | **Publish `docs/design/`** once the privacy scan passes and every personal name is scrubbed. `DECISIONS.md`'s citations then resolve. |
| D4 | Arcs that the header calls `light` | **Detect from content, as a classification rule.** `light` becomes a provisional type that the runner refines to `arc` or `science` from the frame; a manifest still overrides. |

Calls made by the designer (reviewable, not agreed in so many words):

- C1 Grating ids are case-sensitive and spelled exactly as `MAST_common` does; the
  header map stops lower-casing `configuration`.
- C2 The time-zone question stays open; mastpipe makes its assumption **explicit and
  recorded** rather than guessing (§5.2.4).
- C3 Python 3.11 support stays (the fleet pins 3.12 for services; mastpipe is a
  library). Ruff's rule set, line length and pin follow the fleet.
- C4 `FIBPRE` leaves the header map: in MAST a "fiber preset" is the stage position
  (`deepspec` / `highspec`), not a fiber layout. The layout comes from `--layout`
  or the instrument default only.
- C5 Masters enter the idempotence key by checksum; they are no longer rebuilt and
  overwritten on every run when an equal master exists (§5.4).
- C6 Numerical items are limited to the four that are a few lines each and have a
  clear right answer (§5.6). The cosmic-filter question goes to Plan 2.
- C7 (after the builder check) Arc detection is frame-level and runs before the trace
  stage (§5.2.3); product names gain a short checksum suffix (§5.2.1); ingest gates
  on `INSTRUME` (§5.2.1); a dropped calibration frame is a recorded failure with exit
  code 2 (§5.3); the fleet test compares against `common.spec` only, with the config
  lists read as text (§5.1); the lab regression keeps its manifest rows for the flat
  and the background (§6).

## 5. The design

### 5.1 Vocabulary (`mastpipe/vocab.py`)

One module, importable with no dependencies beyond the standard library, restating
what `MAST_common/spec.py` and `config/specs.py` define:

- `FrameType`: `light`, `bias`, `dark`, `flat` — the four the control software writes.
- `GratingName`: `Ca`, `Halpha`, `Mg`, `Future`; `Disperser` adds `Empty`.
- `SpecInstrument`: `deepspec`, `highspec`.
- `LampName`: `ThAr`, `qTh` (the calibration-unit lamps / filter wheels).
- `DeepspecBand`: `I`, `G`, `R`, `U` (for Plan 2; carried now so the mirror is whole).
- `unit_name(number, project="mast") -> "mast01"`, the fleet's unit hostname rule.

mastpipe's own extensions live next to them, clearly labelled as pipeline-only:
`EXPOSURE_TYPES` keeps `arc`, `background`, `unknown` *and gains `light`* (§5.2.3).
Everything that today compares against a literal string for an instrument, grating,
frame type or lamp goes through `vocab`.

A test, `tests/fleet/test_vocab_matches_common.py`, imports `common.spec` (whose
only third-party import is Pydantic) when a `common` package is importable and
asserts member for member that the mirror equals `FrameType`, `Disperser`,
`GratingNames`, `SpecInstruments`, `WheelNames` and `DeepspecBands`. The grating and
lamp lists that `common/config/specs.py` keeps as locals inside its validator are
read from the file's text (AST or regex), not imported, because importing
`common.config` reaches for MongoDB. `unit_name`'s output is checked against
`common.canonical.canonic_unit_name`, the fleet's validator. The `fleet` job
installs mastpipe and nothing of the fleet's pinned environment. It **fails**, not skips, when the
environment variable `MASTPIPE_FLEET_CHECK=1` is set and `common` is not importable
(the fleet guideline's "a green run must mean something"). Locally without the
variable it skips with a message naming the sibling layout.

CI gains a job `fleet` (Linux, Python 3.12): checks this repo out as `mastpipe/` and
`MAST_common` as `common/` next to it, on a branch of the same name as the one under
test if `MAST_common` has one, else `master`; sets `PYTHONPATH` to the workspace root
and `MASTPIPE_FLEET_CHECK=1`; runs the fleet tests only. A vocabulary change in
`MAST_common` turns this job red without anything in mastpipe moving. The job needs
no secrets: both repositories are public.

### 5.2 Ingest contract

#### 5.2.1 Folder shape

`scan_folder` becomes **recursive** (`*.fits`, `*.fit`, `*.fits.fz`, case-insensitive),
skipping hidden entries, `seq.txt`, and any file that is not a FITS image (recorded
as skipped with its reason, never a failure). A frame's `source` is its path
**relative to the night folder**, with forward slashes, so
`acquisition-0012/exposure-003.fits` and `acquisition-0013/exposure-003.fits` are two
sources. Product names gain a short input-checksum suffix (six hex characters) after the
time stamp, so two frames with the same mid-time to the second, or two frames
without a time, never share a name. `identify`'s master lookup matches on the
prefix and tolerates the suffix; the tutorials are executed at docs build and
pick the new names up.

The manifest's `source` column accepts a relative path or a glob pattern over
relative paths (`acquisition-0012/*.fits`); a bare file name still matches when it is
unique in the tree (the laboratory manifests keep working) and is an error when it
is ambiguous. Rows that match nothing are reported, not dropped silently.

The header map gains an `instrument` field (`INSTRUME`, lower-cased: `MAST_spec`
writes `Highspec` on every spectrograph frame; the unit cameras write their own
names). Ingest skips, with a reason, any frame whose instrument is stated and is
not the one requested; a frame with no instrument card is accepted (the laboratory
frames have none). The operator normally points quick-look at `<night>/highspec/`;
pointing at `<night>/` works because the gate drops the unit cameras' and DeepSpec's
frames. Where `MAST_spec` will file HighSpec biases, darks and flats is not yet
defined (its HighSpec path takes only `light` today); the simulator writes them in
acquisition folders of their own under the same root, which is a guess and is
labelled as one.

#### 5.2.2 Grating

With no `GRATING` card and no file-name alias, the configuration is **unresolved**.
Quick-look then requires `--configuration <GratingName>`; without it the run stops
before reducing anything with a message listing the valid names. There is no silent
default. The file-name regex maps the laboratory aliases to the canonical names
through the map option (`mgb518 → Mg`, `halpha → Halpha`, `caiihk → Ca`).

The grating block files are renamed to the canonical names (`gratings/Mg.yaml`) and
the registry refuses a block whose name is not a `GratingName`.

#### 5.2.3 Frame type

The header map maps `light` to **`light`**, a provisional exposure type meaning
"shutter open, source not stated". Ingest keeps classification rules for frames
whose header has no valid type (bias by exposure time, dark by lamp + shutter, arc by
lamp card once one exists), and one rule now applies to `light` as well: a `light`
frame with an exposure time below the bias threshold is a bias. The runner resolves
the remaining `light` frames **before the trace stage**, because on a real night the
full arc is itself `light` and is the trace source when the night has no flat. The
statistic is therefore **frame-level**, on the frame collapsed along the spatial axis
(a 1-D spectrum): an emission-line-dominated spectrum (most of the flux above a
robust continuum sits in narrow peaks) makes the frame `arc`; a continuum-dominated
one makes it `science`; a frame with no signal above the noise is `unknown`. The
exact statistic and thresholds are the builder's; the contract is that every
laboratory Th-Ar focus frame and every simulated arc is `arc`, every simulated
science frame is `science`, the laboratory flat is `science` (a continuum, which
only a manifest or a `flat` header type can tell apart) and the laboratory lamp-off
background is `unknown`. `classify_arc` then says full or calibration-fiber as
today, once the traces exist.

The resolution is recorded in `classification_rule` (`content:arc_like`,
`content:continuum`) and in the product's provenance, and a manifest row's
`exposure_type` still wins over it. A frame the header calls `bias`, `dark` or
`flat` is never re-classified.

#### 5.2.4 Time

Three changes:

- The Newton's `FRAME` card ("Start of Frame Exposure", with milliseconds) is tried
  before `DATE` ("file creation date"). `DATE-OBS` stays first for frames that have
  it. DeepSpec's `T_MID` (UTC, mid-exposure) is declared for Plan 2.
- Every time field in a header map carries a `timezone` option: `utc`, an IANA zone
  name, or `assume-utc`. A stamp with an explicit offset is honoured whatever the
  option says. A naive stamp is interpreted in the named zone and converted to UTC;
  under `assume-utc` it is taken as UTC **and** the run logs one warning per run and
  every product carries `MPTZ = 'assumed-utc'`. HighSpec ships with `assume-utc`
  until the control-software clock is settled (one frame with a known UTC time
  settles it; then the option flips to `Asia/Jerusalem` or `utc` in one line).
- All ISO parsing goes through one helper in `mastpipe.core.times` (there are two
  copies today); it always returns an aware UTC instant. Formatted `mid_time`
  strings stay naive ISO UTC (no `Z` suffix), as in Plan 1, so product names,
  identifiers and stored values do not change form (settled during execution,
  2026-11-04). `make_identifier` normalises offsets to UTC.

#### 5.2.5 Targets per fiber

`FrameMeta.target: str | None` becomes `FrameMeta.targets: Mapping[str, str]`
keyed by telescope name (`mast01 → "HD 1234"`). Sources, in order: a header-map
field spec of a new kind, `per_telescope`, with a keyword template (for example
`OBJ{nn}` or `TARG{nn}`; none is written today, so the template is a seam, not a
promise); a manifest column `targets` written `mast01=HD 1234;mast02=…`; and a
single `OBJECT` card, which applies to every science fiber (the laboratory and
pre-contract case). The combination step groups by target across fibers as before
and now actually sees more than one group on a simulated survey night.

The simulator writes per-fiber targets through the template so the end-to-end
test exercises several target groups per exposure.

#### 5.2.6 Telescope column

`survey45.csv`'s `telescope` column becomes `mast01`…`mast20`; `TraceSet`, products
and QC use the same names. (Which physical trace is `mast01` is the instrument
team's answer and is out of scope; the file's order is a placeholder and says so.)

### 5.3 Runner robustness

- **Per-frame isolation for calibrations.** Master building loads frames one at a
  time; a frame that fails to load or has the wrong shape or binning is recorded in
  `failed` and dropped. The master is built from the survivors; only when none
  survive does that kind count as missing, and then the library is tried as today.
  A group fails only when it has no usable traces at all. A dropped frame is a
  **failure** of that frame (`failed`, and the command's exit code 2, which is what
  a supervisor should see: the night was reduced and something was wrong); a file
  that is not a FITS image is `skipped`. A group's binning is that of its science
  frames (the most common value); a calibration frame with another binning is
  dropped from the group, not the other way round. The end-to-end test plants a
  truncated bias, an odd-binned flat and a corrupt arc and asserts the night still
  yields science products and reports three failures.
- **Atomic writes.** Every product, master, reference, index and QC file is written
  to a temporary name in the same directory and moved into place with `os.replace`.
  A reader never sees a partial file.
- **Index lock.** `CalibrationLibrary` takes an exclusive lock on `index.lock` for
  the duration of a read-modify-write (a small cross-platform helper over `fcntl` /
  `msvcrt`, with a timeout and a clear error). The index is re-read under the lock
  before it is modified, so a concurrent `identify` and `quicklook` both keep their
  entries. A test runs two writers in processes against one library.

### 5.4 Idempotence key

The key of a science product becomes: the frame's checksum; the configuration hash
(which now includes `time_policy` and every `--set` override); the **hash of the
reference solution file**; the hashes of the drift-arc frames; the hash of the layout
CSV; the checksums of the masters used; the time-zone option. `MPREFH`, `MPLAYOUT`
plus a new `MPLAYH` (layout hash) and `MPMKH` (masters hash) carry the parts. A
product with the old key shape is simply reduced again once.

Masters: when the library already holds a master of the same kind, camera, binning,
readout key, configuration and **input checksums**, it is reused rather than rebuilt
and overwritten; `--force` rebuilds.

### 5.5 Fiber keying and the simulator

- `layout_index` is the layout's `trace` column, not the row position; `TraceSet`,
  solutions, drift models and sky pairing key on it everywhere (this is what
  `DECISIONS.md` already claims). The shipped CSVs are unaffected; a test shuffles
  the rows of a layout and asserts identical products.
- The simulator gives every fiber its own small perturbation of the dispersion
  (`x_c`, `B` jitter per trace, seeded) and a drift that varies linearly with `y`, so
  that a swapped or shifted fiber identity changes the truth comparison. End-to-end
  tolerances are re-derived from the new truth; the lab regression is untouched.

### 5.6 Numerics (bounded)

1. Master-frame variances propagate into the calibrated frame's variance.
2. Lines flagged `saturated` are excluded from the wavelength fit (they are already
   excluded from `use_x` in the line table; the fit must honour the flag).
3. The drift model records `n_usable`; with fewer than three usable arc fibers the
   scatter gate cannot judge, so the model is marked `usable` only if `n_usable >= 3`
   or the configured `min_fibers` says otherwise, and the QC report and `MPDUSE`
   say why.
4. A combined spectrum keeps every version that **all** contributors carry; a
   version missing on some contributor is combined from the contributors that have
   it and recorded in provenance (`MPVMISS` lists the fibers without it), never
   dropped silently.

### 5.7 Readout key

`readout_mode` becomes a composite of `READMODE` and `OUTPTAMP` (`image/conventional`,
`image/em`); the library matches on it. The detector config gains an optional
`readouts:` table of gain and read noise per amplifier, falling back to the camera's
defaults (still placeholders, still said so).

### 5.8 Fleet alignment

- `ruff`: the fleet's `select` (`F,E,W,I,N,B,UP,C90,SIM,TID`) plus `NPY`; line length
  125; `ruff==0.16.0` pinned in the `dev` extra; `format --check` before `check`;
  `N806`/`N803` allowed for physics-named variables by a documented per-file ignore
  rather than a blanket one; `max-complexity` set to the worst current function and
  ratcheted down as functions are split.
- CI: a separate `lint` job on Linux; the test matrix keeps three platforms × 3.11,
  3.12; `concurrency` group per ref; the `fleet` job of §5.1; the privacy scan stays
  in the test job.
- `CLAUDE.md`: references `mast-claude-config` and its guidelines (edit needs the
  lead's approval at implementation time).

### 5.9 Documents

- `docs/design/` leaves `.gitignore` and is committed, after the privacy scan and a
  name scrub (the laboratory notebook, the DeepSpec instrument team, the project lead
  are roles, not names). Documents that reference each other by file name keep
  working.
- `CLAUDE.md` and `README.md` stop naming `ext/pypeit_adapter.py` as if it existed:
  "an optional adapter is planned; the `pypeit` extra pre-installs its dependency".
- `docs/running/night-folder.md` describes the real folder shape, the provisional
  `light` type and its resolution, the required `--configuration`, the time-zone
  option, the manifest patterns and the per-fiber targets.
- `DECISIONS.md` gets one entry per decision in §4; `CHANGES.md` one section for the
  release.

## 6. Acceptance

1. A simulated night written in `MAST_spec`'s shape — `acquisition-NNNN/exposure-NNN.fits`,
   `IMAGETYP = 'light'` on arcs and science alike, no `GRATING`, `OBJECT`, `LAMP` or
   `FIBPRE` — reduces with `mastpipe quicklook night/ --out out --instrument highspec
   --configuration Mg --layout survey45` and no manifest, producing the same science
   products as the fully-labelled night.
2. The laboratory regression passes unchanged in its numbers, staging the same four
   calibration files and focus scan as today, with a manifest that says `Mg`, keeps
   `exposure_type` only for the flat and the lamp-off background (both `light` in
   their headers, both continuum-or-nothing by content) and gives the Th-Ar focus
   frames no type: content detection classifies them.
3. The bad-calibration-frame night yields science products.
4. Two processes writing the same library leave a valid index holding both entries.
5. Editing the layout CSV, the reference file, a master or `--time-policy` causes a
   re-reduction; re-running unchanged reduces nothing.
6. Shuffling layout rows changes no product; a deliberately swapped fiber in the
   simulator fails the end-to-end test.
7. The `fleet` job is green against `MAST_common` `master`; changing one literal in
   `mastpipe.vocab` turns it red.
8. `ruff format --check` and `ruff check` pass under the fleet configuration; CI has
   `lint`, `test`, `fleet` and `docs` jobs.
9. `DECISIONS.md`'s citations open on GitHub.

## 7. Left to the builder

- The frame-level arc-likeness statistic (§5.2.3) and its thresholds, as long as the
  contract holds on the lab frames and the simulator.
- The lock helper's implementation and timeout.
- How the per-telescope header template is spelled in the header map.
- Where the `fleet` tests live and how `MASTPIPE_FLEET_CHECK` is wired.
- Task order and splitting (the implementation plan).

## 8. Risks

- Renaming gratings and telescopes touches configuration, tests, tutorials and
  docs; a missed `mgb518` in a notebook breaks the docs build (which is the point of
  `-W`).
- The content-based arc rule is a heuristic; its failure mode is a science frame
  with strong emission lines on every fiber. Only `light` frames are subject to it,
  the manifest overrides it, and the header contract will retire it.
- Changing the simulator's truth may expose tolerances that were tuned, not derived.

## 9. Open with the team (not blocking this plan)

- The header contract: `IMAGETYP` that distinguishes an arc (or a `LAMP` card with
  the lamp name), `GRATING` with the `MAST_common` name, per-unit target cards, a UTC
  mid-time card like DeepSpec's `T_MID`, the amplifier. All of it is known to
  `MAST_spec` at exposure time.
- Whether the Newton's `DATE` / `FRAME` stamps are local or UTC.
- Whether "science exposures carry no arc light" stands against the paper's
  continuous calibration trace.
- Who records the slit order, and where.
- A pip-installable vocabulary package extracted from `MAST_common`, so that the
  mirror of §5.1 can go away.
