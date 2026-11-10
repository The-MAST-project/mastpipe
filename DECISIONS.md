# Decisions

A dated log of design and workflow decisions for mastpipe, newest last. Each
entry records what was decided and, briefly, why; reversals are recorded as new
entries rather than by editing old ones.

## 2026-09-24 — Independent package

- **Build mastpipe as an independent package** seeded from the team's
  HighSpec laboratory reduction notebook, with reduction steps as pure functions on arrays
  and a declarative instrument layer.
- **soxspipe rejected** as a base: echelle-order data model with no fiber
  concept, wavelength solution tied to ESO pinhole frames and predicted-line
  tables, hard-coded instrument branches, frame classification bound to ESO
  header keywords, effectively single-maintainer and GPL-3.
- **PypeIt deferred**: a strong community project whose new `Fiber` pipeline
  has the right shape, but it was unreleased, still carried instrument-specific
  logic and had no concept of arc fibers exposed alongside science fibers. The
  design keeps a future integration open through a single optional, lazily
  imported adapter module; nothing in the first plan depends on it.

## 2026-10-29 — Revision after the team meeting

Recorded from the design revision written after the team meeting with the MAST
group.

- **Reduction engine confirmed**: the independent package seeded from the laboratory
  notebook stays; PypeIt remains a future option behind the optional adapter.
- **Three-cadence wavelength scheme** replaces the daily full arc: a full arc
  every few months fits an independent solution per fiber and becomes the
  stored, versioned reference; a calibration-fiber arc one to three times per
  night (and at every instrument switch) measures a linear drift in pixel space
  against the reference; science exposures carry no arc light and take their
  wavelengths from the reference corrected by the transferred drift. HighSpec
  uses the grating formula, DeepSpec the dispersion equation, behind one
  interface.
- **HighSpec first, with the whole path in Plan 1**: reduced HighSpec spectra
  are needed now, so Plan 1 covers the full HighSpec path to sky-subtracted,
  wavelength-calibrated spectra and within-exposure combination. DeepSpec moves
  to Plan 2, with the October 2026 laboratory package as its reference data.
- **Six stages**: clean, trace and extract, wavelength, background (sky-fiber)
  subtraction, flux calibration, telluric correction; combination follows per
  target. Flux calibration and telluric correction exist as interfaces with
  pass-through implementations until their methods are chosen.
- **Kept versions**: every individual spectrum is kept, with its
  pre-subtraction and pre-calibration versions alongside, as columns of the
  per-fiber spectrum table on one wavelength grid (not separate extensions), so
  the table keeps one row per fiber.
- **Calls made while writing the revision** (each reversible by the user):
  1. Plan 1 scope is the full HighSpec path through sky subtraction and
     within-exposure combination; DeepSpec is Plan 2.
  2. The first HighSpec reference solution is reproduced from the two stored
     anchor pairs through the `identify` command, which is part of Plan 1
     because the survey layout and the other gratings need it anyway.
  3. Without a flat, traces come from any fully lit frame and the relative
     throughput defaults to one with a warning.
  4. Instrument switches need not be observable in the data: drift pairing is
     nearest in time within the configuration, bounded by `drift.max_age`.
  5. Reference validity is time-keyed: an exposure uses the reference valid at
     its mid-time, with an explicit override, so re-reductions reproduce.
  6. Quick-look re-reduces `reference`-status exposures when the night's first
     drift arc arrives; the product key includes the drift-arc set.
  7. Quick-look runs stages 1–4 and within-exposure combination; stages 5–6
     run only when their calibrations exist.
  8. Kept versions are columns in the per-fiber spectrum table on one
     wavelength grid, not separate extensions.

## 2026-10-29 — Instrument definitions as package data

- **One folder per instrument, discovered by the registry**: `config.yaml`
  carries everything except the configuration blocks, which are one YAML file
  each under `gratings/` (HighSpec) or `bands/` (DeepSpec) with `name` equal to
  the file stem; layouts are `layouts/<preset>.csv` for every preset listed in
  `config.yaml`; anchor pairs are optional `anchors/<configuration>_<preset>.yaml`.
  A new grating is therefore a data change with no code edit, which a test
  proves with a temporary extra file.
- **Line tables are generated, not hand-edited**: the Murphy et al. atlas is not
  shipped; `scripts/make_line_table.py` reads it from `MASTPIPE_ATLAS_DIR` and
  writes the committed `lines/<configuration>.csv` (catalogue lines in band with
  the usable flag at the instrument resolution). The header comment records the
  criteria and the counts so a regeneration must reproduce the file unchanged.
- **Layout and anchor coordinates are in the cropped frame** of `config.yaml`,
  the frame every step works in; raw-frame numbers are converted when they are
  entered.

## 2026-10-30 — Trace step

- **Profile noise is the noise of the median, not the profile's own scatter**:
  the notebook's `find_traces` thresholds the collapsed spatial profile at a
  multiple of its own robust sigma, which is a noise estimate only while most
  rows are empty. On the 45-trace survey layout (8.5 px pitch) that sigma
  measures the comb itself (on the simulated flat 7 247 ADU, so the 5σ
  threshold of 36 000 ADU exceeds every peak and nothing is detected), and on
  the three-fiber layout the median of integer-valued pixels is quantised and
  the sigma collapses to 0.1 ADU (noise peaks pass). `mastpipe.steps.trace`
  uses the statistical noise of a median across the collapsed columns
  (`sqrt(pi/2)` times the rows' robust scatter over `sqrt(n)`), which agrees
  with the notebook where the notebook's estimate is a noise estimate. The
  survey layout also needs `trace_halfwidth` below the pitch (4, against the
  HighSpec default of 8, which sets a 16 px minimum separation); the instrument
  configuration is unchanged and the tests pass the value explicitly, so the
  survey preset will need its own trace parameters when it is commissioned.
- **A lit trace is judged by the mean, not the median, of its aperture
  spectrum, against the frame's own brightest trace**: an arc spectrum is
  emission lines on a dark continuum, so its median along x is the level
  between lines (zero in the simulator). The frame first loses a per-column
  running median along y (about four aperture halfwidths), which removes a
  residual background and the smooth wings of lit neighbours; a trace is lit
  when its mean exceeds both `lit_sigma` times the noise of the mean and
  `lit_fraction` (0.1) of the brightest trace's mean, after a three-column
  running median that removes single-pixel cosmic rays (lines are wider than
  one column). An absolute noise floor alone (about 0.4 ADU per column per
  sigma on a 1598-column frame) would call every trace lit after a residual of
  0.3 ADU per pixel and turn calibration-fiber arcs into full arcs.
- **The layout shift search stays below half the smallest layout pitch**
  (`min(20 px, 0.45 x pitch)`, kernel width equal to the matching tolerance): a
  periodic layout matched one pitch off would still pair almost every row with
  the roles one fiber off, so the alias must be outside the search range.
- **The collapse statistic is a parameter** (`TraceParams.collapse`: `median`,
  the notebook's choice and the default; `mean`; `percentile` with
  `collapse_percentile = 90`): a median along x cannot see a trace that carries
  only emission lines (on the simulator's full arc it finds no trace), so when
  the runner falls back to the full arc for trace finding it passes
  `collapse="percentile"`, which recovers all 45 simulated survey traces within
  0.2 px at the detector centre. The detection threshold and the centroid
  threshold use the chosen statistic's own standard error.

## 2026-10-30 — `identify` command

- **The anchor panels rank template lines by their height in the convolved
  template**, not by an `snr` column of the usable-line table: the panels take
  either the atlas (`flag_usable_lines` output) or the shipped line table (no
  `snr`) as the template, and the height at the catalogue wavelength is what
  the eye matches on the figure in both cases. The measured panel keeps the
  notebook's `height` ranking.
- **`run_identify` / `run_solve` take the layout preset name**, not the layout
  table: the name selects the stored anchors and is the reference's
  `fiber_preset`; the shared `_arc_lines` chain, which the runner reuses, takes
  the table. The anchors file carries the trace its pairs were measured on
  (`trace:`), read by `load_anchor_trace`, so a command without `--trace` uses it.
- **Masters are picked up from `out/calib` by kind** (`*_bias.fits`,
  `*_background.fits`, `*_flat.fits`, the latest of each by name); a dark
  master is not used by `identify`, and without masters the frame is only
  given its variance and mask. The measured panel is mirrored when the
  configuration's `dispersion_sign` is negative.

## 2026-10-30 — Calibration library and quick-look

- **The idempotence key of a science product is the input checksum, the
  configuration hash, the reference identifier and the sorted checksums of
  the drift arcs it used** (revision 1 §2.2). A product on disk with the same
  key is skipped; a product whose key changed is reduced again and listed in
  `RunSummary.reduced_again` — this covers the planned case (a `reference`
  product whose exposure now has a drift arc) and equally a `nominal` product
  once a reference exists or any product after a reference change, since the
  key, not the status, is what decides. `--force` reduces regardless.
- **Full arcs are idempotent through the library**: a full arc whose checksum
  is among a stored reference's `source_checksums` is not re-solved (its
  reference is already there), unless forced. Re-solving the same arc gives
  the same identifier (arc time), so a re-solve replaces the library entry
  rather than adding a second one.
- **The index is the library's source of truth** and every reference entry
  carries its camera, because reading a reference back needs the camera's
  pixel size for the grating geometry; `identify` and `solve` register the
  references they write, so a following `quicklook` finds them without a
  manual step. Reference validity is compared on parsed timestamps.
- **Drift arcs are reduced on the arc fibers only** (a trace set restricted
  to the `arc` role goes into the shared `_arc_lines` chain), which is both
  the physics (nothing else is lit) and a five-fold saving in line finding.
- **A night without a bias master** (in the night or the library) is a
  group-level failure rather than a silent reduction without bias
  subtraction; the other masters are optional.
- **A drift model is tied to the reference it was fitted against**: the runner
  keeps the reference identifier with every drift model and pairs a frame only
  with the arcs fitted against its own reference. A full arc solved mid-night
  starts a new reference; under `nearest` or `interpolate` an arc fitted against
  it could otherwise be applied to a frame still using the earlier one. The
  blended model of `interpolate` carries per fiber the summed normalised weight
  of the models holding that fiber in its `WEIGHT` column; the per-arc weights
  are in the science QC JSON (`drift_weights`, aligned with `drift_arcs`).

## 2026-10-30 — Final fix wave of Plan 1

- **Fibers are keyed by the layout's trace id, not by detection order**
  (reverses the trace-index keying the trace step started with). `find_traces`
  still numbers detections by position, but `layout_indexed` re-keys the
  matched traces so that `Trace.index == layout_index` and drops the
  detections that match no layout row; the runner and `identify` both use it,
  so every per-fiber key downstream — `Spectrum.fiber`, the reference's
  `solutions`, the drift model's `per_fiber`, the throughput table — names the
  same physical fiber on every night. A fiber missing from one night's flat is
  absent from its products instead of shifting every fiber above it onto its
  neighbour's solution; an ambient-light detection (the lab flat shows two) is
  dropped, logged and listed in the trace QC with its position. As a guard,
  `apply_drift` withholds the wavelength of a fiber whose trace centre sits
  farther than `trace.layout_tolerance` from the solution's `trace_y_center`
  and the runner records it under `wavelength_skipped`; a reference solved on
  another pixel grid is refused with a `ValueError` naming it. `solve_full_arc`
  and `identify --trace` name fibers by that id as well.
- **Header-map options `missing_values` and `reference`**: a keyword whose raw
  or converted value is a listed sentinel (`TEMP = -999` on the Andor camera)
  counts as absent, and a time field declares the instant it records
  (`start`, `mid` or `end`; HighSpec's `DATE-OBS` / `DATE` is the start).
  Ingest derives `mid_time` from the start (or end) and the exposure time;
  without an exposure time the time is kept as given and `mid_time` is listed
  as unresolved. The simulator writes `DATE-OBS` as the start, so the truth's
  `mid_time` is reached the way the pipeline reaches it.
- **Per-preset parameters** live in `config.yaml` under `layout_params`
  (`preset -> partial StepParams`), layered by `build_effective(layout=...)`
  between the grating block and the command-line overrides; `survey45` narrows
  the trace, line and extraction halfwidths to 4 px for its 8.5 px pitch, and
  its layout table fixes the extraction aperture at 4 px. The tests no longer
  pass those values as overrides: the default path is what runs.
- **A night folder is validated**: a path that is not a directory raises
  `FileNotFoundError` before anything is written (the CLI exits 2), and a
  folder without a readable frame records a night-level failure instead of
  logging "scanned 0 frames" and exiting 0.
- **The layout preset joins the idempotence key** of a science product and is
  stored as `MPLAYOUT` in every product header: two layouts can share an
  effective configuration (and hence a hash) and still name different fibers.

## 2026-10-30 — Rulings made during Plan 1 and recorded at its close

- **`config_hash` covers the values of the effective configuration only**:
  `EffectiveConfig.layers`, which records how the values were reached
  (instrument, configuration block, overrides), is excluded from the hashed
  JSON. The hash names *what* was used, so the same effective configuration
  hashes the same on every platform, and an override that merely restates a
  value (as the command line does when it repeats a default) does not make a
  product look different from one reduced without it.
- **The first reference is reproduced from the notebook's stored anchor pairs, not
  imported from the notebook's solution** (extends call 2 of the revision):
  the two `(wavelength, pixel)` pairs read off the notebook's anchor figure
  are stored with the instrument (`anchors/mgb518_lab3.yaml`, with the trace
  they were measured on) and `identify --stored-anchors` re-solves every fiber
  from them. The reference therefore carries the pipeline's own provenance
  and line list, and the notebook's numbers remain an independent check: the
  laboratory regression compares the two rather than asserting an import.
- **One comb-matching primitive**: `comb_shift` (the Gaussian-kernel
  cross-match of two position combs with a parabolic peak) serves both the
  global shift of `match_layout` (trace rows against the layout) and the
  shift guess of `fit_fiber_drift` (measured line columns against the
  reference's identified lines). The kernel width is a parameter set by each
  caller to the scatter it expects between matching members — the layout
  tolerance in one case, `drift.tol_px` in the other — rather than the grid
  step.
- **The cosmic-ray footprint runs along the dispersion axis** (`(1, size)`
  for an integer `cosmic_size`, a `(rows, columns)` tuple when given): a
  spectral line is extended across rows by the fiber image but is as sharp as
  a cosmic ray along a column, so a square median footprint flagged the cores
  of bright emission lines (the simulated sky lines) while the along-row one
  does not. Tracks of several pixels along a row are a known gap until
  profile-based rejection.
- **The simulator's drift law is anchored at the full arc**: the shift and
  stretch run from their start value at the full arc's mid-time to their end
  value at the last frame's, frames before the arc carrying the start value,
  so with the default `drift_px = (0.0, 0.6)` the full arc itself is
  undrifted. A reference solved from it then has zero drift by construction
  and every later measured drift is relative to it, as it is at the telescope.
- **The usable-line table is generated with the notebook's convolution
  kernel**: a normalised Gaussian truncated at three sigma applied with
  `numpy.convolve`, not `scipy.ndimage.gaussian_filter1d`. The two differ in
  the wings enough to move the 514.6056 nm Th I line across the SNR 20 cut
  (20.4 against 19.996), so the committed `mgb518` table has 31 usable lines
  of 73 and `scripts/make_line_table.py` must reproduce it unchanged.
- **The laboratory regression asserts the layout-matched fibers, with the
  flat calibrated as the notebook does** (bias subtracted, lamp-off background
  kept: `flat_subtract_background` is off). On the real flat the median-noise
  trace estimator also detects two faint ambient-light features (S/N 6–10)
  that the notebook's profile-scatter sigma swallows; `match_layout` sets them
  aside as `extra`, they never solve and they vanish when the background is
  subtracted, so the regression checks the three `arc`-role traces against
  the notebook's centres and leaves the raw count to the trace unit tests.
  The frames are reached only through `MASTPIPE_TEST_DATA` with the layout the
  README documents; nothing about their location is in the repository.

## 2026-11-02 — Public repository

- **The GitHub repository is public** under `The-MAST-project/mastpipe`, like the
  other MAST repositories, replacing the earlier private-for-now intent. The
  privacy rules in `CLAUDE.md` (no paths, hosts, addresses, credentials or real
  data in any tracked file) and the CI privacy scan are what make this safe.
- **No individuals are named in tracked files.** Prose, comments and test
  identifiers credit the HighSpec laboratory reduction notebook as a thing, not
  its author; authorship is "The MAST project" (license, package metadata, docs).
- **`CLAUDE.md` and `DECISIONS.md` are committed**, as in the other MAST
  repositories (nine of which track a `CLAUDE.md`, four a `DECISIONS.md`);
  the design documents under `docs/design/` stay local.

## 2026-11-04 — Plan 1.5: hardening at the boundaries

Recorded from the Plan 1.5 design (`docs/design/2026-11-04-mastpipe-plan-1.5-hardening-design.md`,
§4), whose decisions D1–D4 were taken with the project lead and whose calls C1–C7 the
designer made, and from the rulings made while the plan was executed.

- **Canonical identifiers are `MAST_common`'s spellings** (D1): gratings `Mg`, `Halpha`,
  `Ca` (and `Future`, `Empty` as known but unsupported), instruments `highspec` and
  `deepspec`, telescopes `mast01`…`mast20`. The fleet's control software writes these
  names, so mastpipe's own spellings (`mgb518`, `halpha`, `caiihk`) survive only as
  file-name aliases of the laboratory data.
- **`MAST_common` is mirrored, not imported** (D2): `mastpipe.vocab` restates the
  literals and is what the code uses; a CI job checks `MAST_common` out as a sibling
  (the paired branch, else `master`, as `MAST_control` does) and a test fails when the
  mirror diverges. mastpipe stays installable without the fleet's library.
- **The design documents are published** (D3): `docs/design/` leaves `.gitignore` and is
  committed once the privacy scan passes and every personal name is replaced by a role
  (the HighSpec laboratory reduction notebook and its author, the DeepSpec instrument
  team, the project lead). This reverses the 2026-11-02 entry's "stay local", so that
  the citations in this file resolve. The documents are not pages of the Sphinx site.
- **An arc that the header calls `light` is detected from its content** (D4): `light`
  becomes a provisional exposure type that the runner refines to `arc` or `science`
  from the frame, before the trace stage (on a night without a flat the full arc is the
  trace source); a manifest still overrides. The control software writes `light` for
  every shutter-open exposure and no lamp card, so no header rule can tell them apart.
- **Grating ids are case-sensitive** (C1) and spelled exactly as `MAST_common` does; the
  header map no longer lower-cases `configuration`.
- **The time-zone question stays open, and the assumption is explicit** (C2): every
  header-map time field declares a `timezone` (`utc`, an IANA zone, or `assume-utc`); a
  stamp with an offset is honoured; HighSpec ships with `assume-utc`, which logs one
  warning per run and records `MPTZ = 'assumed-utc'` on every product, until one frame
  with a known UTC time settles the control computer's clock.
- **Python 3.11 support stays** (C3): the fleet pins 3.12 for services, but mastpipe is a
  library. Ruff's rule set, line length (125) and pin (`ruff==0.16.0`) follow the fleet,
  with `NPY` added and physics-named variables allowed by per-file ignores.
- **`FIBPRE` leaves the header map** (C4): in MAST a "fiber preset" is a stage position
  (`deepspec` / `highspec`), not a fiber layout. The layout comes from `--layout`, a
  manifest's `fiber_preset` or the instrument's default.
- **Masters enter the idempotence key by checksum** (C5), and an equal master already in
  the library is reused rather than rebuilt and overwritten on every run; `--force`
  rebuilds.
- **Numerical changes are limited to four bounded items** (C6): master variances in the
  error budget, saturated lines out of the wavelength fit, a drift gate that needs three
  usable arc fibers (`drift.min_fibers`), and versions never dropped in combination
  (`MPVMISS`). The cosmic-ray filter question goes to Plan 2.
- **Calls after the builder check** (C7): arc detection is frame-level and runs before
  the trace stage; product names gain a short input-checksum segment, so two frames with
  the same mid-time to the second never share a name; ingest gates on `INSTRUME`; a
  dropped calibration frame is a recorded failure with exit status 2; the fleet test
  compares against `common.spec` only, reading the configuration lists as text; the
  laboratory regression keeps manifest rows for the flat and the background.
- **Any conditional header-map rule may refine a `light` frame**, not only the bias rule
  the design names: a `light` frame with a lamp card is typed by it, and one with no
  exposure time is `unknown` at ingest rather than resolved from its content. The rules
  already express what a header can prove; content is the fallback for what it cannot.
- **`mid_time` strings stay naive ISO UTC** (`2026-10-01T18:02:30`), not the design's
  `Z`-suffixed form: parsing always yields an aware UTC instant, but the stored and
  printed form is unchanged from Plan 1, so product names, identifiers and tests keep
  their shape. The design document was corrected accordingly.
- **An unopenable FITS file is a failure; a non-FITS file is not scanned**: a file with a
  FITS suffix that astropy cannot open is recorded in `failed` (exit status 2), because
  a supervisor should see it; anything else in the folder (`seq.txt`, notes) is not a
  frame and is ignored. The design's "skipped" means "not scanned".
- **The expected frame shape comes from the detector configuration**, not from the first
  good frame of a batch, which was order-dependent: the crop and shape are unbinned
  pixels, and `expected_shape(detector, binning)` and `load_frame` both scale the crop by
  the binning factors. A group's binning is the most common among its science (and
  unresolved `light`) frames, computed before the masters are built.
- **Master reuse matches the inputs and a calibration-parameter hash, not the
  configuration**: `find_by_inputs` compares kind, camera, input checksums and
  `calib_hash` (the `calib` parameters and the detector), so a changed stacking method
  rebuilds the master and moves the science key. Equal input checksums already imply the
  binning and readout. Since frames that name no grating take the run's
  `--configuration`, the same frames may be reduced under two configurations; stacked
  from the same frames under the same parameters they give the same master (a flat also
  needs the same masters applied, which reuse checks), so ignoring the configuration is
  benign — the reused master merely keeps the name of the configuration it was first
  written under.
- **Frames that name no grating take the run's `--configuration`** before grouping: their
  products, masters and library entries are named and indexed by that grating (`Mg`),
  not `noconf`, so a later night's library lookup and `identify` / `solve` find a site
  night's masters (the control software writes no `GRATING` card today). This resolves
  the Plan 1 deferred item on `noconf` masters.
- **`identify` and `solve` choose masters as the runner does**: through
  `CalibLibrary.find` on the arc's camera, binning, readout key, detector temperature and
  mid-time (the flat on its configuration too), honouring `--library`, instead of the
  latest file name in `OUT/calib`, which could hand an `image/conventional` arc an
  `image/em` bias. A consequence: a master whose temperature is more than 2 K from the
  arc's is not applied. The laboratory full arc (−50 °C) thus no longer takes the
  laboratory bias (−40 °C); its background (no temperature) and flat still apply, and
  the regression's numbers are unchanged (0.016 / 0.033 / 0.012 px, 23 / 25 / 20 lines).
  `identify --library` also adds its reference to that library, as `solve` does.
- **A layout's hash ignores line endings**: `layout_hash` reads CRLF as LF and
  `.gitattributes` keeps CSV files LF, so a Windows checkout does not change the
  idempotence key and reduce every product again.
- **Plan 1 output folders are not migrated**: product names changed, so old products are
  never found or replaced, only left beside the new ones; the upgrade notes ask for a
  fresh `--out` or cleared `spec/` and `qa/` folders, and `quicklook` warns once per run
  when `spec/` holds old-style names. The calibration library may be kept: its old
  entries never match.
- **The missing-configuration message lists the instrument's configurations**, not every
  fleet grating name: `Future` or `Ca` would be refused by `build_effective` anyway.
- **The master variances propagate exactly**: a frame of exposure `t` calibrated with a
  background of exposure `t_bkg` is `raw − (1 − s)·bias − s·background` with
  `s = t / t_bkg` and carries `(1 − s)²·V_bias + s²·V_background`; with a dark, `s = t /
  t_dark` and it carries `(1 − s)²·V_bias + s²·V_dark` (the dark-rate variance is
  `V_dark / t_dark²`, with `t_dark` carried by the master). The bias enters both subtractions, so its variance is
  correlated and must not be counted twice.
- **Frames of another readout than their group are dropped per frame**: the readout key
  is `READMODE`/`OUTPTAMP` (`image/em`), the group follows the most common key among its
  science and not-yet-resolved `light` frames (as for the binning). A calibration frame
  with another key (a bias, dark, background or flat, or a full arc competing as the
  trace source) is left out and recorded in `skipped` with both keys, not as a failure:
  it belongs to another set (an EM bias series beside conventional science), and
  `build_master` refuses mixed amplifiers and would otherwise fail the whole group. A
  science, `light` or arc frame with another key, or another binning, is a failure,
  since it would be reduced with masters and noise of another readout. A night mixing
  frames with and without those cards loses the uncarded ones when the group's key is
  known.
- **The simulator's default jitter and gradient are deliberately small**
  (`dispersion_jitter_px = 0.3`, `drift_gradient_px = 0.05`), below what the end-to-end
  tolerances could detect, so the fixture nights reduce with stable tolerances.
  Sensitivity to fiber-identity and drift-row mistakes is proven by dedicated tests
  instead: the end-to-end swap test at a jitter of 1 px and the drift-fit test at a
  2 px gradient, which a drift flat in `y` fails.

## 2026-11-10 — Tidy-up after Plan 1.5

- **The quick-look runner is split into stage modules**, moved without any change of
  logic: `pipeline.group` (run options and summary, the per-group state `_Group`, the
  folder and policy constants, the logger and the helpers shared by the stages and the
  runner), `pipeline.masters`, `pipeline.traces`, `pipeline.arc_stage` and
  `pipeline.science`; `pipeline.runner` keeps the recipe's description, the group loop
  and `run_quicklook`, and re-exports the public names. **A stage module imports
  `pipeline.group`, never another stage or `pipeline.runner`**: a helper two stages need
  moves to `group`.
- **One length for every short hash**: the configuration, calibration, layout, masters
  and reference hashes keep the first `core.hashing.HASH_LENGTH` (12) hex characters of
  a SHA-256. `core.hashing` imports nothing, so `config`, `instruments`, `steps` and
  `pipeline` all use it without a cycle. (Frame checksums keep their own
  `io.ingest.CHECKSUM_LENGTH`: they identify files, not configurations.)
- **A library flat is reused by the content of its masters, compared by checksum**: the
  flat records `MPMKH`, the hash of the file checksums of the masters it was calibrated
  with (`_masters_hash`, the same hash a science key carries), and is reused only when
  the group's masters hash the same, no longer when their file names merely match. A
  flat written before the card existed reads as unset and is **built again once** rather
  than trusted: its provenance holds only the masters' names, which cannot prove the
  content. The rule that a master built in this run forbids reusing the flat stays (it is
  now redundant but cheap).
- **The flat batch is chosen among the frames that survive**, like the other kinds: every
  exposure-time batch is loaded and calibrated, and the one with the most survivors (the
  longer exposure on a tie) is stacked; the losing flat batches are not stacked.
- **A flat-derived trace set is named with the flat's first input checksum** (the one in
  the flat's own name), through `_master_meta`, rather than `nosum`, so trace sets of one
  stamp built on flats stacked from different frames get different names.
- **`DriftModel.n_usable` stays `int | None`**: an omitted count is derived from
  `per_fiber["usable"]` in `__post_init__`, which test fixtures and the reading of products
  written before `MPDNUSE` rely on; a default of 0 would silently report zero usable
  fibers instead, so the two `or 0` reads stay.
