# Decisions

A dated log of design and workflow decisions for mastpipe, newest last. Each
entry records what was decided and, briefly, why; reversals are recorded as new
entries rather than by editing old ones.

## 2026-09-24 — Independent package

- **Build mastpipe as an independent package** seeded from the team's
  HighSpec reduction notebook, with reduction steps as pure functions on arrays
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
