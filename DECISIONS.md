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
