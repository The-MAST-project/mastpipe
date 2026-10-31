# Overview: six stages and a combination

A science exposure passes through six stages in a fixed order; combination
across telescopes follows per target. Each stage is a pure function on arrays
and small dataclasses; the runner chains them and does all the reading and
writing.

| stage | what it does | where |
| --- | --- | --- |
| 1 Clean | bias and dark (or a matched-exposure background) subtracted; variance from gain and read noise; saturated and non-finite pixels masked; cosmic rays flagged on science frames | `steps.calib`, `steps.cosmics` |
| 2 Trace and extract | traces found on the flat, matched to the fiber layout and keyed by its ids; aperture extraction per fiber with error and mask | `steps.trace`, `steps.extract` |
| 3 Wavelength | the fiber's stored reference solution, corrected by the drift measured on the calibration-fiber arcs nearest in time | `steps.wavecal`, `steps.drift`, `pipeline.arcs` |
| 4 Background | the telescope's sky fiber, scaled by relative throughput and resampled in wavelength, subtracted | `steps.sky` |
| 5 Flux calibration | per-telescope response from a standard star — **method open**; pass-through | `steps.fluxcal` |
| 6 Telluric correction | **method open**; pass-through | `steps.telluric` |
| Combination | the target's science fibers combined on a common grid, every version separately | `steps.combine` |

The quick-look runs stages 1–4 and the combination; stages 5 and 6 run once
their calibrations exist. Every stage keeps what the previous one produced:
the product carries the raw extraction and each later version side by side
([Products policy](products-policy.md)).

## Two run modes, one code path

The observatory's **quick-look** is non-interactive: it uses the most recent
valid calibrations, degrades gracefully (a frame without a paired arc still
gets its reference wavelengths), skips a bad frame and goes on, and is
idempotent so it can run after every exposure. The **manual** reduction is the
same runner with stricter choices: stop at the first error, interpolate the
drift between the arcs before and after the exposure once the night is
complete, redo a stage with `--force`. The products are identical in form.

## What is inherited

The scientific core descends from the HighSpec laboratory
reduction notebook: the calibration arithmetic, trace finding, aperture
extraction, continuum fitting, line finding with deblending, and the
grating-equation wavelength solution anchored on two identified lines and
refined against the Th-Ar atlas of Murphy et al. The ports are one-to-one
where the notebook's assumptions hold; the places where they were changed, and
why, are named in the following pages and in the project's decision log.

## What was added

Everything the notebook did not need in the laboratory: a declarative
instrument layer, the fiber layout as an input table, the three-cadence
wavelength scheme with stored references and drift transfer, sky subtraction,
combination, FITS products with provenance, a calibration library, a
simulator with known truth, and the unattended runner.
