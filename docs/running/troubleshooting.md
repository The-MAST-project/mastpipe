# Troubleshooting

The quick-look never stops at a bad frame unless told to (`--on-error raise`).
Instead it records what it skipped and why, in the console summary, the log
under `out/logs/`, and the QC files. This page lists the messages you will
meet and what to do.

## Frames

"skipped <file>: exposure type unknown (rule:N)" or "(unclassified)"
: The header map could not tell what the frame is: no exposure type keyword
  and no rule matched (a frame with no exposure time is deliberately left
  `unknown`, never scaled by a guess). Add a manifest row for it.

"no readable frame in the night folder" / "night folder … is not a directory"
: The path is wrong, or the folder has no FITS astropy can open. The command
  exits with status 2.

"FAILED <file>: …"
: An exception in that frame's reduction; the traceback is in the log. The
  other frames are unaffected. A failure in the masters or trace finding fails
  the whole camera-and-grating group instead, since nothing downstream can run.

## Arcs and references

"no reference for …; run identify"
: A full arc was found but the library has no reference for its grating and
  layout to seed the solve from. Make one with `identify` (see
  [The first reference solution](first-reference.md)); the next quick-look
  will solve the arc. A calibration-fiber arc with no reference valid at its
  time is skipped the same way ("no reference valid at …; run identify"), and
  an arc whose lit traces match neither pattern is skipped as
  "arc kind unknown (lights N of M traces)".

Status `reference` instead of `transferred`
: No calibration-fiber arc paired with the frame: none in the folder, none of
  the same grating, none fitted against the same reference, none within
  `drift.max_age_hours` (default 12 h), or the arc's drift fit was unusable (too
  few lines or too large a scatter; see its QC). The spectrum is still
  wavelength-calibrated by the reference alone. Early-night frames reduced
  before the first arc get this status and are re-reduced when the arc lands.

Status `nominal`
: No reference solution exists for this grating and layout; the wavelengths
  are from the configured geometry only and are approximate.

A fiber under `wavelength_skipped` in the science QC
: Its trace centre lies farther from the reference solution's than
  `trace.layout_tolerance`, so the solution was not transferred to it. The
  fiber keeps its raw extraction with no wavelength. Usually the layout or the
  reference is stale for this night: make a new reference from a full arc.

`extra` traces in the trace QC
: Detections that match no layout row (in the laboratory flat, two faint
  ambient-light features at S/N 6–10). They are dropped from the products and listed with
  their positions. If a real fiber is among them, the layout table is wrong.

`missing` traces in the trace QC
: Layout rows with no detection within `trace.layout_tolerance`: a dark fiber,
  or a layout whose nominal positions are off by more than the tolerance.
  Missing fibers are absent from the products.

## Sky subtraction

Reason `no sky fiber`
: The layout has no `sky` row for that science fiber's telescope, or the sky
  fiber was not found. The science fiber keeps its raw extraction only.

Reason `sky fiber fully masked`, `no wavelength`, `bad throughput`
: The sky fiber's spectrum is unusable, one side has no wavelength solution,
  or the flat gave a non-positive throughput for either fiber. Same outcome:
  raw extraction only, reason recorded.

## Products

"already reduced" for every frame
: The products exist with the same key. Use `--force` to redo them, or change
  what should change (policy, overrides), which changes the key.

Provenance disagrees with the exposure
: `write_spectra` refuses to write a product whose provenance statuses differ
  from the exposure's. This is a pipeline bug, not a data problem; report it
  with the log.

## Getting more detail

```bash
mastpipe show <product.fits>            # statuses, per-fiber summary, provenance
cat out/qa/<product>_spec.json          # the numbers behind the summary (named like the product)
tail -n 50 out/logs/quicklook_*.log     # the run's log
mastpipe quicklook ... --on-error raise # stop at the first error with a traceback
```
