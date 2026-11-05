# What this release does not do

This release reduces HighSpec end to end on the simulator and on the
laboratory frames. Known limits, and where each is headed:

| topic | state | planned |
| --- | --- | --- |
| DeepSpec | configuration slots exist (model kind, orientation, per-camera grouping) but no DeepSpec definition ships | next release: the four bands, the dispersion-equation model, the krypton-then-Th-Ar bootstrap, the on-sky frames of September 2026 |
| Flux calibration | pass-through; status `uncalibrated` | once the standard-star observing policy is decided |
| Telluric correction | pass-through; status `uncorrected` | with the flux calibration |
| Barycentric corrections | not computed; the header vocabulary has no pointing | next release, stored per telescope, never applied to the grid |
| Optimal extraction | aperture extraction only | Horne extraction with the flat-derived profile, incl. cosmic-ray rejection from profile residuals |
| Combination across exposures | within an exposure only | time grouping: every N frames, a window, a night |
| Flat-fielding within apertures | not implemented | when a properly illuminated flat exists |
| Bad-pixel mask | saturation and non-finite pixels only | from detector parameters and master-frame outliers |
| Cosmic rays along a row | survive the along-dispersion filter | closed by the optimal extraction |
| Sky resampling | linear; a few per cent of a line peak at sub-pixel offsets | a flux-conserving resampler |
| Manual mode | the quick-look runner with `--on-error raise`, `--time-policy interpolate`, `--force` | a `reduce` command with per-stage `--force` |
| Quality pages, watch mode, focus recipe, PypeIt adapter | not present | a later release |
| Arc lamp for the DeepSpec calibration fibers | unknown (HighSpec: Th-Ar) | to be settled with the instrument team |

## Open questions carried from the design

- How standards will be observed decides the flux-calibration recipe; the
  code provides the per-telescope response interface and nothing else.
- Whether a smooth drift model across five calibration fibers reaches the
  few-km/s goal is an empirical question for commissioning; the fiber-to-fiber
  agreement figure is the acceptance test.
- The HighSpec detector's gain and read noise in the instrument file, the
  defaults and the per-amplifier (`conventional`, `em`) values alike, are
  placeholders to be measured.

## Why these limits are safe

Every open stage is an interface with a pass-through default that records its
status in the product, so nothing downstream mistakes an uncalibrated spectrum
for a calibrated one; the deferred stages have slots in the products and the
status vocabulary, so adding them is a change of code, not of products.
