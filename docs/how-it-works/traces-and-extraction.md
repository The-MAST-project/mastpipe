# Cleaning, traces and extraction

## Cleaning

The calibration arithmetic is the notebook's. With a lamp-off background
exposure of the science exposure time, which the laboratory setup uses:

```text
clean = raw − bias − (background − bias) · t / t_background
```

otherwise with a master dark rate `(dark − bias) / t_dark`:

```text
clean = raw − bias − dark_rate · t
```

Masters are median stacks (per exposure time for dark and background, with
optional sigma clipping) carrying validity keys — camera, binning, readout
mode, temperature, time — that the calibration library matches on. The
variance per pixel is the Poisson term of the detected signal plus the read
noise, in ADU². Saturated pixels (in the frame or in a master) and non-finite
values are masked.

**Cosmic rays** are flagged on science frames as residuals above a median
filter. The filter runs *along the dispersion axis* (a 1 × 5 footprint): a
square footprint flags the cores of sky emission lines, which are narrow along
the dispersion but extended across the fiber, as cosmic rays. The price is
that a multi-pixel track lying along a row survives; profile-based rejection
in the optimal extraction of the next release closes that gap.

## Traces

Traces are found on the master flat: a band around the detector centre is
collapsed along the dispersion axis, a running-median baseline is removed, and
peaks above a noise threshold are followed across the detector in blocks of
columns, each trace fitted by a clipped quadratic (the outer fibers curve by
the conical-diffraction signature the notebook measured; a straight line
leaves it in the residuals). Two changes from the notebook:

- the **noise** of the collapsed profile is the statistical noise of the
  median across the band's columns, not the profile's own robust scatter — on
  a dense 45-trace comb the latter *is* the comb and nothing is detected;
- when a night has no flat, the traces are found on a **full arc** with a
  percentile collapse (a median collapse of a line-only frame is dark).

Detected traces are then **matched to the fiber layout**, an input table that
says which trace is which: role (arc, science, sky), telescope, nominal
position. The number of traces is never inferred from the data. Matching
estimates one global shift and assigns layout rows to detections within a
tolerance; the matched traces are re-keyed by the layout's trace ids, unmatched
detections are dropped and listed (the laboratory flat shows two faint
ambient-light features at S/N 6–10), and layout rows without a detection are listed
as missing. On the survey layout the search is bounded to less than half the
8.5 px pitch so a one-pitch alias is impossible.

An arc frame is classified as **full** (every matched trace lit) or
**calibration-fiber** (the arc traces lit, the science and sky traces dark)
from a per-trace level above the noise and above a tenth of the brightest
trace, after a spatial running-median baseline. The flat also gives the
**relative throughput** per fiber used in sky scaling.

## Extraction

Aperture extraction follows each trace with a fixed half-width (9 px for the
laboratory fibers, 4 px for the 8.5 px survey pitch; the layout's `aperture`
column overrides per fiber) and the notebook's three aperture rules: whole
rows, floored rows, or sub-pixel edge weights. Flux is the weighted sum,
error the weighted quadrature sum of the variance, and a column is masked
when its coverage falls below a minimum or any pixel in it is masked. For arc
frames an asymmetrically clipped Chebyshev continuum is fitted and stored; the
continuum-subtracted flux feeds the anchor panels and the cross-correlation
that seeds the other fibers' solutions (the line finder works on the
straightened two-dimensional strip with its own row-wise baseline); the raw
flux is never replaced.

**Line finding** on arc frames is the notebook's two-dimensional method: the
strip around each trace is straightened, its rows continuum-subtracted, peaks
are detected, nearby lines are fitted jointly and deblended by adding Gaussian
components while the residual still peaks above the noise, and every line
carries flags that say why it would or would not be used — isolated, near an
edge, complete, on the trace, fit converged, saturated. The wavelength solve
uses the lines marked `use_x`; the drift fit uses the same column. Two
departures from the notebook: `complete` is judged on the measuring rows
only (so a trace near the detector edge keeps its lines), and a deblended
component pinned at the narrowest allowed width is rejected as a noise spike.

## What the next release changes here

Horne optimal extraction with the flat-derived spatial profile, which also
rejects cosmic rays from profile residuals; and pixel-to-pixel flat-fielding
within the apertures, to be added once a properly illuminated flat exists.
