# Sky subtraction and combination

## Sky subtraction

Each telescope feeds a science fiber and a sky fiber, so every science
spectrum has its own sky measured at the same time through the same optics.
The pairing comes from the fiber layout, never from the data.

The two fibers have different wavelength solutions, so the subtraction is done
**in wavelength space**: the sky fiber's spectrum is scaled by the ratio of the
two fibers' relative throughputs (from the flat), resampled linearly onto the
science fiber's wavelength grid with the error propagated exactly through the
interpolation weights, and subtracted. Pixels outside the sky fiber's
coverage, or next to a masked sky pixel, are masked in the result.

The subtracted sky is kept in the product (`SKY`), and the result is a new
**version** of the spectrum (`FLUX_SKYSUB`, `ERROR_SKYSUB`, `MASK_SKYSUB`);
the raw extraction stays untouched. A science fiber whose sky fiber is
missing, fully masked, without a wavelength solution, or whose throughput is
unusable keeps its raw extraction only, and the reason is recorded.

"Background" here means the sky fiber's light. The host galaxy under a
transient is not removed by this step; that is one reason the raw extraction
is kept.

Linear resampling of a 3-pixel-wide sky line between grids offset by a
fraction of a pixel leaves a few per cent of its peak as residual; a
flux-conserving resampler is a later improvement.

## Combination

After the stages, the science fibers of one target in one exposure (several
telescopes on the same object) are combined:

- a **common wavelength grid** spanning the contributors, log-uniform by
  default, with the step of the finest contributor;
- each contributor's arrays resampled onto it with the same exact error
  propagation as the sky step;
- an **inverse-variance weighted mean** per pixel, error `1/√Σw`, masked where
  no contributor is unmasked, with a one-pass outlier clip where three or more
  contributors overlap;
- **every version** the contributors share combined separately — raw,
  sky-subtracted, and the later ones once they exist — so the combined file
  mirrors the per-exposure one.

Contributors can be chosen: a telescope whitelist or blacklist from the
configuration, and an optional automatic rejection (off by default) on median
signal-to-noise, flux ratio to the group, or the drift measured for the
exposure. Whatever is excluded is **listed with its reason**, not dropped, and
stays in the per-exposure file. (The drift criterion needs the exposure's QC
numbers, which the quick-look does not pass yet, so only the other reasons
can fire today.)

Combination across exposures in time (every N frames, a time window, a
night) is the next release's work; the grouping rule is already recorded in
the product so the two will read the same way.

## Flux calibration and telluric correction

Stages 5 and 6 exist as interfaces with pass-through implementations. How
standard stars will be observed (all telescopes on one standard at once, or
otherwise) is not yet decided, and the telluric method (standard-star division
or a model atmosphere) follows from it. Until then the statuses in every
product read `uncalibrated` and `uncorrected`; the version columns for these
stages exist in the spectra table but are unpopulated (NaN, masked, `HAS_`
false), and the combined table omits them.
