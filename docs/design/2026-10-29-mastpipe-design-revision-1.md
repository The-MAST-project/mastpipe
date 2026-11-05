# mastpipe — design revision 1 (after the team meeting)

*Supersedes the named sections of `2026-09-24-mastpipe-design.md`.*
Written 2026-10-29 (Claude Fable 5.1, from a conversation with the project lead
recording the outcome of the team meeting with the notebook's author and the MAST
group, held after the 2026-10-15 slide deck). Like the September document it is
written to be copied verbatim into the repository: no local paths, hostnames or
credentials.

Everything in the September specification stays in force unless a section
below replaces it. The replaced sections are named at the top of each part.

## 1. Decisions from the meeting

| Topic | September design | Decision now |
|---|---|---|
| Reduction engine | Independent package seeded from the HighSpec laboratory reduction notebook | **Unchanged.** PypeIt stays a future option behind the optional adapter once its `Fiber` pipeline is released; nothing in the first plan depends on it. |
| Instrument priority | HighSpec first, DeepSpec in the third plan | **HighSpec first, and urgently**: reduced HighSpec spectra are needed now, so the **whole HighSpec path to sky-subtracted, wavelength-calibrated spectra moves into Plan 1** (§7). DeepSpec moves to Plan 2, earlier than before, because real DeepSpec data now exist (§4). |
| Wavelength cadence | Full arc daily and on every instrument switch; science exposures carry arc light in five arc fibers | **Replaced** by the three-cadence scheme in §2. Science exposures carry **no** arc light. |
| Wavelength model | Grating equation for both instruments | **Grating formula for HighSpec, dispersion (Hartmann) equation for DeepSpec**, behind one interface (§2.4). |
| Daily solution | Fresh full solve each night | **Reference solution plus linear drift** (§2). |
| Pipeline stages | Ended at sky subtraction and combination; flux and telluric out of scope | **Six stages** ending in flux calibration and telluric correction (§3); their methods are open (§6). |
| Products | Per-exposure spectra with the subtracted sky kept | **Every individual spectrum is kept, with the pre-subtraction and pre-calibration versions alongside** (§5). |
| DeepSpec data | 2024 laboratory frames, not a reference | **The DeepSpec instrument team's October 2026 package** is the reference data set for DeepSpec (§4). |

## 2. Wavelength calibration — replaces §5 items 5–7 of the September design

### 2.1 Three cadences

1. **Full arc, every few months.** Every telescope is fed by a calibration
   lamp, so all 45 traces (and, in the laboratory, however many exist) see arc
   light. From this the pipeline fits a **new, independent wavelength solution
   per fiber**. The previous stored solution serves only as the starting guess
   and as the carrier of line identities; it does not constrain the fit. The
   result replaces the stored **reference solution** for that instrument and
   configuration (grating or band), versioned and dated.
2. **Calibration-fiber arc, one to three times per night and at every
   instrument switch.** Only the five calibration fibers are illuminated, by
   the calibration unit in the instrument room, in a **dedicated arc
   exposure**. From it the pipeline measures the **drift** of the spectrograph
   since the full arc, per arc fiber, against the reference solution. The
   drift is assumed **linear in pixel space** along the dispersion axis
   (a shift and a stretch); the fit reports whether a stretch is significant
   and falls back to shift-only when it is not.
3. **Science exposure.** No arc light. The science and sky fibers take their
   wavelengths from the reference solution corrected by the drift transferred
   from the calibration fibers (§2.2).

An instrument switch is a *trigger* for a calibration-fiber arc, not a
different kind of exposure, and it need not be visible in the data (§8). The
ingest classifies arc exposures by which fibers carry light: a frame whose arc
traces are lit and whose science traces are dark is a calibration-fiber arc;
one where all traces are lit is a full arc. The manifest can override.

**Lamps.** For HighSpec both the full arc and the calibration-unit arc are
**Th-Ar**, so the drift fit reuses the reference solution's line identities
directly. For DeepSpec the calibration-unit lamp, and whether the five
calibration fibers reach DeepSpec at all, is **to be checked** with the
instrument team (§6); until then the drift step is designed and tested for
HighSpec, with the lamp species a per-configuration parameter.

### 2.2 Drift transfer and the time policy

- The drift measured on the five calibration fibers is interpolated across the
  detector by a low-order model in trace position (as in the September §5.6)
  and applied to every science and sky fiber's reference solution.
- A science exposure is paired with calibration-fiber arcs **by time**, within
  the same instrument configuration, among arcs taken since the reference
  solution's full arc and no older than a configured maximum age
  (`drift.max_age`):
  - **quick-look mode**: the most recent arc *before* the exposure (nearest
    before; if none exists within the age limit, the nearest after; if none
    at all, status `reference`).
  - **manual mode**: linear interpolation in time between the arcs before and
    after the exposure when both exist, else nearest.
  The policy is a configuration parameter (`drift.time_policy`:
  `nearest_before` | `nearest` | `interpolate`) with those defaults per mode.
- In quick-look, an exposure reduced with status `reference` because no arc
  of the night existed yet is **re-reduced automatically** when the night's
  first drift arc arrives. The idempotency key of a per-exposure product
  therefore includes the set of drift arcs it used, and the night page shows
  the current status per exposure.
- **Wavelength status** of a science exposure, written to the product header:
  `transferred` (reference plus drift), `reference` (no usable arc since the
  full arc; the reference solution alone), `nominal` (no reference at all;
  geometry only). The status, the reference solution's identifier and hash,
  the drift-arc exposures used and the fitted drift per arc fiber are all
  recorded in the product.
- Quality numbers: per arc fiber, shift, stretch, number of lines and rms;
  across arc fibers, the scatter about the smooth model. Thresholds in
  configuration mark a drift arc as unusable (then the exposure falls through
  the policy above) and flag a drift larger than a configured limit as a
  warning on the night page.

### 2.3 The reference solution on disk

- Stored in the calibration library (September §4.6) as a FITS product per
  instrument, configuration and fiber-stage preset: a per-fiber table of model
  kind and parameters, the identified-lines table (pixel, wavelength, species,
  residual, used flag), fit rms and line count, the anchor range in pixels and
  wavelength (values outside it are extrapolated and flagged), the trace it
  was measured on, and provenance (full-arc exposures, pipeline version,
  configuration hash, date).
- Science products reference it by identifier and hash. A reference solution
  is **never edited in place**; a new full arc creates a new one. Each
  reference is **valid from its full arc's time until the next full arc** of
  the same configuration, and an exposure uses the reference valid at its own
  mid-time, so re-reducing an old night after a newer full arc exists gives
  the same products. A command-line override selects a specific reference.
- The first reference for a new grating or band is made by the manual
  `identify` command (§2.5); after that, full arcs are solved automatically
  from the previous reference. For HighSpec's Mg b grating, the first
  reference for the three laboratory fibers is **reproduced from the notebook's two
  anchor pairs**, stored with the instrument definition, through the same
  `identify` path; the notebook's published numbers are the regression target. The
  survey layout and the other gratings go through `identify` when their first
  full arc is taken.

### 2.4 Two wavelength models, one interface

A `WavelengthModel` has `wavelength(pixel)`, `pixel(wavelength)`, a parameter
vector, a kind, and a fit routine taking (pixel, wavelength, weight) triples
and an initial parameter vector. Two implementations:

- **`grating`** (HighSpec): the notebook's grating-equation solution with the
  geometry constants from the grating block and a low-order distortion
  correction, incidence angle held fixed, exactly as in the September §5.5.
- **`dispersion`** (DeepSpec): the DeepSpec instrument team's Hartmann form,
  `wavelength_air_nm = p0 + p1 / (x - p2)`, with `x` the zero-based column
  along the dispersion axis after the frame is oriented so that dispersion
  runs along columns. Three parameters per fiber. A low-order residual
  polynomial is allowed as an optional extension but is off by default,
  matching the DeepSpec instrument team's fits.

The model kind is set by the instrument configuration block, never inferred.
Both models produce the same `WavelengthSolution` product and the same drift
mechanics; drift is fitted in pixel space and is model-independent.

### 2.5 Bootstrapping a reference — the `identify` command

Replaces the automatic by-eye-anchoring ambition of September §5.5 for the
*first* solution of a configuration, which is a rare, attended job:

- **HighSpec:** Th-Ar lines anchored against the atlas convolved to the
  instrument resolution; the user confirms or supplies a few anchor pairs in a
  small text file, the pipeline fits the grating model and refines with all
  matched lines.
- **DeepSpec (two stages, the DeepSpec instrument team's method):** the resolution is too low to
  solve a Th-Ar lamp blind, so (1) a **krypton** lamp is identified first: few
  strong lines, matched by hand or from the previous reference, pinning the
  dispersion model; (2) the **Th-Ar** lamp is then matched against the Kr
  solution within a pixel tolerance and the model is refined on Kr plus Th-Ar
  with sigma clipping. Both lamps have a lamp-off **background** exposure of
  the same exposure time that is subtracted first. The Kr line catalogue per
  band (wavelength, relative strength, blend-separation limit) and the Th-Ar
  list are configuration data in the band block.

Once a reference exists, the full-arc solve of §2.1 uses it as the initial
guess and carries line identities forward; the manual step is not repeated.

## 3. Pipeline stages — replaces §1.2 scope and §5 ordering

The reduction of one science exposure runs through six stages, in this
order; combination follows per target.

1. **Clean.** Bias, dark and background subtraction, bad-pixel mask,
   cosmic-ray rejection on the 2-D frame, variance from gain and read noise.
   Flat-fielding stays designed-in and off, as before.
2. **Trace and extract.** Trace finding against the layout table, aperture
   extraction (Horne optimal extraction in the second plan), per fiber, with
   variance and mask. Traces are found on the master flat when one exists,
   otherwise on any frame in which every trace is lit (the full arc, in the
   laboratory data); without a flat the relative throughput defaults to one
   and the product carries a warning. On sky, the monthly QTH flats provide
   both.
3. **Wavelength.** Reference solution plus transferred drift (§2). Air
   wavelengths; barycentric corrections computed per telescope and stored,
   never applied.
4. **Background subtraction.** Each science fiber minus the sky fiber of the
   same telescope, scaled by relative throughput, resampled in wavelength
   space. The term "background" here means the sky fiber's spectrum; host-
   galaxy light under a transient is *not* removed by this step (see §5).
5. **Flux calibration.** Per-telescope response curve from a standard star,
   applied to each telescope's spectra. **Method open** (§6); the stage exists
   as an interface with a pass-through implementation that records
   `flux_status = uncalibrated`.
6. **Telluric correction.** **Method open** (§6); same interface pattern,
   `telluric_status = uncorrected`.

**Combination** (September §7) runs after stage 6 per target, with the same
grouping, subset and rejection options as before.

**Quick-look mode** runs stages 1–4 and combination within the exposure;
stages 5–6 run in quick-look only when a response curve and a telluric
recipe for the night already exist in the calibration library. **Manual mode**
runs all six. (Claude's call; override by configuration.)

## 4. DeepSpec — replaces §4.3 and §11 first item, and §5 item 7

Reference data set: the October 2026 package from the DeepSpec instrument team (lamp calibration
work plus raw frames), located through `MASTPIPE_TEST_DATA` and never
committed. What it establishes, and what the configuration must encode:

- **Four cameras, four bands** (U, G, R, I), each a Greateyes 1024 × 1024
  back-illuminated deep-depletion CCD, raw frames 1056 × 1027 pixels
  (overscan to be identified in the detector block), **one single-extension
  FITS file per camera** written by the camera SDK. The September assumption
  of one multi-extension file per exposure is dropped; ingest groups the
  per-camera files of one exposure by time and exposure type instead. The
  header map must support both cases, since the control software may change
  this.
- **Orientation.** Fibers are vertical stripes in the raw frame; the frame is
  transposed so dispersion runs along columns and fibers are stacked along
  rows, with no flips (`raw.T`). Orientation (transpose, flip x, flip y) is a
  detector-block parameter applied once at ingest; all downstream code sees
  dispersion along the second axis for both instruments.
- **Calibration lamps:** krypton, Th-Ar, and a lamp-off background of the
  same exposure time per band; exposure times differ per band (0.1 s to 10 s).
- **Fibers in the package:** four traces per band at slit positions
  −1.65, 0, +1.65, +2.20 mm (U has three), labelled by position and
  explicitly *not* identified with science fibers. The layout table for this
  data set therefore carries roles `unknown` and the millimetre label as the
  trace name; the survey layout is a separate file, as for HighSpec.
- **Wavelength model** `dispersion` (§2.4), reference values per band and
  fiber available from the DeepSpec instrument team's `calibration.json` files for regression
  tests: for G, rms 0.014–0.06 nm per fiber on 7 Kr plus 12–24 Th-Ar lines.
- **Science frames:** one reference star (HD 214200) and one supernova
  (SN 2026aaiv) exposure in G, R and I from 2026-09-15/16, plus a 600 s
  blocked-shutter background in R and an unvalidated 1 ms bias in G. These are
  the end-to-end test for the DeepSpec path once stages 1–4 work; they do not
  come with matched bias, dark or flat frames, so the clean stage must
  tolerate missing masters (quick-look policy) for them.
- **Headers** are camera-SDK headers, not MAST control-software headers:
  `DATE-OBS`, `EXPTIME`, `IMAGETYP` (`light`/`bias`), `SHUTTER`,
  `CCD-TEMP`, `READOUT`, binning. Band and lamp are not in the header; they
  come from the folder and filename, so the DeepSpec header map leans on the
  filename regular expression and the manifest. One frame carries a camera IP
  address in its header: a reminder that raw headers are never copied
  wholesale into products or logs.

What this changes in the build order: the DeepSpec instrument definition
(detector block with orientation, header map, band blocks with the Kr and
Th-Ar catalogues, laboratory layout), the `dispersion` model and the
Kr-then-ThAr `identify` path, tested on the DeepSpec instrument team's lamps against its reference
solutions, open the **second** plan, followed by the DeepSpec science path
(stages 1–4 on the star and supernova frames) and the four-camera product
grouping. Nothing DeepSpec-specific is in the first plan beyond keeping the
interfaces instrument-neutral (model kind, orientation, per-camera grouping
in the data model). HighSpec is the first thing to reach a usable spectrum.

## 5. Products — amends §7 of the September design

- **Every individual spectrum is kept.** Per target the products include
  every per-telescope, per-exposure spectrum, in addition to any combined
  spectrum. Combination never replaces its inputs.
- **Earlier stages are kept alongside later ones.** The per-exposure spectra
  file carries, per fiber, the extracted spectrum **before background
  subtraction**, the background-subtracted spectrum, and (when the stages
  ran) the flux-calibrated and the telluric-corrected spectra, as separate
  columns or extensions with the same wavelength grid. The subtracted sky
  spectrum is kept as before. The combined-spectrum product carries the same
  set, each combined from the corresponding per-exposure version. Rationale
  from the meeting: a transient in a host galaxy whose background was
  mis-subtracted must be recoverable from the raw extraction without a re-run.
- Product headers add `flux_status` and `telluric_status` next to the
  wavelength status, and the reference-solution identifier and drift-arc
  references of §2.2.

## 6. Open questions carried forward

- **Flux calibration.** How standards will be observed is undecided: all
  telescopes on one standard at once (one response per telescope per night)
  or otherwise. The pipeline provides the per-telescope response interface and
  the pass-through default; the recipe is written when the observing policy
  is set.
- **Telluric correction.** Method undecided (standard-star division versus a
  model atmosphere). Interface and pass-through only for now.
- **DeepSpec drift arcs.** Whether the five calibration fibers reach
  DeepSpec, and which lamp (Kr or Th-Ar) the calibration unit feeds them, is
  to be checked with the instrument team. Until then the DeepSpec wavelength
  status after a reference solution is `reference`.
- **Calibration-fiber arc detection.** Classifying full versus
  calibration-fiber arcs from which traces are lit assumes the control
  software does not label them; if a header keyword appears, the header map
  takes precedence.
- **Drift linearity.** The linear-in-pixel drift is the meeting's working
  assumption; the drift fit reports residuals so the assumption is checked on
  real data, and the model order is a parameter.
- **DeepSpec overscan and gain** are not in the package and must come from
  the instrument team.
- Unchanged from September: trace density on HighSpec's 45-trace layout,
  arc-fiber transfer accuracy, Windows CI, PypeIt adapter validity.

## 7. Consequences for the implementation plans

The three-plan roadmap is re-cut so that **Plan 1 delivers usable HighSpec
spectra**: wavelength-calibrated with drift, sky-subtracted, combined within
the exposure, on the synthetic night and on the notebook's laboratory frames. Plan 1
is rewritten against this revision before execution.

**Plan 1 — HighSpec end to end.** Everything in the September Plan 1
(skeleton, CI, privacy scan, data model, configuration, HighSpec definition,
header-map ingest, simulator, detector calibration, traces, aperture
extraction, products, calibration library, runner, `quicklook`) with these
changes and additions:

- Simulator: no arc light in science frames; separate full-arc and
  calibration-fiber-arc frame types with a controllable linear drift between
  them and a controllable time sequence within the night.
- Data model: `WavelengthModel` with the `grating` kind (the `dispersion`
  kind's slot exists, implementation in Plan 2); `ReferenceSolution` product
  with time validity; `DriftModel` keyed to arc exposures and a time policy;
  `flux_status`, `telluric_status` fields; detector orientation and
  per-camera grouping in the model so DeepSpec needs no model change later.
- Ingest: arc sub-classification by lit traces; orientation at load.
- Wavelength: import of the notebook's laboratory solution as the first Mg b
  reference; the full-arc automatic solve from a previous reference (the notebook's
  grating-equation method); HighSpec `identify`; the drift fit per arc fiber,
  the across-detector transfer, the time policy and the quick-look re-reduce.
- Sky subtraction from the sky fiber of the same telescope, in wavelength
  space, with throughput scaling (one without a flat).
- Within-exposure combination across telescopes with the grouping, subset and
  rejection options of September §7.
- Pipeline order and the kept-versions product layout of §5; stage
  interfaces for flux calibration and telluric correction with pass-through
  implementations.
- Regression test against the notebook's laboratory numbers (skipped when
  `MASTPIPE_TEST_DATA` is unset).

**Plan 2 — DeepSpec and better extraction.** DeepSpec instrument definition,
`dispersion` model, Kr-then-ThAr `identify`, regression against the DeepSpec instrument team's
reference solutions, the DeepSpec science path on the star and supernova
frames, four-camera product grouping; Horne optimal extraction with flat
profiles; across-exposure combination (time grouping); the DeepSpec drift
question once answered.

**Plan 3 — Operations.** QA figures and pages, `watch` mode, idempotence
polish, focus recipe, PypeIt adapter, flux calibration and telluric recipes
once their observing policy is set.

## 8. Calls made while writing this revision

Settled here so the builder does not have to guess; each is recorded as a
decision and can be reversed by the user:

- **Plan 1 scope** is the full HighSpec path through sky subtraction and
  within-exposure combination (user decision 2026-10-29); DeepSpec moves to
  Plan 2 (user decision).
- **First HighSpec reference** is reproduced from the notebook's two stored anchor
  pairs through `identify` (nothing to import: the notebook keeps its solution
  in memory only); `identify` is in Plan 1 because the survey layout and the
  other gratings need it anyway.
- **Without a flat**, traces come from any fully lit frame and throughput
  defaults to one with a warning.
- **Instrument switches need not be observable**: drift pairing is nearest in
  time within the configuration, bounded by `drift.max_age`.
- **Reference validity is time-keyed**: an exposure uses the reference valid
  at its mid-time, with an explicit override, so re-reductions reproduce.
- **Quick-look re-reduces** `reference`-status exposures when the night's
  first drift arc arrives; the product key includes the drift-arc set.
- **Quick-look runs stages 1–4** and within-exposure combination; stages 5–6
  run only when their calibrations exist.
- **Kept versions** are columns in the per-fiber spectrum table on one
  wavelength grid (not separate extensions), to keep one row per fiber.
