# The data model

Everything the layers exchange is a dataclass in `mastpipe.core.model`.
Arrays are numpy; masks are boolean with **True = bad**. Each class validates
its shapes and vocabularies on construction, so a wrong-length array or an
unknown status fails at the boundary, not three steps later.

## Frames

`FrameMeta`
: One raw frame described in the pipeline's vocabulary: instrument, camera,
  exposure type, configuration, exposure time, mid-time (UTC, derived from the
  header's start time), binning, readout mode, temperature, lamp state, fiber
  preset, targets (`(telescope, target)` pairs, or one `("*", target)` pair for
  every fiber; `target_for(telescope)` reads them), arc kind (`full`, `calfiber`), checksum, source file name,
  extension, the fields left unresolved, and the classification rule that
  decided the type. Frozen.

`FrameRef`
: Where a frame lives plus its `FrameMeta`; pixels are loaded on demand.

`Frame`
: Data, optional variance and mask, and the meta. After `calibrate_frame` the
  variance and mask are set.

`MasterCalib`, `CalibrationSet`
: A stacked master with validity keys and provenance; the set of masters
  applied to a frame (bias, dark rate, background with its exposure time, gain,
  read noise).

## Traces and spectra

`Trace`, `TraceSet`
: A fitted trace (polynomial, centre, slope, rms, centroid samples with their
  clipping mask, role, telescope, layout index, per-fiber aperture) and the set
  of them with the layout shift, the missing and the extra traces. After
  layout matching `Trace.index` **is** the layout's trace id.

`Spectrum`
: One fiber: pixel grid, raw `flux`/`error`/`mask` (never replaced), role,
  telescope, target, optional `wavelength`, `sky`, `continuum`, and
  `versions`: a dict from `"skysub"`, `"fluxcal"`, `"telluric"` to
  `SpectrumVersion(flux, error, mask)`. `latest(upto)` returns the most
  advanced version present.

`ExposureSpectra`
: All spectra of one exposure with its meta, the three statuses, the
  reference identifier and hash, and the drift arcs used.

`CombinedSpectrum`
: The target's spectra on a common grid, version by version, with contributor
  and exclusion tables, the grouping rule and the barycentric terms.

## Wavelengths

`WavelengthSolution`
: Per fiber: model kind and parameters, pixel count, the identified lines
  (`x, wavelength, species, dx_px`), rms, line count, the anchor range in
  pixels and wavelength, the trace centre, and the geometry constants.

`ReferenceSolution`
: The solutions of one full arc for one instrument, configuration and fiber
  preset, with `created`, `valid_from`, the source checksums, an identifier
  (`{instrument}_{configuration}_{preset}_{YYYYmmddTHHMMSS}`) and a hash over
  the solutions.

`DriftModel`
: One calibration-fiber arc's drift: the arcs' checksums and mid-time, a
  per-fiber table (fiber, trace centre, shift, stretch, lines, rms, usable),
  the smoothed shift and stretch polynomials across the detector, their
  scatter, a usable flag and the pixel count. `shift_at(y)` and `stretch_at(y)`
  evaluate them.

`QCReport`
: A stage's scalar quality numbers, written as JSON.

## Vocabularies

```python
EXPOSURE_TYPES = ("bias", "dark", "flat", "arc", "science", "background", "unknown")
ARC_KINDS = ("full", "calfiber", "unknown")
ROLES = ("arc", "science", "sky", "unknown")
WAVELENGTH_STATUSES = ("none", "nominal", "reference", "transferred", "solved")
FLUX_STATUSES = ("uncalibrated", "calibrated")
TELLURIC_STATUSES = ("uncorrected", "corrected")
VERSIONS = ("skysub", "fluxcal", "telluric")
```

## Configuration models

`mastpipe.config.models` holds the Pydantic models, all frozen and forbidding
unknown keys: `DetectorConfig`, `GeometryConfig` (model kind `grating` or
`dispersion`, groove spacing, incidence angle, order, camera focal length,
central wavelength, dispersion sign), one `*Params` model per step
(`CalibParams`, `TraceParams`, `ExtractParams`, `LineParams`, `WavecalParams`,
`DriftParams`, `SkyParams`, `CombineParams`) gathered in `StepParams`,
`ConfigurationBlock`, `InstrumentConfig` and the merged `EffectiveConfig`.

## Wavelength models

`mastpipe.steps.wavemodel.WavelengthModel` is the protocol every model obeys:
`wavelength(x)`, `pixel(λ)`, `dispersion(x)`, a parameter vector, a kind, and a
dict round trip. `GratingModel` implements it for HighSpec (the grating
equation with a distortion polynomial). The DeepSpec `dispersion` kind has
its slot and raises until the next release. `model_from_solution` rebuilds a
model from a stored `WavelengthSolution` and the grating geometry.
