# Products policy: nothing thrown away

Two decisions from the team meeting shape every file the pipeline writes.

## Every individual spectrum is kept

A target observed by five telescopes has five spectra per exposure, and the
combined spectrum is a sixth file, not a replacement. Combination never
replaces its inputs; an excluded contributor is listed with its reason and
stays in the per-exposure file.

## Earlier stages are kept alongside later ones

A spectrum's product holds, per fiber and on one pixel grid, the raw
extraction (`FLUX`, `ERROR`, `MASK`) and, as separate column sets on the same
row, each later **version**: `SKYSUB`, and `FLUXCAL` and `TELLURIC` once those
stages have calibrations. A presence flag per row says which versions ran.
The combined file carries the same set, each combined from the corresponding
per-exposure version.

The case that motivated this: a supernova in a host galaxy whose sky or host
subtraction went wrong must be recoverable from the raw extraction without a
re-run.

The code enforces it: a `Spectrum`'s raw arrays are never replaced, stages add
versions, and the FITS writer stores them as columns so the table keeps one
row per fiber.

## Provenance

Every product header carries the pipeline version, the configuration hash,
the input frames' checksums, the calibrations used, the reference solution's
identifier and hash, the drift arcs' checksums, the layout preset, the three
statuses (wavelength, flux, telluric), and the steps applied. No raw
instrument keyword is copied into a product; the frame's metadata travels in
the pipeline's own vocabulary. A product can therefore be traced to the full
arc its wavelengths rest on, and two runs that agree on their keys produce the
same bytes.

## Idempotence

A science product's **key** is the input checksum, the configuration hash,
the reference identifier, the sorted drift-arc checksums and the layout
preset. A quick-look run skips a product whose key matches and redoes one
whose key changed — which is how early-night frames reduced before the first
calibration-fiber arc are upgraded from `reference` to `transferred` once the
arc lands. Full arcs are idempotent through the reference that lists them as
its source.

## Names

```text
{instrument}_{configuration}_{camera}_{YYYYmmddTHHMMSS}_{kind}.fits      masters, traces, spec, comb_<target>
{instrument}_{configuration}_{preset}_{YYYYmmddTHHMMSS}_wref.fits        reference solutions
```

Times are the exposure's UTC mid-time, derived by ingest from the header's
start time and the exposure time. Raw file names never appear in product
names.
