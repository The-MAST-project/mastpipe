# Products

Everything is written under the output folder. File names come from the
frame's metadata, never from the raw file name:

```text
{instrument}_{configuration}_{camera}_{YYYYmmddTHHMMSS}_{checksum}_{kind}.fits
```

with the UTC mid-time of the exposure (`notime` without one), `checksum` the
first six hex digits of the input file's SHA-256 (for a master, its first input
frame's), so two exposures with the same mid-time never share a name (`nosum`
on a trace set built from a library flat, which has no input frame), and
`kind` one of `bias`, `dark`,
`background`, `flat`, `traces`, `spec`, `comb_<target>`. A frame whose header
names no grating takes the configuration given with `--configuration`. Reference solutions
are `{instrument}_{configuration}_{preset}_{stamp}_wref.fits`, the stamp being
the full arc's mid-time; the name without `_wref.fits` is the **reference
identifier** quoted in every product that used it.

```text
out/
  calib/   masters, trace sets, references, fiber_agreement.csv, index.json, index.lock
  spec/    *_spec.fits (one per science exposure), *_comb_<target>.fits (one per target)
  qa/      *.json: one per trace set, reference, drift arc and science frame
  logs/    quicklook_<UTC stamp>.log
```

Every product, the QC files and `calib/index.json` are written to a temporary
`<name>.tmp-<pid>` next to the target and moved into place when complete, so an
interrupted run never leaves a partial file. `calib/index.lock` is the lock that
lets several runs (a quick-look and an `identify`) share one `calib/` folder
without losing each other's index entries; leave it in place.

## Output folders of an earlier version

Product names gained the checksum segment and the grating's new name
(`Mg`, not `mgb518`), so a product written by Plan 1 is never found at its new
path. It is not reduced again in place: a new file is written beside it, and
the old `spec/` and `qa/` files stay, so a glob over the folder sees both.
`quicklook` warns once per run when `spec/` holds product names without the
checksum segment. Use a fresh `--out` for such a folder, or clear its `spec/`
and `qa/` folders first. The calibration library in `calib/` may be kept: its
old entries stay in the index but never match (their readout is keyed
`…/unknown`), and the masters are rebuilt beside them.

## Looking inside

```bash
mastpipe show out/spec/highspec_Mg_main_20261001T221500_3f9a1c_spec.fits
mastpipe show out/calib/highspec_Mg_lab3_20260909T113702_wref.fits
```

`show` prints the product kind, the frame and its statuses, a per-fiber
summary (role, telescope, target, median S/N, masked fraction, kept versions,
drift), or the per-fiber solutions of a reference, the traces of a trace set,
the contributors of a combined spectrum or the validity of a master, then the
provenance. The
same information is in the FITS headers and tables, so any FITS reader works.

## Provenance in every header

Every product carries `MP*` cards: the pipeline version (`MPVER`), the
configuration hash (`MPCFGH`), the wavelength, flux and telluric statuses
(`MPWSTAT`, `MPFSTAT`, `MPTSTAT`), the reference identifier and hash
(`MPWREF`, `MPWREFH`), the layout preset (`MPLAYOUT`), and numbered lists of
input checksums, calibrations used, drift arcs used and steps applied. No raw
instrument keyword reaches a product header; the frame's metadata travels as
numbered name/value pairs.

## The spectra table

A science product has one table, `SPECTRA`, with **one row per fiber** on a
common pixel grid.

| columns | content |
| --- | --- |
| `FIBER`, `ROLE`, `TELESCOPE`, `TARGET` | the fiber's identity from the layout; `FIBER` is the layout's trace id |
| `PIXEL`, `FLUX`, `ERROR`, `MASK` | the raw aperture extraction, never overwritten |
| `WAVELENGTH` (`HAS_WAVELENGTH`) | the fiber's wavelengths by the time policy |
| `SKY` | the sky spectrum resampled onto the fiber (`CONTINUUM` is reserved and not populated in this release) |
| `FLUX_SKYSUB`, `ERROR_SKYSUB`, `MASK_SKYSUB` (`HAS_SKYSUB`) | the sky-subtracted version |
| `FLUX_FLUXCAL`, … (`HAS_FLUXCAL`); `FLUX_TELLURIC`, … (`HAS_TELLURIC`) | the later versions, once those stages have calibrations |

Each **kept version** is a set of three columns on the same row. The columns
exist for every version; they are populated only for fibers whose stage ran
(`HAS_<version>` per row, and a presence card on the table), NaN and masked
otherwise. The raw extraction always stays: if the sky or
host subtraction of a transient went wrong, the un-subtracted spectrum is in
the same row, no re-run needed.

A `DRIFT` table is attached when a drift model was applied: per calibration
fiber, the shift and stretch measured, the lines used, the rms, a usable flag
and the blend weight, with the arc's mid-time, the smoothed model's
coefficients, the number of usable arc fibers (`MPDNUSE`) and whether the model
is usable (`MPDUSE`) in its header.

## Wavelength status

| status | meaning |
| --- | --- |
| `transferred` | reference solution corrected by the drift of the paired calibration-fiber arcs |
| `reference` | reference solution alone: no usable calibration-fiber arc of the same grating, fitted against this reference, within `max_age_hours` of the exposure |
| `nominal` | no reference solution exists; wavelengths from the configured grating geometry only |
| `solved` | a reference product: solved from a full arc |

The flux and telluric statuses are `uncalibrated` / `uncorrected` in this
release.

## Combined products

One file per target and exposure: table `COMBINED` with `WAVELENGTH` on a
common grid (log-uniform by default) and `FLUX_<V>`, `ERROR_<V>`, `MASK_<V>`
for `RAW` and every version any contributor carries, each combined from the
contributors that carry it; `CONTRIB` (fiber, telescope, weight) and
`EXCLUDED` (fiber, telescope, reason) tables; the grouping rule in the header,
and `MPVMISS` listing, per version, the contributing fibers without it. Contributors are weighted by inverse variance
with a one-pass outlier clip. Excluded spectra are listed with the reason
(configuration, no wavelength, low S/N, flux ratio) and remain in the
per-exposure file.

## Reference products

`SOLUTIONS`: one row per fiber with the model kind and parameters, the rms and
number of lines, the pixel and wavelength range the lines span (outside it the
solution is extrapolated), the trace centre it was measured on and the
geometry constants. `LINES`: every identified line with its catalogue
wavelength, species and residual. The header carries the identifier, the
hash, the creation and validity times.

## QC files

One JSON per trace set (traces found, matched, missing, extra with positions,
per-trace rms and S/N, layout shift), per reference (per-fiber rms and line
counts), per drift arc (per-fiber shift, stretch, lines, rms; the model's
scatter), and per science frame (per-fiber median S/N and masked fraction,
drift shift and stretch used, the drift arcs and their blend weights, sky
subtraction reasons, `wavelength_skipped` fibers). They are plain JSON for a
dashboard to read.
