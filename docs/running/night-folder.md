# The night folder

`quicklook` reads one night folder of raw FITS frames. It does not care about
file names; it classifies every frame from its header, through the
instrument's header map, and groups the frames by camera and grating.

## The folder shape

At the observatory the spectrograph control software (`MAST_spec`) files every
HighSpec exposure in an acquisition folder of its own:

```text
<night>/
  highspec/
    seq.txt
    acquisition-0001/
      exposure-001.fits
      exposure-002.fits
    acquisition-0002/
      exposure-001.fits
    ...
```

Point `quicklook` at `<night>/highspec/`. Pointing it at `<night>/` works too:
the unit cameras and DeepSpec file their frames under the same night, and
ingest leaves out every frame whose `INSTRUME` card names another instrument
than the one requested (each such frame is listed under *skipped* with the
reason). A frame with no `INSTRUME` card, such as a laboratory frame, is
accepted. Where `MAST_spec` will file the HighSpec biases, darks and flats is
not settled yet; any acquisition folder under the root is read.

The folder is scanned **recursively**. A FITS file is any file whose name ends
in `.fits`, `.fit` or `.fits.fz`, in any case; files and folders whose name
starts with `.` are ignored, and other files (`seq.txt`, logs, notes) are not
scanned at all. A symbolic link to a folder is not descended into (a linked
FITS file is read). A file with a FITS name that cannot be opened is recorded as a
**failure** (`unreadable FITS file`), and the run's exit status is 2.

Every frame is known by its path **relative to the night folder**, with forward
slashes: `acquisition-0012/exposure-003.fits` and
`acquisition-0013/exposure-003.fits` are two frames. That relative path is
what the logs, the QC files and the products' provenance call the frame's
*source*.

## What a night needs

| frames | role |
| --- | --- |
| bias (several) | the pedestal; stacked into a master |
| dark, or a **background** exposure of the science exposure time | the dark current; the notebook's method uses a lamp-off background of matched exposure |
| flat (several) | where the traces are, and the relative throughput of the fibers |
| **full arc** | every fiber lit by the lamp; needed only when a new reference solution is to be made, every few months |
| **calibration-fiber arc** (one to three per night) | the five calibration fibers lit; measures the drift since the full arc |
| science | the targets; no arc light |

A kind the night lacks (say, no bias was taken) is taken from the calibration
library in `out/calib`, matched on camera, binning, readout (the readout mode
and the output amplifier, e.g. `image/em`) and temperature. Dark and
background frames are stacked per exposure time.

A group follows the binning and readout of its science frames and of its
not-yet-resolved `light` frames (the most common among them; all its frames when
it has neither). A calibration frame (bias, dark, background, flat, or a full
arc that is a trace candidate) read out otherwise is left out and listed under
*skipped* with the reason (`readout image/em is not the group's
image/conventional`); it belongs to another set and does not change the exit
status. A calibration frame binned otherwise, and a science, `light` or arc
frame binned or read out otherwise, is dropped and recorded as a failure with
the reason. A night that mixes frames with and without the `READMODE` and
`OUTPTAMP` cards therefore loses the frames without them (their readout is
`unknown/unknown`) whenever the group's own readout is known; give such frames
their readout through the manifest (`readout_mode`, `amplifier`) to keep them.

## What the headers must say

The HighSpec header map reads the Andor Solis keywords and the ones the
control software writes or will write.

### Exposure type

The type comes from `IMAGETYP` or `OBSTYPE`: `bias`/`zero`, `dark`, `flat`,
`arc`/`comp`, `object`/`science`, `background`, and `light`.

`light` is what the control software writes for every shutter-open exposure,
arcs and science frames alike, so it is a **provisional** type: "shutter open,
source not stated". It is resolved in three ways, in order:

1. **By rule, at ingest.** The header map's classification rules that carry a
   condition refine a `light` frame: an exposure time below 1 ms makes it a
   bias, a lamp card (`LAMP = thar`, `flat`, or `off` with the shutter closed)
   makes it an arc, a flat or a dark, and a missing exposure time makes it
   `unknown`. A frame with no valid type at all goes through the same rules and
   falls back to `science`.
2. **By content, in the runner.** Every frame still `light` after ingest is
   calibrated with the night's masters and collapsed along the slit into one
   spectrum before the traces are found (on a night without a flat the full
   arc is itself a `light` frame and is the trace source). A spectrum whose
   fiber light sits mostly in narrow peaks makes the frame an **arc**
   (`content:arc_like`); a continuum makes it **science**
   (`content:continuum`); a frame with no signal above the noise is skipped
   (`content:no_signal`). The rule is written in the frame's QC and in the
   product's provenance.
3. **By the manifest**, which always wins: a frame whose manifest row sets
   `exposure_type` takes that type; a row may set it to `light`, and the frame
   is then resolved from its content as above.

Content cannot tell a flat from a science frame (both are continua), and a
lamp-off background with a little ambient light in it may resolve as `science`.
Declare flats and backgrounds with a header type (`flat`, `background`) or a
manifest row. Whether an arc is a *full* or a *calibration-fiber* arc is told
afterwards from which traces carry light, so no keyword is needed for that.
A frame that ends up `unknown` is skipped and listed with the rule that
decided it (or `unclassified` when none did).

### Grating

The grating comes from the `GRATING` card, spelled as the MAST fleet spells it
(`Mg`, `Halpha`, `Ca`), or from a laboratory file name carrying the old alias
(`mgb518`, `halpha`, `caiihk`), which the header map translates. When neither
names it, pass it on the command line:

```bash
mastpipe quicklook <night>/highspec --out out --instrument highspec --configuration Mg
```

There is no silent default: without `--configuration`, frames whose grating is
unknown are not reduced, the run records a failure that lists the
configurations the instrument defines (`['Mg']` for HighSpec today, not every
fleet grating name), and the exit status is 2. The control software writes no `GRATING` card today,
so a night from the site needs the option.

### Time

The exposure's time comes from `DATE-OBS`, else `FRAME` (the camera driver's
start of frame exposure, with milliseconds), else `DATE` (file creation). The
header map declares the instant the card records (`reference: start`), and
ingest adds half the exposure time to reach the mid-time.

Each time field carries a `timezone` option saying how a stamp without an
offset is read: `utc`, an IANA zone name (`Asia/Jerusalem`), or `assume-utc`.
An IANA name of UTC itself (`UTC`, `Etc/UTC`) is the same as `utc`, and its
frames are labelled `utc`.
A stamp with an explicit offset is honoured whatever the option says. Under
`assume-utc` a naive stamp is taken as UTC, the run logs one warning, and
every product records the assumption in its `MPTZ` card (`assumed-utc`).
HighSpec ships with `assume-utc` until it is established whether the control
computer's clock keeps UTC; the option then changes in one line of
`header_map.yaml`. Mid-times are stored and printed as ISO UTC without a zone
suffix (`2026-10-01T18:02:30`).

### Everything else

- the exposure time from `EXPOSURE` or `EXPTIME`;
- the temperature from `TEMP` or `CCD-TEMP` (−999, Andor's "unstabilised",
  counts as missing);
- the binning from `HBIN`/`VBIN`;
- the readout from `READMODE` and the output amplifier from `OUTPTAMP`
  (`Conventional` → `conventional`, `Electron Multiplying` → `em`), which
  together select the gain and read noise and are matched by the library;
- the targets, one per telescope: the cards `OBJ01`…`OBJ20` when the control
  software writes them (none does yet), else one `OBJECT` card that applies to
  every science fiber. Only a missing card, a blank one or a logical one
  (`OBJECT = F`) names no target (`OBJECT = 0` is the target `0`); a
  calibration frame (bias, dark, flat, arc, background) does not list its
  targets as unresolved, and any other frame without a target does, whether
  its type comes from the header or from a manifest row.

## The manifest

Frames whose headers lack the instrument keywords (laboratory data, or a
header the control software did not write) enter through a **manifest**: a
CSV whose rows name frames by `source` and whose non-empty cells override any
metadata field.

```text
source,exposure_type,configuration,fiber_preset,arc_kind
image_bias.fits,bias,Mg,lab3,
background.fits,background,Mg,lab3,
flat_518.fits,flat,Mg,lab3,
FOCUS_*.fits,,Mg,lab3,
```

```bash
mastpipe quicklook night/ --out out --instrument highspec --manifest night/manifest.csv
```

A row's `source` is one of:

- a **relative path** from the night folder (`acquisition-0012/exposure-003.fits`);
- a **glob pattern** over relative paths (`acquisition-0012/*.fits`;
  `*` and `?` stay within one folder, `**` crosses folders; `[0-9]` and `[!0-9]`
  are character classes, read as Python's `fnmatch` reads them, and a class
  that cannot be read, such as the reversed range `[9-0]`, is an error naming
  the row);
- a **bare file name**, which matches the one file of that name anywhere in the
  tree; a name that several files share is an error that lists them.

A `\` in a source is read as a folder separator, so a manifest written on
Windows (`acquisition-0012\exposure-003.fits`) matches as well.

An `extension` column restricts a row to one HDU of a multi-extension file.
Rows apply in order, so a later row wins over an earlier one for the frames
both name. A row that matches no frame is listed under *skipped* as
`manifest row 3 (exposure-009.fits): matches no frame`, never dropped
silently; a row naming a file that exists but is not a FITS file (a `.txt`, a
`.fts`) is listed there as `not a FITS file`, since the scan never reads it.

Fields you can set: `exposure_type` (`bias`, `dark`, `flat`, `arc`, `light`, `science`,
`background`, `unknown`), `arc_kind` (`full`, `calfiber`), `configuration`,
`fiber_preset`, `targets` (`mast01=HD 1234;mast02=HD 5678`, or a bare name for
every telescope; the older `target` column is still read that way; `;`
separates entries and the first `=` separates the unit from the target, so no
target name can contain `;`, and a bare name cannot contain `=`),
`exposure_time`, `mid_time` (taken as the mid-exposure time; a stamp with an
explicit offset is honoured and converted to UTC, a naive one is UTC),
`binning`, `readout_mode`, `amplifier`, `detector_temperature`, `lamp_state`,
`camera`, `instrument`, `classification_rule`.

An error in the manifest itself (a path that cannot be read, such as a folder; a file
that is not UTF-8 CSV; rows with more cells than the header, as a trailing comma on every
row gives; a missing `source` column or an unknown one; a bad glob; an `extension` that
is not a number; a value that does not parse, such as `targets` naming something that is
not a telescope unit) stops the run before any frame is reduced: `quicklook` prints one
line on stderr naming the file or the row, e.g. `FAILED: manifest row 2: bad glob
'frame[9-0].fits': …`, and exits with status 2. A manifest that does not exist is
reported the same way.

## Layouts

The **fiber layout** says which trace is which: role (`arc`, `science`, `sky`),
telescope (`mast01`…`mast20`), nominal position. It is an input table, never
inferred from the data, and no header card selects it (in MAST a "fiber
preset" is a stage position, not a layout). HighSpec ships two presets:

`lab3`
: the three-fiber laboratory setup (all arc fibers), the default.

`survey45`
: the survey quasi-slit, 45 traces at an 8.5 px pitch: 5 arc, 20 science and
  20 sky, one science and one sky fiber per telescope. Choosing it also selects
  its narrower extraction and trace parameters, so no override is needed. Which
  physical trace belongs to which telescope is a placeholder until the
  instrument team confirms it.

`quicklook` takes the preset from `--layout`, else from a manifest's
`fiber_preset` column, else the instrument's default; `identify` and `solve`
use `--layout` or the default, so pass it for a survey arc. Detections that
match no layout row are dropped and listed in the trace QC; a fiber missing on
a night is simply absent from its products.

## Leftovers of an interrupted run

Every product, master, reference, index and QC file is written under a
temporary name in its own folder (`<name>.tmp-<pid>-<random>`) and moved into place in
one step, so a reader never sees a partial file. A run that is killed may
leave such a temporary behind; it is never read and is safe to delete.
