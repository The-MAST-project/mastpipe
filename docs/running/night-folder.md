# The night folder

`quicklook` reads one folder of raw FITS frames. It does not care about file
names; it classifies every frame from its header, through the instrument's
header map, and groups the frames by camera and grating.

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
and the output amplifier, e.g. `image/em`) and temperature. Dark and background frames are stacked per exposure time.

## What the headers must say

The HighSpec header map reads the Andor Solis keywords and the ones the
control software will write:

- the exposure type from `IMAGETYP` or `OBSTYPE` (`bias`/`zero`, `dark`, `flat`,
  `arc`/`comp`, `light`/`object`/`science`, `background`); when it is missing, from rules on the exposure time
  (≈ 0 s → bias) and the lamp (`LAMP = thar` → arc, `flat` → flat, `off` with the
  shutter closed → dark), else science;
- the exposure time from `EXPOSURE`/`EXPTIME`, the start time from `DATE-OBS`
  (the header map declares it as the *start*; the frame's mid-time is derived by
  adding half the exposure), the temperature from `TEMP` (−999, Andor's
  "unstabilised", counts as missing), binning from `HBIN`/`VBIN`, readout mode
  from `READMODE`, the output amplifier from `OUTPTAMP` (`Conventional` →
  `conventional`, `Electron Multiplying` → `em`), which selects the gain and
  read noise;
- the grating from `GRATING` or the file name, the fiber layout preset from
  `FIBPRE`, each telescope's target from `OBJ01`…`OBJ20` when written, else one
  target for every fiber from `OBJECT`.

An arc is told to be *full* or *calibration-fiber* from which traces carry
light, so no keyword is needed for that. A frame the map cannot classify is
skipped and listed with the rule that matched (or `unclassified` when none did).

## The manifest

Frames whose headers lack the instrument keywords (laboratory data, or a
header the control software did not write) enter through a **manifest**: a
CSV keyed by file name whose non-empty cells override any metadata field.

```text
source,exposure_type,configuration,fiber_preset,arc_kind
image_bias.fits,bias,Mg,lab3,
background.fits,background,Mg,lab3,
flat_518.fits,flat,Mg,lab3,
FOCUS_05_13.30004_mm.fits,arc,Mg,lab3,full
```

```bash
mastpipe quicklook night/ --out out --instrument highspec --manifest night/manifest.csv
```

Fields you can set: `exposure_type` (`bias`, `dark`, `flat`, `arc`, `science`,
`background`, `unknown`), `arc_kind` (`full`, `calfiber`), `configuration`,
`fiber_preset`, `targets` (`mast01=HD 1234;mast02=HD 5678`, or a bare name for
every telescope; the older `target` column is still read that way),
`exposure_time`, `mid_time` (taken as the
mid-exposure time), `binning`, `readout_mode`, `amplifier`, `detector_temperature`,
`lamp_state`, `camera`, `instrument`, `classification_rule`; `extension`
selects an HDU of a multi-extension file.

## Layouts

The **fiber layout** says which trace is which: role (`arc`, `science`, `sky`),
telescope, nominal position. It is an input table, never inferred from the
data. HighSpec ships two presets:

`lab3`
: the three-fiber laboratory setup (all arc fibers), the default.

`survey45`
: the survey quasi-slit, 45 traces at an 8.5 px pitch: 5 arc, 20 science and
  20 sky, one science and one sky fiber per telescope. Choosing it also selects
  its narrower extraction and trace parameters, so no override is needed.

`quicklook` reads the preset from the frames' `FIBPRE` keyword and falls back
to `--layout`, then to the instrument's default; `identify` and `solve` use
`--layout` or the default only, so pass it for a survey arc.
Detections that match no layout row are dropped and listed in the trace QC; a
fiber missing on a night is simply absent from its products.
