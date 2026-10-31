# Adding an instrument

An instrument is a folder under `src/mastpipe/instruments/<name>/`. DeepSpec is
the next one and is the example here; nothing below needs code unless the
instrument needs a wavelength model the package lacks.

## `config.yaml`

```yaml
name: deepspec
cameras:
  g:
    name: g
    shape: [1027, 1056]        # raw ny, nx
    gain: 1.0                  # e-/ADU
    read_noise: 5.0            # e-
    pixel_size_um: 13.0
    saturation: 65000
    crop: null                 # or [x_min, x_max, y_min, y_max], end-exclusive
    dispersion_axis: 1
    transpose: true            # fibers are vertical stripes in the raw frame
    flip_x: false
    flip_y: false
params:                        # step parameters; every field has a default
  trace: {trace_halfwidth: 6, degree: 2}
  extract: {halfwidth: 7.0}
default_layout: lab4
layouts: [lab4, survey45]
layout_params:                 # per-preset overrides, merged when a layout is chosen
  survey45:
    trace: {trace_halfwidth: 4}
    lines: {halfwidth: 4}
    extract: {halfwidth: 4.0}
arc_classification: lit_traces # or `header` when a keyword says full/calfiber
```

Orientation (crop, then transpose, flip x, flip y) is applied once at load, so
every step sees the dispersion along the second axis. Several cameras are
several entries; the configuration blocks are shared by all cameras of the
instrument.

## `header_map.yaml`

The pipeline's vocabulary mapped to the instrument's raw keywords. Each field
lists candidate keywords in priority order and optional transforms:

```yaml
fields:
  exposure_time: {keys: [EXPTIME, EXPOSURE], type: float}
  mid_time: {keys: [DATE-OBS], type: str, reference: start}   # start | mid | end
  detector_temperature: {keys: [CCD-TEMP], type: float, missing_values: [-999]}
  exposure_type: {keys: [IMAGETYP], lower: true, map: {light: science, bias: bias}}
  lamp_state: {keys: [LAMP], lower: true}
  configuration: {keys: [BAND], lower: true, filename_regex: "_(?P<value>[ugri])_"}
  camera: {keys: [CAMERA], lower: true}
  fiber_preset: {keys: [FIBPRE]}
  target: {keys: [OBJECT]}
classification:                # applied in order when the header has no valid type
  - {when: {exposure_time: {lt: 0.001}}, type: bias}
  - {when: {lamp_state: {eq: thar}}, type: arc}
  - {when: {exposure_time: {missing: true}}, type: unknown}
  - {type: science}
```

Nothing downstream of ingest sees a raw keyword, so a change of control
software is a change of this file. Unresolved fields are recorded, never
guessed. The manifest can override any field frame by frame.

## Layouts

`layouts/<preset>.csv` with columns `trace, role, telescope, y_nominal,
slope_nominal, aperture`: which trace is which (`arc`, `science`, `sky`, or
`unknown` for laboratory traces whose identity is not known), which telescope,
its nominal row at the detector centre, its slope, and an optional per-fiber
extraction half-width. The number of traces is never inferred from data;
detections that match no row are dropped.

## Configuration blocks

`bands/<name>.yaml` for a camera-and-filter instrument, `gratings/<name>.yaml`
for a grating one; the registry merges both into `configurations`. The
`geometry.model` selects the wavelength model: `grating` (HighSpec's grating
equation) or `dispersion` (DeepSpec's `λ = p0 + p1/(x − p2)`, to be added with
the instrument). A line table per band (`lines/<name>.csv`) belongs to the block through
its `line_table` key; a second catalogue for a two-stage bootstrap (krypton,
then Th-Ar) will need a new block key when the DeepSpec model is added, since
the models reject unknown keys.

## Anchors

`anchors/<configuration>_<preset>.yaml` with `trace:` and `pairs:` once a first
reference has been made by hand with `identify`.

## Tests

The registry test proves the instrument loads and its layouts validate; a
simulated night gives the end-to-end test before real frames exist (today the
simulator writes HighSpec header cards and uses the grating model, so a new
instrument needs its own header cards and model there first);
real frames enter through `MASTPIPE_TEST_DATA` in `tests/regression/` and skip
without it.
