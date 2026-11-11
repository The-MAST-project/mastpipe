# The first reference solution

Before `quicklook` can give wavelengths, each grating and fiber layout needs a
**reference solution**: an independent wavelength solution per fiber, made
from a *full arc* in which every fiber sees the lamp. It is the rare, attended
step of the pipeline. Once a reference exists, later full arcs are solved
automatically from it (`solve`, or `quicklook` when it meets one), and the
nightly calibration-fiber arcs only measure the drift against it.

The reference is valid from its arc's mid-time until the next full arc of the
same grating and layout, so an old night re-reduced later still uses the
reference that was current then.

## With stored anchors

The HighSpec Mg b grating ships with the two anchor pairs read off the laboratory
notebook for the three-fiber laboratory layout. If your full arc matches that
setup, one command does it:

```bash
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg \
    --stored-anchors
#   reference highspec_Mg_lab3_20260909T113659: 3 fibers solved
#   written: out/calib/highspec_Mg_lab3_20260909T113659_wref.fits
```

The masters of the calibration library (`--library`, default `out/calib`)
that match the arc's camera, binning, readout and temperature (the flat its
configuration too) are applied, the nearest in time, and the log names them
(`masters applied: ...`); a master that matches in all but temperature (more
than 2 °C away) is named in a warning rather than left out silently. Traces
come from the flat master when there is one, otherwise from the arc itself.

`identify` and `solve` use the instrument's default layout (`lab3` for
HighSpec) unless told otherwise: for a survey full arc pass `--layout survey45`.
No header card names the layout; `quicklook` also reads a manifest's
`fiber_preset` column. The
reference is keyed by grating *and* layout, so a reference made with the wrong
preset is not found by `quicklook` later. The
reference is registered in the library and `fiber_agreement.csv` beside it
lists, for every catalogue line identified in two or more fibers, how far the
independent solutions disagree. In the laboratory the median spread is 0.02 px.

## Reading anchors off the panels

For a new grating, or a layout whose stored anchors do not fit, run without
anchors first:

```bash
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg
#   anchor panels: out/anchor_panels.png
#   candidates: out/anchor_candidates.csv
```

```{figure} ../_static/anchor_panels.png
:alt: Anchor panels: the Th-Ar template above, the measured arc spectrum below, strongest lines numbered in each
:width: 100%

The anchor panels for a simulated full arc. Template pair 2 and measured pair 3 are
the same line, as are 4 and 6.
```

The figure shows the Th-Ar template (the Murphy et al. atlas convolved to the
instrument's resolution when `MASTPIPE_ATLAS_DIR` or `--atlas-dir` points at it,
else a synthetic template from the shipped line table) above the measured
spectrum of one trace over the same band, with the strongest well-separated
lines numbered in each. Two well-separated pairs that are clearly the same line
fix the two free parameters of the grating model. Read them off and pass them
back:

```bash
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg \
    --anchors 515.8605:267.22,525.836:1423.22 --trace 0
```

or from a YAML file (`--anchors-file`) with `pairs: [[nm, px], ...]` and
`trace:`. `--trace` is the layout trace id of the fiber the pairs were measured
on; the other fibers are seeded from it by cross-correlation and then solved on
their own.

If the patterns in the two panels do not look like the same spectrum, the
template window is off (wrong configuration) or wavelength runs the other way
across the detector (check the grating block's `dispersion_sign`).

## How a solution is judged

For each fiber the command reports the number of lines used and the rms in
pixels; the reference product stores the identified lines with their
residuals. On the laboratory frame the three fibers solve with 23, 25 and 20
lines at 0.016, 0.033 and 0.012 px rms, against the notebook's 23, 25 and 19
lines at 0.016, 0.038 and 0.010 px. A fiber with
fewer than `wavecal.min_lines` matches (default 6) is not solved and is absent
from the reference; `quicklook` then gives that fiber no wavelength (its
`HAS_WAVELENGTH` flag stays false; nothing else is reported for it).

## Later full arcs

```bash
mastpipe solve night2/arc_full.fits --out out --instrument highspec --configuration Mg
#   reference highspec_Mg_lab3_20261112T010203: 3 fibers solved, seeded by highspec_Mg_lab3_20260909T113659
```

`solve` takes the latest reference of the grating and layout from the library
as its starting point, solves every fiber afresh, and adds the new reference
without replacing the old one. It exits with status 2 and "run identify first"
when the library has none. `quicklook` does the same when it finds a full arc
in a night folder.
