# Running the pipeline

These pages are for the person who has a folder of frames and wants spectra.
They follow the order you will meet things in: install, arrange the night,
make the first wavelength reference for a grating, run the quick-look, read
what came out, and sort out what did not.

```{toctree}
:maxdepth: 1

install
night-folder
first-reference
quicklook
products
troubleshooting
```

## The short version

```bash
pip install "git+<repository>"
# once per grating and fiber layout: the first reference solution from a full arc
mastpipe identify night/arc_full.fits --out out --instrument highspec --configuration Mg --stored-anchors
# every night, from the site's folder (the control software names no grating)
mastpipe quicklook <night>/highspec --out out --instrument highspec --configuration Mg --layout survey45
mastpipe show out/spec/highspec_Mg_main_20261001T221500_3f9a1c_spec.fits
```

A laboratory or simulated night needs neither `--configuration` nor `--layout`
only when its headers (or file names) name the grating and the default `lab3`
layout applies; otherwise give both, or a manifest.

Three words carry most of the meaning:

- a **reference solution** is a per-fiber wavelength solution made from a *full
  arc*, an exposure in which every fiber sees the lamp; it is made every few
  months and stored;
- a **drift model** is the shift and stretch of the spectrum measured on a
  *calibration-fiber arc*, an exposure one to three times a night in which only
  the five calibration fibers see the lamp;
- a science frame's **wavelength status** says what it got: `transferred`
  (reference corrected by drift), `reference` (reference alone, no arc paired),
  or `nominal` (the configured geometry, no reference at all).
