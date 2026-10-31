# How it works

These pages describe the pipeline as a method: what each stage does, which
choices were made and why, and what is still open. They are written from the
code and the decision log; the design documents they summarise are kept with
the project.

```{toctree}
:maxdepth: 1

overview
wavelengths
traces-and-extraction
sky-and-combination
products-policy
limits
```

## The instrument in two sentences

MAST feeds one spectrograph at a time from twenty telescopes through fibers
that form a quasi-slit on the detector: for each telescope a science fiber and
a sky fiber, plus five calibration fibers that only the calibration unit in
the instrument room can light. All telescopes expose together, each on its
own target, so one frame holds twenty independent spectra that are combined
only after every fiber has its own wavelength solution.

HighSpec is the first instrument: one detector, straight nearly parallel
traces, narrow-band gratings (about 14 nm at R ≈ 20 000) with Th-Ar arcs. Its
grating is a configuration file, so adding one needs no code.
