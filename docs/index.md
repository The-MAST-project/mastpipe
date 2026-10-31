# mastpipe

**mastpipe** is the data-reduction pipeline for the MAST fiber-fed spectrographs:
**HighSpec** (R ≈ 20 000, narrow-band gratings, one detector) and **DeepSpec**
(R ≈ 500, four cameras). It takes raw frames to wavelength-calibrated,
sky-subtracted spectra per fiber, combined spectra per target, and
quality-assurance products. It runs unattended at the observatory as a
quick-look and on a laptop as a manual reduction tool.

The reduction steps descend from the HighSpec laboratory
reduction notebook, written as pure functions behind a declarative instrument
layer. This release covers HighSpec end to end; DeepSpec follows.

::::{grid} 1 1 3 3
:gutter: 2

:::{grid-item-card} Running the pipeline
:link: running/index
:link-type: doc
Install, prepare a night folder, make the first reference solution, run
`quicklook`, read the products. For observers.
:::

:::{grid-item-card} How it works
:link: how-it-works/index
:link-type: doc
The stages, the three-cadence wavelength calibration, what is inherited from
the notebook and what is deferred. For the team.
:::

:::{grid-item-card} Extending the pipeline
:link: extending/index
:link-type: doc
The architecture and data model, adding a grating or an instrument, the
simulator and tests, the API reference. For developers.
:::

::::

## In one picture

A night's frames become products in six stages, with combination at the end:

1. **Clean** — bias, dark or background, variance, saturation and cosmic-ray masks.
2. **Trace and extract** — traces matched to the fiber layout, aperture extraction.
3. **Wavelength** — the stored reference solution of each fiber, corrected by the
   drift measured on the calibration-fiber arcs nearest in time.
4. **Background** — the telescope's sky fiber, scaled and resampled, subtracted.
5. **Flux calibration** and 6. **telluric correction** — interfaces whose methods
   are still open; they pass spectra through unchanged for now.

Nothing is thrown away: every fiber's raw extraction stays in the product next to
each later version.

```{toctree}
:hidden:
:maxdepth: 2

running/index
how-it-works/index
extending/index
api/index
changelog
```
