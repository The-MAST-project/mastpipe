# Tutorials

Two notebooks, executed when this documentation is built, so every figure and
number comes from the code as it is. Both run on a **simulated** laboratory night
and need no data; swap the simulated folder for your own and they run on real
frames.

```{toctree}
:maxdepth: 1

api-walkthrough
inspecting-products
```

- **The pipeline step by step** walks through the Python API: from raw frames
  to a wavelength-calibrated, sky-subtracted, combined spectrum, one stage per
  section, with a look at each intermediate object. For the student who will
  extend the code.
- **Inspecting a night's products** opens the folder `quicklook` wrote, reads
  the products back as the pipeline's own objects and plots them, including the
  provenance, the drift and the kept versions. For the observer.

Run them yourself with Jupyter after `pip install -e ".[docs]"`; both are in
`docs/tutorials/` of the repository and are committed without outputs.
