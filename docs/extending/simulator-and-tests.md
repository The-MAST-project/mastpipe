# The simulator and the tests

Survey data did not exist when the pipeline was written, so every step and the
end-to-end path are tested against a **simulator with known truth**. It is the
first thing to reach for when changing a step: simulate a night, run the
change, compare with the truth.

```{figure} ../_static/simulated_frame.png
:alt: A simulated HighSpec science frame with 45 traces along the dispersion axis
:width: 100%

A simulated science frame on the survey layout: 20 science traces with absorption dips, 20 sky traces with emission lines, 5 dark calibration-fiber traces, cosmic rays.
```

## What `simulate_night` writes

`mastpipe.sim.simulator.simulate_night(effective, layout, line_table, out_dir, …)`
writes a night of raw-shaped FITS frames whose headers the HighSpec header map
classifies, plus `truth.npz` and `truth.json`:

- bias, dark, flat and background frames with the detector's noise model;
- one **full arc**: every trace lit, the line table's catalogue lines placed by
  the exact grating equation at the configured geometry, with per-fiber
  throughput; every fiber has its own small perturbation of the dispersion
  solution (`dispersion_jitter_px`, 0.3 px by default: the central wavelength
  moves by up to that much and the chip edges by up to that much again), so no
  two fibers share a pixel-to-wavelength map;
- **calibration-fiber arcs**: only the arc traces lit, the lines shifted and
  stretched by a linear drift law that is zero at the full arc and grows to the
  configured end value at the last frame; the shift also varies linearly with
  the trace's row, by up to `drift_gradient_px` (0.05 px by default) between
  the bottom and the top of the detector at the last frame;
- **science frames**: a continuum with absorption dips in the science fibers,
  a continuum with emission lines in the sky fibers, no arc light, the same
  drift law, cosmic rays on single pixels;
- a frame order in time (bias, dark, flat, full arc, background, then
  alternating arc and science groups at a cadence) so every science frame has
  a calibration-fiber arc before it; `no_cal_arcs=True` writes none;
- `swap_fibers=(a, b)` exchanges what two fibers carry on every frame while the
  truth keeps the unswapped assignment: a pipeline that keys fibers correctly
  then disagrees with the truth on exactly those two, which is how the
  end-to-end test proves the truth comparison can see a fiber-identity bug.

`SimTruth` records the trace polynomials, the nominal and the per-fiber
dispersion solutions, every line's pixel per fiber with an on-chip mask, the
throughputs, the science and sky spectra per fiber, the noise constants, and a
per-frame table with the drift (shift, stretch, row gradient) at each frame's
time; `trace_shifts(file)` gives every fiber's shift on a frame. Compare a
fiber's wavelengths with its own solution (`x_c_per_trace`, `B_per_trace`),
not the nominal one. Everything is deterministic from a seed.

## How the tests use it

`tests/conftest.py` builds one simulated survey night and one laboratory-layout
night per session, plus the calibrated full arc, its line tables and fluxes,
so a test can start at any stage without redoing the earlier ones. Typical
assertions:

- traces within 0.3 px of the truth polynomials at three columns, all 45 found;
- extracted science flux within 2 % of the throughput-scaled truth;
- 45/45 fibers solved from a single anchored model, every wavelength within
  0.3 dispersions of truth; a reference-seeded re-solve equal to the first;
- the drift fit recovering the planted shift to 0.05 px, the transfer landing
  drifted absorption dips within 0.3 dispersions;
- sky lines removed to within 3σ plus the loss of linearly interpolating a
  3 px line between two fibers' grids; combination reducing the noise of four
  copies by two;
- the end-to-end quick-look: every science product `transferred`, wavelengths
  against truth, sky-subtracted versions present, one combined file per
  target, the idempotence and re-reduce behaviour, the `interpolate` blend;
- a night with one fiber darkened: the other 44 keep their own solutions;
- a night with two fibers swapped: exactly those two fail the truth comparison.

## Tiers and markers

All from the root of the clone, where `pyproject.toml` configures pytest:

```bash
pytest -m "not slow"   # unit tests on small arrays and the shared fixtures, about three minutes
pytest                 # plus the simulated-night end-to-end tests, about six minutes
pytest tests/regression  # the laboratory frames, needs MASTPIPE_TEST_DATA
```

Tests that need real data read the root from `MASTPIPE_TEST_DATA` and the
atlas from `MASTPIPE_ATLAS_DIR` and skip when they are unset; a set variable
pointing at an incomplete folder fails loudly rather than skipping. The
regression test stages the laboratory night in a temporary folder with a
manifest and asserts tolerances set from the notebook's numbers (trace
centres within 0.3 px, per-fiber rms at most 0.05 px, at least 19 lines, band
coverage within 0.03 nm, fiber agreement median under 0.03 px).

`pyproject.toml` sets `filterwarnings = ["error"]`: any warning fails the
suite, in CI as well as locally. Keep it that way; a NaN slice or a
deprecated call surfaces at once instead of in production.

## Writing a new test

- Test a step on hand-made arrays with values you computed yourself first;
  then on the simulator's truth.
- Never read real data in `tests/unit` or `tests/e2e`; real frames belong in
  `tests/regression` behind the environment variable.
- Mark an end-to-end quick-look run `@pytest.mark.slow`; unit tests may use
  the shared simulated-night fixtures unmarked.
- Compare arrays with `numpy.testing.assert_array_equal` (NaN-safe) and every
  scalar field of a dataclass when round-tripping a product.
