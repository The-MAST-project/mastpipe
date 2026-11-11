# Architecture: three layers and one rule

```text
src/mastpipe/
  instruments/   declarative configuration: YAML, CSV — no code beyond validation
  steps/         pure functions: arrays and dataclasses in, dataclasses out
  io/            reading frames and writing products
  pipeline/      the runner, the calibration library, time pairing, provenance, commands
  config/        the Pydantic models and the layered loader
  core/          the dataclasses the layers exchange
  sim/           the synthetic-night simulator
  cli.py         argparse entry point; commands register themselves
```

## Instruments: configuration, not code

An instrument is a folder of data files read by `instruments.registry`:

- `config.yaml` — cameras (shape, crop, orientation, gain, read noise,
  saturation), the step parameters, the layout presets and their per-preset
  parameter blocks, how arcs are classified;
- `header_map.yaml` — the pipeline's vocabulary mapped to raw keywords, with
  types, value maps, missing values, time references, a file-name fallback,
  and the classification rules;
- `layouts/*.csv` — one fiber layout per preset;
- `gratings/*.yaml` (HighSpec) or `bands/*.yaml` (DeepSpec) — one block per
  configuration: band, resolving power, geometry, line table;
- `lines/*.csv` — the catalogue lines of a band with a usable flag;
- `anchors/*.yaml` — stored anchor pairs for `identify`.

Adding a grating is one YAML file; the registry discovers it. No code in this
layer runs anything.

## Steps: pure functions

The rule the whole design rests on: a function in `steps/` **never reads a
file, never logs to a file, never looks at a configuration object**. It takes
arrays, dataclasses and plain parameters (the `*Params` models from
`config.models`, which are frozen data) and returns dataclasses. Consequences:

- every step is unit-testable on small synthetic arrays and on the simulator's
  truth;
- the same step serves the quick-look, the manual reduction and `identify`;
- a step could be wrapped by another framework (the PypeIt adapter of a later
  release) without touching it.

The one documented exception is the atlas file reader in `steps.atlas`, used
offline by a script and by `identify`'s panels, never by the runner.

## I/O and pipeline: everything that touches disk

`io.ingest` is the only place, with the header maps, where a raw FITS keyword
name may appear. It turns files into `FrameRef`s whose `FrameMeta` speaks the
pipeline's vocabulary; nothing downstream sees a raw header. `io.products`
writes and reads every product as FITS with provenance cards and gives back
the same dataclasses the steps produce. `io.qa` draws figures.

`pipeline.runner` chains the stages per recipe, decides what to reduce and what
to skip, and writes. The quick-look's stages live in their own modules —
`pipeline.masters`, `pipeline.traces`, `pipeline.arc_stage` and
`pipeline.science` — over the per-group state (masters, traces, drift models)
and shared helpers of `pipeline.group`; a stage module imports `pipeline.group`,
never another stage or the runner. `pipeline.callib` is the calibration library
on disk. `pipeline.arcs` pairs science frames with drift arcs by time.
`pipeline.identify` is the attended bootstrap of a reference.
`pipeline.provenance` builds the header cards. `pipeline.commands` registers the
CLI sub-commands.

## Configuration layering

`config.loader.build_effective` merges, in order, the instrument's step
parameters, the configuration block's parameters, the layout preset's
parameters, and the command-line `--set` overrides, into an `EffectiveConfig`
validated by Pydantic models that forbid unknown keys. Its hash (over the
values only) goes into every product and into the idempotence key.

## Where to add things

| you want to … | touch |
| --- | --- |
| add a grating | `instruments/highspec/gratings/`, `lines/`, `anchors/` — no code |
| add an instrument | a new folder under `instruments/` — see [Adding an instrument](adding-an-instrument.md) |
| add a step parameter | the `*Params` model in `config/models.py` and the step that uses it |
| add a reduction step | a pure function in `steps/`, called from the runner, with its product columns in `io/products.py` |
| add a product | a writer/reader pair in `io/products.py`, a name in `product_name`, a QC report |
| add a command | a registrar in `pipeline/commands.py` appended to `cli.COMMANDS` |
