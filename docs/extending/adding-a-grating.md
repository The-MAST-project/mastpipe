# Adding a HighSpec grating

HighSpec has three gratings today (Mg b, Hα, Ca II H&K) and more are planned.
A grating is a configuration entry, never code: one YAML block, one line
table, and, once the first full arc has been taken, one anchors file. The
registry discovers the block by file name.

## 1. The grating block

`src/mastpipe/instruments/highspec/gratings/<name>.yaml`, where `<name>` is
the configuration name the header map's `GRATING` keyword (or the file-name
pattern) resolves to:

```yaml
name: Halpha
band_nm: [649.0, 664.0]        # the detector's wavelength range, generous
resolving_power: 20000
geometry:
  model: grating               # the wavelength model kind
  groove_spacing_nm: 367.0     # 2724.8 lines/mm
  incidence_angle_deg: 45.0    # held fixed in the fit
  order: 1
  camera_focal_mm: 479.41
  central_wavelength_nm: 656.3
  dispersion_sign: 1           # -1 when wavelength falls with pixel
line_table: lines/Halpha.csv
params: {}                     # per-grating overrides of any step parameter
```

The geometry constants are the optical design's, as in the reduction notebook;
they are never fitted. The `name` must equal the file stem and be one of the
fleet's grating names, `mastpipe.vocab.GRATING_NAMES` (`Ca`, `Halpha`, `Mg`,
`Future`, spelled exactly so); the registry rejects any other.

## 2. The line table

The catalogue lines of the band with a flag saying which survive at the
instrument's resolution. It is generated once from the Th-Ar atlas of Murphy
et al. and committed:

```bash
MASTPIPE_ATLAS_DIR=<atlas folder> python scripts/make_line_table.py \
    --instrument highspec --configuration Halpha --resolving-power 20000
```

The script convolves the atlas to the resolving power and keeps a line when
the centroid of its convolved feature stays within 0.02 FWHM of the catalogue
wavelength with a signal-to-noise above 20, the notebook's criteria. The
output, `lines/<name>.csv`, has columns `wavelength_nm, species, logI, usable`
and a header comment citing the atlas. The simulator uses the same table to
place lines, so a simulated night of the new grating exists as soon as the
table does.

## 3. Add the header mapping

If the control software names the grating differently from `<name>`, add the
spelling to the `configuration` field's `map` in `header_map.yaml`, and to its
`filename_regex` if laboratory files carry the name instead.

## 4. The first reference

Take a full arc, run `identify` without anchors, read two pairs off the
panels, and run it again with `--anchors` (see
[The first reference solution](../running/first-reference.md)). To make those
pairs the stored default for the grating and a layout, write
`anchors/<name>_<preset>.yaml`:

```yaml
trace: 0                              # layout trace id the pairs were measured on
pairs:
  - [656.2852, 812.40]                # [wavelength_nm, pixel]
  - [660.1300, 1241.70]
```

From then on `identify --stored-anchors` works and `quicklook` solves later
full arcs on its own.

## 5. Check

```bash
pytest tests/unit/test_registry.py       # the block loads and validates
mastpipe quicklook sim_night/ --out out --instrument highspec --configuration Halpha
```

`tests/unit/test_registry.py::test_extra_grating_file_is_picked_up` is the
contract: a new block appears in `load_instrument("highspec").configurations`
with no code change.
