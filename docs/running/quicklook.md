# Quick-look: reducing a night

```bash
mastpipe quicklook night/ --out out --instrument highspec --configuration Mg
```

One command reduces everything in the folder it can, records what it could
not, and exits. It is meant to run unattended at the observatory and again on
a laptop; both runs give the same products from the same inputs.

## What happens

Per camera and grating found in the folder:

1. **Masters.** Bias, dark and background frames are stacked (dark and
   background per exposure time) and flats are stacked after bias and dark subtraction;
   each master goes into the library with its validity keys. A kind the night
   lacks is taken from the library.
2. **Traces.** Found on the flat master (or, without a flat, on a full arc),
   matched to the fiber layout and re-keyed by the layout's trace ids;
   unmatched detections are dropped and listed. Relative throughput per fiber
   from the flat.
3. **Arcs.** Each arc frame is classified from which traces carry light. A
   *full* arc is re-solved from the library's reference and the result added
   as a new reference (or skipped with "no reference …; run identify" when
   there is none).
   A *calibration-fiber* arc is fitted for drift against the reference valid at
   its time.
4. **Science frames.** Calibrated, cosmic rays flagged, extracted; wavelengths
   from the reference valid at the frame's mid-time, corrected by the drift of
   the paired arcs; sky subtracted from the telescope's sky fiber; passed
   through the flux and telluric stages; written with provenance and a QC file.
5. **Combination.** The science fibers of each target in the exposure are
   combined on a common wavelength grid, version by version, into one file per
   target.

## Options that matter

`--time-policy nearest_before | nearest | interpolate`
: How a science frame picks its drift arcs. `nearest_before` (the quick-look
  default) uses the latest arc before the frame, falling back to the nearest
  after; `nearest` the closest in time either way; `interpolate` the arcs
  before and after, weighted by time. Only arcs of the same grating, fitted
  against the frame's reference, and within `drift.max_age_hours` (default 12)
  qualify. With none, the frame gets status `reference`.

`--configuration Mg`
: The grating of frames whose header does not name it (no `GRATING` card and
  no grating alias in the file name). HighSpec ships one configuration today,
  `Mg`; the fleet also names `Halpha` and `Ca`, which become valid here once their
  grating files are added under `instruments/highspec/gratings/`. The control
  software writes no such card today, so a night from the site needs this
  option; without it those frames are not reduced and their group,
  `group main/(no configuration)`, fails with the list of the instrument's
  configurations. A frame that names its grating
  (in its header, its file name or a manifest row) keeps it, even a grating
  other than this one; one that does not takes this one, and its products,
  masters and library entries are named by it. Such a frame's product records
  where its grating came from: the provenance step `configuration Mg from the
  run's --configuration (none in the header, file name or manifest)`, and
  `configuration_origin` in its QC file (`run`, or `frame` for a grating the
  frame itself named; in the QC of every science frame, full arc and drift
  arc).

`--manifest file.csv`
: Metadata overrides by relative path, glob pattern or file name; see
  [The night folder](night-folder.md).

`--layout lab3 | survey45`
: The fiber preset (else a manifest's `fiber_preset`, else the instrument's
  default); selects the preset's parameters too.

`--reference ID`
: Pin one reference solution for every frame instead of the one valid at each
  frame's time.

`--on-error continue | raise`
: With `continue` (default) a frame that fails is recorded and the run goes
  on; a failure in the masters or traces fails that camera-and-grating group.
  With `raise` the first error stops the run with its traceback (exit status 1).

`--force`
: Reduce products that already exist, and rebuild masters the library holds.

`--set group.key=value`
: Override any step parameter for this run, e.g. `--set drift.max_age_hours=24`
  or `--set extract.halfwidth=4`. The configuration actually used is hashed
  into every product, so an override changes the idempotence key. An override
  without `=` or with an empty key segment is a usage error (exit status 2).

`--summary-json PATH`
: Write the run's result to `PATH` as one JSON document, for a program that
  starts `quicklook` and needs to know what happened; see
  [The result file](result-file.md).

`--quiet`
: Keep the per-frame progress lines off the console: warnings and errors still go
  to stderr, the printed summary to stdout, and the log file still gets every line.

## Idempotence and re-reduction

A science product on disk is skipped when its **key** matches the run's: the
input frame's checksum, the configuration hash (overrides and `--time-policy`
included), the reference solution's hash, the sorted checksums of the drift
arcs it would use, the layout preset and the hash of its CSV file, a hash of
the masters used and the frame's time zone. A product whose key changed is
reduced again. Masters the calibration library already holds, stacked from
the same frames with the same `calib` parameters and detector (and, for a
flat, calibrated with masters of the same content), are reused rather than
rebuilt; `--force` rebuilds them. The case this is built for: the
first science frames of a night are reduced before any calibration-fiber arc
exists and get status `reference`; when the night's first arc lands and the
quick-look runs again, those frames are reduced once more and become
`transferred`. The run summary lists them under `reduced_again`.

## What the run reports

The console shows one line per product and a summary; the log in `out/logs/`
has everything with timestamps, one new file per run
(`quicklook_<UTC stamp>_<6 hex>.log`, the random suffix keeping two runs started in
the same second apart). When something stops the whole run, the runner writes why to
the log file only (`run failed: ...`, the traceback after it for an unexpected error,
or `run interrupted` for Ctrl-C), and the command prints its own one line on stderr: the
`FAILED: ...` line, or the traceback of an unexpected error. `--quiet` (in Python,
`RunOptions(quiet=True)`) keeps the console to warnings and errors; the log still has
every line, and the run summary's `log` names the file.

The exit status says which kind of thing happened:

| Status | Name | Meaning |
|---|---|---|
| 0 | `complete` | The run finished and nothing failed. |
| 0 | `empty` | The folder holds no frame of this instrument and nothing failed: nothing to do. |
| 1 | `error` | An unexpected error: a defect in mastpipe, or `--on-error raise` re-raising a frame's error. The traceback is on stderr and in the log. |
| 2 | (usage) | The command line is wrong (the argument parser's own status), a malformed `--set` included. |
| 3 | `partial` | The run finished and at least one frame or group failed; the other products are good. |
| 4 | `failed` | The run could not be done, for a reason mastpipe recognizes; nothing was reduced. |

A folder whose only file cannot be read as FITS has no frame but a failure, so it is
`partial` (3), never `empty`. Status 4 is one `FAILED:` line on stderr naming the
folder, the file or the row (see [The manifest](night-folder.md#the-manifest)); its
causes and their `kind` names are listed in [The result file](result-file.md).
Skipped frames and failures are listed with their reasons; see
[Troubleshooting](troubleshooting.md).

## Running on a laptop afterwards

The same command on the same folder and output gives the same products. To
redo a night with a different policy:

```bash
mastpipe quicklook night/ --out out --instrument highspec --configuration Mg --time-policy interpolate --force
```

Products from earlier runs are overwritten only when `--force` is given or
their key changed; the library keeps every reference and master it has seen.
