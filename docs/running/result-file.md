# The result file

A program that starts `quicklook` (the observatory's control software, a cron job, a
batch script) needs to know what happened without reading the log or the printed text.
It gets two things: the **exit status**, and, when it asks for one, a **result file**:

```bash
mastpipe quicklook night/ --out out --instrument highspec --configuration Mg \
    --summary-json out/result.json --quiet
```

The printed output (stdout) stays text for people and is not part of this contract.

## When the file is written

`--summary-json PATH` writes one JSON document to `PATH` as the command's last act, for
every outcome except a usage error (exit status 2, when the command line itself is
wrong) and an interrupt (Ctrl-C). Without the option nothing is written. The folder of
`PATH` is created when missing, and the file is written atomically: a reader sees the
whole document or no file, never half of one.

A run that is interrupted (Ctrl-C) or killed writes no result file, and its exit status
is the interpreter's or the operating system's own (130 or a signal's status, say), not
one of the table below.

If the printed summary cannot be written to stdout (an encoding that cannot spell a file
name, a full device), one line on stderr says `could not print the summary: …`; the exit
status is still the run's and the result file is still written. (On a full device the
interpreter may fail once more to flush stdout as it exits and end with status 120; the
result file then still holds the run's own `exit_code`.)

If the file cannot be written, the exit status is still the run's, and one line on
stderr says `could not write the result file PATH: …`. A caller that asked for a result
and finds none after status 0, 1, 3 or 4 treats the run as `error`. A caller that reuses
one `PATH` across runs **removes the old file first**, since a stale one would otherwise
be read as this run's.

## Exit status

| Status | `status` | Meaning | Result file |
|---|---|---|---|
| 0 | `complete` | The run finished and nothing failed. | yes |
| 0 | `empty` | The folder holds no frame of this instrument and nothing failed. Nothing to do. | yes |
| 1 | `error` | An unexpected error: a defect in mastpipe, or `--on-error raise` re-raising a frame's or group's error (a run error raised inside a group included). Traceback on stderr and in the log. | yes |
| 2 | (usage) | The command line is wrong: the argument parser's own status, a malformed `--set` included. | no |
| 3 | `partial` | The run finished and at least one frame or group failed. The products of the others are good. | yes |
| 4 | `failed` | The run could not be done, for a reason mastpipe recognizes (`run_error.kind` below), found before the calibration library was open. Nothing was reduced. | yes |

A folder with no usable frame in which something did fail (its only file is not
readable FITS) is `partial`, not `empty`, so a folder of corrupt files is never passed
over as empty.

## Fields

The document of schema version 1 always has every field below; a value the run never
reached is `null` (or an empty list).

`schema`
: The schema version, `1`.

`command`
: `"quicklook"`.

`mastpipe_version`
: The version of mastpipe that ran.

`status`, `exit_code`
: The name and number of the exit-status table above; `exit_code` repeats the process's
  status, so the file alone is enough.

`started`, `finished`
: The run's start and end, UTC, `YYYY-MM-DDTHH:MM:SSZ`; set on every run, one stopped
  at once (a night folder that does not exist) included.

`input`, `out`, `library`
: The night folder, the `--out` folder and the calibration library the run used, as
  given (neither resolved nor made absolute) but with forward slashes, as is
  `options.manifest`. `library` is `null` when the run stopped before opening it.

`options`
: The options of the run: `instrument`, `configuration`, `camera`, `layout`,
  `time_policy`, `force`, `on_error`, `manifest` (the path as given, forward slashes), `reference` (the
  `--reference` id) and `set` (the `--set` overrides as given, a list of strings).
  An option not given is `null` (`force` is `false`, `on_error` `"continue"`, `set`
  `[]`).

`log`
: The run's own log file, relative to `out`, with forward slashes
  (`logs/quicklook_<UTC stamp>_<6 hex>.log`); `null` when the run stopped before the log
  existed. Read it from here; never guess the name.

`frames_scanned`
: How many frames the scan found (after the manifest); `null` when the run stopped
  before the scan.

`counts`
: The lengths of the five lists below: `products`, `qc`, `skipped`, `reduced_again`,
  `failed`. They equal the counts in the printed summary.

`products`, `qc`
: Every file the run wrote, as `{"root": ..., "path": ...}`: `root` is `"library"` for a
  file in the calibration library (which may lie inside or outside `--out`), else
  `"out"` for one under the `--out` folder (when both hold a file, the deeper folder is
  its root; which folder holds a file is decided on absolute paths, so `--out .` does not
  hold `../lib`); `path` is relative to that root, with forward slashes.

`skipped`
: What the run left out on purpose: `{"scope", "source", "reason"}`.

`reduced_again`
: The sources (relative to the night folder) of science products reduced again because
  their key changed.

`failed`
: What the run could not do: `{"scope", "source", "kind", "error_type", "message"}`.
  `error_type` names the exception's class, or is `null` when no exception caused it.

`run_error`
: `null`, or, for a `failed` or `error` run, `{"kind", "error_type", "message",
  "retryable"}`: what stopped the run.

On a `failed` or `error` run, the lists hold what the run had gathered when it stopped
(products written before an `--on-error raise` abort included).

### `scope` and `source`

Every `skipped` and `failed` item has both; a `source` is never a sentence.

| `scope` | `source` |
|---|---|
| `frame` | the file's path relative to the night folder, with forward slashes (a target that could not be combined is an item of its exposure, the target named in the text) |
| `group` | `<camera>/<configuration>` |
| `manifest` | the manifest row's own `source` cell; the text names the row |
| `run` | `null` (reserved: nothing uses it yet) |

### `failed[].kind`

`error` unless the cause is one a caller may act on:

| `kind` | Cause |
|---|---|
| `error` | anything else |
| `configuration_unknown` | a group whose frames name no grating, and no `--configuration` was passed |
| a run-error kind | a run error raised inside a group or for a frame, under `--on-error continue`: e.g. `library_busy` (the library's lock was not free in time while a group added to it) or `library_invalid` (the index could not be read while a group wrote to it) |

Any kind of the `run_error.kind` table below may appear on an item, and new kinds may be
added (see Stability).

### `run_error.kind`

A run error stops the run before its first group: before the calibration library was
open, so nothing has been reduced. A run error raised later (under `--on-error raise`,
re-raised from a group) is an `error` run, kind `unexpected`, with its class in
`error_type`.

| `kind` | Cause | `retryable` |
|---|---|---|
| `input_missing` | the night folder is not a directory | false |
| `manifest_missing` | `--manifest` names no file | false |
| `manifest_invalid` | the manifest has an error (the file, a column or a row) | false |
| `instrument_unknown` | no instrument of that name is shipped | false |
| `library_busy` | another process held the calibration library's lock past the timeout when the library was opened | **true** |
| `library_invalid` | the library's `index.json` cannot be read as an index | false |
| `output_unwritable` | `--out`, its log folder or log file, or the library folder cannot be created or locked | false |
| `unexpected` | an `error` run (status 1): any other exception, or a run error raised after the library was open | false |

`retryable` is mastpipe's statement that the same command may succeed later with nothing
changed by a person.

## Stability

Adding a field, or a new `kind`, does not change `schema`. Removing or renaming a field,
or changing what a value means, does. Readers ignore fields they do not know. Both kind
lists grow only by adding names, never by renaming.

## Example

A `partial` run (the committed example a test compares a run's document with; the
placeholders stand for the paths and version of that run):

```{literalinclude} ../../tests/unit/data/summary_schema_1.json
:language: json
```
