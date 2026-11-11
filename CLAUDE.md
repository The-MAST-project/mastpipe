# mastpipe — Claude project notes

**mastpipe** is the data-reduction pipeline for the MAST fiber-fed
spectrographs, **HighSpec** (R ≈ 20 000, narrow-band gratings — three on the
instrument, of which one, `Mg`, is configured here so far; more planned — one
detector) and **DeepSpec** (R ≈ 500, four cameras). It takes raw frames to
wavelength-calibrated, sky-subtracted spectra per fiber, combined spectra per
target, and quality-assurance products. It runs as an unattended quick-look at
the observatory and as a manual reduction tool. Developed by the MAST team; the
scientific core descends from the HighSpec laboratory reduction notebook. The
design is in `docs/design/`.

## Editing this file

`CLAUDE.md` is **read-only by default**. Any change requires explicit user
approval in the conversation before running Edit/Write on it. A deny rule in
`.claude/settings.json` enforces this. Treat the rules here as load-bearing.

## Workflow

- **Plan before substantial work** (a new step, a new recipe, a change to the
  data model or products): explore, propose, list open questions, get explicit
  approval. Then implement.
- **Same model throughout.** Design and implementation are both done with the
  model the user is working with. Do not switch to a smaller model for
  implementation unless the user asks.
- **Test-first.** Every step function gets unit tests on synthetic arrays
  before or with its implementation. Never claim a step works without running
  its tests.
- **Announce model switches**, if any.
- Lightweight tasks (a typo, a single small fix) need no ceremony.
- **Keep `DECISIONS.md` current**: record design decisions (and reversals)
  as part of the change that makes them.
- **Keep the documentation in sync.** Every code change also updates, in the
  same change, the `README.md` and the pages under `docs/` that describe it
  (commands and options, configuration keys, products, behaviour, tutorials)
  and adds a line to `CHANGES.md`. Build the docs with
  `sphinx-build -E -W -b html docs <dir>` (`-E`: the notebook execution cache
  can hide a broken build); a change is not done while it warns.

## Privacy & data safety

- **Observatory data are unpublished research data.** Never send frames,
  spectra, target lists, coordinates, or quality numbers to third-party
  services (pastebins, Gists, image hosts, diagram renderers, external LLM
  tooling) without explicit user approval — *even if the service claims "no
  logging"*.
- **Nothing private in the repository.** No local absolute paths, home
  directories, hostnames, IP addresses, usernames, e-mail addresses, or
  machine-specific configuration in any tracked file, including docs,
  notebooks, tests and CI. Use environment variables (`MASTPIPE_TEST_DATA`
  for real test frames) and relative paths. Strip notebook outputs before
  committing.
- **Test data must be synthetic.** Real frames stay outside the repository and
  are located through `MASTPIPE_TEST_DATA`; tests that need them skip cleanly
  when the variable is unset.
- **Never delete data.** No `rm` on `*.fits`, `*.fits.fz`, `*.csv`, `*.npz`
  data files or anything under a data or products folder. Ask first if a
  destructive action seems necessary.
- **Stay in the repository.** Read and write only inside the repository
  folder; the harness-managed memory and plan directories are the exceptions.

## Secrets & credentials

- **Never commit secrets**: tokens, API keys, `.env`, database URLs with
  passwords, service credentials. `.gitignore` excludes `.env` and `*.local.*`.
  If a secret slips into a tracked file, stop, alert the user, do not push.
- Configuration that differs per machine goes into a git-ignored local file
  or environment variables, documented in a committed example with
  placeholder values only.
- If asked to read or echo a secret in the conversation, confirm first;
  secrets in transcripts are still secrets.

## Network & external services

Allowed without prior approval: read-only documentation fetches (docs sites,
Context7, web search for technical questions); installing packages already
declared in `pyproject.toml`; downloading public reference data (line atlases)
from their canonical sources into `src/mastpipe/data/` when the design calls
for it.

Require explicit user approval: adding a runtime or development dependency;
calling any service that accepts uploads; publishing anything (PyPI,
conda-forge, GitHub Pages); creating or changing GitHub repositories, settings
or collaborators.

## Git & remote safety

- The GitHub repository (`The-MAST-project/mastpipe`) is **public**: every
  tracked file and every commit is world-readable, which is why the privacy
  rules above are strict. Do not change its visibility.
- Commits are fine **when explicitly asked**, in logical chunks; do not bundle
  unrelated changes. Feature branches off `main`.
- **No `git push` without explicit per-push approval.** Never force-push.
  Never `git reset --hard`, `git clean -fd`, `git checkout -- .`, `git
  restore`, or rebase published commits without explicit approval.
- No hook-bypass flags (`--no-verify`, `--no-gpg-sign`). If a hook fails, fix
  the cause. The privacy scan (`python scripts/privacy_scan.py`: paths, e-mails,
  IPs, tokens) must pass; CI runs it on every push.
- Use the machine's existing git identity; do not pass `-c user.*` overrides.

## Stack

Python ≥ 3.11 · numpy · scipy · astropy · pandas · matplotlib · pydantic v2 ·
pyyaml · tzdata (time zones on Windows). Optional extras: `pypeit` (adapter
only), `dev` (pytest, pytest-cov, ruff, nbstripout). ruff is pinned in the
`dev` extra and nowhere else; CI reads the pin from `pyproject.toml`. **No new
top-level dependency without telling the user what and why first.**

## Commands

```bash
pip install -e ".[dev]"          # development install
pytest                           # tiers 1–2, offline, synthetic data only
MASTPIPE_TEST_DATA=<folder> MASTPIPE_ATLAS_DIR=<folder> pytest tests/regression   # tier 3, real lab frames
ruff format --check . && ruff check .
python scripts/privacy_scan.py
mastpipe identify <full-arc.fits> --out <dir> --instrument highspec --configuration Mg --stored-anchors
mastpipe solve <full-arc.fits> --out <dir> --instrument highspec --configuration Mg
mastpipe quicklook <night-folder> --out <dir> --instrument highspec [--layout lab3|survey45] [--manifest frames.csv]
mastpipe quicklook <night-folder> --out <dir> --instrument highspec --configuration Mg --layout survey45   # site form, no manifest
mastpipe show <product.fits>
# quicklook exits 0, or 2 when a frame or group failed; a missing night folder or
# manifest, or a manifest error (`ManifestError`), is one `FAILED:` line on stderr, exit 2.
# `mastpipe reduce` (manual mode, --force <stage>) and `mastpipe focus` are later plans.
```

## Architecture invariants

- **Three layers, strictly separated**: `instruments/` (declarative
  configuration only), `steps/` (pure functions: arrays and dataclasses in,
  dataclasses out; no file access, no global state, no config lookups),
  `pipeline/` + `io/` (ingest, the quick-look runner and its stage modules,
  the calibration library, products, QA, CLI). Do not let a step read a file
  or a config. `core/` (the data model, times, the hash length) is imported
  by every layer and imports none of them.
- **The quick-look recipe is split by stage.** `pipeline/runner.py` chains
  the stages and holds `run_quicklook`; `pipeline/group.py` holds the run
  options, the per-group state and the shared helpers; `masters.py`,
  `traces.py`, `arc_stage.py` and `science.py` hold one stage each. A stage
  module imports `pipeline.group`, never another stage or `pipeline.runner`.
- **No raw header keyword outside `io/ingest.py` and the header maps.**
  Everything downstream uses the pipeline vocabulary in `FrameMeta`.
- **The fiber layout is an input, never inferred.** Trace finding matches
  detections to the layout table and reports discrepancies.
- **A grating or band is configuration, never code.** New HighSpec gratings
  are added as one YAML file under `instruments/highspec/gratings/` (band,
  geometry, atlas, line table); nothing in `steps/` or `pipeline/` may branch
  on a grating name.
- **Every product carries provenance**: pipeline version, configuration hash,
  input checksums, calibration files used, wavelength status
  (`solved` / `transferred` / `reference` / `nominal`).
- **Air wavelengths** internally; vacuum only as an output option.
  Barycentric corrections are stored, never applied to the grid.
- **Cross-platform**: `pathlib` only, no shell-outs, no POSIX-only calls,
  explicit encodings, no assumptions about path separators or line endings.
- **Quick-look never stops on one bad frame; manual mode raises.** Both run
  the same step functions.
- **PypeIt only through a planned optional adapter** (not written yet), lazily
  imported, behind the native step signatures, version-pinned, with its own
  tests; the `pypeit` extra pre-installs its dependency.
- **Vocabulary comes from MAST_common.** Instrument, grating, frame-type, lamp
  and unit names are those of `common.spec`, mirrored in `mastpipe.vocab` and
  checked by the `fleet` CI job; never introduce a name of our own for
  something MAST_common already names.
- **Notebooks contain no logic**; they import the package.

## Testing conventions

- Unit tests build inputs with `mastpipe.sim` or small hand-made arrays,
  fixed seeds, no network, no real data.
- End-to-end tests run the quick-look recipe on a simulated night and compare
  every product with the simulator's truth.
- Regression tests run the laboratory frames and bound the results by the
  laboratory notebook's (trace centres, per-fiber wavelength rms ≤ 0.05 px,
  per-fiber line counts); they skip without real data. Resolving power and
  best focus are not tested: there is no width or focus recipe yet.
- CI runs on Linux, Windows and macOS; a change is not done until the matrix
  is green.

## Fleet conventions

Follow `The-MAST-project/mast-claude-config` (its
`code-validation-and-ci-guidelines.md`), with two deliberate departures: CI
tests Linux, Windows and macOS, and Python 3.11 stays supported, because
mastpipe is a pip-installable library, not a fleet service.

## Out of scope (first release)

Absolute flux calibration, telluric correction, pixel-to-pixel flat-fielding
(present but off), radial-velocity measurement, integration with the
observatory supervisor and dashboard.
If a request lands here, say so rather than adding scope silently.

## Memory

Repo-local memory at `.claude/memory/` with a `MEMORY.md` index, **git-ignored**.
One fact per file with YAML frontmatter (`name`, `description`,
`metadata.type` = user | feedback | project | reference). Read it at the start
of work; do not record what the repository or git history already captures.
