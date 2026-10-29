# mastpipe — Claude project notes

**mastpipe** is the data-reduction pipeline for the MAST fiber-fed
spectrographs, **HighSpec** (R ≈ 20 000, narrow-band gratings — three today,
more planned — one detector) and **DeepSpec** (R ≈ 500, four cameras). It takes raw frames to
wavelength-calibrated, sky-subtracted spectra per fiber, combined spectra per
target, and quality-assurance products. It runs as an unattended quick-look at
the observatory and as a manual reduction tool. Developed by the MAST team;
its scientific core descends from the team's
HighSpec reduction notebook. The design is in `docs/design/`.

## Editing this file

`CLAUDE.md` is **read-only by default**. Any change requires explicit user
approval in the conversation before running Edit/Write on it. A deny rule in
`.claude/settings.json` enforces this. Treat the rules here as load-bearing.

## Workflow

- **Plan before substantial work** (a new step, a new recipe, a change to the
  data model or products): explore, propose, list open questions, get explicit
  approval. Then implement.
- **Same model throughout.** Design and implementation are both done with the
  model the user is working with (Claude Fable at the time of writing). Do not
  switch to a smaller model for implementation unless the user asks.
- **Test-first.** Every step function gets unit tests on synthetic arrays
  before or with its implementation. Never claim a step works without running
  its tests.
- **Announce model switches**, if any.
- Lightweight tasks (a typo, a single small fix) need no ceremony.
- **Keep `DECISIONS.md` current**: record design decisions (and reversals)
  as part of the change that makes them.

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

- The GitHub repository is **private**. Do not change its visibility.
- Commits are fine **when explicitly asked**, in logical chunks; do not bundle
  unrelated changes. Feature branches off `main`.
- **No `git push` without explicit per-push approval.** Never force-push.
  Never `git reset --hard`, `git clean -fd`, `git checkout -- .`, `git
  restore`, or rebase published commits without explicit approval.
- No hook-bypass flags (`--no-verify`, `--no-gpg-sign`). If a hook fails, fix
  the cause. The pre-commit privacy scan (paths, e-mails, IPs, tokens) must
  pass.
- Use the machine's existing git identity; do not pass `-c user.*` overrides.

## Stack

Python ≥ 3.11 · numpy · scipy · astropy · pandas · matplotlib · pydantic v2 ·
pyyaml. Optional extras: `pypeit` (adapter only), `dev` (pytest, pytest-cov,
ruff, nbstripout). **No new top-level dependency without telling the user
what and why first.**

## Commands

```bash
pip install -e ".[dev]"          # development install
pytest                           # tiers 1–2, offline, synthetic data only
MASTPIPE_TEST_DATA=<folder> pytest tests/regression   # tier 3, real lab frames
ruff check src tests && ruff format --check src tests
mastpipe quicklook <night-folder>
mastpipe reduce --config run.yaml --manifest frames.csv
mastpipe focus <scan-folder>
mastpipe show <product.fits>
```

## Architecture invariants

- **Three layers, strictly separated**: `instruments/` (declarative
  configuration only), `steps/` (pure functions: arrays and dataclasses in,
  dataclasses out; no file access, no global state, no config lookups),
  `pipeline/` + `io/` (ingest, runner, products, QA, CLI). Do not let a step
  read a file or a config.
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
- **PypeIt only through `ext/pypeit_adapter.py`**, lazily imported, behind the
  native step signatures, version-pinned, with its own tests.
- **Notebooks contain no logic**; they import the package.

## Testing conventions

- Unit tests build inputs with `mastpipe.sim` or small hand-made arrays,
  fixed seeds, no network, no real data.
- End-to-end tests run the quick-look recipe on a simulated night and compare
  every product with the simulator's truth.
- Regression tests reproduce the published laboratory numbers (per-fiber
  wavelength rms, resolving power, best focus) and skip without real data.
- CI runs on Linux, Windows and macOS; a change is not done until the matrix
  is green.

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
