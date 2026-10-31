# Conventions

The repository's `CLAUDE.md` is the authoritative statement of the working
rules; this page restates the ones a contributor meets first.

## Privacy and data safety

- **Nothing private in the repository.** No local absolute paths, home
  directories, hostnames, IP addresses, usernames, e-mail addresses or
  credentials in any tracked file — code, tests, docs, CI. Real data are
  reached only through `MASTPIPE_TEST_DATA` and `MASTPIPE_ATLAS_DIR`.
  `scripts/privacy_scan.py` enforces it and CI runs it.
- **Test data are synthetic.** Real frames stay outside the repository; tests
  that need them skip cleanly without the variable.
- **Observatory data are unpublished research data.** Frames, spectra, target
  lists and quality numbers never go to third-party services.
- **Never delete data.** No removal of FITS, CSV or product files by code or
  by hand without asking.

## Code

- Python ≥ 3.11; hard dependencies only numpy, scipy, astropy, pandas,
  matplotlib, pydantic v2 and PyYAML. A new top-level dependency is a decision
  to discuss first.
- Cross-platform: `pathlib` only, `encoding="utf-8"` on every text open,
  `newline=""` on CSV writes, no shell-outs, no POSIX-only calls. CI runs on
  Linux, Windows and macOS.
- Steps are pure (see [Architecture](architecture.md)); raw FITS keyword
  names appear only in `io/ingest.py`, `io/headers.py` and the header maps.
- `ruff` for lint and format, line length 100.
- Docstrings in the third person; no first-person voice in code or docs.

## Testing

- Test first: a failing test, then the code, then the test green.
- `filterwarnings = error`: a warning is a failure.
- Mark end-to-end quick-look runs and the laboratory regression `slow`; keep
  the unit tier to a few minutes.

## Products and provenance

- Every product carries the pipeline version, the configuration hash, input
  checksums, calibrations used, reference identifier and hash, drift arcs,
  layout preset and statuses.
- Raw extractions are never replaced; later stages add versions.
- Wavelengths are air wavelengths in nanometres; barycentric corrections are
  stored, never applied.
- The fiber layout is an input table; fibers are keyed by its trace ids.

## Records

- `DECISIONS.md` records design decisions and reversals as part of the change
  that makes them, dated.
- `CHANGES.md` has one entry per change under the version heading.
- Commit messages describe the change in the imperative (`feat(steps): …`,
  `fix(pipeline): …`, `docs: …`); commits are small and logically separate.

## Git

- Local commits on `main`; the GitHub repository is created and pushed only
  with the project owner's explicit approval, as a private repository.
- No force pushes, no history rewriting of published commits, no hook-bypass
  flags.
