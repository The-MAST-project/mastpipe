"""The result of a ``quicklook`` run for a calling program: its outcome and exit code
(:func:`outcome`), and the schema-1 JSON document (:func:`document`, :func:`write`) that
``--summary-json`` writes.

Exit codes: 0 ``complete`` or ``empty``, 1 ``error`` (an unexpected exception), 2 a usage
error (the argument parser's own, never returned here), 3 ``partial`` (some frame or group
failed), 4 ``failed`` (a :class:`~mastpipe.core.errors.RunError` stopped the run before
the calibration library was open, so nothing was reduced; see :func:`stopped_by_run_error`).
The fields and the stability rule are documented on the user-docs page "The result file".
"""

from __future__ import annotations

import json
import os
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from mastpipe.__version__ import __version__
from mastpipe.core.errors import RunError
from mastpipe.io.atomic import atomic_write
from mastpipe.pipeline.group import RunOptions, RunSummary

SCHEMA = 1
"""The result document's schema version: it changes when a field is removed or renamed or
a value changes meaning, never when a field or a ``kind`` is added."""

EXIT_OK, EXIT_ERROR, EXIT_PARTIAL, EXIT_FAILED = 0, 1, 3, 4


def stopped_by_run_error(summary: RunSummary, error: BaseException | None) -> bool:
    """Whether ``error`` is a run error of the result (``failed``, exit 4, with its own
    ``kind`` and ``retryable``): a :class:`~mastpipe.core.errors.RunError` raised before
    the calibration library was open (``summary.library`` still ``None``), when nothing
    has been reduced. Once the library is open, any exception leaving the run -- a
    ``RunError`` re-raised from a group under ``--on-error raise`` included -- is an
    unexpected ``error``."""
    return isinstance(error, RunError) and summary.library is None


def outcome(summary: RunSummary, error: BaseException | None) -> tuple[str, int]:
    """The run's status and exit code: ``("failed", 4)`` when
    :func:`stopped_by_run_error`; ``("error", 1)`` for any other exception;
    ``("partial", 3)`` when something failed; ``("empty", 0)`` when the scan found no frame;
    else ``("complete", 0)``."""
    if stopped_by_run_error(summary, error):
        return "failed", EXIT_FAILED
    if error is not None:
        return "error", EXIT_ERROR
    if summary.failed:
        return "partial", EXIT_PARTIAL
    if summary.frames_scanned == 0:
        return "empty", EXIT_OK
    return "complete", EXIT_OK


def _absolute(path: Path) -> Path:
    """``path`` made absolute and normalised (``..`` removed) without following links."""
    return Path(os.path.abspath(path))


def _written(path: Path, out_dir: Path, library: Path | None) -> dict[str, str]:
    """A written path relative to its root, with forward slashes: of ``out`` and
    ``library``, the deeper root that holds it, decided on absolute paths (so ``--out .``
    does not hold ``../lib``, and a library inside ``--out`` is the master's root)."""
    target = _absolute(path)
    roots = [("out", _absolute(out_dir))]
    if library is not None:
        roots.append(("library", _absolute(library)))
    holding = [(name, root) for name, root in roots if target.is_relative_to(root)]
    if not holding:
        raise ValueError(f"{path} lies outside {out_dir} and the library {library}")
    name, root = max(holding, key=lambda item: len(item[1].parts))
    return {"root": name, "path": target.relative_to(root).as_posix()}


def _posix(path: Path | None) -> str | None:
    return None if path is None else path.as_posix()


def _run_error(summary: RunSummary, error: BaseException | None) -> dict[str, Any] | None:
    if error is None:
        return None
    if stopped_by_run_error(summary, error):
        assert isinstance(error, RunError)
        kind, retryable = error.kind, error.retryable
    else:
        kind, retryable = "unexpected", False
    return {"kind": kind, "error_type": type(error).__name__, "message": str(error), "retryable": retryable}


def document(
    summary: RunSummary,
    options: RunOptions,
    night_dir: Path,
    out_dir: Path,
    overrides: Sequence[str],
    error: BaseException | None,
) -> dict[str, Any]:
    """The schema-1 result document of a run: ``summary`` as the run left it (its facts
    ``None`` and its lists empty when the run stopped before them), the run's options with
    ``overrides`` the ``--set`` texts as given, and ``error`` the exception that ended the
    run, if any."""
    status, code = outcome(summary, error)
    library = summary.library
    return {
        "schema": SCHEMA,
        "command": "quicklook",
        "mastpipe_version": __version__,
        "status": status,
        "exit_code": code,
        "started": summary.started,
        "finished": summary.finished,
        "input": night_dir.as_posix(),
        "out": out_dir.as_posix(),
        "library": _posix(library),
        "options": {
            "instrument": options.instrument,
            "configuration": options.configuration,
            "camera": options.camera,
            "layout": options.layout,
            "time_policy": options.time_policy,
            "force": options.force,
            "on_error": options.on_error,
            "manifest": _posix(options.manifest),
            "reference": options.reference_id,
            "set": list(overrides),
        },
        "log": None if summary.log is None else _absolute(summary.log).relative_to(_absolute(out_dir)).as_posix(),
        "frames_scanned": summary.frames_scanned,
        "counts": {
            "products": len(summary.products),
            "qc": len(summary.qc),
            "skipped": len(summary.skipped),
            "reduced_again": len(summary.reduced_again),
            "failed": len(summary.failed),
        },
        "products": [_written(path, out_dir, library) for path in summary.products],
        "qc": [_written(path, out_dir, library) for path in summary.qc],
        "skipped": [{"scope": item.scope, "source": item.source, "reason": item.reason} for item in summary.skipped],
        "reduced_again": list(summary.reduced_again),
        "failed": [
            {
                "scope": item.scope,
                "source": item.source,
                "kind": item.kind,
                "error_type": item.error_type,
                "message": item.message,
            }
            for item in summary.failed
        ],
        "run_error": _run_error(summary, error),
    }


def write(path: Path, doc: Mapping[str, Any]) -> None:
    """Write ``doc`` to ``path`` as JSON (indent 1, sorted keys, ASCII, a final newline),
    atomically: a reader sees the whole file or none. The folder is created when missing."""
    text = json.dumps(doc, indent=1, sort_keys=True, ensure_ascii=True) + "\n"
    with atomic_write(path) as tmp, tmp.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)
