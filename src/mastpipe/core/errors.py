"""The errors that stop a run: one exception family, :class:`RunError`, whose classes
carry a ``kind`` (a stable name a caller can branch on) and whether trying again later
may succeed (``retryable``).

Each class also keeps the builtin base its cause raised before the family existed
(``FileNotFoundError``, ``KeyError``, ``TimeoutError``, ``ValueError``, ``OSError``), so
code catching that builtin still catches it. The manifest's content error,
:class:`mastpipe.io.ingest.ManifestError` (kind ``manifest_invalid``), belongs to the
family too; it lives with the manifest reader. Like the rest of :mod:`mastpipe.core`, this
module imports no other layer.
"""

from __future__ import annotations

from typing import ClassVar


class RunError(Exception):
    """An error that stops a run. ``kind`` names the cause (``unexpected`` for the base
    class); ``retryable`` is true when the same run may succeed later unchanged."""

    kind: ClassVar[str] = "unexpected"
    retryable: ClassVar[bool] = False


class InputMissingError(RunError, FileNotFoundError):
    """The night folder is not a directory."""

    kind = "input_missing"


class ManifestMissingError(RunError, FileNotFoundError):
    """The manifest file does not exist."""

    kind = "manifest_missing"


class UnknownInstrumentError(RunError, KeyError):
    """No instrument of that name is shipped. Its message is the plain text (a
    ``KeyError`` would quote it)."""

    kind = "instrument_unknown"

    def __str__(self) -> str:
        return Exception.__str__(self)


class LibraryBusyError(RunError, TimeoutError):
    """The calibration library's lock was not free in time; another process holds it."""

    kind = "library_busy"
    retryable = True


class LibraryInvalidError(RunError, ValueError):
    """The calibration library's index cannot be read: not UTF-8, not JSON, or not a list."""

    kind = "library_invalid"


class OutputUnwritableError(RunError, OSError):
    """The output folder, the run's log or the calibration library cannot be created or
    opened."""

    kind = "output_unwritable"
