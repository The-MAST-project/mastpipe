"""Tests for the run errors: one exception family with a ``kind`` and ``retryable``, each
class keeping its builtin base; and for :func:`~mastpipe.pipeline.group.fail`, which
records a failure with the kind of a run error."""

import pytest

from mastpipe.core.errors import (
    InputMissingError,
    LibraryBusyError,
    LibraryInvalidError,
    ManifestMissingError,
    OutputUnwritableError,
    RunError,
    UnknownInstrumentError,
)
from mastpipe.instruments.registry import load_instrument, load_layout
from mastpipe.io.ingest import ManifestError
from mastpipe.pipeline.group import RunSummary, fail


@pytest.mark.parametrize(
    ("cls", "kind", "retryable", "base"),
    [
        (InputMissingError, "input_missing", False, FileNotFoundError),
        (ManifestMissingError, "manifest_missing", False, FileNotFoundError),
        (ManifestError, "manifest_invalid", False, ValueError),
        (UnknownInstrumentError, "instrument_unknown", False, KeyError),
        (LibraryBusyError, "library_busy", True, TimeoutError),
        (LibraryInvalidError, "library_invalid", False, ValueError),
        (OutputUnwritableError, "output_unwritable", False, OSError),
    ],
)
def test_run_error_family(cls, kind, retryable, base):
    error = cls("the message")
    assert isinstance(error, RunError) and isinstance(error, base)
    assert (error.kind, error.retryable) == (kind, retryable)
    assert str(error) == "the message"


def test_run_error_defaults():
    error = RunError("the message")
    assert (error.kind, error.retryable) == ("unexpected", False)


def test_unknown_instrument():
    with pytest.raises(UnknownInstrumentError, match="unknown instrument 'nosuch'"):
        load_instrument("nosuch")


def test_unknown_layout_is_not_a_run_error():
    with pytest.raises(KeyError) as excinfo:
        load_layout("highspec", "nosuch")
    assert not isinstance(excinfo.value, RunError)


def test_fail_records_strings_exceptions_and_kinds():
    summary = RunSummary()
    fail(summary, "frame", "a.fits", "a plain message", prefix="ignored: ")
    fail(summary, "frame", "b.fits", ValueError("bad value"), prefix="masters: ")
    fail(summary, "group", "main/Mg", ValueError(), kind="configuration_unknown")
    fail(summary, "run", None, LibraryBusyError("library busy"))
    fail(summary, "run", None, LibraryBusyError("library busy"), kind="error")
    a, b, c, d, e = summary.failed
    assert (a.message, a.kind, a.error_type) == ("a plain message", "error", None)
    assert (b.message, b.kind, b.error_type) == ("masters: bad value", "error", "ValueError")
    assert (c.message, c.kind, c.error_type) == ("ValueError", "configuration_unknown", "ValueError")
    assert (d.message, d.kind, d.error_type) == ("library busy", "library_busy", "LibraryBusyError")
    assert e.kind == "error"
