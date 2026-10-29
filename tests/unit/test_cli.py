"""Tests for the command-line entry point."""

import pytest

from mastpipe import __version__
from mastpipe.cli import main


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_no_command_prints_help(capsys):
    assert main([]) == 0
    assert "usage: mastpipe" in capsys.readouterr().out
