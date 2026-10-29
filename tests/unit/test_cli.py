"""Tests for the command-line entry point."""

import pytest

from mastpipe import __version__, cli
from mastpipe.cli import main


def register_probe(subparsers):
    parser = subparsers.add_parser("probe")
    parser.add_argument("value")
    parser.set_defaults(func=lambda args: 3 if args.value == "x" else 4)


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_no_command_prints_help(capsys):
    assert main([]) == 0
    assert "usage: mastpipe" in capsys.readouterr().out


def test_registered_command_dispatches(monkeypatch):
    monkeypatch.setattr(cli, "COMMANDS", [register_probe])
    assert main(["probe", "x"]) == 3
    assert main(["probe", "y"]) == 4
