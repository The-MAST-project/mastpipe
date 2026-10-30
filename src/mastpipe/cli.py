"""Command-line entry point.

Sub-commands are added by the registrars listed in :data:`COMMANDS`; each adds
one sub-parser whose ``func`` default runs the command. With no sub-command the
program prints its help and exits successfully.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence

from mastpipe.__version__ import __version__
from mastpipe.pipeline.commands import register_quicklook, register_show, register_solve
from mastpipe.pipeline.identify import register_identify

CommandRegistrar = Callable[[argparse._SubParsersAction], None]
"""Adds one sub-parser to the sub-parsers action and sets its ``func`` default,
a callable taking the parsed ``argparse.Namespace`` and returning the exit code."""

COMMANDS: list[CommandRegistrar] = [
    register_identify,
    register_quicklook,
    register_solve,
    register_show,
]
"""Registrars applied, in order, by :func:`build_parser`."""


def build_parser() -> argparse.ArgumentParser:
    """Return the top-level argument parser with every registered sub-command."""
    parser = argparse.ArgumentParser(
        prog="mastpipe",
        description="Data-reduction pipeline for the MAST fiber-fed spectrographs.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(dest="command", title="commands", metavar="<command>")
    for register in COMMANDS:
        register(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Parse ``argv`` and dispatch; return the process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    return args.func(args)
