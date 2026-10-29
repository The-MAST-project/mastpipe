"""Command-line entry point.

Sub-commands are registered on the parser returned by :func:`build_parser`;
with no sub-command the program prints its help and exits successfully.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from mastpipe.__version__ import __version__


def build_parser() -> argparse.ArgumentParser:
    """Return the top-level argument parser with its sub-parser group."""
    parser = argparse.ArgumentParser(
        prog="mastpipe",
        description="Data-reduction pipeline for the MAST fiber-fed spectrographs.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    parser.add_subparsers(dest="command", title="commands", metavar="<command>")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Parse ``argv`` and dispatch; return the process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
    return 0
