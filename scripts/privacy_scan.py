"""Scan tracked text files for private information that must never be committed.

Checked patterns: local absolute paths, e-mail addresses, IPv4 addresses (loopback and
the unspecified address are allowed) and credential tokens. The scan covers every file
listed by ``git ls-files`` whose suffix marks it as text, skips this script itself,
prints one line per hit and exits with status 1 when anything is found.

Usage: ``python scripts/privacy_scan.py`` from anywhere inside the repository.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

PATTERNS: dict[str, re.Pattern[str]] = {
    "local path": re.compile(r"/mnt/|/home/|/Users/|[A-Za-z]:\\Users\\"),
    "e-mail": re.compile(r"\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}\b"),
    "ipv4": re.compile(r"\b(?:\d{1,3}\.){3}\d{1,3}\b"),
    "token": re.compile(r"\bghp_[A-Za-z0-9]{20,}|\bAKIA[0-9A-Z]{16}\b|\bsk-[A-Za-z0-9_-]{20,}"),
}

ALLOWED_IPV4 = {"0.0.0.0", "127.0.0.1"}

TEXT_SUFFIXES = {
    "",
    ".bat",
    ".cfg",
    ".csv",
    ".dat",
    ".ini",
    ".ipynb",
    ".json",
    ".md",
    ".ps1",
    ".py",
    ".rst",
    ".sh",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}

MAX_SHOWN = 12


def repository_root() -> Path:
    """Return the root of the git repository that contains this script."""
    result = subprocess.run(
        ["git", "rev-parse", "--show-toplevel"],
        cwd=Path(__file__).resolve().parent,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=True,
    )
    return Path(result.stdout.strip())


def tracked_text_files(root: Path) -> list[Path]:
    """Return tracked files under ``root`` whose suffix marks them as text."""
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=True,
    )
    this_script = Path(__file__).resolve()
    files = []
    for name in result.stdout.split("\0"):
        if not name:
            continue
        path = root / name
        if path.suffix.lower() not in TEXT_SUFFIXES or path.resolve() == this_script:
            continue
        if path.is_file():
            files.append(path)
    return files


def redact(match: str) -> str:
    """Shorten a matched string so that a real secret is not echoed in full."""
    if len(match) <= MAX_SHOWN:
        return match
    return match[:MAX_SHOWN] + "..."


def scan_file(path: Path) -> list[tuple[int, str, str]]:
    """Return ``(line number, pattern name, redacted match)`` for every hit in ``path``."""
    hits = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for lineno, line in enumerate(handle, start=1):
            for name, pattern in PATTERNS.items():
                for match in pattern.finditer(line):
                    text = match.group(0)
                    if name == "ipv4" and text in ALLOWED_IPV4:
                        continue
                    hits.append((lineno, name, redact(text)))
    return hits


def main() -> int:
    root = repository_root()
    total = 0
    for path in tracked_text_files(root):
        for lineno, name, text in scan_file(path):
            total += 1
            print(f"{path.relative_to(root).as_posix()}:{lineno}: {name}: {text}")
    if total:
        sys.stdout.flush()
        print(f"privacy scan: {total} hit(s) found", file=sys.stderr)
        return 1
    print("privacy scan: clean")
    return 0


if __name__ == "__main__":
    sys.exit(main())
