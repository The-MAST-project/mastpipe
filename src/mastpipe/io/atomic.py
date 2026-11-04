"""Atomic file replacement and an exclusive advisory file lock.

:func:`atomic_write` lets a writer fill a temporary file next to the target and moves it
into place with :func:`os.replace` only when the writer finishes, so a reader sees either
the old file or the new one, never a partial one. :func:`locked` serialises a
read-modify-write across processes (the calibration library's index). Both work on
POSIX and Windows; the lock is the only platform-specific code, guarded by
``sys.platform``.
"""

from __future__ import annotations

import os
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import IO

LOCK_POLL_S = 0.05
"""Seconds between attempts to take a lock that another holder has."""


@contextmanager
def atomic_write(path: Path | str) -> Iterator[Path]:
    """Yield a temporary path next to ``path``; on normal exit move it onto ``path``.

    The temporary is ``<name>.tmp-<pid>`` in the same folder (so the move is a rename on
    one file system); the folder is created when missing. When the block raises, the
    temporary is removed, ``path`` is left as it was and the error propagates. The writer
    must have closed the temporary before the block ends (Windows refuses to replace a
    file that is open).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp-{os.getpid()}")
    try:
        yield tmp
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


def _try_lock(handle: IO[bytes]) -> bool:
    handle.seek(0)
    try:
        if sys.platform == "win32":
            import msvcrt

            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl

            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        return False
    return True


def _unlock(handle: IO[bytes]) -> None:
    handle.seek(0)
    if sys.platform == "win32":
        import msvcrt

        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
    else:
        import fcntl

        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


@contextmanager
def locked(path: Path | str, timeout_s: float = 30.0) -> Iterator[None]:
    """Hold an exclusive advisory lock on ``path`` for the duration of the block.

    The file is created when missing (with its folder) and never deleted. The lock is
    tried every ``LOCK_POLL_S`` seconds; ``TimeoutError`` when it is not free within
    ``timeout_s``. It is per open file, so a second ``locked`` on the same path waits
    even in the same process. Advisory: it only excludes other users of ``locked``.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+b") as handle:
        deadline = time.monotonic() + timeout_s
        while not _try_lock(handle):
            if time.monotonic() >= deadline:
                raise TimeoutError(f"could not lock {path} within {timeout_s} s")
            time.sleep(LOCK_POLL_S)
        try:
            yield
        finally:
            _unlock(handle)
