"""Tests for the atomic file replacement and the advisory lock of ``mastpipe.io.atomic``."""

import errno
import sys

import pytest

from mastpipe.io.atomic import atomic_write, locked


def test_atomic_write_replaces_the_target(tmp_path):
    target = tmp_path / "sub" / "index.json"
    with atomic_write(target) as tmp:
        assert tmp.parent == target.parent and tmp.name.startswith("index.json.tmp-")
        tmp.write_text("new\n", encoding="utf-8")
        assert not target.exists()
    assert target.read_text(encoding="utf-8") == "new\n"
    with atomic_write(target) as tmp:
        tmp.write_text("newer\n", encoding="utf-8")
    assert target.read_text(encoding="utf-8") == "newer\n"
    assert sorted(p.name for p in target.parent.iterdir()) == ["index.json"]


def test_atomic_write_leaves_no_partial_file(tmp_path):
    target = tmp_path / "product.fits"
    with pytest.raises(RuntimeError, match="interrupted"), atomic_write(target) as tmp:
        tmp.write_text("partial", encoding="utf-8")
        raise RuntimeError("interrupted")
    assert not target.exists()
    assert list(tmp_path.iterdir()) == []
    # an earlier complete file survives an interrupted rewrite untouched
    target.write_text("complete", encoding="utf-8")
    with pytest.raises(RuntimeError), atomic_write(target) as tmp:
        tmp.write_text("partial", encoding="utf-8")
        raise RuntimeError("interrupted")
    assert target.read_text(encoding="utf-8") == "complete"
    assert [p.name for p in tmp_path.iterdir()] == ["product.fits"]


def test_overlapping_writes_to_one_target_do_not_share_a_temporary(tmp_path):
    # Two writers in one process (same pid, e.g. two threads) overlap on one target: each
    # must get its own temporary, so each publishes its own contents and neither move fails.
    target = tmp_path / "index.json"
    with atomic_write(target) as first, atomic_write(target) as second:
        assert first != second
        first.write_text("first\n", encoding="utf-8")
        second.write_text("second\n", encoding="utf-8")
    assert target.read_text(encoding="utf-8") == "first\n"  # the outer block finished last
    assert [p.name for p in tmp_path.iterdir()] == ["index.json"]


def test_lock_times_out_with_a_message(tmp_path):
    with (
        locked(tmp_path / "index.lock"),
        pytest.raises(TimeoutError, match="index.lock"),
        locked(tmp_path / "index.lock", timeout_s=0.2),
    ):
        pass
    # released on exit: it can be taken again, and the lock file stays
    with locked(tmp_path / "index.lock", timeout_s=0.2):
        pass
    assert (tmp_path / "index.lock").is_file()


def test_lock_error_that_is_not_busy_propagates(tmp_path, monkeypatch):
    def unsupported(*args):
        raise OSError(errno.ENOLCK, "no locks available")

    if sys.platform == "win32":
        import msvcrt

        monkeypatch.setattr(msvcrt, "locking", unsupported)
    else:
        import fcntl

        monkeypatch.setattr(fcntl, "flock", unsupported)
    with pytest.raises(OSError) as raised, locked(tmp_path / "index.lock", timeout_s=5.0):
        pass
    assert raised.type is OSError and raised.value.errno == errno.ENOLCK


def test_atomic_write_unlink_failure_does_not_mask_the_error(tmp_path, monkeypatch):
    def failing_unlink(self, missing_ok=False):
        raise PermissionError("cannot remove")

    monkeypatch.setattr(type(tmp_path), "unlink", failing_unlink)
    with pytest.raises(RuntimeError, match="interrupted"), atomic_write(tmp_path / "product.fits") as tmp:
        tmp.write_text("partial", encoding="utf-8")
        raise RuntimeError("interrupted")
