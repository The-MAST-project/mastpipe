"""Tests for the atomic file replacement and the advisory lock of ``mastpipe.io.atomic``."""

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
    with pytest.raises(RuntimeError, match="interrupted"):
        with atomic_write(target) as tmp:
            tmp.write_text("partial", encoding="utf-8")
            raise RuntimeError("interrupted")
    assert not target.exists()
    assert list(tmp_path.iterdir()) == []
    # an earlier complete file survives an interrupted rewrite untouched
    target.write_text("complete", encoding="utf-8")
    with pytest.raises(RuntimeError):
        with atomic_write(target) as tmp:
            tmp.write_text("partial", encoding="utf-8")
            raise RuntimeError("interrupted")
    assert target.read_text(encoding="utf-8") == "complete"
    assert [p.name for p in tmp_path.iterdir()] == ["product.fits"]


def test_lock_times_out_with_a_message(tmp_path):
    with locked(tmp_path / "index.lock"):
        with pytest.raises(TimeoutError, match="index.lock"):
            with locked(tmp_path / "index.lock", timeout_s=0.2):
                pass
    # released on exit: it can be taken again, and the lock file stays
    with locked(tmp_path / "index.lock", timeout_s=0.2):
        pass
    assert (tmp_path / "index.lock").is_file()
