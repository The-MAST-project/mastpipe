"""Unit tests of the quick-look result: the outcome and exit code of a run, and the
schema-1 document, compared with the committed example (the contract)."""

import json
from dataclasses import replace
from pathlib import Path

import pytest

from mastpipe.core.errors import LibraryBusyError
from mastpipe.pipeline.group import Failure, RunOptions, RunSummary, Skipped
from mastpipe.pipeline.summary import document, outcome, write


def summary(**changes) -> RunSummary:
    base = RunSummary(started="2026-01-02T03:04:05Z", finished="2026-01-02T03:04:06Z", frames_scanned=3)
    return replace(base, **changes)


@pytest.mark.parametrize(
    ("changes", "error", "expected"),
    [
        ({}, None, ("complete", 0)),
        ({"frames_scanned": 0}, None, ("empty", 0)),
        ({"failed": [Failure("frame", "a.fits", "boom")]}, None, ("partial", 3)),
        ({"frames_scanned": 0, "failed": [Failure("frame", "a.fits", "unreadable FITS file")]}, None, ("partial", 3)),
        ({}, LibraryBusyError("busy"), ("failed", 4)),
        ({}, RuntimeError("bug"), ("error", 1)),
        # once the library is open, a run error leaving the run is a re-raise: error, 1
        ({"library": Path("lib")}, LibraryBusyError("busy"), ("error", 1)),
    ],
)
def test_outcome(changes, error, expected):
    assert outcome(summary(**changes), error) == expected


def test_document_matches_the_committed_example(tmp_path):
    out, library = tmp_path / "out", tmp_path / "lib"
    run = summary(
        products=[out / "spec" / "p_spec.fits", library / "m_bias.fits"],
        qc=[out / "qa" / "p.json"],
        skipped=[Skipped("group", "main/Mg", "not camera x")],
        failed=[Failure("frame", "acq/e1.fits", "boom", error_type="ValueError")],
        reduced_again=["acq/e2.fits"],
        log=out / "logs" / "quicklook_20260102T030405_abcdef.log",
        library=library,
    )
    doc = document(run, RunOptions("highspec", force=True), Path("night"), out, ["extract.halfwidth=4"], None)
    doc.update(input="<input>", out="<out>", library="<library>", mastpipe_version="<version>")
    example = json.loads((Path(__file__).parent / "data" / "summary_schema_1.json").read_text(encoding="utf-8"))
    assert doc == example
    assert doc["products"][1] == {"root": "library", "path": "m_bias.fits"}
    assert doc["log"] == "logs/quicklook_20260102T030405_abcdef.log"
    assert doc["status"] == "partial" and doc["exit_code"] == 3 and doc["schema"] == 1


def test_run_error_document():
    doc = document(RunSummary(), RunOptions("highspec"), Path("n"), Path("o"), [], LibraryBusyError("busy"))
    assert doc["run_error"] == {
        "kind": "library_busy",
        "error_type": "LibraryBusyError",
        "message": "busy",
        "retryable": True,
    }
    assert doc["frames_scanned"] is None and doc["log"] is None and doc["library"] is None
    assert doc["products"] == [] and doc["status"] == "failed" and doc["exit_code"] == 4
    doc = document(RunSummary(), RunOptions("highspec"), Path("n"), Path("o"), [], RuntimeError("bug"))
    assert doc["run_error"]["kind"] == "unexpected" and doc["run_error"]["retryable"] is False
    assert doc["status"] == "error" and doc["exit_code"] == 1


def test_run_error_after_the_library_opened_is_unexpected():
    run = RunSummary(library=Path("lib"), products=[Path("o/spec/p_spec.fits")])
    doc = document(run, RunOptions("highspec", on_error="raise"), Path("n"), Path("o"), [], LibraryBusyError("busy"))
    assert doc["run_error"] == {
        "kind": "unexpected",
        "error_type": "LibraryBusyError",
        "message": "busy",
        "retryable": False,
    }
    assert doc["status"] == "error" and doc["exit_code"] == 1


def test_paths_as_given_with_forward_slashes():
    options = RunOptions("highspec", manifest=Path("some") / "frames.csv")
    run = RunSummary(library=Path("some") / "lib")
    doc = document(run, options, Path("some") / "night", Path("some") / "out", [], None)
    assert doc["input"] == "some/night" and doc["out"] == "some/out"
    assert doc["library"] == "some/lib" and doc["options"]["manifest"] == "some/frames.csv"


def test_root_is_decided_on_absolute_paths(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "work").mkdir()
    monkeypatch.chdir(tmp_path / "work")
    library = Path("..") / "lib"
    run = summary(products=[library / "m_bias.fits", Path("spec") / "p_spec.fits"], library=library)
    doc = document(run, RunOptions("highspec"), Path("night"), Path("."), [], None)
    assert doc["products"] == [{"root": "library", "path": "m_bias.fits"}, {"root": "out", "path": "spec/p_spec.fits"}]


def test_the_deeper_root_wins():
    out = Path("o")
    run = summary(products=[out / "calib" / "m_bias.fits", out / "spec" / "p_spec.fits"], library=out / "calib")
    doc = document(run, RunOptions("highspec"), Path("night"), out, [], None)
    assert doc["products"] == [{"root": "library", "path": "m_bias.fits"}, {"root": "out", "path": "spec/p_spec.fits"}]


def test_write_is_sorted_ascii_json_with_a_newline(tmp_path):
    path = tmp_path / "new" / "r.json"
    write(path, {"b": "é", "a": [1]})
    text = path.read_bytes().decode("ascii")
    assert text == '{\n "a": [\n  1\n ],\n "b": "\\u00e9"\n}\n'
    assert [p.name for p in path.parent.iterdir()] == ["r.json"]


RESULT = {"schema": 1, "command": "quicklook", "status": "complete"}


@pytest.mark.parametrize(
    "existing",
    [
        b"",  # an empty file (a caller's pre-created placeholder)
        json.dumps(RESULT).encode("utf-8"),  # an earlier result document
        b'\n  {"command": "quicklook", "schema": 99}\n',  # another schema's result
    ],
)
def test_write_replaces_an_empty_file_or_an_earlier_result(tmp_path, existing):
    path = tmp_path / "r.json"
    path.write_bytes(existing)
    write(path, {"schema": 1, "command": "quicklook", "status": "empty"})
    assert json.loads(path.read_text(encoding="utf-8"))["status"] == "empty"
    assert [p.name for p in tmp_path.iterdir()] == ["r.json"]


@pytest.mark.parametrize(
    "existing",
    [
        b'[{"kind": "bias", "path": "m_bias.fits"}]\n',  # the library's index.json
        b"SIMPLE  =                    T" + b" " * 2850,  # a FITS file
        b"source,lamp_state\nx.fits,on\n",  # a manifest
        b'{"command": "identify", "schema": 1}',  # a JSON object of another command
        b'{"command": "quicklook"}',  # no schema
        b'{"schema": 1}',  # no command
        b"\xff\xfe{",  # not UTF-8
        b"{not json",
        b"   ",  # not empty, not JSON
    ],
)
def test_write_refuses_to_overwrite_another_file(tmp_path, existing):
    path = tmp_path / "index.json"
    path.write_bytes(existing)
    with pytest.raises(FileExistsError, match="is not a quick-look result file; not overwritten"):
        write(path, {"schema": 1, "command": "quicklook"})
    assert path.read_bytes() == existing
    assert [p.name for p in tmp_path.iterdir()] == ["index.json"]  # no temporary left


def test_write_refuses_a_folder(tmp_path):
    folder = tmp_path / "r.json"
    folder.mkdir()
    with pytest.raises(OSError):
        write(folder, {"schema": 1, "command": "quicklook"})
    assert folder.is_dir() and list(folder.iterdir()) == []


def test_write_refuses_an_alias_of_another_file(tmp_path):
    """The existing content decides, not the spelling: a link to the library index is refused."""
    index = tmp_path / "calib" / "index.json"
    index.parent.mkdir()
    index.write_text("[]\n", encoding="utf-8")
    link = tmp_path / "r.json"
    try:
        link.symlink_to(index)
    except (OSError, NotImplementedError):
        pytest.skip("cannot create a symbolic link here")
    with pytest.raises(FileExistsError, match="not overwritten"):
        write(link, {"schema": 1, "command": "quicklook"})
    assert link.is_symlink() and index.read_text(encoding="utf-8") == "[]\n"
    with pytest.raises(FileExistsError, match="not overwritten"):
        write(tmp_path / "calib" / ".." / "calib" / "index.json", {"schema": 1, "command": "quicklook"})
    assert index.read_text(encoding="utf-8") == "[]\n"


@pytest.mark.parametrize("name", ["index.json", "index.lock", "sub/r.json"])
def test_write_refuses_a_path_in_the_library_even_when_absent(tmp_path, name):
    """An empty library has no index.json yet: a result written there would take its place."""
    library = tmp_path / "out" / "calib"
    path = library / name
    with pytest.raises(FileExistsError, match="lies in the calibration library .*not written"):
        write(path, {"schema": 1, "command": "quicklook"}, library)
    assert not path.exists()


def test_write_refuses_the_library_under_another_spelling(tmp_path):
    library = tmp_path / "calib"
    library.mkdir()
    alias = tmp_path / "alias"
    try:
        alias.symlink_to(library, target_is_directory=True)
    except (OSError, NotImplementedError):
        pytest.skip("cannot create a symbolic link here")
    with pytest.raises(FileExistsError, match="lies in the calibration library"):
        write(alias / "index.json", {"schema": 1, "command": "quicklook"}, library)
    with pytest.raises(FileExistsError, match="lies in the calibration library"):
        write(tmp_path / "x" / ".." / "calib" / "index.json", {"schema": 1, "command": "quicklook"}, library)
    assert list(library.iterdir()) == []


def test_write_beside_the_library_is_allowed(tmp_path):
    library = tmp_path / "out" / "calib"
    path = tmp_path / "out" / "calib-result.json"  # a sibling whose name starts like the library's
    write(path, {"schema": 1, "command": "quicklook"}, library)
    assert json.loads(path.read_text(encoding="utf-8"))["command"] == "quicklook"
