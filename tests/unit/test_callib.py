"""Tests for the calibration library: masters found by their validity keys, references
found by validity time or identifier, and the ``index.json`` behind both."""

import json
import multiprocessing
import threading

import numpy as np
import pandas as pd
import pytest

from mastpipe.core.model import LINE_COLUMNS, FrameMeta, MasterCalib, WavelengthSolution
from mastpipe.io.atomic import locked
from mastpipe.io.ingest import file_checksum
from mastpipe.io.products import read_master, reference_name
from mastpipe.pipeline.callib import INDEX_FILE, LOCK_FILE, CalibLibrary
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.wavecal import build_reference
from mastpipe.steps.wavemodel import model_from_config

N_PIX = 1598
SHAPE = (8, 10)


def make_meta(kind: str, mid_time: str, **kwargs) -> FrameMeta:
    fields = {
        "instrument": "highspec",
        "camera": "main",
        "exposure_type": kind,
        "configuration": "Mg",
        "exposure_time": 0.0 if kind == "bias" else 300.0,
        "mid_time": mid_time,
        "readout_mode": "image",
        "amplifier": "conventional",
        "detector_temperature": -60.0,
        "checksum": f"{kind}{mid_time[-8:].replace(':', '')}",
        "source": f"{kind}.fits",
    }
    fields.update(kwargs)
    return FrameMeta(**fields)


def make_master(kind: str, mid_time: str, level: float, **kwargs) -> tuple[MasterCalib, FrameMeta]:
    meta = make_meta(kind, mid_time, **kwargs)
    validity = {
        "camera": meta.camera,
        "binning": meta.binning,
        "readout_mode": meta.readout_mode,
        "amplifier": meta.amplifier,
        "detector_temperature": meta.detector_temperature,
        "mid_time": meta.mid_time,
        "exposure_time": meta.exposure_time,
    }
    master = MasterCalib(
        kind=kind,
        data=np.full(SHAPE, level, dtype=np.float32),
        validity=validity,
        n_frames=3,
        provenance={"inputs": [meta.checksum], "n_frames": 3, "method": "median"},
        variance=np.ones(SHAPE, dtype=np.float32),
    )
    return master, meta


CALIB = "calib0123456"
"""A calib hash for the reuse tests."""


def make_provenance() -> Provenance:
    return Provenance("0.0.1", "cfg123", ["in1"], [])


def make_reference(created: str, x_c_offset: float, highspec_effective):
    geom = grating_geometry(highspec_effective.geometry, highspec_effective.detector.pixel_size_um)
    model = model_from_config("grating", geom, N_PIX)
    params = np.array(model.params, dtype=float)
    params[0] += x_c_offset
    x = np.array([200.5, 800.25, 1400.75])
    lines = pd.DataFrame(
        {
            "x": x,
            "wavelength": np.asarray(model.wavelength(x), dtype=float),
            "species": ["Th I", "Ar I", "Th II"],
            "dx_px": np.zeros(3),
        }
    )[list(LINE_COLUMNS)]
    solution = WavelengthSolution(
        fiber=0,
        kind="grating",
        params=params,
        n_pix=N_PIX,
        lines=lines,
        rms_px=0.01,
        n_lines=3,
        anchor_px=(200.5, 1400.75),
        anchor_nm=(float(lines["wavelength"].min()), float(lines["wavelength"].max())),
        trace_y_center=12.0,
        geometry={
            "d_nm": geom.d_nm,
            "alpha_rad": geom.alpha_rad,
            "order": geom.order,
            "B_nom": geom.B_nom,
        },
    )
    return build_reference("highspec", "Mg", "survey45", {0: solution}, created, ("arc",))


# -------------------------------------------------------------------------------- masters


def test_callib_add_writes_master_and_index(tmp_path):
    library = CalibLibrary(tmp_path / "calib")
    master, meta = make_master("bias", "2026-10-01T18:00:00", 500.0)
    path = library.add(master, make_provenance(), meta)
    assert path == tmp_path / "calib" / "highspec_Mg_main_20261001T180000_bias18_bias.fits"
    assert path.is_file()
    back, prov = read_master(path)
    assert back.kind == "bias" and float(back.data[0, 0]) == 500.0
    assert prov.inputs == ["in1"]
    index = json.loads((tmp_path / "calib" / INDEX_FILE).read_text(encoding="utf-8"))
    assert len(index) == 1 and index[0]["kind"] == "bias" and index[0]["path"] == path.name
    # adding the same master again replaces its entry instead of duplicating it
    library.add(master, make_provenance(), meta)
    assert len(CalibLibrary(tmp_path / "calib").list()) == 1
    assert [entry["kind"] for entry in library.list(kind="bias")] == ["bias"]
    assert library.list(kind="flat") == []


def test_find_by_inputs_matches_sorted_checksums(tmp_path):
    library = CalibLibrary(tmp_path / "calib")
    master, meta = make_master("bias", "2026-10-01T18:00:00", 500.0)
    prov = Provenance("0.0.1", "cfg123", ["c3", "c1", "c2"], [])
    path = library.add(master, prov, meta, CALIB)
    entry = library.list(kind="bias")[0]
    assert entry["inputs"] == ["c1", "c2", "c3"]
    assert entry["file_checksum"] == file_checksum(path)
    assert entry["calib_hash"] == CALIB

    # the same frames in any order find it, read back with its path
    found = library.find_by_inputs("bias", ["c2", "c3", "c1"], "main", CALIB)
    assert found is not None
    back, found_path = found
    assert found_path == path and float(back.data[0, 0]) == 500.0
    # the index as written, seen by a fresh instance, finds it too
    assert CalibLibrary(tmp_path / "calib").find_by_inputs("bias", ["c1", "c2", "c3"], "main", CALIB)
    # another set of frames, another kind or another camera does not
    assert library.find_by_inputs("bias", ["c1", "c2"], "main", CALIB) is None
    assert library.find_by_inputs("bias", ["c1", "c2", "c3", "c4"], "main", CALIB) is None
    assert library.find_by_inputs("dark", ["c1", "c2", "c3"], "main", CALIB) is None
    assert library.find_by_inputs("bias", ["c1", "c2", "c3"], "other", CALIB) is None
    # the same frames stacked under other parameters or another detector do not
    assert library.find_by_inputs("bias", ["c1", "c2", "c3"], "main", "calibfedcba98") is None
    # an empty set of frames never matches
    assert library.find_by_inputs("bias", [], "main", CALIB) is None


def test_find_by_inputs_needs_a_recorded_calib_hash(tmp_path):
    """An entry added without a calib hash (as a library written before it existed) is
    never reused: the master is built again once."""
    library = CalibLibrary(tmp_path / "calib")
    master, meta = make_master("bias", "2026-10-01T18:00:00", 500.0)
    library.add(master, Provenance("0.0.1", "cfg123", ["c1"], []), meta)
    assert library.list(kind="bias")[0]["calib_hash"] is None
    assert library.find_by_inputs("bias", ["c1"], "main", CALIB) is None


def test_find_by_inputs_ignores_a_missing_file(tmp_path):
    library = CalibLibrary(tmp_path / "calib")
    master, meta = make_master("bias", "2026-10-01T18:00:00", 500.0)
    path = library.add(master, Provenance("0.0.1", "cfg123", ["c1"], []), meta, CALIB)
    assert library.find_by_inputs("bias", ["c1"], "main", CALIB) is not None
    path.rename(tmp_path / "moved.fits")
    assert library.find_by_inputs("bias", ["c1"], "main", CALIB) is None


def test_callib_find_matches_keys_temperature_and_time(tmp_path):
    library = CalibLibrary(tmp_path / "calib")
    early, early_meta = make_master("bias", "2026-10-01T18:00:00", 1.0)
    late, late_meta = make_master("bias", "2026-10-02T18:00:00", 2.0)
    warm, warm_meta = make_master("bias", "2026-10-03T18:00:00", 3.0, detector_temperature=-40.0)
    binned, binned_meta = make_master("bias", "2026-10-04T18:00:00", 4.0, binning=(2, 2))
    flat, flat_meta = make_master("flat", "2026-10-01T19:00:00", 5.0, configuration="Halpha")
    prov = make_provenance()
    for master, meta in ((early, early_meta), (late, late_meta), (warm, warm_meta)):
        library.add(master, prov, meta)
    library.add(binned, prov, binned_meta)
    library.add(flat, prov, flat_meta)

    # the library is re-read from disk by a new instance
    library = CalibLibrary(tmp_path / "calib")
    found = library.find("bias", "main", (1, 1), "image/conventional")
    assert found is not None
    master, path = found
    assert float(master.data[0, 0]) == 3.0 and path.name.endswith("_bias.fits")  # the latest
    master, _ = library.find("bias", "main", (1, 1), "image/conventional", temperature=-60.0)
    assert float(master.data[0, 0]) == 2.0  # the warm one is outside the tolerance
    master, _ = library.find("bias", "main", (1, 1), "image/conventional", temperature=-41.0)
    assert float(master.data[0, 0]) == 3.0
    master, _ = library.find(
        "bias",
        "main",
        (1, 1),
        "image/conventional",
        temperature=-60.0,
        near_time="2026-10-01T20:00:00",
    )
    assert float(master.data[0, 0]) == 1.0  # nearest in time among the matches
    master, _ = library.find("bias", "main", (2, 2), "image/conventional")
    assert float(master.data[0, 0]) == 4.0
    assert library.find("bias", "main", (1, 1), "fast/conventional") is None
    assert library.find("bias", "other", (1, 1), "image/conventional") is None
    assert library.find("dark", "main", (1, 1), "image/conventional") is None
    assert library.find("flat", "main", (1, 1), "image/conventional", configuration="Mg") is None
    master, _ = library.find("flat", "main", (1, 1), "image/conventional", configuration="Halpha")
    assert float(master.data[0, 0]) == 5.0


def test_callib_find_separates_amplifiers(tmp_path):
    library = CalibLibrary(tmp_path / "calib")
    bias, meta = make_master("bias", "2026-10-01T18:00:00", 1.0)
    library.add(bias, make_provenance(), meta)
    assert library.list("bias")[0]["readout_key"] == "image/conventional"
    assert library.find("bias", "main", (1, 1), "image/em") is None
    assert library.find("bias", "main", (1, 1), "image/conventional") is not None


def test_callib_find_keys_an_older_entry_by_its_readout_mode(tmp_path):
    library = CalibLibrary(tmp_path / "calib")
    bias, meta = make_master("bias", "2026-10-01T18:00:00", 1.0)
    library.add(bias, make_provenance(), meta)
    index = tmp_path / "calib" / INDEX_FILE
    entries = json.loads(index.read_text(encoding="utf-8"))
    for entry in entries:
        del entry["readout_key"]  # as written before the amplifier was recorded
    index.write_text(json.dumps(entries), encoding="utf-8")
    library = CalibLibrary(tmp_path / "calib")
    # its amplifier is unknown, so it never stands for a known one
    assert library.find("bias", "main", (1, 1), "image/conventional") is None
    assert library.find("bias", "main", (1, 1), "image/unknown") is not None


# ------------------------------------------------------------------------------ references


def test_callib_find_reference_by_time_and_identifier(tmp_path, highspec_effective):
    library = CalibLibrary(tmp_path / "calib")
    earlier = make_reference("2026-10-01T18:00:00", 0.0, highspec_effective)
    later = make_reference("2026-10-05T18:00:00", 2.5, highspec_effective)
    prov = make_provenance()
    path = library.add_reference(later, prov, "main")
    assert path == tmp_path / "calib" / reference_name(later) and path.is_file()
    library.add_reference(earlier, prov, "main")
    assert earlier.identifier != later.identifier

    library = CalibLibrary(tmp_path / "calib")
    between = library.find_reference("highspec", "Mg", "survey45", at_time="2026-10-03T00:00:00")
    assert between is not None and between.identifier == earlier.identifier
    assert between.hash == earlier.hash and set(between.solutions) == {0}
    after = library.find_reference("highspec", "Mg", "survey45", at_time="2026-10-06T00:00:00")
    assert after is not None and after.identifier == later.identifier
    # a time-zone-aware time compares as a time, not as a string
    aware = library.find_reference("highspec", "Mg", "survey45", at_time="2026-10-05T19:00:00+02:00")
    assert aware is not None and aware.identifier == earlier.identifier
    latest = library.find_reference("highspec", "Mg", "survey45")
    assert latest is not None and latest.identifier == later.identifier
    assert library.find_reference("highspec", "Mg", "survey45", at_time="2026-09-01T00:00:00") is None
    by_id = library.find_reference(
        "highspec",
        "Mg",
        "survey45",
        at_time="2026-10-06T00:00:00",
        identifier=earlier.identifier,
    )
    assert by_id is not None and by_id.identifier == earlier.identifier
    assert library.find_reference("highspec", "Mg", "survey45", identifier="nope") is None
    assert library.find_reference("highspec", "Mg", "lab3") is None
    assert library.find_reference("highspec", "Halpha", "survey45") is None
    # the index knows the camera and the source checksums of every reference
    entries = library.list(kind="reference")
    assert len(entries) == 2
    assert {entry["camera"] for entry in entries} == {"main"}
    assert library.reference_for_source("highspec", "Mg", "survey45", "arc") is not None
    assert library.reference_for_source("highspec", "Mg", "survey45", "other") is None


def test_callib_reference_with_same_identifier_is_replaced(tmp_path, highspec_effective):
    library = CalibLibrary(tmp_path / "calib")
    first = make_reference("2026-10-01T18:00:00", 0.0, highspec_effective)
    second = make_reference("2026-10-01T18:00:00", 1.0, highspec_effective)
    assert first.identifier == second.identifier and first.hash != second.hash
    library.add_reference(first, make_provenance(), "main")
    library.add_reference(second, make_provenance(), "main")
    assert len(library.list(kind="reference")) == 1
    found = library.find_reference("highspec", "Mg", "survey45")
    assert found is not None and found.hash == second.hash


def test_callib_rejects_unknown_camera_on_read(tmp_path, highspec_effective):
    library = CalibLibrary(tmp_path / "calib")
    ref = make_reference("2026-10-01T18:00:00", 0.0, highspec_effective)
    library.add_reference(ref, make_provenance(), "nocam")
    with pytest.raises(KeyError):
        library.find_reference("highspec", "Mg", "survey45")


# ---------------------------------------------------------------------------- concurrency

N_PER_WORKER = 25


def add_biases(root: str, worker: str, barrier) -> None:
    """Add ``N_PER_WORKER`` distinct bias masters to the library at ``root``; a spawned
    process's target, so it lives at module level. Both workers open the library before
    either adds, so each starts from the same (empty) index."""
    library = CalibLibrary(root)
    barrier.wait(timeout=60)
    for i in range(N_PER_WORKER):
        master, meta = make_master("bias", f"2026-10-01T18:{i:02d}:00", 500.0, checksum=f"{worker}{i:05d}")
        library.add(master, make_provenance(), meta)


def test_two_processes_add_to_one_library(tmp_path):
    context = multiprocessing.get_context("spawn")
    barrier = context.Barrier(2)
    workers = [context.Process(target=add_biases, args=(str(tmp_path), name, barrier)) for name in ("a", "b")]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join(timeout=120)
    assert [worker.exitcode for worker in workers] == [0, 0]
    assert len(CalibLibrary(tmp_path).list("bias")) == 2 * N_PER_WORKER
    index = json.loads((tmp_path / INDEX_FILE).read_text(encoding="utf-8"))
    assert len(index) == 2 * N_PER_WORKER
    assert not [p.name for p in tmp_path.iterdir() if ".tmp-" in p.name]


def test_opening_a_library_waits_for_a_writer(tmp_path):
    root = tmp_path / "calib"
    opened: list[CalibLibrary] = []
    with locked(root / LOCK_FILE):
        reader = threading.Thread(target=lambda: opened.append(CalibLibrary(root)))
        reader.start()
        reader.join(timeout=0.3)
        assert reader.is_alive() and not opened  # blocked while the index is being written
    reader.join(timeout=10)
    assert len(opened) == 1 and opened[0].list() == []
