"""Regression of the HighSpec path on the laboratory frames reduced in the HighSpec laboratory
notebook: the three-fiber Mg b layout (``lab3``), the bias, lamp-off background and flat of
the calibration folder, and the Th-Ar focus-scan frame at best focus as the full arc.

The frames live outside the repository under ``$MASTPIPE_TEST_DATA/highspec/lab3/``
(``calib/`` and ``scan/``, see the README) and the module is skipped when the variable is
unset. The Andor headers carry the exposure time, temperature, binning and readout mode, and
``IMAGETYP`` (``bias`` for the bias, ``light`` for every other frame), but no lamp, grating
or fiber preset. The manifest, as a laboratory run of ``mastpipe quicklook --manifest``
would give it, declares only what the content cannot decide: the grating and the fiber
preset of every frame, and the type of the flat (a continuum, which content calls
science) and of the lamp-off background (no signal). The full arc is left to the runner,
which resolves it from its content. The numbers asserted are the notebook's:
trace centres, per-fiber rms and line counts, band coverage and the fiber-to-fiber
agreement of the reference solution.

One test reads the frames without the manifest, as the control software writes them
(``IMAGETYP = 'light'``), and checks that their content resolves them: every focus frame of
the scan an arc, the flat a continuum, the lamp-off background no signal.

The runs are quick (one arc, three fibers), but ``quicklook`` -> ``identify`` ->
``quicklook`` is one story, so a module-scoped output folder carries it and the tests
read it in order.
"""

import os
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from mastpipe.config.loader import build_effective
from mastpipe.instruments.registry import load_header_map, load_instrument, load_layout
from mastpipe.io.ingest import file_checksum, load_frame, scan_folder
from mastpipe.io.products import read_qc, read_reference
from mastpipe.pipeline.callib import CalibLibrary, registry_geometry
from mastpipe.pipeline.identify import AGREEMENT_FILE, CALIB_FOLDER, _arc_lines, run_identify
from mastpipe.pipeline.runner import RunOptions, run_quicklook
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.classify import classify_light
from mastpipe.steps.drift import fit_fiber_drift
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.trace import find_traces, match_layout
from mastpipe.steps.wavecal import fiber_agreement
from mastpipe.steps.wavemodel import model_from_solution

DATA_VARIABLE = "MASTPIPE_TEST_DATA"
ATLAS_VARIABLE = "MASTPIPE_ATLAS_DIR"

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        not os.environ.get(DATA_VARIABLE),
        reason=f"set {DATA_VARIABLE} to the folder holding highspec/lab3/ (see the README)",
    ),
]

LAB_FOLDER = Path("highspec") / "lab3"
CALIB_FRAMES = {"image_bias.fits": "bias", "background.fits": "background", "flat_518.fits": "flat"}
"""The calibration frames of the laboratory folder and their types."""
FULL_ARC = "FOCUS_05_13.30004_mm.fits"
"""The Th-Ar frame at the best focus of the scan: the full arc."""
MANIFEST_TYPES = {"image_bias.fits": "", "background.fits": "background", "flat_518.fits": "flat", FULL_ARC: ""}
"""The type the manifest gives each staged frame: none for the bias (its header says
``bias``) and the full arc (resolved from its content, as is its arc kind)."""
CONFIGURATION = "Mg"
PRESET = "lab3"

# The notebook's numbers on the same frames.
NOTEBOOK_Y_CENTERS = (17.40, 188.96, 360.72)
"""Trace centres at the detector centre column, cropped frame, px."""
NOTEBOOK_RMS_PX = (0.016, 0.038, 0.010)
NOTEBOOK_N_LINES = (23, 25, 19)
NOTEBOOK_COVERAGE_NM = (513.498, 527.315)

Y_CENTER_TOL_PX = 0.3
TRACE_RMS_MAX_PX = 0.01
SOLUTION_RMS_MAX_PX = 0.05
COVERAGE_NM = (513.50, 527.32)
COVERAGE_TOL_NM = 0.03
MIN_LINES = 19
SPREAD_MEDIAN_MAX_PX = 0.03
SELF_DRIFT_MAX_PX = 0.01


# ------------------------------------------------------------------------------ fixtures


@pytest.fixture(scope="module")
def lab_root() -> Path:
    root = Path(os.environ[DATA_VARIABLE]) / LAB_FOLDER
    missing = [str(Path("calib") / name) for name in CALIB_FRAMES if not (root / "calib" / name).is_file()]
    if not (root / "scan" / FULL_ARC).is_file():
        missing.append(str(Path("scan") / FULL_ARC))
    if missing:
        pytest.fail(f"{DATA_VARIABLE} lacks {missing} under {LAB_FOLDER}")
    return root


@pytest.fixture(scope="module")
def atlas_dir() -> Path | None:
    value = os.environ.get(ATLAS_VARIABLE)
    return Path(value) if value else None


@pytest.fixture(scope="module")
def night(lab_root: Path, tmp_path_factory: pytest.TempPathFactory) -> SimpleNamespace:
    """The staged night: the three calibration frames, the full arc and the manifest."""
    folder = tmp_path_factory.mktemp("lab_night")
    for name in CALIB_FRAMES:
        shutil.copy(lab_root / "calib" / name, folder / name)
    shutil.copy(lab_root / "scan" / FULL_ARC, folder / FULL_ARC)
    rows = ["source,exposure_type,configuration,fiber_preset"]
    rows += [f"{name},{kind},{CONFIGURATION},{PRESET}" for name, kind in MANIFEST_TYPES.items()]
    manifest = folder.parent / f"{folder.name}_manifest.csv"
    manifest.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="")
    return SimpleNamespace(dir=folder, manifest=manifest, arc=folder / FULL_ARC)


@pytest.fixture(scope="module")
def effective():
    return build_effective(load_instrument("highspec"), CONFIGURATION, "main", layout=PRESET)


@pytest.fixture(scope="module")
def out_dir(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("lab_out")


def quicklook(night: SimpleNamespace, out_dir: Path):
    return run_quicklook(night.dir, out_dir, RunOptions("highspec", manifest=night.manifest))


@pytest.fixture(scope="module")
def first_run(night, out_dir):
    """The night reduced before any reference exists: masters and traces only."""
    return quicklook(night, out_dir)


@pytest.fixture(scope="module")
def reference(first_run, night, effective, out_dir, atlas_dir) -> SimpleNamespace:
    """``mastpipe identify --stored-anchors`` on the full arc, with the masters of the run."""
    ref, path = run_identify(night.arc, effective, PRESET, out_dir, anchors=None, anchor_trace=None, atlas_dir=atlas_dir)
    return SimpleNamespace(ref=ref, path=path)


@pytest.fixture(scope="module")
def second_run(reference, night, out_dir) -> SimpleNamespace:
    """The night reduced again once the reference exists, with the calibration folder's
    listing before and after."""
    before = sorted(p.name for p in (out_dir / CALIB_FOLDER).iterdir())
    summary = quicklook(night, out_dir)
    after = sorted(p.name for p in (out_dir / CALIB_FOLDER).iterdir())
    return SimpleNamespace(summary=summary, before=before, after=after)


# ------------------------------------------------------------------------------- helpers


def coverage_nm(solution, geom) -> tuple[float, float]:
    """Wavelengths at the first and last pixel, in increasing order."""
    model = model_from_solution(solution, geom)
    ends = np.asarray(model.wavelength(np.array([0.0, solution.n_pix - 1.0])), dtype=float)
    return float(ends.min()), float(ends.max())


def arc_y_centers(traceset) -> list[float]:
    return [trace.y_center for trace in traceset.by_role("arc")]


# --------------------------------------------------------------------------------- tests


def test_lab_traces_match_notebook(night, effective, first_run, out_dir):
    """Trace finding on the calibrated flat with the HighSpec default ``TraceParams``.

    The flat is calibrated as the runner does it (bias subtracted, the lamp-off background
    kept: ``flat_subtract_background`` is off, as in the notebook). The three fibers are
    found at the notebook's centres with quadratic fits tighter than 0.01 px. The same frame
    also shows two faint features (S/N 6-10 over the noise of the median) where the ambient
    light of the background frame has structure; ``match_layout`` sets them aside as
    ``extra`` and they are absent from the background-subtracted flat, which yields the
    three fibers alone.

    The bias is typed by its header; the full arc, untyped by the manifest, reaches the run
    as ``light`` and the runner resolves it from its content into an arc
    (``content:arc_like``), which the arc stage takes as the full arc (skipped with "run
    identify" before a reference exists).
    """
    detector, params = effective.detector, effective.params
    refs = {
        ref.meta.source: ref
        for ref in scan_folder(night.dir, "highspec", load_header_map("highspec"), manifest=night.manifest)
    }
    assert {name: refs[name].meta.exposure_type for name in CALIB_FRAMES} == CALIB_FRAMES
    assert refs["image_bias.fits"].meta.classification_rule == "header"
    assert (refs[FULL_ARC].meta.exposure_type, refs[FULL_ARC].meta.arc_kind) == ("light", None)
    log = "".join(path.read_text(encoding="utf-8") for path in (out_dir / "logs").glob("quicklook_*.log"))
    assert f"light {FULL_ARC}: arc (content:arc_like)" in log
    assert f"arc {FULL_ARC}: full" in log
    assert any(item.startswith(FULL_ARC) and "run identify" in item for item in first_run.skipped)
    masters = {
        kind: build_master([load_frame(refs[name], detector)], kind, params.calib)
        for name, kind in CALIB_FRAMES.items()
        if kind != "flat"
    }
    calib = build_calibration_set(masters["bias"], None, masters["background"], detector)
    raw_flat = load_frame(refs["flat_518.fits"], detector)
    layout = load_layout("highspec", PRESET)

    flat = calibrate_frame(raw_flat, calib, params.calib.flat_subtract_background, saturation=detector.saturation)
    traceset = match_layout(find_traces(flat.data, params.trace), layout, params.trace.layout_tolerance)
    assert traceset.missing == ()
    arcs = traceset.by_role("arc")
    assert len(arcs) == 3
    assert np.allclose(arc_y_centers(traceset), NOTEBOOK_Y_CENTERS, atol=Y_CENTER_TOL_PX)
    assert all(trace.rms < TRACE_RMS_MAX_PX for trace in arcs), [t.rms for t in arcs]
    assert all(len(trace.coeffs) == params.trace.degree + 1 for trace in arcs)

    subtracted = calibrate_frame(raw_flat, calib, True, saturation=detector.saturation)
    alone = find_traces(subtracted.data, params.trace)
    assert len(alone) == 3
    assert np.allclose([t.y_center for t in alone], arc_y_centers(traceset), atol=0.1)


def test_lab_focus_frames_are_arc_like(lab_root, effective):
    """Without a manifest every laboratory frame but the bias is ``light`` (the control
    software's ``IMAGETYP``, no lamp card), and the content resolves it as the runner
    would: every Th-Ar focus frame of the scan is an arc, whether or not the lamp-off
    background is subtracted (the ambient light it holds is about 30 times the line flux
    and has broad structure of its own); the flat is a continuum, ``science``; the lamp-off
    background, calibrated with the bias alone as on a night where it is not declared, has
    no signal."""
    detector, params = effective.detector, effective.params
    header_map = load_header_map("highspec")
    calib_refs = {ref.meta.source: ref for ref in scan_folder(lab_root / "calib", "highspec", header_map)}
    focus = sorted(scan_folder(lab_root / "scan", "highspec", header_map), key=lambda r: r.meta.source)
    assert len(focus) >= 11
    for ref in [*focus, calib_refs["flat_518.fits"], calib_refs["background.fits"]]:
        assert (ref.meta.exposure_type, ref.meta.classification_rule) == ("light", "header")
    assert calib_refs["image_bias.fits"].meta.exposure_type == "bias"

    def master(name: str, kind: str):
        return build_master([load_frame(calib_refs[name], detector)], kind, params.calib)

    bias, background = master("image_bias.fits", "bias"), master("background.fits", "background")
    bias_only = build_calibration_set(bias, None, None, detector)
    with_background = build_calibration_set(bias, None, background, detector)

    def resolved(ref, calib) -> tuple[str, str]:
        subtract = calib.background is not None
        frame = calibrate_frame(load_frame(ref, detector), calib, subtract, saturation=detector.saturation)
        return classify_light(frame, params.classify)

    for ref in focus:
        for calib in (bias_only, with_background):
            assert resolved(ref, calib) == ("arc", "content:arc_like"), ref.meta.source
    for calib in (bias_only, with_background):
        assert resolved(calib_refs["flat_518.fits"], calib) == ("science", "content:continuum")
    assert resolved(calib_refs["background.fits"], bias_only) == ("unknown", "content:no_signal")


def test_lab_reference_from_stored_anchors(reference, night, effective, out_dir):
    """``run_identify`` with the stored anchor pairs of ``Mg`` / ``lab3`` on the full arc,
    with the masters the quick-look run left in ``out/calib``: three solutions at the
    notebook's quality (rms 0.016 / 0.038 / 0.010 px, 23 / 25 / 19 lines, coverage
    513.498-527.315 nm) and fiber-to-fiber agreement below 0.03 px."""
    ref, path = reference.ref, reference.path
    geom = grating_geometry(effective.geometry, effective.detector.pixel_size_um)
    assert ref.instrument == "highspec" and ref.configuration == CONFIGURATION
    assert ref.fiber_preset == PRESET
    assert ref.source_checksums == (file_checksum(night.arc),)
    assert sorted(ref.solutions) == [0, 1, 2]  # the layout's trace ids
    for fiber, solution in sorted(ref.solutions.items()):
        assert solution.rms_px <= SOLUTION_RMS_MAX_PX, (fiber, solution.rms_px)
        assert solution.n_lines >= MIN_LINES, (fiber, solution.n_lines)
        low, high = coverage_nm(solution, geom)
        assert low == pytest.approx(COVERAGE_NM[0], abs=COVERAGE_TOL_NM), (fiber, low)
        assert high == pytest.approx(COVERAGE_NM[1], abs=COVERAGE_TOL_NM), (fiber, high)
    assert np.allclose(
        sorted(s.trace_y_center for s in ref.solutions.values()),
        NOTEBOOK_Y_CENTERS,
        atol=Y_CENTER_TOL_PX,
    )
    agreement = fiber_agreement(ref.solutions)
    assert len(agreement) >= MIN_LINES
    assert float(agreement["spread_px"].median()) < SPREAD_MEDIAN_MAX_PX
    # the product is on disk, verifies on reading and is in the library
    assert path.is_file() and (out_dir / CALIB_FOLDER / AGREEMENT_FILE).is_file()
    stored, prov = read_reference(path, registry_geometry("main"))
    assert stored.identifier == ref.identifier and stored.hash == ref.hash
    assert prov.wavelength_status == "solved"
    assert any("anchors on trace 0" in step for step in prov.steps)
    library = CalibLibrary(out_dir / CALIB_FOLDER)
    found = library.find_reference("highspec", CONFIGURATION, PRESET)
    assert found is not None and found.identifier == ref.identifier


def test_lab_drift_against_itself_is_zero(reference, night, effective, out_dir):
    """The ``use_x`` line centres of the full arc fitted for drift against the reference
    made from the very same frame: no shift, no stretch, every identified line matched."""
    ref = reference.ref
    lines = _arc_lines(night.arc, effective, load_layout("highspec", PRESET), out_dir)
    # the two ambient-light features of the flat are dropped: layout fibers only, by id
    assert [trace.index for trace in lines.traceset] == [0, 1, 2]
    assert len(lines.traceset.extra) == 2
    fitted = 0
    for trace, table in zip(lines.traceset, lines.tables, strict=True):
        solution = ref.solutions[trace.index]
        x_obs = table.loc[table["use_x"].astype(bool), "x_center"].to_numpy(dtype=float)
        result = fit_fiber_drift(x_obs, solution, effective.params.drift)
        assert abs(result["shift_px"]) < SELF_DRIFT_MAX_PX, (trace.index, result)
        assert result["stretch"] == pytest.approx(0.0, abs=1e-9)
        assert result["n_lines"] == solution.n_lines
        assert result["rms_px"] < SELF_DRIFT_MAX_PX and result["usable"]
        fitted += 1
    assert fitted == 3


def test_lab_quicklook_runs(first_run, reference, second_run, night, out_dir):
    """``run_quicklook`` on the staged night with the manifest: before ``identify`` the
    masters and traces are made and the full arc is skipped with "run identify"; after
    it, the second run finds the arc already solved into the reference, adds no file to
    the calibration folder and fails on nothing."""
    assert first_run.failed == []
    assert first_run.reduced_again == []
    assert any(item.startswith(FULL_ARC) and "run identify" in item for item in first_run.skipped), first_run.skipped
    assert not any(p.name.endswith("_spec.fits") for p in first_run.products)
    library = CalibLibrary(out_dir / CALIB_FOLDER)
    assert {e["kind"] for e in library.list()} >= {"bias", "background", "flat"}
    traces = sorted((out_dir / CALIB_FOLDER).glob("*_traces.fits"))
    assert len(traces) == 1
    qc = read_qc(out_dir / "qa" / traces[0].with_suffix(".json").name)
    assert qc.stage == "trace" and qc.values["source"] == "flat"
    assert qc.values["layout"] == PRESET and qc.values["missing"] == []
    assert qc.values["index"] == [0, 1, 2] and qc.values["n_traces"] == 3
    assert qc.values["n_detected"] == 5 and len(qc.values["extra_y_center"]) == 2
    rows = zip(qc.values["y_center"], qc.values["roles"], strict=True)
    arc_rows = [y for y, role in rows if role == "arc"]
    assert np.allclose(arc_rows, NOTEBOOK_Y_CENTERS, atol=Y_CENTER_TOL_PX)

    summary = second_run.summary
    assert summary.failed == []
    assert summary.reduced_again == []
    assert any(
        item.startswith(FULL_ARC) and f"already solved into reference {reference.ref.identifier}" in item
        for item in summary.skipped
    ), summary.skipped
    assert not any(p.name.endswith("_spec.fits") for p in summary.products)
    assert second_run.after == second_run.before
    assert any(name.endswith("_wref.fits") for name in second_run.after)
    assert len(library.list(kind="reference")) == 1
    # one log per run, unless both runs started within the same second
    logs = sorted((out_dir / "logs").glob("quicklook_*.log"))
    assert 1 <= len(logs) <= 2
    text = "".join(log.read_text(encoding="utf-8") for log in logs)
    assert "run identify" in text and "already solved" in text
