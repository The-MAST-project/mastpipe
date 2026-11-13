"""Unit tests of the quick-look runner's pure helpers: the drift models a frame may be
paired with (those fitted against its own reference), the blend of two drift models
under the ``interpolate`` policy, the group's binning and the per-frame isolation of the
calibration frames."""

import logging
import re
import shutil
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.config.models import DriftParams
from mastpipe.core.errors import LibraryBusyError, LibraryInvalidError, OutputUnwritableError
from mastpipe.core.model import (
    PER_FIBER_COLUMNS,
    CalibrationSet,
    DriftModel,
    FrameMeta,
    FrameRef,
    ReferenceSolution,
)
from mastpipe.instruments.registry import load_instrument
from mastpipe.io.products import read_master, write_master
from mastpipe.pipeline import masters
from mastpipe.pipeline.arcs import pair_arcs
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.group import (
    Failure,
    Skipped,
    _calibrated,
    _configuration_origin,
    _configuration_steps,
    _group_binning,
    _group_readout,
)
from mastpipe.pipeline.masters import _masters_hash, _usable_frames
from mastpipe.pipeline.runner import RunOptions, RunSummary, run_quicklook
from mastpipe.pipeline.science import _arc_weights, _arcs_for, _blend_drifts
from mastpipe.vocab import GRATING_NAMES

N_PIX = 1598


def make_meta(mid_time: str, checksum: str, exposure_type: str = "arc") -> FrameMeta:
    return FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type=exposure_type,
        configuration="Mg",
        mid_time=mid_time,
        checksum=checksum,
        source=f"{checksum}.fits",
    )


def make_reference(identifier: str, valid_from: str) -> ReferenceSolution:
    return ReferenceSolution(
        instrument="highspec",
        configuration="Mg",
        fiber_preset="survey45",
        solutions={},
        created=valid_from,
        valid_from=valid_from,
        source_checksums=(identifier,),
        identifier=identifier,
        hash="h" + identifier,
    )


def make_drift(
    checksum: str,
    mid_time: str,
    shift_coeffs: list[float],
    stretch_coeffs: list[float],
    fibers: dict[int, tuple[float, float, int, float, bool]],
    scatter_px: float = 0.02,
    usable: bool = True,
) -> DriftModel:
    rows = [
        {
            "fiber": fiber,
            "y_center": 10.0 + 8.5 * fiber,
            "shift_px": shift,
            "stretch": stretch,
            "n_lines": n_lines,
            "rms_px": rms,
            "usable": ok,
        }
        for fiber, (shift, stretch, n_lines, rms, ok) in fibers.items()
    ]
    return DriftModel(
        arc_checksums=(checksum,),
        arc_mid_time=mid_time,
        per_fiber=pd.DataFrame(rows, columns=list(PER_FIBER_COLUMNS)),
        shift_coeffs=np.array(shift_coeffs),
        stretch_coeffs=np.array(stretch_coeffs),
        scatter_px=scatter_px,
        usable=usable,
        n_pix=N_PIX,
    )


# ------------------------------------------------------------- arcs tied to their reference


@pytest.mark.parametrize("policy", ["nearest", "interpolate", "nearest_before"])
def test_frame_never_receives_drift_fitted_against_another_reference(policy):
    """A full arc solved mid-night makes reference B; the arcs after it are fitted against
    B. A science frame before that time uses A and must pair with A's arcs only, however
    close B's arcs are in time."""
    ref_a = make_reference("A", "2026-10-01T18:00:00")
    ref_b = make_reference("B", "2026-10-01T21:00:00")
    fibers = {0: (0.3, 0.0, 20, 0.02, True), 44: (0.3, 0.0, 20, 0.02, True)}
    arc_a1 = make_drift("a1", "2026-10-01T19:00:00", [0.3], [0.0], fibers)
    arc_b1 = make_drift("b1", "2026-10-01T21:10:00", [-2.0], [0.0], fibers)
    arc_b2 = make_drift("b2", "2026-10-01T21:30:00", [-2.0], [0.0], fibers)
    arcs = [
        (make_meta("2026-10-01T19:00:00", "a1"), arc_a1, "A"),
        (make_meta("2026-10-01T21:10:00", "b1"), arc_b1, "B"),
        (make_meta("2026-10-01T21:30:00", "b2"), arc_b2, "B"),
    ]
    # the frame sits just before B's validity, between an A arc and the B arcs
    frame = make_meta("2026-10-01T20:55:00", "sci", exposure_type="science")
    params = DriftParams(time_policy=policy)

    mine = _arcs_for(arcs, ref_a)
    assert [meta.checksum for meta, _ in mine] == ["a1"]
    pairs = pair_arcs(frame, mine, ref_a, params)
    assert [model.arc_checksums for model, _ in pairs] == [("a1",)]
    assert pairs[0][1] == pytest.approx(1.0)
    # without the filter the nearest / interpolate policies would have taken the B arc
    # fifteen minutes away (nearest_before alone is safe: no B arc precedes the frame)
    unfiltered = pair_arcs(frame, [(m, d) for m, d, _ in arcs], ref_a, params)
    taken = sorted(c for model, _ in unfiltered for c in model.arc_checksums)
    assert taken == (["a1"] if policy == "nearest_before" else ["a1", "b1"] if policy == "interpolate" else ["b1"])
    # a frame after B's time uses B and only B's arcs
    later = make_meta("2026-10-01T21:20:00", "sci2", exposure_type="science")
    theirs = _arcs_for(arcs, ref_b)
    assert sorted(meta.checksum for meta, _ in theirs) == ["b1", "b2"]
    pairs = pair_arcs(later, theirs, ref_b, params)
    assert pairs and all(c.startswith("b") for model, _ in pairs for c in model.arc_checksums)


# --------------------------------------------------------------------------- the blend


def test_blend_drifts_single_model_is_itself():
    drift = make_drift("x", "2026-10-01T20:00:00", [0.1], [0.0], {0: (0.1, 0.0, 10, 0.01, True)})
    blended, weights = _blend_drifts([(drift, 1.0)])
    assert blended is drift and weights is None


def test_blend_drifts_weights_two_models():
    fibers_a = {0: (0.2, 1e-5, 20, 0.02, True), 20: (0.4, 1e-5, 18, 0.03, True)}
    fibers_b = {0: (1.2, 3e-5, 15, 0.04, True), 44: (1.6, 3e-5, 16, 0.01, False)}
    first = make_drift("aaa", "2026-10-01T20:00:00", [0.2], [1e-5], fibers_a, scatter_px=0.02)
    second = make_drift("bbb", "2026-10-01T21:00:00", [2e-3, 1.2], [0.0, 3e-5], fibers_b, scatter_px=0.05)
    blended, weights = _blend_drifts([(first, 0.3), (second, 0.7)])
    for y in (10.0, 300.0):
        assert blended.shift_at(y) == pytest.approx(0.3 * first.shift_at(y) + 0.7 * second.shift_at(y))
        assert blended.stretch_at(y) == pytest.approx(0.3 * first.stretch_at(y) + 0.7 * second.stretch_at(y))
    # right-aligned coefficients: the constant model contributes to the last slot only
    assert blended.shift_coeffs == pytest.approx([0.7 * 2e-3, 0.3 * 0.2 + 0.7 * 1.2])
    assert blended.stretch_coeffs == pytest.approx([0.0, 0.3 * 1e-5 + 0.7 * 3e-5])
    assert blended.arc_checksums == ("aaa", "bbb")
    assert blended.arc_mid_time == "2026-10-01T20:42:00"  # 0.7 of the hour after 20:00
    assert blended.n_pix == N_PIX
    assert blended.scatter_px == pytest.approx(0.05)
    assert blended.usable is True
    # the fewest usable arc fibers of the blended models, not a count of the merged table
    assert (first.n_usable, second.n_usable, blended.n_usable) == (2, 1, 1)
    table = blended.per_fiber.set_index("fiber")
    assert sorted(table.index) == [0, 20, 44]
    assert table.loc[0, "shift_px"] == pytest.approx(0.3 * 0.2 + 0.7 * 1.2)
    assert table.loc[0, "stretch"] == pytest.approx(0.3 * 1e-5 + 0.7 * 3e-5)
    assert table.loc[0, "n_lines"] == 15 and table.loc[0, "usable"]
    assert table.loc[0, "rms_px"] == pytest.approx(0.3 * 0.02 + 0.7 * 0.04)
    # a fiber only one model holds takes that model's values and its weight
    assert table.loc[20, "shift_px"] == pytest.approx(0.4)
    assert table.loc[44, "shift_px"] == pytest.approx(1.6) and not table.loc[44, "usable"]
    assert weights == {0: pytest.approx(1.0), 20: pytest.approx(0.3), 44: pytest.approx(0.7)}
    # unnormalised weights are normalised; an unusable model makes the blend unusable
    same, _ = _blend_drifts([(first, 3.0), (second, 7.0)])
    assert same.shift_at(100.0) == pytest.approx(blended.shift_at(100.0))
    bad = make_drift("ccc", "2026-10-01T22:00:00", [0.5], [0.0], fibers_a, usable=False)
    assert not _blend_drifts([(first, 0.5), (bad, 0.5)])[0].usable
    assert _arc_weights([(first, 0.3), (second, 0.7)], ["aaa", "bbb"]) == pytest.approx([0.3, 0.7])


@pytest.mark.parametrize(
    ("usable", "scatter", "reason"),
    [
        ([True, True, False], 0.02, "2 usable arc fibers, min_fibers is 3"),
        ([True, True, True], 0.5, "scatter 0.500 px exceeds max_scatter_px 0.05"),
    ],
)
def test_unusable_drift_warning_counts_the_fibers_once(monkeypatch, caplog, usable, scatter, reason):
    """A drift model that is not usable is warned about once per arc: the gate's reason
    (with its own number, usable fibers or scatter) once, then how many arc fibers there
    were, and no number repeated."""
    from mastpipe.pipeline import arc_stage

    fibers = {fiber: (0.2, 1e-5, 20, 0.02, ok) for fiber, ok in zip((0, 20, 44), usable, strict=True)}
    drift = make_drift("aaa", "2026-10-01T20:00:00", [0.2], [1e-5], fibers, scatter_px=scatter, usable=False)
    arc_traces = [SimpleNamespace(index=fiber) for fiber in fibers]
    traceset = SimpleNamespace(by_role=lambda role: arc_traces, shape=(400, N_PIX))
    group = SimpleNamespace(
        effective=None,
        layout=None,
        out_dir=None,
        preset="survey45",
        traces=lambda: traceset,
        injected_masters=dict,
        drift_params=DriftParams(max_scatter_px=0.05),
        arcs=[],
        run_configured=frozenset(),
    )
    monkeypatch.setattr(arc_stage, "_reference_for", lambda group, at_time: make_reference("r1", "2026-10-01T18:00:00"))
    monkeypatch.setattr(arc_stage, "TraceSet", lambda **kwargs: None)
    monkeypatch.setattr(arc_stage, "_arc_lines", lambda *args, **kwargs: SimpleNamespace(tables=[None] * len(fibers)))
    monkeypatch.setattr(arc_stage, "fit_drift", lambda *args: drift)
    monkeypatch.setattr(arc_stage, "_write_qc", lambda *args: None)
    ref = FrameRef(path=Path("aaa.fits"), extension=0, meta=make_meta("2026-10-01T20:00:00", "aaa"))
    with caplog.at_level(logging.WARNING, logger="mastpipe"):
        arc_stage._cal_arc(group, ref)
    (message,) = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert message.startswith("drift arc aaa.fits: model not usable: ")
    assert message.count(reason) == 1 and message.endswith("; 3 arc fibers")
    assert message.count("usable arc fibers") == (1 if "min_fibers" in reason else 0)
    assert message.count("scatter ") == (1 if "max_scatter" in reason else 0)  # the scatter value


@pytest.mark.parametrize("configured", [True, False])
def test_full_arc_records_where_its_configuration_came_from(monkeypatch, configured):
    """A full arc given ``--configuration`` by the run carries the provenance step on the
    reference it is solved into, and ``configuration_origin`` ``run`` in its ``wref`` QC;
    one whose header named its grating carries no step and ``frame``."""
    from mastpipe.pipeline import arc_stage

    meta = make_meta("2026-10-01T20:00:00", "aaa", "arc")
    previous = make_reference("r1", "2026-10-01T18:00:00")
    added, qc = [], []
    library = SimpleNamespace(
        reference_for_source=lambda *args: None,
        find_reference=lambda *args: previous,
        add_reference=lambda solution, prov, camera: added.append(prov) or Path("new_wref.fits"),
    )
    group = SimpleNamespace(
        effective=SimpleNamespace(instrument="highspec", configuration="Mg", camera="main"),
        library=library,
        options=SimpleNamespace(force=False),
        layout=None,
        out_dir=None,
        preset="survey45",
        traces=lambda: None,
        injected_masters=dict,
        summary=RunSummary(),
        run_configured=frozenset({(meta.source, meta.extension)}) if configured else frozenset(),
    )
    lines = SimpleNamespace(traceset=[])
    solution = SimpleNamespace(identifier="r2", valid_from="2026-10-01T20:00:00", solutions={})
    monkeypatch.setattr(arc_stage, "_arc_lines", lambda *args, **kwargs: lines)
    monkeypatch.setattr(arc_stage, "_solve_and_package", lambda *args: solution)
    monkeypatch.setattr(arc_stage, "_provenance", lambda effective, lines, steps, layout: SimpleNamespace(steps=steps))
    monkeypatch.setattr(arc_stage, "_write_qc", lambda group, name, stage, values: qc.append((stage, values)))
    arc_stage._full_arc(group, FrameRef(path=Path("aaa.fits"), extension=0, meta=meta))
    (prov,) = added
    step = "configuration Mg from the run's --configuration (none in the header, file name or manifest)"
    assert (step in prov.steps) is configured
    assert prov.steps[-1] == "wavecal seeded by reference r1"
    ((stage, values),) = qc
    assert stage == "reference" and values["configuration_origin"] == ("run" if configured else "frame")


# ------------------------------------------------------------- no silent configuration


def test_group_without_configuration_fails_naming_the_gratings(sim_night, tmp_path):
    """A frame whose grating is neither in a GRATING card nor in its file name is not
    reduced with the instrument's only grating: the group fails and names the choices,
    the instrument's configurations (not every fleet grating name, some of which it
    cannot reduce)."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    del header["GRATING"]
    fits.PrimaryHDU(data, header).writeto(folder / "frame_0001.fits")
    summary = run_quicklook(folder, tmp_path / "out", RunOptions("highspec"))
    message = "configuration unknown; pass --configuration one of ['Mg']"
    assert summary.failed == [Failure("group", "main/(no configuration)", message, kind="configuration_unknown")]
    assert sorted(load_instrument("highspec").configurations) == ["Mg"]
    assert set(GRATING_NAMES) - {"Mg"}  # names the fleet knows but highspec cannot reduce
    assert summary.products == [] and summary.qc == []


def test_configuration_option_leaves_a_frame_naming_another_grating_alone(sim_night, tmp_path):
    """``--configuration`` names the grating of a frame whose header names none; a frame
    whose GRATING card names another grating keeps its own and is grouped under it (here a
    grating HighSpec does not define, so that group fails on its own)."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    header["GRATING"] = "Ca"
    fits.PrimaryHDU(data, header).writeto(folder / "frame_ca.fits")
    del header["GRATING"]
    fits.PrimaryHDU(data, header).writeto(folder / "frame_none.fits")
    summary = run_quicklook(folder, tmp_path / "out", RunOptions("highspec", configuration="Mg"))
    failed = {item.label: item for item in summary.failed}
    assert sorted(failed) == ["group main/Ca", "group main/Mg"]
    assert "'Ca'" in failed["group main/Ca"].message


def test_configuration_steps_name_the_run_option_only_for_frames_it_configured():
    """The provenance step and the QC origin of a frame given ``--configuration`` by the
    run (:func:`~mastpipe.pipeline.group._configuration_steps`), and of one whose header or
    manifest named its grating."""
    configured = make_meta("2026-10-01T20:00:00", "aaa", "science")
    named = make_meta("2026-10-01T20:00:00", "bbb", "science")
    group = SimpleNamespace(run_configured=frozenset({(configured.source, configured.extension)}))
    assert _configuration_steps(group, configured) == [
        "configuration Mg from the run's --configuration (none in the header, file name or manifest)"
    ]
    assert _configuration_origin(group, configured) == "run"
    assert _configuration_steps(group, named) == []
    assert _configuration_origin(group, named) == "frame"


# ------------------------------------------------------------------ assumed time zone


def test_one_warning_per_run_when_frame_times_carry_no_zone(sim_night, tmp_path):
    """The simulator writes a naive DATE-OBS, which HighSpec's map takes as UTC: the run log
    says so once, however many frames there are, and not at all when every stamp has its
    own offset."""
    sim_dir, _ = sim_night
    warning = "frame times carry no time zone; taken as UTC (header map 'assume-utc')"
    for name, offset in (("naive", ""), ("explicit", "+00:00")):
        folder = tmp_path / name
        folder.mkdir()
        for k, source in enumerate(("sim_bias_00.fits", "sim_bias_01.fits")):
            with fits.open(sim_dir / source) as hdul:
                data = hdul[0].data.copy()
                header = hdul[0].header.copy()
            del header["GRATING"]  # the group then fails fast, after the scan
            header["DATE-OBS"] = header["DATE-OBS"] + offset
            fits.PrimaryHDU(data, header).writeto(folder / f"frame_{k:04d}.fits")
        out = tmp_path / f"out_{name}"
        run_quicklook(folder, out, RunOptions("highspec"))
        (log_file,) = (out / "logs").glob("*.log")
        count = log_file.read_text(encoding="utf-8").count(warning)
        assert count == (1 if name == "naive" else 0), name


def test_one_warning_per_run_for_products_of_an_earlier_version(sim_night, tmp_path):
    """Product names gained a checksum segment, so a product written by an earlier version
    is never found at its new path: a new file is written beside it. The run log says so
    once, naming the ``spec`` folder, however many such files there are, and not at all
    when every product name carries the segment."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data, header = hdul[0].data.copy(), hdul[0].header.copy()
    del header["GRATING"]  # the group then fails fast, after the check
    fits.PrimaryHDU(data, header).writeto(folder / "frame_0001.fits")
    warning = "written by an earlier mastpipe version"
    for name, files in (
        (
            "old",
            [
                "highspec_Mg_main_20261001T180230_spec.fits",
                "highspec_Mg_main_notime_spec.fits",
                "highspec_Mg_main_20261001T180230_comb_target_a.fits",
            ],
        ),
        (
            "new",
            [
                "highspec_Mg_main_20261001T180230_a1b2c3_spec.fits",
                "highspec_Mg_main_notime_nosum_spec.fits",
                "highspec_Mg_main_20261001T180230_a1b2c3_comb_target_a.fits",
            ],
        ),
    ):
        out = tmp_path / f"out_{name}"
        (out / "spec").mkdir(parents=True)
        for file in files:
            (out / "spec" / file).write_bytes(b"")
        run_quicklook(folder, out, RunOptions("highspec"))
        (log_file,) = (out / "logs").glob("*.log")
        lines = [line for line in log_file.read_text(encoding="utf-8").splitlines() if warning in line]
        assert len(lines) == (1 if name == "old" else 0), name
        if lines:
            assert "WARNING" in lines[0] and str(out / "spec") in lines[0] and "3 product" in lines[0]


def test_one_warning_per_run_for_an_amplifier_without_its_noise(sim_night, tmp_path):
    """HighSpec lists gain and read noise for ``conventional`` and ``em``; frames naming
    another amplifier (here a misspelt one) take the camera's defaults, and the run log
    says so once, however many frames there are, and not at all for a listed amplifier."""
    sim_dir, _ = sim_night
    warning = "no gain and read noise listed for amplifier"
    for amplifier in ("Conventionl", "Conventional"):
        folder = tmp_path / amplifier
        folder.mkdir()
        for k, source in enumerate(("sim_bias_00.fits", "sim_bias_01.fits")):
            with fits.open(sim_dir / source) as hdul:
                data = hdul[0].data.copy()
                header = hdul[0].header.copy()
            del header["GRATING"]  # the group then fails fast, after the scan
            header["OUTPTAMP"] = amplifier
            fits.PrimaryHDU(data, header).writeto(folder / f"frame_{k:04d}.fits")
        out = tmp_path / f"out_{amplifier}"
        run_quicklook(folder, out, RunOptions("highspec"))
        (log_file,) = (out / "logs").glob("*.log")
        text = log_file.read_text(encoding="utf-8")
        assert text.count(warning) == (1 if amplifier == "Conventionl" else 0), amplifier
        if amplifier == "Conventionl":
            assert "conventionl" in text


# ---------------------------------------------------- other instruments, unmatched rows


def test_other_instruments_and_unmatched_manifest_rows_are_skipped(sim_night, tmp_path):
    """A nested night holding a DeepSpec frame next to HighSpec's: the DeepSpec frame and
    a manifest row naming no frame are listed in ``skipped`` with their reasons."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    del header["GRATING"]  # the group then fails fast, after the scan
    for name, instrument in (("highspec", "Highspec"), ("deepspec", "Deepspec")):
        (folder / name / "acquisition-0001").mkdir(parents=True)
        header["INSTRUME"] = instrument
        fits.PrimaryHDU(data, header).writeto(folder / name / "acquisition-0001" / "e-001.fits")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("source,target\nnothing.fits,x\n", encoding="utf-8")
    options = RunOptions("highspec", manifest=manifest)
    summary = run_quicklook(folder, tmp_path / "out", options)
    assert summary.skipped == [
        Skipped("frame", "deepspec/acquisition-0001/e-001.fits", "instrument 'deepspec' is not 'highspec'"),
        Skipped("manifest", "nothing.fits", "manifest row 0: matches no frame"),
    ]
    assert [item.label for item in summary.failed] == ["group main/(no configuration)"]


def test_manifest_row_naming_a_non_fits_file_has_its_own_reason(sim_night, tmp_path):
    """A row naming a file that exists but is not a FITS file name is skipped as such, not
    as matching no frame."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    shutil.copy(sim_dir / "sim_bias_00.fits", folder / "sim_bias_00.fits")
    (folder / "notes.txt").write_text("bias taken with the dome closed\n", encoding="utf-8")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("source,lamp_state\nnotes.txt,off\nnothing.fits,off\n", encoding="utf-8")
    summary = run_quicklook(folder, tmp_path / "out", RunOptions("highspec", manifest=manifest))
    assert Skipped("manifest", "notes.txt", "manifest row 0: not a FITS file (.fits, .fit, .fits.fz)") in summary.skipped
    assert Skipped("manifest", "nothing.fits", "manifest row 1: matches no frame") in summary.skipped
    assert not any(row.reason.startswith("manifest row 0") and "matches no frame" in row.reason for row in summary.skipped)


# ------------------------------------------------------------- structured summary items


def test_summary_items_are_structured(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    (empty / "junk.fits").write_bytes(b"this is not a FITS file\n" * 200)
    summary = run_quicklook(empty, tmp_path / "out", RunOptions("highspec"))
    assert summary.failed == [Failure("frame", "junk.fits", "unreadable FITS file")]
    assert summary.frames_scanned == 0
    assert summary.started is not None and summary.finished is not None
    assert summary.started.endswith("Z") and summary.started <= summary.finished


def test_empty_folder_has_no_failure(tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    summary = run_quicklook(empty, tmp_path / "out", RunOptions("highspec"))
    assert summary.failed == [] and summary.skipped == [] and summary.frames_scanned == 0


def test_item_labels():
    assert Skipped("group", "main/Mg", "not camera x").label == "group main/Mg"
    assert Failure("frame", "a/b.fits", "boom").label == "a/b.fits"
    assert Failure("run", None, "boom").label == "run"


# ------------------------------------------------------------ calibration frames, binning


def binned_ref(name: str, exposure_type: str, binning: tuple[int, int], **fields) -> FrameRef:
    meta = FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type=exposure_type,
        configuration="Mg",
        binning=binning,
        exposure_time=1.0,
        source=name,
        **fields,
    )
    return FrameRef(path=Path(name), extension=0, meta=meta)


def test_group_binning_follows_the_science_frames():
    """Many calibration frames at 1x1 do not outvote the science frames at 2x2 (light
    frames, which may be science, count with them); without any, all frames vote."""
    calibrations = [binned_ref(f"bias_{k}.fits", "bias", (1, 1)) for k in range(5)]
    calibrations += [binned_ref(f"flat_{k}.fits", "flat", (1, 1)) for k in range(3)]
    science = [
        binned_ref("sci_0.fits", "science", (2, 2)),
        binned_ref("sci_1.fits", "science", (1, 1)),
        binned_ref("light_0.fits", "light", (2, 2)),
    ]
    assert _group_binning(calibrations + science) == (2, 2)
    odd = [binned_ref("bias_odd.fits", "bias", (2, 1))]
    assert _group_binning(odd + calibrations) == (1, 1)


def test_group_readout_follows_the_science_frames():
    """The group's readout key and amplifier come from the frames that may be science, as
    its binning does; without any, all frames vote."""
    conventional = {"readout_mode": "image", "amplifier": "conventional"}
    em = {"readout_mode": "image", "amplifier": "em"}
    calibrations = [binned_ref(f"bias_{k}.fits", "bias", (1, 1), **conventional) for k in range(5)]
    science = [
        binned_ref("sci_0.fits", "science", (1, 1), **em),
        binned_ref("light_0.fits", "light", (1, 1), **em),
        binned_ref("sci_1.fits", "science", (1, 1), **conventional),
    ]
    assert _group_readout(calibrations + science) == ("image/em", "em")
    assert _group_readout(calibrations) == ("image/conventional", "conventional")
    assert _group_readout([binned_ref("x.fits", "bias", (1, 1))]) == ("unknown/unknown", None)


def usable_group(on_error: str = "continue") -> SimpleNamespace:
    """The parts of a ``_Group`` that :func:`_usable_frames` reads."""
    detector = SimpleNamespace(shape=(8, 10), crop=None, transpose=False, flip_x=False, flip_y=False)
    return SimpleNamespace(
        effective=SimpleNamespace(detector=detector),
        binning=(1, 1),
        readout_key="unknown/unknown",
        options=SimpleNamespace(on_error=on_error),
        summary=SimpleNamespace(failed=[], skipped=[]),
    )


def test_usable_frames_drops_bad_frames_with_their_reason(tmp_path):
    """A frame that cannot be read, has another binning, or another shape than the
    detector's (after crop and orientation) is recorded in ``failed`` and dropped; the
    others are kept in order. Under ``on_error='raise'`` the first bad frame raises."""
    refs = []
    for name, shape, binning in (
        ("a.fits", (8, 10), (1, 1)),
        ("binned.fits", (8, 10), (2, 1)),
        ("b.fits", (8, 10), (1, 1)),
        ("stub.fits", (4, 5), (1, 1)),
        ("c.fits", (8, 10), (1, 1)),
    ):
        fits.PrimaryHDU(np.full(shape, 600, dtype=np.uint16)).writeto(tmp_path / name)
        ref = binned_ref(name, "bias", binning)
        refs.append(FrameRef(path=tmp_path / name, extension=0, meta=ref.meta))
    (tmp_path / "junk.fits").write_bytes(b"this is not a FITS file\n" * 200)
    refs.insert(1, FrameRef(tmp_path / "junk.fits", 0, binned_ref("junk.fits", "bias", (1, 1)).meta))

    group = usable_group()
    frames = _usable_frames(group, refs, calibrate=False)
    assert [frame.meta.source for frame in frames] == ["a.fits", "b.fits", "c.fits"]
    failed = {item.source: item.message for item in group.summary.failed}
    assert list(failed) == ["junk.fits", "binned.fits", "stub.fits"]
    assert "binning 2x1" in failed["binned.fits"] and "1x1" in failed["binned.fits"]
    assert "shape (4, 5)" in failed["stub.fits"] and "(8, 10)" in failed["stub.fits"]

    # a stub sorting first does not set the shape: the good frames after it survive
    group = usable_group()
    frames = _usable_frames(group, [refs[4], refs[0], refs[3]], calibrate=False)
    assert [frame.meta.source for frame in frames] == ["a.fits", "b.fits"]
    assert [item.source for item in group.summary.failed] == ["stub.fits"]

    with pytest.raises(ValueError, match="binning 2x1"):
        _usable_frames(usable_group("raise"), refs[2:], calibrate=False)


def test_usable_frames_skips_a_frame_of_another_readout(tmp_path):
    """A calibration frame read through another amplifier than the group's is recorded in
    ``skipped`` with both readout keys and left out of the master: it is not a failure
    (it may belong to another group's night), and it does not fail the stack or raise
    under ``on_error='raise'``."""
    refs = []
    for name, amplifier in (
        ("a.fits", "conventional"),
        ("em.fits", "em"),
        ("b.fits", "conventional"),
    ):
        fits.PrimaryHDU(np.full((8, 10), 600, dtype=np.uint16)).writeto(tmp_path / name)
        meta = binned_ref(name, "bias", (1, 1), readout_mode="image", amplifier=amplifier).meta
        refs.append(FrameRef(path=tmp_path / name, extension=0, meta=meta))
    group = usable_group()
    group.readout_key = "image/conventional"
    frames = _usable_frames(group, refs, calibrate=False)
    assert [frame.meta.source for frame in frames] == ["a.fits", "b.fits"]
    assert group.summary.failed == []
    assert group.summary.skipped == [Skipped("frame", "em.fits", "readout image/em is not the group's image/conventional")]
    group = usable_group("raise")
    group.readout_key = "image/conventional"
    assert [frame.meta.source for frame in _usable_frames(group, refs, calibrate=False)] == ["a.fits", "b.fits"]
    assert len(group.summary.skipped) == 1


def calibrating_group(t_background: float | None) -> SimpleNamespace:
    """A group whose calibration set holds a bias and a background master, the latter
    with ``t_background`` as its exposure time, and background subtraction switched on."""
    group = usable_group()
    group.effective.detector.saturation = 65000.0
    group.params = SimpleNamespace(
        calib=SimpleNamespace(subtract_background=True, flat_subtract_background=True, cosmic_sigma=8.0, cosmic_size=5)
    )
    group.calib = CalibrationSet(
        bias=np.full((8, 10), 100.0, dtype=np.float32),
        dark_rate=0.0,
        background=np.full((8, 10), 130.0, dtype=np.float32),
        t_background=t_background,
        gain=1.0,
        read_noise=5.0,
        shape=(8, 10),
    )
    return group


def test_background_master_without_exposure_time_is_not_subtracted(tmp_path):
    """A background master whose inputs carried no exposure time cannot be scaled: the
    frames fall back to the bias (dark) path instead of every frame failing with
    'background subtraction requested ...'."""
    fits.PrimaryHDU(np.full((8, 10), 600, dtype=np.uint16)).writeto(tmp_path / "flat.fits")
    meta = binned_ref("flat.fits", "flat", (1, 1)).meta
    meta = FrameMeta(**{**meta.__dict__, "exposure_time": 2.0})
    ref = FrameRef(path=tmp_path / "flat.fits", extension=0, meta=meta)

    group = calibrating_group(t_background=None)
    frames = _usable_frames(group, [ref], calibrate=True)
    assert [frame.meta.source for frame in frames] == ["flat.fits"]
    assert group.summary.failed == []
    assert np.allclose(frames[0].data, 500.0)  # bias only

    frame, steps = _calibrated(group, ref)
    assert steps == ["bias"]
    assert np.allclose(frame.data, 500.0)

    # with an exposure time the background is scaled and subtracted as before
    group = calibrating_group(t_background=4.0)
    frame, steps = _calibrated(group, ref)
    assert steps[0] == "bias" and steps[1].startswith("background x 0.5")
    assert np.allclose(frame.data, 500.0 - 30.0 * 0.5)


# ------------------------------------------------------------------- masters and library

CALIB_OPTIONS = RunOptions("highspec", layout="survey45")


def calibration_night(sim_dir: Path, folder: Path, kinds: tuple[str, ...] = ("bias", "flat")) -> Path:
    """A night folder holding only the simulated night's frames of ``kinds``: the masters
    and the traces are reduced, nothing else."""
    folder.mkdir()
    for kind in kinds:
        for path in sorted(sim_dir.glob(f"sim_{kind}_*.fits")):
            shutil.copy(path, folder / path.name)
    return folder


def library_master(out: Path, kind: str) -> Path:
    (entry,) = CalibLibrary(out / "calib").list(kind=kind)
    return out / "calib" / entry["path"]


def test_library_flat_is_reused_only_with_the_same_masters_content(sim_night, tmp_path):
    """A library flat records the hash of the checksums of the masters it was calibrated
    with (``MPMKH``) and is reused only when the group's masters hash the same: a bias
    rewritten under the same name rebuilds it, and so does, once, a flat written before
    the hash was recorded. The trace set built on the flat carries the checksum of the
    flat's first input frame in its name, as the flat does, not ``nosum``."""
    folder = calibration_night(sim_night[0], tmp_path / "night")
    out = tmp_path / "out"
    first = run_quicklook(folder, out, CALIB_OPTIONS)
    assert first.failed == []
    bias_path, flat_path = library_master(out, "bias"), library_master(out, "flat")
    bias, _ = read_master(bias_path)
    flat, flat_prov = read_master(flat_path)
    assert flat_prov.masters_hash == _masters_hash({"bias": (bias, bias_path)})
    (traces,) = (out / "calib").glob("*_traces.fits")
    assert traces.name.split("_")[4] == flat_path.name.split("_")[4] == flat.provenance["inputs"][0][:6]

    # nothing changed: the flat is reused
    assert flat_path not in run_quicklook(folder, out, CALIB_OPTIONS).products

    # another bias under the same file name (and the same index entry): the bias is
    # reused as it stands, the flat calibrated with the old one is built again
    _, bias_prov = read_master(bias_path)
    write_master(bias_path, replace(bias, data=bias.data + 1.0), bias_prov)
    changed = run_quicklook(folder, out, CALIB_OPTIONS)
    assert changed.failed == []
    assert bias_path not in changed.products and flat_path in changed.products
    assert read_master(flat_path)[1].masters_hash == _masters_hash({"bias": (bias, bias_path)})
    assert flat_path not in run_quicklook(folder, out, CALIB_OPTIONS).products

    # a flat written before the masters hash was recorded is built again once
    with fits.open(flat_path, mode="update") as hdul:
        del hdul[0].header["MPMKH"]
    assert flat_path in run_quicklook(folder, out, CALIB_OPTIONS).products
    assert flat_path not in run_quicklook(folder, out, CALIB_OPTIONS).products
    # every run, rebuild or reuse, wrote the flat's trace set under the one name
    assert list((out / "calib").glob("*_traces.fits")) == [traces]


def test_flat_batch_is_chosen_among_the_frames_that_survive(sim_night, tmp_path, monkeypatch):
    """The flats are batched by exposure time and the batch with the most frames that load
    and calibrate is stacked, as for the other kinds: three 5 s flats of which two fail
    their calibration lose to two 8 s flats. A flat whose calibration raises is recorded
    in ``failed`` with the error."""
    sim_dir, _ = sim_night
    folder = calibration_night(sim_dir, tmp_path / "night")
    for k, source in enumerate(("sim_flat_00.fits", "sim_flat_01.fits")):
        with fits.open(sim_dir / source) as hdul:
            data, header = hdul[0].data.copy(), hdul[0].header.copy()
        header["EXPOSURE"] = 8.0
        header["DATE-OBS"] = f"2026-10-01T20:0{k}:00"
        fits.PrimaryHDU(data, header).writeto(folder / f"sim_flat_long_{k}.fits")
    real = masters.calibrate_frame

    def calibrate(frame, *args, **kwargs):
        if frame.meta.source in ("sim_flat_00.fits", "sim_flat_01.fits"):
            raise RuntimeError("calibration refused")
        return real(frame, *args, **kwargs)

    monkeypatch.setattr(masters, "calibrate_frame", calibrate)
    out = tmp_path / "out"
    summary = run_quicklook(folder, out, CALIB_OPTIONS)
    assert summary.failed == [
        Failure("frame", "sim_flat_00.fits", "calibration refused", error_type="RuntimeError"),
        Failure("frame", "sim_flat_01.fits", "calibration refused", error_type="RuntimeError"),
    ]
    flat, prov = read_master(library_master(out, "flat"))
    assert flat.validity["exposure_time"] == 8.0
    assert flat.n_frames == 2 and len(prov.inputs) == 2


def test_run_error_leaves_the_summary_and_the_log(tmp_path, capsys):
    night = tmp_path / "night"
    night.mkdir()
    out = tmp_path / "out"
    (out / "calib").mkdir(parents=True)
    (out / "calib" / "index.json").write_text("{not json", encoding="utf-8")
    summary = RunSummary()
    with pytest.raises(LibraryInvalidError):
        run_quicklook(night, out, RunOptions("highspec"), summary)
    assert summary.frames_scanned == 0 and summary.finished is not None
    (log_file,) = (out / "logs").glob("quicklook_*.log")
    assert "run failed" in log_file.read_text(encoding="utf-8")
    assert "run failed" not in capsys.readouterr().err  # stderr is the command's


@pytest.mark.parametrize(
    ("options", "error"),
    [
        (RunOptions("highspec"), FileNotFoundError),
        (RunOptions("highspec", on_error="ignore"), ValueError),
        (RunOptions("highspec", time_policy="sometime"), ValueError),
    ],
)
def test_early_failure_is_stamped_and_writes_nothing(tmp_path, options, error):
    """A run stopped by its precondition checks still has ``started`` and ``finished``;
    it creates no ``--out`` and no log."""
    summary = RunSummary()
    with pytest.raises(error):
        run_quicklook(tmp_path / "nowhere", tmp_path / "out", options, summary)
    assert summary.started is not None and summary.finished is not None
    assert summary.started <= summary.finished
    assert summary.log is None and summary.library is None and summary.frames_scanned is None
    assert not (tmp_path / "out").exists()


def test_out_that_is_a_file_is_output_unwritable(tmp_path):
    night = tmp_path / "night"
    night.mkdir()
    (tmp_path / "out").write_text("a file", encoding="utf-8")
    with pytest.raises(OutputUnwritableError):
        run_quicklook(night, tmp_path / "out", RunOptions("highspec"))


def test_two_runs_in_one_second_have_two_logs(tmp_path, monkeypatch):
    night = tmp_path / "night"
    night.mkdir()
    fixed = datetime(2026, 1, 2, 3, 4, 5, tzinfo=UTC)

    class Clock(datetime):
        @classmethod
        def now(cls, tz=None):
            return fixed

    monkeypatch.setattr("mastpipe.pipeline.runner.datetime", Clock)
    first = run_quicklook(night, tmp_path / "out", RunOptions("highspec"))
    second = run_quicklook(night, tmp_path / "out", RunOptions("highspec"))
    assert first.log != second.log and first.log.is_file() and second.log.is_file()
    assert re.fullmatch(r"quicklook_20260102T030405_[0-9a-f]{6}\.log", first.log.name)
    assert first.log.parent == tmp_path / "out" / "logs"


def test_a_drawn_suffix_already_taken_draws_another(tmp_path, monkeypatch):
    night = tmp_path / "night"
    night.mkdir()
    fixed = datetime(2026, 1, 2, 3, 4, 5, tzinfo=UTC)

    class Clock(datetime):
        @classmethod
        def now(cls, tz=None):
            return fixed

    suffixes = iter(["aaaaaa", "aaaaaa", "bbbbbb"])
    monkeypatch.setattr("mastpipe.pipeline.runner.datetime", Clock)
    monkeypatch.setattr("mastpipe.pipeline.runner.secrets.token_hex", lambda n: next(suffixes))
    first = run_quicklook(night, tmp_path / "out", RunOptions("highspec"))
    second = run_quicklook(night, tmp_path / "out", RunOptions("highspec"))
    assert first.log.name == "quicklook_20260102T030405_aaaaaa.log"
    assert second.log.name == "quicklook_20260102T030405_bbbbbb.log"
    assert "quicklook of" in first.log.read_text(encoding="utf-8")


def test_quiet_keeps_progress_off_stderr(tmp_path, capsys):
    night = tmp_path / "night"
    night.mkdir()
    run_quicklook(night, tmp_path / "out", RunOptions("highspec", quiet=True))
    err = capsys.readouterr().err
    assert "INFO" not in err and "WARNING no frames of highspec" in err
    run_quicklook(night, tmp_path / "out", RunOptions("highspec"))
    assert "INFO quicklook of" in capsys.readouterr().err


def test_quiet_run_log_keeps_progress(tmp_path, monkeypatch):
    night = tmp_path / "night"
    night.mkdir()
    monkeypatch.chdir(tmp_path)
    summary = run_quicklook("night", "out", RunOptions("highspec", quiet=True))
    assert summary.log.parent == Path("out", "logs")  # under --out as given
    assert "INFO quicklook of" in summary.log.read_text(encoding="utf-8")


def test_summary_names_the_log_of_a_run_that_raises(tmp_path):
    night = tmp_path / "night"
    night.mkdir()
    out = tmp_path / "out"
    (out / "calib").mkdir(parents=True)
    (out / "calib" / "index.json").write_text("{not json", encoding="utf-8")
    summary = RunSummary()
    with pytest.raises(LibraryInvalidError):
        run_quicklook(night, out, RunOptions("highspec"), summary)
    assert summary.log is not None and summary.log.parent == out / "logs"
    assert "run failed" in summary.log.read_text(encoding="utf-8")


def test_busy_library_stays_library_busy(tmp_path, monkeypatch):
    def busy(path, timeout_s=30.0):
        raise TimeoutError(f"could not lock {path}")

    night = tmp_path / "night"
    night.mkdir()
    monkeypatch.setattr("mastpipe.pipeline.callib.locked", busy)
    with pytest.raises(LibraryBusyError, match="could not lock"):
        run_quicklook(night, tmp_path / "out", RunOptions("highspec"))


def test_unexpected_error_is_logged_with_its_traceback_in_the_log_only(tmp_path, monkeypatch, capsys):
    def defect(*args, **kwargs):
        raise RuntimeError("a defect")

    night = tmp_path / "night"
    night.mkdir()
    monkeypatch.setattr("mastpipe.pipeline.runner.scan_folder", defect)
    summary = RunSummary()
    with pytest.raises(RuntimeError, match="a defect"):
        run_quicklook(night, tmp_path / "out", RunOptions("highspec"), summary)
    text = summary.log.read_text(encoding="utf-8")
    assert "run failed" in text and "Traceback" in text and "RuntimeError: a defect" in text
    err = capsys.readouterr().err
    assert "Traceback" not in err and "a defect" not in err


def test_interrupt_is_one_plain_line_in_the_log(tmp_path, monkeypatch, capsys):
    def interrupt(*args, **kwargs):
        raise KeyboardInterrupt

    night = tmp_path / "night"
    night.mkdir()
    monkeypatch.setattr("mastpipe.pipeline.runner.scan_folder", interrupt)
    summary = RunSummary()
    with pytest.raises(KeyboardInterrupt):
        run_quicklook(night, tmp_path / "out", RunOptions("highspec"), summary)
    text = summary.log.read_text(encoding="utf-8")
    assert "run interrupted" in text and "Traceback" not in text and "run failed" not in text
    assert "run interrupted" not in capsys.readouterr().err
