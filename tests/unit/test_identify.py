"""Tests for the ``identify`` and ``solve`` commands on the simulated night: the anchor
panels and candidate table without anchors, the reference product written from truth
anchors, the seeded re-solve, the anchor-pair parsing and the command line."""

import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml

from mastpipe.cli import main
from mastpipe.config.loader import build_effective
from mastpipe.config.models import CalibParams
from mastpipe.core.model import Frame, ReferenceSolution
from mastpipe.instruments.registry import load_header_map, load_instrument
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.io.products import product_name, read_reference, reference_name, write_master
from mastpipe.io.qa import plot_anchor_panels, strongest_spread
from mastpipe.pipeline import identify
from mastpipe.pipeline.identify import (
    _arc_lines,
    _latest_master,
    parse_anchor_pairs,
    read_anchors_file,
    run_identify,
    run_solve,
)
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.atlas import catalogue_arrays
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.lines import find_lines
from mastpipe.steps.wavecal import wavelength_map
from tests.conftest import truth_model

MIN_SOLVED = 43
MAX_MAP_ERROR_DISP = 0.3
ANCHOR_LOG_INTENSITY = 2.0


# ------------------------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def survey_effective():
    """The HighSpec configuration with the ``survey45`` preset's parameters (the 8.5 px
    pitch narrows the three halfwidths to 4 px)."""
    return build_effective(load_instrument("highspec"), "Mg", "main", layout="survey45")


@pytest.fixture(scope="module")
def arc(sim_night, survey_effective, mg_lines) -> SimpleNamespace:
    """The simulated full arc on disk, its truth anchors and the grating truth."""
    night_dir, truth = sim_night
    row = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    path = night_dir / str(row["file"])
    geom = grating_geometry(survey_effective.geometry, survey_effective.detector.pixel_size_um)
    n_pix = truth.extras["n_pix"]
    cat_wl, _ = catalogue_arrays(mg_lines, usable_only=True)
    usable = np.isin(truth.lines_nm, cat_wl)
    bright = usable & truth.lines_on_chip & (truth.line_logI > ANCHOR_LOG_INTENSITY)
    first, last = np.flatnonzero(bright)[[0, -1]]
    anchors = [
        (float(truth.lines_nm[i]), round(float(truth.line_pixels[0][i]), 1)) for i in (first, last)
    ]
    return SimpleNamespace(
        path=path,
        mid_time=str(row["mid_time"]),
        truth=truth,
        anchors=anchors,
        geom=geom,
        n_pix=n_pix,
        truth_model=lambda fiber: truth_model(truth, fiber, geom, n_pix),
        cat_wl=cat_wl,
    )


@pytest.fixture(scope="module")
def calib_dir(sim_night, survey_effective, tmp_path_factory) -> Path:
    """An output folder whose ``calib/`` holds bias, background and flat masters."""
    night_dir, _ = sim_night
    detector = survey_effective.detector
    out_dir = tmp_path_factory.mktemp("identify_calib")
    refs = scan_folder(night_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        if ref.meta.exposure_type in ("bias", "background", "flat"):
            by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    params = CalibParams()
    bias = build_master(by_type["bias"], "bias", params)
    background = build_master(by_type["background"], "background", params)
    calib = build_calibration_set(bias, None, None, detector)
    flats = [calibrate_frame(f, calib, subtract_background=False) for f in by_type["flat"]]
    flat = build_master(flats, "flat", params)
    prov = Provenance("0.0.1", "abc", [], [])
    for master in (bias, background, flat):
        name = product_name(by_type[master.kind][0].meta, master.kind)
        write_master(out_dir / "calib" / name, master, prov)
    return out_dir


def master_names(calib_dir: Path) -> list[str]:
    return sorted(p.name for p in (calib_dir / "calib").glob("*.fits") if "_wref" not in p.name)


@pytest.fixture(scope="module")
def identified(arc, survey_effective, calib_dir) -> SimpleNamespace:
    ref, path = run_identify(
        arc.path, survey_effective, "survey45", calib_dir, anchors=arc.anchors, anchor_trace=0
    )
    return SimpleNamespace(ref=ref, path=path, out_dir=calib_dir)


def map_error_in_dispersions(solution, arc) -> float:
    """Largest map error against the fiber's own truth model, in dispersions."""
    lam = wavelength_map(solution, arc.geom)
    x = np.arange(arc.n_pix)
    model = arc.truth_model(solution.fiber)
    return float(np.max(np.abs(lam - model.wavelength(x)) / model.dispersion(x)))


# ------------------------------------------------------------------------- anchor panels


def test_strongest_spread_picks_strong_separated_lines_left_to_right():
    positions = np.array([10.0, 12.0, 50.0, 90.0, 95.0, 200.0])
    strength = np.array([5.0, 9.0, 1.0, 7.0, 8.0, 3.0])
    picked = strongest_spread(positions, strength, 3, 0.0, 100.0, 0.06)
    # 12 beats 10 (2 px apart, within 6 % of 100), 95 beats 90, 200 is outside the window,
    # 50 is the third, numbered left to right; a gap equal to the limit is not a conflict
    assert picked == [1, 2, 4]
    assert strongest_spread(positions, strength, 3, 0.0, 100.0, 0.05) == [1, 3, 4]
    assert strongest_spread(positions, strength, 10, 0.0, 300.0, 0.0) == [0, 1, 2, 3, 4, 5]
    assert strongest_spread(positions, strength, 1, 0.0, 300.0, 0.0) == [1]
    assert strongest_spread(np.zeros(0), np.zeros(0), 3, 0.0, 1.0, 0.1) == []


def test_plot_anchor_panels_writes_png_and_numbers_both_panels(tmp_path):
    wave = np.linspace(513.0, 527.5, 4000)
    centres = np.array([514.0, 516.5, 519.0, 521.5, 524.0, 526.5])
    flux = 1.0 + sum(
        h * np.exp(-0.5 * ((wave - c) / 0.02) ** 2)
        for c, h in zip(centres, [5, 9, 2, 7, 8, 4], strict=True)
    )
    usable = pd.DataFrame(
        {"wavelength": np.concatenate([centres, [520.0]]), "usable": [True] * 6 + [False]}
    )
    obs_flux = np.zeros(1000)
    x_lines = np.array([100.0, 300.0, 500.0, 700.0, 900.0])
    for x in x_lines:
        obs_flux += 100.0 * np.exp(-0.5 * ((np.arange(1000) - x) / 1.3) ** 2)
    obs_flux[:5] = np.nan  # masked columns must not break the normalisation
    obs_lines = pd.DataFrame(
        {
            "x_center": np.concatenate([x_lines, [950.0]]),
            "height": [100.0, 110.0, 90.0, 120.0, 80.0, 500.0],
            "use_x": [True] * 5 + [False],
        }
    )
    out = tmp_path / "panels" / "anchor_panels.png"
    template, measured = plot_anchor_panels(
        wave, flux, usable, obs_flux, obs_lines, (513.0, 527.5), n_labels=4, out=out
    )
    assert out.is_file() and out.stat().st_size > 0
    assert list(template.columns) == ["id", "wavelength"]
    assert list(measured.columns) == ["id", "x_center"]
    assert template["id"].tolist() == [1, 2, 3, 4]
    assert measured["id"].tolist() == [1, 2, 3, 4]
    # the strongest four usable lines, numbered left to right; 520.0 is not usable
    assert template["wavelength"].tolist() == [514.0, 516.5, 521.5, 524.0]
    assert measured["x_center"].tolist() == [100.0, 300.0, 500.0, 700.0]
    assert template["wavelength"].is_monotonic_increasing


# ---------------------------------------------------------------------------- identify


def test_identify_without_anchors_writes_panels_and_candidates(arc, survey_effective, tmp_path):
    out_dir = tmp_path / "out"
    ref, png = run_identify(arc.path, survey_effective, "survey45", out_dir, anchor_trace=0)
    assert ref is None
    assert png == out_dir / "anchor_panels.png"
    assert png.is_file() and png.stat().st_size > 0
    table = pd.read_csv(out_dir / "anchor_candidates.csv", encoding="utf-8")
    assert list(table.columns) == ["id", "wavelength_nm", "x_center"]
    template = table.loc[table["id"].str.startswith("template_")]
    measured = table.loc[table["id"].str.startswith("measured_")]
    assert len(template) == 8 and len(measured) == 8
    assert template["id"].tolist() == [f"template_{k}" for k in range(1, 9)]
    assert measured["id"].tolist() == [f"measured_{k}" for k in range(1, 9)]
    # template candidates are usable catalogue lines; measured ones sit on truth line pixels
    assert set(template["wavelength_nm"]) <= set(arc.cat_wl)
    assert template["x_center"].isna().all() and measured["wavelength_nm"].isna().all()
    truth_px = arc.truth.line_pixels[0][arc.truth.lines_on_chip]
    for x in measured["x_center"]:
        assert np.min(np.abs(truth_px - x)) < 0.3
    assert not (out_dir / "calib").exists()


def test_identify_with_anchors_writes_reference(identified, arc, survey_effective):
    ref, path = identified.ref, identified.path
    assert isinstance(ref, ReferenceSolution)
    assert path == identified.out_dir / "calib" / reference_name(ref)
    assert path.is_file()
    back, prov = read_reference(path, lambda instrument, configuration: arc.geom)
    assert back.identifier == ref.identifier and back.hash == ref.hash
    assert back.instrument == "highspec" and back.configuration == "Mg"
    assert back.fiber_preset == "survey45"
    assert back.valid_from == back.created == arc.mid_time
    assert len(back.solutions) >= MIN_SOLVED
    for solution in back.solutions.values():
        assert map_error_in_dispersions(solution, arc) < MAX_MAP_ERROR_DISP
    # provenance: the arc's checksum, the masters used, solved status
    assert back.source_checksums == tuple(prov.inputs) and len(prov.inputs) == 1
    assert prov.wavelength_status == "solved"
    assert sorted(prov.calibrations) == master_names(identified.out_dir)
    assert len(prov.calibrations) == 3
    assert prov.pipeline_version and prov.config_hash
    agreement = pd.read_csv(identified.out_dir / "calib" / "fiber_agreement.csv", encoding="utf-8")
    assert list(agreement.columns) == ["wavelength", "n_fibers", "spread_px"]
    assert agreement["n_fibers"].max() == len(back.solutions)
    assert agreement["spread_px"].max() < MAX_MAP_ERROR_DISP


def test_identify_uses_stored_anchors_and_trace(
    tmp_path, sim_night_lab3, highspec_effective, mg_lines, monkeypatch
):
    """Without an ``anchors`` argument the instrument's stored pairs and trace are used."""
    night_dir, truth = sim_night_lab3
    row = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    cat_wl, _ = catalogue_arrays(mg_lines, usable_only=True)
    bright = (
        np.isin(truth.lines_nm, cat_wl)
        & truth.lines_on_chip
        & (truth.line_logI > ANCHOR_LOG_INTENSITY)
    )
    first, last = np.flatnonzero(bright)[[0, -1]]
    stored = [(float(truth.lines_nm[i]), float(truth.line_pixels[1][i])) for i in (first, last)]
    asked: list[tuple] = []
    seeded: dict = {}
    original = identify.solve_full_arc

    def spy(*args, **kwargs):
        seeded.update(kwargs)
        return original(*args, **kwargs)

    monkeypatch.setattr(identify, "load_anchors", lambda *a: asked.append(a) or stored)
    monkeypatch.setattr(identify, "load_anchor_trace", lambda *a: 1)
    monkeypatch.setattr(identify, "solve_full_arc", spy)
    out_dir = tmp_path / "lab"
    ref, path = run_identify(night_dir / str(row["file"]), highspec_effective, "lab3", out_dir)
    assert asked == [("highspec", "Mg", "lab3")]
    assert seeded["ref_fiber"] == 1
    assert isinstance(ref, ReferenceSolution) and set(ref.solutions) == {0, 1, 2}
    assert ref.fiber_preset == "lab3" and path.is_file()
    assert not (out_dir / "anchor_panels.png").exists()
    # the output folder holds no masters: traces came from the arc, nothing was subtracted
    _, prov = read_reference(path, lambda instrument, configuration: None)
    assert prov.calibrations == []
    assert any(step.startswith("traces arc_full") for step in prov.steps)
    # write_reference=False solves but writes nothing; the path is where it would go
    monkeypatch.setattr(identify, "load_anchor_trace", lambda *a: 0)
    dry, dry_path = run_identify(
        night_dir / str(row["file"]),
        highspec_effective,
        "lab3",
        tmp_path / "dry",
        anchors=stored,
        write_reference=False,
    )
    assert isinstance(dry, ReferenceSolution) and set(dry.solutions) == {0, 1, 2}
    for fiber, solution in ref.solutions.items():  # seeded from trace 0 instead of 1
        assert dry.solutions[fiber].params[0] == pytest.approx(solution.params[0], abs=1e-6)
        assert dry.solutions[fiber].params[1] == pytest.approx(solution.params[1], rel=1e-6)
    assert dry_path == tmp_path / "dry" / "calib" / reference_name(dry)
    assert not (tmp_path / "dry").exists()


def test_latest_master_sorts_by_stamp_within_the_configuration(tmp_path, survey_effective):
    names = [
        "highspec_Mg_main_20260930T235900_bias.fits",
        "highspec_Mg_main_20261001T180000_bias.fits",
        "highspec_Mg_main_notime_bias.fits",  # sorts after every digit stamp by name
        "highspec_noconf_main_20270101T000000_bias.fits",  # another configuration, later
        "highspec_Mg_other_20270101T000000_bias.fits",  # another camera, later
        "highspec_Mg_main_20270101T000000_flat.fits",  # another kind
    ]
    for name in names:
        (tmp_path / name).write_bytes(b"")
    assert sorted(names)[-1] != "highspec_Mg_main_20261001T180000_bias.fits"
    chosen = _latest_master(tmp_path, survey_effective, "bias")
    assert chosen is not None and chosen.name == "highspec_Mg_main_20261001T180000_bias.fits"
    assert _latest_master(tmp_path, survey_effective, "background") is None
    for name in names[:2]:
        (tmp_path / name).unlink()
    chosen = _latest_master(tmp_path, survey_effective, "bias")
    assert chosen is not None and chosen.name == "highspec_Mg_main_notime_bias.fits"


def test_latest_master_reads_the_stamp_before_the_checksum(tmp_path, survey_effective):
    names = [
        "highspec_Mg_main_20260930T235900_ffffff_bias.fits",
        "highspec_Mg_main_20261001T180000_000000_bias.fits",  # latest, lowest checksum
        "highspec_Mg_main_notime_nosum_bias.fits",
        "highspec_Mg_main_20270101T000000_abcdef_flat.fits",  # another kind
    ]
    for name in names:
        (tmp_path / name).write_bytes(b"")
    chosen = _latest_master(tmp_path, survey_effective, "bias")
    assert chosen is not None and chosen.name == "highspec_Mg_main_20261001T180000_000000_bias.fits"
    for name in names[:2]:
        (tmp_path / name).unlink()
    chosen = _latest_master(tmp_path, survey_effective, "bias")
    assert chosen is not None and chosen.name == "highspec_Mg_main_notime_nosum_bias.fits"


def test_arc_lines_lists_only_applied_masters_and_takes_injected_traces(
    sim_night_lab3, highspec_effective, lab3_layout, tmp_path
):
    night_dir, truth = sim_night_lab3
    detector = highspec_effective.detector
    refs = scan_folder(night_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        if ref.meta.exposure_type in ("bias", "background"):
            by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    bias = build_master(by_type["bias"], "bias", CalibParams())
    background = build_master(by_type["background"], "background", CalibParams())
    background.provenance["name"] = "night background"
    row = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    arc = night_dir / str(row["file"])
    masters = {"bias": bias, "background": background}

    lines = _arc_lines(arc, highspec_effective, lab3_layout, tmp_path, masters=masters)
    assert lines.calibrations == [f"bias:{bias.validity['mid_time']}", "night background"]
    assert lines.steps[:2] == ["bias", "background x 0.1 (t 30 s / t_background 300 s)"]
    assert len(lines.traceset) == 3 and lines.steps[2].startswith("traces arc_full (3 found")
    assert abs(float(np.median(lines.frame.data)) - 0.0) < 2.0  # bias and background gone

    # with background subtraction off the background is not an applied calibration
    no_subtract = highspec_effective.model_copy(
        update={
            "params": highspec_effective.params.model_copy(
                update={"calib": CalibParams(subtract_background=False)}
            )
        }
    )
    lines2 = _arc_lines(
        arc, no_subtract, lab3_layout, tmp_path, masters=masters, traceset=lines.traceset
    )
    assert lines2.calibrations == [f"bias:{bias.validity['mid_time']}"]
    assert lines2.steps[0] == "bias"
    assert lines2.steps[1].startswith("background night background present but not subtracted")
    assert lines2.steps[2] == "traces given (3 found, 0 missing)"
    assert lines2.traceset is lines.traceset
    assert not (tmp_path / "calib").exists()


def test_arc_lines_passes_the_calibrated_mask_to_find_lines(
    sim_night_lab3, highspec_effective, lab3_layout, tmp_path, monkeypatch
):
    """The arc's saturation mask reaches the line finder, which flags lines through it."""
    night_dir, truth = sim_night_lab3
    row = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    seen = []

    def spy(*args, **kwargs):
        seen.append(kwargs.get("mask"))
        return find_lines(*args, **kwargs)

    monkeypatch.setattr(identify, "find_lines", spy)
    lines = _arc_lines(night_dir / str(row["file"]), highspec_effective, lab3_layout, tmp_path)
    assert len(seen) == 1 and seen[0] is lines.frame.mask and seen[0] is not None


# ------------------------------------------------------------------------------- solve


def test_run_solve_seeded_by_previous_reference(identified, arc, survey_effective, tmp_path):
    previous = identified.ref
    out_dir = tmp_path / "again"
    # the same masters as the anchor run, so the two solves see the same lines
    (out_dir / "calib").mkdir(parents=True)
    for name in master_names(identified.out_dir):
        shutil.copy(identified.out_dir / "calib" / name, out_dir / "calib" / name)
    again = run_solve(arc.path, survey_effective, "survey45", out_dir, previous)
    path = out_dir / "calib" / reference_name(again)
    assert path.is_file()
    assert set(again.solutions) == set(previous.solutions)
    for fiber, solution in previous.solutions.items():
        other = again.solutions[fiber]
        assert other.lines["wavelength"].tolist() == solution.lines["wavelength"].tolist()
        assert other.params[0] == pytest.approx(solution.params[0], abs=1e-6)
        assert other.params[1] == pytest.approx(solution.params[1], rel=1e-6)
    back, prov = read_reference(path, lambda instrument, configuration: arc.geom)
    assert back.hash == again.hash
    assert prov.reference_id == previous.identifier and prov.reference_hash == previous.hash
    # a reference of another configuration or preset cannot seed this arc
    with pytest.raises(ValueError, match="preset"):
        run_solve(arc.path, survey_effective, "lab3", out_dir, previous)


# --------------------------------------------------------------------------------- CLI


def test_cli_identify_parses_anchor_pairs(tmp_path, sim_night_lab3, capsys):
    pairs = parse_anchor_pairs("515.8605:267.22,525.8360:1423.22")
    assert pairs == [(515.8605, 267.22), (525.836, 1423.22)]
    assert parse_anchor_pairs(" 515.8605:267.22 , 525.8360:1423.22 ") == pairs
    for bad in ("515.8605", "a:b", "515.8605:267.22:3", ""):
        with pytest.raises(ValueError):
            parse_anchor_pairs(bad)
    anchors_file = tmp_path / "anchors.yaml"
    anchors_file.write_text(
        yaml.safe_dump({"trace": 2, "pairs": [[515.8605, 267.22], [525.836, 1423.22]]}),
        encoding="utf-8",
    )
    assert read_anchors_file(anchors_file) == (pairs, 2)
    anchors_file.write_text(yaml.safe_dump({"pairs": [[515.8605, 267.22]]}), encoding="utf-8")
    assert read_anchors_file(anchors_file) == ([(515.8605, 267.22)], None)

    # the no-anchors path end to end on the simulated 3-fiber laboratory night
    night_dir, truth = sim_night_lab3
    row = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    out_dir = tmp_path / "cli"
    argv = [
        "identify",
        str(night_dir / str(row["file"])),
        "--out",
        str(out_dir),
        "--instrument",
        "highspec",
        "--configuration",
        "Mg",
        "--layout",
        "lab3",
        "--set",
        "trace.degree=2",
    ]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert str(out_dir / "anchor_panels.png") in out
    assert (out_dir / "anchor_panels.png").is_file()
    assert (out_dir / "anchor_candidates.csv").is_file()
    # the three anchor sources exclude one another
    with pytest.raises(SystemExit):
        main([*argv, "--anchors", "515.8605:267.22", "--stored-anchors"])
