"""Tests for the instrument registry and the shipped HighSpec definition."""

import shutil

import pytest

from mastpipe.core.model import ROLES
from mastpipe.instruments import registry
from mastpipe.instruments.registry import (
    LAYOUT_COLUMNS,
    LINE_TABLE_COLUMNS,
    instrument_dir,
    list_instruments,
    load_anchor_trace,
    load_anchors,
    load_header_map,
    load_instrument,
    load_layout,
    load_line_table,
)


def test_highspec_is_registered():
    assert "highspec" in list_instruments()
    assert (instrument_dir("highspec") / "config.yaml").is_file()


def test_highspec_config_loads_with_grating_blocks():
    inst = load_instrument("highspec")
    assert inst.name == "highspec"
    assert "main" in inst.cameras
    assert inst.cameras["main"].shape == (400, 1600)
    assert inst.cameras["main"].crop == (1, 1599, 1, 399)
    assert "mgb518" in inst.configurations
    block = inst.configurations["mgb518"]
    assert block.geometry.model == "grating"
    assert block.geometry.central_wavelength_nm == pytest.approx(518.0)
    assert block.band_nm == (513.0, 527.5)
    assert block.line_table == "lines/mgb518.csv"
    assert inst.default_layout == "lab3"
    assert set(inst.layouts) >= {"lab3", "survey45"}
    assert inst.params.lines.sigma_guess == pytest.approx(3.0 / 2.3548)


def test_layout_tables_parse_and_validate():
    lab3 = load_layout("highspec", "lab3")
    assert list(lab3.columns) == list(LAYOUT_COLUMNS)
    assert len(lab3) == 3 and set(lab3.role) <= set(ROLES)
    assert lab3.trace.dtype.kind == "i"
    assert lab3.telescope.isna().all()
    assert lab3.y_nominal.tolist() == pytest.approx([17.4, 189.0, 360.7])
    survey = load_layout("highspec", "survey45")
    assert len(survey) == 45
    assert (survey.role == "arc").sum() == 5
    assert (survey.role == "science").sum() == 20 and (survey.role == "sky").sum() == 20
    assert survey.y_nominal.is_monotonic_increasing
    assert survey.loc[survey.role != "arc", "telescope"].notna().all()
    assert survey.loc[survey.role == "arc", "trace"].tolist() == [0, 11, 22, 33, 44]
    assert sorted(set(survey.loc[survey.role != "arc", "telescope"])) == [
        f"T{i:02d}" for i in range(1, 21)
    ]


def test_unknown_layout_preset():
    with pytest.raises(KeyError):
        load_layout("highspec", "nope")


def test_line_table_loads():
    table = load_line_table("highspec", "mgb518")
    assert list(table.columns) == list(LINE_TABLE_COLUMNS)
    assert len(table) == 73
    assert table.wavelength_nm.is_monotonic_increasing
    assert table.wavelength_nm.between(513.0, 527.5).all()
    assert table.usable.dtype == bool and table.usable.sum() >= 25
    assert "Th I" in set(table.species)
    assert all(isinstance(value, str) for value in table.species)


def test_anchors_load():
    pairs = load_anchors("highspec", "mgb518", "lab3")
    assert len(pairs) == 2
    assert pairs[0][0] == pytest.approx(515.8605)
    assert pairs[0][1] == pytest.approx(267.22)
    assert pairs[1] == pytest.approx((525.8360, 1423.22))
    assert load_anchors("highspec", "mgb518", "survey45") == []


def test_header_map_yaml_has_required_sections():
    hm = load_header_map("highspec")
    assert "fields" in hm and "classification" in hm
    assert "exposure_type" in hm["fields"] and "exposure_time" in hm["fields"]
    assert hm["fields"]["exposure_time"]["keys"] == ["EXPOSURE", "EXPTIME"]
    assert hm["fields"]["exposure_type"]["map"]["zero"] == "bias"
    assert isinstance(hm["classification"], list) and len(hm["classification"]) >= 5


def test_unknown_instrument():
    with pytest.raises(KeyError):
        load_instrument("nope")
    with pytest.raises(KeyError):
        instrument_dir("nope")


def test_extra_grating_file_is_picked_up(tmp_path, monkeypatch):
    """A new HighSpec grating is one YAML file; no code change may be needed."""
    copy = tmp_path / "highspec"
    shutil.copytree(registry.instrument_dir("highspec"), copy)
    (copy / "gratings" / "halpha.yaml").write_text(
        "name: halpha\nband_nm: [649.0, 664.0]\nresolving_power: 20000\n"
        "geometry: {model: grating, groove_spacing_nm: 367.0, incidence_angle_deg: 45.0, "
        "order: 1, camera_focal_mm: 479.41, central_wavelength_nm: 656.3}\nparams: {}\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(registry, "instrument_dir", lambda name: copy)
    inst = registry.load_instrument("highspec")
    assert set(inst.configurations) >= {"mgb518", "halpha"}
    assert inst.configurations["halpha"].geometry.central_wavelength_nm == pytest.approx(656.3)


def test_anchor_trace_loads_with_default():
    assert load_anchor_trace("highspec", "mgb518", "lab3") == 0
    assert load_anchor_trace("highspec", "mgb518", "survey45") == 0  # no file: default
