"""Tests for the configuration models and the layered loader."""

import pytest
import yaml
from pydantic import ValidationError

from mastpipe.config.loader import (
    build_effective,
    config_hash,
    deep_merge,
    load_yaml,
    parse_override,
)
from mastpipe.config.models import (
    ConfigurationBlock,
    DetectorConfig,
    GeometryConfig,
    InstrumentConfig,
    StepParams,
    TraceParams,
)


def make_geometry(**overrides) -> GeometryConfig:
    fields = {
        "model": "grating",
        "groove_spacing_nm": 367.0,
        "incidence_angle_deg": 45.0,
        "order": 1,
        "camera_focal_mm": 479.41,
        "central_wavelength_nm": 518.0,
    }
    fields.update(overrides)
    return GeometryConfig(**fields)


def make_detector(**overrides) -> DetectorConfig:
    fields = {
        "name": "main",
        "shape": (400, 1600),
        "gain": 1.0,
        "read_noise": 5.0,
        "pixel_size_um": 16.0,
        "saturation": 65000.0,
        "crop": (1, 1599, 1, 399),
    }
    fields.update(overrides)
    return DetectorConfig(**fields)


def make_instrument(
    block_params=None, instrument_params=None, layout_params=None
) -> InstrumentConfig:
    if block_params is None:
        block_params = {"trace": {"trace_halfwidth": 4}}
    return InstrumentConfig(
        name="toy",
        cameras={"main": make_detector()},
        configurations={
            "g518": ConfigurationBlock(
                name="g518",
                band_nm=(513.0, 527.5),
                resolving_power=20000.0,
                geometry=make_geometry(),
                line_table="lines/g518.csv",
                params=block_params,
            )
        },
        params=StepParams(**(instrument_params or {})),
        default_layout="lab3",
        layouts=["lab3", "survey45"],
        layout_params=layout_params or {},
    )


def test_deep_merge_override_wins_and_lists_replace():
    base = {"a": {"x": 1, "y": [1, 2], "z": {"p": 1}}, "b": 2, "l": [1, 2, 3]}
    out = deep_merge(base, {"a": {"y": [9], "z": {"q": 2}}, "c": 3, "l": [4]})
    assert out == {"a": {"x": 1, "y": [9], "z": {"p": 1, "q": 2}}, "b": 2, "c": 3, "l": [4]}
    # a scalar in the override replaces a mapping in the base
    assert deep_merge({"a": {"x": 1}}, {"a": 5}) == {"a": 5}
    # neither input is mutated and the result shares nothing with them
    assert base["a"]["y"] == [1, 2] and "c" not in base
    out["a"]["z"]["p"] = 99
    assert base["a"]["z"]["p"] == 1


def test_build_effective_layers_block_over_instrument():
    inst = make_instrument(instrument_params={"trace": {"trace_halfwidth": 10, "x_bin": 30}})
    eff = build_effective(inst, "g518", "main")
    assert eff.params.trace.trace_halfwidth == 4  # block wins over instrument
    assert eff.params.trace.x_bin == 30  # untouched instrument value survives
    assert eff.params.extract.halfwidth == 9.0  # defaults fill the rest
    assert eff.instrument == "toy"
    assert eff.camera == "main"
    assert eff.configuration == "g518"
    assert eff.detector.shape == (400, 1600)
    assert eff.geometry.central_wavelength_nm == 518.0
    assert eff.band_nm == (513.0, 527.5)
    assert eff.resolving_power == 20000.0
    assert eff.line_table == "lines/g518.csv"
    assert eff.layers == ["instrument:toy", "configuration:g518"]

    eff2 = build_effective(
        inst,
        "g518",
        "main",
        overrides=[{"trace": {"trace_halfwidth": 6}}, {"extract": {"method": "mean"}}],
    )
    assert eff2.params.trace.trace_halfwidth == 6  # override wins over block
    assert eff2.params.trace.x_bin == 30
    assert eff2.params.extract.method == "mean"
    assert eff2.layers == [
        "instrument:toy",
        "configuration:g518",
        "override:trace.trace_halfwidth=6",
        'override:extract.method="mean"',
    ]
    # inputs are untouched by layering
    assert inst.params.trace.trace_halfwidth == 10
    assert inst.configurations["g518"].params == {"trace": {"trace_halfwidth": 4}}


def test_build_effective_layers_layout_params_between_block_and_overrides():
    inst = make_instrument(
        block_params={"trace": {"trace_halfwidth": 7}, "extract": {"halfwidth": 9.0}},
        layout_params={"survey45": {"trace": {"trace_halfwidth": 4}, "extract": {"halfwidth": 4}}},
    )
    plain = build_effective(inst, "g518", "main")
    assert plain.params.trace.trace_halfwidth == 7 and plain.params.extract.halfwidth == 9.0
    assert plain.layers == ["instrument:toy", "configuration:g518"]
    lab = build_effective(inst, "g518", "main", layout="lab3")  # no block for this preset
    assert lab.params == plain.params
    assert lab.layers == ["instrument:toy", "configuration:g518", "layout:lab3"]
    survey = build_effective(inst, "g518", "main", layout="survey45")
    assert survey.params.trace.trace_halfwidth == 4 and survey.params.extract.halfwidth == 4.0
    assert survey.layers == ["instrument:toy", "configuration:g518", "layout:survey45"]
    # an override still wins over the preset block
    forced = build_effective(
        inst, "g518", "main", overrides=[{"trace": {"trace_halfwidth": 5}}], layout="survey45"
    )
    assert forced.params.trace.trace_halfwidth == 5 and forced.params.extract.halfwidth == 4.0
    assert config_hash(survey) != config_hash(plain) and config_hash(lab) == config_hash(plain)
    with pytest.raises(KeyError, match="layout"):
        build_effective(inst, "g518", "main", layout="nope")
    # a parameter block for a preset the instrument does not list is rejected
    with pytest.raises(ValidationError, match="layout_params"):
        make_instrument(layout_params={"nope": {"trace": {"trace_halfwidth": 4}}})
    # the shipped HighSpec definition narrows the survey halfwidths to the 8.5 px pitch
    from mastpipe.instruments.registry import load_instrument

    highspec = load_instrument("highspec")
    eff = build_effective(highspec, "Mg", "main", layout="survey45")
    assert eff.params.trace.trace_halfwidth == 4
    assert eff.params.lines.halfwidth == 4 and eff.params.extract.halfwidth == 4.0
    default = build_effective(highspec, "Mg", "main", layout="lab3")
    assert default.params.trace.trace_halfwidth == 8 and default.params.extract.halfwidth == 9.0


def test_build_effective_unknown_configuration_or_camera():
    inst = make_instrument()
    with pytest.raises(KeyError, match="g518"):
        build_effective(inst, "nope", "main")
    with pytest.raises(KeyError, match="main"):
        build_effective(inst, "g518", "nope")


def test_build_effective_rejects_bad_override():
    inst = make_instrument()
    with pytest.raises(ValidationError):
        build_effective(inst, "g518", "main", overrides=[{"trace": {"banana": 1}}])
    with pytest.raises(ValidationError):
        build_effective(inst, "g518", "main", overrides=[{"extract": {"aperture": "triangle"}}])


def test_config_hash_stable_and_sensitive():
    inst = make_instrument()
    a = config_hash(build_effective(inst, "g518", "main"))
    b = config_hash(build_effective(make_instrument(), "g518", "main"))
    assert a == b
    assert len(a) == 12 and int(a, 16) >= 0
    c = config_hash(build_effective(inst, "g518", "main", overrides=[{"trace": {"x_bin": 10}}]))
    assert c != a
    # a detector change is part of the effective configuration too
    other = make_instrument()
    other = other.model_copy(update={"cameras": {"main": make_detector(gain=2.0)}})
    assert config_hash(build_effective(other, "g518", "main")) != a
    # an override that restates the effective value does not change the hash
    same = config_hash(build_effective(inst, "g518", "main", overrides=[{"trace": {"x_bin": 25}}]))
    assert same == a


def test_classify_params_layer_like_any_step():
    inst = make_instrument()
    plain = build_effective(inst, "g518", "main")
    assert plain.params.classify.arc_peak_fraction == 0.5
    tuned = build_effective(inst, "g518", "main", overrides=[{"classify": {"signal_sigma": 8.0}}])
    assert tuned.params.classify.signal_sigma == 8.0
    assert tuned.params.classify.continuum_window == 101
    assert config_hash(tuned) != config_hash(plain)
    with pytest.raises(ValidationError):
        build_effective(inst, "g518", "main", overrides=[{"classify": {"banana": 1}}])


def test_unknown_key_rejected():
    with pytest.raises(ValidationError):
        TraceParams(trace_halfwidth=4, banana=1)
    with pytest.raises(ValidationError):
        StepParams(trace={"banana": 1})
    with pytest.raises(ValidationError):
        make_detector(banana=1)
    with pytest.raises(ValidationError):
        make_geometry(model="prism")
    inst = make_instrument()
    with pytest.raises(ValidationError):
        InstrumentConfig(**inst.model_dump(), banana=1)


def test_parse_override_types():
    assert parse_override("x=1") == {"x": 1}
    assert isinstance(parse_override("x=1")["x"], int)
    assert parse_override("x=1.5") == {"x": 1.5}
    assert parse_override("x=true") == {"x": True}
    assert parse_override("x=abc") == {"x": "abc"}
    assert parse_override("x=[1,2]") == {"x": [1, 2]}
    assert parse_override("x=null") == {"x": None}
    assert parse_override("a.b.c=1.5") == {"a": {"b": {"c": 1.5}}}
    assert parse_override("x=a=b") == {"x": "a=b"}
    assert parse_override("x=") == {"x": ""}
    for bad in ("x", "=1", "a..b=1", ".a=1"):
        with pytest.raises(ValueError):
            parse_override(bad)


def test_load_yaml_reads_mapping(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("name: toy\nparams:\n  trace: {x_bin: 10}\nlayouts: [a, b]\n", encoding="utf-8")
    assert load_yaml(path) == {
        "name": "toy",
        "params": {"trace": {"x_bin": 10}},
        "layouts": ["a", "b"],
    }
    empty = tmp_path / "empty.yaml"
    empty.write_text("", encoding="utf-8")
    assert load_yaml(empty) == {}
    sequence = tmp_path / "list.yaml"
    sequence.write_text("- 1\n- 2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="mapping"):
        load_yaml(sequence)


def test_instrument_round_trips_through_yaml(tmp_path):
    inst = make_instrument()
    path = tmp_path / "toy.yaml"
    path.write_text(yaml.safe_dump(inst.model_dump(mode="json")), encoding="utf-8")
    loaded = InstrumentConfig.model_validate(load_yaml(path))
    assert loaded == inst
    assert loaded.arc_classification == "lit_traces"
    assert loaded.cameras["main"].crop == (1, 1599, 1, 399)
