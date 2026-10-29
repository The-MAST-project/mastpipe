"""Tests for header-map resolution and exposure classification."""

import pytest

from mastpipe.instruments.registry import load_header_map
from mastpipe.io.headers import HeaderMap

SYNTHETIC_MAP = {
    "fields": {
        "exposure_time": {"keys": ["EXPOSURE", "EXPTIME"], "type": "float"},
        "binning_x": {"keys": ["HBIN"], "type": "int", "default": 1},
        "exposure_type": {
            "keys": ["IMAGETYP"],
            "lower": True,
            "map": {"zero": "bias", "object": "science"},
        },
        "readout_mode": {"keys": ["READMODE"], "lower": True},
        "configuration": {
            "keys": ["GRATING"],
            "lower": True,
            "filename_regex": "(?P<value>mgb518|halpha)",
        },
        "target": {"keys": ["OBJECT"]},
    },
    "classification": [
        {"when": {"exposure_time": {"lt": 0.001}}, "type": "bias"},
        {"when": {"exposure_time": {"missing": True}}, "type": "unknown"},
    ],
}


def test_resolve_priority_and_transforms():
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    header = {
        "EXPTIME": "7.5",
        "EXPOSURE": "12.5",
        "IMAGETYP": "OBJECT",
        "READMODE": "Fast",
        "GRATING": "MGB518",
    }
    fields, unresolved = hmap.resolve(header, "frame_0001.fits")
    assert fields["exposure_time"] == pytest.approx(12.5)  # first key in the list wins
    assert isinstance(fields["exposure_time"], float)
    assert fields["exposure_type"] == "science"  # lower-cased, then mapped
    assert fields["readout_mode"] == "fast"  # lower-cased, left as is when unmapped
    assert fields["configuration"] == "mgb518"
    assert fields["binning_x"] == 1  # default used when no key is present
    assert fields["target"] is None
    assert unresolved == ["target"]
    assert set(fields) == set(SYNTHETIC_MAP["fields"])


def test_unmapped_value_is_left_as_is():
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    fields, _ = hmap.resolve({"IMAGETYP": "Comp"}, "x.fits")
    assert fields["exposure_type"] == "comp"


def test_unconvertible_value_is_unresolved():
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    fields, unresolved = hmap.resolve({"EXPOSURE": "n/a"}, "x.fits")
    assert fields["exposure_time"] is None
    assert "exposure_time" in unresolved


def test_filename_regex_fallback():
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    fields, unresolved = hmap.resolve({}, "night01_halpha_0007.fits")
    assert fields["configuration"] == "halpha"
    assert "configuration" not in unresolved
    fields, unresolved = hmap.resolve({}, "night01_0007.fits")
    assert fields["configuration"] is None
    assert "configuration" in unresolved
    # a present keyword beats the file name
    fields, _ = hmap.resolve({"GRATING": "mgb518"}, "night01_halpha_0007.fits")
    assert fields["configuration"] == "mgb518"


def test_classify_header_wins_then_rules():
    hmap = HeaderMap.from_dict(load_header_map("highspec"))
    base = {name: None for name in hmap.fields}
    assert hmap.classify({**base, "exposure_type": "bias", "lamp_state": "thar"}) == (
        "bias",
        "header",
    )
    # an unmapped header value is not a valid type: the rules decide
    assert hmap.classify({**base, "exposure_type": "weird", "exposure_time": 0.0}) == (
        "bias",
        "rule:0",
    )
    assert hmap.classify({**base, "exposure_time": 10.0, "lamp_state": "thar"}) == (
        "arc",
        "rule:2",
    )
    dark = {**base, "exposure_time": 10.0, "lamp_state": "off", "shutter_state": "closed"}
    assert hmap.classify(dark) == ("dark", "rule:1")
    # both conditions of the dark rule must hold
    assert hmap.classify({**base, "exposure_time": 10.0, "lamp_state": "off"}) == (
        "science",
        "rule:5",
    )
    assert hmap.classify({**base, "exposure_time": 10.0, "lamp_state": "flat"}) == (
        "flat",
        "rule:3",
    )


def test_missing_exposure_time_is_unknown():
    hmap = HeaderMap.from_dict(load_header_map("highspec"))
    fields, unresolved = hmap.resolve({"OBJECT": "lab", "LAMP": "thar"}, "lab_0001.fits")
    assert fields["exposure_time"] is None
    assert "exposure_time" in unresolved
    # the lamp rule comes first: a lit lamp is still an arc without an exposure time
    assert hmap.classify(fields) == ("arc", "rule:2")
    fields, _ = hmap.resolve({"OBJECT": "lab"}, "lab_0001.fits")
    assert hmap.classify(fields) == ("unknown", "rule:4")


def test_no_matching_rule_is_unclassified():
    hmap = HeaderMap.from_dict({"fields": SYNTHETIC_MAP["fields"], "classification": []})
    assert hmap.classify({"exposure_time": 5.0}) == ("unknown", "unclassified")
    hmap = HeaderMap.from_dict({"fields": SYNTHETIC_MAP["fields"]})
    assert hmap.classify({}) == ("unknown", "unclassified")


@pytest.mark.parametrize(
    "bad",
    [
        {"fields": {"exposure_time": {"keys": ["EXPOSURE"], "type": "complex"}}},
        {"fields": {"exposure_time": {"keys": ["EXPOSURE"], "typo": 1}}},
        {"fields": {"exposure_time": {}}},
        {"fields": {}, "classification": [{"type": "nonsense"}]},
        {"fields": {}, "classification": [{"when": {"missing_field": {"eq": 1}}, "type": "bias"}]},
        {
            "fields": {"exposure_time": {"keys": ["EXPOSURE"]}},
            "classification": [{"when": {"exposure_time": {"gt": 1}}, "type": "bias"}],
        },
        {"fields": {"configuration": {"filename_regex": "(mgb518)"}}},
    ],
)
def test_invalid_map_is_rejected(bad):
    with pytest.raises(ValueError):
        HeaderMap.from_dict(bad)
