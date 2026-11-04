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
            "filename_regex": "(?P<value>mgb518|halpha)",
            "map": {"mgb518": "Mg", "halpha": "Halpha"},
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
        "GRATING": "Mg",
    }
    fields, unresolved = hmap.resolve(header, "frame_0001.fits")
    assert fields["exposure_time"] == pytest.approx(12.5)  # first key in the list wins
    assert isinstance(fields["exposure_time"], float)
    assert fields["exposure_type"] == "science"  # lower-cased, then mapped
    assert fields["readout_mode"] == "fast"  # lower-cased, left as is when unmapped
    assert fields["configuration"] == "Mg"  # case kept
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
    assert fields["configuration"] == "Halpha"  # the file-name alias, mapped
    assert "configuration" not in unresolved
    fields, unresolved = hmap.resolve({}, "night01_0007.fits")
    assert fields["configuration"] is None
    assert "configuration" in unresolved
    # a present keyword beats the file name
    fields, _ = hmap.resolve({"GRATING": "Mg"}, "night01_halpha_0007.fits")
    assert fields["configuration"] == "Mg"


def test_configuration_keeps_case_and_maps_lab_aliases():
    hm = HeaderMap.from_dict(load_header_map("highspec"))
    assert hm.resolve({"GRATING": "Mg"}, "x.fits")[0]["configuration"] == "Mg"
    assert hm.resolve({}, "FOCUS_00_mgb518.fits")[0]["configuration"] == "Mg"
    assert hm.resolve({}, "x_halpha.fits")[0]["configuration"] == "Halpha"


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
        {"fields": {"configuration": {"filename_regex": "(Mg)"}}},
    ],
)
def test_invalid_map_is_rejected(bad):
    with pytest.raises(ValueError):
        HeaderMap.from_dict(bad)


def test_missing_values_make_a_field_unresolved():
    hmap = HeaderMap.from_dict(
        {
            "fields": {
                "detector_temperature": {
                    "keys": ["TEMP", "CCD-TEMP"],
                    "type": "float",
                    "missing_values": [-999],
                },
                "exposure_time": {"keys": ["EXPOSURE"], "type": "float"},
            }
        }
    )
    fields, unresolved = hmap.resolve({"TEMP": -999, "EXPOSURE": 1.0}, "x.fits")
    assert fields["detector_temperature"] is None and "detector_temperature" in unresolved
    # the sentinel is caught after conversion too, and the next key is tried
    fields, unresolved = hmap.resolve({"TEMP": "-999.0", "CCD-TEMP": -60.5}, "x.fits")
    assert fields["detector_temperature"] == pytest.approx(-60.5)
    assert "detector_temperature" not in unresolved
    fields, _ = hmap.resolve({"TEMP": -60.0}, "x.fits")
    assert fields["detector_temperature"] == pytest.approx(-60.0)
    # the HighSpec map lists the Andor sentinel
    highspec = HeaderMap.from_dict(load_header_map("highspec"))
    fields, unresolved = highspec.resolve({"TEMP": -999, "EXPOSURE": 300.0}, "x.fits")
    assert fields["detector_temperature"] is None and "detector_temperature" in unresolved
    with pytest.raises(ValueError, match="missing_values"):
        HeaderMap.from_dict({"fields": {"t": {"keys": ["T"], "missing_values": 7}}})


def test_time_reference_option_is_validated():
    spec = HeaderMap.from_dict(
        {"fields": {"mid_time": {"keys": ["DATE-OBS"], "type": "str", "reference": "end"}}}
    ).fields["mid_time"]
    assert spec.reference == "end"
    plain = HeaderMap.from_dict({"fields": {"mid_time": {"keys": ["DATE-OBS"]}}})
    assert plain.fields["mid_time"].reference == "start"
    with pytest.raises(ValueError, match="reference"):
        HeaderMap.from_dict({"fields": {"mid_time": {"keys": ["DATE-OBS"], "reference": "late"}}})
    highspec = HeaderMap.from_dict(load_header_map("highspec"))
    assert highspec.fields["mid_time"].reference == "start"


def test_lt_on_a_non_numeric_value_does_not_match():
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    base = {name: None for name in hmap.fields}
    # an untyped string where a number is compared: the rule does not match, nothing raises
    assert hmap.classify({**base, "exposure_time": "short"}) == ("unknown", "unclassified")
    assert hmap.classify({**base, "exposure_time": 0.0}) == ("bias", "rule:0")
    assert hmap.classify({**base, "exposure_time": True}) == ("unknown", "unclassified")
