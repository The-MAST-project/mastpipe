"""Tests for header-map resolution and exposure classification."""

import numpy as np
import pytest
from astropy.io import fits

from mastpipe.instruments.registry import load_header_map
from mastpipe.io.headers import HeaderMap
from mastpipe.io.ingest import scan_folder

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
    spec = HeaderMap.from_dict({"fields": {"mid_time": {"keys": ["DATE-OBS"], "type": "str", "reference": "end"}}}).fields[
        "mid_time"
    ]
    assert spec.reference == "end"
    plain = HeaderMap.from_dict({"fields": {"mid_time": {"keys": ["DATE-OBS"]}}})
    assert plain.fields["mid_time"].reference == "start"
    with pytest.raises(ValueError, match="reference"):
        HeaderMap.from_dict({"fields": {"mid_time": {"keys": ["DATE-OBS"], "reference": "late"}}})
    highspec = HeaderMap.from_dict(load_header_map("highspec"))
    assert highspec.fields["mid_time"].reference == "start"


def test_timezone_option_is_validated():
    def spec(**options):
        raw = {"fields": {"mid_time": {"keys": ["DATE-OBS"], "type": "str", **options}}}
        return HeaderMap.from_dict(raw).fields["mid_time"]

    assert spec().timezone == "utc"
    for zone in ("utc", "assume-utc", "Asia/Jerusalem"):
        assert spec(timezone=zone).timezone == zone
    # an IANA name for UTC is the utc keyword (core.times.check_timezone)
    for zone in ("UTC", "Etc/UTC"):
        assert spec(timezone=zone).timezone == "utc"
    for bad in ("Mars/Olympus", "", "../etc", 3):
        with pytest.raises(ValueError, match="timezone"):
            spec(timezone=bad)
    highspec = HeaderMap.from_dict(load_header_map("highspec"))
    assert highspec.fields["mid_time"].timezone == "assume-utc"
    assert highspec.fields["mid_time"].keys == ("DATE-OBS", "FRAME", "DATE")


def test_lt_on_a_non_numeric_value_does_not_match():
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    base = {name: None for name in hmap.fields}
    # an untyped string where a number is compared: the rule does not match, nothing raises
    assert hmap.classify({**base, "exposure_time": "short"}) == ("unknown", "unclassified")
    assert hmap.classify({**base, "exposure_time": 0.0}) == ("bias", "rule:0")
    assert hmap.classify({**base, "exposure_time": True}) == ("unknown", "unclassified")


# ------------------------------------------------------------- the provisional light type


def test_light_is_provisional_and_rules_refine_it(highspec_map):
    base = {name: None for name in highspec_map.fields}
    light = {**base, "exposure_type": "light"}
    assert highspec_map.classify({**light, "exposure_time": 30.0}) == ("light", "header")
    assert highspec_map.classify({**light, "exposure_time": 0.0}) == ("bias", "rule:0")
    # a rule with conditions refines it; the unconditional fallback never does
    lamp = {**light, "exposure_time": 30.0, "lamp_state": "thar"}
    assert highspec_map.classify(lamp) == ("arc", "rule:2")
    # without a valid header type the rules and the fallback apply as before
    assert highspec_map.classify({**base, "exposure_time": 30.0}) == ("science", "rule:5")


def test_light_with_zero_exposure_is_bias(tmp_path, highspec_map):
    cards = {"IMAGETYP": "Light", "DATE-OBS": "2026-10-01T18:00:00", "INSTRUME": "Highspec"}
    for name, exposure in (("zero.fits", 0.0), ("long.fits", 30.0)):
        hdu = fits.PrimaryHDU(np.zeros((4, 4), dtype=np.uint16))
        hdu.header.update({**cards, "EXPOSURE": exposure})
        hdu.writeto(tmp_path / name)
    refs = {ref.meta.source: ref.meta for ref in scan_folder(tmp_path, "highspec", highspec_map)}
    assert (refs["zero.fits"].exposure_type, refs["zero.fits"].classification_rule) == (
        "bias",
        "rule:0",
    )
    assert (refs["long.fits"].exposure_type, refs["long.fits"].classification_rule) == (
        "light",
        "header",
    )


def test_flat_header_type_is_never_refined(highspec_map):
    base = {name: None for name in highspec_map.fields}
    for kind in ("flat", "bias", "dark", "arc", "science"):
        fields = {**base, "exposure_type": kind, "exposure_time": 0.0, "lamp_state": "thar"}
        assert highspec_map.classify(fields) == (kind, "header")


def test_a_map_without_conditional_rules_leaves_light_as_is():
    raw = {
        "fields": {
            "exposure_type": {"keys": ["IMAGETYP"]},
            "exposure_time": {"keys": ["EXPOSURE"]},
        },
        "classification": [{"type": "science"}],
    }
    hmap = HeaderMap.from_dict(raw)
    assert hmap.classify({"exposure_type": "light", "exposure_time": 0.0}) == ("light", "header")
    assert hmap.classify({"exposure_type": None, "exposure_time": 0.0}) == ("science", "rule:0")


# ------------------------------------------------------------------------- per_unit fields

PER_UNIT_MAP = {
    "fields": {
        "targets": {"per_unit": "OBJ{nn}", "missing_values": ["none"]},
        "target": {"keys": ["OBJECT"]},
    }
}


def test_per_unit_field_collects_present_units_in_unit_order():
    hmap = HeaderMap.from_dict(PER_UNIT_MAP)
    header = {"OBJ03": "HD 3", "OBJ01": "HD 1", "OBJ02": "none", "OBJ04": "", "OBJ21": "x"}
    fields, unresolved = hmap.resolve(header, "frame.fits")
    assert fields["targets"] == (("mast01", "HD 1"), ("mast03", "HD 3"))
    assert unresolved == ["target"]
    fields, unresolved = hmap.resolve({"OBJECT": "star"}, "frame.fits")
    assert fields["targets"] is None and fields["target"] == "star"
    assert unresolved == ["targets"]


def test_per_unit_spec_is_validated():
    with pytest.raises(ValueError, match="nn"):
        HeaderMap.from_dict({"fields": {"targets": {"per_unit": "OBJ"}}})
    for other in ({"keys": ["OBJECT"]}, {"filename_regex": "(?P<value>x)"}, {"default": "x"}):
        with pytest.raises(ValueError, match="per_unit"):
            HeaderMap.from_dict({"fields": {"targets": {"per_unit": "OBJ{nn}", **other}}})


def test_highspec_map_reads_targets_per_unit(highspec_map):
    assert highspec_map.fields["targets"].per_unit == "OBJ{nn}"
    assert highspec_map.fields["target"].keys == ("OBJECT",)


def test_amplifier_is_read_and_keyed(tmp_path, highspec_map):
    cards = {
        "IMAGETYP": "Light",
        "EXPOSURE": 30.0,
        "DATE-OBS": "2026-10-01T18:00:00",
        "INSTRUME": "Highspec",
        "READMODE": "Image",
    }
    for name, amplifier in (("em.fits", "Electron Multiplying"), ("conv.fits", "Conventional")):
        hdu = fits.PrimaryHDU(np.zeros((4, 4), dtype=np.uint16))
        hdu.header.update({**cards, "OUTPTAMP": amplifier})
        hdu.writeto(tmp_path / name)
    hdu = fits.PrimaryHDU(np.zeros((4, 4), dtype=np.uint16))
    hdu.header.update(cards)
    hdu.writeto(tmp_path / "bare.fits")
    refs = {ref.meta.source: ref.meta for ref in scan_folder(tmp_path, "highspec", highspec_map)}
    assert refs["em.fits"].amplifier == "em"
    assert refs["em.fits"].readout_mode == "image"  # OUTPTAMP is no longer a readout-mode key
    assert refs["em.fits"].readout_key == "image/em"
    assert refs["conv.fits"].readout_key == "image/conventional"
    assert refs["bare.fits"].amplifier is None
    assert refs["bare.fits"].readout_key == "image/unknown"


def test_an_undefined_card_counts_as_absent():
    """A FITS card with no value holds astropy's ``Undefined`` (an astropy ``Header`` hands
    it out as ``None``, a mapping built from the cards as is): either way it is absent, the
    next keyword is tried, and a field with nothing else is unresolved."""
    hmap = HeaderMap.from_dict(SYNTHETIC_MAP)
    exposure = fits.Card.fromstring("EXPOSURE=".ljust(80))
    blank_object = fits.Card.fromstring("OBJECT  =".ljust(80))
    assert isinstance(exposure.value, fits.card.Undefined)
    cards = {card.keyword: card.value for card in (exposure, blank_object)}
    header = fits.Header([exposure, blank_object, ("EXPTIME", 5.0)])
    for mapping in ({**cards, "EXPTIME": 5.0}, header):
        fields, unresolved = hmap.resolve(mapping, "x.fits")
        assert fields["exposure_time"] == 5.0
        assert fields["target"] is None and "target" in unresolved


def test_targets_must_be_read_per_unit():
    """Ingest reads ``targets`` as ``(unit, target)`` pairs, so a map that declares it with
    ``keys`` (one value for the frame) is rejected when it is loaded; ``target`` is the
    field for that."""
    with pytest.raises(ValueError, match="targets.*per_unit"):
        HeaderMap.from_dict({"fields": {"targets": {"keys": ["OBJECT"]}}})
    with pytest.raises(ValueError, match="targets.*per_unit"):
        HeaderMap.from_dict({"fields": {"targets": {"default": "sky"}}})
    assert HeaderMap.from_dict({"fields": {"targets": {"per_unit": "OBJ{nn}"}}}).fields["targets"].per_unit
