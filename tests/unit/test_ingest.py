"""Tests for folder scanning, frame loading, manifests and checksums."""

import csv
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.config.models import DetectorConfig
from mastpipe.core.model import FrameMeta, FrameRef
from mastpipe.instruments.registry import load_header_map, load_instrument
from mastpipe.io.headers import HeaderMap
from mastpipe.io.ingest import (
    apply_manifest,
    expected_shape,
    file_checksum,
    load_frame,
    read_manifest,
    scan_folder,
)


def write_fits(path: Path, data: np.ndarray | None, cards: dict) -> None:
    hdu = fits.PrimaryHDU(data)
    hdu.header.update(cards)
    hdu.writeto(path)


def write_manifest(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_detector(**overrides) -> DetectorConfig:
    fields = {
        "name": "main",
        "shape": (4, 6),
        "gain": 1.0,
        "read_noise": 1.0,
        "pixel_size_um": 16.0,
        "saturation": 65535.0,
        "crop": None,
    }
    fields.update(overrides)
    return DetectorConfig(**fields)


@pytest.fixture
def highspec_map() -> HeaderMap:
    return HeaderMap.from_dict(load_header_map("highspec"))


@pytest.fixture
def night(tmp_path: Path) -> Path:
    data = np.arange(12, dtype=np.uint16).reshape(3, 4)
    write_fits(
        tmp_path / "a_bias.fits",
        data,
        {"IMAGETYP": "Bias", "EXPOSURE": 0.0, "DATE-OBS": "2026-10-01T18:00:00", "HBIN": 2},
    )
    write_fits(
        tmp_path / "b_science.fits",
        data + 100,
        {
            "EXPOSURE": 120.0,
            "DATE-OBS": "2026-10-01T18:05:00",
            "OBJECT": "lab target",
            "HBIN": 2,
            "VBIN": 2,
            "TEMP": -60.5,
            "READMODE": "Normal",
            "GRATING": "Mg",
        },
    )
    (tmp_path / "notes.txt").write_text("not a frame\n", encoding="utf-8")
    return tmp_path


def test_file_checksum_is_short_sha256(tmp_path: Path):
    path = tmp_path / "blob.bin"
    path.write_bytes(bytes(range(256)) * 10)
    expected = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    assert file_checksum(path) == expected
    assert len(file_checksum(path)) == 16


def test_scan_folder_builds_refs_with_checksums(night: Path, highspec_map: HeaderMap):
    refs = scan_folder(night, "highspec", highspec_map)
    assert [ref.meta.source for ref in refs] == ["a_bias.fits", "b_science.fits"]
    assert all(isinstance(ref, FrameRef) for ref in refs)
    assert all(ref.extension == 0 and ref.meta.extension == 0 for ref in refs)
    assert all(ref.meta.instrument == "highspec" and ref.meta.camera == "main" for ref in refs)
    assert [ref.meta.checksum for ref in refs] == [file_checksum(ref.path) for ref in refs]
    assert refs[0].meta.checksum != refs[1].meta.checksum
    assert refs[0].path == night / "a_bias.fits"

    bias = refs[0].meta
    assert bias.exposure_type == "bias" and bias.classification_rule == "header"
    assert bias.exposure_time == 0.0
    assert bias.mid_time == "2026-10-01T18:00:00"  # a zero-second exposure: start is mid
    assert bias.binning == (2, 1)  # VBIN absent: the missing axis defaults to 1
    assert bias.targets == () and bias.arc_kind is None
    assert "targets" in bias.unresolved and "exposure_time" not in bias.unresolved
    assert "target" not in bias.unresolved

    sci = refs[1].meta
    assert sci.exposure_type == "science" and sci.classification_rule == "rule:5"
    assert sci.exposure_time == pytest.approx(120.0)
    assert sci.binning == (2, 2)
    assert sci.detector_temperature == pytest.approx(-60.5)
    assert sci.readout_mode == "normal"
    assert sci.configuration == "Mg"
    assert sci.targets == (("*", "lab target"),)  # a single OBJECT applies to every fiber
    assert "targets" not in sci.unresolved and "target" not in sci.unresolved


def test_scan_folder_accepts_raw_mapping_and_pattern(night: Path):
    refs = scan_folder(night, "highspec", load_header_map("highspec"), pattern="a_*.fits")
    assert [ref.meta.source for ref in refs] == ["a_bias.fits"]


def test_missing_exposure_time_is_unknown(tmp_path: Path, highspec_map: HeaderMap):
    write_fits(tmp_path / "lab_0001.fits", np.zeros((2, 2), dtype=np.uint16), {"OBJECT": "lab"})
    (ref,) = scan_folder(tmp_path, "highspec", highspec_map)
    assert ref.meta.exposure_type == "unknown"
    assert ref.meta.classification_rule == "rule:4"
    assert ref.meta.exposure_time is None
    assert "exposure_time" in ref.meta.unresolved


def test_mid_time_is_derived_from_the_start_time(tmp_path: Path, highspec_map: HeaderMap):
    """HighSpec writes the exposure start; a 300 s frame's mid-time is 150 s later. Without
    an exposure time the time stays as given and ``unresolved`` says so."""
    cards = {"IMAGETYP": "light", "EXPOSURE": 300.0, "DATE-OBS": "2026-10-01T20:00:00"}
    write_fits(tmp_path / "sci.fits", np.zeros((2, 2), dtype=np.uint16), cards)
    write_fits(
        tmp_path / "fraction.fits",
        np.zeros((2, 2), dtype=np.uint16),
        {"IMAGETYP": "flat", "EXPOSURE": 5.0, "DATE-OBS": "2026-10-01T20:00:00.5"},
    )
    write_fits(
        tmp_path / "notime.fits",
        np.zeros((2, 2), dtype=np.uint16),
        {"LAMP": "thar", "DATE-OBS": "2026-10-01T21:00:00"},
    )
    refs = {ref.meta.source: ref.meta for ref in scan_folder(tmp_path, "highspec", highspec_map)}
    assert refs["sci.fits"].mid_time == "2026-10-01T20:02:30"
    assert "mid_time" not in refs["sci.fits"].unresolved
    assert refs["fraction.fits"].mid_time == "2026-10-01T20:00:03"
    assert refs["notime.fits"].exposure_time is None
    assert refs["notime.fits"].mid_time == "2026-10-01T21:00:00"
    assert "mid_time" in refs["notime.fits"].unresolved

    # the reference of the time field decides the arithmetic
    for reference, expected in (("mid", "2026-10-01T20:00:00"), ("end", "2026-10-01T19:57:30")):
        raw = load_header_map("highspec")
        raw["fields"]["mid_time"] = {**raw["fields"]["mid_time"], "reference": reference}
        (ref,) = scan_folder(tmp_path, "highspec", HeaderMap.from_dict(raw), pattern="sci.fits")
        assert ref.meta.mid_time == expected, reference


def test_frame_meta_records_the_timezone_used(tmp_path: Path, highspec_map: HeaderMap):
    """A naive FRAME stamp under HighSpec's ``assume-utc`` is kept as UTC and labelled
    ``assumed-utc``; a stamp with its own offset is converted and labelled ``utc``; a
    zone named in the map converts a naive stamp and is the label."""
    blank = np.zeros((2, 2), dtype=np.uint16)
    write_fits(
        tmp_path / "frame.fits",
        blank,
        {"IMAGETYP": "light", "EXPOSURE": 60.0, "FRAME": "2026-09-09T11:35:24.500"},
    )
    write_fits(
        tmp_path / "offset.fits",
        blank,
        {"IMAGETYP": "light", "EXPOSURE": 60.0, "DATE-OBS": "2026-11-04T20:00:00+02:00"},
    )
    write_fits(tmp_path / "untimed.fits", blank, {"IMAGETYP": "light", "EXPOSURE": 60.0})
    refs = {ref.meta.source: ref.meta for ref in scan_folder(tmp_path, "highspec", highspec_map)}
    assert refs["frame.fits"].mid_time == "2026-09-09T11:35:54"
    assert refs["frame.fits"].timezone == "assumed-utc"
    assert refs["offset.fits"].mid_time == "2026-11-04T18:00:30"
    assert refs["offset.fits"].timezone == "utc"
    assert refs["untimed.fits"].mid_time is None and refs["untimed.fits"].timezone is None

    raw = load_header_map("highspec")
    raw["fields"]["mid_time"] = {**raw["fields"]["mid_time"], "timezone": "Asia/Jerusalem"}
    local = {ref.meta.source: ref.meta for ref in scan_folder(tmp_path, "highspec", HeaderMap.from_dict(raw))}
    assert local["frame.fits"].mid_time == "2026-09-09T08:35:54"
    assert local["frame.fits"].timezone == "Asia/Jerusalem"
    assert local["offset.fits"].timezone == "utc"

    # a manifest mid_time is taken as given, in UTC
    manifest = tmp_path / "frames.csv"
    write_manifest(manifest, [{"source": "frame.fits", "mid_time": "2026-09-09T12:00:00"}])
    (ref,) = scan_folder(tmp_path, "highspec", highspec_map, pattern="frame.fits", manifest=manifest)
    assert ref.meta.mid_time == "2026-09-09T12:00:00" and ref.meta.timezone == "utc"


def test_scan_folder_rejects_a_missing_folder(tmp_path: Path, highspec_map: HeaderMap):
    with pytest.raises(FileNotFoundError, match="nope"):
        scan_folder(tmp_path / "nope", "highspec", highspec_map)
    (tmp_path / "file.fits").write_bytes(b"")
    with pytest.raises(FileNotFoundError):
        scan_folder(tmp_path / "file.fits", "highspec", highspec_map)
    assert scan_folder(tmp_path, "highspec", highspec_map, unreadable=[]) == []


def test_unreadable_file_collected_not_raised(night: Path, highspec_map: HeaderMap):
    (night / "c_broken.fits").write_text("this is not a FITS file\n", encoding="utf-8")
    with pytest.raises(OSError):
        scan_folder(night, "highspec", highspec_map)
    unreadable: list[Path] = []
    refs = scan_folder(night, "highspec", highspec_map, unreadable=unreadable)
    assert unreadable == [night / "c_broken.fits"]
    assert [ref.meta.source for ref in refs] == ["a_bias.fits", "b_science.fits"]


def test_load_frame_crops_and_orients(tmp_path: Path):
    data = np.arange(24, dtype=np.int32).reshape(4, 6)
    write_fits(tmp_path / "f.fits", data, {})
    meta = FrameMeta(instrument="highspec", camera="main", exposure_type="flat", source="f.fits")
    ref = FrameRef(tmp_path / "f.fits", 0, meta)

    plain = load_frame(ref, make_detector())
    assert plain.data.dtype == np.float32
    assert plain.variance is None and plain.mask is None
    assert plain.meta is meta
    np.testing.assert_array_equal(plain.data, data)

    detector = make_detector(crop=(1, 5, 0, 3), transpose=True, flip_x=True)
    frame = load_frame(ref, detector)
    expected = np.array([[13, 7, 1], [14, 8, 2], [15, 9, 3], [16, 10, 4]], dtype=np.float32)
    np.testing.assert_array_equal(frame.data, expected)
    assert frame.data.dtype == np.float32
    # the loaded shape is a function of the detector alone
    assert expected_shape(detector) == frame.data.shape
    assert expected_shape(make_detector()) == plain.data.shape

    flipped = load_frame(ref, make_detector(flip_y=True))
    np.testing.assert_array_equal(flipped.data, data[::-1, :])


def test_expected_shape_scales_with_binning():
    """``shape`` and ``crop`` are unbinned pixels: at 1x1 the loaded shape is the crop's (or
    the frame's), oriented; at 2x2 each extent halves; x binning divides the raw x extent
    and y binning the raw y extent, before the transpose."""
    highspec = load_instrument("highspec").cameras["main"]
    assert expected_shape(highspec) == (398, 1598)
    assert expected_shape(highspec, (1, 1)) == (398, 1598)
    assert expected_shape(highspec, (2, 2)) == (199, 799)
    assert expected_shape(make_detector()) == (4, 6)
    assert expected_shape(make_detector(), (2, 2)) == (2, 3)
    rotated = make_detector(shape=(8, 12), crop=(2, 10, 0, 8), transpose=True)
    assert expected_shape(rotated) == (8, 8)
    assert expected_shape(rotated, (2, 1)) == (4, 8)  # x halves, then becomes the rows


def test_load_frame_of_a_binned_frame(tmp_path: Path):
    """A frame binned 2x2 (``HBIN = VBIN = 2``, half the raw size on each axis) is cropped
    at the halved bounds and loads to ``expected_shape(detector, (2, 2))``."""
    detector = load_instrument("highspec").cameras["main"]
    raw_y, raw_x = detector.shape
    data = np.full((raw_y // 2, raw_x // 2), 600, dtype=np.uint16)
    write_fits(tmp_path / "b.fits", data, {"HBIN": 2, "VBIN": 2})
    meta = FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type="bias",
        binning=(2, 2),
        source="b.fits",
    )
    frame = load_frame(FrameRef(tmp_path / "b.fits", 0, meta), detector)
    assert frame.data.shape == expected_shape(detector, (2, 2)) == (199, 799)


def test_manifest_overrides_fields(night: Path, highspec_map: HeaderMap, tmp_path: Path):
    path = tmp_path / "manifest.csv"
    write_manifest(
        path,
        [
            {
                "source": "b_science.fits",
                "exposure_type": "arc",
                "target": "",
                "arc_kind": "full",
                "exposure_time": "",
            },
            {
                "source": "a_bias.fits",
                "exposure_type": "",
                "target": "renamed",
                "arc_kind": "",
                "exposure_time": "0.5",
            },
        ],
    )
    manifest = read_manifest(path)
    assert list(manifest.columns) == [
        "source",
        "exposure_type",
        "target",
        "arc_kind",
        "exposure_time",
    ]

    refs = scan_folder(night, "highspec", highspec_map)
    new = apply_manifest(refs, manifest)
    assert len(new) == 2 and all(isinstance(ref, FrameRef) for ref in new)
    assert new[0].meta.exposure_type == "bias" and new[0].meta.classification_rule == "header"
    assert new[0].meta.targets == (("*", "renamed"),)
    assert new[0].meta.exposure_time == pytest.approx(0.5)
    assert new[1].meta.exposure_type == "arc" and new[1].meta.classification_rule == "manifest"
    assert new[1].meta.arc_kind == "full"
    assert new[1].meta.targets == (("*", "lab target"),)  # an empty cell is not an override
    assert new[1].meta.exposure_time == pytest.approx(120.0)
    assert new[1].meta.checksum == refs[1].meta.checksum
    assert refs[1].meta.exposure_type == "science"  # the input refs are untouched

    # the same through scan_folder, from the path
    again = scan_folder(night, "highspec", highspec_map, manifest=path)
    assert [ref.meta.exposure_type for ref in again] == ["bias", "arc"]


def test_manifest_rejects_bad_values_and_columns(night: Path, highspec_map: HeaderMap, tmp_path):
    refs = scan_folder(night, "highspec", highspec_map)
    bad_type = tmp_path / "bad_type.csv"
    write_manifest(bad_type, [{"source": "a_bias.fits", "exposure_type": "lamp"}])
    with pytest.raises(ValueError, match="exposure_type"):
        apply_manifest(refs, read_manifest(bad_type))
    bad_kind = tmp_path / "bad_kind.csv"
    write_manifest(bad_kind, [{"source": "a_bias.fits", "arc_kind": "partial"}])
    with pytest.raises(ValueError, match="arc_kind"):
        apply_manifest(refs, read_manifest(bad_kind))
    bad_column = tmp_path / "bad_column.csv"
    write_manifest(bad_column, [{"source": "a_bias.fits", "exposure": "1"}])
    with pytest.raises(ValueError, match="exposure"):
        read_manifest(bad_column)
    no_source = tmp_path / "no_source.csv"
    write_manifest(no_source, [{"file": "a_bias.fits", "target": "x"}])
    with pytest.raises(ValueError, match="source"):
        read_manifest(no_source)


def test_mef_extensions_become_separate_refs(tmp_path: Path):
    hmap = HeaderMap.from_dict(
        {
            "fields": {
                "exposure_time": {"keys": ["EXPOSURE"], "type": "float"},
                "camera": {"keys": ["EXTNAME"], "lower": True},
                "target": {"keys": ["OBJECT"]},
            },
            "classification": [{"type": "science"}],
        }
    )
    primary = fits.PrimaryHDU()
    primary.header["EXPOSURE"] = 30.0
    primary.header["OBJECT"] = "shared"
    blue = fits.ImageHDU(np.ones((2, 3), dtype=np.uint16), name="BLUE")
    red = fits.ImageHDU(np.ones((2, 3), dtype=np.uint16), name="RED")
    red.header["OBJECT"] = "red only"
    fits.HDUList([primary, blue, red]).writeto(tmp_path / "mef.fits")

    refs = scan_folder(tmp_path, "deepspec", hmap)
    assert [ref.extension for ref in refs] == [1, 2]
    assert [ref.meta.extension for ref in refs] == [1, 2]
    assert [ref.meta.camera for ref in refs] == ["blue", "red"]
    assert all(ref.meta.source == "mef.fits" for ref in refs)
    assert all(ref.meta.instrument == "deepspec" for ref in refs)
    assert [ref.meta.exposure_time for ref in refs] == [30.0, 30.0]  # inherited from primary
    assert [ref.meta.target_for("mast01") for ref in refs] == ["shared", "red only"]
    assert len({ref.meta.checksum for ref in refs}) == 1

    frame = load_frame(refs[1], make_detector(shape=(2, 3)))
    assert frame.data.shape == (2, 3)

    # a manifest row with an extension applies to that extension only
    path = tmp_path / "manifest.csv"
    write_manifest(path, [{"source": "mef.fits", "extension": "2", "exposure_type": "arc"}])
    new = apply_manifest(refs, read_manifest(path))
    assert [ref.meta.exposure_type for ref in new] == ["science", "arc"]
    write_manifest(path, [{"source": "mef.fits", "extension": "", "exposure_type": "flat"}])
    new = apply_manifest(refs, read_manifest(path))
    assert [ref.meta.exposure_type for ref in new] == ["flat", "flat"]


# ------------------------------------------- nested night folders, instrument gate, manifest


BLANK = np.zeros((2, 2), dtype=np.uint16)


@pytest.fixture
def nested(tmp_path: Path) -> Path:
    """A night filed the way the control software does it: one folder per acquisition,
    ``seq.txt`` counters in the parents, ``exposure-001.fits`` in every acquisition."""
    night = tmp_path / "night"
    for acquisition, kind in (("acquisition-0001", "Arc"), ("acquisition-0002", "Light")):
        folder = night / "highspec" / acquisition
        folder.mkdir(parents=True)
        cards = {"IMAGETYP": kind, "EXPOSURE": 60.0, "INSTRUME": "Highspec"}
        write_fits(folder / "exposure-001.fits", BLANK, cards)
    write_fits(
        night / "highspec" / "acquisition-0002" / "exposure-002.fits",
        BLANK,
        {"IMAGETYP": "Flat", "EXPOSURE": 5.0},
    )
    (night / "seq.txt").write_text("2\n", encoding="utf-8")
    (night / "highspec" / "seq.txt").write_text("2\n", encoding="utf-8")
    return night


def test_fibpre_card_is_ignored(tmp_path: Path, highspec_map: HeaderMap):
    """In MAST a fiber preset is a stage position, not a fiber layout: FIBPRE is not read."""
    write_fits(tmp_path / "frame.fits", BLANK, {"IMAGETYP": "Light", "FIBPRE": "survey45"})
    (ref,) = scan_folder(tmp_path, "highspec", highspec_map)
    assert ref.meta.fiber_preset is None


def test_scan_is_recursive_with_relative_sources(nested: Path, highspec_map: HeaderMap):
    refs = scan_folder(nested, "highspec", highspec_map)
    assert [ref.meta.source for ref in refs] == [
        "highspec/acquisition-0001/exposure-001.fits",
        "highspec/acquisition-0002/exposure-001.fits",
        "highspec/acquisition-0002/exposure-002.fits",
    ]
    assert refs[0].path == nested / "highspec" / "acquisition-0001" / "exposure-001.fits"
    assert [ref.meta.exposure_type for ref in refs] == ["arc", "light", "flat"]
    # a flat scan keeps to the folder itself
    assert scan_folder(nested, "highspec", highspec_map, recursive=False) == []


def test_fits_suffixes_are_matched_case_insensitively(tmp_path: Path, highspec_map: HeaderMap):
    for name in ("a.fits", "b.FITS", "c.fit", "d.Fit"):
        write_fits(tmp_path / name, BLANK, {"IMAGETYP": "Light"})
    fits.CompImageHDU(BLANK).writeto(tmp_path / "e.fits.fz")
    write_fits(tmp_path / "f.fits.gz.bak", BLANK, {"IMAGETYP": "Light"})
    sources = [ref.meta.source for ref in scan_folder(tmp_path, "highspec", highspec_map)]
    assert sources == ["a.fits", "b.FITS", "c.fit", "d.Fit", "e.fits.fz"]


def test_non_fits_and_hidden_files_are_not_scanned(nested: Path, highspec_map: HeaderMap):
    (nested / "plan.json").write_text("{}\n", encoding="utf-8")
    (nested / "highspec" / ".hidden.fits").write_text("not a FITS file\n", encoding="utf-8")
    hidden = nested / ".sync" / "acquisition-0009"
    hidden.mkdir(parents=True)
    (hidden / "exposure-001.fits").write_text("not a FITS file\n", encoding="utf-8")
    unreadable: list[Path] = []
    refs = scan_folder(nested, "highspec", highspec_map, unreadable=unreadable)
    assert unreadable == []
    assert len(refs) == 3
    assert not any("hidden" in ref.meta.source or ".sync" in ref.meta.source for ref in refs)


def test_instrument_gate_skips_other_instruments(nested: Path, highspec_map: HeaderMap):
    deepspec = nested / "deepspec" / "acquisition-0001"
    deepspec.mkdir(parents=True)
    write_fits(deepspec / "exposure-001.fits", BLANK, {"IMAGETYP": "Light", "INSTRUME": "Deepspec"})
    skipped: list[tuple[str, str]] = []
    refs = scan_folder(nested, "highspec", highspec_map, skipped=skipped)
    sources = [ref.meta.source for ref in refs]
    assert "highspec/acquisition-0001/exposure-001.fits" in sources  # INSTRUME = 'Highspec'
    assert "highspec/acquisition-0002/exposure-002.fits" in sources  # no INSTRUME
    assert len(refs) == 3
    assert skipped == [("deepspec/acquisition-0001/exposure-001.fits", "instrument 'deepspec' is not 'highspec'")]
    assert all("instrument" not in ref.meta.unresolved for ref in refs)
    # without a list the other instrument's frames are still left out
    assert len(scan_folder(nested, "highspec", highspec_map)) == 3


def test_manifest_matches_relative_path_glob_and_unique_name(nested: Path, highspec_map: HeaderMap, tmp_path: Path):
    refs = scan_folder(nested, "highspec", highspec_map)
    path = tmp_path / "manifest.csv"
    write_manifest(
        path,
        [
            {"source": "highspec/acquisition-0001/*.fits", "target": "", "arc_kind": "full"},
            {"source": "exposure-002.fits", "target": "flat lamp", "arc_kind": ""},
            {
                "source": "highspec/acquisition-0002/exposure-001.fits",
                "target": "star",
                "arc_kind": "",
            },
            {"source": "*.fits", "target": "top level only", "arc_kind": ""},
            {"source": "**/exposure-00[12].fits", "target": "", "arc_kind": ""},
        ],
    )
    unmatched: list[str] = []
    new = apply_manifest(refs, read_manifest(path), unmatched)
    assert [ref.meta.arc_kind for ref in new] == ["full", None, None]
    assert [ref.meta.target_for(None) for ref in new] == [None, "star", "flat lamp"]
    # "*" does not cross "/": the fourth row matches nothing at the top level
    assert unmatched == ["manifest row 3: *.fits"]


def test_bare_name_ambiguous_is_an_error(nested: Path, highspec_map: HeaderMap, tmp_path: Path):
    refs = scan_folder(nested, "highspec", highspec_map)
    path = tmp_path / "manifest.csv"
    write_manifest(path, [{"source": "exposure-001.fits", "exposure_type": "flat"}])
    with pytest.raises(ValueError, match="acquisition-0001.*acquisition-0002"):
        apply_manifest(refs, read_manifest(path))


def test_unmatched_manifest_rows_are_reported(night: Path, highspec_map: HeaderMap, tmp_path):
    refs = scan_folder(night, "highspec", highspec_map)
    path = tmp_path / "manifest.csv"
    write_manifest(path, [{"source": "nothing.fits", "target": "x"}])
    unmatched: list[str] = []
    apply_manifest(refs, read_manifest(path), unmatched)
    assert unmatched == ["manifest row 0: nothing.fits"]


def test_manifest_mid_time_is_normalised_to_utc(night: Path, highspec_map: HeaderMap, tmp_path):
    path = tmp_path / "manifest.csv"
    write_manifest(path, [{"source": "a_bias.fits", "mid_time": "2026-10-01T20:00:00.4+02:00"}])
    (ref, _) = scan_folder(night, "highspec", highspec_map, manifest=path)
    assert ref.meta.mid_time == "2026-10-01T18:00:00" and ref.meta.timezone == "utc"
    write_manifest(path, [{"source": "a_bias.fits", "mid_time": "yesterday-ish"}])
    with pytest.raises(ValueError, match="manifest row 0: bad mid_time"):
        scan_folder(night, "highspec", highspec_map, manifest=path)


# --------------------------------------------------------------------------- targets per unit


def test_per_unit_cards_become_targets(tmp_path: Path, highspec_map: HeaderMap):
    cards = {"EXPOSURE": 300.0, "OBJECT": "target_a+target_c", "OBJ01": "target_a"}
    cards["OBJ03"] = "target_c"
    write_fits(tmp_path / "science.fits", np.ones((2, 3), dtype=np.uint16), cards)
    (ref,) = scan_folder(tmp_path, "highspec", highspec_map)
    assert ref.meta.targets == (("mast01", "target_a"), ("mast03", "target_c"))
    assert ref.meta.target_for("mast03") == "target_c" and ref.meta.target_for("mast02") is None
    assert "targets" not in ref.meta.unresolved and "target" not in ref.meta.unresolved


def test_single_object_applies_to_all(tmp_path: Path, highspec_map: HeaderMap):
    cards = {"EXPOSURE": 300.0, "OBJECT": "target_a"}
    write_fits(tmp_path / "science.fits", np.ones((2, 3), dtype=np.uint16), cards)
    (ref,) = scan_folder(tmp_path, "highspec", highspec_map)
    assert ref.meta.targets == (("*", "target_a"),)
    assert ref.meta.target_for("mast17") == "target_a"


@pytest.mark.parametrize("blank", ["", "   "])
def test_empty_object_is_no_target(tmp_path: Path, highspec_map: HeaderMap, blank: str):
    """An empty or blank ``OBJECT`` card names no target: the frame's targets are
    unresolved, as without the card, rather than one target ``''`` for every fiber."""
    cards = {"EXPOSURE": 300.0, "OBJECT": blank}
    write_fits(tmp_path / "science.fits", np.ones((2, 3), dtype=np.uint16), cards)
    (ref,) = scan_folder(tmp_path, "highspec", highspec_map)
    assert ref.meta.targets == ()
    assert ref.meta.target_for("mast01") is None
    assert "targets" in ref.meta.unresolved


def test_manifest_targets_column(night: Path, highspec_map: HeaderMap, tmp_path: Path):
    path = tmp_path / "manifest.csv"
    write_manifest(
        path,
        [
            {"source": "b_science.fits", "targets": "mast01=HD 1; mast02 = HD 2;"},
            {"source": "a_bias.fits", "targets": "target_b"},
        ],
    )
    new = scan_folder(night, "highspec", highspec_map, manifest=path)
    assert new[1].meta.targets == (("mast01", "HD 1"), ("mast02", "HD 2"))
    assert new[0].meta.targets == (("*", "target_b"),)
    # the targets column wins over a legacy target column of the same row
    both = tmp_path / "both.csv"
    write_manifest(both, [{"source": "b_science.fits", "target": "old", "targets": "mast05=new"}])
    again = scan_folder(night, "highspec", highspec_map, manifest=both)
    assert again[1].meta.targets == (("mast05", "new"),)
    for bad in ("mast01=", "=HD 1", "mast01=a;mast01=b"):
        broken = tmp_path / "broken.csv"
        write_manifest(broken, [{"source": "b_science.fits", "targets": bad}])
        with pytest.raises(ValueError, match="targets"):
            apply_manifest(scan_folder(night, "highspec", highspec_map), read_manifest(broken))


def test_manifest_overrides_clear_unresolved(tmp_path: Path):
    """A field the manifest sets is no longer unresolved; the others stay listed. A
    ``target`` cell (the legacy column) resolves ``targets``."""
    meta = FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type="light",
        source="x.fits",
        unresolved=("mid_time", "targets", "arc_kind", "configuration"),
    )
    ref = FrameRef(path=tmp_path / "x.fits", extension=0, meta=meta)
    manifest = pd.DataFrame(
        [{"source": "x.fits", "mid_time": "2026-10-01T18:00:00", "target": "HD 1", "arc_kind": "full"}],
        dtype=str,
    )
    new = apply_manifest([ref], manifest)[0].meta
    assert new.mid_time == "2026-10-01T18:00:00"
    assert new.targets == (("*", "HD 1"),)
    assert new.arc_kind == "full"
    assert new.unresolved == ("configuration",)
    assert ref.meta.unresolved == ("mid_time", "targets", "arc_kind", "configuration")

    # a row that overrides nothing of the unresolved fields leaves the list as it is
    other = pd.DataFrame([{"source": "x.fits", "exposure_time": "3"}], dtype=str)
    assert apply_manifest([ref], other)[0].meta.unresolved == ref.meta.unresolved
