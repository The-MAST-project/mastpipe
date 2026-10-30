"""Tests for folder scanning, frame loading, manifests and checksums."""

import csv
import hashlib
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from mastpipe.config.models import DetectorConfig
from mastpipe.core.model import FrameMeta, FrameRef
from mastpipe.instruments.registry import load_header_map
from mastpipe.io.headers import HeaderMap
from mastpipe.io.ingest import (
    apply_manifest,
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
            "GRATING": "MgB518",
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
    assert bias.target is None and bias.arc_kind is None
    assert "target" in bias.unresolved and "exposure_time" not in bias.unresolved

    sci = refs[1].meta
    assert sci.exposure_type == "science" and sci.classification_rule == "rule:5"
    assert sci.exposure_time == pytest.approx(120.0)
    assert sci.binning == (2, 2)
    assert sci.detector_temperature == pytest.approx(-60.5)
    assert sci.readout_mode == "normal"
    assert sci.configuration == "mgb518"
    assert sci.target == "lab target"


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

    flipped = load_frame(ref, make_detector(flip_y=True))
    np.testing.assert_array_equal(flipped.data, data[::-1, :])


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
    assert new[0].meta.target == "renamed"
    assert new[0].meta.exposure_time == pytest.approx(0.5)
    assert new[1].meta.exposure_type == "arc" and new[1].meta.classification_rule == "manifest"
    assert new[1].meta.arc_kind == "full"
    assert new[1].meta.target == "lab target"  # an empty cell is not an override
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
    assert [ref.meta.target for ref in refs] == ["shared", "red only"]  # extension wins
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
