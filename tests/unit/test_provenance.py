"""Tests for the provenance cards that carry a science product's idempotence key."""

from astropy.io import fits

from mastpipe.pipeline.provenance import Provenance


def header_of(prov: Provenance) -> fits.Header:
    header = fits.Header()
    for key, value, comment in prov.to_header_cards():
        header[key] = (value, comment)
    return header


def test_provenance_round_trips_new_cards():
    prov = Provenance(
        pipeline_version="0.0.1",
        config_hash="abcdef012345",
        inputs=["0123456789abcdef"],
        calibrations=["bias_a"],
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
        layout="survey45",
        layout_hash="a1b2c3d4e5f6",
        masters_hash="0f1e2d3c4b5a",
        timezone="assumed-utc",
    )
    header = header_of(prov)
    assert header["MPLAYH"] == "a1b2c3d4e5f6"
    assert header["MPMKH"] == "0f1e2d3c4b5a"
    assert header["MPTZ"] == "assumed-utc"
    assert Provenance.from_header(header) == prov

    # unset values are written as "none" and read back as None
    bare = Provenance("0.0.1", "abcdef012345", [], [])
    header = header_of(bare)
    assert header["MPLAYH"] == "none" and header["MPMKH"] == "none"
    back = Provenance.from_header(header)
    assert back.layout_hash is None and back.masters_hash is None and back.timezone is None

    # a product written before these cards existed reads them as unset
    header = header_of(prov)
    for key in ("MPLAYH", "MPMKH", "MPTZ"):
        del header[key]
    back = Provenance.from_header(header)
    assert back.layout_hash is None and back.masters_hash is None and back.timezone is None
    assert back.layout == "survey45" and back.reference_hash == "0011223344aa"
