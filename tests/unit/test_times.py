"""Tests for the one time helper: zones, offsets and the ``mid_time`` format."""

import pandas as pd
import pytest

from mastpipe.core.times import is_explicit, to_utc, utc_iso


def test_naive_is_read_in_the_named_zone():
    assert utc_iso(to_utc("2026-09-09T11:35:24", "Asia/Jerusalem")) == "2026-09-09T08:35:24"


def test_assume_utc_keeps_the_clock():
    assert utc_iso(to_utc("2026-09-09T11:35:24", "assume-utc")) == "2026-09-09T11:35:24"
    assert utc_iso(to_utc("2026-09-09T11:35:24")) == "2026-09-09T11:35:24"


def test_explicit_offset_is_honoured_under_assume_utc():
    assert utc_iso(to_utc("2026-11-04T20:00:00+02:00", "assume-utc")) == "2026-11-04T18:00:00"
    assert utc_iso(to_utc("2026-11-04T20:00:00+02:00", "Asia/Tokyo")) == "2026-11-04T18:00:00"
    assert is_explicit("2026-11-04T20:00:00+02:00") and not is_explicit("2026-11-04T20:00:00")
    assert is_explicit("2026-11-04T20:00:00Z")
    assert not is_explicit("not a time")


def test_unknown_zone_is_an_error():
    with pytest.raises(ValueError):
        to_utc("2026-01-01T00:00:00", "Mars/Olympus")


def test_unparsable_value_is_an_error():
    for bad in ("", "now", "yesterday at noon"):
        with pytest.raises(ValueError):
            to_utc(bad)


def test_result_is_aware_utc_and_the_iso_form_drops_fractions():
    stamp = to_utc("2026-09-09T11:35:24.567", "assume-utc")
    assert stamp.tzinfo is not None and stamp.utcoffset() == pd.Timedelta(0)
    assert stamp.microsecond == 567000
    assert utc_iso(stamp) == "2026-09-09T11:35:24"
    assert utc_iso(pd.Timestamp("2026-09-09T11:35:24")) == "2026-09-09T11:35:24"
