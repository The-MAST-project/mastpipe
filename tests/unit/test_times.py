"""Tests for the one time helper: zones, offsets and the ``mid_time`` format."""

import pandas as pd
import pytest

from mastpipe.core.times import ASSUMED_UTC, check_timezone, is_explicit, read_time, to_utc, utc_iso, zone_label


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


@pytest.mark.parametrize("zone", ["UTC", "Etc/UTC"])
def test_utc_zone_names_are_the_utc_keyword(zone):
    """An IANA name for UTC means what ``utc`` means, so it is normalised to it here, the
    one place a zone is checked: frames read under either label the same way."""
    assert check_timezone(zone) == "utc"
    assert zone_label(check_timezone(zone), explicit=False) == "utc"


def test_zone_labels():
    assert zone_label("utc", explicit=False) == "utc"
    assert zone_label("assume-utc", explicit=False) == ASSUMED_UTC == "assumed-utc"
    assert zone_label("Asia/Jerusalem", explicit=False) == "Asia/Jerusalem"
    for zone in ("utc", "assume-utc", "Asia/Jerusalem"):
        assert zone_label(zone, explicit=True) == "utc"


def test_read_time_parses_once_and_says_whether_the_offset_was_given():
    stamp, explicit = read_time("2026-11-04T20:00:00+02:00", "assume-utc")
    assert explicit and utc_iso(stamp) == "2026-11-04T18:00:00"
    stamp, explicit = read_time("2026-09-09T11:35:24", "Asia/Jerusalem")
    assert not explicit and utc_iso(stamp) == "2026-09-09T08:35:24"
    with pytest.raises(ValueError):
        read_time("now")


def test_repeated_and_skipped_local_times_take_the_first_offset():
    """A wall-clock time that a daylight-saving change repeats (autumn) or skips (spring)
    is read with its first (fold 0) offset, never an error."""
    # 2026-10-25 03:00 CEST -> 02:00 CET: 02:30 happens twice, the first at +02:00
    assert utc_iso(to_utc("2026-10-25T02:30:00", "Europe/Berlin")) == "2026-10-25T00:30:00"
    # 2026-03-29 02:00 CET -> 03:00 CEST: 02:30 never happens, read at the earlier +01:00
    assert utc_iso(to_utc("2026-03-29T02:30:00", "Europe/Berlin")) == "2026-03-29T01:30:00"
