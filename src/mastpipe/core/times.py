"""Times: the one place where a time string is parsed and turned into UTC.

A time is ISO-8601 (``2026-10-01T18:02:30``, optionally with fractional seconds and an
offset such as ``+02:00`` or ``Z``). A value that carries its own offset is *explicit* and
is converted to UTC whatever else is said. A *naive* value (no offset) is read according
to ``timezone``:

``"utc"``
    the value is UTC (the pipeline's own strings: every ``mid_time`` it writes is UTC);
``"assume-utc"``
    the value is taken as UTC although nobody has established that its clock is (a header
    map's way of saying so out loud);
an IANA zone name (``"Asia/Jerusalem"``)
    the value is a wall-clock time in that zone.

:func:`utc_iso` writes the form every ``mid_time`` uses: naive ISO-8601 UTC, whole
seconds (fractions are dropped). Nothing here reads a file or a configuration.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

import pandas as pd

UTC = "utc"
ASSUME_UTC = "assume-utc"
KEYWORDS = (UTC, ASSUME_UTC)
"""The ``timezone`` values that are not IANA zone names."""


def check_timezone(timezone: Any) -> str:
    """``timezone`` when it is ``"utc"``, ``"assume-utc"`` or a known IANA zone name;
    ``ValueError`` otherwise."""
    if timezone in KEYWORDS:
        return timezone
    if not isinstance(timezone, str) or not timezone:
        raise ValueError(f"timezone must be one of {KEYWORDS} or an IANA zone, got {timezone!r}")
    try:
        ZoneInfo(timezone)
    except (ZoneInfoNotFoundError, ValueError) as error:
        raise ValueError(f"unknown timezone {timezone!r}") from error
    return timezone


def _parse(value: str) -> datetime:
    try:
        return datetime.fromisoformat(str(value).strip())
    except ValueError as error:
        raise ValueError(f"not an ISO-8601 time: {value!r}") from error


def is_explicit(value: str) -> bool:
    """Whether ``value`` carries its own UTC offset (false for an unparsable value)."""
    try:
        return _parse(value).tzinfo is not None
    except ValueError:
        return False


def to_utc(value: str, timezone: str = UTC) -> pd.Timestamp:
    """``value`` as an aware UTC timestamp, a naive value read according to ``timezone``
    (see the module docstring). ``ValueError`` for an unparsable value or an unknown zone."""
    timezone = check_timezone(timezone)
    parsed = _parse(value)
    if parsed.tzinfo is None:
        # replace() never raises: a wall-clock time repeated or skipped by a daylight-saving
        # change is read with its first (fold 0) offset
        parsed = parsed.replace(tzinfo=ZoneInfo("UTC" if timezone in KEYWORDS else timezone))
    return pd.Timestamp(parsed).tz_convert("UTC")


def utc_iso(stamp: pd.Timestamp | datetime) -> str:
    """``stamp`` as naive ISO-8601 UTC to the second (``2026-10-01T18:02:30``); a naive
    ``stamp`` is taken to be UTC."""
    stamp = pd.Timestamp(stamp)
    if stamp.tzinfo is not None:
        stamp = stamp.tz_convert("UTC").tz_localize(None)
    return stamp.isoformat(timespec="seconds")
