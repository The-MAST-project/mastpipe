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
    the value is a wall-clock time in that zone. A name of UTC itself (``"UTC"``,
    ``"Etc/UTC"`` and their aliases, :data:`UTC_NAMES`) is the ``"utc"`` keyword:
    :func:`check_timezone` normalises it, so both spellings read and label a time alike.

:func:`read_time` parses a value once and says whether it carried its own offset;
:func:`zone_label` turns that into the ``FrameMeta.timezone`` label: ``"utc"`` for an
explicit value or under ``"utc"``, ``"assumed-utc"`` (:data:`ASSUMED_UTC`) under
``"assume-utc"``, else the zone name.

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
ASSUMED_UTC = "assumed-utc"
"""The ``FrameMeta.timezone`` label of a stamp without an offset taken as UTC under
``assume-utc``."""
UTC_NAMES = frozenset({"UTC", "Etc/UTC", "Etc/UCT", "UCT", "Etc/Universal", "Universal", "Etc/Zulu", "Zulu"})
"""IANA names of UTC itself, read as the ``"utc"`` keyword."""


def check_timezone(timezone: Any) -> str:
    """``timezone`` when it is ``"utc"``, ``"assume-utc"`` or a known IANA zone name, a name
    of UTC itself (:data:`UTC_NAMES`) as ``"utc"``; ``ValueError`` otherwise."""
    if timezone in KEYWORDS:
        return timezone
    if not isinstance(timezone, str) or not timezone:
        raise ValueError(f"timezone must be one of {KEYWORDS} or an IANA zone, got {timezone!r}")
    if timezone in UTC_NAMES:
        return UTC
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


def read_time(value: str, timezone: str = UTC) -> tuple[pd.Timestamp, bool]:
    """``(stamp, explicit)``: ``value`` parsed once as an aware UTC timestamp, a naive value
    read according to ``timezone`` (see the module docstring), and whether it carried its
    own offset. ``ValueError`` for an unparsable value or an unknown zone."""
    timezone = check_timezone(timezone)
    parsed = _parse(value)
    explicit = parsed.tzinfo is not None
    if not explicit:
        # replace() never raises: a wall-clock time repeated or skipped by a daylight-saving
        # change is read with its first (fold 0) offset
        parsed = parsed.replace(tzinfo=ZoneInfo("UTC" if timezone in KEYWORDS else timezone))
    return pd.Timestamp(parsed).tz_convert("UTC"), explicit


def to_utc(value: str, timezone: str = UTC) -> pd.Timestamp:
    """``value`` as an aware UTC timestamp, a naive value read according to ``timezone``
    (see the module docstring). ``ValueError`` for an unparsable value or an unknown zone."""
    return read_time(value, timezone)[0]


def zone_label(timezone: str, explicit: bool) -> str:
    """The ``FrameMeta.timezone`` label of a time read under ``timezone`` (a value checked
    by :func:`check_timezone`): ``"utc"`` when the value carried its own offset or under
    ``"utc"``, :data:`ASSUMED_UTC` under ``"assume-utc"``, else the zone name."""
    if explicit:
        return UTC
    return ASSUMED_UTC if timezone == ASSUME_UTC else timezone


def utc_iso(stamp: pd.Timestamp | datetime) -> str:
    """``stamp`` as naive ISO-8601 UTC to the second (``2026-10-01T18:02:30``); a naive
    ``stamp`` is taken to be UTC."""
    stamp = pd.Timestamp(stamp)
    if stamp.tzinfo is not None:
        stamp = stamp.tz_convert("UTC").tz_localize(None)
    return stamp.isoformat(timespec="seconds")
