"""Header maps: translate raw header keywords into pipeline vocabulary and classify frames.

A header map is the ``header_map.yaml`` of an instrument, interpreted here. It has two
sections:

``fields``
    pipeline name -> ``{keys, type, lower, map, default, filename_regex, missing_values,
    reference, timezone}``. ``keys`` are raw keywords tried in order (the first present wins); the
    value is then converted to ``type`` (``float``, ``int`` or ``str``), lower-cased when
    ``lower`` is true and translated through ``map`` (an unmapped value is left as is). A
    keyword whose raw or converted value is listed in ``missing_values`` (a sentinel such
    as the ``-999`` of an unread sensor) counts as absent and the next keyword is tried.
    When no keyword is present ``filename_regex`` is searched in the file name (its
    ``value`` group is the value, transformed the same way), and ``default`` is used when
    that fails too. A field that still has no value is *unresolved*: ``None`` in the result
    and listed by name. A value that cannot be converted counts as unresolved as well.
    ``reference`` (``start``, ``mid`` or ``end``; default ``start``) says which instant of
    the exposure a time field records; ingest derives the mid-time from it. ``timezone``
    (``utc``, ``assume-utc`` or an IANA zone name; default ``utc``) says how a time
    field's stamp without an offset is read (see :mod:`mastpipe.core.times`); ``assume-utc``
    reads it as UTC but records that this is an assumption.

``classification``
    an ordered list of rules ``{when: {<field>: {lt: v} | {eq: v} | {missing: true}}, type}``
    applied when the header carries no valid ``exposure_type``. Every condition of a rule
    must hold; a rule without ``when`` matches unconditionally and serves as the fallback.

Raw keyword names live only in the YAML files; this module never names one.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from astropy.io.fits.card import Undefined

from mastpipe.core.model import EXPOSURE_TYPES
from mastpipe.core.times import UTC, check_timezone

TYPES: dict[str, type] = {"float": float, "int": int, "str": str}
FIELD_OPTIONS = (
    "keys",
    "type",
    "lower",
    "map",
    "default",
    "filename_regex",
    "missing_values",
    "reference",
    "timezone",
)
OPERATORS = ("lt", "eq", "missing")
TIME_REFERENCES = ("start", "mid", "end")
"""Which instant of the exposure a time field records (the ``reference`` option)."""
UNCLASSIFIED = "unknown"

_NO_DEFAULT = object()


@dataclass(frozen=True)
class FieldSpec:
    """How one pipeline field is read from a header or a file name."""

    name: str
    keys: tuple[str, ...] = ()
    type: str | None = None
    lower: bool = False
    map: Mapping[Any, Any] | None = None
    default: Any = _NO_DEFAULT
    filename_regex: re.Pattern[str] | None = None
    missing_values: tuple[Any, ...] = ()
    reference: str = "start"
    timezone: str = UTC

    @property
    def has_default(self) -> bool:
        return self.default is not _NO_DEFAULT

    def is_missing(self, value: Any) -> bool:
        """Whether ``value`` is one of the field's ``missing_values`` sentinels."""
        try:
            return any(value == sentinel for sentinel in self.missing_values)
        except (TypeError, ValueError):  # incomparable value: not a sentinel
            return False

    def transform(self, value: Any) -> Any:
        """Apply ``type``, ``lower`` and ``map`` in that order; ``ValueError`` when the
        conversion fails."""
        if self.type is not None:
            try:
                value = TYPES[self.type](value)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{self.name}: cannot convert {value!r} to {self.type}") from error
        if self.lower and isinstance(value, str):
            value = value.lower()
        if self.map is not None:
            try:
                value = self.map.get(value, value)
            except TypeError:  # unhashable value: nothing to map
                pass
        return value


@dataclass(frozen=True)
class Rule:
    """One classification rule: ``conditions`` are ``(field, operator, value)`` triples."""

    conditions: tuple[tuple[str, str, Any], ...]
    exposure_type: str

    def matches(self, fields: Mapping[str, Any]) -> bool:
        return all(_holds(fields.get(name), op, value) for name, op, value in self.conditions)


def _numeric(value: Any) -> bool:
    return isinstance(value, int | float) and not isinstance(value, bool)


def _holds(actual: Any, op: str, expected: Any) -> bool:
    if op == "missing":
        return (actual is None) == bool(expected)
    if actual is None:
        return False
    if op == "eq":
        return actual == expected
    # "lt": an untyped string (or anything else that is not a number) never matches
    if not _numeric(actual) or not _numeric(expected):
        return False
    return actual < expected


def _present(value: Any) -> bool:
    """A header value counts as present unless it is ``None`` or a FITS undefined card."""
    return value is not None and not isinstance(value, Undefined)


@dataclass(frozen=True)
class HeaderMap:
    """The parsed header map of an instrument: field specs and classification rules."""

    fields: Mapping[str, FieldSpec]
    rules: tuple[Rule, ...] = ()

    @classmethod
    def from_dict(cls, d: Mapping[str, Any]) -> HeaderMap:
        """Parse the raw YAML mapping; ``ValueError`` for an unknown option, type or operator,
        a rule referring to a field the map does not define, or an invalid exposure type."""
        raw_fields = d.get("fields") or {}
        if not isinstance(raw_fields, Mapping):
            raise ValueError("header map: 'fields' must be a mapping")
        fields = {name: _parse_field(name, spec) for name, spec in raw_fields.items()}
        raw_rules = d.get("classification") or []
        if not isinstance(raw_rules, list):
            raise ValueError("header map: 'classification' must be a list of rules")
        rules = tuple(_parse_rule(i, rule, fields) for i, rule in enumerate(raw_rules))
        return cls(fields=fields, rules=rules)

    def resolve(self, header: Mapping[str, Any], filename: str) -> tuple[dict[str, Any], list[str]]:
        """Every field of the map read from ``header`` (any mapping, an astropy ``Header``
        included) with ``filename`` as the regex fallback; unresolved fields are ``None`` and
        listed by name in map order."""
        values: dict[str, Any] = {}
        unresolved: list[str] = []
        for name, spec in self.fields.items():
            value = _resolve_field(spec, header, filename)
            values[name] = value
            if value is None:
                unresolved.append(name)
        return values, unresolved

    def classify(self, fields: Mapping[str, Any]) -> tuple[str, str]:
        """``(exposure_type, rule name)``: the header's own valid type (``"header"``), else
        the first matching rule (``"rule:<i>"``), else ``"unknown"`` (``"unclassified"``)."""
        from_header = fields.get("exposure_type")
        if from_header in EXPOSURE_TYPES:
            return from_header, "header"
        for i, rule in enumerate(self.rules):
            if rule.matches(fields):
                return rule.exposure_type, f"rule:{i}"
        return UNCLASSIFIED, "unclassified"


def _parse_field(name: str, spec: Any) -> FieldSpec:
    if not isinstance(spec, Mapping):
        raise ValueError(f"header map field {name!r}: must be a mapping")
    unknown = sorted(set(spec) - set(FIELD_OPTIONS))
    if unknown:
        raise ValueError(f"header map field {name!r}: unknown options {unknown}")
    keys = spec.get("keys") or ()
    if isinstance(keys, str):
        keys = (keys,)
    keys = tuple(str(key) for key in keys)
    type_name = spec.get("type")
    if type_name is not None and type_name not in TYPES:
        raise ValueError(f"header map field {name!r}: type must be one of {sorted(TYPES)}")
    regex = None
    if spec.get("filename_regex") is not None:
        regex = re.compile(spec["filename_regex"])
        if "value" not in regex.groupindex:
            raise ValueError(f"header map field {name!r}: filename_regex needs a 'value' group")
    if not keys and regex is None and "default" not in spec:
        raise ValueError(f"header map field {name!r}: needs keys, filename_regex or default")
    mapping = spec.get("map")
    if mapping is not None and not isinstance(mapping, Mapping):
        raise ValueError(f"header map field {name!r}: map must be a mapping")
    missing = spec.get("missing_values") or ()
    if isinstance(missing, str) or not isinstance(missing, list | tuple):
        raise ValueError(f"header map field {name!r}: missing_values must be a list")
    reference = spec.get("reference", "start")
    if reference not in TIME_REFERENCES:
        raise ValueError(
            f"header map field {name!r}: reference must be one of {TIME_REFERENCES}, "
            f"got {reference!r}"
        )
    try:
        timezone = check_timezone(spec.get("timezone", UTC))
    except ValueError as error:
        raise ValueError(f"header map field {name!r}: timezone: {error}") from error
    return FieldSpec(
        name=name,
        keys=keys,
        type=type_name,
        lower=bool(spec.get("lower", False)),
        map=dict(mapping) if mapping is not None else None,
        default=spec["default"] if "default" in spec else _NO_DEFAULT,
        filename_regex=regex,
        missing_values=tuple(missing),
        reference=str(reference),
        timezone=timezone,
    )


def _parse_rule(index: int, rule: Any, fields: Mapping[str, FieldSpec]) -> Rule:
    label = f"header map rule {index}"
    if not isinstance(rule, Mapping) or "type" not in rule:
        raise ValueError(f"{label}: must be a mapping with a 'type'")
    unknown = sorted(set(rule) - {"when", "type"})
    if unknown:
        raise ValueError(f"{label}: unknown keys {unknown}")
    exposure_type = rule["type"]
    if exposure_type not in EXPOSURE_TYPES:
        raise ValueError(f"{label}: type must be one of {EXPOSURE_TYPES}, got {exposure_type!r}")
    when = rule.get("when") or {}
    if not isinstance(when, Mapping):
        raise ValueError(f"{label}: 'when' must be a mapping of field -> condition")
    conditions: list[tuple[str, str, Any]] = []
    for name, condition in when.items():
        if name not in fields:
            raise ValueError(f"{label}: unknown field {name!r}")
        if not isinstance(condition, Mapping) or not condition:
            raise ValueError(
                f"{label}: condition on {name!r} must be a mapping of operator -> value"
            )
        for op, value in condition.items():
            if op not in OPERATORS:
                raise ValueError(f"{label}: operator must be one of {OPERATORS}, got {op!r}")
            conditions.append((name, op, value))
    return Rule(conditions=tuple(conditions), exposure_type=exposure_type)


def _resolve_field(spec: FieldSpec, header: Mapping[str, Any], filename: str) -> Any:
    for key in spec.keys:
        if key in header:
            value = header[key]
            if not _present(value) or spec.is_missing(value):
                continue
            converted = _convert(spec, value)
            if converted is not None and spec.is_missing(converted):
                continue
            return converted
    if spec.filename_regex is not None:
        match = spec.filename_regex.search(filename)
        if match is not None:
            return _convert(spec, match.group("value"))
    if spec.has_default:
        return spec.default
    return None


def _convert(spec: FieldSpec, value: Any) -> Any:
    try:
        return spec.transform(value)
    except ValueError:
        return None
