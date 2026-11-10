"""Header maps: translate raw header keywords into pipeline vocabulary and classify frames.

A header map is the ``header_map.yaml`` of an instrument, interpreted here. It has two
sections:

``fields``
    pipeline name -> ``{keys, type, lower, map, default, filename_regex, missing_values,
    reference, timezone, per_unit}``. ``keys`` are raw keywords tried in order (the first
    present wins); the value is then converted to ``type`` (``float``, ``int`` or
    ``str``), lower-cased when ``lower`` is true and translated through ``map`` (an
    unmapped value is left as is). A keyword whose raw or converted value is listed in
    ``missing_values`` (a sentinel such as the ``-999`` of an unread sensor) counts as
    absent and the next keyword is tried.
    When no keyword is present ``filename_regex`` is searched in the file name (its
    ``value`` group is the value, transformed the same way), and ``default`` is used when
    that fails too. A field that still has no value is *unresolved*: ``None`` in the result
    and listed by name. A value that cannot be converted counts as unresolved as well.
    ``reference`` (``start``, ``mid`` or ``end``; default ``start``) says which instant of
    the exposure a time field records; ingest derives the mid-time from it. ``timezone``
    (``utc``, ``assume-utc`` or an IANA zone name; default ``utc``) says how a time
    field's stamp without an offset is read (see :mod:`mastpipe.core.times`); ``assume-utc``
    reads it as UTC but records that this is an assumption.

    A field may instead declare ``per_unit``, a keyword template containing ``{nn}``, in
    place of ``keys``, ``filename_regex`` and ``default``: the template is expanded with
    each telescope unit's two-digit number, and the value is the tuple of ``(unit name,
    value)`` pairs of the units whose keyword is present (each value transformed as above;
    a blank or ``missing_values`` value counts as absent), in unit order, or unresolved
    when no unit's keyword is present. The fields of :data:`PER_UNIT_FIELDS` (``targets``)
    are read this way only: ingest takes their value as ``(unit, value)`` pairs.

``classification``
    an ordered list of rules ``{when: {<field>: {lt: v} | {eq: v} | {missing: true}}, type}``
    applied when the header carries no valid ``exposure_type``. Every condition of a rule
    must hold; a rule without ``when`` matches unconditionally and serves as the fallback.

A valid header type wins, except ``light``: it is provisional ("shutter open, source not
stated", what the control software writes for arcs and science frames alike), so the rules
**with conditions** are tried on it and the first that matches wins (a ``light`` frame with
a zero exposure time is a bias). With no such match it stays ``light``, for the runner to
resolve from the frame's content; the unconditional fallback applies only to a frame
without a valid header type.

Raw keyword names live only in the YAML files; this module never names one.
"""

from __future__ import annotations

import contextlib
import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from astropy.io.fits.card import Undefined

from mastpipe.core.model import EXPOSURE_TYPES
from mastpipe.core.times import UTC, check_timezone
from mastpipe.vocab import UNIT_COUNT, unit_name

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
    "per_unit",
)
OPERATORS = ("lt", "eq", "missing")
TIME_REFERENCES = ("start", "mid", "end")
"""Which instant of the exposure a time field records (the ``reference`` option)."""
UNIT_PLACEHOLDER = "{nn}"
"""Where a ``per_unit`` template takes the two-digit unit number."""
PER_UNIT_FIELDS = ("targets",)
"""Pipeline fields whose value is one entry per telescope unit: they must declare
``per_unit`` (``target`` is the single-value field for every fiber)."""
UNCLASSIFIED = "unknown"
PROVISIONAL = "light"
"""The header type that the rules with conditions may still refine."""

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
    per_unit: str | None = None

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
            with contextlib.suppress(TypeError):  # unhashable value: nothing to map
                value = self.map.get(value, value)
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
        a field of :data:`PER_UNIT_FIELDS` without ``per_unit``, a rule referring to a field
        the map does not define, or an invalid exposure type."""
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
        the first matching rule (``"rule:<i>"``), else ``"unknown"`` (``"unclassified"``).
        A ``light`` header type is refined by the first matching rule that has conditions
        and stays ``("light", "header")`` when none matches."""
        from_header = fields.get("exposure_type")
        provisional = from_header == PROVISIONAL
        if from_header in EXPOSURE_TYPES and not provisional:
            return from_header, "header"
        for i, rule in enumerate(self.rules):
            if provisional and not rule.conditions:
                continue
            if rule.matches(fields):
                return rule.exposure_type, f"rule:{i}"
        if provisional:
            return PROVISIONAL, "header"
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
    per_unit = spec.get("per_unit")
    if per_unit is not None:
        if not isinstance(per_unit, str) or UNIT_PLACEHOLDER not in per_unit:
            raise ValueError(
                f"header map field {name!r}: per_unit must be a keyword template "
                f"containing {UNIT_PLACEHOLDER}, got {per_unit!r}"
            )
        clashing = [option for option in ("keys", "filename_regex", "default") if option in spec]
        if clashing:
            raise ValueError(f"header map field {name!r}: per_unit excludes {', '.join(clashing)}")
    elif name in PER_UNIT_FIELDS:
        raise ValueError(
            f"header map field {name!r}: must be read per_unit (a keyword template with "
            f"{UNIT_PLACEHOLDER}); a single value for every fiber is the 'target' field"
        )
    elif not keys and regex is None and "default" not in spec:
        raise ValueError(f"header map field {name!r}: needs keys, filename_regex or default")
    mapping = spec.get("map")
    if mapping is not None and not isinstance(mapping, Mapping):
        raise ValueError(f"header map field {name!r}: map must be a mapping")
    missing = spec.get("missing_values") or ()
    if isinstance(missing, str) or not isinstance(missing, list | tuple):
        raise ValueError(f"header map field {name!r}: missing_values must be a list")
    reference = spec.get("reference", "start")
    if reference not in TIME_REFERENCES:
        raise ValueError(f"header map field {name!r}: reference must be one of {TIME_REFERENCES}, got {reference!r}")
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
        default=spec.get("default", _NO_DEFAULT),
        filename_regex=regex,
        missing_values=tuple(missing),
        reference=str(reference),
        timezone=timezone,
        per_unit=per_unit,
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
            raise ValueError(f"{label}: condition on {name!r} must be a mapping of operator -> value")
        for op, value in condition.items():
            if op not in OPERATORS:
                raise ValueError(f"{label}: operator must be one of {OPERATORS}, got {op!r}")
            conditions.append((name, op, value))
    return Rule(conditions=tuple(conditions), exposure_type=exposure_type)


def _resolve_field(spec: FieldSpec, header: Mapping[str, Any], filename: str) -> Any:
    if spec.per_unit is not None:
        return _resolve_per_unit(spec, spec.per_unit, header)
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


def _resolve_per_unit(spec: FieldSpec, template: str, header: Mapping[str, Any]) -> tuple[tuple[str, Any], ...] | None:
    """``(unit name, value)`` for every unit whose expanded ``template`` keyword is present;
    ``None`` when none is."""
    pairs: list[tuple[str, Any]] = []
    for number in range(1, UNIT_COUNT + 1):
        key = template.replace(UNIT_PLACEHOLDER, f"{number:02d}")
        if key not in header:
            continue
        value = header[key]
        if not _present(value) or spec.is_missing(value):
            continue
        if isinstance(value, str) and not value.strip():
            continue
        converted = _convert(spec, value)
        if converted is None or spec.is_missing(converted):
            continue
        pairs.append((unit_name(number), converted))
    return tuple(pairs) or None


def _convert(spec: FieldSpec, value: Any) -> Any:
    try:
        return spec.transform(value)
    except ValueError:
        return None
