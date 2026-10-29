"""Layered configuration: YAML loading, deep merge, effective config and its hash.

The effective configuration of a run is built per (camera, configuration) by layering,
in order, the instrument's step parameters, the configuration block's partial
parameters and any command-line overrides; later layers win key by key. The result is
validated as :class:`~mastpipe.config.models.StepParams` and hashed for provenance.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import Any

import yaml

from mastpipe.config.models import EffectiveConfig, InstrumentConfig, StepParams


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Return a new mapping: override wins; nested mappings merge; lists and scalars replace.

    Neither input is modified and the result shares no containers with them.
    """
    out = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def load_yaml(path: Path | str) -> dict[str, Any]:
    """Read a YAML file whose top level is a mapping; an empty file is an empty mapping."""
    with Path(path).open(encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"{path}: top level must be a mapping, got {type(data).__name__}")
    return data


def _flatten(mapping: dict[str, Any], prefix: str = "") -> Iterator[tuple[str, Any]]:
    """Yield ``(dotted.key, leaf value)`` pairs of a nested mapping, in key order."""
    for key, value in mapping.items():
        path = f"{prefix}{key}"
        if isinstance(value, dict) and value:
            yield from _flatten(value, f"{path}.")
        else:
            yield path, value


def build_effective(
    instrument: InstrumentConfig,
    configuration: str,
    camera: str,
    overrides: Iterable[dict[str, Any]] = (),
) -> EffectiveConfig:
    """Layer instrument params -> configuration-block params -> overrides, in that order.

    Raises ``KeyError`` for an unknown configuration or camera and a Pydantic
    ``ValidationError`` when the merged parameters contain an unknown key or a bad value.
    """
    if configuration not in instrument.configurations:
        known = sorted(instrument.configurations)
        raise KeyError(f"unknown configuration {configuration!r} for {instrument.name}: {known}")
    if camera not in instrument.cameras:
        known = sorted(instrument.cameras)
        raise KeyError(f"unknown camera {camera!r} for {instrument.name}: {known}")
    block = instrument.configurations[configuration]
    layers = [f"instrument:{instrument.name}", f"configuration:{configuration}"]
    merged = deep_merge(instrument.params.model_dump(), block.params)
    for override in overrides:
        merged = deep_merge(merged, override)
        layers.extend(f"override:{key}={json.dumps(value)}" for key, value in _flatten(override))
    return EffectiveConfig(
        instrument=instrument.name,
        camera=camera,
        configuration=configuration,
        detector=instrument.cameras[camera],
        geometry=block.geometry,
        band_nm=block.band_nm,
        resolving_power=block.resolving_power,
        line_table=block.line_table,
        params=StepParams.model_validate(merged),
        layers=layers,
    )


def config_hash(effective: EffectiveConfig) -> str:
    """First 12 hex digits of the sha256 of the effective configuration's canonical JSON.

    The hash covers the values only (``layers``, which records how they were reached, is
    excluded), so the same effective configuration hashes the same on every platform and
    whether or not an override merely restates a value.
    """
    payload = effective.model_dump(mode="json", exclude={"layers"})
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:12]


def parse_override(text: str) -> dict[str, Any]:
    """Parse one ``key.path=value`` override into a nested mapping.

    The value is read as JSON when it parses (``1`` -> int, ``1.5`` -> float, ``true`` ->
    bool, ``null`` -> None, ``[1,2]`` -> list) and kept as the raw string otherwise.
    """
    key, separator, raw = text.partition("=")
    if not separator:
        raise ValueError(f"override {text!r} must have the form key.path=value")
    parts = key.split(".")
    if any(not part for part in parts):
        raise ValueError(f"override {text!r} has an empty key segment")
    try:
        value: Any = json.loads(raw)
    except json.JSONDecodeError:
        value = raw
    nested: dict[str, Any] = {parts[-1]: value}
    for part in reversed(parts[:-1]):
        nested = {part: nested}
    return nested
