"""Locate and load instrument definitions shipped as package data.

An instrument is a folder next to this module. ``config.yaml`` carries everything except the
configuration blocks, which are one YAML file each under ``gratings/`` (HighSpec) or
``bands/`` (DeepSpec) and are merged into ``configurations`` here, so adding a grating is a
data change with no code edit. Layouts are ``layouts/<preset>.csv`` for every preset named
in ``config.yaml``; line tables are referenced by the configuration block; anchor pairs are
``anchors/<configuration>_<preset>.yaml`` and are optional.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import pandas as pd

from mastpipe.config.loader import load_yaml
from mastpipe.config.models import ConfigurationBlock, InstrumentConfig
from mastpipe.core.model import ROLES
from mastpipe.vocab import GRATING_NAMES

LAYOUT_COLUMNS = ("trace", "role", "telescope", "y_nominal", "slope_nominal", "aperture")
LINE_TABLE_COLUMNS = ("wavelength_nm", "species", "logI", "usable")
BLOCK_FOLDERS = ("gratings", "bands")
LAYOUT_HASH_LENGTH = 12


def _root() -> Path:
    return Path(__file__).resolve().parent


def list_instruments() -> list[str]:
    """Names of the instruments shipped with the package, sorted."""
    return sorted(path.name for path in _root().iterdir() if (path / "config.yaml").is_file())


def instrument_dir(name: str) -> Path:
    """Folder of the instrument ``name``; ``KeyError`` when no such instrument is shipped."""
    path = _root() / name
    if not (path / "config.yaml").is_file():
        raise KeyError(f"unknown instrument {name!r}; known: {list_instruments()}")
    return path


def load_instrument(name: str) -> InstrumentConfig:
    """Validate ``config.yaml`` with every ``gratings/*.yaml`` and ``bands/*.yaml`` merged in.

    Each block file is validated on its own as a :class:`ConfigurationBlock`; its ``name``
    must equal the file stem, so that the configuration name, the file and the line table
    it references cannot drift apart. A grating block's ``name`` must be one of
    :data:`~mastpipe.vocab.GRATING_NAMES`, the names the control software writes.
    """
    folder = instrument_dir(name)
    raw = load_yaml(folder / "config.yaml")
    blocks: dict[str, Any] = dict(raw.get("configurations") or {})
    for sub in BLOCK_FOLDERS:
        for path in sorted((folder / sub).glob("*.yaml")):
            block = ConfigurationBlock.model_validate(load_yaml(path))
            if block.name != path.stem:
                raise ValueError(f"{path.name}: block name {block.name!r} must equal the file stem")
            if block.name in blocks:
                raise ValueError(f"{path.name}: configuration {block.name!r} defined twice")
            if sub == "gratings" and block.name not in GRATING_NAMES:
                raise ValueError(f"grating block {block.name!r} is not a MAST_common grating name {GRATING_NAMES}")
            blocks[block.name] = block.model_dump()
    raw["configurations"] = blocks
    return InstrumentConfig.model_validate(raw)


def load_header_map(name: str) -> dict[str, Any]:
    """The raw ``header_map.yaml`` mapping; :mod:`mastpipe.io.headers` interprets it."""
    return load_yaml(instrument_dir(name) / "header_map.yaml")


def load_layout(name: str, preset: str) -> pd.DataFrame:
    """The fiber layout ``preset`` as a table with :data:`LAYOUT_COLUMNS`, sorted by ``y_nominal``.

    ``trace`` is int, ``telescope`` is a string or ``None`` (arc fibers have none), the
    remaining columns are float (``aperture`` is NaN when the layout leaves it to the step
    parameters). Roles must be in :data:`~mastpipe.core.model.ROLES` and trace ids unique.
    """
    instrument = load_instrument(name)
    if preset not in instrument.layouts:
        raise KeyError(f"unknown layout preset {preset!r} for {name}; known: {instrument.layouts}")
    path = instrument_dir(name) / "layouts" / f"{preset}.csv"
    table = pd.read_csv(path, dtype={"telescope": "string"}, encoding="utf-8")
    missing = [column for column in LAYOUT_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"{path.name}: missing layout columns {missing}")
    table = table.loc[:, list(LAYOUT_COLUMNS)].copy()
    table["trace"] = table["trace"].astype(int)
    table["telescope"] = table["telescope"].astype(object).where(table["telescope"].notna(), None)
    for column in ("y_nominal", "slope_nominal", "aperture"):
        table[column] = table[column].astype(float)
    bad_roles = set(table["role"]) - set(ROLES)
    if bad_roles:
        raise ValueError(f"{path.name}: unknown roles {sorted(bad_roles)}; allowed: {ROLES}")
    if table["trace"].duplicated().any():
        raise ValueError(f"{path.name}: duplicate trace ids")
    return table.sort_values("y_nominal", kind="stable").reset_index(drop=True)


def layout_hash(name: str, preset: str) -> str:
    """The first :data:`LAYOUT_HASH_LENGTH` hex characters of the SHA-256 of the layout
    ``preset``'s CSV file bytes with CRLF line endings read as LF (part of a science
    product's idempotence key: an edited layout reduces the night again, a checkout with
    other line endings does not). ``KeyError`` for an unknown preset."""
    instrument = load_instrument(name)
    if preset not in instrument.layouts:
        raise KeyError(f"unknown layout preset {preset!r} for {name}; known: {instrument.layouts}")
    data = (instrument_dir(name) / "layouts" / f"{preset}.csv").read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()[:LAYOUT_HASH_LENGTH]


def load_line_table(name: str, configuration: str) -> pd.DataFrame:
    """The catalogue line table of ``configuration`` with :data:`LINE_TABLE_COLUMNS`.

    Sorted by wavelength (air nm); ``usable`` is bool. ``KeyError`` for an unknown
    configuration, ``ValueError`` when the block declares no line table.
    """
    instrument = load_instrument(name)
    if configuration not in instrument.configurations:
        known = sorted(instrument.configurations)
        raise KeyError(f"unknown configuration {configuration!r} for {name}; known: {known}")
    relative = instrument.configurations[configuration].line_table
    if relative is None:
        raise ValueError(f"configuration {configuration!r} of {name} has no line table")
    path = instrument_dir(name) / Path(relative)
    table = pd.read_csv(path, comment="#", encoding="utf-8")
    missing = [column for column in LINE_TABLE_COLUMNS if column not in table.columns]
    if missing:
        raise ValueError(f"{path.name}: missing line-table columns {missing}")
    table = table.loc[:, list(LINE_TABLE_COLUMNS)].copy()
    table["wavelength_nm"] = table["wavelength_nm"].astype(float)
    table["species"] = table["species"].astype(str)
    table["logI"] = table["logI"].astype(float)
    table["usable"] = table["usable"].astype(bool)
    return table.sort_values("wavelength_nm", kind="stable").reset_index(drop=True)


def _anchors_path(name: str, configuration: str, preset: str) -> Path:
    return instrument_dir(name) / "anchors" / f"{configuration}_{preset}.yaml"


def load_anchors(name: str, configuration: str, preset: str) -> list[tuple[float, float]]:
    """Stored ``(wavelength nm, pixel)`` anchor pairs for a configuration and layout preset.

    Read from ``anchors/<configuration>_<preset>.yaml``; an empty list when the file does
    not exist. Pixels are in the cropped frame of the layout, measured on the trace that
    :func:`load_anchor_trace` names.
    """
    path = _anchors_path(name, configuration, preset)
    if not path.is_file():
        return []
    data = load_yaml(path)
    pairs = data.get("pairs") or []
    anchors: list[tuple[float, float]] = []
    for pair in pairs:
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError(f"{path.name}: each anchor pair must be [wavelength_nm, pixel]")
        anchors.append((float(pair[0]), float(pair[1])))
    return anchors


def load_anchor_trace(name: str, configuration: str, preset: str) -> int:
    """The trace index the stored anchor pairs were measured on (the file's ``trace:``).

    ``0`` when the anchors file does not exist or does not name a trace.
    """
    path = _anchors_path(name, configuration, preset)
    if not path.is_file():
        return 0
    trace = load_yaml(path).get("trace")
    return 0 if trace is None else int(trace)
