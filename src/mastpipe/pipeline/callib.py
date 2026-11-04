"""The calibration library: master calibrations and reference solutions under one folder,
found by the keys that decide where they apply.

The library is a folder (``out/calib`` by default) holding the FITS products written by
:mod:`mastpipe.io.products` and an ``index.json`` listing them, one entry per product with
the validity keys that :meth:`CalibLibrary.find` and :meth:`CalibLibrary.find_reference`
match on. The index is the source of truth: a product dropped into the folder by hand is
not seen until it is added. Entries carry the ``camera`` their product belongs to, because
reading a reference back needs the camera's pixel size to rebuild its grating geometry
(:func:`mastpipe.steps.geometry.grating_geometry` from the instrument definition).

References are time-keyed: :meth:`CalibLibrary.find_reference` returns the reference with
the latest ``valid_from`` at or before the requested time (the latest overall without a
time), or the one named by its identifier. Times are compared as parsed timestamps (a
naive time is UTC), never as strings. Adding a product whose file name or identifier is
already indexed replaces that entry and overwrites the file. This is pipeline-layer code:
it reads and writes files.

Several processes may share one library (a quick-look and an ``identify`` run against the
same folder). Products and the index are written atomically, and every change to the
index re-reads it under an exclusive lock on ``index.lock``
(:func:`mastpipe.io.atomic.locked`), so concurrent writers keep each other's entries.
Lookups read the index as last loaded or written by this instance, without the lock: an
entry another process added since is not seen until this instance next writes or is
re-created, which is acceptable for finding calibrations.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from mastpipe.core.model import MasterCalib, ReferenceSolution
from mastpipe.core.times import to_utc
from mastpipe.instruments.registry import load_instrument
from mastpipe.io.atomic import atomic_write, locked
from mastpipe.io.products import (
    product_name,
    read_master,
    read_reference,
    reference_name,
    write_master,
    write_reference,
)
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.geometry import GratingGeometry, grating_geometry

INDEX_FILE = "index.json"
LOCK_FILE = "index.lock"
"""The advisory lock file next to the index, taken for every read-modify-write of it."""
REFERENCE_KIND = "reference"
"""The ``kind`` of a reference entry; masters carry their ``MasterCalib.kind``."""


def registry_geometry(camera: str) -> Any:
    """A ``geometry_lookup`` for :func:`~mastpipe.io.products.read_reference` through the
    instrument registry, with the pixel size of ``camera``."""

    def lookup(instrument: str, configuration: str) -> GratingGeometry:
        definition = load_instrument(instrument)
        if camera not in definition.cameras:
            raise KeyError(
                f"unknown camera {camera!r} for {instrument}: {sorted(definition.cameras)}"
            )
        block = definition.configurations[configuration]
        return grating_geometry(block.geometry, definition.cameras[camera].pixel_size_um)

    return lookup


class CalibLibrary:
    """Masters and references under ``root``, indexed in ``root/index.json``."""

    def __init__(self, root: Path | str) -> None:
        self.root = Path(root)
        # Read under the lock: on Windows ``os.replace`` onto an index another process
        # has open fails, so a reader must not hold it while a writer saves. Readers wait
        # at most one write. The lock helper creates ``root`` and the lock file, so an
        # absent or empty library still opens (with no entries).
        with locked(self.root / LOCK_FILE):
            self._entries: list[dict[str, Any]] = self._load()
        self._references: dict[str, ReferenceSolution] = {}

    # ------------------------------------------------------------------------ the index

    def _load(self) -> list[dict[str, Any]]:
        path = self.root / INDEX_FILE
        if not path.is_file():
            return []
        entries = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(entries, list):
            raise ValueError(f"{path}: the index must be a JSON list of entries")
        return entries

    def _save(self) -> None:
        """Write the index atomically; callers hold the lock (see :meth:`_replace`)."""
        text = json.dumps(self._entries, indent=1, sort_keys=True, ensure_ascii=True)
        with atomic_write(self.root / INDEX_FILE) as tmp:
            with tmp.open("w", encoding="utf-8", newline="") as handle:
                handle.write(text + "\n")

    def _replace(self, key: str, value: Any, entry: dict[str, Any]) -> None:
        """Append ``entry``, dropping any entry whose ``key`` already equals ``value``.

        Under the lock on ``index.lock`` the index is re-read from disk first, so entries
        another process added since this instance loaded it are kept; ``TimeoutError``
        from :func:`~mastpipe.io.atomic.locked` when the lock is not free in time.
        """
        with locked(self.root / LOCK_FILE):
            self._entries = [e for e in self._load() if e.get(key) != value]
            self._entries.append(entry)
            self._save()

    def list(self, kind: str | None = None) -> list[dict[str, Any]]:
        """Copies of the index entries, those of ``kind`` only when given."""
        return [dict(e) for e in self._entries if kind is None or e.get("kind") == kind]

    # --------------------------------------------------------------------------- masters

    def add(self, master: MasterCalib, prov: Provenance, meta_for_name: Any) -> Path:
        """Write ``master`` as ``root/<product_name(meta_for_name, kind)>`` and index it.

        The entry records the master's validity (``camera``, ``binning``, ``readout_mode``,
        ``detector_temperature``, ``mid_time``, ``exposure_time``), the instrument and
        configuration of ``meta_for_name`` and the number of frames. A master of the same
        file name replaces the earlier entry.
        """
        name = product_name(meta_for_name, master.kind)
        path = self.root / name
        write_master(path, master, prov)
        validity = master.validity
        binning = validity.get("binning") or meta_for_name.binning
        entry = {
            "kind": master.kind,
            "path": name,
            "instrument": meta_for_name.instrument,
            "configuration": meta_for_name.configuration,
            "camera": validity.get("camera") or meta_for_name.camera,
            "binning": [int(b) for b in binning],
            "readout_mode": validity.get("readout_mode"),
            "detector_temperature": validity.get("detector_temperature"),
            "mid_time": validity.get("mid_time"),
            "exposure_time": validity.get("exposure_time"),
            "n_frames": int(master.n_frames),
        }
        self._replace("path", name, entry)
        return path

    def find(
        self,
        kind: str,
        camera: str,
        binning: tuple[int, int],
        readout_mode: str | None,
        configuration: str | None = None,
        temperature: float | None = None,
        temperature_tolerance: float = 2.0,
        near_time: str | None = None,
    ) -> tuple[MasterCalib, Path] | None:
        """The master of ``kind`` that applies, read back, with its path; ``None`` without one.

        ``camera``, ``binning`` and ``readout_mode`` must match exactly; ``configuration``
        must match when given (a flat belongs to its grating); the entry's temperature must
        be within ``temperature_tolerance`` of ``temperature`` when both are known. Among
        several matches the nearest ``mid_time`` to ``near_time`` wins (an entry without a
        time counts as farthest), the latest ``mid_time`` without ``near_time``. Reads this
        instance's view of the index without the lock; a master another process added
        since may be missed (a stale read is acceptable here).
        """
        matches = []
        for entry in self._entries:
            if entry.get("kind") != kind or entry.get("camera") != camera:
                continue
            if tuple(entry.get("binning") or ()) != tuple(int(b) for b in binning):
                continue
            if entry.get("readout_mode") != readout_mode:
                continue
            if configuration is not None and entry.get("configuration") != configuration:
                continue
            entry_temperature = entry.get("detector_temperature")
            if (
                temperature is not None
                and entry_temperature is not None
                and abs(float(entry_temperature) - float(temperature)) > temperature_tolerance
            ):
                continue
            matches.append(entry)
        if not matches:
            return None
        if near_time is not None:
            target = to_utc(near_time)

            def distance(entry: dict[str, Any]) -> float:
                time = entry.get("mid_time")
                if time is None:
                    return float("inf")
                return abs((to_utc(time) - target).total_seconds())

            chosen = min(matches, key=distance)
        else:
            chosen = max(
                matches,
                key=lambda e: (
                    to_utc(e["mid_time"])
                    if e.get("mid_time")
                    else pd.Timestamp.min.tz_localize("UTC")
                ),
            )
        path = self.root / str(chosen["path"])
        return read_master(path)[0], path

    # ------------------------------------------------------------------------ references

    def add_reference(self, ref: ReferenceSolution, prov: Provenance, camera: str) -> Path:
        """Write ``ref`` as ``root/<identifier>_wref.fits`` and index it with ``camera``.

        A reference of the same identifier (the same instrument, configuration, preset and
        arc time) replaces the earlier entry and file.
        """
        name = reference_name(ref)
        path = self.root / name
        write_reference(path, ref, prov)
        entry = {
            "kind": REFERENCE_KIND,
            "path": name,
            "identifier": ref.identifier,
            "instrument": ref.instrument,
            "configuration": ref.configuration,
            "fiber_preset": ref.fiber_preset,
            "camera": camera,
            "created": ref.created,
            "valid_from": ref.valid_from,
            "hash": ref.hash,
            "source_checksums": list(ref.source_checksums),
        }
        self._references.pop(name, None)
        self._replace("identifier", ref.identifier, entry)
        return path

    def _reference_entries(
        self, instrument: str, configuration: str, fiber_preset: str
    ) -> list[dict[str, Any]]:
        return [
            e
            for e in self._entries
            if e.get("kind") == REFERENCE_KIND
            and e.get("instrument") == instrument
            and e.get("configuration") == configuration
            and e.get("fiber_preset") == fiber_preset
        ]

    def _read_reference(self, entry: dict[str, Any]) -> ReferenceSolution:
        name = str(entry["path"])
        cached = self._references.get(name)
        if cached is None:
            cached, _ = read_reference(self.root / name, registry_geometry(str(entry["camera"])))
            self._references[name] = cached
        return cached

    def find_reference(
        self,
        instrument: str,
        configuration: str,
        fiber_preset: str,
        at_time: str | None = None,
        identifier: str | None = None,
    ) -> ReferenceSolution | None:
        """The reference of ``instrument`` / ``configuration`` / ``fiber_preset`` to use.

        With ``identifier`` that reference (``None`` when the library does not hold it);
        otherwise the one with the latest ``valid_from`` at or before ``at_time``, or the
        latest overall when ``at_time`` is ``None``. ``None`` when no reference qualifies.
        The file is read once per library instance (``ValueError`` from
        :func:`~mastpipe.io.products.read_reference` when it does not verify).
        """
        entries = self._reference_entries(instrument, configuration, fiber_preset)
        if identifier is not None:
            entries = [e for e in entries if e.get("identifier") == identifier]
        elif at_time is not None:
            limit = to_utc(at_time)
            entries = [e for e in entries if to_utc(str(e["valid_from"])) <= limit]
        if not entries:
            return None
        chosen = max(entries, key=lambda e: (to_utc(str(e["valid_from"])), str(e["identifier"])))
        return self._read_reference(chosen)

    def reference_for_source(
        self, instrument: str, configuration: str, fiber_preset: str, checksum: str
    ) -> dict[str, Any] | None:
        """The index entry of the reference solved from the frame with ``checksum``, if any.

        Reads this instance's view of the index without the lock; a reference another
        process added since may be missed (a stale read is acceptable here).
        """
        for entry in self._reference_entries(instrument, configuration, fiber_preset):
            if checksum in (entry.get("source_checksums") or []):
                return dict(entry)
        return None
