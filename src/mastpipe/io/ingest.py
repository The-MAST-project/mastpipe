"""Ingest: scan a folder of raw FITS files into :class:`FrameRef` objects and load frames.

:func:`scan_folder` opens every file matching a pattern, resolves each image HDU's header
through the instrument's :class:`~mastpipe.io.headers.HeaderMap` and classifies it,
producing one :class:`~mastpipe.core.model.FrameRef` per image HDU (the primary HDU
when it carries data, else every image extension, each read with the primary header
underneath its own). A manifest CSV can override any metadata field frame by frame; it
is how laboratory data without instrument keywords enter the pipeline.
:func:`load_frame` reads one HDU as ``float32``, crops it and applies the detector's
orientation. Downstream code sees only the pipeline vocabulary carried by
:class:`~mastpipe.core.model.FrameMeta`.
"""

from __future__ import annotations

import dataclasses
import hashlib
from collections import ChainMap
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.config.models import DetectorConfig
from mastpipe.core.model import ARC_KINDS, EXPOSURE_TYPES, Frame, FrameMeta, FrameRef
from mastpipe.io.headers import HeaderMap

CHECKSUM_LENGTH = 16
CHUNK_SIZE = 1 << 20
DEFAULT_CAMERA = "main"
MANIFEST_RULE = "manifest"

#: Header-map fields copied verbatim into :class:`FrameMeta` when they resolve.
PASSTHROUGH_FIELDS = (
    "configuration",
    "exposure_time",
    "mid_time",
    "readout_mode",
    "detector_temperature",
    "lamp_state",
    "fiber_preset",
    "target",
)


def _parse_binning(text: str) -> tuple[int, int]:
    """``"2x2"`` (or ``"2,2"`` / ``"2 2"``) to ``(2, 2)``."""
    parts = text.lower().replace("x", " ").replace(",", " ").split()
    if len(parts) != 2:
        raise ValueError(f"binning must be '<x>x<y>', got {text!r}")
    return int(parts[0]), int(parts[1])


#: Metadata fields a manifest may override, with the parser applied to the CSV cell.
MANIFEST_FIELDS: dict[str, Callable[[str], Any]] = {
    "instrument": str,
    "camera": str,
    "exposure_type": str,
    "configuration": str,
    "exposure_time": float,
    "mid_time": str,
    "binning": _parse_binning,
    "readout_mode": str,
    "detector_temperature": float,
    "lamp_state": str,
    "fiber_preset": str,
    "target": str,
    "arc_kind": str,
    "classification_rule": str,
}
MANIFEST_KEYS = ("source", "extension")


def file_checksum(path: Path | str) -> str:
    """The first 16 hex characters of the file's SHA-256, read in chunks."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()[:CHECKSUM_LENGTH]


def scan_folder(
    folder: Path | str,
    instrument: str,
    header_map: HeaderMap | Mapping[str, Any],
    pattern: str = "*.fit*",
    manifest: pd.DataFrame | Path | str | None = None,
    unreadable: list[Path] | None = None,
) -> list[FrameRef]:
    """One :class:`FrameRef` per image HDU of every file matching ``pattern``, in name order.

    ``header_map`` is a parsed :class:`HeaderMap` or the raw mapping it is built from.
    A file astropy cannot open raises ``OSError``, unless ``unreadable`` is a list: the
    path is appended to it and scanning continues. ``manifest`` (a table from
    :func:`read_manifest` or its path) is applied to the result.
    """
    if not isinstance(header_map, HeaderMap):
        header_map = HeaderMap.from_dict(header_map)
    refs: list[FrameRef] = []
    for path in sorted(Path(folder).glob(pattern)):
        if not path.is_file():
            continue
        try:
            with fits.open(path) as hdul:
                hdus = _image_hdus(hdul)
        except OSError:
            if unreadable is None:
                raise
            unreadable.append(path)
            continue
        checksum = file_checksum(path)
        for index, header in hdus:
            meta = _build_meta(header_map, header, path.name, instrument, index, checksum)
            refs.append(FrameRef(path=path, extension=index, meta=meta))
    if manifest is not None:
        if not isinstance(manifest, pd.DataFrame):
            manifest = read_manifest(manifest)
        refs = apply_manifest(refs, manifest)
    return refs


def _image_hdus(hdul: fits.HDUList) -> list[tuple[int, Mapping[str, Any]]]:
    """``(index, header)`` of the HDUs that count as frames: the primary when it has
    data, else every image extension with the primary header merged underneath."""
    primary = hdul[0].header
    if primary.get("NAXIS", 0) > 0:
        return [(0, fits.Header(primary))]
    found: list[tuple[int, Mapping[str, Any]]] = []
    for index, hdu in enumerate(hdul):
        if index > 0 and getattr(hdu, "is_image", False) and hdu.header.get("NAXIS", 0) > 0:
            found.append((index, ChainMap(fits.Header(hdu.header), fits.Header(primary))))
    return found


def _build_meta(
    header_map: HeaderMap,
    header: Mapping[str, Any],
    filename: str,
    instrument: str,
    extension: int,
    checksum: str,
) -> FrameMeta:
    fields, unresolved = header_map.resolve(header, filename)
    exposure_type, rule = header_map.classify(fields)
    arc_kind = fields.get("arc_kind")
    if arc_kind is not None and arc_kind not in ARC_KINDS:
        arc_kind = None
        unresolved.append("arc_kind")
    camera = fields.get("camera")
    passthrough = {name: fields.get(name) for name in PASSTHROUGH_FIELDS}
    return FrameMeta(
        instrument=instrument,
        camera=DEFAULT_CAMERA if camera is None else str(camera),
        exposure_type=exposure_type,
        binning=(_axis_binning(fields.get("binning_x")), _axis_binning(fields.get("binning_y"))),
        arc_kind=arc_kind,
        checksum=checksum,
        source=filename,
        extension=extension,
        unresolved=tuple(unresolved),
        classification_rule=rule,
        **passthrough,
    )


def _axis_binning(value: Any) -> int:
    return 1 if value is None else int(value)


def load_frame(ref: FrameRef, detector: DetectorConfig) -> Frame:
    """The HDU of ``ref`` as a ``float32`` :class:`Frame`: cropped to ``detector.crop``
    (``x_min, x_max, y_min, y_max``, end-exclusive) when set, then transposed, flipped
    in x and flipped in y as the detector asks, in that order. No variance or mask yet."""
    with fits.open(ref.path) as hdul:
        data = np.asarray(hdul[ref.extension].data, dtype=np.float32)
    if data.ndim != 2:
        raise ValueError(
            f"{ref.meta.source}[{ref.extension}]: expected a 2-D image, got {data.ndim}-D"
        )
    if detector.crop is not None:
        x_min, x_max, y_min, y_max = detector.crop
        data = data[y_min:y_max, x_min:x_max]
    if detector.transpose:
        data = data.T
    if detector.flip_x:
        data = data[:, ::-1]
    if detector.flip_y:
        data = data[::-1, :]
    return Frame(data=np.ascontiguousarray(data), meta=ref.meta)


def read_manifest(path: Path | str) -> pd.DataFrame:
    """A manifest CSV as a table of strings: ``source`` (file name), optional ``extension``,
    and any column of :data:`MANIFEST_FIELDS`; empty cells mean "no override"."""
    path = Path(path)
    table = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
    if "source" not in table.columns:
        raise ValueError(f"{path.name}: a manifest needs a 'source' column")
    unknown = [
        column for column in table.columns if column not in (*MANIFEST_KEYS, *MANIFEST_FIELDS)
    ]
    if unknown:
        raise ValueError(f"{path.name}: unknown manifest columns {unknown}")
    return table


def apply_manifest(refs: Iterable[FrameRef], manifest: pd.DataFrame) -> list[FrameRef]:
    """New refs with the manifest's non-empty cells written over their metadata.

    A row matches a ref by ``source``; when the row has a non-empty ``extension`` it
    matches that HDU only. Rows apply in order (later rows win), rows matching no ref are
    ignored. Overriding ``exposure_type`` records ``"manifest"`` as the classification
    rule unless the row sets one. ``ValueError`` for a value outside
    :data:`~mastpipe.core.model.EXPOSURE_TYPES` or :data:`~mastpipe.core.model.ARC_KINDS`.
    """
    overrides = [_row_overrides(row, index) for index, row in manifest.iterrows()]
    out: list[FrameRef] = []
    for ref in refs:
        meta = ref.meta
        for source, extension, values in overrides:
            if source == meta.source and (extension is None or extension == ref.extension):
                meta = dataclasses.replace(meta, **values)
        out.append(dataclasses.replace(ref, meta=meta) if meta is not ref.meta else ref)
    return out


def _row_overrides(row: pd.Series, index: Any) -> tuple[str, int | None, dict[str, Any]]:
    label = f"manifest row {index}"
    source = str(row["source"]).strip()
    if not source:
        raise ValueError(f"{label}: empty source")
    extension_text = str(row["extension"]).strip() if "extension" in row else ""
    extension = int(extension_text) if extension_text else None
    values: dict[str, Any] = {}
    for name, parse in MANIFEST_FIELDS.items():
        if name not in row:
            continue
        text = str(row[name]).strip()
        if not text:
            continue
        try:
            values[name] = parse(text)
        except ValueError as error:
            raise ValueError(f"{label}: bad {name} value {text!r}") from error
    if "exposure_type" in values and values["exposure_type"] not in EXPOSURE_TYPES:
        raise ValueError(f"{label}: exposure_type must be one of {EXPOSURE_TYPES}")
    if "arc_kind" in values and values["arc_kind"] not in ARC_KINDS:
        raise ValueError(f"{label}: arc_kind must be one of {ARC_KINDS}")
    if "exposure_type" in values:
        values.setdefault("classification_rule", MANIFEST_RULE)
    return source, extension, values
