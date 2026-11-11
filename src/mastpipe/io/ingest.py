"""Ingest: scan a folder of raw FITS files into :class:`FrameRef` objects and load frames.

:func:`scan_folder` walks a night folder and its subfolders (the control software files
each acquisition in a folder of its own), opens every FITS file, resolves each image HDU's
header through the instrument's :class:`~mastpipe.io.headers.HeaderMap` and classifies it,
producing one :class:`~mastpipe.core.model.FrameRef` per image HDU (the primary HDU
when it carries data, else every image extension, each read with the primary header
underneath its own); a symbolic link to a folder is not descended into. A frame's
``source`` is its path relative to the night folder, with forward slashes; a frame whose
header names another instrument is left out. The ``mid_time`` of a frame is derived from
the header's time field and the exposure time according to the field's ``reference``
(``start``: plus half the exposure; ``end``: minus half; ``mid``: as given) and converted
to UTC, a stamp without an offset read according to the field's ``timezone``
(:mod:`mastpipe.core.times`); ``FrameMeta.timezone`` records how. When the exposure time
is unknown the time is not shifted and ``unresolved`` lists ``mid_time``. A frame's
``targets`` come from the map's per-unit ``targets`` field (one target per telescope
unit), else from its single ``target`` field as ``(("*", target),)``, which applies to
every science fiber (a missing card, a blank string or a logical card names none); when
neither resolves ``unresolved`` lists ``targets``, except for the
:data:`CALIBRATION_TYPES`, which have no target. A manifest CSV can
override any metadata field frame by frame (its ``mid_time`` is taken as the mid-time,
UTC unless it carries an offset, which is honoured); a row names its frames by relative
path, glob pattern or unique file name. It is how laboratory data without instrument
keywords enter the pipeline.
:func:`load_frame` reads one HDU as ``float32``, crops it and applies the detector's
orientation. Downstream code sees only the pipeline vocabulary carried by
:class:`~mastpipe.core.model.FrameMeta`.
"""

from __future__ import annotations

import dataclasses
import fnmatch
import hashlib
import re
from collections import ChainMap
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path, PureWindowsPath
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.config.models import DetectorConfig
from mastpipe.core.model import (
    ALL_TELESCOPES,
    ARC_KINDS,
    EXPOSURE_TYPES,
    Frame,
    FrameMeta,
    FrameRef,
)
from mastpipe.core.times import ASSUMED_UTC as ASSUMED_UTC  # lives in core.times; importable here too
from mastpipe.core.times import UTC, read_time, to_utc, utc_iso, zone_label
from mastpipe.io.headers import HeaderMap
from mastpipe.vocab import is_unit_name

CHECKSUM_LENGTH = 16
CHUNK_SIZE = 1 << 20
DEFAULT_CAMERA = "main"
MANIFEST_RULE = "manifest"
TIME_FIELD = "mid_time"
"""The header-map field carrying the exposure's time, whose ``reference`` says which instant."""
INSTRUMENT_FIELD = "instrument"
"""The header-map field naming the instrument that wrote a frame: frames of another are skipped."""
DEFAULT_PATTERN = "*.fit*"
"""``scan_folder``'s default ``pattern``: any file with a FITS suffix (:data:`FITS_SUFFIXES`)."""
FITS_SUFFIXES = (".fits", ".fit", ".fits.fz")
"""File-name endings of FITS files, compared case-insensitively."""
GLOB_CHARS = frozenset("*?[")
TARGETS_FIELD = "targets"
"""The header-map field (``per_unit``) naming each telescope unit's target."""
TARGET_FIELD = "target"
"""The header-map field naming one target for every fiber, the fallback of ``targets``."""
CALIBRATION_TYPES = ("bias", "dark", "flat", "arc", "background")
"""Exposure types with no target: ``targets`` is never listed as unresolved for them."""

#: Header-map fields copied verbatim into :class:`FrameMeta` when they resolve.
PASSTHROUGH_FIELDS = (
    "configuration",
    "exposure_time",
    "mid_time",
    "timezone",
    "readout_mode",
    "amplifier",
    "detector_temperature",
    "lamp_state",
    "fiber_preset",
)


def _manifest_time(text: str) -> str:
    """A manifest ``mid_time`` cell, read as UTC unless it carries an offset, in the
    pipeline's form; ``ValueError`` when it does not parse."""
    return utc_iso(to_utc(text, UTC))


def _parse_binning(text: str) -> tuple[int, int]:
    """``"2x2"`` (or ``"2,2"`` / ``"2 2"``) to ``(2, 2)``."""
    parts = text.lower().replace("x", " ").replace(",", " ").split()
    if len(parts) != 2:
        raise ValueError(f"binning must be '<x>x<y>', got {text!r}")
    return int(parts[0]), int(parts[1])


def _parse_targets(text: str) -> tuple[tuple[str, str], ...]:
    """``"mast01=HD 1234;mast02=HD 5678"`` to ``(("mast01", "HD 1234"), ("mast02", "HD 5678"))``;
    an entry without ``=`` (such as a bare ``"HD 1234"``) applies to every telescope, as
    ``("*", "HD 1234")``. Empty entries are ignored; ``ValueError`` for an empty name or
    target, a name that is not a unit name, or a telescope named twice."""
    pairs: list[tuple[str, str]] = []
    for entry in text.split(";"):
        if not entry.strip():
            continue
        if "=" in entry:
            telescope, _, target = (part.strip() for part in entry.partition("="))
            if not is_unit_name(telescope):
                raise ValueError(f"{telescope!r} is not a telescope unit name")
        else:
            telescope, target = ALL_TELESCOPES, entry.strip()
        if not target:
            raise ValueError(f"no target for {telescope!r}")
        if any(telescope == seen for seen, _ in pairs):
            raise ValueError(f"{telescope!r} is named twice")
        pairs.append((telescope, target))
    return tuple(pairs)


def _parse_single_target(text: str) -> tuple[tuple[str, str], ...]:
    """A legacy ``target`` cell: one target for every telescope."""
    return ((ALL_TELESCOPES, text),)


#: Metadata fields a manifest may override, with the parser applied to the CSV cell.
MANIFEST_FIELDS: dict[str, Callable[[str], Any]] = {
    "instrument": str,
    "camera": str,
    "exposure_type": str,
    "configuration": str,
    "exposure_time": float,
    "mid_time": _manifest_time,
    "binning": _parse_binning,
    "readout_mode": str,
    "amplifier": str,
    "detector_temperature": float,
    "lamp_state": str,
    "fiber_preset": str,
    "targets": _parse_targets,
    "target": _parse_single_target,
    "arc_kind": str,
    "classification_rule": str,
}
MANIFEST_KEYS = ("source", "extension")
LEGACY_MANIFEST_FIELDS = {TARGET_FIELD: TARGETS_FIELD}
"""Manifest columns kept for older manifests, with the field they set (a ``target`` cell
is one target for every telescope); the new column wins when a row fills both."""


class ManifestError(ValueError):
    """An error in a manifest's content, naming the file or the row: a file pandas cannot
    read as UTF-8 CSV, a missing ``source`` column or an unknown one, or a row whose
    source, extension or values are bad (:func:`read_manifest`, :func:`apply_manifest`).
    A ``ValueError``, so callers catching that still catch it; a command catches this
    class to report the manifest as a user error, apart from the run's own errors."""


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
    pattern: str = DEFAULT_PATTERN,
    manifest: pd.DataFrame | Path | str | None = None,
    unreadable: list[Path] | None = None,
    recursive: bool = True,
    skipped: list[tuple[str, str]] | None = None,
    unmatched: list[str] | None = None,
) -> list[FrameRef]:
    """One :class:`FrameRef` per image HDU of every FITS file under ``folder``, in path order.

    ``recursive`` walks the subfolders too; ``source`` is the file's path relative to
    ``folder`` with forward slashes. Files and folders whose name starts with ``.`` are
    ignored. With the default ``pattern`` a FITS file is one whose name ends in
    :data:`FITS_SUFFIXES` (any case); another ``pattern`` is matched against the file name
    instead (``identify`` scans a single file that way). ``header_map`` is a parsed
    :class:`HeaderMap` or the raw mapping it is built from. A frame whose ``instrument``
    field resolves to another instrument than ``instrument`` (compared case-insensitively)
    is left out, and appended to ``skipped`` as ``(source, reason)`` when it is a list (once
    per file, however many of its HDUs are frames).
    ``FileNotFoundError`` when ``folder`` is not a directory; an empty list when it holds
    no matching file. A file astropy cannot open raises ``OSError``, unless ``unreadable``
    is a list: the path is appended to it and scanning continues. ``manifest`` (a table
    from :func:`read_manifest` or its path) is applied to the result, its rows that match
    no frame appended to ``unmatched`` when it is a list (see :func:`apply_manifest`).
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise FileNotFoundError(f"{folder} is not a directory")
    if not isinstance(header_map, HeaderMap):
        header_map = HeaderMap.from_dict(header_map)
    refs: list[FrameRef] = []
    for path, source in _fits_files(folder, pattern, recursive):
        try:
            with fits.open(path) as hdul:
                hdus = _image_hdus(hdul)
        except OSError:
            if unreadable is None:
                raise
            unreadable.append(path)
            continue
        checksum: str | None = None
        other: str | None = None
        for index, header in hdus:
            fields, unresolved = header_map.resolve(header, path.name)
            written = fields.get(INSTRUMENT_FIELD)
            if written is not None and str(written).lower() != instrument.lower():
                reason = f"instrument {written!r} is not {instrument!r}"
                if skipped is not None and reason != other:
                    skipped.append((source, reason))
                other = reason
                continue
            if checksum is None:
                checksum = file_checksum(path)
            meta = _build_meta(header_map, fields, unresolved, source, instrument, index, checksum)
            refs.append(FrameRef(path=path, extension=index, meta=meta))
    if manifest is not None:
        if not isinstance(manifest, pd.DataFrame):
            manifest = read_manifest(manifest)
        refs = apply_manifest(refs, manifest, unmatched)
    return refs


def _fits_files(folder: Path, pattern: str, recursive: bool) -> list[tuple[Path, str]]:
    """``(path, source)`` of the files to scan under ``folder``, sorted by ``source``."""
    found: list[tuple[Path, str]] = []
    for path in folder.rglob("*") if recursive else folder.iterdir():
        relative = path.relative_to(folder)
        if any(part.startswith(".") for part in relative.parts) or not path.is_file():
            continue
        if pattern == DEFAULT_PATTERN:
            if not path.name.lower().endswith(FITS_SUFFIXES):
                continue
        elif not fnmatch.fnmatch(path.name, pattern):
            continue
        found.append((path, relative.as_posix()))
    return sorted(found, key=lambda item: item[1])


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
    fields: dict[str, Any],
    unresolved: list[str],
    source: str,
    instrument: str,
    extension: int,
    checksum: str,
) -> FrameMeta:
    """The :class:`FrameMeta` of one HDU from its resolved header-map ``fields``.

    ``instrument`` is the instrument being reduced, so a header without an ``instrument``
    field leaves nothing unresolved."""
    unresolved = [name for name in unresolved if name != INSTRUMENT_FIELD]
    exposure_type, rule = header_map.classify(fields)
    spec = header_map.fields.get(TIME_FIELD)
    reference = "mid" if spec is None else spec.reference
    timezone = UTC if spec is None else spec.timezone
    mid_time, derived, zone = _mid_time(fields.get(TIME_FIELD), fields.get("exposure_time"), reference, timezone)
    fields[TIME_FIELD] = mid_time
    fields["timezone"] = zone
    if not derived:
        unresolved.append(TIME_FIELD)
    arc_kind = fields.get("arc_kind")
    if arc_kind is not None and arc_kind not in ARC_KINDS:
        arc_kind = None
        unresolved.append("arc_kind")
    # a logical card (OBJECT = F) names no target; OBJECT = 0 is the target "0"
    targets = tuple(pair for pair in fields.get(TARGETS_FIELD) or () if _names_target(pair[1]))
    if not targets:
        single = fields.get(TARGET_FIELD)
        targets = ((ALL_TELESCOPES, single),) if _names_target(single) else ()
    targets = tuple((str(telescope), str(target)) for telescope, target in targets)
    unresolved = [name for name in unresolved if name not in (TARGETS_FIELD, TARGET_FIELD)]
    if not targets and exposure_type not in CALIBRATION_TYPES:
        unresolved.append(TARGETS_FIELD)
    camera = fields.get("camera")
    passthrough = {name: fields.get(name) for name in PASSTHROUGH_FIELDS}
    return FrameMeta(
        instrument=instrument,
        camera=DEFAULT_CAMERA if camera is None else str(camera),
        exposure_type=exposure_type,
        binning=(_axis_binning(fields.get("binning_x")), _axis_binning(fields.get("binning_y"))),
        targets=targets,
        arc_kind=arc_kind,
        checksum=checksum,
        source=source,
        extension=extension,
        unresolved=tuple(unresolved),
        classification_rule=rule,
        **passthrough,
    )


def _names_target(value: Any) -> bool:
    """Whether a header value names a target: not a missing card, a blank string or a
    logical card."""
    return value is not None and not isinstance(value, (bool, np.bool_)) and bool(str(value).strip())


def _axis_binning(value: Any) -> int:
    return 1 if value is None else int(value)


def _mid_time(value: Any, exposure_time: Any, reference: str, timezone: str = UTC) -> tuple[str | None, bool, str | None]:
    """``(mid_time, derived, zone)``: the exposure's UTC mid-time from the header's time
    ``value`` (an ISO-8601 string) and ``exposure_time`` in seconds, according to
    ``reference``; a ``value`` without an offset is read according to ``timezone``.

    ``zone`` is the :attr:`FrameMeta.timezone` label (:func:`~mastpipe.core.times.zone_label`):
    ``"utc"`` for a value with its own offset (or under ``timezone="utc"``),
    ``"assumed-utc"`` under ``"assume-utc"``, else the zone name. ``derived`` is false
    when the time is given but could not be turned into a mid-time: without an exposure
    time it is converted but not shifted; a time that does not parse is passed on as it
    came, with no zone. ``None`` stays ``None`` (and counts as derived: there is nothing
    to note twice).
    """
    if value is None:
        return None, True, None
    text = str(value)
    try:
        stamp, explicit = read_time(text, timezone)
    except ValueError:
        return text, False, None
    zone = zone_label(timezone, explicit)
    if reference == "mid":
        return utc_iso(stamp), True, zone
    if exposure_time is None:
        return utc_iso(stamp), False, zone
    try:
        half = pd.Timedelta(seconds=float(exposure_time) / 2.0)
    except (TypeError, ValueError):
        return utc_iso(stamp), False, zone
    mid = stamp + half if reference == "start" else stamp - half
    return utc_iso(mid), True, zone


def _binned_crop(detector: DetectorConfig, binning: tuple[int, int]) -> tuple[int, int, int, int] | None:
    """``detector.crop`` (unbinned pixels) in the pixels of a frame binned ``binning``
    (``(x, y)`` factors): each bound integer-divided by its axis's factor."""
    if detector.crop is None:
        return None
    (bx, by), (x_min, x_max, y_min, y_max) = binning, detector.crop
    return x_min // bx, x_max // bx, y_min // by, y_max // by


def expected_shape(detector: DetectorConfig, binning: tuple[int, int] = (1, 1)) -> tuple[int, int]:
    """Shape of a frame of ``detector`` binned ``binning`` (``(x, y)`` factors) after
    :func:`load_frame`'s crop and orientation. ``detector.shape`` and ``crop`` are in
    unbinned pixels; a binned frame is cropped and sized in binned pixels, each axis divided
    by its factor before the transpose."""
    bx, by = binning
    crop = _binned_crop(detector, binning)
    if crop is None:
        ny, nx = detector.shape[0] // by, detector.shape[1] // bx
    else:
        x_min, x_max, y_min, y_max = crop
        ny, nx = y_max - y_min, x_max - x_min
    return (nx, ny) if detector.transpose else (ny, nx)


def load_frame(ref: FrameRef, detector: DetectorConfig) -> Frame:
    """The HDU of ``ref`` as a ``float32`` :class:`Frame`: cropped to ``detector.crop``
    (``x_min, x_max, y_min, y_max``, end-exclusive, in unbinned pixels; a frame binned
    ``ref.meta.binning`` is cropped at the bounds divided by its factors) when set, then
    transposed, flipped in x and flipped in y as the detector asks, in that order. A frame
    of the detector's size loads to :func:`expected_shape` of its binning. No variance or
    mask yet."""
    with fits.open(ref.path) as hdul:
        data = np.asarray(hdul[ref.extension].data, dtype=np.float32)
    if data.ndim != 2:
        raise ValueError(f"{ref.meta.source}[{ref.extension}]: expected a 2-D image, got {data.ndim}-D")
    crop = _binned_crop(detector, ref.meta.binning)
    if crop is not None:
        x_min, x_max, y_min, y_max = crop
        data = data[y_min:y_max, x_min:x_max]
    if detector.transpose:
        data = data.T
    if detector.flip_x:
        data = data[:, ::-1]
    if detector.flip_y:
        data = data[::-1, :]
    return Frame(data=np.ascontiguousarray(data), meta=ref.meta)


def read_manifest(path: Path | str) -> pd.DataFrame:
    """A manifest CSV as a table of strings: ``source`` (a relative path, a glob pattern or
    a file name, see :func:`apply_manifest`), optional ``extension``, and any column of
    :data:`MANIFEST_FIELDS`; empty cells mean "no override". :class:`ManifestError` naming
    the file when it is not UTF-8 CSV (pandas' parser errors and ``UnicodeDecodeError`` are
    ``ValueError`` subclasses), cannot be opened (a folder, say), has rows longer than its
    header, lacks ``source`` or has an unknown column; ``FileNotFoundError`` when it does
    not exist."""
    path = Path(path)
    try:
        table = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8")
    except FileNotFoundError:
        raise
    except OSError as error:
        raise ManifestError(f"{path.name}: cannot be read: {error}") from error
    except ValueError as error:
        raise ManifestError(f"{path.name}: not a readable UTF-8 CSV file: {error}") from error
    if not isinstance(table.index, pd.RangeIndex):
        # every row one cell longer than the header: pandas made the first column the
        # index, shifting the others under the wrong names
        raise ManifestError(f"{path.name}: rows have more cells than the header")
    if "source" not in table.columns:
        raise ManifestError(f"{path.name}: a manifest needs a 'source' column")
    unknown = [column for column in table.columns if column not in (*MANIFEST_KEYS, *MANIFEST_FIELDS)]
    if unknown:
        raise ManifestError(f"{path.name}: unknown manifest columns {unknown}")
    return table


def apply_manifest(
    refs: Iterable[FrameRef],
    manifest: pd.DataFrame,
    unmatched: list[str] | None = None,
) -> list[FrameRef]:
    """New refs with the manifest's non-empty cells written over their metadata.

    A ``\\`` in a row's ``source`` is a folder separator, read as ``/`` (a manifest written
    on Windows). A row's ``source`` matches a ref whose relative path it is; else, when it
    holds glob characters, the refs whose relative path it matches (``*`` and ``?`` stay
    within a folder, ``**`` crosses folders); else a bare file name matches the one file of
    that name, and ``ValueError`` names the candidates when several files have it. When the
    row has a non-empty ``extension`` it matches that HDU only. Rows apply in order (later
    rows win). A field a row sets is removed from the frame's ``unresolved``, and a
    ``target`` cell resolves ``targets`` (as does a row that makes the frame one of the
    :data:`CALIBRATION_TYPES`, which have no target); a row that makes it another type
    while it has no target lists ``targets``, as ingest does for such a header. A row
    matching no ref is appended to ``unmatched`` as ``"manifest row {index}: {source}"``
    when it is a list (:func:`names_non_fits_file` tells whether its source names a file
    the scan never read). Reads no file. Overriding ``exposure_type`` records
    ``"manifest"`` as the classification rule unless the row sets one.
    :class:`ManifestError` (a ``ValueError``) naming the row for an empty ``source``, a bad
    glob, a file name several files have, an ``extension`` that is not an integer, a value
    outside :data:`~mastpipe.core.model.EXPOSURE_TYPES` or
    :data:`~mastpipe.core.model.ARC_KINDS`, or a cell its parser rejects (a ``mid_time``
    or ``targets`` that does not parse, say).
    """
    refs = list(refs)
    sources = sorted({ref.meta.source for ref in refs})
    rows: list[tuple[set[str], int | None, dict[str, Any]]] = []
    for index, row in manifest.iterrows():
        source, extension, values = _row_overrides(row, index)
        matched = _match_sources(source, sources, f"manifest row {index}")
        rows.append((matched, extension, values))
        if unmatched is not None and not any(_applies(ref, matched, extension) for ref in refs):
            # the runner splits this label on its first ": " to get the source back
            unmatched.append(f"manifest row {index}: {source}")
    out: list[FrameRef] = []
    for ref in refs:
        meta = ref.meta
        for matched, extension, values in rows:
            if _applies(ref, matched, extension):
                cleared = set(values)
                exposure_type = values.get("exposure_type")
                if exposure_type in CALIBRATION_TYPES:
                    cleared.add(TARGETS_FIELD)
                unresolved = tuple(name for name in meta.unresolved if name not in cleared)
                meta = dataclasses.replace(meta, **values, unresolved=unresolved)
                # a type with a target but none named: listed, as ingest lists it for a header
                lacks_target = not meta.targets and TARGETS_FIELD not in meta.unresolved
                if exposure_type not in (None, *CALIBRATION_TYPES) and lacks_target:
                    meta = dataclasses.replace(meta, unresolved=(*meta.unresolved, TARGETS_FIELD))
        out.append(dataclasses.replace(ref, meta=meta) if meta is not ref.meta else ref)
    return out


def _applies(ref: FrameRef, matched: set[str], extension: int | None) -> bool:
    return ref.meta.source in matched and (extension is None or extension == ref.extension)


def names_non_fits_file(folder: Path | str, source: str) -> bool:
    """Whether a manifest row's ``source`` names an existing file under the night
    ``folder`` whose name is not a FITS file name (:data:`FITS_SUFFIXES`), so the scan never
    read it: a relative path, or a bare file name anywhere below ``folder``; a ``\\`` is a
    folder separator, as in :func:`apply_manifest`. A glob, a rooted source (``/x``) and one
    naming a drive (``C:x``) name no such file, and hidden files and folders are not looked
    at, as the scan does not. A runner uses it to say why a row in :func:`apply_manifest`'s
    ``unmatched`` matched no frame.
    """
    folder = Path(folder)
    source = source.strip().replace("\\", "/")
    # a drive is tested the Windows way on every platform (rglob refuses one on Windows)
    if not source or GLOB_CHARS & set(source) or source.startswith("/") or PureWindowsPath(source).drive:
        return False
    candidates = [folder.joinpath(*source.split("/"))] if "/" in source else list(folder.rglob(source))
    for path in candidates:
        # relative_to is lexical: it keeps a ".." part (rejected by the dot test below) and
        # raises only for a path on another drive (Windows "C:dir/x.txt", refused above)
        try:
            relative = path.relative_to(folder)
        except ValueError:
            continue
        if any(part.startswith(".") for part in relative.parts) or not path.is_file():
            continue
        if not path.name.lower().endswith(FITS_SUFFIXES):
            return True
    return False


def _match_sources(source: str, sources: list[str], label: str) -> set[str]:
    """The relative ``sources`` a manifest row's ``source`` names (see :func:`apply_manifest`)."""
    if source in sources:
        return {source}
    if GLOB_CHARS & set(source):
        try:
            regex = _glob_regex(source)
        except re.error as error:  # a reversed range, say: reported like the other row errors
            raise ManifestError(f"{label}: bad glob {source!r}: {error}") from error
        return {candidate for candidate in sources if regex.fullmatch(candidate)}
    if "/" in source:
        return set()
    named = [candidate for candidate in sources if candidate.rsplit("/", 1)[-1] == source]
    if len(named) > 1:
        raise ManifestError(f"{label}: {source!r} names several files: {', '.join(named)}")
    return set(named)


def _glob_regex(pattern: str) -> re.Pattern[str]:
    """A glob over relative POSIX paths as a regex: ``*`` and ``?`` do not cross ``/``,
    ``**`` does (``**/`` also matches no folder at all); ``[...]`` is a character class
    (``[!...]`` its complement), never matching ``/``. As in :mod:`fnmatch`, a ``]`` right
    after ``[`` or ``[!`` is a member, ``a-z`` is a range, and every other character of a
    class is literal (``[``, ``&``, ``|``, ``~``, ``^``, a ``-`` that starts or ends it); a
    ``[`` without a closing ``]`` is literal. ``re.error`` for a reversed range
    (``[9-0]``). ``pattern`` holds no ``\\`` (a manifest source has its backslashes read
    as ``/`` before it gets here)."""
    parts: list[str] = []
    i = 0
    while i < len(pattern):
        close = _class_end(pattern, i) if pattern[i] == "[" else -1
        if pattern.startswith("**/", i):
            parts.append("(?:.*/)?")
            i += 3
        elif pattern.startswith("**", i):
            parts.append(".*")
            i += 2
        elif pattern[i] == "*":
            parts.append("[^/]*")
            i += 1
        elif pattern[i] == "?":
            parts.append("[^/]")
            i += 1
        elif close != -1:
            body = pattern[i + 1 : close]
            negate = body.startswith("!")
            members = _class_members(body[1:] if negate else body)
            parts.append(f"[^/{members}]" if negate else f"(?!/)[{members}]")
            i = close + 1
        else:
            parts.append(re.escape(pattern[i]))
            i += 1
    return re.compile("".join(parts))


def _class_end(pattern: str, start: int) -> int:
    """Index of the ``]`` closing the class opened at ``pattern[start]``, ``-1`` without
    one; a ``]`` first in the class (after ``[`` or ``[!``) is a member, as in fnmatch."""
    j = start + 1
    if pattern.startswith("!", j):
        j += 1
    if pattern.startswith("]", j):
        j += 1
    return pattern.find("]", j)


def _class_members(body: str) -> str:
    """The body of a glob class as the body of a regex class: ranges ``x-y`` kept, every
    character escaped, so nothing in it is special to ``re`` (no nested set, no ``&&``)."""
    members: list[str] = []
    k = 0
    while k < len(body):
        if k + 2 < len(body) and body[k + 1] == "-":
            members.append(f"{re.escape(body[k])}-{re.escape(body[k + 2])}")
            k += 3
        else:
            members.append(re.escape(body[k]))
            k += 1
    return "".join(members)


def _row_overrides(row: pd.Series, index: Any) -> tuple[str, int | None, dict[str, Any]]:
    label = f"manifest row {index}"
    source = str(row["source"]).strip().replace("\\", "/")  # a Windows separator is a folder separator
    if not source:
        raise ManifestError(f"{label}: empty source")
    extension_text = str(row["extension"]).strip() if "extension" in row else ""
    try:
        extension = int(extension_text) if extension_text else None
    except ValueError as error:
        raise ManifestError(f"{label}: bad extension {extension_text!r}") from error
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
            raise ManifestError(f"{label}: bad {name} value {text!r}: {error}") from error
    for legacy, field in LEGACY_MANIFEST_FIELDS.items():
        if legacy in values:
            values.setdefault(field, values.pop(legacy))
    if "exposure_type" in values and values["exposure_type"] not in EXPOSURE_TYPES:
        raise ManifestError(f"{label}: exposure_type must be one of {EXPOSURE_TYPES}")
    if "arc_kind" in values and values["arc_kind"] not in ARC_KINDS:
        raise ManifestError(f"{label}: arc_kind must be one of {ARC_KINDS}")
    if "exposure_type" in values:
        values.setdefault("classification_rule", MANIFEST_RULE)
    if "mid_time" in values:
        values["timezone"] = UTC  # a manifest's mid_time is taken as given, in UTC
    return source, extension, values
