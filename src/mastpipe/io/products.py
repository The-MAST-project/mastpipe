"""FITS products: every dataclass the steps produce, written with its provenance and read
back as the same dataclass.

Layout rules shared by every product:

- The primary header carries the product-level cards: the
  :class:`~mastpipe.pipeline.provenance.Provenance` cards, the frame's
  :class:`~mastpipe.core.model.FrameMeta` as numbered ``MPNMK`` / ``MPMKnn`` (field name) /
  ``MPMVnn`` (JSON value) cards, and the product's own scalars. Extension headers carry only
  table-level cards (which columns are populated, the drift scalars).
- Keywords are at most eight characters and all start with ``MP``; mappings are numbered
  key/value schemes with a count card (``MPNVAL`` / ``MPVALKnn`` / ``MPVALVnn`` for a
  master's validity, ``MPNBC`` / ``MPBCKnn`` / ``MPBCVnn`` for barycentric terms, ``MPNXP``
  / ``MPXPKnn`` / ``MPXPVnn`` for the free-form provenance dicts of
  :class:`~mastpipe.core.model.MasterCalib` and
  :class:`~mastpipe.core.model.ExposureSpectra`). Values are JSON; a JSON list of a
  validity value or of the ``binning`` / ``unresolved`` fields reads back as a tuple.
- No raw instrument keyword ever reaches a product header: frame metadata travels only as
  the ``MPMKnn`` / ``MPMVnn`` pairs. An unset scalar is the string ``"none"``.
- Per-fiber tables hold one row per fiber with vector columns on the fiber's pixel grid.
  Optional arrays are NaN-filled (masks ``True``-filled) and flagged by a per-row ``HAS_*``
  boolean column plus a table-level presence card: ``MPHASSKY`` / ``MPHASFLX`` /
  ``MPHASTEL`` for the kept ``skysub`` / ``fluxcal`` / ``telluric`` versions and ``MPGRDWAV``
  / ``MPGRDSKY`` / ``MPGRDCON`` for the wavelength, sky and continuum grids.
- Tables from ``DataFrame``s (``DRIFT``, ``CONTRIB``, ``EXCLUDED``) keep the frame's
  columns, upper-cased on disk; ``None`` strings become empty strings.
- Files are written with ``overwrite=True``; arrays read back are copies in native byte
  order.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Callable, Iterable
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

from mastpipe.core.model import (
    COMBINED_VERSIONS,
    LINE_COLUMNS,
    VERSIONS,
    CombinedSpectrum,
    DriftModel,
    ExposureSpectra,
    FrameMeta,
    MasterCalib,
    QCReport,
    ReferenceSolution,
    Spectrum,
    SpectrumVersion,
    Trace,
    TraceSet,
    WavelengthSolution,
)
from mastpipe.pipeline.provenance import NONE, Provenance
from mastpipe.steps.geometry import GratingGeometry
from mastpipe.steps.wavecal import solutions_hash
from mastpipe.steps.wavemodel import model_from_solution

MAX_TRACE_DEGREE = 3
"""Highest trace polynomial degree the ``TRACES`` table stores (``COEFFS`` has 4 slots)."""
MAX_PARAMS = 12
"""Slots of the ``PARAMS`` column of the ``SOLUTIONS`` table."""

VERSION_CARDS = {"skysub": "MPHASSKY", "fluxcal": "MPHASFLX", "telluric": "MPHASTEL"}
"""Presence card of each kept version on the ``SPECTRA`` table."""
GRID_CARDS = {"wavelength": "MPGRDWAV", "sky": "MPGRDSKY", "continuum": "MPGRDCON"}
"""Presence card of each optional per-fiber grid on the ``SPECTRA`` table."""

META_SCHEME = ("MPNMK", "MPMK", "MPMV")
VALIDITY_SCHEME = ("MPNVAL", "MPVALK", "MPVALV")
BARYCENTRIC_SCHEME = ("MPNBC", "MPBCK", "MPBCV")
EXTRA_SCHEME = ("MPNXP", "MPXPK", "MPXPV")
"""``(count card, key prefix, value prefix)`` of the numbered mapping schemes (2 digits)."""
SOURCE_SCHEME = ("MPNSRC", "MPSRC", 3)
"""Count card, prefix and digits of a reference's source checksums."""

META_TUPLE_FIELDS = ("binning", "unresolved")
SAFE_TARGET = re.compile(r"[^A-Za-z0-9_-]")


# ------------------------------------------------------------------------------ names


def product_name(meta: FrameMeta, kind: str, ext: str = ".fits") -> str:
    """``{instrument}_{configuration|noconf}_{camera}_{YYYYmmddTHHMMSS|notime}_{kind}{ext}``.

    The stamp is the frame's UTC mid-time; a time-zone-aware ``mid_time`` is converted.
    """
    stamp = "notime"
    if meta.mid_time is not None:
        time = pd.Timestamp(meta.mid_time)
        if time.tzinfo is not None:
            time = time.tz_convert("UTC")
        stamp = time.strftime("%Y%m%dT%H%M%S")
    configuration = meta.configuration or "noconf"
    return f"{meta.instrument}_{configuration}_{meta.camera}_{stamp}_{kind}{ext}"


def reference_name(ref: ReferenceSolution) -> str:
    """``{identifier}_wref.fits``."""
    return f"{ref.identifier}_wref.fits"


def combined_name(meta: FrameMeta, target: str) -> str:
    """:func:`product_name` of kind ``comb_<target>``, the target reduced to ``[A-Za-z0-9_-]``."""
    return product_name(meta, f"comb_{SAFE_TARGET.sub('_', target)}")


# ---------------------------------------------------------------------------- helpers


def _json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=True, separators=(",", ":"))


CARD_LENGTH = 80
LONG_STRING = 68
"""Longest string value astropy writes on one card; longer ones go to CONTINUE cards."""


def _set(header: fits.Header, key: str, value: Any, comment: str) -> None:
    """``header[key] = (value, comment)`` with the comment trimmed to the room the card has.

    astropy warns and truncates when keyword, value and comment overrun 80 columns, and
    only moves to CONTINUE cards when the string value alone exceeds :data:`LONG_STRING`
    characters; trimming the comment here keeps every write warning-free.
    """
    if isinstance(value, str) and len(value) <= LONG_STRING:
        rendered = max(len(value) + value.count("'") + 2, 10)
    elif isinstance(value, str):
        rendered = 0  # CONTINUE cards carry the comment along
    else:
        rendered = 20
    if rendered:
        room = CARD_LENGTH - 10 - rendered - 3
        comment = comment[: max(room, 0)]
    header[key] = (value, comment)


def _add_cards(header: fits.Header, cards: Iterable[tuple[str, Any, str]]) -> None:
    for key, value, comment in cards:
        _set(header, key, value, comment)


def _add_mapping(header: fits.Header, scheme: tuple[str, str, str], mapping: dict) -> None:
    """Write ``mapping`` as the numbered key / JSON-value cards of ``scheme``."""
    count_key, key_prefix, value_prefix = scheme
    if len(mapping) >= 100:
        raise ValueError(f"{count_key}: {len(mapping)} entries exceed the {key_prefix} cards")
    _set(header, count_key, len(mapping), f"number of {key_prefix} entries")
    for n, (key, value) in enumerate(mapping.items(), 1):
        _set(header, f"{key_prefix}{n:02d}", str(key), f"{key_prefix} key {n}")
        _set(header, f"{value_prefix}{n:02d}", _json(value), f"{key_prefix} value {n} (JSON)")


def _read_mapping(header: fits.Header, scheme: tuple[str, str, str]) -> dict[str, Any] | None:
    """The inverse of :func:`_add_mapping`; ``None`` when the scheme's count card is absent."""
    count_key, key_prefix, value_prefix = scheme
    if count_key not in header:
        return None
    count = int(header[count_key])
    return {
        str(header[f"{key_prefix}{n:02d}"]): json.loads(header[f"{value_prefix}{n:02d}"])
        for n in range(1, count + 1)
    }


def _tuplify(value: Any) -> Any:
    return tuple(value) if isinstance(value, list) else value


def _add_meta(header: fits.Header, meta: FrameMeta) -> None:
    _add_mapping(header, META_SCHEME, asdict(meta))


def _read_meta(header: fits.Header) -> FrameMeta:
    fields = _read_mapping(header, META_SCHEME) or {}
    for name in META_TUPLE_FIELDS:
        if name in fields:
            fields[name] = tuple(fields[name])
    return FrameMeta(**fields)


def _add_provenance(header: fits.Header, prov: Provenance) -> None:
    _add_cards(header, prov.to_header_cards())


def _optional_card(value: float | bool | None) -> float | bool | str:
    return NONE if value is None else value


def _optional_float(value: Any) -> float | None:
    return None if isinstance(value, str) and value == NONE else float(value)


def _native(array: Any) -> np.ndarray:
    """A copy of a FITS array in native byte order, detached from the file."""
    array = np.asarray(array)
    return array.astype(array.dtype.type)


def _strings(column: Any) -> list[str]:
    return [str(value) for value in column]


def _none_if_empty(text: str) -> str | None:
    return text if text else None


def _string_column(name: str, values: Iterable[str | None]) -> fits.Column:
    texts = ["" if value is None else str(value) for value in values]
    width = max([1, *(len(text) for text in texts)])
    return fits.Column(name=name, format=f"A{width}", array=np.array(texts, dtype=f"U{width}"))


def _column(name: str, fmt: str, array: Any) -> fits.Column:
    return fits.Column(name=name, format=fmt, array=np.asarray(array))


def _vector_format(n: int, code: str) -> str:
    return f"{n}{code}"


def _table(name: str, columns: list[fits.Column]) -> fits.BinTableHDU:
    return fits.BinTableHDU.from_columns(columns, name=name)


def _table_from_frame(name: str, frame: pd.DataFrame) -> fits.BinTableHDU:
    """A ``DataFrame`` as a table, columns upper-cased: bool -> ``L``, int -> ``K``, float ->
    ``D``, anything else as strings (``None`` -> ``""``)."""
    columns = []
    for label in frame.columns:
        series = frame[label]
        kind = series.dtype.kind
        if kind == "b":
            columns.append(_column(str(label).upper(), "L", series.to_numpy(dtype=bool)))
        elif kind in "iu":
            columns.append(_column(str(label).upper(), "K", series.to_numpy(dtype=np.int64)))
        elif kind == "f":
            columns.append(_column(str(label).upper(), "D", series.to_numpy(dtype=float)))
        else:
            columns.append(_string_column(str(label).upper(), series.tolist()))
    return _table(name, columns)


def _frame_from_table(hdu: fits.BinTableHDU) -> pd.DataFrame:
    """The inverse of :func:`_table_from_frame`, columns lower-cased."""
    data = hdu.data
    out: dict[str, Any] = {}
    for label in hdu.columns.names:
        column = data[label]
        if column.dtype.kind in "SU":
            out[label.lower()] = _strings(column)
        else:
            out[label.lower()] = _native(column)
    return pd.DataFrame(out)


def _write(path: Path | str, hdus: list[Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return path


# ---------------------------------------------------------------------------- masters


def write_master(path: Path | str, master: MasterCalib, prov: Provenance) -> Path:
    """The master image as the primary HDU with its variance as the ``VARIANCE`` image
    (absent when ``master.variance is None``).

    Cards: ``MPKIND``, ``MPNFRM``, the validity as ``MPVALKnn`` / ``MPVALVnn``, the master's
    own provenance dict as ``MPXPKnn`` / ``MPXPVnn`` and the provenance cards.
    """
    primary = fits.PrimaryHDU(np.asarray(master.data))
    header = primary.header
    _set(header, "MPKIND", master.kind, "master calibration kind")
    _set(header, "MPNFRM", int(master.n_frames), "number of frames stacked")
    _add_mapping(header, VALIDITY_SCHEME, master.validity)
    _add_mapping(header, EXTRA_SCHEME, master.provenance)
    _add_provenance(header, prov)
    hdus: list[Any] = [primary]
    if master.variance is not None:
        hdus.append(fits.ImageHDU(np.asarray(master.variance), name="VARIANCE"))
    return _write(path, hdus)


def read_master(path: Path | str) -> tuple[MasterCalib, Provenance]:
    """The inverse of :func:`write_master`; validity lists come back as tuples."""
    with fits.open(path) as hdul:
        header = hdul[0].header
        validity = {
            k: _tuplify(v) for k, v in (_read_mapping(header, VALIDITY_SCHEME) or {}).items()
        }
        master = MasterCalib(
            kind=str(header["MPKIND"]),
            data=_native(hdul[0].data),
            validity=validity,
            n_frames=int(header["MPNFRM"]),
            provenance=_read_mapping(header, EXTRA_SCHEME) or {},
            variance=_native(hdul["VARIANCE"].data) if "VARIANCE" in hdul else None,
        )
        return master, Provenance.from_header(header)


# ----------------------------------------------------------------------------- traces


def write_traceset(path: Path | str, traceset: TraceSet, prov: Provenance) -> Path:
    """Tables ``TRACES`` (one row per trace: ``INDEX, DEGREE, COEFFS(4), YCENTER, SLOPE,
    RMS, NPOINTS, SNR, ROLE, TELESCOPE, LAYOUT, APERTURE``) and ``SAMPLES`` (one row per
    centroid sample: ``INDEX, X, Y, KEPT``).

    ``COEFFS`` holds the ``polyval`` coefficients right-aligned to degree
    :data:`MAX_TRACE_DEGREE`, ``DEGREE`` saying how many are real; a higher degree is a
    ``ValueError``. ``LAYOUT`` is -1 and ``APERTURE`` NaN when unset; empty ``ROLE`` /
    ``TELESCOPE`` mean ``None``. Primary cards: ``MPSHIFT`` (``"none"`` when unmatched),
    ``MPMISS`` and ``MPEXTRA`` (comma-joined layout indices, empty when none), ``MPNY``,
    ``MPNX`` and the provenance cards.
    """
    traces = traceset.traces
    n = len(traces)
    coeffs = np.zeros((n, MAX_TRACE_DEGREE + 1))
    degree = np.zeros(n, dtype=np.int32)
    for i, trace in enumerate(traces):
        c = np.asarray(trace.coeffs, dtype=float)
        if c.ndim != 1 or c.size == 0 or c.size > MAX_TRACE_DEGREE + 1:
            raise ValueError(
                f"trace {trace.index}: polynomial degree {c.size - 1} is not storable "
                f"(the TRACES table holds degree 0 to {MAX_TRACE_DEGREE})"
            )
        degree[i] = c.size - 1
        coeffs[i, -c.size :] = c
    traces_hdu = _table(
        "TRACES",
        [
            _column("INDEX", "J", [t.index for t in traces]),
            _column("DEGREE", "J", degree),
            _column("COEFFS", _vector_format(MAX_TRACE_DEGREE + 1, "D"), coeffs),
            _column("YCENTER", "D", [t.y_center for t in traces]),
            _column("SLOPE", "D", [t.slope for t in traces]),
            _column("RMS", "D", [t.rms for t in traces]),
            _column("NPOINTS", "J", [t.n_points for t in traces]),
            _column("SNR", "D", [t.snr for t in traces]),
            _string_column("ROLE", [t.role for t in traces]),
            _string_column("TELESCOPE", [t.telescope for t in traces]),
            _column(
                "LAYOUT", "J", [-1 if t.layout_index is None else t.layout_index for t in traces]
            ),
            _column(
                "APERTURE", "D", [np.nan if t.aperture is None else t.aperture for t in traces]
            ),
        ],
    )
    index = np.concatenate([np.full(len(t.x_samples), t.index) for t in traces] or [[]])
    samples_hdu = _table(
        "SAMPLES",
        [
            _column("INDEX", "J", index.astype(np.int32)),
            _column("X", "D", np.concatenate([t.x_samples for t in traces] or [[]])),
            _column("Y", "D", np.concatenate([t.y_samples for t in traces] or [[]])),
            _column("KEPT", "L", np.concatenate([t.kept for t in traces] or [[]]).astype(bool)),
        ],
    )
    primary = fits.PrimaryHDU()
    header = primary.header
    _set(header, "MPSHIFT", _optional_card(traceset.layout_shift), "layout shift, px")
    _set(header, "MPMISS", ",".join(str(k) for k in traceset.missing), "layout traces not found")
    _set(header, "MPEXTRA", ",".join(str(k) for k in traceset.extra), "traces not in the layout")
    _set(header, "MPNY", int(traceset.shape[0]), "frame rows")
    _set(header, "MPNX", int(traceset.shape[1]), "frame columns")
    _add_provenance(header, prov)
    return _write(path, [primary, traces_hdu, samples_hdu])


def _parse_indices(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split(",") if part.strip())


def read_traceset(path: Path | str) -> tuple[TraceSet, Provenance]:
    """The inverse of :func:`write_traceset`."""
    with fits.open(path) as hdul:
        header = hdul[0].header
        rows = hdul["TRACES"].data
        samples = hdul["SAMPLES"].data
        sample_index = _native(samples["INDEX"])
        x_samples = _native(samples["X"])
        y_samples = _native(samples["Y"])
        kept = _native(samples["KEPT"]).astype(bool)
        coeffs_all = _native(rows["COEFFS"])
        roles = _strings(rows["ROLE"])
        telescopes = _strings(rows["TELESCOPE"])
        traces = []
        for i in range(len(rows)):
            index = int(rows["INDEX"][i])
            degree = int(rows["DEGREE"][i])
            selected = sample_index == index
            layout = int(rows["LAYOUT"][i])
            aperture = float(rows["APERTURE"][i])
            traces.append(
                Trace(
                    index=index,
                    coeffs=coeffs_all[i, MAX_TRACE_DEGREE - degree :].copy(),
                    y_center=float(rows["YCENTER"][i]),
                    slope=float(rows["SLOPE"][i]),
                    rms=float(rows["RMS"][i]),
                    n_points=int(rows["NPOINTS"][i]),
                    snr=float(rows["SNR"][i]),
                    x_samples=x_samples[selected],
                    y_samples=y_samples[selected],
                    kept=kept[selected],
                    role=_none_if_empty(roles[i]),
                    telescope=_none_if_empty(telescopes[i]),
                    layout_index=None if layout < 0 else layout,
                    aperture=None if math.isnan(aperture) else aperture,
                )
            )
        traceset = TraceSet(
            traces=traces,
            shape=(int(header["MPNY"]), int(header["MPNX"])),
            layout_shift=_optional_float(header["MPSHIFT"]),
            missing=_parse_indices(str(header["MPMISS"])),
            extra=_parse_indices(str(header["MPEXTRA"])),
        )
        return traceset, Provenance.from_header(header)


# ---------------------------------------------------------------------------- spectra


DUPLICATED_FIELDS = (
    "wavelength_status",
    "flux_status",
    "telluric_status",
    "reference_id",
    "reference_hash",
    "drift_arcs",
)
"""Fields :class:`~mastpipe.core.model.ExposureSpectra` and :class:`Provenance` both carry."""


def _unset(value: Any) -> bool:
    return value is None or value == NONE or (isinstance(value, list | tuple) and not value)


def _exposure_provenance(prov: Provenance, es: ExposureSpectra) -> Provenance:
    """``prov`` with its unset duplicated fields filled from the exposure.

    A field of :data:`DUPLICATED_FIELDS` that ``prov`` leaves unset (``None``, ``"none"`` or
    empty) is taken from ``es``; one that is set must agree with ``es``, else
    ``ValueError`` naming it: a provenance built with the wrong status is a bug upstream
    that must not be papered over on disk.
    """
    filled: dict[str, Any] = {}
    for name in DUPLICATED_FIELDS:
        theirs = getattr(es, name)
        ours = getattr(prov, name)
        if name == "drift_arcs":
            theirs, ours = list(theirs), list(ours)
        if _unset(ours):
            filled[name] = theirs
        elif ours != theirs:
            raise ValueError(f"provenance {name}={ours!r} disagrees with the exposure's {theirs!r}")
    return replace(prov, **filled)


DRIFT_SHIFT_SCHEME = ("MPDNSH", "MPDSH", 2)
DRIFT_STRETCH_SCHEME = ("MPDNST", "MPDST", 2)
"""Count card, prefix and digits of the drift shift / stretch polynomial coefficients."""


def _add_coefficients(header: fits.Header, scheme: tuple[str, str, int], coeffs: Any) -> None:
    count_key, prefix, digits = scheme
    values = np.asarray(coeffs, dtype=float).reshape(-1)
    if values.size >= 10**digits:
        raise ValueError(f"{values.size} coefficients exceed the {prefix} cards")
    _set(header, count_key, int(values.size), f"number of {prefix} coefficients")
    for n, value in enumerate(values, 1):
        _set(header, f"{prefix}{n:0{digits}d}", float(value), f"{prefix} coefficient {n}")


def _read_coefficients(header: fits.Header, scheme: tuple[str, str, int]) -> np.ndarray:
    count_key, prefix, digits = scheme
    count = int(header[count_key])
    return np.array([float(header[f"{prefix}{n:0{digits}d}"]) for n in range(1, count + 1)])


def _drift_hdu(drift: DriftModel, weights: dict[int, float] | None) -> fits.BinTableHDU:
    """``DriftModel.per_fiber`` plus a ``WEIGHT`` column (1.0 for a fiber ``weights`` lacks),
    with the model's scalars and polynomials as cards."""
    table = drift.per_fiber.copy()
    fibers = table["fiber"].to_numpy()
    table["weight"] = [1.0 if weights is None else float(weights.get(int(f), 1.0)) for f in fibers]
    hdu = _table_from_frame("DRIFT", table)
    header = hdu.header
    _set(header, "MPDTIME", drift.arc_mid_time, "drift arc mid-time, UTC")
    _set(header, "MPDSCAT", float(drift.scatter_px), "drift scatter, px")
    _set(header, "MPDUSE", bool(drift.usable), "drift model usable")
    _set(header, "MPDNPIX", int(drift.n_pix), "pixels along dispersion")
    _add_coefficients(header, DRIFT_SHIFT_SCHEME, drift.shift_coeffs)
    _add_coefficients(header, DRIFT_STRETCH_SCHEME, drift.stretch_coeffs)
    return hdu


def _read_drift(hdu: fits.BinTableHDU, arc_checksums: tuple[str, ...]) -> DriftModel:
    """The inverse of :func:`_drift_hdu`; the checksums come from the provenance cards."""
    header = hdu.header
    per_fiber = _frame_from_table(hdu).drop(columns=["weight"])
    return DriftModel(
        arc_checksums=arc_checksums,
        arc_mid_time=str(header["MPDTIME"]),
        per_fiber=per_fiber,
        shift_coeffs=_read_coefficients(header, DRIFT_SHIFT_SCHEME),
        stretch_coeffs=_read_coefficients(header, DRIFT_STRETCH_SCHEME),
        scatter_px=float(header["MPDSCAT"]),
        usable=bool(header["MPDUSE"]),
        n_pix=int(header["MPDNPIX"]),
    )


def _stack(spectra: list[Spectrum], get: Callable[[Spectrum], Any], fill: Any, dtype: Any):
    """Rows of ``get(spectrum)`` (``fill`` where ``None``) and the per-row presence flags."""
    n_pix = len(spectra[0].pixel)
    out = np.full((len(spectra), n_pix), fill, dtype=dtype)
    present = np.zeros(len(spectra), dtype=bool)
    for i, spectrum in enumerate(spectra):
        array = get(spectrum)
        if array is not None:
            out[i] = array
            present[i] = True
    return out, present


def write_spectra(
    path: Path | str,
    es: ExposureSpectra,
    prov: Provenance,
    drift: DriftModel | None = None,
    weights: dict[int, float] | None = None,
) -> Path:
    """Table ``SPECTRA``, one row per fiber: ``FIBER, ROLE, TELESCOPE, TARGET`` and the
    vector columns ``PIXEL, FLUX, ERROR, MASK``, the optional grids ``WAVELENGTH, SKY,
    CONTINUUM`` and ``FLUX_<V>, ERROR_<V>, MASK_<V>`` for every kept version ``V`` in
    :data:`~mastpipe.core.model.VERSIONS`, each optional set NaN- / ``True``-filled with a
    ``HAS_<name>`` column per row and the presence cards :data:`VERSION_CARDS` /
    :data:`GRID_CARDS` on the table.

    The primary header carries the frame metadata (``MPMKnn`` / ``MPMVnn``), the exposure's
    free-form provenance dict (``MPXPKnn`` / ``MPXPVnn``, absent when ``None``) and the
    provenance cards. The statuses, reference and drift arcs ``prov`` leaves unset are
    filled from ``es``; set ones must agree with ``es`` (``ValueError`` naming the field).

    ``drift``, when given, is written as the ``DRIFT`` table: ``DriftModel.per_fiber`` plus
    a ``WEIGHT`` column (``weights`` by fiber, 1.0 when absent) with the cards ``MPDTIME``
    (``arc_mid_time``), ``MPDSCAT``, ``MPDUSE``, ``MPDNPIX`` and the shift / stretch
    coefficients as ``MPDNSH`` / ``MPDSHnn`` and ``MPDNST`` / ``MPDSTnn``; the arc checksums
    are the provenance's ``MPDARCnn`` cards.

    ``ValueError`` without spectra or when the spectra differ in pixel-grid length.
    """
    spectra = es.spectra
    if not spectra:
        raise ValueError("write_spectra needs at least one spectrum")
    lengths = {len(spectrum.pixel) for spectrum in spectra}
    if len(lengths) != 1:
        raise ValueError(f"spectra of one exposure must share a pixel grid, got lengths {lengths}")
    n_pix = lengths.pop()
    if n_pix == 0:
        raise ValueError("spectra have an empty pixel grid")
    integral = all(np.issubdtype(np.asarray(s.pixel).dtype, np.integer) for s in spectra)
    columns = [
        _column("FIBER", "J", [s.fiber for s in spectra]),
        _string_column("ROLE", [s.role for s in spectra]),
        _string_column("TELESCOPE", [s.telescope for s in spectra]),
        _string_column("TARGET", [s.target for s in spectra]),
        _column(
            "PIXEL",
            _vector_format(n_pix, "K" if integral else "D"),
            _stack(spectra, lambda s: s.pixel, 0, np.int64 if integral else float)[0],
        ),
        _column(
            "FLUX", _vector_format(n_pix, "D"), _stack(spectra, lambda s: s.flux, np.nan, float)[0]
        ),
        _column(
            "ERROR",
            _vector_format(n_pix, "D"),
            _stack(spectra, lambda s: s.error, np.nan, float)[0],
        ),
        _column(
            "MASK", _vector_format(n_pix, "L"), _stack(spectra, lambda s: s.mask, True, bool)[0]
        ),
    ]
    presence: list[tuple[str, bool, str]] = []
    for name, card in GRID_CARDS.items():
        values, present = _stack(spectra, lambda s, n=name: getattr(s, n), np.nan, float)
        columns.append(_column(f"HAS_{name.upper()}", "L", present))
        columns.append(_column(name.upper(), _vector_format(n_pix, "D"), values))
        presence.append((card, bool(present.any()), f"{name} column populated"))
    for version in VERSIONS:
        tag = version.upper()
        flux, present = _stack(
            spectra, lambda s, v=version: _version_part(s, v, "flux"), np.nan, float
        )
        error, _ = _stack(spectra, lambda s, v=version: _version_part(s, v, "error"), np.nan, float)
        mask, _ = _stack(spectra, lambda s, v=version: _version_part(s, v, "mask"), True, bool)
        columns.append(_column(f"HAS_{tag}", "L", present))
        columns.append(_column(f"FLUX_{tag}", _vector_format(n_pix, "D"), flux))
        columns.append(_column(f"ERROR_{tag}", _vector_format(n_pix, "D"), error))
        columns.append(_column(f"MASK_{tag}", _vector_format(n_pix, "L"), mask))
        presence.append((VERSION_CARDS[version], bool(present.any()), f"{version} version kept"))
    spectra_hdu = _table("SPECTRA", columns)
    _add_cards(spectra_hdu.header, presence)

    primary = fits.PrimaryHDU()
    _add_meta(primary.header, es.meta)
    if es.provenance is not None:
        _add_mapping(primary.header, EXTRA_SCHEME, es.provenance)
    _add_provenance(primary.header, _exposure_provenance(prov, es))
    hdus: list[Any] = [primary, spectra_hdu]
    if drift is not None:
        hdus.append(_drift_hdu(drift, weights))
    return _write(path, hdus)


def _version_part(spectrum: Spectrum, version: str, part: str) -> np.ndarray | None:
    kept = spectrum.versions.get(version)
    return None if kept is None else getattr(kept, part)


def read_spectra(
    path: Path | str,
) -> tuple[ExposureSpectra, Provenance, DriftModel | None]:
    """The inverse of :func:`write_spectra`: ``(exposure, provenance, drift)``.

    ``versions`` and the optional grids are populated only where their presence card and
    the row's ``HAS_*`` flag are set; ``drift`` is ``None`` without a ``DRIFT`` table and
    otherwise the :class:`~mastpipe.core.model.DriftModel` rebuilt from it (its
    ``arc_checksums`` are the provenance's drift arcs; the ``WEIGHT`` column stays in the
    file). The exposure's statuses, reference and drift arcs are the provenance's.
    """
    with fits.open(path) as hdul:
        header = hdul[0].header
        table = hdul["SPECTRA"]
        data = table.data
        n = len(data)
        prov = Provenance.from_header(header)
        fibers = _native(data["FIBER"])
        roles = _strings(data["ROLE"])
        telescopes = _strings(data["TELESCOPE"])
        targets = _strings(data["TARGET"])
        pixel = _native(data["PIXEL"])
        flux = _native(data["FLUX"])
        error = _native(data["ERROR"])
        mask = _native(data["MASK"]).astype(bool)
        grids: dict[str, tuple[np.ndarray, np.ndarray] | None] = {}
        for name, card in GRID_CARDS.items():
            if bool(table.header[card]):
                grids[name] = (_native(data[name.upper()]), _native(data[f"HAS_{name.upper()}"]))
            else:
                grids[name] = None
        kept: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        for version in VERSIONS:
            if bool(table.header[VERSION_CARDS[version]]):
                tag = version.upper()
                kept[version] = (
                    _native(data[f"FLUX_{tag}"]),
                    _native(data[f"ERROR_{tag}"]),
                    _native(data[f"MASK_{tag}"]).astype(bool),
                    _native(data[f"HAS_{tag}"]).astype(bool),
                )
        spectra = []
        for i in range(n):
            versions = {
                version: SpectrumVersion(f[i], e[i], m[i])
                for version, (f, e, m, present) in kept.items()
                if present[i]
            }
            optional = {
                name: (values[0][i] if values is not None and values[1][i] else None)
                for name, values in grids.items()
            }
            spectra.append(
                Spectrum(
                    fiber=int(fibers[i]),
                    pixel=pixel[i],
                    flux=flux[i],
                    error=error[i],
                    mask=mask[i],
                    role=_none_if_empty(roles[i]),
                    telescope=_none_if_empty(telescopes[i]),
                    target=_none_if_empty(targets[i]),
                    wavelength=optional["wavelength"],
                    sky=optional["sky"],
                    continuum=optional["continuum"],
                    versions=versions,
                )
            )
        es = ExposureSpectra(
            meta=_read_meta(header),
            spectra=spectra,
            wavelength_status=prov.wavelength_status,
            flux_status=prov.flux_status,
            telluric_status=prov.telluric_status,
            reference_id=prov.reference_id,
            reference_hash=prov.reference_hash,
            drift_arcs=tuple(prov.drift_arcs),
            provenance=_read_mapping(header, EXTRA_SCHEME),
        )
        drift = _read_drift(hdul["DRIFT"], tuple(prov.drift_arcs)) if "DRIFT" in hdul else None
        return es, prov, drift


# --------------------------------------------------------------------------- reference


def write_reference(path: Path | str, ref: ReferenceSolution, prov: Provenance) -> Path:
    """Tables ``SOLUTIONS`` (one row per fiber: ``FIBER, KIND, NPARAMS, PARAMS(12), NPIX,
    RMSPX, NLINES, ANCHPX0, ANCHPX1, ANCHNM0, ANCHNM1, YCENTER, GEOMETRY``) and ``LINES``
    (one row per identified line: ``FIBER, X, WAVELENGTH, SPECIES, DXPX``).

    ``PARAMS`` is zero-padded to :data:`MAX_PARAMS` slots with ``NPARAMS`` real ones (more
    is a ``ValueError``); ``GEOMETRY`` is the solution's provenance dict as JSON. Primary
    cards: ``MPREFID``, ``MPREFH``, ``MPCREAT``, ``MPVALID``, ``MPINST``, ``MPCONF``,
    ``MPPRESET``, the source checksums as ``MPNSRC`` / ``MPSRCnnn`` and the provenance
    cards.
    """
    fibers = sorted(ref.solutions)
    solutions = [ref.solutions[fiber] for fiber in fibers]
    params = np.zeros((len(solutions), MAX_PARAMS))
    n_params = np.zeros(len(solutions), dtype=np.int32)
    for i, solution in enumerate(solutions):
        p = np.asarray(solution.params, dtype=float)
        if p.ndim != 1 or p.size > MAX_PARAMS:
            raise ValueError(
                f"fiber {solution.fiber}: {p.size} parameters exceed the {MAX_PARAMS} PARAMS slots"
            )
        n_params[i] = p.size
        params[i, : p.size] = p
    solutions_hdu = _table(
        "SOLUTIONS",
        [
            _column("FIBER", "J", fibers),
            _string_column("KIND", [s.kind for s in solutions]),
            _column("NPARAMS", "J", n_params),
            _column("PARAMS", _vector_format(MAX_PARAMS, "D"), params),
            _column("NPIX", "J", [s.n_pix for s in solutions]),
            _column("RMSPX", "D", [s.rms_px for s in solutions]),
            _column("NLINES", "J", [s.n_lines for s in solutions]),
            _column("ANCHPX0", "D", [s.anchor_px[0] for s in solutions]),
            _column("ANCHPX1", "D", [s.anchor_px[1] for s in solutions]),
            _column("ANCHNM0", "D", [s.anchor_nm[0] for s in solutions]),
            _column("ANCHNM1", "D", [s.anchor_nm[1] for s in solutions]),
            _column("YCENTER", "D", [s.trace_y_center for s in solutions]),
            _string_column("GEOMETRY", [_json(s.geometry) for s in solutions]),
        ],
    )
    line_fiber = np.concatenate([np.full(len(s.lines), s.fiber) for s in solutions] or [[]])
    lines_hdu = _table(
        "LINES",
        [
            _column("FIBER", "J", line_fiber.astype(np.int32)),
            _column("X", "D", _concat_lines(solutions, "x")),
            _column("WAVELENGTH", "D", _concat_lines(solutions, "wavelength")),
            _string_column("SPECIES", [sp for s in solutions for sp in s.lines["species"]]),
            _column("DXPX", "D", _concat_lines(solutions, "dx_px")),
        ],
    )
    primary = fits.PrimaryHDU()
    header = primary.header
    _set(header, "MPREFID", ref.identifier, "reference solution identifier")
    _set(header, "MPREFH", ref.hash, "hash of the solutions")
    _set(header, "MPCREAT", ref.created, "created, UTC")
    _set(header, "MPVALID", ref.valid_from, "valid from, UTC")
    _set(header, "MPINST", ref.instrument, "instrument")
    _set(header, "MPCONF", ref.configuration, "configuration")
    _set(header, "MPPRESET", ref.fiber_preset, "fiber layout preset")
    count_key, prefix, digits = SOURCE_SCHEME
    if len(ref.source_checksums) >= 10**digits:
        raise ValueError(f"{len(ref.source_checksums)} source checksums exceed the {prefix} cards")
    _set(header, count_key, len(ref.source_checksums), "number of source frames")
    for n, checksum in enumerate(ref.source_checksums, 1):
        _set(header, f"{prefix}{n:0{digits}d}", checksum, f"source frame checksum {n}")
    _add_provenance(header, prov)
    return _write(path, [primary, solutions_hdu, lines_hdu])


def _concat_lines(solutions: list[WavelengthSolution], column: str) -> np.ndarray:
    parts = [s.lines[column].to_numpy(dtype=float) for s in solutions]
    return np.concatenate(parts) if parts else np.zeros(0)


def _check_geometry(solution: WavelengthSolution, geom: GratingGeometry, identifier: str) -> None:
    """``ValueError`` when a stored grating constant differs from the configuration's."""
    for key, stored in solution.geometry.items():
        if not hasattr(geom, key):
            continue
        current = float(getattr(geom, key))
        if not math.isclose(float(stored), current, rel_tol=1e-9, abs_tol=1e-12):
            raise ValueError(
                f"reference {identifier}, fiber {solution.fiber}: stored geometry {key}={stored!r} "
                f"differs from the configuration's {current!r}"
            )


def read_reference(
    path: Path | str, geometry_lookup: Callable[[str, str], GratingGeometry]
) -> tuple[ReferenceSolution, Provenance]:
    """The inverse of :func:`write_reference`.

    ``geometry_lookup(instrument, configuration)`` supplies the grating geometry of the
    stored ``MPINST`` / ``MPCONF``; every solution's stored geometry is checked against
    it, its model is rebuilt with :func:`~mastpipe.steps.wavemodel.model_from_solution`
    and the solutions are re-hashed against ``MPREFH`` (when set), so an edited or
    mismatched reference is a ``ValueError`` here rather than a wrong wavelength later.
    """
    with fits.open(path) as hdul:
        header = hdul[0].header
        identifier = str(header["MPREFID"])
        instrument = str(header["MPINST"])
        configuration = str(header["MPCONF"])
        geom = geometry_lookup(instrument, configuration)
        rows = hdul["SOLUTIONS"].data
        lines = hdul["LINES"].data
        line_fiber = _native(lines["FIBER"])
        line_x = _native(lines["X"])
        line_wl = _native(lines["WAVELENGTH"])
        line_dx = _native(lines["DXPX"])
        line_species = np.array(_strings(lines["SPECIES"]), dtype=object)
        params_all = _native(rows["PARAMS"])
        kinds = _strings(rows["KIND"])
        geometries = _strings(rows["GEOMETRY"])
        solutions: dict[int, WavelengthSolution] = {}
        for i in range(len(rows)):
            fiber = int(rows["FIBER"][i])
            n_params = int(rows["NPARAMS"][i])
            selected = line_fiber == fiber
            table = pd.DataFrame(
                {
                    "x": line_x[selected],
                    "wavelength": line_wl[selected],
                    "species": list(line_species[selected]),
                    "dx_px": line_dx[selected],
                }
            )[list(LINE_COLUMNS)]
            solution = WavelengthSolution(
                fiber=fiber,
                kind=kinds[i],
                params=params_all[i, :n_params].copy(),
                n_pix=int(rows["NPIX"][i]),
                lines=table,
                rms_px=float(rows["RMSPX"][i]),
                n_lines=int(rows["NLINES"][i]),
                anchor_px=(float(rows["ANCHPX0"][i]), float(rows["ANCHPX1"][i])),
                anchor_nm=(float(rows["ANCHNM0"][i]), float(rows["ANCHNM1"][i])),
                trace_y_center=float(rows["YCENTER"][i]),
                geometry=json.loads(geometries[i]),
            )
            _check_geometry(solution, geom, identifier)
            model_from_solution(solution, geom)
            solutions[fiber] = solution
        stored_hash = str(header["MPREFH"])
        if stored_hash and solutions_hash(solutions) != stored_hash:
            raise ValueError(f"reference {identifier}: solutions do not match the stored hash")
        count_key, prefix, digits = SOURCE_SCHEME
        checksums = tuple(
            str(header[f"{prefix}{n:0{digits}d}"]) for n in range(1, int(header[count_key]) + 1)
        )
        ref = ReferenceSolution(
            instrument=instrument,
            configuration=configuration,
            fiber_preset=str(header["MPPRESET"]),
            solutions=solutions,
            created=str(header["MPCREAT"]),
            valid_from=str(header["MPVALID"]),
            source_checksums=checksums,
            identifier=identifier,
            hash=stored_hash,
        )
        return ref, Provenance.from_header(header)


# ---------------------------------------------------------------------------- combined


def write_combined(
    path: Path | str, cs: CombinedSpectrum, meta: FrameMeta, prov: Provenance
) -> Path:
    """Table ``COMBINED``, one row per wavelength: ``WAVELENGTH`` and ``FLUX_<V>, ERROR_<V>,
    MASK_<V>`` for every version present in ``cs.versions`` (``RAW`` included), plus the
    tables ``CONTRIB`` and ``EXCLUDED`` from the contributor and excluded frames.

    Primary cards: ``MPTARGET``, ``MPGROUP``, the barycentric terms as ``MPNBC`` /
    ``MPBCKnn`` / ``MPBCVnn``, the frame metadata and the provenance cards.
    """
    columns = [_column("WAVELENGTH", "D", np.asarray(cs.wavelength, dtype=float))]
    for version in COMBINED_VERSIONS:
        kept = cs.versions.get(version)
        if kept is None:
            continue
        tag = version.upper()
        columns.append(_column(f"FLUX_{tag}", "D", np.asarray(kept.flux, dtype=float)))
        columns.append(_column(f"ERROR_{tag}", "D", np.asarray(kept.error, dtype=float)))
        columns.append(_column(f"MASK_{tag}", "L", np.asarray(kept.mask, dtype=bool)))
    primary = fits.PrimaryHDU()
    header = primary.header
    _set(header, "MPTARGET", cs.target, "target")
    _set(header, "MPGROUP", cs.grouping, "grouping of the combined spectra")
    _add_mapping(header, BARYCENTRIC_SCHEME, cs.barycentric)
    _add_meta(header, meta)
    _add_provenance(header, prov)
    return _write(
        path,
        [
            primary,
            _table("COMBINED", columns),
            _table_from_frame("CONTRIB", cs.contributors),
            _table_from_frame("EXCLUDED", cs.excluded),
        ],
    )


def read_combined(path: Path | str) -> tuple[CombinedSpectrum, FrameMeta, Provenance]:
    """The inverse of :func:`write_combined`: ``(combined, meta, provenance)``."""
    with fits.open(path) as hdul:
        header = hdul[0].header
        table = hdul["COMBINED"]
        names = set(table.columns.names)
        versions = {}
        for version in COMBINED_VERSIONS:
            tag = version.upper()
            if f"FLUX_{tag}" in names:
                versions[version] = SpectrumVersion(
                    _native(table.data[f"FLUX_{tag}"]),
                    _native(table.data[f"ERROR_{tag}"]),
                    _native(table.data[f"MASK_{tag}"]).astype(bool),
                )
        cs = CombinedSpectrum(
            target=str(header["MPTARGET"]),
            wavelength=_native(table.data["WAVELENGTH"]),
            versions=versions,
            contributors=_frame_from_table(hdul["CONTRIB"]),
            excluded=_frame_from_table(hdul["EXCLUDED"]),
            grouping=str(header["MPGROUP"]),
            barycentric=_read_mapping(header, BARYCENTRIC_SCHEME) or {},
        )
        return cs, _read_meta(header), Provenance.from_header(header)


# ---------------------------------------------------------------------------------- qc


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"QC value of type {type(value).__name__} is not JSON-serialisable")


def write_qc(path: Path | str, report: QCReport) -> Path:
    """``{"stage": ..., "values": {...}}`` as JSON (``indent=1``, sorted keys, numpy scalars
    and arrays converted)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {"stage": report.stage, "values": report.values}
    text = json.dumps(payload, indent=1, sort_keys=True, default=_jsonable)
    path.write_text(text + "\n", encoding="utf-8")
    return path


def read_qc(path: Path | str) -> QCReport:
    """The inverse of :func:`write_qc`."""
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return QCReport(stage=str(payload["stage"]), values=dict(payload["values"]))
