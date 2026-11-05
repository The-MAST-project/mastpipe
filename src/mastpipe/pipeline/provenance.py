"""Provenance: what went into a product, what was applied, and how far it got.

Every product header carries one :class:`Provenance` as FITS cards. Keywords are at most
eight characters, so lists are numbered schemes with a count card: ``MPNIN`` /
``MPIN0001``… (input checksums), ``MPNCAL`` / ``MPCAL001``… (calibrations used), ``MPNDARC``
/ ``MPDARC01``… (drift-arc checksums) and ``MPNSTEP`` / ``MPSTEP01``… (steps applied). The
scalar cards are ``MPVER`` (pipeline version), ``MPCFGH`` (effective-configuration hash),
``MPWSTAT`` / ``MPFSTAT`` / ``MPTSTAT`` (wavelength, flux and telluric status), ``MPWREF``
and ``MPWREFH`` (reference-solution identifier and hash), ``MPLAYOUT`` (the fiber layout
preset the product was reduced with), ``MPLAYH`` (the hash of that preset's CSV file),
``MPMKH`` (the hash of the masters' file checksums) and ``MPTZ`` (the time zone of the
frame's clock, as ingest recorded it: ``utc``, ``assumed-utc`` or an IANA zone name).
Counts are written as integers, everything else as strings; an unset value is the string
``"none"``.

:meth:`Provenance.to_header_cards` and :meth:`Provenance.from_header` are each other's
inverse. Nothing here touches a file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from mastpipe.core.model import FLUX_STATUSES, TELLURIC_STATUSES, WAVELENGTH_STATUSES

NONE = "none"
"""Card value standing for an unset field."""

NUMBERED = (
    ("MPNIN", "MPIN", 4, "inputs", "input checksum"),
    ("MPNCAL", "MPCAL", 3, "calibrations", "calibration used"),
    ("MPNDARC", "MPDARC", 2, "drift_arcs", "drift-arc checksum"),
    ("MPNSTEP", "MPSTEP", 2, "steps", "step applied"),
)
"""``(count card, item prefix, digits, attribute, comment)`` of every numbered scheme."""


def _check_choice(name: str, value: Any, choices: tuple[str, ...]) -> None:
    if value not in choices:
        raise ValueError(f"{name} must be one of {choices}, got {value!r}")


def _to_card(value: str | None) -> str:
    return NONE if value is None else str(value)


def _from_card(value: Any) -> str | None:
    text = str(value)
    return None if text == NONE else text


def _optional(header: Any, key: str) -> str | None:
    """A scalar card added after the first products were written: unset when absent."""
    return _from_card(header[key]) if key in header else None


@dataclass
class Provenance:
    """What a product was made from and what was applied to it.

    ``inputs`` are the checksums of the raw frames, ``calibrations`` the identifiers of the
    masters used, ``drift_arcs`` the checksums of the calibration-fiber arcs whose drift was
    transferred, ``steps`` the operations applied, in order, ``layout`` the fiber
    layout preset (part of a science product's idempotence key) and ``timezone`` how the
    input frames' times were read (``FrameMeta.timezone``). The three statuses follow
    the vocabularies of :mod:`mastpipe.core.model`. Lists are stored as lists, so a tuple
    given here compares equal after a round trip.
    """

    pipeline_version: str
    config_hash: str
    inputs: list[str]
    calibrations: list[str]
    wavelength_status: str = "none"
    flux_status: str = "uncalibrated"
    telluric_status: str = "uncorrected"
    reference_id: str | None = None
    reference_hash: str | None = None
    drift_arcs: list[str] = field(default_factory=list)
    steps: list[str] = field(default_factory=list)
    layout: str | None = None
    timezone: str | None = None
    layout_hash: str | None = None
    masters_hash: str | None = None

    def __post_init__(self) -> None:
        _check_choice("wavelength_status", self.wavelength_status, WAVELENGTH_STATUSES)
        _check_choice("flux_status", self.flux_status, FLUX_STATUSES)
        _check_choice("telluric_status", self.telluric_status, TELLURIC_STATUSES)
        for name in ("inputs", "calibrations", "drift_arcs", "steps"):
            setattr(self, name, [str(item) for item in getattr(self, name)])

    def to_header_cards(self) -> list[tuple[str, str | int, str]]:
        """``(keyword, value, comment)`` triples, counts as integers, the rest as strings.

        ``ValueError`` when a list has more entries than its numbered scheme can hold.
        """
        cards: list[tuple[str, str | int, str]] = [
            ("MPVER", self.pipeline_version, "mastpipe version"),
            ("MPCFGH", self.config_hash, "effective configuration hash"),
            ("MPWSTAT", self.wavelength_status, "wavelength status"),
            ("MPFSTAT", self.flux_status, "flux status"),
            ("MPTSTAT", self.telluric_status, "telluric status"),
            ("MPWREF", _to_card(self.reference_id), "reference solution identifier"),
            ("MPWREFH", _to_card(self.reference_hash), "reference solution hash"),
            ("MPLAYOUT", _to_card(self.layout), "fiber layout preset"),
            ("MPTZ", _to_card(self.timezone), "time zone of the frame's clock"),
            ("MPLAYH", _to_card(self.layout_hash), "fiber layout file hash"),
            ("MPMKH", _to_card(self.masters_hash), "hash of the masters' checksums"),
        ]
        for count_key, prefix, digits, attribute, comment in NUMBERED:
            items = getattr(self, attribute)
            if len(items) >= 10**digits:
                raise ValueError(f"{attribute}: {len(items)} entries exceed the {prefix} cards")
            cards.append((count_key, len(items), f"number of {comment}s"))
            for n, item in enumerate(items, 1):
                cards.append((f"{prefix}{n:0{digits}d}", item, f"{comment} {n}"))
        return cards

    @classmethod
    def from_header(cls, header: Any) -> Provenance:
        """The inverse of :meth:`to_header_cards` on a FITS header (or any mapping).

        ``KeyError`` when a scalar card or a numbered card the count announces is missing;
        ``MPLAYOUT``, ``MPTZ``, ``MPLAYH`` and ``MPMKH``, added later than the others, read
        as unset when absent.
        """
        lists: dict[str, list[str]] = {}
        for count_key, prefix, digits, attribute, _ in NUMBERED:
            count = int(header[count_key])
            lists[attribute] = [str(header[f"{prefix}{n:0{digits}d}"]) for n in range(1, count + 1)]
        return cls(
            pipeline_version=str(header["MPVER"]),
            config_hash=str(header["MPCFGH"]),
            inputs=lists["inputs"],
            calibrations=lists["calibrations"],
            wavelength_status=str(header["MPWSTAT"]),
            flux_status=str(header["MPFSTAT"]),
            telluric_status=str(header["MPTSTAT"]),
            reference_id=_from_card(header["MPWREF"]),
            reference_hash=_from_card(header["MPWREFH"]),
            drift_arcs=lists["drift_arcs"],
            steps=lists["steps"],
            layout=_optional(header, "MPLAYOUT"),
            timezone=_optional(header, "MPTZ"),
            layout_hash=_optional(header, "MPLAYH"),
            masters_hash=_optional(header, "MPMKH"),
        )
