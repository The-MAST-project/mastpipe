"""The spectrograph vocabulary of the MAST fleet, mirrored from ``MAST_common``.

The fleet's shared library ``MAST_common`` (imported as the package ``common``)
defines the names the spectrograph control software writes into frames and
configuration. ``mastpipe`` must speak the same vocabulary but cannot depend on
``common``, which needs the observatory's environment to import, so this module
keeps a standard-library-only copy. A fleet test (``tests/fleet``) checks every
constant against the real definitions; CI runs it with ``MAST_common`` checked out
next to this repository.

Each constant mirrors one name in ``MAST_common``:

- :data:`FRAME_TYPES` -- the values of ``common/spec.py: FrameType``.
- :data:`GRATING_NAMES` -- ``common/spec.py: GratingNames``.
- :data:`DISPERSERS` -- ``common/spec.py: Disperser``.
- :data:`SPEC_INSTRUMENTS` -- ``common/spec.py: SpecInstruments``.
- :data:`LAMP_NAMES` -- ``common/spec.py: WheelNames``.
- :data:`DEEPSPEC_BANDS` -- ``common/spec.py: DeepspecBands``.
- :data:`UNIT_COUNT` and :func:`unit_name` -- the unit-naming rule of
  ``common/config/site.py`` (project name followed by the two-digit unit number).

``mastpipe``'s own exposure types (``arc``, ``background``, ``unknown``,
``science``) are pipeline-only and live in :mod:`mastpipe.core.model`.
"""

from __future__ import annotations

import re

FRAME_TYPES: tuple[str, ...] = ("light", "bias", "dark", "flat")
GRATING_NAMES: tuple[str, ...] = ("Ca", "Halpha", "Mg", "Future")
DISPERSERS: tuple[str, ...] = ("Ca", "Mg", "Halpha", "Empty")
SPEC_INSTRUMENTS: tuple[str, ...] = ("deepspec", "highspec")
LAMP_NAMES: tuple[str, ...] = ("ThAr", "qTh")
DEEPSPEC_BANDS: tuple[str, ...] = ("I", "G", "R", "U")
UNIT_COUNT = 20

_UNIT_NAME = re.compile(r"mast(\d{2})")


def unit_name(number: int, project: str = "mast") -> str:
    """Return the canonical name of a telescope unit, e.g. ``mast01``.

    Raises ``ValueError`` when ``number`` is outside ``1..UNIT_COUNT``.
    """
    if not 1 <= number <= UNIT_COUNT:
        raise ValueError(f"unit number {number} is outside 1..{UNIT_COUNT}")
    return f"{project}{number:02d}"


def is_unit_name(text: str) -> bool:
    """Whether ``text`` is a canonical unit name, ``mast01`` to ``mast20``."""
    match = _UNIT_NAME.fullmatch(text)
    return match is not None and 1 <= int(match.group(1)) <= UNIT_COUNT


def unit_number(name: str) -> int:
    """Return the number of a canonical unit name (``mast07`` gives 7).

    Raises ``ValueError`` when ``name`` is not a canonical unit name.
    """
    if not is_unit_name(name):
        raise ValueError(f"{name!r} is not a unit name (mast01..mast{UNIT_COUNT:02d})")
    return int(name[len("mast") :])
