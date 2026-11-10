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
- :data:`UNIT_COUNT`, :data:`PROJECT` and :func:`unit_name` -- the unit-naming rule of
  ``common/config/site.py`` (project name followed by the two-digit unit number);
  :func:`is_unit_name` and :func:`unit_number` read it back, all three with the same
  optional ``project`` argument.

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

PROJECT = "mast"
"""The project prefix of the unit names (``common/config/site.py``'s project name)."""


def _unit_pattern(project: str) -> re.Pattern[str]:
    # [0-9], not \d: \d accepts any Unicode digit
    return re.compile(re.escape(project) + r"([0-9]{2})")


def unit_name(number: int, project: str = PROJECT) -> str:
    """Return the canonical name of a telescope unit, e.g. ``mast01``.

    Raises ``ValueError`` when ``number`` is outside ``1..UNIT_COUNT``.
    """
    if not 1 <= number <= UNIT_COUNT:
        raise ValueError(f"unit number {number} is outside 1..{UNIT_COUNT}")
    return f"{project}{number:02d}"


def is_unit_name(text: str, project: str = PROJECT) -> bool:
    """Whether ``text`` is a canonical unit name of ``project``, ``mast01`` to ``mast20``
    (two ASCII digits)."""
    match = _unit_pattern(project).fullmatch(text)
    return match is not None and 1 <= int(match.group(1)) <= UNIT_COUNT


def unit_number(name: str, project: str = PROJECT) -> int:
    """Return the number of a canonical unit name of ``project`` (``mast07`` gives 7).

    Raises ``ValueError`` when ``name`` is not a canonical unit name.
    """
    if not is_unit_name(name, project):
        raise ValueError(f"{name!r} is not a unit name ({project}01..{project}{UNIT_COUNT:02d})")
    return int(name[len(project) :])
