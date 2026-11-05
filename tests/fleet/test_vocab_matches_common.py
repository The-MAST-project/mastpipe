"""Check :mod:`mastpipe.vocab` against the fleet's shared library ``MAST_common``.

``MAST_common`` is imported as the package ``common``: check it out as a folder
named ``common`` next to this repository and put their parent folder on
``PYTHONPATH``. Without it these tests skip; with ``MASTPIPE_FLEET_CHECK=1`` (as in
CI's ``fleet`` job) a missing ``common`` is an error instead, so the check cannot
pass having proved nothing.
"""

import importlib
import os
import re
import typing
from pathlib import Path

import pytest

from mastpipe.vocab import (
    DEEPSPEC_BANDS,
    DISPERSERS,
    FRAME_TYPES,
    GRATING_NAMES,
    LAMP_NAMES,
    SPEC_INSTRUMENTS,
    unit_name,
)

REQUIRED = os.environ.get("MASTPIPE_FLEET_CHECK") == "1"
SIBLING_LAYOUT = (
    "MAST_common is not importable as 'common': check it out as a sibling folder "
    "named 'common' and put the parent folder on PYTHONPATH"
)
common_spec = (
    importlib.import_module("common.spec") if REQUIRED else pytest.importorskip("common.spec", reason=SIBLING_LAYOUT)
)


def _literal_values(alias) -> tuple[str, ...]:
    """The values of a ``typing.Literal`` alias, in declaration order."""
    return typing.get_args(alias)


def _config_list(text: str, name: str) -> set[str]:
    """The string items of a ``name = [...]`` list literal found in ``text``."""
    match = re.search(rf"{name} = \[(.*?)\]", text)
    assert match is not None, f"{name} not found in common/config/specs.py"
    return set(match.group(1).replace('"', "").replace(" ", "").split(","))


def test_frame_types_match():
    assert tuple(m.value for m in common_spec.FrameType) == FRAME_TYPES


def test_grating_names_match():
    assert _literal_values(common_spec.GratingNames) == GRATING_NAMES


def test_dispersers_match():
    assert _literal_values(common_spec.Disperser) == DISPERSERS


def test_instruments_match():
    assert _literal_values(common_spec.SpecInstruments) == SPEC_INSTRUMENTS


def test_lamps_match():
    assert _literal_values(common_spec.WheelNames) == LAMP_NAMES


def test_bands_match():
    assert _literal_values(common_spec.DeepspecBands) == DEEPSPEC_BANDS


def test_config_lists_match_as_text():
    # common/config/specs.py keeps its lists as locals inside a validator, and
    # common.config needs a database driver to import, so read the file as text.
    path = Path(common_spec.__file__).with_name("config").joinpath("specs.py")
    text = path.read_text(encoding="utf-8")
    assert _config_list(text, "valid_grating_names") == set(GRATING_NAMES)
    assert _config_list(text, "valid_wheel_names") == set(LAMP_NAMES)


def test_unit_names_are_canonical():
    canonic = importlib.import_module("common.utils").canonic_unit_name
    assert all(canonic(unit_name(n)) == unit_name(n) for n in range(1, 21))
