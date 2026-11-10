"""Tests for the spectrograph vocabulary mirrored from the fleet's shared library."""

import pytest

from mastpipe.vocab import (
    DEEPSPEC_BANDS,
    DISPERSERS,
    FRAME_TYPES,
    GRATING_NAMES,
    LAMP_NAMES,
    SPEC_INSTRUMENTS,
    is_unit_name,
    unit_name,
    unit_number,
)


def test_unit_name_is_two_digit():
    assert unit_name(1) == "mast01" and unit_name(20) == "mast20"


@pytest.mark.parametrize("n", [0, 21])
def test_unit_name_rejects_out_of_range(n):
    with pytest.raises(ValueError):
        unit_name(n)


def test_unit_number_round_trips():
    assert all(unit_number(unit_name(n)) == n for n in range(1, 21))
    assert not is_unit_name("T01")


def test_vocabulary_values():
    assert GRATING_NAMES == ("Ca", "Halpha", "Mg", "Future")
    assert FRAME_TYPES == ("light", "bias", "dark", "flat")
    assert LAMP_NAMES == ("ThAr", "qTh")
    assert DISPERSERS == ("Ca", "Mg", "Halpha", "Empty")
    assert SPEC_INSTRUMENTS == ("deepspec", "highspec")
    assert DEEPSPEC_BANDS == ("I", "G", "R", "U")


def test_unit_digits_are_ascii():
    """``\\d`` would accept any Unicode digit (``mast\\u0660\\u0661`` is Arabic-Indic 01)."""
    assert not is_unit_name("mast\u0660\u0661")
    with pytest.raises(ValueError):
        unit_number("mast\u0660\u0661")
    assert not is_unit_name("mast1") and not is_unit_name("mast001") and not is_unit_name("Mast01")


def test_the_project_prefix_is_an_argument_of_all_three():
    """``unit_name``, ``is_unit_name`` and ``unit_number`` take the same ``project``
    argument, with the same default."""
    assert unit_name(7, "lab") == "lab07"
    assert is_unit_name("lab07", "lab") and not is_unit_name("lab07")
    assert unit_number("lab07", "lab") == 7
    with pytest.raises(ValueError, match="lab01"):
        unit_number("mast07", "lab")
    assert all(unit_number(unit_name(n, "lab"), "lab") == n for n in range(1, 21))
