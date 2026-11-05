"""Tests for the spectrograph vocabulary mirrored from the fleet's shared library."""

import pytest

from mastpipe.vocab import (
    FRAME_TYPES,
    GRATING_NAMES,
    LAMP_NAMES,
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
