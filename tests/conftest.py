"""Shared pytest configuration for the mastpipe test suite.

The session-scoped simulated nights are built once per session (seed 0) and shared by
every test that needs frames on disk; tests must treat them as read-only.
"""

from pathlib import Path

import pandas as pd
import pytest

from mastpipe.config.loader import build_effective
from mastpipe.config.models import EffectiveConfig
from mastpipe.instruments.registry import load_instrument, load_layout, load_line_table
from mastpipe.sim.simulator import SimTruth, simulate_night


@pytest.fixture(scope="session")
def highspec_effective() -> EffectiveConfig:
    return build_effective(load_instrument("highspec"), "mgb518", "main")


@pytest.fixture(scope="session")
def survey_layout() -> pd.DataFrame:
    return load_layout("highspec", "survey45")


@pytest.fixture(scope="session")
def lab3_layout() -> pd.DataFrame:
    return load_layout("highspec", "lab3")


@pytest.fixture(scope="session")
def mgb518_lines() -> pd.DataFrame:
    return load_line_table("highspec", "mgb518")


@pytest.fixture(scope="session")
def sim_night(
    tmp_path_factory: pytest.TempPathFactory,
    highspec_effective: EffectiveConfig,
    survey_layout: pd.DataFrame,
    mgb518_lines: pd.DataFrame,
) -> tuple[Path, SimTruth]:
    out_dir = tmp_path_factory.mktemp("sim_night")
    truth = simulate_night(
        highspec_effective, survey_layout, mgb518_lines, out_dir, seed=0, n_science=4
    )
    return out_dir, truth


@pytest.fixture(scope="session")
def sim_night_lab3(
    tmp_path_factory: pytest.TempPathFactory,
    highspec_effective: EffectiveConfig,
    lab3_layout: pd.DataFrame,
    mgb518_lines: pd.DataFrame,
) -> tuple[Path, SimTruth]:
    out_dir = tmp_path_factory.mktemp("sim_night_lab3")
    truth = simulate_night(
        highspec_effective,
        lab3_layout,
        mgb518_lines,
        out_dir,
        seed=0,
        n_science=4,
        fiber_preset="lab3",
    )
    return out_dir, truth
