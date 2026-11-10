"""Shared pytest configuration for the mastpipe test suite.

The session-scoped simulated nights are built once per session (seed 0) and shared by
every test that needs frames on disk; tests must treat them as read-only. The calibrated
full arc of the survey night, its truth-built traces, its line tables and its
continuum-subtracted fluxes are session fixtures as well, shared by the line-finding and
wavelength-solution tests.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.loader import build_effective
from mastpipe.config.models import CalibParams, EffectiveConfig, ExtractParams, LineParams
from mastpipe.core.model import Frame, Trace, TraceSet
from mastpipe.instruments.registry import (
    load_header_map,
    load_instrument,
    load_layout,
    load_line_table,
)
from mastpipe.io.headers import HeaderMap
from mastpipe.io.ingest import load_frame, scan_folder
from mastpipe.sim.simulator import FWHM_PX, GAUSSIAN_FWHM, SimTruth, simulate_night
from mastpipe.steps.calib import build_calibration_set, build_master, calibrate_frame
from mastpipe.steps.extract import continuum_subtract_all, extract_aperture
from mastpipe.steps.geometry import GratingGeometry, grating_geometry
from mastpipe.steps.lines import find_lines
from mastpipe.steps.wavemodel import GratingModel

ARC_STRIP_HALFWIDTH = 4
"""Rows on either side of a simulated trace (3 px FWHM image, 8.5 px pitch): the strip holds
the whole profile and stays 4.5 rows clear of the neighbouring fiber's centre."""
ARC_SIGMA_PX = FWHM_PX / GAUSSIAN_FWHM
"""Gaussian sigma of the simulated line image, px."""


@pytest.fixture(scope="session")
def highspec_effective() -> EffectiveConfig:
    return build_effective(load_instrument("highspec"), "Mg", "main")


@pytest.fixture
def highspec_map() -> HeaderMap:
    """HighSpec's shipped header map, parsed."""
    return HeaderMap.from_dict(load_header_map("highspec"))


@pytest.fixture(scope="session")
def survey_layout() -> pd.DataFrame:
    return load_layout("highspec", "survey45")


@pytest.fixture(scope="session")
def lab3_layout() -> pd.DataFrame:
    return load_layout("highspec", "lab3")


@pytest.fixture(scope="session")
def mg_lines() -> pd.DataFrame:
    return load_line_table("highspec", "Mg")


@pytest.fixture(scope="session")
def sim_night(
    tmp_path_factory: pytest.TempPathFactory,
    highspec_effective: EffectiveConfig,
    survey_layout: pd.DataFrame,
    mg_lines: pd.DataFrame,
) -> tuple[Path, SimTruth]:
    out_dir = tmp_path_factory.mktemp("sim_night")
    truth = simulate_night(highspec_effective, survey_layout, mg_lines, out_dir, seed=0, n_science=4)
    return out_dir, truth


@pytest.fixture(scope="session")
def sim_night_lab3(
    tmp_path_factory: pytest.TempPathFactory,
    highspec_effective: EffectiveConfig,
    lab3_layout: pd.DataFrame,
    mg_lines: pd.DataFrame,
) -> tuple[Path, SimTruth]:
    out_dir = tmp_path_factory.mktemp("sim_night_lab3")
    truth = simulate_night(
        highspec_effective,
        lab3_layout,
        mg_lines,
        out_dir,
        seed=0,
        n_science=4,
        fiber_preset="lab3",
    )
    return out_dir, truth


# ------------------------------------------------------------ the calibrated full arc


def truth_model(truth: SimTruth, fiber: int, geom: GratingGeometry, n_pix: int) -> GratingModel:
    """The grating model the simulator drew fiber ``fiber`` (a layout trace id) from: every
    trace has its own ``(x_c, B)``."""
    row = truth.row_of(fiber)
    x_c, B = truth.x_c_per_trace[row], truth.B_per_trace[row]
    return GratingModel(float(x_c), float(B), geom.alpha_rad, None, n_pix, geom)


def truth_traceset(truth: SimTruth, shape: tuple[int, int]) -> TraceSet:
    """Traces built from the simulator's polynomial coefficients, with the layout's roles and
    telescopes (``None`` for the arc fibers)."""
    traces = [
        Trace(
            index=i,
            coeffs=np.asarray(coeffs, dtype=float),
            y_center=float(np.polyval(coeffs, shape[1] // 2)),
            slope=0.0,
            rms=0.0,
            n_points=10,
            snr=100.0,
            x_samples=np.arange(10.0),
            y_samples=np.full(10, float(np.polyval(coeffs, shape[1] // 2))),
            kept=np.ones(10, dtype=bool),
            role=row["role"],
            telescope=None if isinstance(row["telescope"], float) else row["telescope"],
        )
        for i, (coeffs, (_, row)) in enumerate(zip(truth.trace_coeffs, truth.layout.iterrows(), strict=True))
    ]
    return TraceSet(traces=traces, shape=shape)


@pytest.fixture(scope="session")
def sim_full_arc(sim_night, highspec_effective) -> SimpleNamespace:
    """The simulated full arc, calibrated, with traces built from the truth.

    ``truth``, ``arc`` (calibrated :class:`Frame`), ``traces``, ``params`` (the
    :class:`LineParams` of the simulated line image), ``geom`` and ``n_pix``.
    """
    out_dir, truth = sim_night
    detector = highspec_effective.detector
    refs = scan_folder(out_dir, "highspec", load_header_map("highspec"))
    by_type: dict[str, list[Frame]] = {}
    for ref in refs:
        by_type.setdefault(ref.meta.exposure_type, []).append(load_frame(ref, detector))
    params = CalibParams()
    masters = [build_master(by_type[kind], kind, params) for kind in ("bias", "dark", "background")]
    calib = build_calibration_set(*masters, detector)
    kind_of = dict(zip(truth.frames["file"], truth.frames["arc_kind"], strict=True))
    full = next(frame for frame in by_type["arc"] if kind_of[frame.meta.source] == "full")
    arc = calibrate_frame(full, calib, saturation=detector.saturation)
    return SimpleNamespace(
        truth=truth,
        arc=arc,
        traces=truth_traceset(truth, arc.shape),
        params=LineParams(halfwidth=ARC_STRIP_HALFWIDTH, sigma_guess=ARC_SIGMA_PX),
        geom=grating_geometry(highspec_effective.geometry, detector.pixel_size_um),
        n_pix=arc.shape[1],
    )


@pytest.fixture(scope="session")
def sim_arc_lines(sim_full_arc) -> tuple[list[pd.DataFrame], list[np.ndarray]]:
    """Line tables and strips of every trace of the simulated full arc."""
    return find_lines(sim_full_arc.arc.data, sim_full_arc.arc.variance, sim_full_arc.traces, sim_full_arc.params)


@pytest.fixture(scope="session")
def sim_arc_fluxes(sim_full_arc) -> list[np.ndarray]:
    """Continuum-subtracted extracted fluxes of every trace of the simulated full arc."""
    spectra = extract_aperture(sim_full_arc.arc, sim_full_arc.traces, ExtractParams(halfwidth=ARC_STRIP_HALFWIDTH))
    spectra = continuum_subtract_all(spectra, ExtractParams())
    return [spectrum.flux - spectrum.continuum for spectrum in spectra]
