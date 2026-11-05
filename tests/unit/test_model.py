"""Tests for the core data model."""

import numpy as np
import pandas as pd
import pytest

from mastpipe.core.model import (
    VERSIONS,
    CombinedSpectrum,
    DriftModel,
    ExposureSpectra,
    Frame,
    FrameMeta,
    Spectrum,
    SpectrumVersion,
    Trace,
    TraceSet,
    WavelengthSolution,
    make_identifier,
)


def make_meta(**overrides) -> FrameMeta:
    fields = {"instrument": "highspec", "camera": "main", "exposure_type": "science"}
    fields.update(overrides)
    return FrameMeta(**fields)


def make_trace(index, role=None, telescope=None, y_center=0.0) -> Trace:
    x_samples = np.arange(5.0)
    return Trace(
        index=index,
        coeffs=np.array([0.0, y_center]),
        y_center=y_center,
        slope=0.0,
        rms=0.1,
        n_points=5,
        snr=50.0,
        x_samples=x_samples,
        y_samples=np.full(5, y_center),
        kept=np.ones(5, dtype=bool),
        role=role,
        telescope=telescope,
    )


def make_spectrum(n=6, **overrides) -> Spectrum:
    fields = {
        "fiber": 0,
        "pixel": np.arange(n),
        "flux": np.ones(n),
        "error": np.full(n, 0.1),
        "mask": np.zeros(n, dtype=bool),
    }
    fields.update(overrides)
    return Spectrum(**fields)


def make_version(n=6, level=2.0) -> SpectrumVersion:
    return SpectrumVersion(np.full(n, level), np.full(n, 0.1), np.zeros(n, dtype=bool))


def test_frame_meta_rejects_unknown_exposure_type():
    with pytest.raises(ValueError, match="exposure_type"):
        make_meta(exposure_type="banana")
    assert make_meta(exposure_type="bias").exposure_type == "bias"


def test_frame_meta_rejects_unknown_arc_kind():
    with pytest.raises(ValueError, match="arc_kind"):
        make_meta(exposure_type="arc", arc_kind="banana")
    meta = make_meta(exposure_type="arc", arc_kind="calfiber")
    assert meta.arc_kind == "calfiber"
    assert meta.binning == (1, 1)
    with pytest.raises(AttributeError):
        meta.camera = "other"  # frozen


def test_target_for_prefers_exact_then_star():
    meta = make_meta(targets=(("*", "all"), ("mast02", "two")))
    assert meta.target_for("mast02") == "two"
    assert meta.target_for("mast07") == "all"
    assert meta.target_for(None) == "all"
    assert make_meta().target_for("mast01") is None
    only = make_meta(targets=(("mast01", "one"),))
    assert only.target_for("mast01") == "one" and only.target_for("mast02") is None
    assert only.target_for(None) is None


def test_frame_meta_rejects_malformed_targets():
    with pytest.raises(ValueError, match="targets"):
        make_meta(targets=(("mast01",),))
    with pytest.raises(ValueError, match="targets"):
        make_meta(targets=(("mast01", 3),))


def test_frame_rejects_shape_mismatch():
    data = np.zeros((4, 5), dtype=np.float32)
    with pytest.raises(ValueError, match="variance"):
        Frame(data, make_meta(), variance=np.zeros((4, 4)))
    with pytest.raises(ValueError, match="mask"):
        Frame(data, make_meta(), mask=np.zeros((5, 4), dtype=bool))
    with pytest.raises(ValueError, match="bool"):
        Frame(data, make_meta(), mask=np.zeros((4, 5), dtype=int))
    frame = Frame(data, make_meta(), mask=np.zeros((4, 5), dtype=bool))
    assert frame.shape == (4, 5)
    assert frame.variance is None


def test_trace_evaluates_polynomial():
    trace = make_trace(0)
    trace.coeffs = np.array([0.001, 0.5, 10.0])
    assert trace.y(2.0) == pytest.approx(11.004)
    np.testing.assert_allclose(trace.y(np.array([0.0, 2.0])), [10.0, 11.004])


def test_traceset_by_role_and_telescope():
    traces = [
        make_trace(0, role="arc", y_center=10.0),
        make_trace(1, role="science", telescope="mast01", y_center=20.0),
        make_trace(2, role="sky", telescope="mast01", y_center=30.0),
        make_trace(3, role="science", telescope="mast02", y_center=40.0),
    ]
    traceset = TraceSet(traces=traces, shape=(100, 50))
    assert len(traceset) == 4
    assert [trace.index for trace in traceset] == [0, 1, 2, 3]
    assert [trace.index for trace in traceset.by_role("science")] == [1, 3]
    assert traceset.by_role("unknown") == []
    by_t1 = traceset.by_telescope("mast01")
    assert set(by_t1) == {"science", "sky"}
    assert by_t1["sky"].index == 2
    assert traceset.by_telescope("mast09") == {}
    assert traceset.missing == () and traceset.extra == ()


def test_spectrum_versions_validate_and_latest():
    assert VERSIONS == ("skysub", "fluxcal", "telluric")
    raw = make_spectrum()
    latest = raw.latest(upto=None)
    assert isinstance(latest, SpectrumVersion)
    assert latest.flux is raw.flux and latest.mask is raw.mask
    assert raw.versions == {}

    skysub = make_version(level=2.0)
    fluxcal = make_version(level=3.0)
    one = make_spectrum(versions={"skysub": skysub})
    assert one.latest() is skysub
    assert one.latest(upto="skysub") is skysub
    assert one.latest(upto="telluric") is skysub
    two = make_spectrum(versions={"skysub": skysub, "fluxcal": fluxcal})
    assert two.latest() is fluxcal
    assert two.latest(upto="skysub") is skysub
    only_fluxcal = make_spectrum(versions={"fluxcal": fluxcal})
    assert only_fluxcal.latest(upto="skysub").flux is only_fluxcal.flux

    with pytest.raises(ValueError, match="banana"):
        make_spectrum(versions={"banana": skysub})
    with pytest.raises(ValueError, match="skysub"):
        make_spectrum(versions={"skysub": make_version(n=7)})
    with pytest.raises(ValueError, match="banana"):
        one.latest(upto="banana")
    with pytest.raises(ValueError, match="flux"):
        make_spectrum(flux=np.ones(7))
    with pytest.raises(ValueError, match="error"):
        SpectrumVersion(np.ones(6), np.ones(5), np.zeros(6, dtype=bool))


def test_exposure_spectra_statuses():
    exposure = ExposureSpectra(meta=make_meta(), spectra=[make_spectrum()])
    assert exposure.wavelength_status == "none"
    assert exposure.flux_status == "uncalibrated"
    assert exposure.telluric_status == "uncorrected"
    assert exposure.drift_arcs == ()
    with pytest.raises(ValueError, match="flux_status"):
        ExposureSpectra(meta=make_meta(), spectra=[], flux_status="banana")
    with pytest.raises(ValueError, match="wavelength_status"):
        ExposureSpectra(meta=make_meta(), spectra=[], wavelength_status="banana")
    with pytest.raises(ValueError, match="telluric_status"):
        ExposureSpectra(meta=make_meta(), spectra=[], telluric_status="banana")


def test_reference_identifier_format():
    assert make_identifier("highspec", "Mg", "lab3", "2026-09-09T12:34:56") == "highspec_Mg_lab3_20260909T123456"


def make_per_fiber(**overrides) -> pd.DataFrame:
    table = {
        "fiber": [0, 1],
        "y_center": [10.0, 90.0],
        "shift_px": [0.5, 1.3],
        "stretch": [0.0, 1e-4],
        "n_lines": [20, 18],
        "rms_px": [0.05, 0.06],
        "usable": [True, True],
    }
    table.update(overrides)
    return pd.DataFrame(table)


def make_drift(per_fiber=None) -> DriftModel:
    return DriftModel(
        arc_checksums=("0123456789abcdef",),
        arc_mid_time="2026-09-09T12:00:00",
        per_fiber=make_per_fiber() if per_fiber is None else per_fiber,
        shift_coeffs=np.array([0.01, 0.4]),
        stretch_coeffs=np.array([1e-4]),
        scatter_px=0.02,
        usable=True,
        n_pix=100,
    )


def test_drift_model_evaluates_at_y():
    model = make_drift()
    assert model.shift_at(10.0) == pytest.approx(0.5)
    np.testing.assert_allclose(model.shift_at(np.array([0.0, 90.0])), [0.4, 1.3])
    assert model.stretch_at(50.0) == pytest.approx(1e-4)
    np.testing.assert_allclose(model.stretch_at(np.array([0.0, 90.0])), [1e-4, 1e-4])
    with pytest.raises(ValueError, match="per_fiber"):
        make_drift(per_fiber=make_per_fiber().drop(columns="usable"))


def test_tables_require_named_columns():
    lines = pd.DataFrame({"x": [10.0], "wavelength": [518.36], "species": ["Th"], "dx_px": [0.0]})
    solution_fields = {
        "fiber": 0,
        "kind": "grating",
        "params": np.zeros(4),
        "n_pix": 100,
        "rms_px": 0.05,
        "n_lines": 1,
        "anchor_px": (10.0, 90.0),
        "anchor_nm": (517.0, 519.0),
        "trace_y_center": 42.0,
        "geometry": {},
    }
    assert WavelengthSolution(lines=lines, **solution_fields).lines is lines
    with pytest.raises(ValueError, match="lines"):
        WavelengthSolution(lines=lines.drop(columns="species"), **solution_fields)

    wavelength = np.linspace(517.0, 519.0, 6)
    contributors = pd.DataFrame({"fiber": [1], "telescope": ["mast01"], "source": ["a.fits"], "weight": [1.0]})
    excluded = pd.DataFrame({"fiber": [3], "telescope": ["mast02"], "source": ["a.fits"], "reason": ["low snr"]})
    combined_fields = {
        "target": "star",
        "wavelength": wavelength,
        "versions": {"raw": make_version(), "skysub": make_version()},
        "grouping": "exposure",
        "barycentric": {"mast01": 0.1},
    }
    combined = CombinedSpectrum(contributors=contributors, excluded=excluded, **combined_fields)
    assert set(combined.versions) == {"raw", "skysub"}
    with pytest.raises(ValueError, match="contributors"):
        CombinedSpectrum(contributors=contributors.drop(columns="weight"), excluded=excluded, **combined_fields)
    with pytest.raises(ValueError, match="excluded"):
        CombinedSpectrum(contributors=contributors, excluded=excluded.drop(columns="reason"), **combined_fields)
    with pytest.raises(ValueError, match="banana"):
        CombinedSpectrum(
            contributors=contributors,
            excluded=excluded,
            **{**combined_fields, "versions": {"banana": make_version()}},
        )
    with pytest.raises(ValueError, match="raw"):
        CombinedSpectrum(
            contributors=contributors,
            excluded=excluded,
            **{**combined_fields, "versions": {"raw": make_version(n=7)}},
        )


def test_readout_key_joins_mode_and_amplifier():
    assert make_meta(readout_mode="image", amplifier="em").readout_key == "image/em"
    assert make_meta(readout_mode="image").readout_key == "image/unknown"
    assert make_meta().readout_key == "unknown/unknown"
