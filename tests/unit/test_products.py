"""Tests for the FITS products and their provenance: every product kind is written to a
temporary folder and read back, and every array and scalar field is compared. The
reference round trip also rebuilds the grating model and checks the solutions hash; the
last test opens every product and asserts that no raw instrument keyword leaked into any
header."""

import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.core.model import (
    LINE_COLUMNS,
    CombinedSpectrum,
    DriftModel,
    ExposureSpectra,
    FrameMeta,
    MasterCalib,
    QCReport,
    Spectrum,
    SpectrumVersion,
    Trace,
    TraceSet,
    WavelengthSolution,
)
from mastpipe.io.products import (
    combined_name,
    product_name,
    read_combined,
    read_master,
    read_qc,
    read_reference,
    read_spectra,
    read_traceset,
    reference_name,
    write_combined,
    write_master,
    write_qc,
    write_reference,
    write_spectra,
    write_traceset,
)
from mastpipe.pipeline.provenance import Provenance
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.wavecal import build_reference, solutions_hash
from mastpipe.steps.wavemodel import GratingModel, model_from_config, model_from_solution

RAW_KEYWORDS = ("EXPOSURE", "IMAGETYP", "OBJECT", "DATE-OBS", "LAMP")
N_PIX = 1600
N_SPEC = 50
CHECKSUMS = ("0123456789abcdef", "fedcba9876543210", "00ff00ff00ff00ff")

META = FrameMeta(
    instrument="highspec",
    camera="main",
    exposure_type="science",
    configuration="Mg",
    exposure_time=300.0,
    mid_time="2026-10-01T20:30:15",
    binning=(1, 1),
    readout_mode="fast",
    detector_temperature=-60.5,
    lamp_state="off",
    fiber_preset="survey45",
    targets=(("mast01", "HD 1234"), ("*", "target_b")),
    arc_kind=None,
    checksum=CHECKSUMS[0],
    source="sim_science_01.fits",
    extension=0,
    unresolved=("lamp_state",),
    classification_rule="header",
)


def make_provenance() -> Provenance:
    return Provenance(
        pipeline_version="0.0.1",
        config_hash="abcdef012345",
        inputs=[CHECKSUMS[0], CHECKSUMS[1]],
        calibrations=["bias_a", "dark_b", "background_c"],
        wavelength_status="transferred",
        flux_status="uncalibrated",
        telluric_status="uncorrected",
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
        drift_arcs=[CHECKSUMS[1], CHECKSUMS[2]],
        steps=["bias", "background", "extract"],
        layout="survey45",
        timezone="assumed-utc",
    )


def make_spectrum(fiber: int, rng: np.random.Generator, with_versions: bool) -> Spectrum:
    flux = rng.normal(100.0, 5.0, N_SPEC)
    error = np.sqrt(np.abs(flux))
    mask = np.zeros(N_SPEC, dtype=bool)
    mask[fiber] = True
    versions = {}
    if with_versions:
        versions["skysub"] = SpectrumVersion(
            flux - 20.0, error * 1.1, mask | (np.arange(N_SPEC) == 7)
        )
    return Spectrum(
        fiber=fiber,
        pixel=np.arange(N_SPEC),
        flux=flux,
        error=error,
        mask=mask,
        role=("science", "sky", None)[fiber % 3],
        telescope=("mast01", "mast01", None)[fiber % 3],
        target=("HD 1234", None, None)[fiber % 3],
        wavelength=517.0 + 0.00866 * np.arange(N_SPEC) if with_versions else None,
        sky=np.full(N_SPEC, 20.0) if with_versions else None,
        continuum=np.full(N_SPEC, 100.0) if with_versions else None,
        versions=versions,
    )


def assert_arrays_equal(actual: np.ndarray | None, expected: np.ndarray | None) -> None:
    if expected is None:
        assert actual is None
    else:
        assert actual is not None
        np.testing.assert_array_equal(actual, expected)


def assert_version_equal(actual: SpectrumVersion, expected: SpectrumVersion) -> None:
    assert_arrays_equal(actual.flux, expected.flux)
    assert_arrays_equal(actual.error, expected.error)
    assert_arrays_equal(actual.mask, expected.mask)
    assert actual.mask.dtype == np.bool_


def assert_spectrum_equal(actual: Spectrum, expected: Spectrum) -> None:
    assert actual.fiber == expected.fiber
    for name in ("role", "telescope", "target"):
        assert getattr(actual, name) == getattr(expected, name), name
    for name in ("pixel", "flux", "error", "mask", "wavelength", "sky", "continuum"):
        assert_arrays_equal(getattr(actual, name), getattr(expected, name))
    assert actual.pixel.dtype.kind == expected.pixel.dtype.kind
    assert set(actual.versions) == set(expected.versions)
    for key, version in expected.versions.items():
        assert_version_equal(actual.versions[key], version)


def all_headers(path) -> list[fits.Header]:
    with fits.open(path) as hdul:
        return [hdu.header.copy() for hdu in hdul]


# ----------------------------------------------------------------------------- provenance


def test_provenance_cards_roundtrip():
    prov = make_provenance()
    cards = prov.to_header_cards()
    header = fits.Header()
    for key, value, comment in cards:
        header[key] = (value, comment)
    keys = [key for key, _, _ in cards]
    for key in ("MPVER", "MPCFGH", "MPWSTAT", "MPFSTAT", "MPTSTAT", "MPWREF", "MPWREFH"):
        assert key in keys
    assert header["MPNIN"] == 2 and header["MPIN0001"] == CHECKSUMS[0]
    assert header["MPNCAL"] == 3 and header["MPCAL003"] == "background_c"
    assert header["MPNDARC"] == 2 and header["MPDARC02"] == CHECKSUMS[2]
    assert header["MPSTEP03"] == "extract"
    assert header["MPWREF"] == prov.reference_id and header["MPWREFH"] == prov.reference_hash
    assert header["MPLAYOUT"] == "survey45"
    assert all(len(key) <= 8 for key in keys)
    assert Provenance.from_header(header) == prov
    # a product written before the layout card existed reads back with no layout
    del header["MPLAYOUT"]
    assert Provenance.from_header(header).layout is None
    assert header["MPTZ"] == "assumed-utc"
    # nor does one written before the time-zone card existed
    del header["MPTZ"]
    assert Provenance.from_header(header).timezone is None

    bare = Provenance(
        pipeline_version="0.0.1", config_hash="abcdef012345", inputs=[], calibrations=[]
    )
    header = fits.Header()
    for key, value, comment in bare.to_header_cards():
        header[key] = (value, comment)
    assert header["MPWREF"] == "none" and header["MPWREFH"] == "none"
    assert header["MPLAYOUT"] == "none"
    assert header["MPTZ"] == "none"
    assert header["MPNIN"] == 0 and "MPIN0001" not in header
    back = Provenance.from_header(header)
    assert back == bare
    assert back.reference_id is None and back.drift_arcs == [] and back.steps == []
    assert back.layout is None and back.timezone is None


def test_product_name_format():
    assert product_name(META, "spec") == "highspec_Mg_main_20261001T203015_012345_spec.fits"
    assert product_name(META, "qc", ext=".json") == (
        "highspec_Mg_main_20261001T203015_012345_qc.json"
    )
    bare = replace(META, configuration=None, mid_time=None, checksum="")
    assert product_name(bare, "trace") == "highspec_noconf_main_notime_nosum_trace.fits"
    aware = replace(META, mid_time="2026-10-01T22:30:15+02:00")
    assert product_name(aware, "spec") == "highspec_Mg_main_20261001T203015_012345_spec.fits"
    assert combined_name(META, "HD 1234/b") == (
        "highspec_Mg_main_20261001T203015_012345_comb_HD_1234_b.fits"
    )


def test_product_name_differs_for_equal_times():
    meta = replace(META, mid_time="2026-10-01T18:02:30")
    a = replace(meta, checksum="abcdef0123456789")
    b = replace(meta, checksum="123456abcdef0000")
    assert product_name(a, "spec") != product_name(b, "spec")
    assert product_name(a, "spec") == "highspec_Mg_main_20261001T180230_abcdef_spec.fits"
    undated = replace(meta, mid_time=None)
    assert product_name(replace(undated, checksum=a.checksum), "spec") != product_name(
        replace(undated, checksum=b.checksum), "spec"
    )


# --------------------------------------------------------------------------------- masters


def test_master_roundtrip(tmp_path):
    rng = np.random.default_rng(1)
    data = rng.normal(500.0, 3.0, (8, 12))
    variance = rng.uniform(1.0, 2.0, (8, 12))
    validity = {
        "camera": "main",
        "binning": (1, 1),
        "readout_mode": "fast",
        "detector_temperature": None,
        "mid_time": "2026-10-01T18:05:00",
        "exposure_time": 0.0,
    }
    master = MasterCalib(
        kind="bias",
        data=data,
        validity=validity,
        n_frames=5,
        provenance={"inputs": list(CHECKSUMS), "n_frames": 5, "method": "median"},
        variance=variance,
    )
    prov = make_provenance()
    path = tmp_path / "bias.fits"
    write_master(path, master, prov)
    back, prov_back = read_master(path)
    assert_arrays_equal(back.data, data)
    assert_arrays_equal(back.variance, variance)
    assert back.kind == "bias" and back.n_frames == 5
    assert back.validity == validity
    assert back.provenance == master.provenance
    assert prov_back == prov
    header = fits.getheader(path)
    assert header["MPKIND"] == "bias" and header["MPNFRM"] == 5
    assert header["MPVALK01"] == "camera" and json.loads(header["MPVALV01"]) == "main"

    plain = MasterCalib(kind="dark", data=data, validity=validity, n_frames=3)
    write_master(tmp_path / "dark.fits", plain, prov)
    back, _ = read_master(tmp_path / "dark.fits")
    assert back.variance is None and back.provenance == {}
    with fits.open(tmp_path / "dark.fits") as hdul:
        assert "VARIANCE" not in hdul


# ---------------------------------------------------------------------------------- traces


def make_trace(index: int, coeffs: list[float], **kwargs) -> Trace:
    x = np.linspace(0.0, N_PIX - 1, 7)
    kept = np.ones(7, dtype=bool)
    kept[index % 7] = False
    return Trace(
        index=index,
        coeffs=np.asarray(coeffs),
        y_center=float(np.polyval(coeffs, N_PIX / 2)),
        slope=float(coeffs[-2]),
        rms=0.01 * (index + 1),
        n_points=7,
        snr=50.0 + index,
        x_samples=x,
        y_samples=np.polyval(coeffs, x),
        kept=kept,
        **kwargs,
    )


def test_traceset_roundtrip_degrees(tmp_path):
    traces = [
        make_trace(0, [0.002, 10.0], role="arc", telescope=None, layout_index=0, aperture=None),
        make_trace(1, [1e-6, 0.002, 20.0], role="science", telescope="mast01", layout_index=1),
        make_trace(2, [1e-9, 1e-6, 0.002, 30.0], aperture=3.5),
    ]
    traceset = TraceSet(traces, shape=(400, N_PIX), layout_shift=0.37, missing=(3,), extra=(44,))
    prov = make_provenance()
    path = tmp_path / "traces.fits"
    write_traceset(path, traceset, prov)
    back, prov_back = read_traceset(path)
    assert prov_back == prov
    assert back.shape == (400, N_PIX)
    assert back.layout_shift == pytest.approx(0.37)
    assert back.missing == (3,) and back.extra == (44,)
    assert len(back) == 3
    for actual, expected in zip(back.traces, traces, strict=True):
        assert actual.index == expected.index
        assert_arrays_equal(actual.coeffs, expected.coeffs)
        assert len(actual.coeffs) == len(expected.coeffs)
        for name in ("y_center", "slope", "rms", "snr", "aperture"):
            assert getattr(actual, name) == pytest.approx(getattr(expected, name)), name
        for name in ("n_points", "role", "telescope", "layout_index"):
            assert getattr(actual, name) == getattr(expected, name), name
        assert_arrays_equal(actual.x_samples, expected.x_samples)
        assert_arrays_equal(actual.y_samples, expected.y_samples)
        assert_arrays_equal(actual.kept, expected.kept)
        assert actual.kept.dtype == np.bool_
    with fits.open(path) as hdul:
        header = hdul[0].header
        assert header["MPSHIFT"] == pytest.approx(0.37)
        assert header["MPMISS"] == "3" and header["MPEXTRA"] == "44"
        assert header["MPNY"] == 400 and header["MPNX"] == N_PIX
        assert hdul["TRACES"].data["COEFFS"].shape == (3, 4)
        assert list(hdul["TRACES"].data["DEGREE"]) == [1, 2, 3]
        assert len(hdul["SAMPLES"].data) == 21

    bare = TraceSet(traces[:1], shape=(400, N_PIX))
    write_traceset(tmp_path / "bare.fits", bare, prov)
    back, _ = read_traceset(tmp_path / "bare.fits")
    assert back.layout_shift is None and back.missing == () and back.extra == ()
    header = fits.getheader(tmp_path / "bare.fits")
    assert header["MPMISS"] == "" and header["MPEXTRA"] == ""

    quartic = TraceSet([make_trace(0, [1.0, 1.0, 1.0, 1.0, 1.0])], shape=(400, N_PIX))
    with pytest.raises(ValueError, match="degree"):
        write_traceset(tmp_path / "quartic.fits", quartic, prov)


# --------------------------------------------------------------------------------- spectra


def test_meta_targets_round_trip_and_legacy_target_card(tmp_path):
    es = ExposureSpectra(
        meta=META,
        spectra=[make_spectrum(0, np.random.default_rng(3), with_versions=False)],
        wavelength_status="transferred",
        flux_status="uncalibrated",
        telluric_status="uncorrected",
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
        drift_arcs=(CHECKSUMS[1], CHECKSUMS[2]),
    )
    path = tmp_path / product_name(META, "spec")
    write_spectra(path, es, make_provenance())
    back = read_spectra(path)[0].meta
    assert back.targets == META.targets and isinstance(back.targets[0], tuple)
    # a product written before targets per fiber carries a single target string
    with fits.open(path, mode="update") as hdul:
        header = hdul[0].header
        n = next(k for k in range(1, 100) if header[f"MPMK{k:02d}"] == "targets")
        header[f"MPMK{n:02d}"] = "target"
        header[f"MPMV{n:02d}"] = json.dumps("target_a")
    assert read_spectra(path)[0].meta.targets == (("*", "target_a"),)


def test_spectra_roundtrip_with_versions_and_drift(tmp_path):
    rng = np.random.default_rng(2)
    spectra = [make_spectrum(k, rng, with_versions=True) for k in range(3)]
    es = ExposureSpectra(
        meta=META,
        spectra=spectra,
        wavelength_status="transferred",
        flux_status="uncalibrated",
        telluric_status="uncorrected",
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
        drift_arcs=(CHECKSUMS[1], CHECKSUMS[2]),
        provenance={"idempotence_key": "abc", "n_arcs": 2},
    )
    drift = make_drift()
    prov = make_provenance()
    path = tmp_path / product_name(META, "spec")
    write_spectra(path, es, prov, drift, weights={0: 0.6, 44: 0.4})
    back, prov_back, drift_back = read_spectra(path)
    assert prov_back == prov
    assert back.meta == META
    assert back.wavelength_status == "transferred"
    assert back.flux_status == "uncalibrated" and back.telluric_status == "uncorrected"
    assert back.reference_id == es.reference_id and back.reference_hash == es.reference_hash
    assert back.drift_arcs == es.drift_arcs
    assert back.provenance == es.provenance
    assert len(back.spectra) == 3
    for actual, expected in zip(back.spectra, spectra, strict=True):
        assert_spectrum_equal(actual, expected)
        assert set(actual.versions) == {"skysub"}
    assert drift_back is not None
    assert_drift_equal(drift_back, drift)
    with fits.open(path) as hdul:
        header = hdul["SPECTRA"].header
        assert header["MPHASSKY"] is True
        assert header["MPHASFLX"] is False and header["MPHASTEL"] is False
        table = hdul["SPECTRA"].data
        assert np.all(np.isnan(table["FLUX_FLUXCAL"])) and np.all(table["MASK_FLUXCAL"])
        assert table["FLUX_SKYSUB"].shape == (3, N_SPEC)
        assert "DRIFT" in hdul
        drift_header = hdul["DRIFT"].header
        assert drift_header["MPDSCAT"] == pytest.approx(0.021)
        assert drift_header["MPDUSE"] is True
        assert drift_header["MPDTIME"] == drift.arc_mid_time
        assert drift_header["MPDNPIX"] == N_PIX
        assert drift_header["MPDNSH"] == 2 and drift_header["MPDSH01"] == pytest.approx(1e-4)
        assert drift_header["MPDNST"] == 1 and drift_header["MPDST01"] == pytest.approx(2e-6)
        assert list(hdul["DRIFT"].data["WEIGHT"]) == [0.6, 0.4]
        primary = hdul[0].header
        assert primary["MPNMK"] == 20
        meta_keys = {primary[f"MPMK{n:02d}"] for n in range(1, 21)}
        assert {"instrument", "camera", "exposure_type", "binning", "unresolved"} <= meta_keys
        assert "amplifier" in meta_keys
        assert "timezone" in meta_keys
        assert json.loads(primary["MPMV08"]) == [1, 1]


def test_drift_n_usable_round_trips_and_mpduse_says_why(tmp_path):
    rng = np.random.default_rng(3)
    es = ExposureSpectra(
        meta=META,
        spectra=[make_spectrum(0, rng, with_versions=False)],
        wavelength_status="transferred",
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
        drift_arcs=(CHECKSUMS[1], CHECKSUMS[2]),
    )
    drift = replace(make_drift(), usable=False, n_usable=2)
    path = tmp_path / product_name(META, "spec")
    write_spectra(path, es, make_provenance(), drift)
    _, _, back = read_spectra(path)
    assert back is not None and back.n_usable == 2 and back.usable is False
    with fits.open(path, mode="update") as hdul:
        header = hdul["DRIFT"].header
        assert header["MPDNUSE"] == 2 and header["MPDUSE"] is False
        assert "2 usable arc fibers" in header.comments["MPDUSE"]
        del header["MPDNUSE"]
    # a product written before the card existed derives the count from the fiber table
    _, _, back = read_spectra(path)
    assert back is not None and back.n_usable == 2


def make_drift() -> DriftModel:
    return DriftModel(
        arc_checksums=(CHECKSUMS[1], CHECKSUMS[2]),
        arc_mid_time="2026-10-01T20:18:00",
        per_fiber=pd.DataFrame(
            {
                "fiber": [0, 44],
                "y_center": [10.5, 390.2],
                "shift_px": [0.31, 0.35],
                "stretch": [0.0, 1e-5],
                "n_lines": [25, 24],
                "rms_px": [0.02, 0.03],
                "usable": [True, True],
            }
        ),
        shift_coeffs=np.array([1e-4, 0.3]),
        stretch_coeffs=np.array([2e-6]),
        scatter_px=0.021,
        usable=True,
        n_pix=N_PIX,
    )


def assert_drift_equal(actual: DriftModel, expected: DriftModel) -> None:
    assert actual.arc_checksums == expected.arc_checksums
    assert actual.arc_mid_time == expected.arc_mid_time
    pd.testing.assert_frame_equal(actual.per_fiber, expected.per_fiber, check_dtype=False)
    assert_arrays_equal(actual.shift_coeffs, expected.shift_coeffs)
    assert_arrays_equal(actual.stretch_coeffs, expected.stretch_coeffs)
    assert actual.scatter_px == pytest.approx(expected.scatter_px)
    assert actual.usable is expected.usable
    assert actual.n_usable == expected.n_usable
    assert actual.n_pix == expected.n_pix
    assert actual.shift_at(100.0) == pytest.approx(expected.shift_at(100.0))


def test_spectra_provenance_disagreement_raises(tmp_path):
    rng = np.random.default_rng(7)
    es = ExposureSpectra(
        meta=META,
        spectra=[make_spectrum(0, rng, with_versions=False)],
        wavelength_status="reference",
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
    )
    prov = Provenance("0.0.1", "abcdef012345", [CHECKSUMS[0]], [], wavelength_status="nominal")
    with pytest.raises(ValueError, match="wavelength_status"):
        write_spectra(tmp_path / "spec.fits", es, prov)
    prov = replace(prov, wavelength_status="reference", reference_hash="deadbeef0000")
    with pytest.raises(ValueError, match="reference_hash"):
        write_spectra(tmp_path / "spec.fits", es, prov)
    prov = replace(prov, reference_hash=None, drift_arcs=[CHECKSUMS[2]])
    with pytest.raises(ValueError, match="drift_arcs"):
        write_spectra(tmp_path / "spec.fits", es, prov)


def test_spectra_bare_provenance_filled_from_exposure(tmp_path):
    rng = np.random.default_rng(8)
    es = ExposureSpectra(
        meta=META,
        spectra=[make_spectrum(0, rng, with_versions=False)],
        wavelength_status="transferred",
        reference_id="highspec_Mg_survey45_20261001T180000",
        reference_hash="0011223344aa",
        drift_arcs=(CHECKSUMS[1],),
    )
    bare = Provenance("0.0.1", "abcdef012345", [CHECKSUMS[0]], ["bias_a"])
    path = tmp_path / "spec.fits"
    write_spectra(path, es, bare)
    back, prov_back, _ = read_spectra(path)
    assert prov_back.wavelength_status == "transferred"
    assert prov_back.reference_id == es.reference_id
    assert prov_back.reference_hash == es.reference_hash
    assert prov_back.drift_arcs == [CHECKSUMS[1]]
    assert prov_back.inputs == [CHECKSUMS[0]] and prov_back.calibrations == ["bias_a"]
    assert back.wavelength_status == "transferred" and back.drift_arcs == (CHECKSUMS[1],)


def test_spectra_roundtrip_without_versions(tmp_path):
    rng = np.random.default_rng(3)
    spectra = [make_spectrum(k, rng, with_versions=False) for k in range(2)]
    es = ExposureSpectra(meta=META, spectra=spectra)
    prov = Provenance("0.0.1", "abcdef012345", [CHECKSUMS[0]], [])
    path = tmp_path / "spec.fits"
    write_spectra(path, es, prov)
    back, prov_back, drift_back = read_spectra(path)
    assert drift_back is None
    assert back.wavelength_status == "none"
    assert back.reference_id is None and back.reference_hash is None
    assert back.drift_arcs == () and back.provenance is None
    assert prov_back == prov
    for actual, expected in zip(back.spectra, spectra, strict=True):
        assert_spectrum_equal(actual, expected)
        assert actual.versions == {}
        assert actual.wavelength is None and actual.sky is None and actual.continuum is None
    with fits.open(path) as hdul:
        assert "DRIFT" not in hdul
        header = hdul["SPECTRA"].header
        for key in ("MPHASSKY", "MPHASFLX", "MPHASTEL", "MPGRDWAV", "MPGRDSKY", "MPGRDCON"):
            assert header[key] is False, key


def test_spectra_partial_versions_per_fiber(tmp_path):
    rng = np.random.default_rng(4)
    spectra = [make_spectrum(k, rng, with_versions=k == 1) for k in range(3)]
    es = ExposureSpectra(meta=META, spectra=spectra, wavelength_status="reference")
    prov = Provenance("0.0.1", "abcdef012345", [CHECKSUMS[0]], [])
    path = tmp_path / "spec.fits"
    write_spectra(path, es, prov)
    back, _, _ = read_spectra(path)
    assert fits.getheader(path, "SPECTRA")["MPHASSKY"] is True
    for actual, expected in zip(back.spectra, spectra, strict=True):
        assert_spectrum_equal(actual, expected)
    assert back.spectra[0].versions == {} and set(back.spectra[1].versions) == {"skysub"}


# ------------------------------------------------------------------------------- reference


def make_lines(x: list[float], wavelength: list[float]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "x": np.asarray(x, dtype=float),
            "wavelength": np.asarray(wavelength, dtype=float),
            "species": ["Th I", "Ar I", "Th II"][: len(x)],
            "dx_px": np.linspace(-0.02, 0.02, len(x)),
        }
    )[list(LINE_COLUMNS)]


def make_solution(fiber: int, model: GratingModel, lines: pd.DataFrame) -> WavelengthSolution:
    return WavelengthSolution(
        fiber=fiber,
        kind=model.kind,
        params=model.params,
        n_pix=model.n_pix,
        lines=lines,
        rms_px=0.013 + fiber * 0.001,
        n_lines=len(lines),
        anchor_px=(float(lines["x"].min()), float(lines["x"].max())),
        anchor_nm=(float(lines["wavelength"].min()), float(lines["wavelength"].max())),
        trace_y_center=10.0 + 8.5 * fiber,
        geometry={
            "d_nm": model.geom.d_nm,
            "alpha_rad": model.geom.alpha_rad,
            "order": model.geom.order,
            "B_nom": model.geom.B_nom,
        },
    )


@pytest.fixture
def reference(highspec_effective):
    geom = grating_geometry(highspec_effective.geometry, highspec_effective.detector.pixel_size_um)
    nominal = model_from_config("grating", geom, N_PIX)
    assert isinstance(nominal, GratingModel)
    with_poly = replace(nominal, x_c=nominal.x_c + 3.25, poly=np.array([0.4, 0.0, 0.0, 0.0]))
    x = [200.5, 800.25, 1400.75]
    solutions = {
        0: make_solution(0, nominal, make_lines(x, nominal.wavelength(np.array(x)))),
        2: make_solution(2, with_poly, make_lines(x[:2], with_poly.wavelength(np.array(x[:2])))),
    }
    ref = build_reference(
        "highspec", "Mg", "survey45", solutions, "2026-10-01T18:00:00", CHECKSUMS[:2]
    )
    return ref, geom


def test_reference_roundtrip(tmp_path, reference):
    ref, geom = reference
    prov = make_provenance()
    path = tmp_path / reference_name(ref)
    assert path.name == "highspec_Mg_survey45_20261001T180000_wref.fits"
    write_reference(path, ref, prov)
    back, prov_back = read_reference(path, lambda instrument, configuration: geom)
    assert prov_back == prov
    for name in ("instrument", "configuration", "fiber_preset", "created", "valid_from"):
        assert getattr(back, name) == getattr(ref, name), name
    assert back.identifier == ref.identifier
    assert back.source_checksums == ref.source_checksums
    assert back.hash == ref.hash
    assert solutions_hash(back.solutions) == ref.hash
    assert set(back.solutions) == {0, 2}
    grid = np.arange(N_PIX, dtype=float)
    for fiber, expected in ref.solutions.items():
        actual = back.solutions[fiber]
        assert actual.fiber == fiber and actual.kind == "grating"
        assert_arrays_equal(actual.params, expected.params)
        assert len(actual.params) == len(expected.params)
        assert actual.n_pix == expected.n_pix and actual.n_lines == expected.n_lines
        assert actual.rms_px == pytest.approx(expected.rms_px)
        assert actual.anchor_px == pytest.approx(expected.anchor_px)
        assert actual.anchor_nm == pytest.approx(expected.anchor_nm)
        assert actual.trace_y_center == pytest.approx(expected.trace_y_center)
        assert actual.geometry == expected.geometry
        pd.testing.assert_frame_equal(actual.lines, expected.lines, check_dtype=False)
        model = model_from_solution(actual, geom)
        np.testing.assert_allclose(
            model.wavelength(grid), expected_model(expected, geom), rtol=0, atol=1e-12
        )
    with fits.open(path) as hdul:
        header = hdul[0].header
        assert header["MPREFID"] == ref.identifier and header["MPREFH"] == ref.hash
        assert header["MPCREAT"] == ref.created and header["MPVALID"] == ref.valid_from
        assert header["MPINST"] == "highspec" and header["MPCONF"] == "Mg"
        assert header["MPPRESET"] == "survey45"
        assert hdul["SOLUTIONS"].data["PARAMS"].shape == (2, 12)
        assert list(hdul["SOLUTIONS"].data["NPARAMS"]) == [3, 7]
        assert len(hdul["LINES"].data) == 5

    other = replace(geom, order=geom.order + 1)
    with pytest.raises(ValueError, match="geometry"):
        read_reference(path, lambda instrument, configuration: other)


def expected_model(solution: WavelengthSolution, geom) -> np.ndarray:
    return model_from_solution(solution, geom).wavelength(np.arange(solution.n_pix, dtype=float))


# -------------------------------------------------------------------------------- combined


def test_combined_roundtrip(tmp_path):
    rng = np.random.default_rng(5)
    wavelength = 515.0 + 0.004 * np.arange(N_SPEC)
    raw = SpectrumVersion(
        rng.normal(1.0, 0.1, N_SPEC), np.full(N_SPEC, 0.1), np.zeros(N_SPEC, bool)
    )
    skysub = SpectrumVersion(raw.flux - 0.1, raw.error, np.arange(N_SPEC) == 3)
    contributors = pd.DataFrame(
        {
            "fiber": [1, 4],
            "telescope": ["mast01", "mast02"],
            "source": ["sim_science_01.fits", "sim_science_01.fits"],
            "weight": [0.6, 0.4],
        }
    )
    excluded = pd.DataFrame(
        {
            "fiber": [7],
            "telescope": ["mast03"],
            "source": ["sim_science_01.fits"],
            "reason": ["low snr"],
        }
    )
    cs = CombinedSpectrum(
        target="HD 1234",
        wavelength=wavelength,
        versions={"raw": raw, "skysub": skysub},
        contributors=contributors,
        excluded=excluded,
        grouping="exposure",
        barycentric={"bjd_tdb": 2461000.5, "rv_kms": 12.25},
        version_missing={"skysub": [4]},
    )
    prov = make_provenance()
    path = tmp_path / combined_name(META, cs.target)
    write_combined(path, cs, META, prov)
    back, meta_back, prov_back = read_combined(path)
    assert prov_back == prov and meta_back == META
    assert back.target == "HD 1234" and back.grouping == "exposure"
    assert back.barycentric == cs.barycentric
    assert back.version_missing == {"skysub": [4]}
    assert_arrays_equal(back.wavelength, wavelength)
    assert set(back.versions) == {"raw", "skysub"}
    assert_version_equal(back.versions["raw"], raw)
    assert_version_equal(back.versions["skysub"], skysub)
    pd.testing.assert_frame_equal(back.contributors, contributors, check_dtype=False)
    pd.testing.assert_frame_equal(back.excluded, excluded, check_dtype=False)
    with fits.open(path) as hdul:
        columns = hdul["COMBINED"].columns.names
        assert {"WAVELENGTH", "FLUX_RAW", "ERROR_RAW", "MASK_RAW", "FLUX_SKYSUB"} <= set(columns)
        assert "FLUX_FLUXCAL" not in columns
        assert len(hdul["COMBINED"].data) == N_SPEC
        header = hdul[0].header
        assert header["MPTARGET"] == "HD 1234" and header["MPGROUP"] == "exposure"
        assert header["MPBCK01"] == "bjd_tdb" and json.loads(header["MPBCV01"]) == 2461000.5
        assert len(hdul["CONTRIB"].data) == 2 and len(hdul["EXCLUDED"].data) == 1
        assert json.loads(header["MPVMISS"]) == {"skysub": [4]}

    empty = replace(cs, excluded=excluded.iloc[:0], versions={"raw": raw}, version_missing={})
    write_combined(tmp_path / "empty.fits", empty, META, prov)
    back, _, _ = read_combined(tmp_path / "empty.fits")
    assert len(back.excluded) == 0 and list(back.excluded.columns) == list(excluded.columns)
    assert set(back.versions) == {"raw"} and back.version_missing == {}
    # a product written before the card existed reads as nothing missing
    with fits.open(path, mode="update") as hdul:
        del hdul[0].header["MPVMISS"]
    assert read_combined(path)[0].version_missing == {}


# -------------------------------------------------------------------------------------- qc


def test_qc_json_roundtrip(tmp_path):
    report = QCReport(
        stage="extract",
        values={"median_snr": np.float64(12.5), "n_masked": np.int64(3), "ok": np.bool_(True)},
    )
    path = tmp_path / "qc.json"
    write_qc(path, report)
    text = path.read_text(encoding="utf-8")
    assert json.loads(text) == {
        "stage": "extract",
        "values": {"median_snr": 12.5, "n_masked": 3, "ok": True},
    }
    assert text.index('"stage"') < text.index('"values"')
    back = read_qc(path)
    assert back.stage == "extract" and back.values == {
        "median_snr": 12.5,
        "n_masked": 3,
        "ok": True,
    }


# ------------------------------------------------------------------------ no raw keywords


def test_no_raw_keywords_in_products(tmp_path, reference):
    rng = np.random.default_rng(6)
    prov = make_provenance()
    ref, geom = reference
    master = MasterCalib("bias", rng.normal(size=(4, 6)), {"camera": "main"}, 2, variance=None)
    traceset = TraceSet([make_trace(0, [0.002, 10.0], role="arc")], shape=(400, N_PIX))
    es = ExposureSpectra(
        META,
        [make_spectrum(0, rng, with_versions=True)],
        wavelength_status=prov.wavelength_status,
        reference_id=prov.reference_id,
        reference_hash=prov.reference_hash,
        drift_arcs=tuple(prov.drift_arcs),
    )
    cs = CombinedSpectrum(
        "HD 1234",
        np.arange(N_SPEC, dtype=float),
        {"raw": SpectrumVersion(np.ones(N_SPEC), np.ones(N_SPEC), np.zeros(N_SPEC, bool))},
        pd.DataFrame({"fiber": [1], "telescope": ["mast01"], "source": ["s"], "weight": [1.0]}),
        pd.DataFrame({"fiber": [], "telescope": [], "source": [], "reason": []}),
        "exposure",
        {},
    )
    paths = {
        "master": tmp_path / "master.fits",
        "traces": tmp_path / "traces.fits",
        "spectra": tmp_path / "spectra.fits",
        "reference": tmp_path / "reference.fits",
        "combined": tmp_path / "combined.fits",
    }
    write_master(paths["master"], master, prov)
    write_traceset(paths["traces"], traceset, prov)
    write_spectra(paths["spectra"], es, prov)
    write_reference(paths["reference"], ref, prov)
    write_combined(paths["combined"], cs, META, prov)
    read_master(paths["master"])
    read_traceset(paths["traces"])
    read_spectra(paths["spectra"])
    read_reference(paths["reference"], lambda instrument, configuration: geom)
    read_combined(paths["combined"])
    for kind, path in paths.items():
        for header in all_headers(path):
            for raw in RAW_KEYWORDS:
                assert raw not in header, f"{kind}: raw keyword {raw} leaked"
            for key in header:
                assert len(key) <= 8 and key == key.upper(), f"{kind}: non-standard card {key}"
