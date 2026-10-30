"""End-to-end quick-look on the simulated night: nominal wavelengths without a reference,
transferred wavelengths against the simulator truth once ``identify`` has made one, the
reference-only night, idempotence with the re-reduce rule, the error policy, and the
``quicklook`` / ``solve`` / ``show`` commands.

The runs are slow (line finding on every arc), so one module-scoped output folder carries
the nominal run, the reference made by ``run_identify`` and the transferred run, and the
tests share them read-only.
"""

import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.cli import main
from mastpipe.config.loader import build_effective
from mastpipe.instruments.registry import load_instrument
from mastpipe.io.ingest import file_checksum
from mastpipe.io.products import (
    combined_name,
    read_combined,
    read_qc,
    read_reference,
    read_spectra,
)
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.identify import run_identify
from mastpipe.pipeline.runner import RunOptions, RunSummary, run_quicklook
from mastpipe.sim.simulator import simulate_night, undrift_x
from mastpipe.steps.atlas import catalogue_arrays
from mastpipe.steps.combine import group_by_target
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.wavemodel import GratingModel

pytestmark = pytest.mark.slow

ANCHOR_LOG_INTENSITY = 2.0
MAX_MAP_ERROR_DISP = 0.3
SURVEY_OVERRIDES = {
    "trace": {"trace_halfwidth": 4},
    "extract": {"halfwidth": 4},
    "lines": {"halfwidth": 4},
}
"""The 45-trace survey layout has an 8.5 px pitch; the HighSpec defaults span the
neighbours, so the survey runs narrow the three halfwidths (as the identify tests do)."""
SURVEY_SET = ["trace.trace_halfwidth=4", "extract.halfwidth=4", "lines.halfwidth=4"]
"""The same overrides as ``--set`` arguments, so the command line hashes the same config."""
SHORT_AGE = {"drift": {"max_age_hours": 0.3}}
"""An age limit of 18 minutes: a copied arc pairs with the science frames next to it only."""
N_SCIENCE = 4


# ------------------------------------------------------------------------------- helpers


def options(*extra_overrides: dict, **kwargs) -> RunOptions:
    return RunOptions(
        instrument="highspec", overrides=[SURVEY_OVERRIDES, *extra_overrides], **kwargs
    )


def truth_anchors(truth, cat_wl) -> list[tuple[float, float]]:
    usable = np.isin(truth.lines_nm, cat_wl)
    bright = usable & truth.lines_on_chip & (truth.line_logI > ANCHOR_LOG_INTENSITY)
    first, last = np.flatnonzero(bright)[[0, -1]]
    return [
        (float(truth.lines_nm[i]), round(float(truth.line_pixels[0][i]), 1)) for i in (first, last)
    ]


def science_files(truth) -> list[str]:
    return truth.frames.loc[truth.frames["exposure_type"] == "science", "file"].tolist()


def spec_products(out_dir: Path) -> list[Path]:
    return sorted((out_dir / "spec").glob("*_spec.fits"))


def expected_arc_before(truth, source: str) -> str:
    """The file of the latest calibration-fiber arc at or before the frame's mid-time."""
    frames = truth.frames
    t = pd.Timestamp(frames.loc[frames["file"] == source, "mid_time"].iloc[0])
    arcs = frames.loc[frames["arc_kind"] == "calfiber"]
    before = arcs.loc[pd.to_datetime(arcs["mid_time"]) <= t]
    return str(before.iloc[-1]["file"])


def map_error_in_dispersions(lam: np.ndarray, night, shift: float, stretch: float) -> float:
    x = np.arange(night.n_pix, dtype=float)
    lam_true = night.truth_model.wavelength(undrift_x(x, shift, stretch, night.n_pix))
    return float(np.max(np.abs(lam - lam_true) / night.truth_model.dispersion(x)))


# ------------------------------------------------------------------------------ fixtures


@pytest.fixture(scope="module")
def night(sim_night, mgb518_lines) -> SimpleNamespace:
    night_dir, truth = sim_night
    effective = build_effective(load_instrument("highspec"), "mgb518", "main", [SURVEY_OVERRIDES])
    geom = grating_geometry(effective.geometry, effective.detector.pixel_size_um)
    n_pix = truth.extras["n_pix"]
    cat_wl, _ = catalogue_arrays(mgb518_lines, usable_only=True)
    full = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    return SimpleNamespace(
        dir=night_dir,
        truth=truth,
        effective=effective,
        geom=geom,
        n_pix=n_pix,
        anchors=truth_anchors(truth, cat_wl),
        truth_model=GratingModel(truth.x_c, truth.B, geom.alpha_rad, None, n_pix, geom),
        full_arc=night_dir / str(full["file"]),
        cal_arcs=truth.frames.loc[truth.frames["arc_kind"] == "calfiber", "file"].tolist(),
    )


@pytest.fixture(scope="module")
def out_dir(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("quicklook")


@pytest.fixture(scope="module")
def nominal_run(night, out_dir) -> RunSummary:
    """The night reduced before any reference exists."""
    return run_quicklook(night.dir, out_dir, options())


@pytest.fixture(scope="module")
def identified(nominal_run, night, out_dir) -> SimpleNamespace:
    """The reference made by ``identify`` from truth anchors, into the same output folder."""
    ref, path = run_identify(
        night.full_arc, night.effective, "survey45", out_dir, anchors=night.anchors, anchor_trace=0
    )
    _, prov = read_reference(path, lambda instrument, configuration: night.geom)
    return SimpleNamespace(ref=ref, path=path, prov=prov)


@pytest.fixture(scope="module")
def transferred_run(identified, night, out_dir) -> RunSummary:
    """The night reduced again, now with the reference and its drift arcs."""
    return run_quicklook(night.dir, out_dir, options())


def no_cal_night(tmp_path: Path, night, survey_layout, mgb518_lines, identified) -> Path:
    """A night without calibration-fiber arcs, its library seeded with the reference."""
    folder = tmp_path / "night"
    simulate_night(
        night.effective,
        survey_layout,
        mgb518_lines,
        folder,
        seed=0,
        n_science=N_SCIENCE,
        no_cal_arcs=True,
    )
    out = tmp_path / "out"
    CalibLibrary(out / "calib").add_reference(identified.ref, identified.prov, "main")
    return folder


# --------------------------------------------------------------------- no reference yet


def test_quicklook_without_reference_gives_nominal(nominal_run, night, out_dir):
    summary = nominal_run
    assert summary.failed == []
    assert summary.reduced_again == []
    sources = science_files(night.truth)
    products = spec_products(out_dir)
    assert len(products) == N_SCIENCE
    for path in products:
        es, prov, drift = read_spectra(path)
        assert es.meta.source in sources
        assert es.wavelength_status == prov.wavelength_status == "nominal"
        assert prov.reference_id is None and prov.drift_arcs == []
        assert drift is None
        assert all(spectrum.wavelength is not None for spectrum in es.spectra)
        assert prov.inputs == [es.meta.checksum]
    # the full arc cannot be solved and the drift arcs have nothing to drift against
    assert any(night.full_arc.name in item and "run identify" in item for item in summary.skipped)
    for arc in night.cal_arcs:
        assert any(arc in item and "no reference" in item for item in summary.skipped)
    # masters, traces and the log are in place
    library = CalibLibrary(out_dir / "calib")
    assert {entry["kind"] for entry in library.list()} == {"bias", "dark", "background", "flat"}
    assert library.find("bias", "main", (1, 1), "image") is not None
    assert len(list((out_dir / "calib").glob("*_traces.fits"))) == 1
    assert any(path.name.endswith("_traces.json") for path in summary.qc)
    logs = sorted((out_dir / "logs").glob("quicklook_*.log"))
    assert logs
    assert "run identify" in logs[0].read_text(encoding="utf-8")


# ------------------------------------------------------------ with reference and arcs


def test_quicklook_end_to_end_on_simulated_night(transferred_run, identified, night, out_dir):
    summary = transferred_run
    truth = night.truth
    assert summary.failed == []
    sources = science_files(truth)
    # the nominal products were reduced again, now against the reference and its arcs
    assert sorted(summary.reduced_again) == sorted(sources)
    assert any(
        night.full_arc.name in item and identified.ref.identifier in item
        for item in summary.skipped
    )
    products = spec_products(out_dir)
    assert len(products) == N_SCIENCE
    checksums = {name: file_checksum(night.dir / name) for name in night.cal_arcs}
    science_traces = 0
    for path in products:
        es, prov, drift = read_spectra(path)
        row = truth.frames.loc[truth.frames["file"] == es.meta.source].iloc[0]
        assert es.wavelength_status == prov.wavelength_status == "transferred"
        assert prov.reference_id == identified.ref.identifier
        assert prov.reference_hash == identified.ref.hash
        expected_arc = expected_arc_before(truth, es.meta.source)
        assert prov.drift_arcs == [checksums[expected_arc]]
        assert drift is not None and drift.usable
        assert drift.arc_checksums == (checksums[expected_arc],)
        assert drift.per_fiber["usable"].all()
        shift, stretch = float(row["shift_px"]), float(row["stretch"])
        for spectrum in es.spectra:
            if spectrum.role != "science":
                continue
            science_traces += 1
            assert spectrum.wavelength is not None
            error = map_error_in_dispersions(spectrum.wavelength, night, shift, stretch)
            assert error < MAX_MAP_ERROR_DISP, (es.meta.source, spectrum.fiber, error)
            assert "skysub" in spectrum.versions
            assert spectrum.sky is not None
            assert np.isfinite(
                spectrum.versions["skysub"].flux[~spectrum.versions["skysub"].mask]
            ).all()
        # QC lists the drift shift that was used and the per-fiber numbers
        qc = read_qc(out_dir / "qa" / path.with_suffix(".json").name)
        assert qc.stage == "science"
        assert qc.values["wavelength_status"] == "transferred"
        assert qc.values["drift_arcs"] == [checksums[expected_arc]]
        shifts = [s for s in qc.values["drift_shift_px"] if s is not None]
        assert shifts and np.allclose(shifts, drift.shift_at(0.0), atol=0.2)
        assert len(qc.values["median_snr"]) == len(es.spectra)
        assert len(qc.values["masked_fraction"]) == len(es.spectra)
        assert qc.values["sky"] == []
        # a combined product per target of the exposure
        for target in group_by_target(es.spectra):
            combined = out_dir / "spec" / combined_name(es.meta, target)
            assert combined.is_file(), combined.name
            cs, meta, comb_prov = read_combined(combined)
            assert cs.target == target and meta.source == es.meta.source
            assert {"raw", "skysub"} <= set(cs.versions)
            assert comb_prov.wavelength_status == "transferred"
            assert len(cs.contributors) == sum(
                1 for s in es.spectra if s.role == "science" and s.target == target
            )
            assert cs.excluded.empty
    assert science_traces == N_SCIENCE * len(night.truth.layout.query("role == 'science'"))
    # each drift arc has its QC file with the fitted drift
    drift_values = [read_qc(p).values for p in summary.qc if p.name.endswith("_drift.json")]
    assert len(drift_values) == len(night.cal_arcs)
    for values in drift_values:
        assert values["usable"] is True
        assert values["reference_id"] == identified.ref.identifier
        assert len(values["per_fiber"]) == len(night.truth.layout.query("role == 'arc'"))


def test_quicklook_without_drift_arcs_gives_reference_status(
    tmp_path, night, survey_layout, mgb518_lines, identified
):
    folder = no_cal_night(tmp_path, night, survey_layout, mgb518_lines, identified)
    summary = run_quicklook(folder, tmp_path / "out", options())
    assert summary.failed == []
    products = spec_products(tmp_path / "out")
    assert len(products) == N_SCIENCE
    for path in products:
        es, prov, drift = read_spectra(path)
        assert es.wavelength_status == prov.wavelength_status == "reference"
        assert prov.reference_id == identified.ref.identifier
        assert prov.drift_arcs == [] and drift is None
        assert all(s.wavelength is not None for s in es.spectra)
        qc = read_qc(tmp_path / "out" / "qa" / path.with_suffix(".json").name)
        assert qc.values["wavelength_status"] == "reference"
        assert all(s is None for s in qc.values["drift_shift_px"])


def test_quicklook_idempotent_then_reduced_again(
    tmp_path, night, survey_layout, mgb518_lines, identified
):
    folder = no_cal_night(tmp_path, night, survey_layout, mgb518_lines, identified)
    out = tmp_path / "out"
    first = run_quicklook(folder, out, options(SHORT_AGE))
    assert first.failed == [] and first.reduced_again == []
    sources = sorted(p.name for p in spec_products(out))
    assert len(sources) == N_SCIENCE

    # nothing changed: every science product is skipped
    again = run_quicklook(folder, out, options(SHORT_AGE))
    assert again.failed == [] and again.reduced_again == []
    skipped = [item for item in again.skipped if "already reduced" in item]
    assert len(skipped) == N_SCIENCE
    assert not any(p.name.endswith("_spec.fits") for p in again.products)

    # a calibration-fiber arc arrives: the frames within its age limit are reduced again
    arc = night.cal_arcs[0]
    shutil.copy(night.dir / arc, folder / arc)
    third = run_quicklook(folder, out, options(SHORT_AGE))
    assert third.failed == []
    t_arc = pd.Timestamp(
        night.truth.frames.loc[night.truth.frames["file"] == arc, "mid_time"].iloc[0]
    )
    statuses = {}
    for path in spec_products(out):
        es, prov, _ = read_spectra(path)
        statuses[es.meta.source] = (prov.wavelength_status, pd.Timestamp(es.meta.mid_time))
    near = sorted(s for s, (_, t) in statuses.items() if abs(t - t_arc) <= pd.Timedelta("18min"))
    far = sorted(s for s, (_, t) in statuses.items() if abs(t - t_arc) > pd.Timedelta("18min"))
    assert near and far
    assert sorted(third.reduced_again) == near
    assert all(statuses[s][0] == "transferred" for s in near)
    assert all(statuses[s][0] == "reference" for s in far)
    assert all(any(s in item and "already reduced" in item for item in third.skipped) for s in far)
    # and with --force everything is reduced once more
    forced = run_quicklook(folder, out, options(SHORT_AGE, force=True))
    assert forced.failed == []
    assert not any("already reduced" in item for item in forced.skipped)


def test_quicklook_interpolate_blends_the_two_arcs(
    tmp_path, transferred_run, identified, night, out_dir
):
    """Under ``interpolate`` a science frame between the two drift arcs carries both: the
    product's DRIFT table reads back as one blended model with both checksums and the
    QC lists the two weights. The run shares the module library but writes elsewhere."""
    out = tmp_path / "interp"
    summary = run_quicklook(
        night.dir, out, options(time_policy="interpolate", library=out_dir / "calib")
    )
    assert summary.failed == []
    truth = night.truth
    arcs = truth.frames.loc[truth.frames["arc_kind"] == "calfiber"]
    t_first, t_last = (
        pd.Timestamp(arcs.iloc[0]["mid_time"]),
        pd.Timestamp(arcs.iloc[-1]["mid_time"]),
    )
    checksums = {file_checksum(night.dir / name) for name in night.cal_arcs}
    between = 0
    for path in spec_products(out):
        es, prov, drift = read_spectra(path)
        assert es.wavelength_status == "transferred" and drift is not None
        t = pd.Timestamp(es.meta.mid_time)
        qc = read_qc(out / "qa" / path.with_suffix(".json").name)
        if t_first < t < t_last:
            between += 1
            assert set(prov.drift_arcs) == checksums and len(prov.drift_arcs) == 2
            assert drift.arc_checksums == tuple(sorted(checksums))
            weights = qc.values["drift_weights"]
            assert len(weights) == 2 and sum(weights) == pytest.approx(1.0)
            assert all(0.0 < w < 1.0 for w in weights)
            assert t_first < pd.Timestamp(drift.arc_mid_time) < t_last
            assert len(drift.per_fiber) == len(truth.layout.query("role == 'arc'"))
            row = truth.frames.loc[truth.frames["file"] == es.meta.source].iloc[0]
            for spectrum in es.spectra:
                if spectrum.role == "science":
                    error = map_error_in_dispersions(
                        spectrum.wavelength, night, float(row["shift_px"]), float(row["stretch"])
                    )
                    assert error < MAX_MAP_ERROR_DISP
        else:
            assert len(prov.drift_arcs) == 1 and qc.values["drift_weights"] == [1.0]
    assert between == 2


# ------------------------------------------------------------------------ error policy


def write_frame(path: Path, data: np.ndarray, **cards) -> None:
    hdu = fits.PrimaryHDU(data)
    hdu.header.update(cards)
    hdu.writeto(path, overwrite=True)


def test_unknown_frames_are_skipped_with_reason(tmp_path):
    folder = tmp_path / "night"
    folder.mkdir()
    # no exposure time and no image type: the header map's "unknown" rule
    write_frame(
        folder / "junk.fits",
        np.zeros((4, 4), dtype=np.uint16),
        **{"DATE-OBS": "2026-10-01T20:00:00"},
    )
    # a frame the manifest declares unknown
    write_frame(
        folder / "weird.fits",
        np.zeros((4, 4), dtype=np.uint16),
        IMAGETYP="light",
        EXPOSURE=10.0,
        **{"DATE-OBS": "2026-10-01T20:12:00"},
    )
    # and one science frame, so that the group exists but has nothing to calibrate with
    write_frame(
        folder / "sci.fits",
        np.zeros((4, 4), dtype=np.uint16),
        IMAGETYP="light",
        EXPOSURE=300.0,
        **{"DATE-OBS": "2026-10-01T20:24:00"},
    )
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("source,exposure_type\nweird.fits,unknown\n", encoding="utf-8", newline="")
    summary = run_quicklook(
        folder, tmp_path / "out", RunOptions("highspec", configuration="mgb518", manifest=manifest)
    )
    assert any("junk.fits" in item and "rule:4" in item for item in summary.skipped)
    assert any("weird.fits" in item and "manifest" in item for item in summary.skipped)
    assert not any("sci.fits" in item for item in summary.skipped)
    assert summary.products == []
    # the group has nothing to calibrate with, which is recorded, not raised
    assert len(summary.failed) == 1 and "bias" in summary.failed[0][1]
    with pytest.raises(RuntimeError, match="bias"):
        run_quicklook(
            folder,
            tmp_path / "out",
            RunOptions("highspec", configuration="mgb518", manifest=manifest, on_error="raise"),
        )


def test_quicklook_continues_after_bad_frame(tmp_path, night):
    folder = tmp_path / "night"
    folder.mkdir()
    truth = night.truth
    keep = truth.frames.loc[
        truth.frames["exposure_type"].isin(["bias", "dark", "flat", "background"])
    ]
    for name in keep["file"]:
        shutil.copy(night.dir / name, folder / name)
    good = science_files(truth)[:2]
    for name in good:
        shutil.copy(night.dir / name, folder / name)
    # a science frame whose image is cut down to a stub: the calibration cannot apply
    header = fits.getheader(night.dir / good[0])
    cards = {key: header[key] for key in ("IMAGETYP", "EXPOSURE", "GRATING", "FIBPRE", "OBJECT")}
    cards["DATE-OBS"] = "2026-10-01T23:00:00"
    write_frame(folder / "sim_science_bad.fits", np.full((20, 30), 600, dtype=np.uint16), **cards)
    # and a file that is not FITS at all
    (folder / "sim_science_garbage.fits").write_bytes(b"this is not a FITS file\n" * 200)

    summary = run_quicklook(folder, tmp_path / "out", options())
    failed = dict(summary.failed)
    assert "sim_science_bad.fits" in failed and "shape" in failed["sim_science_bad.fits"]
    assert "sim_science_garbage.fits" in failed
    assert set(failed) == {"sim_science_bad.fits", "sim_science_garbage.fits"}
    produced = {read_spectra(p)[0].meta.source for p in spec_products(tmp_path / "out")}
    assert produced == set(good)
    assert not (
        tmp_path / "out" / "spec" / "highspec_mgb518_main_20261001T230000_spec.fits"
    ).exists()

    with pytest.raises(ValueError, match="shape"):
        run_quicklook(folder, tmp_path / "out2", options(on_error="raise"))


# ------------------------------------------------------------------------------- the CLI


def test_cli_quicklook_and_show(transferred_run, identified, night, out_dir, capsys):
    argv = ["quicklook", str(night.dir), "--out", str(out_dir), "--instrument", "highspec"]
    for item in SURVEY_SET:
        argv += ["--set", item]
    assert main(argv) == 0
    out = capsys.readouterr().out
    assert "skipped" in out and "0 failed" in out
    # the command line hashes the same configuration: everything is already reduced
    assert len(spec_products(out_dir)) == N_SCIENCE
    assert "already reduced" in out

    spec = spec_products(out_dir)[0]
    assert main(["show", str(spec)]) == 0
    out = capsys.readouterr().out
    assert "spectra" in out and "transferred" in out and identified.ref.identifier in out
    assert "science" in out and "sky" in out and "arc" in out
    assert main(["show", str(identified.path)]) == 0
    out = capsys.readouterr().out
    assert "reference" in out and identified.ref.identifier in out and "solved" in out
    traces = next((out_dir / "calib").glob("*_traces.fits"))
    assert main(["show", str(traces)]) == 0
    assert "traces" in capsys.readouterr().out
    bias = next((out_dir / "calib").glob("*_bias.fits"))
    assert main(["show", str(bias)]) == 0
    assert "master" in capsys.readouterr().out
    combined = next((out_dir / "spec").glob("*_comb_*.fits"))
    assert main(["show", str(combined)]) == 0
    assert "combined" in capsys.readouterr().out


def test_cli_solve_uses_the_library(sim_night_lab3, mgb518_lines, tmp_path, capsys):
    night_dir, truth = sim_night_lab3
    full = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    arc = str(night_dir / str(full["file"]))
    out = tmp_path / "out"
    base = ["--out", str(out), "--instrument", "highspec", "--configuration", "mgb518"]
    assert main(["solve", arc, *base]) == 2
    assert "run identify first" in capsys.readouterr().out
    cat_wl, _ = catalogue_arrays(mgb518_lines, usable_only=True)
    pairs = ",".join(f"{wl}:{px}" for wl, px in truth_anchors(truth, cat_wl))
    assert main(["identify", arc, *base, "--anchors", pairs, "--trace", "0"]) == 0
    capsys.readouterr()
    library = CalibLibrary(out / "calib")
    entries = library.list(kind="reference")
    assert len(entries) == 1 and entries[0]["camera"] == "main"
    # an empty library elsewhere has nothing to seed from
    other = tmp_path / "library"
    assert main(["solve", arc, *base, "--library", str(other)]) == 2
    assert "run identify first" in capsys.readouterr().out
    # the seeded solve of the same arc replaces the reference instead of adding a second
    assert main(["solve", arc, *base]) == 0
    text = capsys.readouterr().out
    assert "seeded by" in text and entries[0]["identifier"] in text
    assert len(CalibLibrary(out / "calib").list(kind="reference")) == 1
    assert main(["show", str(out / "calib" / f"{entries[0]['identifier']}_wref.fits")]) == 0
    assert "seeded by" in capsys.readouterr().out
