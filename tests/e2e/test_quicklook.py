"""End-to-end quick-look on the simulated night: nominal wavelengths without a reference,
transferred wavelengths against the simulator truth once ``identify`` has made one, the
reference-only night, idempotence with the re-reduce rule, fibers keyed by layout id when
a fiber is missing or an extra detection appears, the error policy, and the ``quicklook``
/ ``solve`` / ``show`` commands.

The runs are slow (line finding on every arc), so one module-scoped output folder carries
the nominal run, the reference made by ``run_identify`` and the transferred run, and the
tests share them read-only.
"""

import shutil
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.cli import main
from mastpipe.config.loader import build_effective
from mastpipe.instruments import registry
from mastpipe.instruments.registry import layout_hash, load_instrument
from mastpipe.io.ingest import file_checksum
from mastpipe.io.products import (
    combined_name,
    read_combined,
    read_master,
    read_qc,
    read_reference,
    read_spectra,
    read_traceset,
)
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.identify import run_identify
from mastpipe.pipeline.runner import NO_FRAMES, RunOptions, RunSummary, run_quicklook
from mastpipe.sim.simulator import BIAS_LEVEL_ADU, simulate_night, undrift_x
from mastpipe.steps.atlas import catalogue_arrays
from mastpipe.steps.combine import group_by_target
from mastpipe.steps.geometry import grating_geometry
from mastpipe.steps.wavecal import build_reference
from tests.conftest import truth_model

pytestmark = pytest.mark.slow

ANCHOR_LOG_INTENSITY = 2.0
MAX_MAP_ERROR_DISP = 0.3
SHORT_AGE = {"drift": {"max_age_hours": 0.3}}
"""An age limit of 18 minutes: a copied arc pairs with the science frames next to it only."""
N_SCIENCE = 4
DROPPED_FIBER = 16
"""The survey fiber (a science fiber of mast08) whose light the missing-fiber night lacks."""
DARK_ROWS = 4
"""Rows on either side of the dropped trace set to the bias level in the flats."""
EXTRA_ROW = 100.0
"""Cropped row of the spurious ridge added between fibers 0 and 1 of the lab3 night."""


# ------------------------------------------------------------------------------- helpers


def options(*extra_overrides: dict, **kwargs) -> RunOptions:
    """The survey night runs on the instrument's own per-preset parameters; its frames
    carry no layout, so the survey45 layout is named."""
    kwargs.setdefault("layout", "survey45")
    return RunOptions(instrument="highspec", overrides=list(extra_overrides), **kwargs)


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


def map_error_in_dispersions(lam: np.ndarray, night, fiber: int, source: str) -> float:
    """Largest ``|lam - lam_true|`` of fiber ``fiber`` on frame ``source``, in dispersions:
    the truth is the fiber's own grating model, drifted by the frame's shift at its row."""
    truth = night.truth
    x = np.arange(night.n_pix, dtype=float)
    model = truth_model(truth, fiber, night.geom, night.n_pix)
    shift = truth.trace_shifts(source)[fiber]
    stretch = float(truth.frames.loc[truth.frames["file"] == source, "stretch"].iloc[0])
    lam_true = model.wavelength(undrift_x(x, shift, stretch, night.n_pix))
    return float(np.max(np.abs(lam - lam_true) / model.dispersion(x)))


def map_errors(out_dir: Path, night) -> dict[str, dict[int, float]]:
    """Every product's map error per fiber against the truth: source -> fiber -> error."""
    errors = {}
    for path in spec_products(out_dir):
        es, _, _ = read_spectra(path)
        errors[es.meta.source] = {
            s.fiber: map_error_in_dispersions(s.wavelength, night, s.fiber, es.meta.source)
            for s in es.spectra
        }
    return errors


# ------------------------------------------------------------------------------ fixtures


@pytest.fixture(scope="module")
def night(sim_night, mg_lines) -> SimpleNamespace:
    night_dir, truth = sim_night
    effective = build_effective(load_instrument("highspec"), "Mg", "main", layout="survey45")
    geom = grating_geometry(effective.geometry, effective.detector.pixel_size_um)
    n_pix = truth.extras["n_pix"]
    cat_wl, _ = catalogue_arrays(mg_lines, usable_only=True)
    full = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    return SimpleNamespace(
        dir=night_dir,
        truth=truth,
        effective=effective,
        geom=geom,
        n_pix=n_pix,
        anchors=truth_anchors(truth, cat_wl),
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


def no_cal_night(tmp_path: Path, night, survey_layout, mg_lines, identified) -> Path:
    """A night without calibration-fiber arcs, its library seeded with the reference."""
    folder = tmp_path / "night"
    simulate_night(
        night.effective,
        survey_layout,
        mg_lines,
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
        assert es.wavelength_status == prov.wavelength_status == "transferred"
        assert prov.reference_id == identified.ref.identifier
        assert prov.reference_hash == identified.ref.hash
        expected_arc = expected_arc_before(truth, es.meta.source)
        assert prov.drift_arcs == [checksums[expected_arc]]
        assert drift is not None and drift.usable
        assert drift.arc_checksums == (checksums[expected_arc],)
        assert drift.per_fiber["usable"].all()
        for spectrum in es.spectra:
            if spectrum.role != "science":
                continue
            science_traces += 1
            assert spectrum.wavelength is not None
            error = map_error_in_dispersions(
                spectrum.wavelength, night, spectrum.fiber, es.meta.source
            )
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
        assert len(qc.values["targets"]) == 20  # one target per telescope unit
        # a combined product per target of the exposure: four targets, five telescopes each
        groups = group_by_target(es.spectra)
        assert sorted(groups) == sorted(truth.extras["targets"]) and len(groups) == 4
        for target in groups:
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
    tmp_path, night, survey_layout, mg_lines, identified
):
    folder = no_cal_night(tmp_path, night, survey_layout, mg_lines, identified)
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
    tmp_path, night, survey_layout, mg_lines, identified
):
    folder = no_cal_night(tmp_path, night, survey_layout, mg_lines, identified)
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


MASTERS = ("bias", "dark", "background", "flat")


def master_files(out: Path) -> dict[str, int]:
    """The master files of the library with their modification times (ns)."""
    calib = out / "calib"
    names = [entry["path"] for entry in CalibLibrary(calib).list() if entry["kind"] in MASTERS]
    return {name: (calib / name).stat().st_mtime_ns for name in names}


def config_hashes(out: Path) -> set[str]:
    return {read_spectra(path)[1].config_hash for path in spec_products(out)}


def test_quicklook_rereduces_when_layout_reference_master_or_policy_changes(
    tmp_path, night, survey_layout, mg_lines, identified, monkeypatch
):
    """Every part of the idempotence key reduces the night again when it changes: the time
    policy (in the configuration hash), the layout CSV, the reference under the same
    identifier and the masters (rebuilt by ``force`` or under other stacking parameters);
    otherwise masters are reused."""
    folder = no_cal_night(tmp_path, night, survey_layout, mg_lines, identified)
    out = tmp_path / "out"
    first = run_quicklook(folder, out, options())
    assert first.failed == [] and first.reduced_again == []
    sources = sorted(read_spectra(path)[0].meta.source for path in spec_products(out))
    assert len(sources) == N_SCIENCE
    index = (out / "calib" / "index.json").read_text(encoding="utf-8")
    masters = master_files(out)
    assert {name.rsplit("_", 1)[1] for name in masters} >= {"bias.fits", "flat.fits"}
    for path in spec_products(out):
        prov = read_spectra(path)[1]
        assert prov.layout_hash == layout_hash("highspec", "survey45")
        assert prov.masters_hash is not None and len(prov.masters_hash) == 12
        assert prov.reference_hash == identified.ref.hash

    # (a) nothing changed: nothing is reduced again and the masters are reused untouched
    again = run_quicklook(folder, out, options())
    assert again.failed == [] and again.reduced_again == []
    assert sum("already reduced" in item for item in again.skipped) == N_SCIENCE
    assert (out / "calib" / "index.json").read_text(encoding="utf-8") == index
    assert master_files(out) == masters
    assert not any(path.name in masters for path in again.products)

    # (b) another time policy is another configuration: every science frame again
    plain = config_hashes(out)
    nearest = run_quicklook(folder, out, options(time_policy="nearest"))
    assert nearest.failed == [] and sorted(nearest.reduced_again) == sources
    assert len(config_hashes(out)) == 1 and config_hashes(out) != plain
    back = run_quicklook(folder, out, options())
    assert sorted(back.reduced_again) == sources and config_hashes(out) == plain

    # (c) a byte appended to the layout CSV
    copy = tmp_path / "instrument" / "highspec"
    shutil.copytree(registry.instrument_dir("highspec"), copy)
    with (copy / "layouts" / "survey45.csv").open("ab") as handle:
        handle.write(b"\n")
    with monkeypatch.context() as patch:
        patch.setattr(registry, "instrument_dir", lambda name: copy)
        edited = run_quicklook(folder, out, options())
        assert edited.failed == [] and sorted(edited.reduced_again) == sources
        assert all(
            read_spectra(path)[1].layout_hash == layout_hash("highspec", "survey45")
            for path in spec_products(out)
        )
    assert sorted(run_quicklook(folder, out, options()).reduced_again) == sources

    # (d) the reference re-solved under the same identifier, with another hash
    ref = identified.ref
    solutions = {}
    for fiber, solution in ref.solutions.items():
        params = np.array(solution.params, dtype=float)
        params[0] += 0.01
        solutions[fiber] = replace(solution, params=params)
    resolved = build_reference(
        ref.instrument,
        ref.configuration,
        ref.fiber_preset,
        solutions,
        ref.created,
        ref.source_checksums,
    )
    assert resolved.identifier == ref.identifier and resolved.hash != ref.hash
    CalibLibrary(out / "calib").add_reference(resolved, identified.prov, "main")
    rehashed = run_quicklook(folder, out, options())
    assert rehashed.failed == [] and sorted(rehashed.reduced_again) == sources
    assert all(read_spectra(path)[1].reference_hash == resolved.hash for path in spec_products(out))
    index = (out / "calib" / "index.json").read_text(encoding="utf-8")
    assert master_files(out) == masters

    # (e) --force rebuilds the masters (and so reduces everything again)
    forced = run_quicklook(folder, out, options(force=True))
    assert forced.failed == [] and sorted(forced.reduced_again) == sources
    rebuilt = master_files(out)
    assert set(rebuilt) == set(masters)
    assert all(rebuilt[name] != masters[name] for name in masters)
    assert {path.name for path in forced.products} >= set(masters)
    # and the next plain run reuses the rebuilt masters
    index = (out / "calib" / "index.json").read_text(encoding="utf-8")
    last = run_quicklook(folder, out, options())
    assert last.failed == [] and last.reduced_again == []
    assert master_files(out) == rebuilt
    assert (out / "calib" / "index.json").read_text(encoding="utf-8") == index

    # (f) other stacking parameters: the masters are built again (not reused under the
    # old ones) and every science frame is reduced again with the new masters
    hashes = {read_spectra(path)[1].masters_hash for path in spec_products(out)}
    mean = run_quicklook(folder, out, options({"calib": {"stack_method": "mean"}}))
    assert mean.failed == [] and sorted(mean.reduced_again) == sources
    restacked = master_files(out)
    assert set(restacked) == set(masters)
    assert all(restacked[name] != rebuilt[name] for name in masters)
    assert {path.name for path in mean.products} >= set(masters)
    assert all(
        read_master(out / "calib" / name)[0].provenance["method"] == "mean" for name in masters
    )
    new = {read_spectra(path)[1].masters_hash for path in spec_products(out)}
    assert len(new) == 1 and new != hashes


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
            for spectrum in es.spectra:
                if spectrum.role == "science":
                    error = map_error_in_dispersions(
                        spectrum.wavelength, night, spectrum.fiber, es.meta.source
                    )
                    assert error < MAX_MAP_ERROR_DISP
        else:
            assert len(prov.drift_arcs) == 1 and qc.values["drift_weights"] == [1.0]
    assert between == 2


# ------------------------------------------------------------------ fibers keyed by layout


def copy_frames(night_dir: Path, truth, folder: Path, types: list[str], arc_kinds=()) -> None:
    frames = truth.frames
    keep = frames["exposure_type"].isin(types) | frames["arc_kind"].isin(list(arc_kinds))
    for name in frames.loc[keep, "file"]:
        shutil.copy(night_dir / name, folder / name)


def edit_flats(folder: Path, truth, change) -> None:
    """Apply ``change(data, x_min, y_min)`` in place to every raw flat of the copied night."""
    x_min, _, y_min, _ = truth.crop
    for name in truth.frames.loc[truth.frames["exposure_type"] == "flat", "file"]:
        with fits.open(folder / name, mode="update") as hdul:
            change(hdul[0].data, x_min, y_min)


def layout_rows(truth) -> dict[int, tuple[str, str | None]]:
    """Layout trace id -> ``(role, telescope)`` of the simulated night."""
    rows = {}
    for _, row in truth.layout.iterrows():
        telescope = row["telescope"]
        if telescope is None or (isinstance(telescope, float) and np.isnan(telescope)):
            telescope = None
        rows[int(row["trace"])] = (str(row["role"]), telescope)
    return rows


SWAPPED = (1, 3)
"""Two science fibers (mast01, mast02) whose light the swapped night exchanges."""
SWAP_JITTER_PX = 1.0
"""Per-trace dispersion jitter of the swapped night. At the default 0.3 px the maps of
fibers 1 and 3 differ by at most 0.43 px on the chip (seed 0) and the swap shows as a map
error of 0.36 and 0.74 dispersions against 0.26 at most for the other fibers: detected, but
0.06 above ``MAX_MAP_ERROR_DISP``. At 1 px the swapped fibers are off by 1.5 and 1.9
dispersions and the other fibers stay at 0.26 at most, so the test does not hang on a
few hundredths."""


def test_swapped_fibers_are_detected(tmp_path, night, survey_layout, mg_lines):
    """A night whose fibers 1 and 3 carry each other's spectra, its truth unswapped: the
    pipeline keys fibers correctly, so its wavelengths are those of the light each trace
    carries and the truth comparison fails for exactly the two swapped fibers. This is what
    makes the end-to-end comparison sensitive to a fiber-identity bug. The identification
    seeds every fiber but the anchored one by cross-correlation, so the jitter is exercised
    there too."""
    folder = tmp_path / "night"
    truth = simulate_night(
        night.effective,
        survey_layout,
        mg_lines,
        folder,
        seed=0,
        n_bias=2,
        n_dark=1,
        n_flat=1,
        n_cal_arc=1,
        n_science=1,
        dispersion_jitter_px=SWAP_JITTER_PX,
        swap_fibers=SWAPPED,
    )
    swapped = SimpleNamespace(truth=truth, geom=night.geom, n_pix=night.n_pix)
    out = tmp_path / "out"
    full = folder / str(truth.frames.loc[truth.frames["arc_kind"] == "full", "file"].iloc[0])
    cat_wl, _ = catalogue_arrays(mg_lines, usable_only=True)
    anchors = truth_anchors(truth, cat_wl)  # fiber 0 is not swapped
    run_identify(full, night.effective, "survey45", out, anchors=anchors, anchor_trace=0)
    summary = run_quicklook(folder, out, options())
    assert summary.failed == []
    errors = map_errors(out, swapped)
    assert len(errors) == 1
    for per_fiber in errors.values():
        assert len(per_fiber) == len(truth.layout)
        wrong = {fiber for fiber, error in per_fiber.items() if error > MAX_MAP_ERROR_DISP}
        assert wrong == set(SWAPPED), per_fiber


def test_quicklook_keys_fibers_by_layout_id_when_one_is_missing(
    tmp_path, transferred_run, identified, night, out_dir
):
    """A night whose flat lacks one fiber's light, reduced against the reference of the
    complete night: the traces are keyed by the layout's trace id, so every other fiber
    keeps its own solution (wavelengths match the truth) and the missing fiber is simply
    absent from the products, with nothing shifted or skipped."""
    truth = night.truth
    folder = tmp_path / "night"
    folder.mkdir()
    copy_frames(
        night.dir, truth, folder, ["bias", "dark", "flat", "background", "science"], ["calfiber"]
    )
    coeffs = truth.trace_coeffs[DROPPED_FIBER]
    x = np.arange(int(truth.extras["n_pix"]))

    def darken(data, x_min, y_min):
        rows = np.rint(np.polyval(coeffs, x)).astype(int) + y_min
        for offset in range(-DARK_ROWS, DARK_ROWS + 1):
            data[rows + offset, x + x_min] = BIAS_LEVEL_ADU

    edit_flats(folder, truth, darken)
    out = tmp_path / "out"
    summary = run_quicklook(folder, out, options(library=out_dir / "calib"))
    assert summary.failed == []

    trace_qc = read_qc(next(p for p in summary.qc if p.name.endswith("_traces.json")))
    assert trace_qc.values["missing"] == [DROPPED_FIBER]
    assert trace_qc.values["n_traces"] == len(truth.layout) - 1
    assert DROPPED_FIBER not in trace_qc.values["index"]
    assert trace_qc.values["extra"] == [] and trace_qc.values["extra_y_center"] == []
    expected = set(truth.layout["trace"].astype(int)) - {DROPPED_FIBER}
    rows = layout_rows(truth)
    products = spec_products(out)
    assert len(products) == N_SCIENCE
    for path in products:
        es, prov, drift = read_spectra(path)
        assert es.wavelength_status == "transferred" and drift is not None
        assert prov.reference_id == identified.ref.identifier and prov.layout == "survey45"
        assert {s.fiber for s in es.spectra} == expected
        for spectrum in es.spectra:
            assert (spectrum.role, spectrum.telescope) == rows[spectrum.fiber]
            assert spectrum.wavelength is not None
            error = map_error_in_dispersions(
                spectrum.wavelength, night, spectrum.fiber, es.meta.source
            )
            assert error < MAX_MAP_ERROR_DISP, (es.meta.source, spectrum.fiber, error)
        qc = read_qc(out / "qa" / path.with_suffix(".json").name)
        assert qc.values["wavelength_skipped"] == []
        assert DROPPED_FIBER not in qc.values["fibers"]
        # the dropped science fiber's telescope keeps its sky fiber, which has nothing to serve
        assert qc.values["sky"] == []


def test_quicklook_drops_detections_outside_the_layout(tmp_path, sim_night_lab3):
    """A bright ridge between two fibers of the three-fiber night is detected on the flat
    but matches no layout row: it is dropped from the traces and the products, logged, and
    listed in the trace QC with its position.

    The layout has arc fibers only, so the simulated "science" frames carry no fiber light:
    read as the simulator writes them (``IMAGETYP = 'light'``) their content has no signal
    and they are skipped; declared ``science`` by a manifest they are reduced, and their
    products carry the layout's fibers."""
    night_dir, truth = sim_night_lab3
    folder = tmp_path / "night"
    folder.mkdir()
    copy_frames(night_dir, truth, folder, ["bias", "dark", "flat", "background", "science"])
    sigma = truth.fwhm_px / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    def ridge(data, x_min, y_min):
        rows = np.arange(-5, 6)
        profile = 15000.0 * np.exp(-0.5 * (rows / sigma) ** 2)
        centre = int(round(EXTRA_ROW)) + y_min
        data[centre - 5 : centre + 6, x_min : data.shape[1] - x_min] += profile[:, None].astype(
            data.dtype
        )

    edit_flats(folder, truth, ridge)
    sources = science_files(truth)
    plain = run_quicklook(folder, tmp_path / "out_light", RunOptions("highspec"))
    assert plain.failed == [] and spec_products(tmp_path / "out_light") == []
    for source in sources:
        assert any(
            item.startswith(source) and "content:no_signal" in item for item in plain.skipped
        ), (source, plain.skipped)

    manifest = tmp_path / "manifest.csv"
    rows = ["source,exposure_type", *(f"{source},science" for source in sources)]
    manifest.write_text("\n".join(rows) + "\n", encoding="utf-8", newline="")
    out = tmp_path / "out"
    summary = run_quicklook(folder, out, RunOptions("highspec", manifest=manifest))
    assert summary.failed == []
    trace_qc = read_qc(next(p for p in summary.qc if p.name.endswith("_traces.json")))
    assert trace_qc.values["layout"] == "lab3"
    assert trace_qc.values["n_detected"] == 4 and trace_qc.values["n_traces"] == 3
    assert trace_qc.values["index"] == [0, 1, 2] and trace_qc.values["missing"] == []
    assert len(trace_qc.values["extra"]) == 1
    assert trace_qc.values["extra_y_center"] == [pytest.approx(EXTRA_ROW, abs=1.0)]
    assert len(trace_qc.values["extra_snr"]) == 1
    products = spec_products(out)
    assert len(products) == len(sources) > 0
    for path in products:
        es, prov, _ = read_spectra(path)
        assert [s.fiber for s in es.spectra] == [0, 1, 2] and prov.layout == "lab3"
    log = next((out / "logs").glob("quicklook_*.log")).read_text(encoding="utf-8")
    assert "match no row of layout lab3" in log and f"y = {EXTRA_ROW:.1f}" in log


def test_reference_on_another_pixel_grid_is_refused(
    tmp_path, night, survey_layout, mg_lines, identified
):
    """A reference solved on another detector crop names itself and the frame width in the
    recorded failure instead of breaking inside the product writer."""
    folder = no_cal_night(tmp_path, night, survey_layout, mg_lines, identified)
    ref = identified.ref
    narrow = build_reference(
        ref.instrument,
        ref.configuration,
        ref.fiber_preset,
        {fiber: replace(s, n_pix=1000) for fiber, s in ref.solutions.items()},
        "2026-10-01T17:00:00",
        ref.source_checksums,
    )
    CalibLibrary(tmp_path / "out" / "calib").add_reference(narrow, identified.prov, "main")
    summary = run_quicklook(folder, tmp_path / "out", options(reference_id=narrow.identifier))
    assert len(summary.failed) == N_SCIENCE
    for source, reason in summary.failed:
        assert source.startswith("sim_science_")
        assert "ValueError" in reason and narrow.identifier in reason
        assert "1000" in reason and str(night.n_pix) in reason
    assert spec_products(tmp_path / "out") == []


# ------------------------------------------------------------------------ error policy


def write_frame(path: Path, data: np.ndarray, **cards) -> None:
    hdu = fits.PrimaryHDU(data)
    hdu.header.update(cards)
    hdu.writeto(path, overwrite=True)


def test_missing_or_empty_night_folder_is_not_a_success(tmp_path, capsys):
    nowhere = tmp_path / "nowhere"
    with pytest.raises(FileNotFoundError, match="nowhere"):
        run_quicklook(nowhere, tmp_path / "out", RunOptions("highspec"))
    assert not (tmp_path / "out").exists()  # nothing is written for a mistyped path
    argv = ["quicklook", str(nowhere), "--out", str(tmp_path / "out"), "--instrument", "highspec"]
    assert main(argv) == 2
    assert "nowhere" in capsys.readouterr().out

    empty = tmp_path / "empty"
    empty.mkdir()
    summary = run_quicklook(empty, tmp_path / "out", RunOptions("highspec"))
    assert summary.products == [] and summary.skipped == []
    assert summary.failed == [("night empty", NO_FRAMES)]
    argv = ["quicklook", str(empty), "--out", str(tmp_path / "out"), "--instrument", "highspec"]
    assert main(argv) == 2
    assert NO_FRAMES in capsys.readouterr().out
    # a folder whose only file cannot be read is as empty
    (empty / "junk.fits").write_bytes(b"this is not a FITS file\n" * 200)
    summary = run_quicklook(empty, tmp_path / "out", RunOptions("highspec"))
    assert [source for source, _ in summary.failed] == ["junk.fits", "night empty"]


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
        folder, tmp_path / "out", RunOptions("highspec", configuration="Mg", manifest=manifest)
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
            RunOptions("highspec", configuration="Mg", manifest=manifest, on_error="raise"),
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
    cards = {key: header[key] for key in ("IMAGETYP", "EXPOSURE", "GRATING", "OBJECT")}
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
    # the stub frame's product (its 300 s start at 23:00 gives a mid-time of 23:02:30)
    assert not list(
        (tmp_path / "out" / "spec").glob("highspec_Mg_main_20261001T230230_*_spec.fits")
    )

    with pytest.raises(ValueError, match="shape"):
        run_quicklook(folder, tmp_path / "out2", options(on_error="raise"))


def test_quicklook_survives_bad_calibration_frames(tmp_path, night):
    """A truncated bias, a flat binned 2x1 and a corrupt arc are each dropped and recorded
    as a failure of that frame; the masters are stacked from the other frames and the night
    still yields every science product. The command reports the failures (exit code 2)."""
    folder = tmp_path / "night"
    shutil.copytree(night.dir, folder)
    truth = night.truth
    bias, flat = (
        str(truth.frames.loc[truth.frames["exposure_type"] == kind, "file"].iloc[0])
        for kind in ("bias", "flat")
    )
    # a bias cut down to a stub, with the cards of a real one
    header = fits.getheader(night.dir / bias)
    stub = fits.PrimaryHDU(np.full((20, 30), 600, dtype=np.uint16), header)
    stub.writeto(folder / "sim_bias_bad.fits")
    # a real flat's image under another binning
    with fits.open(night.dir / flat) as hdul:
        data, header = hdul[0].data.copy(), hdul[0].header.copy()
    header["HBIN"] = 2
    fits.PrimaryHDU(data, header).writeto(folder / "sim_flat_bad.fits")
    # an arc that is not FITS at all: the scan reports it as an unreadable file
    (folder / "sim_arc_bad.fits").write_bytes(b"this is not a FITS file\n" * 200)

    out = tmp_path / "out"
    summary = run_quicklook(folder, out, options())
    failed = dict(summary.failed)
    assert set(failed) == {"sim_bias_bad.fits", "sim_flat_bad.fits", "sim_arc_bad.fits"}
    assert "shape" in failed["sim_bias_bad.fits"]
    assert "binning 2x1" in failed["sim_flat_bad.fits"]
    assert failed["sim_arc_bad.fits"] == "unreadable FITS file"
    assert len(spec_products(out)) == N_SCIENCE
    # the masters were stacked from the good frames only
    n_bias, n_flat = (
        int((truth.frames["exposure_type"] == kind).sum()) for kind in ("bias", "flat")
    )
    for kind, n in (("bias", n_bias), ("flat", n_flat)):
        _, prov = read_master(next((out / "calib").glob(f"*_{kind}.fits")))
        assert len(prov.inputs) == n, kind

    argv = ["quicklook", str(folder), "--out", str(out), "--instrument", "highspec"]
    assert main(argv + ["--layout", "survey45"]) == 2


def test_quicklook_drops_a_bad_full_arc_among_the_trace_candidates(tmp_path, night):
    """A night without a flat finds its traces on a full arc. A second full arc binned 2x1
    is a candidate that cannot be calibrated: it is recorded once in ``failed`` and left
    out, the traces still come from the good full arc and every science frame is reduced."""
    folder = tmp_path / "night"
    folder.mkdir()
    truth = night.truth
    copy_frames(night.dir, truth, folder, ["bias", "dark", "background", "arc", "science"])
    with fits.open(night.full_arc) as hdul:
        data, header = hdul[0].data.copy(), hdul[0].header.copy()
    header["HBIN"] = 2
    fits.PrimaryHDU(data, header).writeto(folder / "sim_arc_full_bad.fits")

    out = tmp_path / "out"
    summary = run_quicklook(folder, out, options())
    assert [source for source, _ in summary.failed] == ["sim_arc_full_bad.fits"]
    trace_qc = read_qc(next(p for p in summary.qc if p.name.endswith("_traces.json")))
    assert trace_qc.values["source"] == "arc_full"
    _, prov = read_traceset(next((out / "calib").glob("*_traces.fits")))
    assert prov.inputs == [file_checksum(night.full_arc)]
    assert len(spec_products(out)) == N_SCIENCE


# --------------------------------------------------------------- light frames by content


def test_resolve_light_sorts_frames_by_content(tmp_path, night):
    """The night as the control software writes it: every arc and science frame is
    ``IMAGETYP = 'light'`` with no lamp card, and there is no flat. The runner resolves the
    light frames from their content before the trace stage, so the traces come from the
    full arc it found, the arcs go to the arc stage and the science frames are reduced,
    each with its classification rule in the QC. A manifest row's type wins over the
    content: that frame never reaches the resolution."""
    folder = tmp_path / "night"
    folder.mkdir()
    truth = night.truth
    copy_frames(night.dir, truth, folder, ["bias", "dark", "background", "arc", "science"])
    for name in truth.frames.loc[truth.frames["exposure_type"].isin(["arc", "science"]), "file"]:
        with fits.open(folder / name, mode="update") as hdul:
            hdul[0].header["IMAGETYP"] = "light"
            del hdul[0].header["LAMP"]
    sources = science_files(truth)
    declared = sources[0]
    manifest = tmp_path / "manifest.csv"
    manifest.write_text(f"source,exposure_type\n{declared},science\n", encoding="utf-8", newline="")
    out = tmp_path / "out"
    summary = run_quicklook(folder, out, options(manifest=manifest))
    assert summary.failed == []

    trace_qc = read_qc(next(p for p in summary.qc if p.name.endswith("_traces.json")))
    assert trace_qc.values["source"] == "arc_full"
    # the full arc reached the arc stage (no reference yet), and so did the drift arcs
    assert any(night.full_arc.name in item and "run identify" in item for item in summary.skipped)
    for arc in night.cal_arcs:
        assert any(arc in item and "no reference" in item for item in summary.skipped)
    products = spec_products(out)
    assert sorted(read_spectra(p)[0].meta.source for p in products) == sorted(sources)
    rules = {}
    for path in products:
        es, prov, _ = read_spectra(path)
        assert es.wavelength_status == "nominal"
        qc = read_qc(out / "qa" / path.with_suffix(".json").name)
        rules[es.meta.source] = qc.values["classification_rule"]
        if es.meta.source != declared:
            assert any("content:continuum" in step for step in prov.steps)
    assert rules.pop(declared) == "manifest"
    assert set(rules.values()) == {"content:continuum"}
    log = next((out / "logs").glob("quicklook_*.log")).read_text(encoding="utf-8")
    assert "content:arc_like" in log


# ------------------------------------------------------------------------------- the CLI


def test_cli_quicklook_and_show(transferred_run, identified, night, out_dir, capsys):
    argv = ["quicklook", str(night.dir), "--out", str(out_dir), "--instrument", "highspec"]
    assert main(argv + ["--layout", "survey45"]) == 0
    out = capsys.readouterr().out
    assert "skipped" in out and "0 failed" in out
    # the command line reaches the same configuration (--layout picks the survey
    # parameters) and the same layout: everything is already reduced
    assert len(spec_products(out_dir)) == N_SCIENCE
    assert "already reduced" in out

    spec = spec_products(out_dir)[0]
    assert main(["show", str(spec)]) == 0
    out = capsys.readouterr().out
    assert "spectra" in out and "transferred" in out and identified.ref.identifier in out
    assert "science" in out and "sky" in out and "arc" in out
    first = night.truth.layout.dropna(subset=["telescope", "target"]).iloc[0]
    assert f"{first.telescope}={first.target}" in out  # the targets per telescope
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


def test_cli_solve_uses_the_library(sim_night_lab3, mg_lines, tmp_path, capsys):
    night_dir, truth = sim_night_lab3
    full = truth.frames.loc[truth.frames["arc_kind"] == "full"].iloc[0]
    arc = str(night_dir / str(full["file"]))
    out = tmp_path / "out"
    base = ["--out", str(out), "--instrument", "highspec", "--configuration", "Mg"]
    assert main(["solve", arc, *base]) == 2
    assert "run identify first" in capsys.readouterr().out
    cat_wl, _ = catalogue_arrays(mg_lines, usable_only=True)
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
