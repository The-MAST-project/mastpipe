"""Unit tests of the quick-look runner's pure helpers: the drift models a frame may be
paired with (those fitted against its own reference), the blend of two drift models
under the ``interpolate`` policy, the group's binning and the per-frame isolation of the
calibration frames."""

import shutil
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from mastpipe.config.models import DriftParams
from mastpipe.core.model import (
    PER_FIBER_COLUMNS,
    CalibrationSet,
    DriftModel,
    FrameMeta,
    FrameRef,
    ReferenceSolution,
)
from mastpipe.instruments.registry import load_instrument
from mastpipe.io.products import read_master, write_master
from mastpipe.pipeline import masters
from mastpipe.pipeline.arcs import pair_arcs
from mastpipe.pipeline.callib import CalibLibrary
from mastpipe.pipeline.group import _calibrated, _group_binning, _group_readout
from mastpipe.pipeline.masters import _masters_hash, _usable_frames
from mastpipe.pipeline.runner import RunOptions, run_quicklook
from mastpipe.pipeline.science import _arc_weights, _arcs_for, _blend_drifts
from mastpipe.vocab import GRATING_NAMES

N_PIX = 1598


def make_meta(mid_time: str, checksum: str, exposure_type: str = "arc") -> FrameMeta:
    return FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type=exposure_type,
        configuration="Mg",
        mid_time=mid_time,
        checksum=checksum,
        source=f"{checksum}.fits",
    )


def make_reference(identifier: str, valid_from: str) -> ReferenceSolution:
    return ReferenceSolution(
        instrument="highspec",
        configuration="Mg",
        fiber_preset="survey45",
        solutions={},
        created=valid_from,
        valid_from=valid_from,
        source_checksums=(identifier,),
        identifier=identifier,
        hash="h" + identifier,
    )


def make_drift(
    checksum: str,
    mid_time: str,
    shift_coeffs: list[float],
    stretch_coeffs: list[float],
    fibers: dict[int, tuple[float, float, int, float, bool]],
    scatter_px: float = 0.02,
    usable: bool = True,
) -> DriftModel:
    rows = [
        {
            "fiber": fiber,
            "y_center": 10.0 + 8.5 * fiber,
            "shift_px": shift,
            "stretch": stretch,
            "n_lines": n_lines,
            "rms_px": rms,
            "usable": ok,
        }
        for fiber, (shift, stretch, n_lines, rms, ok) in fibers.items()
    ]
    return DriftModel(
        arc_checksums=(checksum,),
        arc_mid_time=mid_time,
        per_fiber=pd.DataFrame(rows, columns=list(PER_FIBER_COLUMNS)),
        shift_coeffs=np.array(shift_coeffs),
        stretch_coeffs=np.array(stretch_coeffs),
        scatter_px=scatter_px,
        usable=usable,
        n_pix=N_PIX,
    )


# ------------------------------------------------------------- arcs tied to their reference


@pytest.mark.parametrize("policy", ["nearest", "interpolate", "nearest_before"])
def test_frame_never_receives_drift_fitted_against_another_reference(policy):
    """A full arc solved mid-night makes reference B; the arcs after it are fitted against
    B. A science frame before that time uses A and must pair with A's arcs only, however
    close B's arcs are in time."""
    ref_a = make_reference("A", "2026-10-01T18:00:00")
    ref_b = make_reference("B", "2026-10-01T21:00:00")
    fibers = {0: (0.3, 0.0, 20, 0.02, True), 44: (0.3, 0.0, 20, 0.02, True)}
    arc_a1 = make_drift("a1", "2026-10-01T19:00:00", [0.3], [0.0], fibers)
    arc_b1 = make_drift("b1", "2026-10-01T21:10:00", [-2.0], [0.0], fibers)
    arc_b2 = make_drift("b2", "2026-10-01T21:30:00", [-2.0], [0.0], fibers)
    arcs = [
        (make_meta("2026-10-01T19:00:00", "a1"), arc_a1, "A"),
        (make_meta("2026-10-01T21:10:00", "b1"), arc_b1, "B"),
        (make_meta("2026-10-01T21:30:00", "b2"), arc_b2, "B"),
    ]
    # the frame sits just before B's validity, between an A arc and the B arcs
    frame = make_meta("2026-10-01T20:55:00", "sci", exposure_type="science")
    params = DriftParams(time_policy=policy)

    mine = _arcs_for(arcs, ref_a)
    assert [meta.checksum for meta, _ in mine] == ["a1"]
    pairs = pair_arcs(frame, mine, ref_a, params)
    assert [model.arc_checksums for model, _ in pairs] == [("a1",)]
    assert pairs[0][1] == pytest.approx(1.0)
    # without the filter the nearest / interpolate policies would have taken the B arc
    # fifteen minutes away (nearest_before alone is safe: no B arc precedes the frame)
    unfiltered = pair_arcs(frame, [(m, d) for m, d, _ in arcs], ref_a, params)
    taken = sorted(c for model, _ in unfiltered for c in model.arc_checksums)
    assert taken == (["a1"] if policy == "nearest_before" else ["a1", "b1"] if policy == "interpolate" else ["b1"])
    # a frame after B's time uses B and only B's arcs
    later = make_meta("2026-10-01T21:20:00", "sci2", exposure_type="science")
    theirs = _arcs_for(arcs, ref_b)
    assert sorted(meta.checksum for meta, _ in theirs) == ["b1", "b2"]
    pairs = pair_arcs(later, theirs, ref_b, params)
    assert pairs and all(c.startswith("b") for model, _ in pairs for c in model.arc_checksums)


# --------------------------------------------------------------------------- the blend


def test_blend_drifts_single_model_is_itself():
    drift = make_drift("x", "2026-10-01T20:00:00", [0.1], [0.0], {0: (0.1, 0.0, 10, 0.01, True)})
    blended, weights = _blend_drifts([(drift, 1.0)])
    assert blended is drift and weights is None


def test_blend_drifts_weights_two_models():
    fibers_a = {0: (0.2, 1e-5, 20, 0.02, True), 20: (0.4, 1e-5, 18, 0.03, True)}
    fibers_b = {0: (1.2, 3e-5, 15, 0.04, True), 44: (1.6, 3e-5, 16, 0.01, False)}
    first = make_drift("aaa", "2026-10-01T20:00:00", [0.2], [1e-5], fibers_a, scatter_px=0.02)
    second = make_drift("bbb", "2026-10-01T21:00:00", [2e-3, 1.2], [0.0, 3e-5], fibers_b, scatter_px=0.05)
    blended, weights = _blend_drifts([(first, 0.3), (second, 0.7)])
    for y in (10.0, 300.0):
        assert blended.shift_at(y) == pytest.approx(0.3 * first.shift_at(y) + 0.7 * second.shift_at(y))
        assert blended.stretch_at(y) == pytest.approx(0.3 * first.stretch_at(y) + 0.7 * second.stretch_at(y))
    # right-aligned coefficients: the constant model contributes to the last slot only
    assert blended.shift_coeffs == pytest.approx([0.7 * 2e-3, 0.3 * 0.2 + 0.7 * 1.2])
    assert blended.stretch_coeffs == pytest.approx([0.0, 0.3 * 1e-5 + 0.7 * 3e-5])
    assert blended.arc_checksums == ("aaa", "bbb")
    assert blended.arc_mid_time == "2026-10-01T20:42:00"  # 0.7 of the hour after 20:00
    assert blended.n_pix == N_PIX
    assert blended.scatter_px == pytest.approx(0.05)
    assert blended.usable is True
    # the fewest usable arc fibers of the blended models, not a count of the merged table
    assert (first.n_usable, second.n_usable, blended.n_usable) == (2, 1, 1)
    table = blended.per_fiber.set_index("fiber")
    assert sorted(table.index) == [0, 20, 44]
    assert table.loc[0, "shift_px"] == pytest.approx(0.3 * 0.2 + 0.7 * 1.2)
    assert table.loc[0, "stretch"] == pytest.approx(0.3 * 1e-5 + 0.7 * 3e-5)
    assert table.loc[0, "n_lines"] == 15 and table.loc[0, "usable"]
    assert table.loc[0, "rms_px"] == pytest.approx(0.3 * 0.02 + 0.7 * 0.04)
    # a fiber only one model holds takes that model's values and its weight
    assert table.loc[20, "shift_px"] == pytest.approx(0.4)
    assert table.loc[44, "shift_px"] == pytest.approx(1.6) and not table.loc[44, "usable"]
    assert weights == {0: pytest.approx(1.0), 20: pytest.approx(0.3), 44: pytest.approx(0.7)}
    # unnormalised weights are normalised; an unusable model makes the blend unusable
    same, _ = _blend_drifts([(first, 3.0), (second, 7.0)])
    assert same.shift_at(100.0) == pytest.approx(blended.shift_at(100.0))
    bad = make_drift("ccc", "2026-10-01T22:00:00", [0.5], [0.0], fibers_a, usable=False)
    assert not _blend_drifts([(first, 0.5), (bad, 0.5)])[0].usable
    assert _arc_weights([(first, 0.3), (second, 0.7)], ["aaa", "bbb"]) == pytest.approx([0.3, 0.7])


# ------------------------------------------------------------- no silent configuration


def test_group_without_configuration_fails_naming_the_gratings(sim_night, tmp_path):
    """A frame whose grating is neither in a GRATING card nor in its file name is not
    reduced with the instrument's only grating: the group fails and names the choices,
    the instrument's configurations (not every fleet grating name, some of which it
    cannot reduce)."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    del header["GRATING"]
    fits.PrimaryHDU(data, header).writeto(folder / "frame_0001.fits")
    summary = run_quicklook(folder, tmp_path / "out", RunOptions("highspec"))
    message = "configuration unknown; pass --configuration one of ['Mg']"
    assert summary.failed == [("group main/None", message)]
    assert sorted(load_instrument("highspec").configurations) == ["Mg"]
    assert set(GRATING_NAMES) - {"Mg"}  # names the fleet knows but highspec cannot reduce
    assert summary.products == [] and summary.qc == []


# ------------------------------------------------------------------ assumed time zone


def test_one_warning_per_run_when_frame_times_carry_no_zone(sim_night, tmp_path):
    """The simulator writes a naive DATE-OBS, which HighSpec's map takes as UTC: the run log
    says so once, however many frames there are, and not at all when every stamp has its
    own offset."""
    sim_dir, _ = sim_night
    warning = "frame times carry no time zone; taken as UTC (header map 'assume-utc')"
    for name, offset in (("naive", ""), ("explicit", "+00:00")):
        folder = tmp_path / name
        folder.mkdir()
        for k, source in enumerate(("sim_bias_00.fits", "sim_bias_01.fits")):
            with fits.open(sim_dir / source) as hdul:
                data = hdul[0].data.copy()
                header = hdul[0].header.copy()
            del header["GRATING"]  # the group then fails fast, after the scan
            header["DATE-OBS"] = header["DATE-OBS"] + offset
            fits.PrimaryHDU(data, header).writeto(folder / f"frame_{k:04d}.fits")
        out = tmp_path / f"out_{name}"
        run_quicklook(folder, out, RunOptions("highspec"))
        (log_file,) = (out / "logs").glob("*.log")
        count = log_file.read_text(encoding="utf-8").count(warning)
        assert count == (1 if name == "naive" else 0), name


def test_one_warning_per_run_for_products_of_an_earlier_version(sim_night, tmp_path):
    """Product names gained a checksum segment, so a product written by an earlier version
    is never found at its new path: a new file is written beside it. The run log says so
    once, naming the ``spec`` folder, however many such files there are, and not at all
    when every product name carries the segment."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data, header = hdul[0].data.copy(), hdul[0].header.copy()
    del header["GRATING"]  # the group then fails fast, after the check
    fits.PrimaryHDU(data, header).writeto(folder / "frame_0001.fits")
    warning = "written by an earlier mastpipe version"
    for name, files in (
        (
            "old",
            [
                "highspec_Mg_main_20261001T180230_spec.fits",
                "highspec_Mg_main_notime_spec.fits",
                "highspec_Mg_main_20261001T180230_comb_target_a.fits",
            ],
        ),
        (
            "new",
            [
                "highspec_Mg_main_20261001T180230_a1b2c3_spec.fits",
                "highspec_Mg_main_notime_nosum_spec.fits",
                "highspec_Mg_main_20261001T180230_a1b2c3_comb_target_a.fits",
            ],
        ),
    ):
        out = tmp_path / f"out_{name}"
        (out / "spec").mkdir(parents=True)
        for file in files:
            (out / "spec" / file).write_bytes(b"")
        run_quicklook(folder, out, RunOptions("highspec"))
        (log_file,) = (out / "logs").glob("*.log")
        lines = [line for line in log_file.read_text(encoding="utf-8").splitlines() if warning in line]
        assert len(lines) == (1 if name == "old" else 0), name
        if lines:
            assert "WARNING" in lines[0] and str(out / "spec") in lines[0] and "3 product" in lines[0]


def test_one_warning_per_run_for_an_amplifier_without_its_noise(sim_night, tmp_path):
    """HighSpec lists gain and read noise for ``conventional`` and ``em``; frames naming
    another amplifier (here a misspelt one) take the camera's defaults, and the run log
    says so once, however many frames there are, and not at all for a listed amplifier."""
    sim_dir, _ = sim_night
    warning = "no gain and read noise listed for amplifier"
    for amplifier in ("Conventionl", "Conventional"):
        folder = tmp_path / amplifier
        folder.mkdir()
        for k, source in enumerate(("sim_bias_00.fits", "sim_bias_01.fits")):
            with fits.open(sim_dir / source) as hdul:
                data = hdul[0].data.copy()
                header = hdul[0].header.copy()
            del header["GRATING"]  # the group then fails fast, after the scan
            header["OUTPTAMP"] = amplifier
            fits.PrimaryHDU(data, header).writeto(folder / f"frame_{k:04d}.fits")
        out = tmp_path / f"out_{amplifier}"
        run_quicklook(folder, out, RunOptions("highspec"))
        (log_file,) = (out / "logs").glob("*.log")
        text = log_file.read_text(encoding="utf-8")
        assert text.count(warning) == (1 if amplifier == "Conventionl" else 0), amplifier
        if amplifier == "Conventionl":
            assert "conventionl" in text


# ---------------------------------------------------- other instruments, unmatched rows


def test_other_instruments_and_unmatched_manifest_rows_are_skipped(sim_night, tmp_path):
    """A nested night holding a DeepSpec frame next to HighSpec's: the DeepSpec frame and
    a manifest row naming no frame are listed in ``skipped`` with their reasons."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    with fits.open(sim_dir / "sim_bias_00.fits") as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    del header["GRATING"]  # the group then fails fast, after the scan
    for name, instrument in (("highspec", "Highspec"), ("deepspec", "Deepspec")):
        (folder / name / "acquisition-0001").mkdir(parents=True)
        header["INSTRUME"] = instrument
        fits.PrimaryHDU(data, header).writeto(folder / name / "acquisition-0001" / "e-001.fits")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("source,target\nnothing.fits,x\n", encoding="utf-8")
    options = RunOptions("highspec", manifest=manifest)
    summary = run_quicklook(folder, tmp_path / "out", options)
    assert summary.skipped == [
        "deepspec/acquisition-0001/e-001.fits: instrument 'deepspec' is not 'highspec'",
        "manifest row 0 (nothing.fits): matches no frame",
    ]
    assert [source for source, _ in summary.failed] == ["group main/None"]


def test_manifest_row_naming_a_non_fits_file_has_its_own_reason(sim_night, tmp_path):
    """A row naming a file that exists but is not a FITS file name is skipped as such, not
    as matching no frame."""
    sim_dir, _ = sim_night
    folder = tmp_path / "night"
    folder.mkdir()
    shutil.copy(sim_dir / "sim_bias_00.fits", folder / "sim_bias_00.fits")
    (folder / "notes.txt").write_text("bias taken with the dome closed\n", encoding="utf-8")
    manifest = tmp_path / "manifest.csv"
    manifest.write_text("source,lamp_state\nnotes.txt,off\nnothing.fits,off\n", encoding="utf-8")
    summary = run_quicklook(folder, tmp_path / "out", RunOptions("highspec", manifest=manifest))
    assert "manifest row 0 (notes.txt): not a FITS file (.fits, .fit, .fits.fz)" in summary.skipped
    assert "manifest row 1 (nothing.fits): matches no frame" in summary.skipped
    assert not any(row.startswith("manifest row 0") and "matches no frame" in row for row in summary.skipped)


# ------------------------------------------------------------ calibration frames, binning


def binned_ref(name: str, exposure_type: str, binning: tuple[int, int], **fields) -> FrameRef:
    meta = FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type=exposure_type,
        configuration="Mg",
        binning=binning,
        exposure_time=1.0,
        source=name,
        **fields,
    )
    return FrameRef(path=Path(name), extension=0, meta=meta)


def test_group_binning_follows_the_science_frames():
    """Many calibration frames at 1x1 do not outvote the science frames at 2x2 (light
    frames, which may be science, count with them); without any, all frames vote."""
    calibrations = [binned_ref(f"bias_{k}.fits", "bias", (1, 1)) for k in range(5)]
    calibrations += [binned_ref(f"flat_{k}.fits", "flat", (1, 1)) for k in range(3)]
    science = [
        binned_ref("sci_0.fits", "science", (2, 2)),
        binned_ref("sci_1.fits", "science", (1, 1)),
        binned_ref("light_0.fits", "light", (2, 2)),
    ]
    assert _group_binning(calibrations + science) == (2, 2)
    odd = [binned_ref("bias_odd.fits", "bias", (2, 1))]
    assert _group_binning(odd + calibrations) == (1, 1)


def test_group_readout_follows_the_science_frames():
    """The group's readout key and amplifier come from the frames that may be science, as
    its binning does; without any, all frames vote."""
    conventional = {"readout_mode": "image", "amplifier": "conventional"}
    em = {"readout_mode": "image", "amplifier": "em"}
    calibrations = [binned_ref(f"bias_{k}.fits", "bias", (1, 1), **conventional) for k in range(5)]
    science = [
        binned_ref("sci_0.fits", "science", (1, 1), **em),
        binned_ref("light_0.fits", "light", (1, 1), **em),
        binned_ref("sci_1.fits", "science", (1, 1), **conventional),
    ]
    assert _group_readout(calibrations + science) == ("image/em", "em")
    assert _group_readout(calibrations) == ("image/conventional", "conventional")
    assert _group_readout([binned_ref("x.fits", "bias", (1, 1))]) == ("unknown/unknown", None)


def usable_group(on_error: str = "continue") -> SimpleNamespace:
    """The parts of a ``_Group`` that :func:`_usable_frames` reads."""
    detector = SimpleNamespace(shape=(8, 10), crop=None, transpose=False, flip_x=False, flip_y=False)
    return SimpleNamespace(
        effective=SimpleNamespace(detector=detector),
        binning=(1, 1),
        readout_key="unknown/unknown",
        options=SimpleNamespace(on_error=on_error),
        summary=SimpleNamespace(failed=[], skipped=[]),
    )


def test_usable_frames_drops_bad_frames_with_their_reason(tmp_path):
    """A frame that cannot be read, has another binning, or another shape than the
    detector's (after crop and orientation) is recorded in ``failed`` and dropped; the
    others are kept in order. Under ``on_error='raise'`` the first bad frame raises."""
    refs = []
    for name, shape, binning in (
        ("a.fits", (8, 10), (1, 1)),
        ("binned.fits", (8, 10), (2, 1)),
        ("b.fits", (8, 10), (1, 1)),
        ("stub.fits", (4, 5), (1, 1)),
        ("c.fits", (8, 10), (1, 1)),
    ):
        fits.PrimaryHDU(np.full(shape, 600, dtype=np.uint16)).writeto(tmp_path / name)
        ref = binned_ref(name, "bias", binning)
        refs.append(FrameRef(path=tmp_path / name, extension=0, meta=ref.meta))
    (tmp_path / "junk.fits").write_bytes(b"this is not a FITS file\n" * 200)
    refs.insert(1, FrameRef(tmp_path / "junk.fits", 0, binned_ref("junk.fits", "bias", (1, 1)).meta))

    group = usable_group()
    frames = _usable_frames(group, refs, calibrate=False)
    assert [frame.meta.source for frame in frames] == ["a.fits", "b.fits", "c.fits"]
    failed = dict(group.summary.failed)
    assert list(failed) == ["junk.fits", "binned.fits", "stub.fits"]
    assert "binning 2x1" in failed["binned.fits"] and "1x1" in failed["binned.fits"]
    assert "shape (4, 5)" in failed["stub.fits"] and "(8, 10)" in failed["stub.fits"]

    # a stub sorting first does not set the shape: the good frames after it survive
    group = usable_group()
    frames = _usable_frames(group, [refs[4], refs[0], refs[3]], calibrate=False)
    assert [frame.meta.source for frame in frames] == ["a.fits", "b.fits"]
    assert [source for source, _ in group.summary.failed] == ["stub.fits"]

    with pytest.raises(ValueError, match="binning 2x1"):
        _usable_frames(usable_group("raise"), refs[2:], calibrate=False)


def test_usable_frames_skips_a_frame_of_another_readout(tmp_path):
    """A calibration frame read through another amplifier than the group's is recorded in
    ``skipped`` with both readout keys and left out of the master: it is not a failure
    (it may belong to another group's night), and it does not fail the stack or raise
    under ``on_error='raise'``."""
    refs = []
    for name, amplifier in (
        ("a.fits", "conventional"),
        ("em.fits", "em"),
        ("b.fits", "conventional"),
    ):
        fits.PrimaryHDU(np.full((8, 10), 600, dtype=np.uint16)).writeto(tmp_path / name)
        meta = binned_ref(name, "bias", (1, 1), readout_mode="image", amplifier=amplifier).meta
        refs.append(FrameRef(path=tmp_path / name, extension=0, meta=meta))
    group = usable_group()
    group.readout_key = "image/conventional"
    frames = _usable_frames(group, refs, calibrate=False)
    assert [frame.meta.source for frame in frames] == ["a.fits", "b.fits"]
    assert group.summary.failed == []
    assert group.summary.skipped == ["em.fits: readout image/em is not the group's image/conventional"]
    group = usable_group("raise")
    group.readout_key = "image/conventional"
    assert [frame.meta.source for frame in _usable_frames(group, refs, calibrate=False)] == ["a.fits", "b.fits"]
    assert len(group.summary.skipped) == 1


def calibrating_group(t_background: float | None) -> SimpleNamespace:
    """A group whose calibration set holds a bias and a background master, the latter
    with ``t_background`` as its exposure time, and background subtraction switched on."""
    group = usable_group()
    group.effective.detector.saturation = 65000.0
    group.params = SimpleNamespace(
        calib=SimpleNamespace(subtract_background=True, flat_subtract_background=True, cosmic_sigma=8.0, cosmic_size=5)
    )
    group.calib = CalibrationSet(
        bias=np.full((8, 10), 100.0, dtype=np.float32),
        dark_rate=0.0,
        background=np.full((8, 10), 130.0, dtype=np.float32),
        t_background=t_background,
        gain=1.0,
        read_noise=5.0,
        shape=(8, 10),
    )
    return group


def test_background_master_without_exposure_time_is_not_subtracted(tmp_path):
    """A background master whose inputs carried no exposure time cannot be scaled: the
    frames fall back to the bias (dark) path instead of every frame failing with
    'background subtraction requested ...'."""
    fits.PrimaryHDU(np.full((8, 10), 600, dtype=np.uint16)).writeto(tmp_path / "flat.fits")
    meta = binned_ref("flat.fits", "flat", (1, 1)).meta
    meta = FrameMeta(**{**meta.__dict__, "exposure_time": 2.0})
    ref = FrameRef(path=tmp_path / "flat.fits", extension=0, meta=meta)

    group = calibrating_group(t_background=None)
    frames = _usable_frames(group, [ref], calibrate=True)
    assert [frame.meta.source for frame in frames] == ["flat.fits"]
    assert group.summary.failed == []
    assert np.allclose(frames[0].data, 500.0)  # bias only

    frame, steps = _calibrated(group, ref)
    assert steps == ["bias"]
    assert np.allclose(frame.data, 500.0)

    # with an exposure time the background is scaled and subtracted as before
    group = calibrating_group(t_background=4.0)
    frame, steps = _calibrated(group, ref)
    assert steps[0] == "bias" and steps[1].startswith("background x 0.5")
    assert np.allclose(frame.data, 500.0 - 30.0 * 0.5)


# ------------------------------------------------------------------- masters and library

CALIB_OPTIONS = RunOptions("highspec", layout="survey45")


def calibration_night(sim_dir: Path, folder: Path, kinds: tuple[str, ...] = ("bias", "flat")) -> Path:
    """A night folder holding only the simulated night's frames of ``kinds``: the masters
    and the traces are reduced, nothing else."""
    folder.mkdir()
    for kind in kinds:
        for path in sorted(sim_dir.glob(f"sim_{kind}_*.fits")):
            shutil.copy(path, folder / path.name)
    return folder


def library_master(out: Path, kind: str) -> Path:
    (entry,) = CalibLibrary(out / "calib").list(kind=kind)
    return out / "calib" / entry["path"]


def test_library_flat_is_reused_only_with_the_same_masters_content(sim_night, tmp_path):
    """A library flat records the hash of the checksums of the masters it was calibrated
    with (``MPMKH``) and is reused only when the group's masters hash the same: a bias
    rewritten under the same name rebuilds it, and so does, once, a flat written before
    the hash was recorded. The trace set built on the flat carries the checksum of the
    flat's first input frame in its name, as the flat does, not ``nosum``."""
    folder = calibration_night(sim_night[0], tmp_path / "night")
    out = tmp_path / "out"
    first = run_quicklook(folder, out, CALIB_OPTIONS)
    assert first.failed == []
    bias_path, flat_path = library_master(out, "bias"), library_master(out, "flat")
    bias, _ = read_master(bias_path)
    flat, flat_prov = read_master(flat_path)
    assert flat_prov.masters_hash == _masters_hash({"bias": (bias, bias_path)})
    (traces,) = (out / "calib").glob("*_traces.fits")
    assert traces.name.split("_")[4] == flat_path.name.split("_")[4] == flat.provenance["inputs"][0][:6]

    # nothing changed: the flat is reused
    assert flat_path not in run_quicklook(folder, out, CALIB_OPTIONS).products

    # another bias under the same file name (and the same index entry): the bias is
    # reused as it stands, the flat calibrated with the old one is built again
    _, bias_prov = read_master(bias_path)
    write_master(bias_path, replace(bias, data=bias.data + 1.0), bias_prov)
    changed = run_quicklook(folder, out, CALIB_OPTIONS)
    assert changed.failed == []
    assert bias_path not in changed.products and flat_path in changed.products
    assert read_master(flat_path)[1].masters_hash == _masters_hash({"bias": (bias, bias_path)})
    assert flat_path not in run_quicklook(folder, out, CALIB_OPTIONS).products

    # a flat written before the masters hash was recorded is built again once
    with fits.open(flat_path, mode="update") as hdul:
        del hdul[0].header["MPMKH"]
    assert flat_path in run_quicklook(folder, out, CALIB_OPTIONS).products
    assert flat_path not in run_quicklook(folder, out, CALIB_OPTIONS).products
    # every run, rebuild or reuse, wrote the flat's trace set under the one name
    assert list((out / "calib").glob("*_traces.fits")) == [traces]


def test_flat_batch_is_chosen_among_the_frames_that_survive(sim_night, tmp_path, monkeypatch):
    """The flats are batched by exposure time and the batch with the most frames that load
    and calibrate is stacked, as for the other kinds: three 5 s flats of which two fail
    their calibration lose to two 8 s flats. A flat whose calibration raises is recorded
    in ``failed`` with the error."""
    sim_dir, _ = sim_night
    folder = calibration_night(sim_dir, tmp_path / "night")
    for k, source in enumerate(("sim_flat_00.fits", "sim_flat_01.fits")):
        with fits.open(sim_dir / source) as hdul:
            data, header = hdul[0].data.copy(), hdul[0].header.copy()
        header["EXPOSURE"] = 8.0
        header["DATE-OBS"] = f"2026-10-01T20:0{k}:00"
        fits.PrimaryHDU(data, header).writeto(folder / f"sim_flat_long_{k}.fits")
    real = masters.calibrate_frame

    def calibrate(frame, *args, **kwargs):
        if frame.meta.source in ("sim_flat_00.fits", "sim_flat_01.fits"):
            raise RuntimeError("calibration refused")
        return real(frame, *args, **kwargs)

    monkeypatch.setattr(masters, "calibrate_frame", calibrate)
    out = tmp_path / "out"
    summary = run_quicklook(folder, out, CALIB_OPTIONS)
    assert summary.failed == [
        ("sim_flat_00.fits", "RuntimeError: calibration refused"),
        ("sim_flat_01.fits", "RuntimeError: calibration refused"),
    ]
    flat, prov = read_master(library_master(out, "flat"))
    assert flat.validity["exposure_time"] == 8.0
    assert flat.n_frames == 2 and len(prov.inputs) == 2
