"""Tests for the pairing of science exposures with calibration-fiber arcs by time."""

import numpy as np
import pandas as pd
import pytest

from mastpipe.config.models import DriftParams
from mastpipe.core.model import PER_FIBER_COLUMNS, DriftModel, FrameMeta, ReferenceSolution
from mastpipe.pipeline.arcs import pair_arcs, status_for

T0 = pd.Timestamp("2026-10-01T22:00:00")


def meta(offset_minutes: float, exposure_type: str = "arc", configuration: str = "mgb518"):
    return FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type=exposure_type,
        configuration=configuration,
        mid_time=(T0 + pd.Timedelta(minutes=offset_minutes)).isoformat(timespec="seconds"),
        arc_kind="calfiber" if exposure_type == "arc" else None,
        checksum=f"sum{offset_minutes:+.0f}",
    )


def model(offset_minutes: float, usable: bool = True) -> DriftModel:
    return DriftModel(
        arc_checksums=(f"sum{offset_minutes:+.0f}",),
        arc_mid_time=(T0 + pd.Timedelta(minutes=offset_minutes)).isoformat(timespec="seconds"),
        per_fiber=pd.DataFrame({column: [] for column in PER_FIBER_COLUMNS}),
        shift_coeffs=np.array([offset_minutes / 100.0]),
        stretch_coeffs=np.array([0.0]),
        scatter_px=0.01,
        usable=usable,
        n_pix=1000,
    )


def arc(offset_minutes: float, **kwargs) -> tuple[FrameMeta, DriftModel]:
    configuration = kwargs.pop("configuration", "mgb518")
    return meta(offset_minutes, configuration=configuration), model(offset_minutes, **kwargs)


def reference(valid_from_minutes: float = -600.0) -> ReferenceSolution:
    valid_from = (T0 + pd.Timedelta(minutes=valid_from_minutes)).isoformat(timespec="seconds")
    return ReferenceSolution(
        instrument="highspec",
        configuration="mgb518",
        fiber_preset="survey45",
        solutions={},
        created=valid_from,
        valid_from=valid_from,
        source_checksums=(),
        identifier="highspec_mgb518_survey45_x",
    )


def mid_times(pairs) -> list[str]:
    return [drift.arc_mid_time for drift, _ in pairs]


def test_pairing_policies():
    exposure = meta(0.0, "science")
    arcs = [arc(-120.0), arc(-10.0), arc(15.0)]
    before = arcs[1][1].arc_mid_time
    after = arcs[2][1].arc_mid_time

    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="nearest_before"))
    assert mid_times(pairs) == [before] and pairs[0][1] == 1.0

    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="nearest"))
    assert mid_times(pairs) == [before] and pairs[0][1] == 1.0

    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="interpolate"))
    assert mid_times(pairs) == [before, after]
    assert [weight for _, weight in pairs] == pytest.approx([0.6, 0.4])
    assert status_for(pairs) == "transferred"


def test_nearest_before_falls_back_to_the_nearest_after():
    exposure = meta(0.0, "science")
    arcs = [arc(15.0), arc(40.0)]
    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="nearest_before"))
    assert mid_times(pairs) == [arcs[0][1].arc_mid_time] and pairs[0][1] == 1.0
    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="interpolate"))
    assert mid_times(pairs) == [arcs[0][1].arc_mid_time] and pairs[0][1] == 1.0


def test_nearest_prefers_the_closer_arc_after():
    exposure = meta(0.0, "science")
    arcs = [arc(-30.0), arc(5.0)]
    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="nearest"))
    assert mid_times(pairs) == [arcs[1][1].arc_mid_time]
    pairs = pair_arcs(exposure, arcs, reference(), DriftParams(time_policy="nearest_before"))
    assert mid_times(pairs) == [arcs[0][1].arc_mid_time]


def test_pairing_ignores_other_configuration_and_stale_arcs():
    exposure = meta(0.0, "science")
    params = DriftParams(time_policy="nearest", max_age_hours=3.0)
    good = arc(-100.0)
    arcs = [
        arc(-5.0, configuration="other"),  # another grating
        arc(-3.0, usable=False),  # the drift fit failed
        arc(-200.0),  # older than max_age_hours
        arc(250.0),  # too far after as well
        arc(-150.0),  # before the reference's valid_from
        good,
    ]
    pairs = pair_arcs(exposure, arcs, reference(valid_from_minutes=-120.0), params)
    assert mid_times(pairs) == [good[1].arc_mid_time]
    # an arc taken exactly at valid_from qualifies, one just before it does not
    pairs = pair_arcs(exposure, arcs, reference(valid_from_minutes=-100.0), params)
    assert mid_times(pairs) == [good[1].arc_mid_time]
    pairs = pair_arcs(exposure, arcs, reference(valid_from_minutes=-99.0), params)
    assert pairs == []


def test_no_candidates_gives_reference_status():
    exposure = meta(0.0, "science")
    assert pair_arcs(exposure, [], reference(), DriftParams()) == []
    assert status_for([]) == "reference"
    stale = [arc(-20.0 * 60.0)]
    assert pair_arcs(exposure, stale, reference(), DriftParams()) == []
    # an exposure without a mid-time cannot be paired
    no_time = FrameMeta(
        instrument="highspec", camera="main", exposure_type="science", configuration="mgb518"
    )
    assert pair_arcs(no_time, [arc(-10.0)], reference(), DriftParams()) == []


def test_times_are_compared_in_utc():
    exposure = meta(0.0, "science")
    aware = FrameMeta(
        instrument="highspec",
        camera="main",
        exposure_type="arc",
        configuration="mgb518",
        mid_time="2026-10-01T23:50:00+02:00",  # 21:50 UTC, ten minutes before
        arc_kind="calfiber",
    )
    pairs = pair_arcs(exposure, [(aware, model(-10.0))], reference(), DriftParams())
    assert len(pairs) == 1
