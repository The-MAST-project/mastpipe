"""Drift: the shift and stretch of the calibration fibers against a reference solution,
smoothed across the detector and transferred to every fiber.

Between full arcs only the calibration fibers see the lamp, in dedicated exposures a few
times a night. Each such arc gives, per arc fiber, a linear drift of the measured lines
against the reference's identified lines (:func:`fit_fiber_drift`): a global shift guess
by :func:`~mastpipe.steps.trace.comb_shift`, a one-to-one nearest match within
``tol_px``, and a clipped least-squares fit of

    x_obs - x_ref = shift + stretch (x_ref - n_pix / 2),

the stretch kept only when it is significant against its own least-squares error. Across
the arc fibers a polynomial in the trace centre ``y`` models the shift and the stretch
(:func:`fit_drift`), so that the drift of a science or sky fiber is read off at its own
``y``. :func:`apply_drift` undoes that drift: the observed column ``x`` of a fiber maps back
to the reference column ``x_ref = (x - s + e n / 2) / (1 + e)``, where the reference model
of the fiber gives the wavelength. Several drift models enter with weights, so that a
science exposure between two arcs can interpolate them. Fibers are keyed by their layout
id throughout (``Trace.index`` after :func:`~mastpipe.steps.trace.layout_indexed`); as a
guard against a mis-keyed fiber, a trace whose centre sits farther from the reference's
``trace_y_center`` than the caller's tolerance is not transferred but reported.

Everything here is pure: arrays, tables, dataclasses and a parameter model in, dataclasses
and arrays out; the time pairing of exposures with arcs lives in
:mod:`mastpipe.pipeline.arcs`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from mastpipe.config.models import DriftParams
from mastpipe.core.model import (
    PER_FIBER_COLUMNS,
    DriftModel,
    FrameMeta,
    ReferenceSolution,
    TraceSet,
    WavelengthSolution,
)
from mastpipe.steps.geometry import GratingGeometry
from mastpipe.steps.trace import comb_shift, robust_sigma
from mastpipe.steps.wavemodel import model_from_solution

CLIP_ITERATIONS = 5
"""Upper bound on the fit-and-clip rounds of :func:`fit_fiber_drift`."""
SIGMA_FLOOR_PX = 1e-6
"""Smallest robust sigma used for clipping, so a perfect fit does not clip everything."""


# ------------------------------------------------------------------------------ per fiber


def _match(x_obs: np.ndarray, x_ref: np.ndarray, guess: float, tol: float):
    """One-to-one nearest pairs ``(x_obs, x_ref)`` with ``|x_obs - guess - x_ref| <= tol``,
    closest pairs first."""
    distance = np.abs((x_obs[:, None] - guess) - x_ref[None, :])
    i_obs, i_ref = np.nonzero(distance <= tol)
    order = np.argsort(distance[i_obs, i_ref], kind="stable")
    used_obs = np.zeros(x_obs.size, dtype=bool)
    used_ref = np.zeros(x_ref.size, dtype=bool)
    pairs: list[tuple[int, int]] = []
    for k in order:
        i, j = int(i_obs[k]), int(i_ref[k])
        if used_obs[i] or used_ref[j]:
            continue
        used_obs[i] = used_ref[j] = True
        pairs.append((i, j))
    if not pairs:
        return np.array([]), np.array([])
    i, j = np.array(pairs).T
    return x_obs[i], x_ref[j]


def _linear_fit(t: np.ndarray, d: np.ndarray) -> tuple[float, float, float]:
    """Least squares of ``d = shift + stretch t`` on the kept points: ``(shift, stretch,
    sigma_stretch)``. One point gives its offset alone; ``sigma_stretch`` is infinite when
    the residual degrees of freedom are zero (two points fix the line exactly)."""
    if t.size < 2:
        return float(np.mean(d)), 0.0, np.inf
    design = np.vstack([np.ones_like(t), t]).T
    coeffs, *_ = np.linalg.lstsq(design, d, rcond=None)
    dof = t.size - 2
    if dof <= 0:
        return float(coeffs[0]), float(coeffs[1]), np.inf
    residual = d - design @ coeffs
    covariance = np.linalg.pinv(design.T @ design) * (float(residual @ residual) / dof)
    return float(coeffs[0]), float(coeffs[1]), float(np.sqrt(covariance[1, 1]))


def _result(shift: float, stretch: float, n_lines: int, rms: float, params: DriftParams):
    usable = n_lines >= params.min_lines and np.isfinite(rms) and rms <= params.max_scatter_px
    return {
        "shift_px": float(shift),
        "stretch": float(stretch),
        "n_lines": int(n_lines),
        "rms_px": float(rms),
        "usable": bool(usable),
    }


def fit_fiber_drift(
    x_obs: np.ndarray, reference: WavelengthSolution, params: DriftParams
) -> dict[str, Any]:
    """The shift and stretch of one fiber's measured lines against its reference solution.

    ``x_obs`` are the measured line columns of the fiber on a calibration-fiber arc and the
    predicted columns are the reference's identified lines ``reference.lines.x``. A global
    shift guess comes from :func:`~mastpipe.steps.trace.comb_shift` of the reference comb
    onto the observed one within ``params.max_shift_px``; measured and predicted lines are
    then matched one-to-one by distance within ``params.tol_px`` and ``x_obs - x_ref = shift
    + stretch (x_ref - n/2)`` is fitted by least squares, dropping residuals beyond
    ``params.clip_sigma`` robust sigmas and repeating until the set is stable. The stretch
    is set to zero (and the shift refitted as a plain offset) when ``|stretch| <
    params.stretch_significance`` times its least-squares error.

    Returns ``shift_px`` (positive when the lines moved to larger ``x``), ``stretch``,
    ``n_lines`` (the matched lines kept by the clipping), ``rms_px`` (their residual rms,
    ``nan`` without any) and ``usable`` (``n_lines >= params.min_lines`` and
    ``rms_px <= params.max_scatter_px``).
    """
    x_obs = np.asarray(x_obs, dtype=float).ravel()
    x_obs = x_obs[np.isfinite(x_obs)]
    x_ref = reference.lines["x"].to_numpy(dtype=float)
    x_ref = x_ref[np.isfinite(x_ref)]
    if x_obs.size == 0 or x_ref.size == 0:
        return _result(0.0, 0.0, 0, float("nan"), params)
    # the shift that lays the reference comb on the observed one: positive when the
    # lines moved to larger x, the sense of the fitted shift
    guess = comb_shift(x_ref, x_obs, params.max_shift_px, width=params.tol_px)
    matched_obs, matched_ref = _match(x_obs, x_ref, guess, params.tol_px)
    if matched_obs.size == 0:
        return _result(0.0, 0.0, 0, float("nan"), params)

    t = matched_ref - reference.n_pix / 2.0
    d = matched_obs - matched_ref
    keep = np.ones(d.size, dtype=bool)
    for _ in range(CLIP_ITERATIONS):
        shift, stretch, sigma_stretch = _linear_fit(t[keep], d[keep])
        residual = d - (shift + stretch * t)
        sigma = max(robust_sigma(residual[keep]), SIGMA_FLOOR_PX)
        new_keep = np.abs(residual) <= params.clip_sigma * sigma
        if new_keep.sum() == 0 or np.array_equal(new_keep, keep):
            break
        keep = new_keep
    else:
        shift, stretch, sigma_stretch = _linear_fit(t[keep], d[keep])

    if abs(stretch) < params.stretch_significance * sigma_stretch:
        stretch = 0.0
        shift = float(np.mean(d[keep]))
    residual = d[keep] - (shift + stretch * t[keep])
    rms = float(np.sqrt(np.mean(residual**2)))
    return _result(shift, stretch, int(keep.sum()), rms, params)


# ------------------------------------------------------------------- across the detector


def _y_centers(traceset: TraceSet) -> dict[int, float]:
    return {int(trace.index): float(trace.y_center) for trace in traceset.traces}


def _poly_across(y: np.ndarray, values: np.ndarray, degree: int) -> np.ndarray:
    """Polynomial of ``degree`` through ``values(y)`` in ``polyval`` order; a constant (no
    fit) for one point or degree zero, so a single fiber never reaches ``polyfit``."""
    if y.size == 0:
        return np.array([0.0])
    if y.size == 1 or degree <= 0:
        return np.array([float(np.mean(values))])
    return np.polyfit(y, values, degree)


def fit_drift(
    line_tables_by_fiber: dict[int, pd.DataFrame],
    reference: ReferenceSolution,
    traceset: TraceSet,
    arc_meta: FrameMeta,
    params: DriftParams,
) -> DriftModel:
    """The drift of one calibration-fiber arc against ``reference``, modelled across ``y``.

    ``line_tables_by_fiber`` holds the :func:`~mastpipe.steps.lines.find_lines` table of
    every lit arc fiber keyed by fiber index; the ``x_center`` of its ``use_x`` rows go
    through :func:`fit_fiber_drift` against the fiber's reference solution (a fiber the
    reference lacks cannot be fitted and is left out). The trace centre ``y_center`` of
    each fiber comes from ``traceset`` (the reference solution's when the trace is absent).
    Over the usable fibers a polynomial of degree ``min(params.degree, n_usable - 1)`` in
    ``y_center`` is fitted to the shifts and another to the stretches; ``scatter_px`` is
    the rms of the usable shifts about that model and the model is ``usable`` when at least
    one fiber is and ``scatter_px <= params.max_scatter_px``. Without any usable fiber both
    polynomials are zero.
    """
    y_centers = _y_centers(traceset)
    rows: list[dict[str, Any]] = []
    for fiber in sorted(line_tables_by_fiber):
        solution = reference.solutions.get(int(fiber))
        if solution is None:
            continue
        table = line_tables_by_fiber[fiber]
        use = table["use_x"].to_numpy(dtype=bool)
        x_obs = table.loc[use, "x_center"].to_numpy(dtype=float)
        result = fit_fiber_drift(x_obs, solution, params)
        rows.append(
            {
                "fiber": int(fiber),
                "y_center": y_centers.get(int(fiber), float(solution.trace_y_center)),
                **result,
            }
        )
    per_fiber = pd.DataFrame(rows, columns=list(PER_FIBER_COLUMNS))
    per_fiber = per_fiber.astype(
        {
            "fiber": int,
            "y_center": float,
            "shift_px": float,
            "stretch": float,
            "n_lines": int,
            "rms_px": float,
            "usable": bool,
        }
    )

    usable = per_fiber.loc[per_fiber["usable"]]
    n_usable = len(usable)
    y = usable["y_center"].to_numpy(dtype=float)
    degree = min(params.degree, n_usable - 1)
    shift_coeffs = _poly_across(y, usable["shift_px"].to_numpy(dtype=float), degree)
    stretch_coeffs = _poly_across(y, usable["stretch"].to_numpy(dtype=float), degree)
    if n_usable:
        about_model = usable["shift_px"].to_numpy(dtype=float) - np.polyval(shift_coeffs, y)
        scatter = float(np.sqrt(np.mean(about_model**2)))
    else:
        scatter = float("nan")
    n_pix = next(
        (solution.n_pix for solution in reference.solutions.values()), int(traceset.shape[1])
    )
    return DriftModel(
        arc_checksums=(arc_meta.checksum,),
        arc_mid_time=arc_meta.mid_time or "",
        per_fiber=per_fiber,
        shift_coeffs=np.asarray(shift_coeffs, dtype=float),
        stretch_coeffs=np.asarray(stretch_coeffs, dtype=float),
        scatter_px=scatter,
        usable=bool(n_usable >= 1 and scatter <= params.max_scatter_px),
        n_pix=int(n_pix),
    )


# ------------------------------------------------------------------------------- transfer


def apply_drift(
    reference: ReferenceSolution,
    drifts: list[tuple[DriftModel, float]],
    traceset: TraceSet,
    geom: GratingGeometry,
    tolerance: float | None = None,
    n_pix: int | None = None,
) -> tuple[dict[int, np.ndarray], dict[int, str]]:
    """The wavelength of every column of every fiber the reference solves, after the drift.

    ``drifts`` pairs drift models with weights (the time pairing of
    :func:`mastpipe.pipeline.arcs.pair_arcs`); at a fiber's trace centre ``y`` the shift is
    ``s = sum w_k shift_k(y)`` and the stretch ``e = sum w_k stretch_k(y)``, the observed
    column ``x`` maps back to the reference column ``x_ref = (x - s + e n / 2) / (1 + e)``
    and the fiber's reference model, rebuilt with ``geom``, gives ``wavelength(x_ref)``.
    Without any drift the result is the reference map itself. Returns ``(wavelengths,
    skipped)``: one array per transferred fiber of the reference, on the fiber's pixel
    grid ``0 .. n_pix - 1``, and the reason for every fiber that was not transferred.

    A fiber is skipped when ``tolerance`` is given and its trace in ``traceset`` sits more
    than ``tolerance`` rows from the solution's ``trace_y_center`` (the fiber is not the
    one the solution was made for); a fiber the trace set lacks is mapped at the
    reference's own centre. ``ValueError`` when ``n_pix`` is given and a solution's pixel
    grid differs from it (a reference from another detector crop), naming the reference.
    """
    y_centers = _y_centers(traceset)
    out: dict[int, np.ndarray] = {}
    skipped: dict[int, str] = {}
    for fiber, solution in reference.solutions.items():
        if n_pix is not None and int(solution.n_pix) != int(n_pix):
            raise ValueError(
                f"reference {reference.identifier} solves fiber {fiber} on {solution.n_pix} "
                f"pixels, the frame has {n_pix}: the reference is for another detector crop"
            )
        y = y_centers.get(int(fiber), float(solution.trace_y_center))
        offset = abs(y - float(solution.trace_y_center))
        if tolerance is not None and offset > tolerance:
            skipped[int(fiber)] = (
                f"trace centre y = {y:.2f} is {offset:.1f} px from the reference's "
                f"{solution.trace_y_center:.2f} (tolerance {tolerance:.1f} px)"
            )
            continue
        s = sum(weight * float(model.shift_at(y)) for model, weight in drifts)
        e = sum(weight * float(model.stretch_at(y)) for model, weight in drifts)
        n = solution.n_pix
        x = np.arange(n, dtype=float)
        x_ref = (x - s + e * n / 2.0) / (1.0 + e)
        wavelength = np.asarray(model_from_solution(solution, geom).wavelength(x_ref), float)
        out[int(fiber)] = wavelength
    return out, skipped
