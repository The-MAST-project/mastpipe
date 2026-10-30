"""Full-arc wavelength solve: a port of the notebook's ``solve_trace`` and its helpers.

Every few months all fibers see the lamp and each gets a fresh, independent wavelength
solution. Per fiber the solve is the notebook's identify -> fit -> clip loop
(:func:`solve_fiber`): the current grating model predicts the wavelength of every measured
line (:func:`identify_lines`, nearest catalogue line within a tolerance, rejected when the
runner-up is close), :func:`~mastpipe.steps.wavemodel.fit_model` refits ``x_c``, ``B``
[, ``alpha``] and the distortion term to the matches, residuals beyond ``clip_sigma``
robust sigmas are dropped and the fit repeated, and the loop stops when the wavelength map
moves by less than ``converge_px`` anywhere on the detector and the line set no longer
changes, or after ``n_passes``. The identification always runs on the grating model, the
only one that can say where a line should be before any line is identified.

:func:`solve_full_arc` runs it on every fiber of a full arc. The first time, a single
anchored model of one fiber seeds the others after sliding it by the cross-correlation
shift of their continuum-subtracted spectra (:func:`shift_to_reference`: the fibers see the
same lamp through the same grating, so their spectra are the same up to a shift). Later, a
previous :class:`~mastpipe.core.model.ReferenceSolution` seeds every fiber from its own
earlier model. The result is a :class:`~mastpipe.core.model.WavelengthSolution` per fiber,
packaged by :func:`build_reference`; :func:`fiber_agreement` is the notebook's external
accuracy check (the same line, fiber to fiber) and :func:`wavelength_map` evaluates a
stored solution on its pixel grid.

Everything here is pure: arrays, tables and parameter models in, dataclasses and tables
out; nothing prints or touches a file.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

from mastpipe.config.models import WavecalParams
from mastpipe.core.model import (
    LINE_COLUMNS,
    ReferenceSolution,
    TraceSet,
    WavelengthSolution,
    make_identifier,
)
from mastpipe.steps.geometry import GratingGeometry
from mastpipe.steps.trace import robust_sigma
from mastpipe.steps.wavemodel import GratingModel, WavelengthModel, fit_model, model_from_solution

FINAL_MODELS = ("grating", "poly")
"""What :func:`solve_fiber` may fit to the identified lines; only ``"grating"`` is implemented."""
SIGMA_FLOOR_PX = 1e-6
"""Smallest robust sigma used for clipping, so a perfect fit does not clip everything."""
LINE_SET_DECIMALS = 6
"""Wavelength rounding (nm) when comparing the identified line sets of two passes."""
HASH_DECIMALS = 9
"""Rounding of the parameters and line positions that enter the reference hash."""
HASH_LENGTH = 12
"""Hex digits of the reference hash."""
AGREEMENT_COLUMNS = ("wavelength", "n_fibers", "spread_px")
GEOMETRY_KEYS = ("d_nm", "alpha_rad", "order", "B_nom")
"""The grating constants a :class:`WavelengthSolution` records for provenance."""


def _empty_lines() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "x": pd.Series(dtype=float),
            "wavelength": pd.Series(dtype=float),
            "species": pd.Series(dtype=object),
            "dx_px": pd.Series(dtype=float),
        }
    )


# ------------------------------------------------------------------------ identification


def identify_lines(
    x_obs, model: WavelengthModel, cat_wl, cat_sp, tol_px: float, unique_tol_px: float
) -> pd.DataFrame:
    """Match measured line columns to catalogue wavelengths through ``model``.

    Every column of ``x_obs`` is mapped to a wavelength and matched to the nearest line of
    ``cat_wl`` (sorted, air nm; ``cat_sp`` its species). The match is kept when the offset
    ``dx_px = (lam(x) - cat) / dispersion(x)`` is within ``tol_px`` and the runner-up
    catalogue line is more than ``unique_tol_px`` away, so that a blend of two catalogue
    lines never enters the fit. Returns a :data:`~mastpipe.core.model.LINE_COLUMNS` table
    (``x, wavelength, species, dx_px``) in the order of ``x_obs``.
    """
    x_obs = np.asarray(x_obs, dtype=float).ravel()
    cat_wl = np.asarray(cat_wl, dtype=float).ravel()
    cat_sp = np.asarray(cat_sp, dtype=object).ravel()
    if cat_wl.shape != cat_sp.shape:
        raise ValueError(f"cat_wl and cat_sp differ in length: {cat_wl.size} vs {cat_sp.size}")
    if x_obs.size == 0 or cat_wl.size == 0:
        return _empty_lines()
    lam = np.atleast_1d(np.asarray(model.wavelength(x_obs), dtype=float))
    disp = np.atleast_1d(np.asarray(model.dispersion(x_obs), dtype=float))
    if cat_wl.size == 1:
        take = np.zeros(lam.size, dtype=int)
        runner_up = np.full(lam.size, np.inf)
    else:
        j = np.clip(np.searchsorted(cat_wl, lam), 1, cat_wl.size - 1)
        take = np.where(np.abs(lam - cat_wl[j - 1]) < np.abs(lam - cat_wl[j]), j - 1, j)
        # a match is only trustworthy if the runner-up is far away
        distances = np.sort(np.abs(lam[:, None] - cat_wl[None, :]) / disp[:, None], axis=1)
        runner_up = distances[:, 1]
    d_px = (lam - cat_wl[take]) / disp
    ok = (np.abs(d_px) < tol_px) & (runner_up > unique_tol_px)
    return pd.DataFrame(
        {
            "x": x_obs[ok],
            "wavelength": cat_wl[take][ok],
            "species": cat_sp[take][ok],
            "dx_px": d_px[ok],
        }
    )


def shift_to_reference(flux, flux_ref, max_shift: int = 200) -> float:
    """Pixel shift of ``flux`` relative to ``flux_ref`` by cross-correlation.

    Both spectra are NaN-filled with zero and normalised to zero mean and unit scatter;
    the lag in ``[-max_shift, max_shift]`` that maximises their dot product is refined by a
    parabola through the peak and its neighbours. Positive when the features of ``flux``
    sit at larger columns than in ``flux_ref``, so the shift adds to the reference's
    ``x_c``.
    """
    a = np.nan_to_num(np.asarray(flux_ref, dtype=float).ravel())
    b = np.nan_to_num(np.asarray(flux, dtype=float).ravel())
    if a.size != b.size:
        raise ValueError(f"flux and flux_ref differ in length: {b.size} vs {a.size}")
    a = (a - a.mean()) / (a.std() or 1.0)
    b = (b - b.mean()) / (b.std() or 1.0)
    lags = np.arange(-int(max_shift), int(max_shift) + 1)
    cc = np.array([np.dot(np.roll(b, -int(lag)), a) for lag in lags], dtype=float)
    i = int(np.argmax(cc))
    peak = float(i)
    if 0 < i < cc.size - 1:
        y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
        curvature = y0 - 2.0 * y1 + y2
        if curvature != 0.0:
            peak = i + 0.5 * (y0 - y2) / curvature
    return float(np.interp(peak, np.arange(lags.size), lags))


# ------------------------------------------------------------------------- one fiber


def solve_fiber(
    x_obs,
    model0: WavelengthModel,
    cat_wl,
    cat_sp,
    params: WavecalParams,
    final_model: str = "grating",
) -> tuple[GratingModel | None, pd.DataFrame]:
    """identify -> fit -> clip from ``model0``, repeated until it stops moving.

    ``x_obs`` are the measured columns of one fiber's usable lines, ``model0`` the grating
    model that seeds the first identification (tolerance ``params.tol_first_px``; later
    passes use ``tol_later_px``). Each pass fits the matches with
    :func:`~mastpipe.steps.wavemodel.fit_model` (``params.fit_alpha``, ``params.poly_deg``;
    the distortion term starts at zero every pass), drops the lines whose residual is more
    than ``clip_sigma`` robust sigmas from the median residual (when at least ``min_lines``
    remain) and refits. A pass that moves ``lambda(x)`` by less than ``converge_px`` (in
    pixels of the central dispersion) anywhere on the detector and changes no line has
    converged; otherwise ``n_passes`` are run. Fewer than ``min_lines`` matches at any pass
    give up: ``(None, matches)``.

    Returns the fitted model and the identified lines, whose ``dx_px`` column holds the
    final fit residual ``x_model(wavelength) - x`` in pixels (the notebook's ``resid_px``).
    ``final_model`` is ``"grating"``; ``"poly"`` raises ``NotImplementedError``.
    """
    if final_model not in FINAL_MODELS:
        raise ValueError(f"final_model must be one of {FINAL_MODELS}, got {final_model!r}")
    if final_model == "poly":
        raise NotImplementedError("the polynomial final model is not part of this plan")
    if not isinstance(model0, GratingModel):
        raise NotImplementedError("solve_fiber identifies and fits on the grating model only")
    x_obs = np.asarray(x_obs, dtype=float).ravel()
    model = model0
    n_pix = model.n_pix
    xs = np.arange(n_pix, dtype=float)
    ident = _empty_lines()
    resid = np.zeros(0)
    prev_lam: np.ndarray | None = None
    prev_set: tuple[float, ...] | None = None
    for p in range(params.n_passes):
        tol = params.tol_first_px if p == 0 else params.tol_later_px
        ident = identify_lines(x_obs, model, cat_wl, cat_sp, tol, params.unique_tol_px)
        if len(ident) < params.min_lines:
            return None, ident
        new, resid = fit_model(
            model,
            ident["x"].to_numpy(),
            ident["wavelength"].to_numpy(),
            fit_alpha=params.fit_alpha,
            poly_deg=params.poly_deg,
        )
        sigma = robust_sigma(resid)
        keep = np.abs(resid - np.median(resid)) < params.clip_sigma * max(sigma, SIGMA_FLOOR_PX)
        if params.min_lines <= keep.sum() < len(ident):
            ident = ident.loc[keep].reset_index(drop=True)
            new, resid = fit_model(
                model,
                ident["x"].to_numpy(),
                ident["wavelength"].to_numpy(),
                fit_alpha=params.fit_alpha,
                poly_deg=params.poly_deg,
            )
        model = new
        lam_now = np.asarray(model.wavelength(xs), dtype=float)
        line_set = tuple(np.round(ident["wavelength"].to_numpy(dtype=float), LINE_SET_DECIMALS))
        if prev_lam is not None:
            moved = float(np.abs(lam_now - prev_lam).max() / model.dispersion(n_pix // 2))
            if moved < params.converge_px and line_set == prev_set:
                break
        prev_lam, prev_set = lam_now, line_set
    ident = ident.copy()
    ident["dx_px"] = resid
    return model, ident


# ------------------------------------------------------------------------ the full arc


def _usable_columns(table: pd.DataFrame) -> np.ndarray:
    if len(table) == 0 or "use_x" not in table.columns:
        return np.zeros(0)
    use = table["use_x"].to_numpy(dtype=bool)
    return table.loc[use, "x_center"].to_numpy(dtype=float)


def _warm_start(solution: WavelengthSolution, geom: GratingGeometry, n_pix: int) -> GratingModel:
    """``x_c, B, alpha`` of a stored solution as a seed; the distortion term is left to the fit."""
    if solution.n_pix != n_pix:
        raise ValueError(
            f"fiber {solution.fiber}: the previous solution has n_pix {solution.n_pix}, "
            f"the spectrum {n_pix}"
        )
    previous = model_from_solution(solution, geom)
    if not isinstance(previous, GratingModel):
        raise NotImplementedError("solve_full_arc seeds from grating solutions only")
    return GratingModel(previous.x_c, previous.B, previous.alpha, None, n_pix, geom)


def _seed_models(
    start: ReferenceSolution | WavelengthModel,
    fibers: list[int],
    flux_list: list[np.ndarray],
    traceset: TraceSet,
    ref_fiber: int,
    geom: GratingGeometry | None,
) -> dict[int, GratingModel]:
    """The starting model of every requested fiber; ``fibers`` and ``ref_fiber`` are trace
    indices (``Trace.index``), the lists are in trace order."""
    position = _positions(traceset)
    seeds: dict[int, GratingModel] = {}
    if isinstance(start, ReferenceSolution):
        if geom is None:
            raise ValueError("solve_full_arc needs geom to rebuild the models of a reference")
        if not start.solutions:
            raise ValueError("the previous reference holds no solutions to seed from")
        previous = start.solutions
        centres = {fiber: solution.trace_y_center for fiber, solution in previous.items()}
        for fiber in fibers:
            k = position[fiber]
            n_pix = len(flux_list[k])
            if fiber in previous:
                seeds[fiber] = _warm_start(previous[fiber], geom, n_pix)
                continue
            # absent from the previous reference: the nearest solved fiber, slid
            y = traceset.traces[k].y_center
            j = min(centres, key=lambda other: abs(centres[other] - y))
            seed = _warm_start(previous[j], geom, n_pix)
            shift = (
                shift_to_reference(flux_list[k], flux_list[position[j]]) if j in position else 0.0
            )
            seeds[fiber] = GratingModel(seed.x_c + shift, seed.B, seed.alpha, None, n_pix, geom)
        return seeds
    if not isinstance(start, GratingModel):
        raise NotImplementedError("solve_full_arc seeds from a grating model only")
    if ref_fiber not in position:
        raise ValueError(f"ref_fiber {ref_fiber} is not among the traces {sorted(position)}")
    reference_flux = flux_list[position[ref_fiber]]
    for fiber in fibers:
        k = position[fiber]
        n_pix = len(flux_list[k])
        shift = 0.0 if fiber == ref_fiber else shift_to_reference(flux_list[k], reference_flux)
        seeds[fiber] = GratingModel(
            start.x_c + shift, start.B, start.alpha, None, n_pix, start.geom
        )
    return seeds


def _positions(traceset: TraceSet) -> dict[int, int]:
    """Trace index -> position in ``traceset`` (``ValueError`` when an index repeats)."""
    position = {int(trace.index): k for k, trace in enumerate(traceset.traces)}
    if len(position) != len(traceset.traces):
        raise ValueError("solve_full_arc needs unique trace indices")
    return position


def solve_full_arc(
    line_tables: list[pd.DataFrame],
    flux_list: list[np.ndarray],
    traceset: TraceSet,
    start: ReferenceSolution | WavelengthModel,
    cat_wl,
    cat_sp,
    params: WavecalParams,
    fibers: list[int] | None = None,
    ref_fiber: int = 0,
    geom: GratingGeometry | None = None,
) -> dict[int, WavelengthSolution]:
    """Solve every fiber of a full arc independently with :func:`solve_fiber`.

    ``line_tables`` are the ``find_lines`` tables and ``flux_list`` the continuum-subtracted
    fluxes (``flux - continuum``), one per trace of ``traceset`` in trace order; a fiber's
    measured columns are the ``x_center`` of its ``use_x`` rows. Fibers are named by their
    trace index (``Trace.index``, the layout id after
    :func:`~mastpipe.steps.trace.layout_indexed`): the result is keyed by it, ``fibers``
    (default: every trace) restricts the solve to those indices and ``ref_fiber`` is the
    index of the anchored trace.

    Seeding: when ``start`` is a single grating model it belongs to fiber ``ref_fiber`` and
    every other fiber starts from it slid by :func:`shift_to_reference` of its flux against
    the reference fiber's. When ``start`` is a previous
    :class:`~mastpipe.core.model.ReferenceSolution`, each fiber warm-starts from its own
    earlier ``x_c, B, alpha`` (rebuilt with ``geom``, which is then required); a fiber the
    previous reference lacks starts from the solved fiber nearest in ``trace_y_center``,
    slid by :func:`shift_to_reference`.

    Fibers without usable lines or with fewer than ``params.min_lines`` matches are absent
    from the result. Each :class:`~mastpipe.core.model.WavelengthSolution` records the
    fitted parameters, the identified lines (``dx_px`` the fit residual), the rms residual,
    the range of the identified lines as ``anchor_px``/``anchor_nm``, the trace's
    ``y_center`` and the grating constants :data:`GEOMETRY_KEYS` for provenance.
    """
    n = len(traceset)
    if len(line_tables) != n or len(flux_list) != n:
        raise ValueError(
            f"line_tables ({len(line_tables)}) and flux_list ({len(flux_list)}) must have one "
            f"entry per trace ({n})"
        )
    position = _positions(traceset)
    fibers = sorted(position) if fibers is None else [int(k) for k in fibers]
    for fiber in fibers:
        if fiber not in position:
            raise ValueError(f"fiber {fiber} is not among the traces {sorted(position)}")
    seeds = _seed_models(start, fibers, flux_list, traceset, ref_fiber, geom)
    solutions: dict[int, WavelengthSolution] = {}
    for fiber in fibers:
        k = position[fiber]
        x_obs = _usable_columns(line_tables[k])
        if x_obs.size == 0:
            continue
        model, ident = solve_fiber(x_obs, seeds[fiber], cat_wl, cat_sp, params)
        if model is None:
            continue
        resid = ident["dx_px"].to_numpy(dtype=float)
        g = model.geom
        solutions[fiber] = WavelengthSolution(
            fiber=fiber,
            kind=model.kind,
            params=model.params,
            n_pix=model.n_pix,
            lines=ident[list(LINE_COLUMNS)],
            rms_px=float(np.sqrt(np.mean(resid**2))),
            n_lines=len(ident),
            anchor_px=(float(ident["x"].min()), float(ident["x"].max())),
            anchor_nm=(float(ident["wavelength"].min()), float(ident["wavelength"].max())),
            trace_y_center=float(traceset.traces[k].y_center),
            geometry={key: getattr(g, key) for key in GEOMETRY_KEYS},
        )
    return solutions


# ------------------------------------------------------------------- reference products


def solutions_hash(solutions: dict[int, WavelengthSolution]) -> str:
    """SHA-256 (:data:`HASH_LENGTH` hex) over every solution's fiber, parameters and line
    positions and wavelengths, rounded to :data:`HASH_DECIMALS` decimals, in fiber order."""
    digest = hashlib.sha256()
    for fiber in sorted(solutions):
        solution = solutions[fiber]
        digest.update(f"fiber {fiber}\n".encode())
        params = ",".join(f"{v:.{HASH_DECIMALS}f}" for v in np.asarray(solution.params, float))
        digest.update(f"params {params}\n".encode())
        for x, wl in solution.lines[["x", "wavelength"]].to_numpy(dtype=float):
            digest.update(f"line {x:.{HASH_DECIMALS}f} {wl:.{HASH_DECIMALS}f}\n".encode())
    return digest.hexdigest()[:HASH_LENGTH]


def build_reference(
    instrument: str,
    configuration: str,
    fiber_preset: str,
    solutions: dict[int, WavelengthSolution],
    created: str,
    source_checksums: tuple[str, ...],
) -> ReferenceSolution:
    """Package the solutions of one full arc as a :class:`~mastpipe.core.model.ReferenceSolution`.

    ``identifier`` comes from :func:`~mastpipe.core.model.make_identifier`, ``valid_from``
    is ``created`` (the arc's mid-time) and ``hash`` is :func:`solutions_hash`.
    """
    return ReferenceSolution(
        instrument=instrument,
        configuration=configuration,
        fiber_preset=fiber_preset,
        solutions=solutions,
        created=created,
        valid_from=created,
        source_checksums=tuple(source_checksums),
        identifier=make_identifier(instrument, configuration, fiber_preset, created),
        hash=solutions_hash(solutions),
    )


def fiber_agreement(solutions: dict[int, WavelengthSolution]) -> pd.DataFrame:
    """The same line, fiber to fiber: the notebook's external accuracy check.

    Each fiber is solved on its own, so for every catalogue wavelength identified in at
    least two fibers the spread of its residual ``dx_px`` across them, ``spread_px =
    max - min``, measures how far the independent solutions disagree, free of any fit
    statistic. Columns :data:`AGREEMENT_COLUMNS` (``wavelength, n_fibers, spread_px``),
    sorted by wavelength.
    """
    frames = [
        solution.lines[["wavelength", "dx_px"]].assign(fiber=fiber)
        for fiber, solution in solutions.items()
        if len(solution.lines)
    ]
    if not frames:
        return pd.DataFrame({column: pd.Series(dtype=float) for column in AGREEMENT_COLUMNS})
    table = pd.concat(frames, ignore_index=True)
    grouped = table.groupby("wavelength", sort=True)["dx_px"]
    out = pd.DataFrame(
        {
            "wavelength": grouped.size().index.to_numpy(dtype=float),
            "n_fibers": grouped.size().to_numpy(dtype=int),
            "spread_px": (grouped.max() - grouped.min()).to_numpy(dtype=float),
        }
    )
    return out.loc[out["n_fibers"] >= 2].reset_index(drop=True)


def wavelength_map(solution: WavelengthSolution, geom: GratingGeometry) -> np.ndarray:
    """Wavelength (nm) of every pixel ``0 .. n_pix - 1`` of a stored solution.

    The model is rebuilt with :func:`~mastpipe.steps.wavemodel.model_from_solution`, which
    needs the configuration's ``geom``.
    """
    model = model_from_solution(solution, geom)
    return np.asarray(model.wavelength(np.arange(solution.n_pix, dtype=float)), dtype=float)
