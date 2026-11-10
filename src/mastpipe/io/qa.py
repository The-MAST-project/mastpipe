"""Quality-assurance figures. The only module that draws; the steps never do.

This plan's one figure is the anchor panels of ``mastpipe identify``, a port of the
notebook's ``plot_anchor_panels``: the catalogue atlas convolved to the instrument's
resolution above, the measured continuum-subtracted arc spectrum of one trace below, both
normalised to their peak and numbered on their strongest well-separated lines, so that a
few template-measured pairs can be read off by pattern and typed back as anchors.

Figures are drawn with matplotlib's object-oriented interface on an ``Agg`` canvas and
written to file (atomically, through :func:`~mastpipe.io.atomic.atomic_write`), so the
command runs unattended and on machines without a display, never touches the global
pyplot state, and leaves a notebook's own backend alone.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from mastpipe.io.atomic import atomic_write

FIGURE_SIZE = (13.0, 7.2)
"""Inches, the notebook's ``figwidth`` and height."""
DPI = 150
TEMPLATE_COLOR = "black"
TEMPLATE_LABEL_COLOR = "tab:blue"
MEASURED_COLOR = "tab:red"


def strongest_spread(positions, strength, n: int, x0: float, x1: float, min_sep: float) -> list[int]:
    """Indices of the ``n`` strongest entries within ``[x0, x1]``, no two closer than
    ``min_sep`` of the window's width, sorted by position (numbered left to right).

    Entries are visited in decreasing ``strength``; one is kept when it lies inside the
    window and is at least ``min_sep (x1 - x0)`` from every entry already kept.
    """
    positions = np.asarray(positions, dtype=float).ravel()
    strength = np.asarray(strength, dtype=float).ravel()
    if positions.shape != strength.shape:
        raise ValueError(f"positions and strength differ in length: {positions.size} vs {strength.size}")
    kept: list[int] = []
    gap = min_sep * (x1 - x0)
    for k in np.argsort(strength, kind="stable")[::-1]:
        p = positions[k]
        if not x0 <= p <= x1 or any(abs(p - positions[j]) < gap for j in kept):
            continue
        kept.append(int(k))
        if len(kept) >= n:
            break
    return sorted(kept, key=lambda k: positions[k])


def _label_peaks(ax, positions, heights, ids, color: str) -> None:
    for p, h, i in zip(positions, heights, ids, strict=True):
        ax.annotate(
            str(i),
            xy=(p, h),
            xytext=(0, 6),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8,
            color="white",
            weight="bold",
            bbox={"boxstyle": "round,pad=0.18", "fc": color, "ec": "none"},
        )


def _normalised(flux) -> np.ndarray:
    """``flux`` over its largest finite value (unchanged when there is none)."""
    flux = np.asarray(flux, dtype=float)
    finite = flux[np.isfinite(flux)]
    peak = float(finite.max()) if finite.size else 0.0
    return flux / peak if peak > 0 else flux


def _interp(x, grid, values) -> list[float]:
    good = np.isfinite(values)
    return [float(np.interp(v, grid[good], values[good])) for v in x] if good.any() else []


def plot_anchor_panels(
    atlas_wave,
    atlas_flux_conv,
    usable: pd.DataFrame,
    obs_flux,
    obs_lines: pd.DataFrame,
    band: tuple[float, float],
    n_labels: int = 8,
    min_sep: float = 0.04,
    flip_x: bool = False,
    *,
    out: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Template above, measurement below, at the same scale, both numbered; written to ``out``.

    The template is ``atlas_flux_conv`` on ``atlas_wave`` over ``band`` (air nm), the
    atlas convolved to the instrument's resolution; ``usable`` lists the catalogue lines
    (``wavelength``, bool ``usable``) and the ``n_labels`` usable lines with the highest
    template flux at their wavelength, no two closer than ``min_sep`` of the band, are
    numbered left to right. The measurement is ``obs_flux`` (one trace's
    continuum-subtracted spectrum, NaN where masked) against pixel, and its ``n_labels``
    tallest ``use_x`` lines of ``obs_lines`` (``x_center``, ``height``) are numbered the
    same way; ``flip_x`` mirrors that panel for a configuration whose wavelength falls
    with pixel. Match by pattern, not by height: the relative intensities of the atlas
    lamp and the instrument's differ.

    Returns ``(template, measured)``: ``id, wavelength`` and ``id, x_center`` tables, ids
    from 1 in each. The PNG is written to a temporary next to ``out`` and moved into place
    (:func:`~mastpipe.io.atomic.atomic_write`), so a failed drawing leaves no partial file;
    the parent folder is created.
    """
    out = Path(out)
    atlas_wave = np.asarray(atlas_wave, dtype=float).ravel()
    atlas_flux_conv = np.asarray(atlas_flux_conv, dtype=float).ravel()
    obs_flux = np.asarray(obs_flux, dtype=float).ravel()
    n_pix = obs_flux.size
    fig = Figure(figsize=FIGURE_SIZE, layout="constrained")
    FigureCanvasAgg(fig)
    axes = fig.subplots(2, 1)

    ax = axes[0]
    inside = (atlas_wave >= band[0]) & (atlas_wave <= band[1])
    tw, tf = atlas_wave[inside], _normalised(atlas_flux_conv[inside])
    ax.plot(tw, tf, color=TEMPLATE_COLOR, lw=0.8)
    candidates = usable.loc[usable["usable"].to_numpy(dtype=bool)].reset_index(drop=True)
    waves = candidates["wavelength"].to_numpy(dtype=float)
    strength = _interp(waves, tw, tf) if tw.size else [0.0] * waves.size
    pick = strongest_spread(waves, strength, n_labels, band[0], band[1], min_sep)
    template = pd.DataFrame({"id": np.arange(1, len(pick) + 1), "wavelength": waves[pick]})
    _label_peaks(
        ax,
        template["wavelength"],
        [strength[k] for k in pick],
        template["id"],
        TEMPLATE_LABEL_COLOR,
    )
    ax.set_xlim(*band)
    ax.set_ylim(-0.05, 1.25)
    ax.set_xlabel("wavelength [nm]")
    ax.set_ylabel("normalised flux")
    ax.set_title("TEMPLATE - catalogue atlas at the instrument resolution", loc="left")
    ax.grid(alpha=0.25)

    ax = axes[1]
    x_pix = np.arange(n_pix, dtype=float)
    of = _normalised(obs_flux)
    ax.plot(x_pix, of, color=MEASURED_COLOR, lw=0.8)
    lines = obs_lines.loc[obs_lines["use_x"].to_numpy(dtype=bool)].reset_index(drop=True)
    x_lines = lines["x_center"].to_numpy(dtype=float)
    pick = strongest_spread(x_lines, lines["height"].to_numpy(dtype=float), n_labels, 0.0, n_pix - 1.0, min_sep)
    measured = pd.DataFrame({"id": np.arange(1, len(pick) + 1), "x_center": x_lines[pick]})
    _label_peaks(
        ax,
        measured["x_center"],
        _interp(measured["x_center"], x_pix, of),
        measured["id"],
        MEASURED_COLOR,
    )
    ax.set_xlim((n_pix - 1.0, 0.0) if flip_x else (0.0, n_pix - 1.0))
    ax.set_ylim(-0.05, 1.25)
    ax.set_xlabel("pixel" + (" (mirrored)" if flip_x else ""))
    ax.set_ylabel("normalised flux")
    ax.set_title("MEASURED - continuum-subtracted arc spectrum of the anchor trace", loc="left")
    ax.grid(alpha=0.25)

    with atomic_write(out) as tmp:
        fig.savefig(tmp, dpi=DPI, format="png")  # the temporary's name has no .png suffix
    return template, measured
