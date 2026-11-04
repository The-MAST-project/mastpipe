"""Produce the figures embedded in the documentation from a simulated night.

Run from the repository root; writes PNG files under ``docs/_static/``. The figures are
committed so that the documentation build stays offline and light. Nothing here reads real
data.

Run: python scripts/make_doc_figures.py
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from astropy.io import fits  # noqa: E402

from mastpipe.config.loader import build_effective  # noqa: E402
from mastpipe.instruments.registry import load_instrument, load_layout, load_line_table  # noqa: E402
from mastpipe.pipeline.identify import run_identify  # noqa: E402
from mastpipe.sim.simulator import simulate_night  # noqa: E402

STATIC = Path(__file__).resolve().parent.parent / "docs" / "_static"


def main() -> int:
    STATIC.mkdir(parents=True, exist_ok=True)
    instrument = load_instrument("highspec")
    effective = build_effective(instrument, "Mg", "main", layout="survey45")
    layout = load_layout("highspec", "survey45")
    lines = load_line_table("highspec", "Mg")
    with tempfile.TemporaryDirectory() as tmp:
        night = Path(tmp) / "night"
        truth = simulate_night(effective, layout, lines, night, seed=0, n_science=2)
        frames = truth.frames
        full_arc = night / frames.loc[frames["arc_kind"] == "full", "file"].iloc[0]
        science = night / frames.loc[frames["exposure_type"] == "science", "file"].iloc[0]

        # 1. the anchor panels, as identify draws them without anchors
        out = Path(tmp) / "out"
        _, png = run_identify(
            full_arc, effective, "survey45", out, anchors=None, write_reference=False
        )
        (STATIC / "anchor_panels.png").write_bytes(Path(png).read_bytes())

        # 2. a simulated science frame: 45 traces, dispersion along x
        data = fits.getdata(science).astype(float)
        fig, ax = plt.subplots(figsize=(9, 3.2), dpi=120)
        lo, hi = np.percentile(data, [5, 99.5])
        ax.imshow(data, origin="lower", cmap="gray", vmin=lo, vmax=hi, aspect="auto")
        ax.set_xlabel("column (dispersion)")
        ax.set_ylabel("row (fibers)")
        ax.set_title(
            "Simulated HighSpec science frame, survey layout: 20 science, 20 sky, 5 dark arc traces"
        )
        fig.tight_layout()
        fig.savefig(STATIC / "simulated_frame.png")
        plt.close(fig)
    print(f"wrote {STATIC / 'anchor_panels.png'} and {STATIC / 'simulated_frame.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
