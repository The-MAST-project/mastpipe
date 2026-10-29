"""Generate the provisional 45-trace survey layout of HighSpec.

Placeholder geometry until the survey fiber head is measured: 45 traces at a pitch of
8.5 px from y = 12.0 (cropped coordinates), the five arc fibers at slots 0, 11, 22, 33 and
44, and the remaining 40 slots alternating science and sky per telescope T01..T20. Writes
``src/mastpipe/instruments/highspec/layouts/survey45.csv``; run it again after editing
the numbers below and commit the result.
"""

from __future__ import annotations

import csv
from pathlib import Path

N_TRACES = 45
Y_START = 12.0
PITCH = 8.5
SLOPE = 0.0020
ARC_SLOTS = (0, 11, 22, 33, 44)
OUTPUT = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "mastpipe"
    / "instruments"
    / "highspec"
    / "layouts"
    / "survey45.csv"
)


def layout_rows() -> list[tuple[str, ...]]:
    rows: list[tuple[str, ...]] = []
    telescope = 1
    next_role = "science"
    for index in range(N_TRACES):
        y_nominal = f"{Y_START + PITCH * index:.1f}"
        if index in ARC_SLOTS:
            rows.append((str(index), "arc", "", y_nominal, f"{SLOPE:.4f}", ""))
            continue
        rows.append((str(index), next_role, f"T{telescope:02d}", y_nominal, f"{SLOPE:.4f}", ""))
        if next_role == "science":
            next_role = "sky"
        else:
            next_role = "science"
            telescope += 1
    return rows


def main() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["trace", "role", "telescope", "y_nominal", "slope_nominal", "aperture"])
        writer.writerows(layout_rows())
    print(OUTPUT.relative_to(Path(__file__).resolve().parents[1]).as_posix())


if __name__ == "__main__":
    main()
