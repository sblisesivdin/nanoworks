#!/usr/bin/env python3
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Audit archived Si QE/GPAW RPA x-direction spectra without DFT dependencies.

Usage:
    python check_silicon_rpa.py silicon-qe-optics-reference.zip silicon-gpaw-rpa-reference.zip

Uses only original seven-column AllData tables inside the existing tutorial ZIPs.
Does not alter or normalize input spectra. All values are finite-window checks,
not sum rules, converged material predictions, or automated linewidth fits.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import math
from pathlib import Path
from zipfile import ZipFile


# expected_points distinguish staged and combined runs without assuming
# their parent directory names inside the reference archives.
CASES = (
    ("staged", "QE", "", 1001, 3.700815, 451.089106, 136.650181),
    ("staged", "GPAW", "NLFC", 157, 3.730863, 194.888208, 115.712327),
    ("staged", "GPAW", "LFC", 157, 3.730863, 146.409241, 105.878612),
    ("combined", "QE", "", 101, 2.600740, 581.993527, 314.481151),
    ("combined", "GPAW", "NLFC", 79, 2.690322, 256.268452, 262.542969),
    ("combined", "GPAW", "LFC", 79, 2.690322, 276.260615, 257.672204),
)


def parse_table(data: bytes) -> list[tuple[float, ...]]:
    """Read unmodified seven-column Nanoworks output; reject malformed data."""
    lines = io.StringIO(data.decode("utf-8-sig"))
    header = next(lines, "")
    if "Energy(eV)" not in header or "Eps_img(1)" not in header:
        raise ValueError("Missing expected optical table header")
    rows: list[tuple[float, ...]] = []
    for number, line in enumerate(lines, start=2):
        if not line.strip():
            continue
        parts = line.split()
        if len(parts) != 7:
            raise ValueError(f"Line {number}: expected seven columns")
        row = tuple(float(x) for x in parts)
        if not all(map(math.isfinite, row)):
            raise ValueError(f"Line {number}: non-finite optical value")
        if rows and row[0] <= rows[-1][0]:
            raise ValueError(f"Line {number}: non-increasing photon energy")
        rows.append(row)
    if not rows:
        raise ValueError("Empty optical data table")
    return rows


def read_case(z: ZipFile, engine: str, kind: str, points: int):
    suffix = f"-OPTICAL-{engine}-Result-Calculation-RPA-"
    if kind:
        suffix += kind + "-"
    suffix += "AllData_xdirection.dat"
    matches = []
    for item in z.infolist():
        if item.is_dir() or not item.filename.endswith(suffix):
            continue
        data = z.read(item)
        table = parse_table(data)
        if len(table) == points:
            matches.append((item.filename, data, table))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {engine} {kind or 'RPA'} table with "
            f"{points} rows; found {len(matches)} (suffix {suffix})"
        )
    return matches[0]


def eps2_at(rows, energy: float) -> float:
    if energy < rows[0][0] or energy > rows[-1][0]:
        raise ValueError("Integration bounds outside original energy grid")
    for prev, nxt in zip(rows, rows[1:]):
        if prev[0] <= energy <= nxt[0]:
            ratio = (energy - prev[0]) / (nxt[0] - prev[0])
            return prev[2] + ratio * (nxt[2] - prev[2])
    return rows[-1][2]


def spectral_area(rows, lower: float = 0.001, upper: float = 8.0) -> float:
    """Trapezoid on native grid; interpolate interval endpoints only."""
    samples = [(lower, eps2_at(rows, lower))]
    samples.extend((r[0], r[2]) for r in rows if lower < r[0] < upper)
    samples.append((upper, eps2_at(rows, upper)))
    return sum(
        (right[0] - left[0]) * (right[1] + left[1]) / 2
        for left, right in zip(samples, samples[1:])
    )


def check_close(actual: float, expected: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=0.001):
        raise ValueError(f"{label}: {actual:.9f} != archived {expected:.9f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qe_zip", type=Path, help="QE silicon optics reference ZIP")
    parser.add_argument("gpaw_zip", type=Path, help="GPAW RPA silicon reference ZIP")
    args = parser.parse_args()
    results = {}
    with ZipFile(args.qe_zip) as qe, ZipFile(args.gpaw_zip) as gpaw:
        print("# Original x-direction RPA data; integration range 0.001--8.0 eV")
        print("# case,engine,treatment,points,peak_energy_eV,peak_epsilon2_1,integral_epsilon2_dE_eV,source_sha256")
        for case, engine, kind, count, ex, ey, ea in CASES:
            filename, source, rows = read_case(
                qe if engine == "QE" else gpaw, engine, kind, count
            )
            peak = max(rows, key=lambda r: r[2])
            area = spectral_area(rows)
            check_close(peak[0], ex, f"{case} {engine} peak energy")
            check_close(peak[2], ey, f"{case} {engine} peak epsilon2")
            check_close(area, ea, f"{case} {engine} area")
            results[(case, engine, kind)] = (peak[2], area)
            print(
                f"{case},{engine},{kind or 'epsilon.x'},{count},"
                f"{peak[0]:.6f},{peak[2]:.6f},{area:.6f},"
                f"{hashlib.sha256(source).hexdigest()}"
            )
    print("# peak_ratio_QE_to_GPAW_NLFC,area_ratio_QE_to_GPAW_NLFC")
    for case in ("staged", "combined"):
        qe_peak, qe_area = results[(case, "QE", "")]
        gp_peak, gp_area = results[(case, "GPAW", "NLFC")]
        print(f"# {case}: {qe_peak / gp_peak:.6f}, {qe_area / gp_area:.6f}")


if __name__ == "__main__":
    main()
