#!/usr/bin/env python3
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Compare new QE epsilon.x and GPAW *NLFC* Si RPA tables without scaling.

Example:
  python compare_silicon_rpa_runs.py QE-AllData_xdirection.dat GPAW-NLFC-AllData_xdirection.dat

Supply the actual structure-prefixed filenames in your results directories.
The defaults integrate the common 0.001--8 eV range and use the input
pair's nominal 0.10 eV isolated-line FWHM. No DFT library is required.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path


def read_table(path: Path) -> list[tuple[float, float]]:
    """Read the energy and imaginary-dielectric columns of a Nanoworks DAT."""
    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8-sig") as file:
        header = file.readline()
        if "Energy(eV)" not in header or "Eps_img(1)" not in header:
            raise ValueError(f"{path}: unsupported optical DAT header")
        for number, line in enumerate(file, start=2):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 7:
                raise ValueError(f"{path}:{number}: expected 7 optical columns")
            values = tuple(float(item) for item in fields)
            if not all(math.isfinite(item) for item in values):
                raise ValueError(f"{path}:{number}: non-finite optical value")
            energy, eps2 = values[0], values[2]
            if rows and energy <= rows[-1][0]:
                raise ValueError(f"{path}:{number}: energy grid not increasing")
            rows.append((energy, eps2))
    if len(rows) < 3:
        raise ValueError(f"{path}: expected at least three photon-energy samples")
    return rows


def value_at(rows: list[tuple[float, float]], energy: float) -> float:
    """Linear interpolation at endpoints only; never extrapolate."""
    if energy < rows[0][0] - 1e-10 or energy > rows[-1][0] + 1e-10:
        raise ValueError("Requested integration window exceeds a data grid")
    for (left_e, left_y), (right_e, right_y) in zip(rows, rows[1:]):
        if left_e - 1e-10 <= energy <= right_e + 1e-10:
            fraction = (energy - left_e) / (right_e - left_e)
            return left_y + (right_y - left_y) * fraction
    return rows[-1][1]


def moments(rows: list[tuple[float, float]], lower: float, upper: float) -> tuple[float, float]:
    """Return integral epsilon2 dE [eV], integral E*epsilon2 dE [eV^2]."""
    if not 0 <= lower < upper:
        raise ValueError("Integration interval must obey 0 <= lower < upper")
    included = [(lower, value_at(rows, lower))]
    included.extend(pair for pair in rows if lower < pair[0] < upper)
    included.append((upper, value_at(rows, upper)))
    area = weighted_area = 0.0
    for (e0, y0), (e1, y1) in zip(included, included[1:]):
        delta = e1 - e0
        area += 0.5 * delta * (y0 + y1)
        weighted_area += 0.5 * delta * (e0 * y0 + e1 * y1)
    return area, weighted_area


def inspect(path: Path, lower: float, upper: float, width: float) -> dict:
    rows = read_table(path)
    peak_index = max(range(len(rows)), key=lambda index: rows[index][1])
    if peak_index in (0, len(rows) - 1):
        raise ValueError(f"{path}: peak at grid boundary, extend frequency range")
    spacing = max(
        rows[peak_index][0] - rows[peak_index - 1][0],
        rows[peak_index + 1][0] - rows[peak_index][0],
    )
    area, weighted = moments(rows, lower, upper)
    return {
        "points": len(rows),
        "peak_eV": rows[peak_index][0],
        "peak_epsilon2": rows[peak_index][1],
        "integral_epsilon2_eV": area,
        "integral_E_epsilon2_eV2": weighted,
        "peak_step_eV": spacing,
        "nominal_FWHM_to_step": width / spacing,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qe_dat", type=Path, help="QE RPA AllData_xdirection.dat")
    parser.add_argument("gpaw_nlfc_dat", type=Path, help="GPAW RPA NLFC AllData_xdirection.dat")
    parser.add_argument("--emin", type=float, default=0.001)
    parser.add_argument("--emax", type=float, default=8.0)
    parser.add_argument("--fwhm", type=float, default=0.10,
                        help="Nominal isolated-line target FWHM in eV")
    args = parser.parse_args()
    if not math.isfinite(args.fwhm) or args.fwhm <= 0.0:
        parser.error("--fwhm must be a positive finite number")
    for value in (args.emin, args.emax):
        if not math.isfinite(value):
            parser.error("Integration bounds must be finite")

    stats = {
        "QE": inspect(args.qe_dat, args.emin, args.emax, args.fwhm),
        "GPAW_NLFC": inspect(args.gpaw_nlfc_dat, args.emin, args.emax, args.fwhm),
    }
    print("# Observed peaks and native-grid trapezoid integrals; no amplitude scaling")
    print("# Nominal FWHM is an isolated-line *estimate*, not a fitted width")
    print("engine,points,peak_energy_eV,peak_epsilon2_1,integral_epsilon2_dE_eV,"
          "integral_E_epsilon2_dE_eV2,local_step_eV,nominal_FWHM_over_step,source_sha256")
    for engine, data in stats.items():
        print(
            f"{engine},{data['points']},{data['peak_eV']:.6f},"
            f"{data['peak_epsilon2']:.6f},{data['integral_epsilon2_eV']:.6f},"
            f"{data['integral_E_epsilon2_eV2']:.6f},"
            f"{data['peak_step_eV']:.6f},{data['nominal_FWHM_to_step']:.3f},"
            f"{data['sha256']}"
        )
    qe, gpaw = stats["QE"], stats["GPAW_NLFC"]
    print(f"# QE/GPAW_NLFC peak ratio: {qe['peak_epsilon2'] / gpaw['peak_epsilon2']:.6f}"
          if gpaw["peak_epsilon2"] else "# QE/GPAW_NLFC peak ratio: undefined")
    print(
        f"# QE/GPAW_NLFC [{args.emin:.6f}, {args.emax:.6f}] eV area ratio: "
        f"{qe['integral_epsilon2_eV'] / gpaw['integral_epsilon2_eV']:.6f}"
        if gpaw["integral_epsilon2_eV"] else
        "# QE/GPAW_NLFC area ratio: undefined"
    )


if __name__ == "__main__":
    main()
