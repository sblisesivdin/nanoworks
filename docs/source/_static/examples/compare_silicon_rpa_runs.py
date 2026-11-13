#!/usr/bin/env python3
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Compare new QE epsilon.x and GPAW *NLFC* Si RPA tables without scaling.

Example:
  python compare_silicon_rpa_runs.py QE-AllData_xdirection.dat GPAW-NLFC-AllData_xdirection.dat

Supply the actual structure-prefixed filenames in your results directories.
The defaults integrate the common 0.001--8 eV range and use the input
pair's nominal 0.10 eV isolated-line FWHM. CSV rows go to stdout;
notes, ratios, and warnings go to stderr. No DFT library is required.
"""

from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path
import sys


OPTICAL_COLUMNS = (
    "Energy(eV)", "Eps_real(1)", "Eps_img(1)", "Refractive_Index(1)",
    "Extinction_Index(1)", "Absorption(1/cm)", "Reflectivity(1)",
)


def read_table(path: Path) -> list[tuple[float, float]]:
    """Read the energy and imaginary-dielectric columns of a Nanoworks DAT."""
    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8-sig") as file:
        header = file.readline()
        if tuple(header.split()) != OPTICAL_COLUMNS:
            raise ValueError(
                f"{path}: unexpected optical DAT header or column order; "
                "expected the seven named Nanoworks optical columns"
            )
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
            if energy < 0.0:
                raise ValueError(f"{path}:{number}: negative photon energy")
            if rows and energy <= rows[-1][0]:
                raise ValueError(f"{path}:{number}: energy grid not increasing")
            rows.append((energy, eps2))
    if len(rows) < 3:
        raise ValueError(f"{path}: expected at least three photon-energy samples")
    return rows


def value_at(rows: list[tuple[float, float]], energy: float) -> float:
    """Linearly interpolate within the provided energy grid; never extrapolate."""
    if not math.isfinite(energy) or energy < rows[0][0] or energy > rows[-1][0]:
        raise ValueError("Requested integration window exceeds a data grid")
    for (left_e, left_y), (right_e, right_y) in zip(rows, rows[1:]):
        if left_e <= energy <= right_e:
            fraction = (energy - left_e) / (right_e - left_e)
            return left_y + (right_y - left_y) * fraction
    return rows[-1][1]


def samples_in_window(
    rows: list[tuple[float, float]], lower: float, upper: float
) -> list[tuple[float, float]]:
    """Clip the original grid to a covered window without regridding."""
    if not (math.isfinite(lower) and math.isfinite(upper) and
            0.0 <= lower < upper):
        raise ValueError("Integration interval must obey 0 <= lower < upper")
    samples = [(lower, value_at(rows, lower))]
    samples.extend(pair for pair in rows if lower < pair[0] < upper)
    samples.append((upper, value_at(rows, upper)))
    return samples


def moments(rows: list[tuple[float, float]], lower: float, upper: float) -> tuple[float, float]:
    """Integrate a piecewise-linear epsilon2 spectrum exactly.

    epsilon2 dE is trapezoidal; E*epsilon2 dE uses the exact primitive
    of E times each linearly interpolated dielectric segment.
    Neither estimate recovers unresolved peaks between samples.
    """
    included = samples_in_window(rows, lower, upper)
    area = weighted_area = 0.0
    for (e0, y0), (e1, y1) in zip(included, included[1:]):
        delta = e1 - e0
        area += 0.5 * delta * (y0 + y1)
        weighted_area += delta * (
            (2.0 * e0 + e1) * y0 + (e0 + 2.0 * e1) * y1
        ) / 6.0
    return area, weighted_area


def window_max_step(rows, lower: float, upper: float) -> float:
    """Largest adjacent interval on the original grid inside the window."""
    selected = samples_in_window(rows, lower, upper)
    return max(b[0] - a[0] for a, b in zip(selected, selected[1:]))


def inspect(path: Path, lower: float, upper: float, width: float) -> dict:
    rows = read_table(path)
    # Report the peak in the SAME window used for the finite-window moments.
    # Otherwise a large high-energy maximum outside --emax may be mistaken
    # for a peak belonging to the compared energy interval.
    selected = samples_in_window(rows, lower, upper)
    peak_index = max(range(len(selected)), key=lambda index: selected[index][1])
    if peak_index in (0, len(selected) - 1):
        raise ValueError(
            f"{path}: peak at selected window boundary; widen --emin/--emax"
        )
    spacing = max(
        selected[peak_index][0] - selected[peak_index - 1][0],
        selected[peak_index + 1][0] - selected[peak_index][0],
    )
    area, weighted = moments(rows, lower, upper)
    max_step = window_max_step(rows, lower, upper)
    return {
        "points": len(rows),
        "peak_eV": selected[peak_index][0],
        "peak_epsilon2": selected[peak_index][1],
        "integral_epsilon2_eV": area,
        "integral_E_epsilon2_eV2": weighted,
        "peak_step_eV": spacing,
        "nominal_FWHM_to_step": width / spacing,
        "max_step_in_window_eV": max_step,
        "min_nominal_FWHM_to_step": width / max_step,
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
    parser.add_argument("--min-grid-ratio", type=float, default=5.0,
                        help="Heuristic minimum FWHM / grid step throughout the "
                             "analysis window (default: 5; not convergence proof)")
    args = parser.parse_args()
    if not math.isfinite(args.fwhm) or args.fwhm <= 0.0:
        parser.error("--fwhm must be a positive finite number")
    if not math.isfinite(args.min_grid_ratio) or args.min_grid_ratio <= 0.0:
        parser.error("--min-grid-ratio must be a positive finite number")
    for value in (args.emin, args.emax):
        if not math.isfinite(value):
            parser.error("Integration bounds must be finite")

    stats = {
        "QE": inspect(args.qe_dat, args.emin, args.emax, args.fwhm),
        "GPAW_NLFC": inspect(args.gpaw_nlfc_dat, args.emin, args.emax, args.fwhm),
    }
    print("Peaks and native-grid integrals use the same selected energy window; no amplitude scaling.", file=sys.stderr)
    print("Nominal FWHM is an isolated-line estimate, not a fitted width.", file=sys.stderr)
    print("E*epsilon2 integrates the piecewise-linear response exactly.", file=sys.stderr)
    print("Grid warnings are heuristics, not convergence tests.", file=sys.stderr)
    print("engine,points,peak_energy_eV,peak_epsilon2_1,integral_epsilon2_dE_eV,"
          "integral_E_epsilon2_dE_eV2,local_step_eV,nominal_FWHM_over_step,"
          "max_window_step_eV,min_FWHM_over_window_step,resolution_warning,source_sha256")
    for engine, data in stats.items():
        sparse = data["min_nominal_FWHM_to_step"] < args.min_grid_ratio
        print(
            f"{engine},{data['points']},{data['peak_eV']:.6f},"
            f"{data['peak_epsilon2']:.6f},{data['integral_epsilon2_eV']:.6f},"
            f"{data['integral_E_epsilon2_eV2']:.6f},"
            f"{data['peak_step_eV']:.6f},{data['nominal_FWHM_to_step']:.3f},"
            f"{data['max_step_in_window_eV']:.6f},"
            f"{data['min_nominal_FWHM_to_step']:.3f},"
            f"{'YES' if sparse else 'NO'},{data['sha256']}"
        )
        if sparse:
            print(
                f"# WARNING: {engine} is below {args.min_grid_ratio:g} "
                "nominal-FWHM intervals per largest step in the window.",
                file=sys.stderr,
            )
    qe, gpaw = stats["QE"], stats["GPAW_NLFC"]
    print(
        f"QE/GPAW_NLFC peak ratio: "
        + (f"{qe['peak_epsilon2'] / gpaw['peak_epsilon2']:.6f}"
           if gpaw["peak_epsilon2"] else "undefined"),
        file=sys.stderr,
    )
    print(
        f"QE/GPAW_NLFC [{args.emin:.6f}, {args.emax:.6f}] eV area ratio: "
        + (f"{qe['integral_epsilon2_eV'] / gpaw['integral_epsilon2_eV']:.6f}"
           if gpaw["integral_epsilon2_eV"] else "undefined"),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
