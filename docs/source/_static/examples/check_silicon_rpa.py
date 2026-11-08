#!/usr/bin/env python3
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Audit archived Si QE/GPAW RPA x-direction spectra without DFT dependencies.

Usage:
    python check_silicon_rpa.py silicon-qe-optics-reference.zip silicon-gpaw-rpa-reference.zip

Uses only original seven-column AllData tables inside the existing tutorial ZIPs.
Does not alter or normalize input spectra. CSV rows are written to stdout;
notes, ratios, and warnings to stderr. All values are finite-window checks,
not sum rules, converged material predictions, or automated linewidth fits.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import math
from pathlib import Path, PurePosixPath
import sys
from zipfile import ZipFile


# These records refer to the actual staged/ and combined/ reference ZIP
# layouts. QE tables are named optical-x.dat in its archive; GPAW tables
# retain the full structure-prefixed AllData_xdirection.dat names.
# The row counts and numerical checks validate each independently.
CASES = (
    ("staged", "QE", "", 1001, 3.700815, 451.089106, 136.650181, 475.153835, 0.019999),
    ("staged", "GPAW", "NLFC", 157, 3.730863, 194.888208, 115.712327, 398.876659, 0.136864),
    ("staged", "GPAW", "LFC", 157, 3.730863, 146.409241, 105.878612, 376.085037, 0.136864),
    ("combined", "QE", "", 101, 2.600740, 581.993527, 314.481151, 979.918073, 0.099990),
    ("combined", "GPAW", "NLFC", 79, 2.690322, 256.268452, 262.542969, 818.145815, 0.271871),
    ("combined", "GPAW", "LFC", 79, 2.690322, 276.260615, 257.672204, 814.208467, 0.271871),
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


def read_case(z: ZipFile, case: str, engine: str, kind: str, points: int):
    """Match the *published* ZIP structure, then independently check row count."""
    if engine == "QE":
        basename = "optical-x.dat"
    else:
        basename = (
            f"Si_mp-149_primitive_Example-OPTICAL-GPAW-"
            f"Result-Calculation-RPA-{kind}-AllData_xdirection.dat"
        )
    matches = []
    for item in z.infolist():
        if item.is_dir():
            continue
        parts = PurePosixPath(item.filename).parts
        if case not in parts or parts[-1] != basename:
            continue
        data = z.read(item)
        table = parse_table(data)
        matches.append((item.filename, data, table))
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {case} {engine} {kind or 'epsilon.x'} "
            f"table ({basename}); found {len(matches)}"
        )
    filename, data, table = matches[0]
    if len(table) != points:
        raise ValueError(
            f"{filename}: expected {points} rows, got {len(table)}"
        )
    return filename, data, table


def eps2_at(rows, energy: float) -> float:
    if energy < rows[0][0] or energy > rows[-1][0]:
        raise ValueError("Integration bounds outside original energy grid")
    for prev, nxt in zip(rows, rows[1:]):
        if prev[0] <= energy <= nxt[0]:
            ratio = (energy - prev[0]) / (nxt[0] - prev[0])
            return prev[2] + ratio * (nxt[2] - prev[2])
    return rows[-1][2]


def window_samples(rows, lower: float, upper: float):
    """Clip a native grid to a finite integration window, without regridding."""
    if not (math.isfinite(lower) and math.isfinite(upper)
            and 0 <= lower < upper):
        raise ValueError("Integration interval must obey 0 <= lower < upper")
    samples = [(lower, eps2_at(rows, lower))]
    samples.extend((r[0], r[2]) for r in rows if lower < r[0] < upper)
    samples.append((upper, eps2_at(rows, upper)))
    return samples


def spectral_moments(rows, lower: float = 0.001, upper: float = 8.0):
    """Exact zeroth/first moments of piecewise-linear epsilon2 samples."""
    samples = window_samples(rows, lower, upper)
    zeroth = first = 0.0
    for (e0, y0), (e1, y1) in zip(samples, samples[1:]):
        delta = e1 - e0
        zeroth += 0.5 * delta * (y0 + y1)
        first += delta * (
            (2.0 * e0 + e1) * y0 + (e0 + 2.0 * e1) * y1
        ) / 6.0
    return zeroth, first


def spectral_area(rows, lower: float = 0.001, upper: float = 8.0) -> float:
    """Trapezoidal area of the original epsilon2 samples."""
    return spectral_moments(rows, lower, upper)[0]


def window_max_step(rows, lower: float = 0.001, upper: float = 8.0):
    """Largest frequency step across the complete selected window."""
    samples = window_samples(rows, lower, upper)
    return max(b[0] - a[0] for a, b in zip(samples, samples[1:]))


def peak_spacing(rows, peak_index: int) -> float:
    """Largest immediate energy gap adjacent to the sampled peak."""
    if not 0 < peak_index < len(rows) - 1:
        raise ValueError("Peak at the boundary: sampling cannot be assessed")
    return max(
        rows[peak_index][0] - rows[peak_index - 1][0],
        rows[peak_index + 1][0] - rows[peak_index][0],
    )


def nominal_fwhm(case: str, engine: str) -> float:
    """Isolated-transition *estimate*, not fitted width of the Si spectrum."""
    input_eta = 0.05 if case == "staged" else 0.10
    return input_eta * (2.0 if engine == "GPAW" else 1.0)


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
        print("Original x-direction RPA data; integration range 0.001--8.0 eV.", file=sys.stderr)
        print("Nominal isolated-line FWHM is only a sampling diagnostic.", file=sys.stderr)
        print("case,engine,treatment,points,peak_energy_eV,peak_epsilon2_1,integral_epsilon2_dE_eV,integral_E_epsilon2_dE_eV2,adjacent_step_eV,max_window_step_eV,nominal_FWHM_eV,min_FWHM_over_window_step,source_sha256")
        for case, engine, kind, count, ex, ey, ea, em, es in CASES:
            filename, source, rows = read_case(
                qe if engine == "QE" else gpaw, case, engine, kind, count
            )
            peak_index = max(range(len(rows)), key=lambda i: rows[i][2])
            peak = rows[peak_index]
            area, first_moment = spectral_moments(rows)
            step = peak_spacing(rows, peak_index)
            max_step = window_max_step(rows)
            width = nominal_fwhm(case, engine)
            check_close(peak[0], ex, f"{case} {engine} peak energy")
            check_close(peak[2], ey, f"{case} {engine} peak epsilon2")
            check_close(area, ea, f"{case} {engine} area")
            check_close(first_moment, em, f"{case} {engine} first moment")
            check_close(max_step, es, f"{case} {engine} max grid step")
            results[(case, engine, kind)] = (peak[2], area, first_moment)
            print(
                f"{case},{engine},{kind or 'epsilon.x'},{count},"
                f"{peak[0]:.6f},{peak[2]:.6f},{area:.6f},"
                f"{first_moment:.6f},{step:.6f},{max_step:.6f},"
                f"{width:.6f},{width / max_step:.3f},"
                f"{hashlib.sha256(source).hexdigest()}"
            )
    print("QE/GPAW_NLFC peak and finite-window area ratios:", file=sys.stderr)
    for case in ("staged", "combined"):
        qe_peak, qe_area, qe_moment = results[(case, "QE", "")]
        gp_peak, gp_area, gp_moment = results[(case, "GPAW", "NLFC")]
        print(
            f"{case}: peak={qe_peak / gp_peak:.6f}, "
            f"area={qe_area / gp_area:.6f}, "
            f"first_moment={qe_moment / gp_moment:.6f}",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
