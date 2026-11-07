# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Offline tests for the downloadable silicon cross-engine optics audit."""

import importlib.util
import io
from pathlib import Path
import unittest
from zipfile import ZipFile


_SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "docs/source/_static/examples/check_silicon_rpa.py"
)
_spec = importlib.util.spec_from_file_location("silicon_rpa_audit", _SCRIPT)
assert _spec is not None and _spec.loader is not None
_audit = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_audit)


def _dat(rows=None):
    if rows is None:
        rows = ((0.001, 0.0), (4.0, 2.0), (8.0, 0.0))
    lines = [
        "Energy(eV) Eps_real(1) Eps_img(1) "
        "Refractive_Index(1) Extinction_Index(1) Absorption(1/cm) Reflectivity(1)"
    ]
    lines += [
        f"{energy} 1.0 {imag} 1.0 0.0 0.0 0.0"
        for energy, imag in rows
    ]
    return ("\n".join(lines) + "\n").encode("utf-8")


def _zip(entries):
    stream = io.BytesIO()
    with ZipFile(stream, "w") as archive:
        for name, data in entries:
            archive.writestr(name, data)
    stream.seek(0)
    return ZipFile(stream)


class SiliconRPAAuditTests(unittest.TestCase):
    def test_original_qe_zip_names_are_selected_by_case(self):
        with _zip([
            ("staged/optical-x.dat", _dat()),
            ("combined/optical-x.dat", _dat()),
            ("staged/optical-y.dat", _dat()),
        ]) as archive:
            name, _, rows = _audit.read_case(
                archive, "staged", "QE", "", 3
            )
            self.assertEqual(name, "staged/optical-x.dat")
            self.assertEqual(len(rows), 3)
            other, _, _ = _audit.read_case(
                archive, "combined", "QE", "", 3
            )
            self.assertEqual(other, "combined/optical-x.dat")

    def test_original_gpaw_zip_name_is_selected(self):
        basename = (
            "Si_mp-149_primitive_Example-OPTICAL-GPAW-"
            "Result-Calculation-RPA-NLFC-AllData_xdirection.dat"
        )
        with _zip([
            ("silicon-gpaw-rpa/staged/" + basename, _dat()),
            ("silicon-gpaw-rpa/combined/" + basename, _dat()),
        ]) as archive:
            name, _, rows = _audit.read_case(
                archive, "staged", "GPAW", "NLFC", 3
            )
            self.assertEqual(name, "silicon-gpaw-rpa/staged/" + basename)
            self.assertEqual(len(rows), 3)

    def test_missing_case_and_wrong_point_count_fail(self):
        with _zip([("staged/optical-x.dat", _dat())]) as archive:
            with self.assertRaisesRegex(ValueError, "Expected exactly one"):
                _audit.read_case(archive, "combined", "QE", "", 3)
            with self.assertRaisesRegex(ValueError, "expected 4 rows"):
                _audit.read_case(archive, "staged", "QE", "", 4)

    def test_duplicate_matches_fail_instead_of_selecting_first(self):
        with _zip([
            ("staged/optical-x.dat", _dat()),
            ("extra/staged/optical-x.dat", _dat()),
        ]) as archive:
            with self.assertRaisesRegex(ValueError, "found 2"):
                _audit.read_case(archive, "staged", "QE", "", 3)

    def test_invalid_numeric_tables_fail(self):
        with self.assertRaisesRegex(ValueError, "non-increasing"):
            _audit.parse_table(_dat(
                ((0.001, 0.0), (4.0, 2.0), (3.0, 0.0))
            ))
        with self.assertRaisesRegex(ValueError, "non-finite"):
            _audit.parse_table(_dat(
                ((0.001, 0.0), (4.0, float("nan")), (8.0, 0.0))
            ))
        with self.assertRaisesRegex(ValueError, "seven columns"):
            _audit.parse_table(
                b"Energy(eV) Eps_img(1)\n0.0 1.0\n"
            )

    def test_analytic_trapezoid_and_sampling_diagnostic(self):
        rows = _audit.parse_table(_dat())
        self.assertAlmostEqual(_audit.spectral_area(rows), 7.999)
        self.assertAlmostEqual(_audit.peak_spacing(rows, 1), 4.0)
        self.assertAlmostEqual(_audit.nominal_fwhm("staged", "QE"), 0.05)
        self.assertAlmostEqual(_audit.nominal_fwhm("staged", "GPAW"), 0.10)
        self.assertAlmostEqual(_audit.nominal_fwhm("combined", "QE"), 0.10)
        self.assertAlmostEqual(_audit.nominal_fwhm("combined", "GPAW"), 0.20)
        with self.assertRaisesRegex(ValueError, "boundary"):
            _audit.peak_spacing(rows, 0)


if __name__ == "__main__":
    unittest.main()
