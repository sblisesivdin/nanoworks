# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""Offline checks for the Si QE/GPAW matched-width diagnostic inputs."""

import ast
import contextlib
import csv
import importlib.util
import io
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


_ROOT = Path(__file__).resolve().parents[1]
_EXAMPLE = _ROOT / "nanoworks/examples/Si-2atoms-optical"
_SCRIPT = _ROOT / "docs/source/_static/examples/compare_silicon_rpa_runs.py"
_spec = importlib.util.spec_from_file_location("si_rpa_comparison", _SCRIPT)
assert _spec is not None and _spec.loader is not None
_comp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_comp)


def _input_settings(name):
    path = _EXAMPLE / name
    settings = {}
    for node in ast.parse(path.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name):
                settings[target.id] = ast.literal_eval(node.value)
    return settings


def _dat(energies, images):
    header = (
        "Energy(eV) Eps_real(1) Eps_img(1) Refractive_Index(1) "
        "Extinction_Index(1) Absorption(1/cm) Reflectivity(1)"
    )
    return header + "\n" + "".join(
        f"{e} 1.0 {y} 1.0 0.0 0.0 0.0\n"
        for e, y in zip(energies, images)
    )


class MatchedWidthInputTests(unittest.TestCase):
    def test_engine_and_ground_state_isolation(self):
        q = _input_settings("Si-QE-RPA-matched-width.py")
        g = _input_settings("Si-GPAW-RPA-matched-width.py")
        self.assertEqual((q["Engine"], g["Engine"]), ("QE", "GPAW"))
        self.assertNotEqual(q["Outdirname"], g["Outdirname"])
        for c in (q, g):
            self.assertTrue(c["Ground_calc"])
            self.assertTrue(c["Optical_calc"])
            for feature in (
                "Geo_optim", "Elastic_calc", "DOS_calc", "Band_calc",
                "Density_calc", "Phonon_calc", "SOC_calc", "Spin_calc",
            ):
                self.assertFalse(c[feature], feature)
            self.assertEqual(c["XC_calc"], "PBE")
            self.assertEqual(c["Wavefunction_cutoff"], 340)
            self.assertEqual(c["Opt_calc_type"], "RPA")
            self.assertEqual(c["Opt_num_of_bands"], 8)
            self.assertTrue(c["Gamma"])
            self.assertTrue(c["Opt_gamma"])
            for axis in "xyz":
                self.assertEqual(c[f"Ground_kpts_{axis}"], 2)
                self.assertEqual(c[f"Opt_kpts_{axis}"], 2)
            self.assertEqual(c["Occupation_scheme"], "fermi-dirac")
            self.assertEqual(c["Smearing_width"], 0.05)
            self.assertEqual(c["Opt_FD_smearing"], 0.05)

    def test_matched_isolated_linewidth_and_grid_resolution(self):
        q = _input_settings("Si-QE-RPA-matched-width.py")
        g = _input_settings("Si-GPAW-RPA-matched-width.py")
        target = 0.10
        self.assertAlmostEqual(q["Opt_eta"], target)
        self.assertAlmostEqual(2 * g["Opt_eta"], target)
        self.assertEqual(q["Opt_num_of_data"], 1001)
        self.assertEqual(q["Opt_max_en"], 10.0)
        self.assertEqual(g["Opt_cut_of_energy"], 30)
        self.assertIsNone(g["Opt_nblocks"])
        self.assertAlmostEqual(g["Opt_domega0"], 0.01)
        self.assertEqual(g["Opt_omega2"], 20.0)
        qe_step = 10.0 / (q["Opt_num_of_data"] - 1)
        self.assertGreater(target / qe_step, 7.0)
        # GPAW documented non-linear grid, local derivative d(omega)/di.
        gpaw_step_at_8 = (
            g["Opt_domega0"]
            * (1 + (math.sqrt(2) - 1) * 8 / g["Opt_omega2"]) ** 2
        )
        self.assertGreater(target / gpaw_step_at_8, 7.0)


class MatchedWidthComparatorTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.dat = Path(self.temporary.name) / "optical.dat"
        self.dat.write_text(
            _dat((0, 1, 2, 3), (0, 2, 2, 0)), encoding="utf-8"
        )

    def test_valid_native_grid_moments_and_peak(self):
        rows = _comp.read_table(self.dat)
        self.assertEqual(rows, [(0, 0), (1, 2), (2, 2), (3, 0)])
        self.assertAlmostEqual(_comp.value_at(rows, 0.5), 1.0)
        self.assertEqual(_comp.moments(rows, 0.5, 2.5), (3.5, 5.25))
        stats = _comp.inspect(self.dat, 0.5, 2.5, 0.10)
        self.assertEqual(stats["points"], 4)
        self.assertEqual(stats["peak_eV"], 1.0)
        self.assertAlmostEqual(stats["nominal_FWHM_to_step"], 0.1)
        self.assertAlmostEqual(stats["max_step_in_window_eV"], 1.0)
        self.assertAlmostEqual(stats["min_nominal_FWHM_to_step"], 0.1)
        self.assertEqual(len(stats["sha256"]), 64)

    def test_first_moment_uses_exact_piecewise_linear_integral(self):
        self.dat.write_text(
            _dat((0, 1, 2), (0, 1, 2)), encoding="utf-8"
        )
        rows = _comp.read_table(self.dat)
        area, first_moment = _comp.moments(rows, 0.0, 2.0)
        self.assertAlmostEqual(area, 2.0)
        self.assertAlmostEqual(first_moment, 8.0 / 3.0)
        # Plain trapezoids of E*epsilon2 would incorrectly give 3.0.
        self.assertNotAlmostEqual(first_moment, 3.0)

    def test_window_max_step_uses_all_intervals_not_only_peak(self):
        self.dat.write_text(
            _dat((0, 0.2, 0.7, 2.0), (0, 2, 3, 0)),
            encoding="utf-8",
        )
        rows = _comp.read_table(self.dat)
        self.assertAlmostEqual(_comp.window_max_step(rows, 0.1, 1.4), 0.7)
        stats = _comp.inspect(self.dat, 0.1, 1.4, 0.1)
        self.assertAlmostEqual(stats["max_step_in_window_eV"], 0.7)
        self.assertAlmostEqual(stats["min_nominal_FWHM_to_step"], 1/7)

    def test_peak_is_selected_within_the_moment_window(self):
        # The global maximum is at 4 eV, outside the requested 0--3 eV
        # comparison window. The reported peak must be the 1 eV maximum.
        self.dat.write_text(
            _dat((0, 1, 2, 3, 4, 5), (0, 4, 1, 1, 20, 0)),
            encoding="utf-8",
        )
        stats = _comp.inspect(self.dat, 0.0, 3.0, 0.10)
        self.assertAlmostEqual(stats["peak_eV"], 1.0)
        self.assertAlmostEqual(stats["peak_epsilon2"], 4.0)
        self.assertEqual(stats["points"], 6)
        with self.assertRaisesRegex(ValueError, "selected window boundary"):
            _comp.inspect(self.dat, 0.0, 0.5, 0.10)

    def test_invalid_tables_and_uncovered_window_fail(self):
        with self.assertRaisesRegex(ValueError, "exceeds a data grid"):
            _comp.moments(_comp.read_table(self.dat), 0.0, 8.0)
        with self.assertRaisesRegex(ValueError, "interval"):
            _comp.moments(_comp.read_table(self.dat), 2.0, 1.0)
        for text, pattern in (
            (_dat((0, 1, 1), (1, 2, 3)), "not increasing"),
            (_dat((0, 1, 2), (1, float("nan"), 3)), "non-finite"),
            ("Energy(eV) Eps_img(1)\n0.0 1.0\n", "expected 7"),
        ):
            self.dat.write_text(text, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, pattern):
                _comp.read_table(self.dat)

    def test_manual_csv_output_and_user_chosen_window(self):
        capture = io.StringIO()
        diagnostics = io.StringIO()
        argv = [
            "compare_silicon_rpa_runs.py", str(self.dat), str(self.dat),
            "--emin=0.5", "--emax=2.5", "--fwhm=0.10",
        ]
        with (
            patch.object(sys, "argv", argv),
            contextlib.redirect_stdout(capture),
            contextlib.redirect_stderr(diagnostics),
        ):
            _comp.main()
        output = capture.getvalue()
        rows = list(csv.DictReader(io.StringIO(output)))
        self.assertEqual(len(rows), 2)
        self.assertEqual([row["engine"] for row in rows], ["QE", "GPAW_NLFC"])
        self.assertAlmostEqual(float(rows[0]["integral_epsilon2_dE_eV"]), 3.5)
        self.assertAlmostEqual(float(rows[0]["integral_E_epsilon2_dE_eV2"]), 5.25)
        self.assertEqual(rows[1]["resolution_warning"], "YES")
        self.assertEqual(len(rows[0]["source_sha256"]), 64)
        self.assertNotIn("#", output)
        self.assertIn(
            "[0.500000, 2.500000] eV area ratio: 1.000000",
            diagnostics.getvalue(),
        )
        self.assertIn("WARNING: GPAW_NLFC is below 5", diagnostics.getvalue())


if __name__ == "__main__":
    unittest.main()
