# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression coverage for phonon tables, diagnostics and required exports."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

import numpy as np

from nanoworks.phonon_results import write_gpaw_phonon_results


class TestGPAWPhononResults(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.prefix = str(Path(temporary.name) / 'Ni-PHONON-GPAW')
        self.phonon = Mock()
        self.phonon.get_frequencies.return_value = [-.02, .1, .2]
        self.phonon.mesh = SimpleNamespace(qpoints=np.array([[0, 0, 0], [.5, 0, 0]]),
            weights=np.array([1, 3]), frequencies=np.array([[-.02, .1, .2], [-.3, .4, .5]]))
        self.phonon.total_dos = SimpleNamespace(
            frequency_points=np.array([-.3, 0., .5]), dos=np.array([.1, 1., .2]))
        self.phonon.get_band_structure_dict.return_value = {
            'distances': [np.array([0., .5])],
            'frequencies': [np.array([[-.02, .1, .2], [-.3, .4, .5]])]}
        self.band_path = ([np.array([[0, 0, 0], [.5, 0, 0]])], ['G', 'X'], [False])
        self.phonon.get_thermal_properties_dict.return_value = {
            'temperatures': [0., 100.], 'free_energy': [.1, .2],
            'entropy': [0., 1.], 'heat_capacity': [0., 2.]}

    def test_signed_tables_and_summary_use_completed_mesh_and_band_data(self):
        report = write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path,
                                          [2, 2, 2], (0, 100, 100))
        self.phonon.run_mesh.assert_called_once_with([2, 2, 2])
        calls = [call[0] for call in self.phonon.method_calls]
        self.assertLess(calls.index('run_mesh'), calls.index('run_total_dos'))
        self.phonon.run_total_dos.assert_called_once()
        self.assertEqual(report['status'], 'complete')
        self.assertEqual(report['gamma_frequencies_thz'], [-.02, .1, .2])
        self.assertEqual(report['minimum_mesh_qpoint'], [.5, 0, 0])
        self.assertEqual(report['negative_mesh_mode_count'], 2)
        self.assertEqual(report['mesh_modes_below_reporting_threshold'], 1)
        self.assertAlmostEqual(report['weighted_mesh_fraction_below_reporting_threshold'], .25)
        self.assertAlmostEqual(report['minimum_band_frequency_thz'], -.3)
        np.testing.assert_allclose(np.loadtxt(report['mesh_data_file'])[:, 4:],
                                   self.phonon.mesh.frequencies)
        np.testing.assert_allclose(np.loadtxt(report['dos_data_file'])[:, 0], [-.3, 0., .5])
        np.testing.assert_allclose(np.loadtxt(report['band_data_file'], skiprows=1)[:, 1:],
                                   self.phonon.get_band_structure_dict.return_value['frequencies'][0])
        thermal = np.loadtxt(report['thermal_data_file'], delimiter=',', skiprows=1)
        np.testing.assert_array_equal(thermal[:, 0], [0, 100])
        self.phonon.write_yaml_band_structure.assert_called_once()
        self.phonon.write_yaml_thermal_properties.assert_called_once()
        self.phonon.plot_band_structure_and_dos.return_value.savefig.assert_called_once()

    def test_required_yaml_export_failure_is_not_reported_as_success(self):
        self.phonon.write_yaml_band_structure.side_effect = OSError('disk full')
        with self.assertRaisesRegex(OSError, 'disk full'):
            write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2])
        self.phonon.plot_band_structure_and_dos.assert_not_called()

    def test_invalid_dos_stops_before_exporting_bands(self):
        self.phonon.total_dos.dos = np.array([np.nan, 1., .2])
        with self.assertRaisesRegex(ValueError, 'finite frequency and DOS'):
            write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2])
        self.phonon.run_band_structure.assert_not_called()

    def test_thermal_output_is_optional(self):
        report = write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2])
        self.assertIsNone(report['thermal_data_file'])
        self.assertIsNone(report['thermal_yaml_file'])
        self.phonon.run_thermal_properties.assert_not_called()
