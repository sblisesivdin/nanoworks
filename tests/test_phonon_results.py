# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression coverage for phonon tables, diagnostics and required exports."""

import tempfile
import json
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from nanoworks.phonon_results import (
    qpoint_frequencies, write_gpaw_phonon_results, prepare_gpaw_postprocess_plan, postprocess_gpaw_plan,
)


class TestQpointFrequencies(unittest.TestCase):
    def test_batch_preserves_order_negative_modes_and_detaches_results(self):
        modes = np.array([[-.3, .2, .4], [.1, .5, .9]])
        phonon = Mock()
        phonon.run_qpoints.return_value = SimpleNamespace(frequencies=modes)
        points = [[.5, 0, 0], [0, 0, 0]]
        result = qpoint_frequencies(phonon, points)
        np.testing.assert_array_equal(result, modes)
        phonon.run_qpoints.assert_called_once()
        np.testing.assert_array_equal(phonon.run_qpoints.call_args.args[0], points)
        modes[:] = 99
        np.testing.assert_array_equal(result, [[-.3, .2, .4], [.1, .5, .9]])
        phonon.get_frequencies.assert_not_called()

    def test_earlier_api_reads_stored_qpoint_result(self):
        phonon = Mock()
        phonon.run_qpoints.return_value = None
        phonon.qpoints = SimpleNamespace(frequencies=[[-.02, .1, .2]])
        np.testing.assert_array_equal(qpoint_frequencies(phonon, [[0, 0, 0]]),
                                      [[-.02, .1, .2]])
        phonon.get_frequencies.assert_not_called()

    def test_invalid_qpoints_stop_before_frequency_calculation(self):
        for points in ([], [0, 0, 0], [[0, 0]], [[0, np.nan, 0]]):
            with self.subTest(points=points):
                phonon = Mock()
                with self.assertRaisesRegex(ValueError, 'q-points'):
                    qpoint_frequencies(phonon, points)
                phonon.run_qpoints.assert_not_called()

    def test_incomplete_or_nonfinite_modes_are_rejected(self):
        for modes in ([.1, .2, .3], [[.1, .2, .3], [.4, .5, .6]],
                      [[]], [[np.nan, .1, .2]], [[np.inf, .1, .2]]):
            with self.subTest(modes=modes):
                phonon = Mock()
                phonon.run_qpoints.return_value = SimpleNamespace(frequencies=modes)
                with self.assertRaisesRegex(ValueError, 'matching finite mode'):
                    qpoint_frequencies(phonon, [[0, 0, 0]])

    def test_missing_stored_result_is_rejected(self):
        phonon = Mock()
        phonon.run_qpoints.return_value = None
        phonon.qpoints = None
        with self.assertRaisesRegex(ValueError, 'did not produce'):
            qpoint_frequencies(phonon, [[0, 0, 0]])


class TestGPAWPhononResults(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.prefix = str(Path(temporary.name) / 'Ni-PHONON-GPAW')
        self.phonon = Mock()
        self.phonon.run_qpoints.return_value = SimpleNamespace(
            frequencies=np.array([[-.02, .1, .2]]))
        self.phonon.get_frequencies.side_effect = AssertionError("Deprecated frequency API used")
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

    def test_invalid_band_path_stops_before_frequency_work_and_exports(self):
        invalid = (None, [], ([], [], []), ([[[0, 0, 0]]], [], []),
                   ([[[0, 0, 0], [0.5, 0]]], [], []),
                   ([[[0, 0, 0], [np.nan, 0, 0]]], [], []))
        for path in invalid:
            with self.subTest(path=path):
                with self.assertRaisesRegex(ValueError, 'band'):
                    write_gpaw_phonon_results(self.phonon, self.prefix, path, [2, 2, 2])
                self.phonon.run_qpoints.assert_not_called()
                self.phonon.run_mesh.assert_not_called()
                self.assertEqual(list(Path(self.prefix).parent.iterdir()), [])

    def test_mesh_mode_count_must_match_gamma_before_mesh_export(self):
        self.phonon.mesh.frequencies = np.ones((2, 2))
        with self.assertRaisesRegex(ValueError, 'mesh mode count'):
            write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2])
        self.assertFalse(Path(self.prefix + '-Result-Mesh-THz.dat').exists())
        self.phonon.run_total_dos.assert_not_called()

    def test_invalid_band_results_preserve_previous_band_table(self):
        target = Path(self.prefix + '-Result-Band.dat')
        target.write_text('previous completed band table')
        self.phonon.get_band_structure_dict.return_value['frequencies'][0] = np.ones((2, 2))
        with self.assertRaisesRegex(ValueError, 'matching finite distances and modes'):
            write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2])
        self.assertEqual(target.read_text(), 'previous completed band table')
        self.phonon.write_yaml_band_structure.assert_not_called()

    def test_two_dimensional_thermal_column_is_rejected_before_thermal_csv(self):
        self.phonon.get_thermal_properties_dict.return_value['free_energy'] = [[.1, .2], [.3, .4]]
        with self.assertRaisesRegex(ValueError, 'four matching finite arrays'):
            write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2], [0, 100, 100])
        self.assertFalse(Path(self.prefix + '-Result-Thermal-Properties.csv').exists())
        self.phonon.write_yaml_thermal_properties.assert_not_called()

    def test_thermal_output_is_optional(self):
        report = write_gpaw_phonon_results(self.phonon, self.prefix, self.band_path, [2, 2, 2])
        self.assertIsNone(report['thermal_data_file'])
        self.assertIsNone(report['thermal_yaml_file'])
        self.phonon.run_thermal_properties.assert_not_called()


class TestGPAWPostprocessPlan(unittest.TestCase):
    def setUp(self):
        try:
            from phonopy import Phonopy
            from phonopy.structure.atoms import PhonopyAtoms
        except ModuleNotFoundError:
            self.skipTest('phonopy is optional')
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.prefix = str(Path(temporary.name) / 'Ni-PHONON-GPAW')
        cell = PhonopyAtoms(symbols=['Ni', 'Ni'], cell=np.eye(3) * 4,
            scaled_positions=[[0, 0, 0], [.5, .5, .5]], magnetic_moments=[2, -2], masses=[60, 60])
        self.matrix = np.diag([2, 1, 1])
        self.phonon = Phonopy(cell, self.matrix)
        count = len(self.phonon.supercell)
        constants = np.zeros((count, count, 3, 3))
        for index in range(count):
            constants[index, index] = np.eye(3) * 2
        self.phonon.force_constants = constants
        np.save(self.prefix + '-Result-Force-Constants.npy', constants)
        self.path = ([np.array([[0, 0, 0], [.5, 0, 0]])], ['G', 'X'], [False])
        self.plan = prepare_gpaw_postprocess_plan(self.phonon, self.prefix, self.matrix,
            self.path, [2, 2, 2], None, {'xc': 'PBE', 'spin_polarized': True,
                                       'reused_force_constants': False})

    def report(self):
        return json.loads(Path(self.prefix + '-Result-Summary.json').read_text())

    def test_archived_geometry_and_constants_reconstruct_without_gpaw(self):
        filename = self.prefix + '-Input-Postprocess.json'
        stored = json.loads(Path(filename).read_text())
        with patch.dict('sys.modules', {'gpaw': None}):
            with patch('nanoworks.phonon_results.write_gpaw_phonon_results',
                       return_value={'status': 'complete', 'engine': 'GPAW'}) as export:
                result = postprocess_gpaw_plan(stored)
        reconstructed = export.call_args.args[0]
        np.testing.assert_array_equal(reconstructed.force_constants, self.phonon.force_constants)
        np.testing.assert_allclose(qpoint_frequencies(reconstructed, [[.2, 0, 0]]),
                                   qpoint_frequencies(self.phonon, [[.2, 0, 0]]))
        np.testing.assert_array_equal(reconstructed.unitcell.masses, self.phonon.unitcell.masses)
        self.assertTrue(result['analysis_only'])
        self.assertTrue(result['reused_force_constants'])
        self.assertEqual(self.report()['status'], 'complete')

    def test_changed_force_constants_are_rejected_before_exports(self):
        np.save(self.plan['force_constants_file'], self.phonon.force_constants * 2)
        with patch('nanoworks.phonon_results.write_gpaw_phonon_results') as export:
            with self.assertRaisesRegex(ValueError, 'force constants changed'):
                postprocess_gpaw_plan(self.plan)
        export.assert_not_called()
        self.assertEqual(self.report()['status'], 'failed')

    def test_changed_geometry_is_rejected(self):
        self.plan['unitcell']['cell'][0][0] += 1
        with self.assertRaisesRegex(ValueError, 'geometry or provenance changed'):
            postprocess_gpaw_plan(self.plan)

    def test_changed_phonopy_version_is_rejected(self):
        with patch('phonopy.__version__', 'different-version'):
            with self.assertRaisesRegex(ValueError, 'Phonopy version changed'):
                postprocess_gpaw_plan(self.plan)

    def test_analysis_mesh_can_change_without_altering_physical_snapshot(self):
        self.plan['dos_mesh'] = [4, 4, 4]
        with patch('nanoworks.phonon_results.write_gpaw_phonon_results',
                   return_value={'status': 'complete', 'engine': 'GPAW'}) as export:
            postprocess_gpaw_plan(self.plan)
        self.assertEqual(export.call_args.args[3], [4, 4, 4])

    def test_invalid_analysis_mesh_is_rejected_without_exports(self):
        self.plan['dos_mesh'] = [2.5, 2, 2]
        with patch('nanoworks.phonon_results.write_gpaw_phonon_results') as export:
            with self.assertRaisesRegex(ValueError, 'Phonon_qpts'):
                postprocess_gpaw_plan(self.plan)
        export.assert_not_called()

    def test_invalid_edited_band_path_is_rejected_before_exports(self):
        self.plan['band_path'] = ([[[0, 0, 0], [float('inf'), 0, 0]]], ['G', 'X'], [False])
        with patch('nanoworks.phonon_results.write_gpaw_phonon_results') as export:
            with self.assertRaisesRegex(ValueError, 'band segment'):
                postprocess_gpaw_plan(self.plan)
        export.assert_not_called()
        self.assertEqual(self.report()['status'], 'failed')
        self.assertTrue(Path(self.plan['force_constants_file']).is_file())

    def test_edited_band_path_can_change_without_altering_physical_snapshot(self):
        self.plan['band_path'] = ([[[0, 0, 0], [.25, .25, 0], [.5, .5, 0]]], ['G', 'M'], [False])
        with patch('nanoworks.phonon_results.write_gpaw_phonon_results',
                   return_value={'status': 'complete', 'engine': 'GPAW'}) as export:
            postprocess_gpaw_plan(self.plan)
        np.testing.assert_array_equal(export.call_args.args[2][0][0],
                                      [[0, 0, 0], [.25, .25, 0], [.5, .5, 0]])
        self.assertEqual(self.report()['status'], 'complete')

    def test_export_failure_preserves_the_archived_retry_plan(self):
        with patch('nanoworks.phonon_results.write_gpaw_phonon_results', side_effect=OSError('disk full')):
            with self.assertRaisesRegex(OSError, 'disk full'):
                postprocess_gpaw_plan(self.plan)
        self.assertEqual(self.report()['status'], 'failed')
        self.assertTrue(Path(self.prefix + '-Input-Postprocess.json').is_file())
        self.assertTrue(Path(self.plan['force_constants_file']).is_file())


class TestPhononNumericalTables(unittest.TestCase):
    def test_projected_dos_keeps_signed_grid_and_atom_rows(self):
        from nanoworks.phonon_results import validate_projected_dos
        grid, values = validate_projected_dos([-.2, 0, .2], [[1, 2, 3], [4, 5, 6]], 2)
        np.testing.assert_array_equal(grid, [-.2, 0, .2])
        np.testing.assert_array_equal(values.sum(axis=0), [5, 7, 9])

    def test_invalid_projected_dos_cannot_be_exported_as_complete(self):
        from nanoworks.phonon_results import validate_projected_dos
        invalid = (([0, 1], [[1, 2]], 2), ([0, 1], [[1, 2, 3]], 1),
                   ([0, 0], [[1, 2]], 1), ([1, 0], [[1, 2]], 1),
                   ([0, np.inf], [[1, 2]], 1), ([0, 1], [[1, np.nan]], 1))
        for grid, values, count in invalid:
            with self.subTest(grid=grid, values=values), self.assertRaisesRegex(ValueError, 'projected DOS'):
                validate_projected_dos(grid, values, count)

    def test_thermal_columns_require_matching_finite_one_dimensional_arrays(self):
        from nanoworks.phonon_results import validate_thermal_table
        valid = dict(temperatures=[0, 100], free_energy=[.1, .2], entropy=[0, 1], heat_capacity=[0, 2])
        self.assertEqual(validate_thermal_table(valid).shape, (2, 4))
        for changes in (dict(free_energy=[[1, 2], [3, 4]]), dict(entropy=[0]),
                        dict(heat_capacity=[0, np.nan]), dict(temperatures=[100, 0]),
                        dict(temperatures=[0, 0]), dict(temperatures=[-1, 100])):
            with self.subTest(changes=changes), self.assertRaisesRegex(ValueError, 'thermal properties'):
                validate_thermal_table({**valid, **changes})
