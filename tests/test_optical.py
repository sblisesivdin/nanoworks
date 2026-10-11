# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import unittest
from types import SimpleNamespace

from nanoworks.optical import validate_optical_settings


class TestOpticalSettings(unittest.TestCase):
    def config(self, **changes):
        values = dict(Engine='GPAW', Opt_calc_type='RPA', Opt_num_of_bands=16,
            Opt_eta=.05, Opt_FD_smearing=.05, Opt_min_en=0., Opt_max_en=20.,
            Opt_num_of_data=401, Opt_shift_en=0., Opt_cut_of_energy=100.,
            Opt_domega0=.05, Opt_omega2=5.,
            Opt_BSE_valence=range(0, 3), Opt_BSE_conduction=range(4, 7))
        values.update(changes)
        return SimpleNamespace(**values)

    def test_gpaw_rejects_unsupported_mode_and_spin_bse(self):
        for changes, message in ((dict(Mode='FD'), 'PW mode'),
                                 (dict(Mode='LCAO'), 'PW mode'),
                                 (dict(Opt_calc_type='BSE', Spin_calc=True), 'Spin_calc')):
            with self.subTest(changes=changes), self.assertRaisesRegex(ValueError, message):
                validate_optical_settings(self.config(**changes))
        validate_optical_settings(self.config(Spin_calc=True))
        validate_optical_settings(self.config(Engine='QE', Spin_calc=True))

    def test_shift_is_validated_only_when_used(self):
        for engine, method in (('GPAW', 'BSE'), ('QE', 'RPA')):
            for shift in (float('nan'), float('inf'), True):
                with self.subTest(engine=engine, shift=shift), self.assertRaisesRegex(ValueError, 'Opt_shift_en'):
                    validate_optical_settings(self.config(Engine=engine, Opt_calc_type=method, Opt_shift_en=shift))
            for shift in (-.5, 0., 1.25):
                self.assertEqual(validate_optical_settings(self.config(
                    Engine=engine, Opt_calc_type=method, Opt_shift_en=shift))['Opt_shift_en'], shift)
        validate_optical_settings(self.config(Opt_shift_en=float('nan')))

    def test_bse_rejects_invalid_band_selections(self):
        import numpy as np
        invalid = ([], [0, 0], [-1], [16], [1.5], [True],
                   [np.bool_(False)], [[1]], np.array([[1]]),
                   np.array([1.0]), '1', 1, {1}, {1: 2})
        for name in ('Opt_BSE_valence', 'Opt_BSE_conduction'):
            for bands in invalid:
                with self.subTest(name=name, bands=bands), self.assertRaisesRegex(ValueError, name):
                    validate_optical_settings(self.config(Opt_calc_type='BSE', **{name: bands}))

    def test_bse_preserves_explicit_band_order_and_allows_shared_bands(self):
        import numpy as np
        result = validate_optical_settings(self.config(Opt_calc_type='BSE',
            Opt_BSE_valence=np.array([2, 0], dtype=np.int64),
            Opt_BSE_conduction=(2, 3)))
        self.assertEqual(result['Opt_BSE_valence'], [2, 0])
        self.assertEqual(result['Opt_BSE_conduction'], [2, 3])
        validate_optical_settings(self.config(Opt_calc_type='BSE',
            Opt_BSE_valence=range(3), Opt_BSE_conduction=range(4, 16)))

    def test_unused_bse_selections_do_not_block_rpa(self):
        for engine in ('GPAW', 'QE'):
            validate_optical_settings(self.config(Engine=engine,
                Opt_BSE_valence=[], Opt_BSE_conduction='unused'))

    def test_invalid_common_parameters(self):
        for engine in ('GPAW', 'QE'):
            for change in (dict(Opt_eta=float('nan')), dict(Opt_eta=-1),
                           dict(Opt_num_of_bands=2.5), dict(Opt_num_of_bands=True),
                           dict(Opt_FD_smearing=-.1)):
                with self.subTest(engine=engine, change=change), self.assertRaises(ValueError):
                    validate_optical_settings(self.config(Engine=engine, **change))

    def test_grid_validation_applies_to_qe_and_gpaw_bse(self):
        for engine, method in (('QE', 'RPA'), ('GPAW', 'BSE')):
            for change in (dict(Opt_num_of_data=2.5), dict(Opt_num_of_data=True),
                           dict(Opt_num_of_data=1), dict(Opt_min_en=20),
                           dict(Opt_max_en=float('inf'))):
                with self.subTest(engine=engine, change=change), self.assertRaises(ValueError):
                    validate_optical_settings(self.config(Engine=engine, Opt_calc_type=method, **change))

    def test_unused_grid_and_backend_parameters_do_not_block_execution(self):
        validate_optical_settings(self.config(Opt_num_of_data=0, Opt_min_en=30, Opt_max_en=0))
        validate_optical_settings(self.config(Engine='QE', Opt_domega0=0, Opt_cut_of_energy=0))
        for name in ('Opt_domega0', 'Opt_omega2', 'Opt_cut_of_energy'):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, name):
                validate_optical_settings(self.config(**{name: 0}))

    def test_method_normalization_and_zero_smearing_semantics(self):
        self.assertEqual(validate_optical_settings(self.config(Opt_calc_type=' rpa '))['Opt_calc_type'], 'RPA')
        validate_optical_settings(self.config(Opt_FD_smearing=0))
        with self.assertRaisesRegex(ValueError, 'Opt_FD_smearing'):
            validate_optical_settings(self.config(Engine='QE', Opt_FD_smearing=0))

    def test_qe_broadening_must_match_epsilon_input_requirements(self):
        with self.assertRaisesRegex(ValueError, 'Opt_eta'):
            validate_optical_settings(self.config(Engine='QE', Opt_eta=0))


class TestOpticalSummary(unittest.TestCase):
    def test_summary_records_settings_actual_grid_and_table_hash(self):
        import hashlib
        import json
        import tempfile
        from pathlib import Path
        from nanoworks.optical import write_optical_summary
        with tempfile.TemporaryDirectory() as directory:
            table = Path(directory) / 'x.dat'
            table.write_text('columns\n0 3 4 2 1 0 .2\n1 3 4 2 1 100 .2\n')
            config = TestOpticalSettings().config(Engine='QE', XC_calc='PBE', Spin_calc=False)
            path = Path(directory) / 'summary.json'
            self.assertEqual(write_optical_summary(path, config, {'x': table}), str(path))
            report = json.loads(path.read_text())
            self.assertEqual(report['schema_version'], 2)
            self.assertEqual(report['tables']['x']['file'], 'x.dat')
            self.assertEqual(report['engine'], 'QE')
            self.assertEqual(report['method'], 'RPA')
            self.assertEqual(report['tables']['x']['npoints'], 2)
            self.assertEqual(report['tables']['x']['energy_max_ev'], 1.)
            self.assertEqual(report['tables']['x']['sha256'], hashlib.sha256(table.read_bytes()).hexdigest())
            self.assertEqual(report['settings']['Opt_num_of_data'], 401)
            self.assertEqual(report['columns'][5], 'Absorption(1/cm)')

    def test_invalid_table_retains_previous_summary(self):
        import tempfile
        from pathlib import Path
        from nanoworks.optical import write_optical_summary
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / 'summary.json'
            summary.write_text('previous completed summary')
            table = Path(directory) / 'invalid.dat'
            table.write_text('columns\n0 1 2\n1 2 3\n')
            config = TestOpticalSettings().config(Engine='QE', XC_calc='PBE', Spin_calc=False)
            with self.assertRaisesRegex(ValueError, 'Optical tables'):
                write_optical_summary(summary, config, {'x': table})
            self.assertEqual(summary.read_text(), 'previous completed summary')
            with self.assertRaisesRegex(ValueError, 'at least one'):
                write_optical_summary(summary, config, {})


class TestOpticalSummaryLoading(unittest.TestCase):
    def make_results(self, directory):
        from pathlib import Path
        from nanoworks.optical import write_optical_summary
        directory = Path(directory)
        table = directory / 'x.dat'
        table.write_text('columns\n0 3 4 2 1 0 .2\n1 3 4 2 1 100 .2\n')
        summary = directory / 'summary.json'
        config = TestOpticalSettings().config(Engine='QE', XC_calc='PBE', Spin_calc=False)
        write_optical_summary(summary, config, {'x': table})
        return summary, table

    def test_moved_results_use_paths_relative_to_summary(self):
        import tempfile
        from pathlib import Path
        from nanoworks.optical import load_optical_summary
        with tempfile.TemporaryDirectory() as directory:
            original = Path(directory) / 'original'; original.mkdir()
            summary, _ = self.make_results(original)
            moved = Path(directory) / 'moved'
            original.rename(moved)
            report, tables = load_optical_summary(moved / summary.name)
            self.assertEqual(report['schema_version'], 2)
            self.assertEqual(tables['x'].shape, (2, 7))

    def test_legacy_summary_finds_relocated_sibling_table(self):
        import json
        import tempfile
        from nanoworks.optical import load_optical_summary
        with tempfile.TemporaryDirectory() as directory:
            summary, _ = self.make_results(directory)
            report = json.loads(summary.read_text())
            report['schema_version'] = 1
            report['tables']['x']['file'] = '/old-computer/results/x.dat'
            summary.write_text(json.dumps(report))
            _, tables = load_optical_summary(summary)
            self.assertIn('x', tables)

    def test_changed_or_missing_table_is_rejected(self):
        import tempfile
        from nanoworks.optical import load_optical_summary
        with tempfile.TemporaryDirectory() as directory:
            summary, table = self.make_results(directory)
            table.write_text(table.read_text().replace('100', '200'))
            with self.assertRaisesRegex(ValueError, 'content differs'):
                load_optical_summary(summary)
            table.unlink()
            with self.assertRaisesRegex(FileNotFoundError, 'missing'):
                load_optical_summary(summary)

    def test_malformed_schema_and_grid_metadata_are_rejected(self):
        import json
        import tempfile
        from nanoworks.optical import load_optical_summary
        with tempfile.TemporaryDirectory() as directory:
            summary, _ = self.make_results(directory)
            original = json.loads(summary.read_text())
            for field, value in (('schema_version', True), ('schema_version', 99),
                                 ('status', 'failed'), ('columns', []), ('tables', {})):
                with self.subTest(field=field):
                    report = dict(original); report[field] = value
                    summary.write_text(json.dumps(report))
                    with self.assertRaises(ValueError):
                        load_optical_summary(summary)
            original['tables']['x']['npoints'] = 3
            summary.write_text(json.dumps(original))
            with self.assertRaisesRegex(ValueError, 'grid differs'):
                load_optical_summary(summary)


class TestOpticalTables(unittest.TestCase):
    def test_finite_negative_dielectric_values_are_preserved(self):
        import numpy as np
        from nanoworks.optical import validate_optical_table
        data = np.array([[0, -2, 1, .3, 1.4, 0, .5], [1, -1, 1, .4, 1.1, 2, .3]])
        result = validate_optical_table(data)
        np.testing.assert_array_equal(result, data)
        data[0, 1] = 99
        self.assertEqual(result[0, 1], -2)

    def test_malformed_or_nonfinite_spectra_are_rejected(self):
        import numpy as np
        from nanoworks.optical import validate_optical_table
        valid = np.ones((2, 7)); valid[:, 0] = [0, 1]
        cases = [valid[:, :6], valid[:1], valid.tolist()[0],
                 valid.copy(), valid.copy(), valid.copy()]
        cases[3][1, 1] = np.nan
        cases[4][:, 0] = [1, 0]
        cases[5][:, 0] = [-1, 0]
        for table in cases:
            with self.subTest(table=table), self.assertRaisesRegex(ValueError, 'Optical tables'):
                validate_optical_table(table)


class TestOpticalExportMPI(unittest.TestCase):
    def test_root_failure_is_broadcast(self):
        from unittest.mock import Mock, patch
        from nanoworks.optical import run_optical_exports
        callback = Mock(side_effect=OSError('disk full'))
        with patch('ase.parallel.world', SimpleNamespace(rank=0)):
            with patch('ase.parallel.broadcast', side_effect=lambda value, **kwargs: value) as broadcast:
                with self.assertRaisesRegex(RuntimeError, 'Optical export failed.*disk full'):
                    run_optical_exports(callback)
        callback.assert_called_once()
        self.assertEqual(broadcast.call_args.args[0][0], 'OSError: disk full')

    def test_nonroot_receives_failure_without_opening_files(self):
        from unittest.mock import Mock, patch
        from nanoworks.optical import run_optical_exports
        callback = Mock()
        with patch('ase.parallel.world', SimpleNamespace(rank=1)):
            with patch('ase.parallel.broadcast', return_value=('OSError: disk full', None)):
                with self.assertRaisesRegex(RuntimeError, 'disk full'):
                    run_optical_exports(callback)
        callback.assert_not_called()


class TestDerivedOpticalProperties(unittest.TestCase):
    def test_analytic_dielectric_response_and_units(self):
        import numpy as np
        from nanoworks.optical import derive_optical_table
        table = derive_optical_table([0, 1], [3, 3], [4, 4])
        np.testing.assert_array_equal(table[:, 3], [2, 2])
        np.testing.assert_array_equal(table[:, 4], [1, 1])
        np.testing.assert_allclose(table[:, 6], [.2, .2])
        self.assertAlmostEqual(table[1, 5], 101354.61433096675)
        self.assertEqual(table[0, 5], 0)

    def test_weak_absorption_survives_dielectric_subtraction_cancellation(self):
        import numpy as np
        from nanoworks.optical import derive_optical_table
        table = derive_optical_table([0, 1], [4, 4], [1e-12, 1e-12])
        np.testing.assert_allclose(table[:, 4], [2.5e-13, 2.5e-13], rtol=1e-14, atol=0)
        self.assertGreater(table[1, 5], 0)
        np.testing.assert_allclose(table[:, 3] ** 2 - table[:, 4] ** 2, [4, 4])
        np.testing.assert_allclose(2 * table[:, 3] * table[:, 4], [1e-12, 1e-12], atol=0)

    def test_negative_real_dielectric_and_sign_convention(self):
        import numpy as np
        from nanoworks.optical import derive_optical_table
        table = derive_optical_table([0, 1], [-4, -4], [0, 0])
        np.testing.assert_array_equal(table[:, 3], [0, 0])
        np.testing.assert_array_equal(table[:, 4], [2, 2])
        np.testing.assert_array_equal(table[:, 6], [1, 1])
        signed = derive_optical_table([0, 1], [3, 3], [-4, -4])
        np.testing.assert_array_equal(signed[:, 2], [-4, -4])
        np.testing.assert_array_equal(signed[:, 4], [1, 1])

    def test_nonfinite_or_mismatched_dielectric_data_is_rejected(self):
        from nanoworks.optical import derive_optical_table
        for values in (([0, 1], [3], [4, 4]), ([0, 1], [3, float('nan')], [4, 4])):
            with self.subTest(values=values), self.assertRaisesRegex(ValueError, 'Dielectric spectra'):
                derive_optical_table(*values)
