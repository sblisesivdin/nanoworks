# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from nanoworks.result_plotting import (
    load_dos_table, plot_dos_result, plot_saved_result, load_band_table, plot_band_result,
    load_pdos_table, plot_pdos_result,
)


class TestSavedResultPlotting(unittest.TestCase):
    def test_new_and_legacy_pdos_headers_preserve_all_orbital_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-DOS-QE-Result-PDOS.dat'
            path.write_text('Energy(eV) s-total(states/eV) p_j0.5(states/eV) TOTAL(states/eV)\n'
                            '-1 1 2 3\n1 2 3 5\n')
            csv = path.with_suffix('.csv')
            csv.write_text('Energy,s-total,p_j0.5,TOTAL\n-1,1,2,3\n1,2,3,5\n')
            columns, data = load_pdos_table(path)
            legacy_columns, legacy_data = load_pdos_table(csv)
            self.assertEqual(columns, legacy_columns)
            np.testing.assert_array_equal(data, legacy_data)

    def test_spin_pdos_column_mismatch_preserves_existing_figure(self):
        with tempfile.TemporaryDirectory() as directory:
            up = Path(directory) / 'Si-DOS-GPAW-Result-PDOS-Up.dat'
            down = up.with_name('Si-DOS-GPAW-Result-PDOS-Down.dat')
            up.write_text('Energy(eV) s-total(states/eV)\n-1 1\n1 2\n')
            down.write_text('Energy(eV) p-total(states/eV)\n-1 1\n1 2\n')
            output = up.with_name('Si-DOS-GPAW-Graph-PDOS.png')
            output.write_bytes(b'previous figure')
            with patch('nanoworks.result_plotting._save_figure') as save:
                with self.assertRaisesRegex(ValueError, 'identical columns'):
                    plot_pdos_result(up)
            save.assert_not_called()
            self.assertEqual(output.read_bytes(), b'previous figure')

    def test_malformed_pdos_header_and_grid_are_rejected(self):
        cases = ['Energy(meV) s-total(states/eV)\n0 1\n1 2\n',
                 'Energy(eV) s-total(states/eV)\n0 1\n0 2\n',
                 'Energy(eV) s-total(states/eV)\n0 nan\n1 2\n']
        for content in cases:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'Si-DOS-QE-Result-PDOS.dat'
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_pdos_table(path)

    def test_soc_dos_uses_one_saved_channel(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'WSe2-DOS-GPAW-Result-DOS-SOC.dat'
            path.write_text('Energy(eV) DOS(states/eV)\n-1 1\n1 2\n')
            with patch('nanoworks.result_plotting._save_figure') as save:
                plot_saved_result(path)
            self.assertEqual(save.call_args.args[1], path.with_name('WSe2-DOS-GPAW-Graph-DOS-SOC.png'))

    def test_long_and_wide_band_tables_preserve_saved_energies(self):
        with tempfile.TemporaryDirectory() as directory:
            long = Path(directory) / 'Si-BAND-GPAW-Result-Band.dat'
            wide = long.with_name('Si-BAND-GPAW-Result-Band-XYYY.dat')
            long.write_text('KPointIndex(1) Energy(eV)\n0 -2\n1 -1\n\n0 3\n1 4\n')
            wide.write_text('KPointIndex(1) Band_1(eV) Band_2(eV)\n0 -2 3\n1 -1 4\n')
            for path in (long, wide):
                indices, energies = load_band_table(path)
                np.testing.assert_array_equal(indices, [0, 1])
                np.testing.assert_array_equal(energies, [[-2, 3], [-1, 4]])
            long.write_text('0 -2\n1 -1\n\n0 3\n1 4\n')
            np.testing.assert_array_equal(load_band_table(long)[1], [[-2, 3], [-1, 4]])

    def test_inconsistent_or_nonfinite_band_blocks_are_rejected(self):
        cases = ['0 1\n1 2\n\n0 3\n2 4\n', '0 1\n1 nan\n',
                 '0.5 1\n1.5 2\n', '1 1\n0 2\n']
        for content in cases:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'Si-BAND-QE-Result-Band.dat'
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_band_table(path)

    def test_band_save_failure_preserves_existing_png_and_closes_figure(self):
        import matplotlib.pyplot as plt
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-BAND-QE-Result-Band.dat'
            path.write_text('KPointIndex(1) Energy(eV)\n0 -1\n1 1\n')
            output = path.with_name('Si-BAND-QE-Graph-Band-Index.png')
            output.write_bytes(b'previous figure')
            before = set(plt.get_fignums())
            with patch('matplotlib.figure.Figure.savefig', side_effect=OSError('disk full')):
                with self.assertRaises(OSError):
                    plot_band_result(path)
            self.assertEqual(output.read_bytes(), b'previous figure')
            self.assertEqual(set(plt.get_fignums()), before)
            self.assertFalse(Path(str(output) + '.tmp').exists())

    def test_missing_spin_band_partner_is_rejected_before_drawing(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-BAND-GPAW-Result-Band-Up.dat'
            path.write_text('KPointIndex(1) Energy(eV)\n0 -1\n1 1\n')
            with patch('nanoworks.result_plotting._save_figure') as save:
                with self.assertRaises(OSError):
                    plot_saved_result(path)
            save.assert_not_called()

    def test_dos_energy_reference_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-DOS-QE-Result-DOS.csv'
            data = np.array([[-3., 1.], [2., 4.]])
            np.savetxt(path, data, delimiter=',')
            np.testing.assert_array_equal(load_dos_table(path), data)

    def test_invalid_dos_tables_are_rejected(self):
        for content in ('0,1\n', '0,1,2\n1,2,3\n', '1,2\n0,3\n', '0,nan\n1,2\n'):
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'dos.csv'
                path.write_text(content)
                with self.assertRaises(ValueError):
                    load_dos_table(path)

    def test_missing_or_mismatched_spin_partner_does_not_replace_figure(self):
        with tempfile.TemporaryDirectory() as directory:
            up = Path(directory) / 'Si-DOS-GPAW-Result-DOS-Up.csv'
            down = up.with_name('Si-DOS-GPAW-Result-DOS-Down.csv')
            output = up.with_name('Si-DOS-GPAW-Graph-DOS.png')
            output.write_bytes(b'previous figure')
            np.savetxt(up, [[0, 1], [1, 2]], delimiter=',')
            with self.assertRaises(OSError):
                plot_dos_result(up)
            np.savetxt(down, [[0, 1], [2, 2]], delimiter=',')
            with self.assertRaisesRegex(ValueError, 'identical energy grids'):
                plot_dos_result(up)
            self.assertEqual(output.read_bytes(), b'previous figure')

    def test_dos_plot_is_written_beside_table_without_engine_loading(self):
        import matplotlib.pyplot as plt
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-DOS-QE-Result-DOS.csv'
            np.savetxt(path, [[-1, 1], [1, 2]], delimiter=',')
            before = set(plt.get_fignums())
            with patch('nanoworks.engine.load_engine_module') as engine:
                result = plot_saved_result(path)
            engine.assert_not_called()
            output = Path(result['figure_files']['DOS'][0])
            self.assertEqual(output.parent, path.parent)
            self.assertTrue(output.is_file())
            self.assertEqual(set(plt.get_fignums()), before)

    def test_common_command_still_accepts_optional_optical_summary(self):
        with patch('nanoworks.optical.plot_optical_summary', return_value={}) as optical:
            self.assertEqual(plot_saved_result('summary.json'), {})
            optical.assert_called_once_with(Path('summary.json').resolve())

    def test_optical_exports_one_dat_with_units_and_no_csv_duplicate(self):
        from nanoworks.optical import write_optical_table, OPTICAL_TABLE_HEADER
        from nanoworks.result_plotting import load_optical_table
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-OPTICAL-QE-Result-AllData.dat'
            data = np.array([[0, 3, 4, 2, 1, 0, .2], [1, 3, 4, 2, 1, 100, .2]])
            write_optical_table(path, data)
            self.assertFalse(path.with_suffix('.csv').exists())
            self.assertEqual(path.read_text().splitlines()[0].split(), OPTICAL_TABLE_HEADER.split())
            np.testing.assert_array_equal(load_optical_table(path), data)

    def test_new_dos_dat_header_is_skipped_and_reference_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-DOS-QE-Result-DOS.dat'
            path.write_text('Energy(eV) DOS(states/eV)\n-3 1\n2 4\n')
            np.testing.assert_array_equal(load_dos_table(path), [[-3, 1], [2, 4]])

    def test_standalone_table_plots_without_summary_and_preserves_values(self):
        from nanoworks.optical import OPTICAL_TABLE_HEADER
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'manually-edited.csv'
            data = np.array([[0, -3, 4, 2, 1, 0, .2], [1, 3, 4, 2, 1, 100, .2]])
            np.savetxt(path, data, delimiter=',', header=','.join(OPTICAL_TABLE_HEADER.split()), comments='')
            original = path.read_bytes()
            with patch('nanoworks.optical.generate_optical_figures', return_value=['figure.png']) as plot:
                with patch('nanoworks.optical.load_optical_summary') as summary:
                    result = plot_saved_result(path)
            summary.assert_not_called()
            np.testing.assert_array_equal(plot.call_args.args[0], data)
            self.assertEqual(plot.call_args.args[1], str(path.with_suffix('')))
            self.assertEqual(result['output_directory'], directory)
            self.assertEqual(path.read_bytes(), original)

    def test_invalid_optical_header_or_grid_does_not_plot(self):
        from nanoworks.optical import OPTICAL_TABLE_HEADER
        cases = ['kpoint energy\n0 1\n1 2\n',
                 OPTICAL_TABLE_HEADER + '\n0 3 4 2 1 0 .2\n0 3 4 2 1 100 .2\n']
        for content in cases:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'table.dat'
                path.write_text(content)
                with patch('nanoworks.optical.generate_optical_figures') as plot:
                    with self.assertRaises(ValueError):
                        plot_saved_result(path)
                plot.assert_not_called()


if __name__ == '__main__':
    unittest.main()
