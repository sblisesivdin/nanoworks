# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from nanoworks.result_plotting import load_dos_table, plot_dos_result, plot_saved_result


class TestSavedResultPlotting(unittest.TestCase):
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
