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

    def test_common_command_dispatches_optical_and_rejects_band_dat(self):
        with patch('nanoworks.optical.plot_optical_summary', return_value={}) as optical:
            self.assertEqual(plot_saved_result('summary.json'), {})
            optical.assert_called_once_with(Path('summary.json').resolve())
        with self.assertRaisesRegex(ValueError, 'Supported saved results'):
            plot_saved_result('Si-BAND-QE-Result-Band.dat')


if __name__ == '__main__':
    unittest.main()
