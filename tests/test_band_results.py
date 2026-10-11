# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from nanoworks.band_results import (
    BAND_PATH_COLUMNS, load_band_path_table, write_band_path_table,
    write_band_energy_tables,
)
from nanoworks.result_plotting import load_band_table


class TestBandEnergyTables(unittest.TestCase):
    def test_long_and_wide_tables_preserve_every_spin_energy(self):
        energies = np.array([[[-2., 1.], [-1., 2.], [0., 3.]],
                             [[-3., 4.], [-2., 5.], [-1., 6.]]])
        with tempfile.TemporaryDirectory() as directory:
            paths = write_band_energy_tables(Path(directory) / 'Ni-BAND-GPAW',
                                             energies, ['Down', 'Up'])
            self.assertEqual(len(paths), 4)
            for path, expected in zip(paths, [energies[0], energies[0], energies[1], energies[1]]):
                indices, values = load_band_table(path)
                np.testing.assert_array_equal(indices, np.arange(3))
                np.testing.assert_array_equal(values, expected)

    def test_staging_failure_preserves_all_existing_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'Si-BAND-GPAW'
            paths = write_band_energy_tables(prefix, np.ones((1, 3, 2)))
            previous = {path: Path(path).read_bytes() for path in paths}
            original = Path.open

            def fail_wide(path, *args, **kwargs):
                if str(path).endswith('-XYYY.dat.tmp'):
                    raise OSError('injected write failure')
                return original(path, *args, **kwargs)

            with patch.object(Path, 'open', fail_wide):
                with self.assertRaises(OSError):
                    write_band_energy_tables(prefix, np.zeros((1, 3, 2)))
            self.assertEqual(previous, {path: Path(path).read_bytes() for path in paths})
            self.assertFalse(list(Path(directory).glob('*.tmp')))

    def test_invalid_energies_and_labels_do_not_replace_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'Si-BAND-GPAW'
            paths = write_band_energy_tables(prefix, np.ones((1, 3, 2)))
            previous = {path: Path(path).read_bytes() for path in paths}
            for energies, labels in [(np.full((1, 3, 2), np.nan), None),
                                     (np.ones((2, 3, 2)), ['Up', 'Up']),
                                     (np.ones((1, 1, 2)), None)]:
                with self.assertRaises(ValueError):
                    write_band_energy_tables(prefix, energies, labels)
            self.assertEqual(previous, {path: Path(path).read_bytes() for path in paths})


class TestBandPathTables(unittest.TestCase):
    def test_labels_and_disconnected_segments_round_trip_as_four_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'Si-BAND-QE-Result-KPath.dat'
            write_band_path_table(path, [0, .5, .5, 1], [0, .5, .5, 1], ['G', 'X', 'W', 'K'])
            rows = path.read_text().splitlines()
            self.assertEqual(rows[0].split(), BAND_PATH_COLUMNS)
            self.assertTrue(all(len(row.split()) == 4 for row in rows))
            result = load_band_path_table(path, np.arange(4))
            np.testing.assert_array_equal(result['distances'], [0, .5, .5, 1])
            np.testing.assert_array_equal(result['breaks'], [False, False, True, False])
            self.assertEqual(result['labels'], ['G', 'X,W', 'K'])

    def test_invalid_axis_or_unsampled_ticks_do_not_replace_previous_table(self):
        cases = [([0, np.nan, 1], [0, 1], ['G', 'X']),
                 ([0, 1, .5], [0, .5], ['G', 'X']),
                 ([0, .5, 1], [.25], ['X'])]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'path.dat'
            for distances, ticks, labels in cases:
                with self.subTest(distances=distances):
                    path.write_bytes(b'previous axis')
                    with self.assertRaises(ValueError):
                        write_band_path_table(path, distances, ticks, labels)
                    self.assertEqual(path.read_bytes(), b'previous axis')

    def test_path_grid_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'path.dat'
            write_band_path_table(path, [0, .5, 1], [0, 1], ['G', 'X'])
            with self.assertRaisesRegex(ValueError, 'do not match'):
                load_band_path_table(path, np.arange(2))


if __name__ == '__main__':
    unittest.main()
