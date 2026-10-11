# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import tempfile
import unittest
from pathlib import Path

import numpy as np

from nanoworks.band_results import BAND_PATH_COLUMNS, load_band_path_table, write_band_path_table


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
