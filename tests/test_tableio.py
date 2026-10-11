# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import io
import unittest
from urllib.parse import unquote

from nanoworks.tableio import DatWriter, table_header, csv_line_to_dat


class TestResultTableContract(unittest.TestCase):
    def test_header_names_do_not_create_extra_columns(self):
        header = table_header(['Energy', 'Projected DOS(states/eV)', 'Step'])
        self.assertEqual(header.split(), ['Energy(eV)', 'Projected_DOS(states/eV)', 'Step(1)'])

    def test_missing_and_text_fields_keep_column_positions(self):
        stream = io.StringIO()
        writer = DatWriter(stream)
        writer.writerow(['Energy', 'label', 'weight'])
        writer.writerow([-1.25, 'W d, j=5/2', None])
        header, row = stream.getvalue().splitlines()
        self.assertEqual(len(header.split()), len(row.split()))
        self.assertEqual(unquote(row.split()[1]), 'W d, j=5/2')
        self.assertEqual(row.split()[2], 'nan')
        with self.assertRaises(ValueError):
            writer.writerow([1, 2])

    def test_duplicate_and_empty_headers_are_rejected(self):
        for labels in (['x', 'x'], ['x', '']):
            with self.subTest(labels=labels), self.assertRaises(ValueError):
                table_header(labels)

    def test_md_formatted_values_keep_precision_and_missing_fields(self):
        self.assertEqual(csv_line_to_dat('1,1.200000,,2.0\n'), '1 1.200000 nan 2.0\n')


if __name__ == '__main__':
    unittest.main()
