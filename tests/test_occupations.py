# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import unittest

from nanoworks.occupations import (
    resolve_engine_occupation,
    resolve_gpaw_occupation,
    resolve_qe_occupation_input,
    validate_occupation_settings,
)


class TestPortableOccupations(unittest.TestCase):

    def test_aliases_are_normalized(self):
        self.assertEqual(
            validate_occupation_settings('cold', 0.2),
            {'scheme': 'marzari-vanderbilt', 'width': 0.2},
        )

    def test_fixed_ignores_width(self):
        self.assertEqual(
            validate_occupation_settings('fixed', -1.0),
            {'scheme': 'fixed', 'width': None},
        )

    def test_smearing_requires_positive_finite_width(self):
        for width in (None, 0.0, -0.1, float('inf'), True):
            with self.subTest(width=width):
                with self.assertRaisesRegex(ValueError, 'Smearing_width'):
                    validate_occupation_settings('fermi-dirac', width)

    def test_unknown_scheme_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Occupation_scheme'):
            validate_occupation_settings('tetrahedra', 0.05)

    def test_gpaw_mapping_uses_fixed_uniform(self):
        self.assertEqual(
            resolve_gpaw_occupation(scheme='fixed', width=None),
            {'name': 'fixed-uniform'},
        )

    def test_qe_mapping_uses_fixed(self):
        self.assertEqual(
            resolve_qe_occupation_input(scheme='fixed', width=None),
            'fixed',
        )

    def test_smearing_mapping_is_engine_neutral(self):
        expected = {'name': 'methfessel-paxton', 'width': 0.1}
        self.assertEqual(
            resolve_engine_occupation(
                'GPAW',
                scheme='methfessel-paxton',
                width=0.1,
            ),
            expected,
        )
        self.assertEqual(
            resolve_engine_occupation(
                'QE',
                scheme='methfessel-paxton',
                width=0.1,
            ),
            expected,
        )

    def test_unknown_engine_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported DFT engine'):
            resolve_engine_occupation(
                'other',
                scheme='fixed',
                width=None,
            )


if __name__ == '__main__':
    unittest.main()
