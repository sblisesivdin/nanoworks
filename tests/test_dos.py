# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import unittest

from nanoworks.dos import resolve_dos_settings, validate_dos_settings


class TestPortableDOSSettings(unittest.TestCase):

    def test_smearing_requires_positive_finite_width(self):
        for width in (None, 0.0, -0.1, float('inf'), True):
            with self.subTest(width=width):
                with self.assertRaisesRegex(ValueError, 'DOS_width'):
                    validate_dos_settings('smearing', width)

    def test_tetrahedron_normalizes_width_to_zero(self):
        self.assertEqual(
            validate_dos_settings('tetrahedra', 0.2),
            {'integration': 'tetrahedron', 'width_ev': 0.0},
        )

    def test_unknown_integration_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'DOS_integration'):
            validate_dos_settings('optimized-tetrahedron', 0.1)

    def test_gpaw_uses_ground_occupation_and_portable_width(self):
        ground = {'name': 'fermi-dirac', 'width': 0.05}
        settings = resolve_dos_settings(
            'GPAW',
            'smearing',
            0.15,
            ground,
        )

        self.assertIs(settings['electronic_occupation'], ground)
        self.assertEqual(settings['width_ev'], 0.15)
        self.assertIsNone(settings['bz_sum'])

    def test_qe_tetrahedron_sets_electronic_and_dos_integration(self):
        settings = resolve_dos_settings(
            'QE',
            'tetrahedron',
            0.1,
            {'name': 'fermi-dirac', 'width': 0.05},
        )

        self.assertEqual(settings['electronic_occupation'], 'tetrahedra')
        self.assertEqual(settings['bz_sum'], 'tetrahedra')
        self.assertIsNone(settings['degauss_ev'])

    def test_qe_smearing_uses_explicit_gaussian_broadening(self):
        ground = {'name': 'fermi-dirac', 'width': 0.05}
        settings = resolve_dos_settings(
            'QE',
            'smearing',
            0.2,
            ground,
        )

        self.assertIs(settings['electronic_occupation'], ground)
        self.assertEqual(settings['bz_sum'], 'smearing')
        self.assertEqual(settings['degauss_ev'], 0.2)
        self.assertEqual(settings['ngauss'], 0)

    def test_unknown_engine_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported DFT engine'):
            resolve_dos_settings('other', 'smearing', 0.1, None)


if __name__ == '__main__':
    unittest.main()
