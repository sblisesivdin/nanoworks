"""Tests for engine-neutral Hubbard-U settings."""

import unittest

from nanoworks.hubbard import (
    normalize_hubbard_u,
    resolve_gpaw_hubbard,
)


class TestHubbardSettings(unittest.TestCase):
    def test_explicit_manifolds_and_energies_are_normalized(self):
        self.assertEqual(
            normalize_hubbard_u({
                ' O-02p ': '7',
                'Zn-3d': 10,
            }),
            {
                'O-2p': 7.0,
                'Zn-3d': 10.0,
            },
        )

    def test_non_dictionary_is_rejected(self):
        with self.assertRaisesRegex(TypeError, 'dictionary'):
            normalize_hubbard_u([('O-2p', 7.0)])

    def test_implicit_orbital_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'explicit'):
            normalize_hubbard_u({'O-p': 7.0})

    def test_non_numeric_energy_is_rejected(self):
        for value in (True, 'large'):
            with self.subTest(value=value):
                with self.assertRaisesRegex(TypeError, 'numeric'):
                    normalize_hubbard_u({'O-2p': value})

    def test_nonpositive_or_nonfinite_energy_is_rejected(self):
        for value in (0.0, -1.0, float('inf'), float('nan')):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, 'finite'):
                    normalize_hubbard_u({'O-2p': value})

    def test_gpaw_translation_drops_principal_number(self):
        self.assertEqual(
            resolve_gpaw_hubbard({
                'O-2p': 7.0,
                'Zn-3d': 10.0,
            }),
            {
                'O': ':p,7.0',
                'Zn': ':d,10.0',
            },
        )

    def test_gpaw_rejects_multiple_corrections_per_element(self):
        with self.assertRaisesRegex(ValueError, 'one Hubbard_U'):
            resolve_gpaw_hubbard({
                'Fe-3d': 4.0,
                'Fe-4s': 1.0,
            })


if __name__ == '__main__':
    unittest.main()
