# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import unittest

import numpy as np
from ase import Atoms

from nanoworks.elasticity import (
    analyze_elastic_stability,
    calculate_2d_elastic_properties,
    normalize_elastic_dimensionality,
    resolve_elastic_dimensionality,
)


class TestElasticity(unittest.TestCase):

    def test_dimensionality_values_are_normalized(self):
        self.assertEqual(normalize_elastic_dimensionality('Auto'), 'auto')
        self.assertEqual(normalize_elastic_dimensionality('2d'), '2D')
        self.assertEqual(normalize_elastic_dimensionality('3D'), '3D')

    def test_invalid_dimensionality_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Elastic_dimensionality'):
            normalize_elastic_dimensionality('slab')

    def test_auto_detects_z_normal_slab(self):
        atoms = Atoms(
            'C2',
            scaled_positions=[(0.0, 0.0, 0.5), (0.5, 0.5, 0.5)],
            cell=[2.5, 2.5, 20.0],
            pbc=True,
        )

        result = resolve_elastic_dimensionality(atoms)

        self.assertEqual(result['resolved'], '2D')
        self.assertEqual(result['normal_axis'], 'z')
        self.assertAlmostEqual(result['vacuum_gap_angstrom'], 20.0)

    def test_auto_keeps_bulk_as_3d(self):
        atoms = Atoms(
            'Si2',
            scaled_positions=[(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)],
            cell=[5.43, 5.43, 5.43],
            pbc=True,
        )

        result = resolve_elastic_dimensionality(atoms)

        self.assertEqual(result['resolved'], '3D')

    def test_auto_does_not_treat_one_atom_cubic_cell_as_slab(self):
        atoms = Atoms(
            'Al',
            scaled_positions=[(0.0, 0.0, 0.0)],
            cell=[6.0, 6.0, 6.0],
            pbc=True,
        )

        result = resolve_elastic_dimensionality(atoms)

        self.assertEqual(result['resolved'], '3D')

    def test_explicit_dimensionality_overrides_detection(self):
        atoms = Atoms(
            'Si',
            scaled_positions=[(0.0, 0.0, 0.0)],
            cell=[5.43, 5.43, 5.43],
            pbc=True,
        )

        result = resolve_elastic_dimensionality(
            atoms,
            dimensionality='2D',
        )

        self.assertEqual(result['resolved'], '2D')
        self.assertFalse(result['automatic'])

    def test_2d_properties_remove_supercell_vacuum_scaling(self):
        tensor = np.zeros((6, 6))
        tensor[0, 0] = 100.0
        tensor[1, 1] = 100.0
        tensor[0, 1] = tensor[1, 0] = 20.0
        tensor[5, 5] = 40.0

        result = calculate_2d_elastic_properties(
            tensor,
            normal_cell_length_angstrom=20.0,
        )

        np.testing.assert_allclose(
            result['stiffness_n_per_m'],
            [
                [200.0, 40.0, 0.0],
                [40.0, 200.0, 0.0],
                [0.0, 0.0, 80.0],
            ],
        )
        self.assertAlmostEqual(
            result['young_modulus_first_n_per_m'],
            192.0,
        )
        self.assertAlmostEqual(
            result['young_modulus_second_n_per_m'],
            192.0,
        )
        self.assertAlmostEqual(result['shear_modulus_n_per_m'], 80.0)
        self.assertAlmostEqual(result['poisson_ratio_first_second'], 0.2)
        self.assertAlmostEqual(result['poisson_ratio_second_first'], 0.2)

    def test_stable_3d_tensor_is_positive_definite(self):
        tensor = np.array([
            [200, 120, 120, 0, 0, 0],
            [120, 200, 120, 0, 0, 0],
            [120, 120, 200, 0, 0, 0],
            [0, 0, 0, 80, 0, 0],
            [0, 0, 0, 0, 80, 0],
            [0, 0, 0, 0, 0, 80],
        ], dtype=float)

        result = analyze_elastic_stability(tensor)

        self.assertTrue(result['tensor_symmetric'])
        self.assertTrue(result['mechanically_stable'])
        self.assertAlmostEqual(result['minimum_eigenvalue'], 80.0)
        self.assertEqual(result['units'], 'GPa')

    def test_negative_stiffness_eigenvalue_is_unstable(self):
        tensor = np.eye(6) * 50.0
        tensor[5, 5] = -1.0

        result = analyze_elastic_stability(tensor)

        self.assertFalse(result['mechanically_stable'])
        self.assertAlmostEqual(result['minimum_eigenvalue'], -1.0)

    def test_tensor_asymmetry_is_reported_separately(self):
        tensor = np.eye(6) * 100.0
        tensor[0, 1] = 10.0
        tensor[1, 0] = 8.0

        result = analyze_elastic_stability(tensor)

        self.assertFalse(result['tensor_symmetric'])
        self.assertAlmostEqual(
            result['maximum_tensor_asymmetry_gpa'],
            2.0,
        )
        self.assertTrue(result['mechanically_stable'])

    def test_2d_stability_uses_in_plane_stiffness(self):
        tensor = np.eye(6)
        stiffness = np.array([
            [200.0, 40.0, 0.0],
            [40.0, 200.0, 0.0],
            [0.0, 0.0, 80.0],
        ])

        result = analyze_elastic_stability(
            tensor,
            dimensionality='2D',
            stiffness_2d_n_per_m=stiffness,
        )

        self.assertTrue(result['mechanically_stable'])
        self.assertEqual(result['units'], 'N/m')
        self.assertEqual(len(result['eigenvalues']), 3)


if __name__ == '__main__':
    unittest.main()
