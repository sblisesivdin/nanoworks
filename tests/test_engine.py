import unittest
from unittest.mock import patch
from ase import Atoms

from nanoworks.engine import (
    normalize_engine_name,
    resolve_stage_kpoint_settings,
    resolve_stage_occupation,
    load_engine_module,
    resolve_initial_magnetic_moments,
)


class TestEngine(unittest.TestCase):

    def test_engine_name_is_normalized(self):
        self.assertEqual(
            normalize_engine_name('gpaw'),
            'GPAW',
        )
        self.assertEqual(
            normalize_engine_name(' GPAW '),
            'GPAW',
        )
    
    def test_stage_kpoints_fall_back_to_ground(self):
        density, size, gamma = resolve_stage_kpoint_settings(
            stage_density=None,
            stage_size=(None, None, None),
            stage_gamma=None,
            ground_density=4.0,
            ground_size=(6, 6, 4),
            ground_gamma=True,
        )

        self.assertEqual(density, 4.0)
        self.assertEqual(size, (6, 6, 4))
        self.assertTrue(gamma)
    
    def test_stage_mesh_overrides_ground_density(self):
        density, size, gamma = resolve_stage_kpoint_settings(
            stage_density=None,
            stage_size=(16, 16, 8),
            stage_gamma=False,
            ground_density=4.0,
            ground_size=(6, 6, 4),
            ground_gamma=True,
        )

        self.assertIsNone(density)
        self.assertEqual(size, (16, 16, 8))
        self.assertFalse(gamma)
    
    def test_stage_occupation_falls_back_to_ground(self):
        ground = {
            'name': 'fermi-dirac',
            'width': 0.05,
        }

        self.assertIs(
            resolve_stage_occupation(None, ground),
            ground,
        )

        dos = {
            'name': 'tetrahedron-method',
        }

        self.assertIs(
            resolve_stage_occupation(dos, ground),
            dos,
        )
    
    @patch('nanoworks.engine.import_module')
    def test_engine_module_is_loaded_lazily(self, import_module):
        expected = object()
        import_module.return_value = expected

        result = load_engine_module('gpaw')

        self.assertIs(result, expected)
        import_module.assert_called_once_with(
            'nanoworks.engine.gpaw'
        )
    
    @patch('nanoworks.engine.import_module')
    def test_qe_engine_module_is_loaded_lazily(self, import_module):
        expected = object()
        import_module.return_value = expected

        result = load_engine_module('qe')

        self.assertIs(result, expected)
        import_module.assert_called_once_with(
            'nanoworks.engine.qe'
        )


    def test_unknown_engine_is_rejected(self):
        with self.assertRaises(ValueError):
            load_engine_module('unknown')

    def test_initial_magnetic_moments_support_scalar(self):
        atoms = Atoms(
            'FeO',
        )

        moments = resolve_initial_magnetic_moments(
            atoms,
            magmom_per_atom=2.0,
        )

        self.assertEqual(
            moments,
            [2.0, 2.0],
        )

    def test_initial_magnetic_moments_support_species_mapping(self):
        atoms = Atoms(
            'Fe2O3',
        )

        moments = resolve_initial_magnetic_moments(
            atoms,
            magmom_per_atom={
                'Fe': 4.0,
                'O': 0.0,
            },
        )

        self.assertEqual(
            moments,
            [
                4.0,
                4.0,
                0.0,
                0.0,
                0.0,
            ],
        )

    def test_initial_magnetic_moments_default_missing_species_to_zero(self):
        atoms = Atoms(
            'FeO',
        )

        moments = resolve_initial_magnetic_moments(
            atoms,
            magmom_per_atom={
                'Fe': 2.2,
            },
        )

        self.assertEqual(
            moments,
            [2.2, 0.0],
        )

    def test_initial_magnetic_moments_support_per_atom_sequence(self):
        atoms = Atoms(
            'Fe2O3',
        )

        moments = resolve_initial_magnetic_moments(
            atoms,
            magmom_per_atom=[
                4.0,
                -4.0,
                0.0,
                0.0,
                0.0,
            ],
        )

        self.assertEqual(
            moments,
            [
                4.0,
                -4.0,
                0.0,
                0.0,
                0.0,
            ],
        )

    def test_single_atom_moment_overrides_resolved_moment(self):
        atoms = Atoms(
            'Fe2O3',
        )

        moments = resolve_initial_magnetic_moments(
            atoms,
            magmom_per_atom={
                'Fe': 4.0,
                'O': 0.0,
            },
            magmom_single_atom=[
                1,
                -4.0,
            ],
        )

        self.assertEqual(
            moments,
            [
                4.0,
                -4.0,
                0.0,
                0.0,
                0.0,
            ],
        )

if __name__ == '__main__':
    unittest.main()
