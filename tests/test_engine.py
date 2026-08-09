import unittest

from nanoworks.engine import (
    normalize_engine_name,
    resolve_stage_kpoint_settings,
    resolve_stage_occupation,
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


if __name__ == '__main__':
    unittest.main()
