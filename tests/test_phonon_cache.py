# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression cases for interrupted phonon force caches and MPI cache I/O."""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from nanoworks.phonon_cache import (
    collective_cache_call, force_signature, load_force_constants,
    load_verified_force, save_verified_force,
)


class TestPhononCache(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / 'force.npy'
        self.forces = np.array([[.1, 0, 0], [-.1, 0, 0]])

    def test_completed_displacement_survives_missing_global_completion(self):
        save_verified_force(self.path, self.forces, 'settings-A', 2)
        # No global force-constant or completion file exists yet.
        np.testing.assert_array_equal(load_verified_force(self.path, 'settings-A', 2), self.forces)
        self.assertIsNone(load_verified_force(self.path, 'settings-B', 2))

    def test_unverified_legacy_and_truncated_records_are_cache_misses(self):
        np.save(self.path, self.forces)
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))
        save_verified_force(self.path, self.forces, 'settings-A', 2)
        self.path.write_bytes(b'\x93NUMPY truncated')
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))
        Path(str(self.path) + '.json').write_text('{ interrupted metadata')
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))

    def test_altered_force_contents_and_wrong_shapes_are_rejected(self):
        save_verified_force(self.path, self.forces, 'settings-A', 2)
        np.save(self.path, self.forces * 2)
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))
        np.save(self.path, np.zeros((1, 3)))
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))
        np.save(self.path, np.full((2, 3), np.nan))
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))

    def test_invalid_new_forces_do_not_replace_previous_valid_record(self):
        save_verified_force(self.path, self.forces, 'settings-A', 2)
        with self.assertRaisesRegex(ValueError, 'finite natoms-by-3'):
            save_verified_force(self.path, np.full((2, 3), np.inf), 'settings-B', 2)
        np.testing.assert_array_equal(load_verified_force(self.path, 'settings-A', 2), self.forces)

    def test_interruption_between_array_and_metadata_never_certifies_new_forces(self):
        save_verified_force(self.path, self.forces, 'settings-A', 2)
        with patch('nanoworks.phonon_cache.write_json_atomic', side_effect=OSError('disk full')):
            with self.assertRaises(OSError):
                save_verified_force(self.path, self.forces * 2, 'settings-B', 2)
        self.assertIsNone(load_verified_force(self.path, 'settings-A', 2))
        self.assertIsNone(load_verified_force(self.path, 'settings-B', 2))

    def test_geometry_magnetism_and_settings_change_displacement_signature(self):
        cell = SimpleNamespace(symbols=['Ni', 'Ni'], cell=np.eye(3),
            positions=np.array([[0., 0, 0], [.5, .5, .5]]), magnetic_moments=np.array([2., -2.]))
        original = force_signature({'cutoff': 500, 'hubbard_u': 6}, cell)
        self.assertNotEqual(original, force_signature({'cutoff': 600, 'hubbard_u': 6}, cell))
        cell.positions[0, 0] += .01
        self.assertNotEqual(original, force_signature({'cutoff': 500, 'hubbard_u': 6}, cell))
        cell.positions[0, 0] = 0
        cell.magnetic_moments[:] = 0
        self.assertNotEqual(original, force_signature({'cutoff': 500, 'hubbard_u': 6}, cell))

    def test_force_constant_cache_accepts_compact_and_full_shapes_only(self):
        for shape in ((2, 4, 3, 3), (4, 4, 3, 3)):
            np.save(self.path, np.zeros(shape))
            self.assertEqual(load_force_constants(self.path, 2, 4).shape, shape)
        np.save(self.path, np.zeros((3, 4, 3, 3)))
        self.assertIsNone(load_force_constants(self.path, 2, 4))
        np.save(self.path, np.full((2, 4, 3, 3), np.nan))
        self.assertIsNone(load_force_constants(self.path, 2, 4))

    def test_root_io_failure_is_shared_with_all_ranks(self):
        callback = Mock(side_effect=OSError('disk full'))
        with patch('ase.parallel.world', SimpleNamespace(rank=0)):
            with patch('ase.parallel.broadcast', side_effect=lambda value, **kwargs: value) as broadcast:
                with self.assertRaisesRegex(RuntimeError, 'disk full'):
                    collective_cache_call(callback)
        callback.assert_called_once()
        self.assertEqual(broadcast.call_args.args[0][0], 'OSError: disk full')

    def test_nonroot_rank_receives_force_result_without_reading_files(self):
        callback = Mock()
        with patch('ase.parallel.world', SimpleNamespace(rank=1)):
            with patch('ase.parallel.broadcast', return_value=(None, self.forces)):
                result = collective_cache_call(callback)
        callback.assert_not_called()
        np.testing.assert_array_equal(result, self.forces)
