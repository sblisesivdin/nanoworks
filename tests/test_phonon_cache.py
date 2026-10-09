# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression cases for interrupted phonon force caches and MPI cache I/O."""

import tempfile
import json
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from nanoworks.phonon_cache import (
    collective_cache_call, force_signature, load_force_constants,
    load_verified_force, save_verified_force,
    load_verified_force_constants, save_verified_force_constants,
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

    def test_verified_force_constants_accept_full_and_compact_records(self):
        metadata = self.path.with_suffix('.json')
        settings = {'cutoff': 500, 'masses': [60, 60]}
        for shape in ((2, 4, 3, 3), (4, 4, 3, 3)):
            with self.subTest(shape=shape):
                constants = np.arange(np.prod(shape), dtype=float).reshape(shape) / 100
                save_verified_force_constants(self.path, metadata, constants, settings, 2, 4)
                np.testing.assert_array_equal(load_verified_force_constants(
                    self.path, metadata, settings, 2, 4), constants)
                self.assertIsNone(load_verified_force_constants(
                    self.path, metadata, dict(settings, cutoff=600), 2, 4))
                self.assertIsNone(load_verified_force_constants(
                    self.path, metadata, dict(settings, masses=[62, 62]), 2, 4))

    def test_altered_finite_constants_do_not_invalidate_displacement_forces(self):
        metadata = self.path.with_suffix('.json')
        constants = np.ones((2, 4, 3, 3))
        force_path = self.path.with_name('displacement.npy')
        save_verified_force(force_path, self.forces, 'settings-A', 2)
        save_verified_force_constants(self.path, metadata, constants, {'cutoff': 500}, 2, 4)
        constants[0, 0, 0, 0] += .1
        np.save(self.path, constants)
        self.assertIsNone(load_verified_force_constants(
            self.path, metadata, {'cutoff': 500}, 2, 4))
        np.testing.assert_array_equal(load_verified_force(force_path, 'settings-A', 2), self.forces)

    def test_unhashed_or_incomplete_force_constant_completion_is_rejected(self):
        metadata = self.path.with_suffix('.json')
        constants = np.ones((2, 4, 3, 3))
        settings = {'cutoff': 500}
        for record in (settings, {}, {'schema': 1, 'units': 'eV/Angstrom^2', 'settings': settings}):
            np.save(self.path, constants)
            metadata.write_text(json.dumps(record))
            self.assertIsNone(load_verified_force_constants(self.path, metadata, settings, 2, 4))
        metadata.write_text('{ interrupted metadata')
        self.assertIsNone(load_verified_force_constants(self.path, metadata, settings, 2, 4))

    def test_force_constant_write_failure_never_certifies_changed_array(self):
        metadata = self.path.with_suffix('.json')
        constants = np.ones((2, 4, 3, 3))
        settings = {'cutoff': 500}
        save_verified_force_constants(self.path, metadata, constants, settings, 2, 4)
        with patch('nanoworks.phonon_cache.write_json_atomic', side_effect=OSError('disk full')):
            with self.assertRaises(OSError):
                save_verified_force_constants(self.path, metadata, constants * 2, settings, 2, 4)
        self.assertIsNone(load_verified_force_constants(self.path, metadata, settings, 2, 4))

    def test_invalid_new_constants_preserve_previous_verified_cache(self):
        metadata = self.path.with_suffix('.json')
        constants = np.ones((2, 4, 3, 3))
        settings = {'cutoff': 500}
        save_verified_force_constants(self.path, metadata, constants, settings, 2, 4)
        for invalid in (np.zeros((3, 4, 3, 3)), np.full((2, 4, 3, 3), np.nan)):
            with self.subTest(shape=invalid.shape), self.assertRaisesRegex(ValueError, 'compact/full'):
                save_verified_force_constants(self.path, metadata, invalid, settings, 2, 4)
        np.testing.assert_array_equal(load_verified_force_constants(
            self.path, metadata, settings, 2, 4), constants)

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
