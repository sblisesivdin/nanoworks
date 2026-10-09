# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Reject invalid phonon parameters before force calculations or postprocessing."""

import unittest
from types import SimpleNamespace

import numpy as np

from nanoworks.phonon_settings import validate_phonon_settings


class TestPhononSettings(unittest.TestCase):
    def config(self, **changes):
        values = dict(Engine='GPAW', Hubbard_U=None, Phonon_supercell=[2, 2, 1],
            Phonon_qpts_x=8, Phonon_qpts_y=8, Phonon_qpts_z=1, Phonon_npoints=61,
            Phonon_kpts_x=None, Phonon_kpts_y=None, Phonon_kpts_z=1,
            Phonon_PW_cutoff=None, Phonon_displacement=.01, Phonon_thermal_calc=True,
            Phonon_T_min=0., Phonon_T_max=300., Phonon_T_step=10., bulk_configuration=[0, 1])
        values.update(changes)
        return SimpleNamespace(**values)

    def test_vector_normalization_and_supercell_atom_count(self):
        values, report = validate_phonon_settings(self.config())
        np.testing.assert_array_equal(values['Phonon_supercell'], np.diag([2, 2, 1]))
        self.assertEqual(report['supercell_atoms'], 8)
        self.assertEqual(report['dos_mesh'], [8, 8, 1])

    def test_nondiagonal_matrix_is_supported_only_by_finite_displacements(self):
        matrix = [[2, 1, 0], [0, 2, 0], [0, 0, 1]]
        for engine, hubbard in (('GPAW', None), ('QE', {'Ni-3d': 6})):
            _, report = validate_phonon_settings(self.config(
                Engine=engine, Hubbard_U=hubbard, Phonon_supercell=matrix))
            self.assertEqual(report['supercell_multiplier'], 4)
        with self.assertRaisesRegex(ValueError, 'Native QE DFPT'):
            validate_phonon_settings(self.config(Engine='QE', Phonon_supercell=matrix))

    def test_counts_never_truncate_fractional_values_or_accept_booleans(self):
        for name in ('Phonon_kpts_x', 'Phonon_qpts_z', 'Phonon_npoints'):
            for value in (2.5, True, 0, '2'):
                with self.subTest(name=name, value=value), self.assertRaisesRegex(ValueError, name):
                    validate_phonon_settings(self.config(**{name: value}))

    def test_invalid_supercell_rejected(self):
        for matrix in ([2, 0, 1], [2, -1, 1], [[1, 0, 0], [0, 0, 0], [0, 0, 1]],
                       [[1, 0, 0], [0, 1, 0], [0, 0, -1]], [2.5, 2, 1], [True, 2, 1]):
            with self.subTest(matrix=matrix), self.assertRaisesRegex(ValueError, 'Phonon_supercell'):
                validate_phonon_settings(self.config(Phonon_supercell=matrix))

    def test_invalid_thermal_ranges_displacements_and_cutoffs(self):
        for changes in ({'Phonon_T_step': 0}, {'Phonon_T_min': -1},
                        {'Phonon_T_max': -1}, {'Phonon_T_max': float('inf')},
                        {'Phonon_T_min': 400, 'Phonon_T_max': 300},
                        {'Phonon_displacement': float('nan')},
                        {'Phonon_displacement': 0}, {'Phonon_PW_cutoff': -1}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                validate_phonon_settings(self.config(**changes))

    def test_mass_validation_is_deferred_until_final_structure_is_known(self):
        from ase import Atoms
        atoms = Atoms('Ni', cell=[4, 4, 4], pbc=True, masses=[0])
        config = self.config(bulk_configuration=atoms)
        _, report = validate_phonon_settings(config, validate_structure=False)
        self.assertNotIn('atomic_masses_amu', report)
        with self.assertRaisesRegex(ValueError, 'atomic masses'):
            validate_phonon_settings(config)
        # Native QE DFPT does not use the ASE/Phonopy finite-displacement mass path.
        _, report = validate_phonon_settings(self.config(Engine='QE', bulk_configuration=atoms))
        self.assertNotIn('atomic_masses_amu', report)

    def test_unused_thermal_and_native_dfpt_displacement_do_not_block(self):
        _, report = validate_phonon_settings(self.config(Engine='QE', Phonon_thermal_calc=False,
            Phonon_T_step=0, Phonon_displacement=-1))
        self.assertEqual(report['dfpt_qpoint_grid'], [2, 2, 1])
        self.assertNotIn('supercell_atoms', report)
        self.assertNotIn('temperature_range_kelvin', report)


class TestPhononAtomicMasses(unittest.TestCase):
    def test_invalid_mass_arrays_are_rejected(self):
        from nanoworks.phonon_settings import validate_atomic_masses
        for masses in ([0, 60], [-1, 60], [np.nan, 60], [np.inf, 60], [60], [[60, 60]]):
            with self.subTest(masses=masses), self.assertRaisesRegex(ValueError, 'atomic masses'):
                validate_atomic_masses(masses, 2)


class TestEditableAnalysisSettings(unittest.TestCase):
    def test_mesh_rejects_fractional_boolean_and_wrong_dimensions(self):
        from nanoworks.phonon_settings import validate_dos_mesh
        self.assertEqual(validate_dos_mesh([2, 3, 4]), [2, 3, 4])
        for mesh in (None, [], [2, 2], [2, 2, 0], [2, 2.5, 2], [True, 2, 2], '222'):
            with self.subTest(mesh=mesh), self.assertRaises(ValueError):
                validate_dos_mesh(mesh)

    def test_thermal_range_rejects_disabled_like_strings_and_invalid_ranges(self):
        from nanoworks.phonon_settings import validate_temperature_range
        self.assertIsNone(validate_temperature_range(None))
        self.assertEqual(validate_temperature_range([0, 300, 100]), [0., 300., 100.])
        for temperature in ([], [0, 300], [0, 300, 0], [300, 0, 10],
                            [0, float('nan'), 10], [False, 300, 10], '123'):
            with self.subTest(temperature=temperature), self.assertRaises(ValueError):
                validate_temperature_range(temperature)
