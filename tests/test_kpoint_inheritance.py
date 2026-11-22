# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from ase import Atoms

from nanoworks.engine import resolve_stage_kpoint_settings
from nanoworks.dftsolve import DFTConfig, dftsolve as DFTSolver, prepare_qe_dry_run


class TestKpointInheritance(unittest.TestCase):
    def setUp(self):
        self.atoms = Atoms('Si2', cell=[4, 6, 20], pbc=True)

    def resolve(self, **changes):
        values = dict(stage_density=None, stage_size=(None, None, 2), stage_gamma=None,
            ground_density=3., ground_size=(17, 18, 19), ground_gamma=True, atoms=self.atoms)
        values.update(changes)
        return resolve_stage_kpoint_settings(**values)

    def test_partial_mesh_inherits_the_density_mesh_instead_of_unused_axis_defaults(self):
        self.assertEqual(self.resolve(), (None, (5, 4, 2), True))
        self.assertEqual(self.resolve(stage_size=(8, None, None)), (None, (8, 4, 1), True))
        self.assertEqual(self.resolve(stage_gamma=False), (None, (5, 4, 2), False))

    def test_inheritance_is_resolved_again_for_the_final_cell(self):
        self.atoms.set_cell([8, 6, 20])
        self.assertEqual(self.resolve(), (None, (3, 4, 2), True))
        self.atoms.pbc = [True, False, False]
        self.assertEqual(self.resolve(), (None, (3, 1, 2), True))

    def test_explicit_density_full_mesh_and_plain_ground_mesh_keep_their_priority(self):
        self.assertEqual(self.resolve(stage_density=7., atoms=None), (7., (17, 18, 19), True))
        self.assertEqual(self.resolve(stage_size=(9, 8, 1), atoms=None), (None, (9, 8, 1), True))
        self.assertEqual(self.resolve(ground_density=None, atoms=None), (None, (17, 18, 2), True))
        self.assertEqual(self.resolve(stage_size=(None, None, None), atoms=None),
                         (3., (17, 18, 19), True))

    def test_partial_density_mesh_rejects_missing_structure_or_invalid_density(self):
        with self.assertRaisesRegex(ValueError, 'require the structure'):
            self.resolve(atoms=None)
        for density in (0., -1., float('nan'), float('inf')):
            with self.subTest(density=density), self.assertRaisesRegex(ValueError, 'Ground_kpts_density'):
                self.resolve(ground_density=density)

    def test_qe_dry_run_uses_same_partial_mesh_for_dos_optics_and_elastic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / 'Si.upf').write_text(
                '<UPF><PP_HEADER pseudo_type="NC" z_valence="4"/></UPF>', encoding='utf-8')
            config = DFTConfig(Engine='QE', Ground_calc=True, DOS_calc=True,
                Optical_calc=True, Elastic_calc=True, bulk_configuration=self.atoms,
                Ground_kpts_density=3., DOS_kpts_z=1, Opt_kpts_z=1, Elastic_kpts_z=1,
                Pseudo_dir=str(root), Pseudopotentials={'Si': 'Si.upf'}, DOS_gamma=False)
            plan = prepare_qe_dry_run(config, root / 'sample')
            jobs = [job for job in plan['jobs'] if job['id'] in ('dos-nscf', 'optical-nscf', 'elastic')]
            self.assertEqual(len(jobs), 3)
            for job in jobs:
                text = Path(job['input_file']).read_text()
                shifts = '0 1 0' if job['id'] == 'dos-nscf' else '0 0 0'
                self.assertIn('5 4 1 ' + shifts, text)

    def test_gpaw_optical_preparation_receives_resolved_partial_mesh(self):
        config = DFTConfig(Engine='GPAW', Optical_calc=True, Opt_calc_type='RPA',
            bulk_configuration=self.atoms, Ground_kpts_density=3., Opt_kpts_z=1)
        solver = object.__new__(DFTSolver)
        solver.__dict__.update(vars(config))
        solver.struct = 'sample'
        solver.engine = SimpleNamespace(is_hybrid=Mock(return_value=False),
            prepare_optical_calc=Mock(side_effect=RuntimeError('preparation reached')))
        with self.assertRaisesRegex(RuntimeError, 'preparation reached'):
            solver._opticalcalc_gpaw()
        arguments = solver.engine.prepare_optical_calc.call_args.kwargs
        self.assertIsNone(arguments['kpoint_density'])
        self.assertEqual(arguments['kpoint_size'], (5, 4, 1))
