# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression coverage for spin/U finite-displacement QE phonons."""

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.units import Bohr
from nanoworks.engine import qe
from nanoworks.qe_phonon import (
    prepare_force_plan, run_force_plan, supercell_kpoints, make_phonon, _save_force, postprocess,
)


class TestQEForceParser(unittest.TestCase):
    def test_last_force_block_units_and_fortran_exponents(self):
        text = '''Forces acting on atoms (cartesian axes, Ry/au):
        atom 1 type 1 force = 99 99 99
        Forces acting on atoms (cartesian axes, Ry/au):
        atom 1 type 1 force = 1.0D-3 -.2D-3 0
        atom 2 type 1 force = -1.0D-3 .2D-3 0'''
        actual = qe.parse_pw_forces(text, 2)
        expected = np.array([[.001, -.0002, 0], [-.001, .0002, 0]])
        np.testing.assert_allclose(actual, expected * qe.EV_PER_RYDBERG / Bohr)

    def test_incomplete_last_block_is_not_replaced_with_old_forces(self):
        text = '''Forces acting on atoms (cartesian axes, Ry/au):
        atom 1 type 1 force = 1 0 0
        atom 2 type 1 force = -1 0 0
        Forces acting on atoms (cartesian axes, Ry/au):
        atom 1 type 1 force = 2 0 0'''
        with self.assertRaisesRegex(ValueError, 'incomplete'):
            qe.parse_pw_forces(text, 2)

    def test_unconverged_scf_with_job_done_cannot_supply_forces(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / 'pw.out'
            path.write_text('''Program PWSCF v.7.4.1 starts
            ! total energy = -10 Ry
            convergence NOT achieved
            Forces acting on atoms (cartesian axes, Ry/au):
            atom 1 type 1 force = 0 0 0
            JOB DONE.''', encoding='utf-8')
            with self.assertRaisesRegex(RuntimeError, 'did not converge'):
                qe.read_pw_force_result(path, 1)

    def test_supercell_mesh_retains_ground_resolution(self):
        atoms = Atoms('Ni', cell=[4, 4, 4], pbc=True)
        self.assertEqual(supercell_kpoints(atoms, atoms.repeat((2, 1, 1)), (4, 4, 4)), (2, 4, 4))


class TestQEFiniteDisplacements(unittest.TestCase):
    def setUp(self):
        try:
            import phonopy
        except ModuleNotFoundError:
            self.skipTest('phonopy is optional')
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.pseudo = self.root / 'Ni.upf'
        self.pseudo.write_text('''<UPF version="2.0.1"><PP_HEADER element="Ni" z_valence="10.0"/>
        <PP_PSWFC><PP_CHI.1 label="3D"/></PP_PSWFC></UPF>''', encoding='utf-8')
        self.atoms = Atoms('Ni2', positions=[[0, 0, 0], [2, 2, 2]], cell=[4, 4, 4], pbc=True)
        self.config = SimpleNamespace(Spin_calc=True, Magmom_per_atom=[2, -2], Magmom_single_atom=None,
            Phonon_supercell=np.diag([2, 1, 1]), Phonon_displacement=.01,
            Ground_kpts_density=None, Ground_kpts_x=4, Ground_kpts_y=4, Ground_kpts_z=4,
            Phonon_kpts_x=None, Phonon_kpts_y=None, Phonon_kpts_z=None,
            Occupation_scheme='fermi-dirac', Smearing_width=.05, SCF_accuracy='tight',
            SCF_max_steps=None, SCF_mixing=None, Electronic_solver='default',
            Ground_num_of_bands=20, Wavefunction_cutoff=500, Phonon_PW_cutoff=None,
            Density_cutoff_ratio=4, Gamma=True, Ground_gamma=None, Total_charge=1,
            Hubbard_U={'Ni-3d': 6.0}, XC_calc='PBE', Pseudo_xc='pbe',
            Electrostatic_boundary='periodic', Electrostatic_normal_axis='z',
            Dipole_correction=False, vdW_calc='None', Phonon_path='GX', Phonon_npoints=5,
            Phonon_qpts_x=2, Phonon_qpts_y=2, Phonon_qpts_z=2,
            Phonon_acoustic_sum_rule=True, Phonon_thermal_calc=False,
            Phonon_T_min=0, Phonon_T_max=300, Phonon_T_step=100)

    def plan(self):
        return prepare_force_plan(self.config, self.atoms, self.root / 'Ni',
                                  self.root, {'Ni': 'Ni.upf'})

    def test_force_inputs_preserve_spin_u_and_extensive_settings(self):
        plan = self.plan()
        self.assertGreater(len(plan['jobs']), 1)
        self.assertEqual(plan['electronic_kpoints'], [2, 4, 4])
        for job in plan['jobs']:
            text = job['input_text']
            self.assertIn('tprnfor = .true.', text)
            self.assertIn('nspin = 2', text)
            self.assertIn('HUBBARD (ortho-atomic)', text)
            self.assertIn('U Ni1-3d 6', text)
            self.assertIn('U Ni2-3d 6', text)
            self.assertIn('tot_charge = 2', text)
            self.assertIn('nbnd = 40', text)
            self.assertEqual(job['natoms'], 4)
        self.assertEqual(len({job['state_dir'] for job in plan['jobs']}), len(plan['jobs']))

    def test_resume_invalidates_forces_when_pseudo_content_changes(self):
        plan = self.plan()
        with patch('nanoworks.qe_phonon.qe.run_pw_forces', return_value=np.zeros((4, 3))) as run:
            with patch('nanoworks.qe_phonon.postprocess', return_value={}):
                run_force_plan(plan)
                count = run.call_count
                self.assertEqual(count, len(plan['jobs']))
                run_force_plan(plan)
                self.assertEqual(run.call_count, count)
                self.pseudo.write_text(self.pseudo.read_text() + '\n<!-- revised UPF -->\n')
                updated = self.plan()
                run_force_plan(updated)
                self.assertEqual(run.call_count, count * 2)

    def test_normalized_force_cache_is_readable_json(self):
        plan = self.plan()
        with patch('nanoworks.qe_phonon.qe.run_pw_forces', return_value=np.zeros((4, 3))):
            with patch('nanoworks.qe_phonon.postprocess', return_value={}):
                run_force_plan(plan)
        for job in plan['jobs']:
            record = json.loads(Path(job['cache_file']).read_text())
            self.assertEqual(record['signature'], job['signature'])
            self.assertEqual(len(record['forces_ev_angstrom']), 4)

    def test_harmonic_forces_produce_band_dos_and_force_constants(self):
        plan = self.plan()
        phonon = make_phonon(plan['unitcell'], plan['supercell'], plan['displacement'])
        cells = [phonon.supercell, *phonon.supercells_with_displacements]
        positions = cells[0].positions
        residual = np.array([[.01, 0, 0], [-.01, 0, 0], [.01, 0, 0], [-.01, 0, 0]])
        for job, cell in zip(plan['jobs'], cells):
            _save_force(job, residual - 5.0 * (cell.positions - positions))
        with patch('nanoworks.qe_phonon._plot', return_value=Path('phonon.png')):
            result = postprocess(plan)
        self.assertEqual(result['frequencies']['nmodes'], 6)
        self.assertEqual(result['dos']['natoms'], 2)
        self.assertAlmostEqual(result['residual_force_max_ev_angstrom'], .01)
        self.assertGreater(np.max(result['frequencies']['frequencies_thz']), 0)
        self.assertTrue(Path(result['band_data_file']).is_file())
        self.assertTrue(Path(result['dos_data_file']).is_file())
        self.assertGreater(np.linalg.norm(np.load(result['force_constants_file'])), 0)
