# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression coverage for spin/U finite-displacement QE phonons."""

import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from ase import Atoms
from ase.units import Bohr
from nanoworks.engine import qe
from nanoworks.phonon_cache import force_digest
from nanoworks.phonon_results import qpoint_frequencies
from nanoworks.qe_phonon import (
    prepare_force_plan, run_force_plan, supercell_kpoints, make_phonon, _save_force, postprocess,
    record_force_result, write_mesh_data, has_verified_force,
    pw_executable_identity, begin_force_plan,
)


class TestQEForceParser(unittest.TestCase):
    def test_executable_identity_tracks_content_not_path_or_timestamp(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            first, second = Path(tmpdir) / 'pw.x', Path(tmpdir) / 'pw-copy.x'
            first.write_bytes(b'QE build A')
            second.write_bytes(first.read_bytes())
            stamp = first.stat().st_mtime_ns
            with patch('nanoworks.qe_phonon.shutil.which', return_value=str(first)):
                initial = pw_executable_identity()
                first.write_bytes(b'QE build B')
                os.utime(first, ns=(stamp, stamp))
                changed = pw_executable_identity()
            with patch('nanoworks.qe_phonon.shutil.which', return_value=str(second)):
                copied = pw_executable_identity()
            self.assertEqual(initial['sha256'], copied['sha256'])
            self.assertNotEqual(initial['sha256'], changed['sha256'])

    def test_mesh_report_retains_negative_modes_and_weights(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            frequencies = np.array([[-.05, .3, 1], [-.2, -.11, 2]])
            path, report = write_mesh_data(str(Path(tmpdir) / 'Ni'),
                [[0, 0, 0], [.5, 0, 0]], [1, 3], frequencies)
            self.assertEqual(report['minimum_mesh_qpoint'], [.5, 0, 0])
            self.assertEqual(report['minimum_mesh_mode_index'], 1)
            self.assertEqual(report['negative_mesh_mode_count'], 3)
            self.assertEqual(report['mesh_modes_below_reporting_threshold'], 2)
            self.assertAlmostEqual(report['weighted_mesh_fraction_below_reporting_threshold'], .5)
            np.testing.assert_allclose(np.loadtxt(path)[:, 4:], frequencies)

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
        self.binary = {'path': '/mock/qe/pw.x', 'sha256': 'mock-QE-build-A'}
        identity = patch('nanoworks.qe_phonon.pw_executable_identity', return_value=self.binary)
        identity.start()
        self.addCleanup(identity.stop)
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

    def test_changed_physical_plan_is_rejected_before_force_or_analysis_work(self):
        changes = (
            lambda plan: plan['unitcell']['scaled_positions'][0].__setitem__(0, .1),
            lambda plan: plan['unitcell']['masses'].__setitem__(0, 240),
            lambda plan: plan['unitcell']['magnetic_moments'].__setitem__(0, -2),
            lambda plan: plan.__setitem__('displacement', .02),
            lambda plan: plan['jobs'].reverse(),
            lambda plan: plan['jobs'][0].__setitem__('signature', 'altered'),
        )
        for change in changes:
            with self.subTest(change=change):
                plan = self.plan()
                change(plan)
                with patch('nanoworks.qe_phonon.qe.run_pw_forces') as run:
                    with self.assertRaisesRegex(ValueError, 'geometry or force-job provenance'):
                        run_force_plan(plan)
                run.assert_not_called()
                with patch('nanoworks.qe_phonon.make_phonon') as build:
                    with self.assertRaisesRegex(ValueError, 'geometry or force-job provenance'):
                        postprocess(plan)
                build.assert_not_called()

    def test_old_unsigned_plan_requires_regeneration(self):
        plan = self.plan()
        plan['schema'] = 2
        plan.pop('physical_signature')
        with patch('nanoworks.qe_phonon.qe.run_pw_forces') as run:
            with self.assertRaisesRegex(ValueError, 'lacks physical provenance'):
                run_force_plan(plan)
        run.assert_not_called()

    def test_analysis_choices_do_not_change_physical_plan_signature(self):
        from nanoworks.qe_phonon import _physical_signature, _validate_plan_resources
        plan = self.plan()
        original = plan['physical_signature']
        plan['dos_mesh'] = [4, 4, 4]
        plan['acoustic_sum_rule'] = False
        plan['thermal'] = True
        plan['temperature'] = [0, 200, 50]
        plan['band_path']['kpoints'][0] = [.1, 0, 0]
        self.assertEqual(_physical_signature(plan), original)
        _validate_plan_resources(plan)

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

    def test_force_plan_preserves_custom_masses_without_repeating_force_scfs(self):
        self.atoms.set_masses([60, 60])
        first = self.plan()
        self.assertEqual(first['unitcell']['masses'], [60, 60])
        phonon = make_phonon(first['unitcell'], first['supercell'], first['displacement'])
        np.testing.assert_array_equal(phonon.unitcell.masses, [60, 60])
        np.testing.assert_array_equal(phonon.supercell.masses, [60] * 4)
        self.atoms.set_masses([62, 62])
        second = self.plan()
        self.assertEqual(second['unitcell']['masses'], [62, 62])
        # Isotope mass changes the dynamical analysis, not Born-Oppenheimer forces.
        self.assertEqual([job['signature'] for job in first['jobs']],
                         [job['signature'] for job in second['jobs']])

    def test_mass_scaled_frequencies_retain_inverse_square_root_relation(self):
        self.atoms.set_masses([60, 60])
        plan = self.plan()
        first = make_phonon(plan['unitcell'], plan['supercell'], plan['displacement'])
        count = len(first.supercell)
        constants = np.zeros((count, count, 3, 3))
        for index in range(count):
            constants[index, index] = np.eye(3)
        first.force_constants = constants
        heavier = dict(plan['unitcell'], masses=[240, 240])
        second = make_phonon(heavier, plan['supercell'], plan['displacement'])
        second.force_constants = constants
        np.testing.assert_allclose(qpoint_frequencies(second, [[.2, 0, 0]]),
                                   qpoint_frequencies(first, [[.2, 0, 0]]) / 2)

    def test_legacy_force_plan_without_masses_uses_phonopy_defaults(self):
        plan = self.plan()
        legacy = dict(plan['unitcell'])
        del legacy['masses']
        phonon = make_phonon(legacy, plan['supercell'], plan['displacement'])
        self.assertTrue(np.isfinite(phonon.unitcell.masses).all())
        self.assertTrue(np.all(phonon.unitcell.masses > 0))

    def test_force_plan_does_not_truncate_fractional_electronic_mesh(self):
        self.config.Phonon_kpts_x = 2.5
        with patch('nanoworks.qe_phonon.qe.render_pw_input') as render:
            with self.assertRaisesRegex(ValueError, 'Phonon_kpts_x'):
                self.plan()
        render.assert_not_called()

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
            self.assertEqual(record['schema'], 1)
            self.assertEqual(record['units'], 'eV/Angstrom')
            self.assertEqual(record['force_sha256'], force_digest(record['forces_ev_angstrom']))

    def test_postprocess_rejects_unverified_output_without_running_qe(self):
        plan = self.plan()
        Path(plan['jobs'][0]['output_file']).write_text('old SCF output')
        with patch('nanoworks.qe_phonon.qe.read_pw_force_result') as read:
            with patch('nanoworks.qe_phonon.qe.run_pw_forces') as run:
                with self.assertRaisesRegex(ValueError, 'No matching verified forces'):
                    postprocess(plan)
        read.assert_not_called()
        run.assert_not_called()
        report = json.loads(Path(plan['struct'] + '-PHONON-QE-Result-Summary.json').read_text())
        self.assertEqual(report['status'], 'failed')
        self.assertEqual(report['completed_force_jobs'], 0)

    def test_shell_resume_checks_signature_shape_and_current_input(self):
        plan = self.plan()
        job = plan['jobs'][0]
        Path(job['input_file']).write_text(job['input_text'])
        self.assertFalse(has_verified_force(plan, job['id']))
        _save_force(job, np.zeros((4, 3)), binary_identity=self.binary)
        self.assertTrue(has_verified_force(plan, job['id']))
        cache = Path(job['cache_file'])
        record = json.loads(cache.read_text())
        record['signature'] = 'stale'
        cache.write_text(json.dumps(record))
        self.assertFalse(has_verified_force(plan, job['id']))
        record['signature'] = job['signature']
        record['forces_ev_angstrom'] = [[0, 0, 0]]
        cache.write_text(json.dumps(record))
        self.assertFalse(has_verified_force(plan, job['id']))
        cache.write_text('{ interrupted write')
        self.assertFalse(has_verified_force(plan, job['id']))
        _save_force(job, np.zeros((4, 3)), binary_identity=self.binary)
        Path(job['input_file']).write_text(job['input_text'] + '\n! altered input\n')
        with self.assertRaisesRegex(ValueError, 'input changed'):
            has_verified_force(plan, job['id'])

    def test_changed_finite_force_is_recomputed_individually(self):
        plan = self.plan()
        expected = np.arange(12, dtype=float).reshape(4, 3) / 10
        with patch('nanoworks.qe_phonon.qe.run_pw_forces', return_value=expected) as run:
            with patch('nanoworks.qe_phonon.postprocess', return_value={}):
                run_force_plan(plan)
                count = run.call_count
                job = plan['jobs'][0]
                Path(job['input_file']).write_text(job['input_text'])
                cache = Path(job['cache_file'])
                record = json.loads(cache.read_text())
                # JSON formatting has no effect on the numerical content hash.
                cache.write_text(json.dumps(record, sort_keys=True, separators=(',', ':')))
                self.assertTrue(has_verified_force(plan, job['id']))
                record['forces_ev_angstrom'][0][0] += .1
                cache.write_text(json.dumps(record))
                self.assertFalse(has_verified_force(plan, job['id']))
                run_force_plan(plan)
                self.assertEqual(run.call_count, count + 1)
                self.assertEqual(run.call_args.args[0], job['input_file'])
                self.assertTrue(has_verified_force(plan, job['id']))
                np.testing.assert_array_equal(
                    json.loads(cache.read_text())['forces_ev_angstrom'], expected)

    def test_postprocess_rejects_altered_force_content_without_running_qe(self):
        plan = self.plan()
        for job in plan['jobs']:
            _save_force(job, np.zeros((4, 3)), binary_identity=self.binary)
        cache = Path(plan['jobs'][0]['cache_file'])
        record = json.loads(cache.read_text())
        record['forces_ev_angstrom'][0][0] = 1
        cache.write_text(json.dumps(record))
        with patch('nanoworks.qe_phonon.qe.run_pw_forces') as run:
            with self.assertRaisesRegex(ValueError, 'No matching verified forces'):
                postprocess(plan)
        run.assert_not_called()
        report = json.loads(Path(plan['struct'] + '-PHONON-QE-Result-Summary.json').read_text())
        self.assertEqual(report['status'], 'failed')
        self.assertEqual(report['completed_force_jobs'], len(plan['jobs']) - 1)

    def test_shell_resume_rejects_old_or_invalid_integrity_metadata(self):
        plan = self.plan()
        job = plan['jobs'][0]
        Path(job['input_file']).write_text(job['input_text'])
        cache = Path(job['cache_file'])
        for field, value in (('force_sha256', None), ('schema', 99), ('units', 'Ry/Bohr')):
            with self.subTest(field=field):
                _save_force(job, np.zeros((4, 3)), binary_identity=self.binary)
                record = json.loads(cache.read_text())
                if value is None:
                    del record[field]
                else:
                    record[field] = value
                cache.write_text(json.dumps(record))
                self.assertFalse(has_verified_force(plan, job['id']))

    def test_shell_resume_refuses_changed_pseudopotential(self):
        plan = self.plan()
        job = plan['jobs'][0]
        Path(job['input_file']).write_text(job['input_text'])
        _save_force(job, np.zeros((4, 3)), binary_identity=self.binary)
        self.pseudo.write_text(self.pseudo.read_text() + '\n<!-- changed -->\n')
        with self.assertRaisesRegex(ValueError, 'Pseudopotential content changed'):
            has_verified_force(plan, job['id'])

    def test_changed_executable_stops_before_reusing_or_launching_forces(self):
        plan = self.plan()
        other = {'path': '/mock/qe/pw.x', 'sha256': 'mock-QE-build-B'}
        with patch('nanoworks.qe_phonon.pw_executable_identity', return_value=other):
            with patch('nanoworks.qe_phonon.qe.run_pw_forces') as run:
                with self.assertRaisesRegex(ValueError, 'pw.x content changed'):
                    run_force_plan(plan)
        run.assert_not_called()

    def test_deck_prepared_without_qe_binds_on_execution_host(self):
        with patch('nanoworks.qe_phonon.pw_executable_identity', return_value=None):
            plan = self.plan()
        self.assertIsNone(plan['pw_executable'])
        begin_force_plan(plan)
        stored = json.loads(Path(plan['manifest_file']).read_text())
        self.assertEqual(stored['pw_executable'], self.binary)
        with patch('nanoworks.qe_phonon.pw_executable_identity',
                   return_value={'path': '/other/pw.x', 'sha256': 'other-build'}):
            with self.assertRaisesRegex(ValueError, 'pw.x content changed'):
                begin_force_plan(stored)

    def test_executable_change_during_scf_never_records_result(self):
        plan = self.plan()
        other = {'path': '/mock/qe/pw.x', 'sha256': 'mock-QE-build-B'}
        with patch('nanoworks.qe_phonon.pw_executable_identity',
                   side_effect=[self.binary, self.binary, other]):
            with patch('nanoworks.qe_phonon.qe.run_pw_forces', return_value=np.zeros((4, 3))):
                with self.assertRaisesRegex(ValueError, 'changed during'):
                    run_force_plan(plan)
        self.assertFalse(Path(plan['jobs'][0]['cache_file']).exists())

    def test_mixed_executable_forces_are_rejected_without_requiring_qe(self):
        plan = self.plan()
        for index, job in enumerate(plan['jobs']):
            binary = self.binary if index == 0 else {'sha256': 'other-QE-build'}
            _save_force(job, np.zeros((4, 3)), binary_identity=binary)
        with patch('nanoworks.qe_phonon.pw_executable_identity', return_value=None) as identity:
            with self.assertRaisesRegex(ValueError, 'inconsistent pw.x'):
                postprocess(plan)
        identity.assert_not_called()

    def test_record_requires_current_input_and_fresh_converged_output(self):
        plan = self.plan()
        job = plan['jobs'][0]
        Path(job['input_file']).write_text(job['input_text'])
        output = Path(job['output_file'])
        output.write_text('Program PWSCF v.7.4.1 starts\n'
            '! total energy = -10 Ry\nconvergence has been achieved\n'
            'Forces acting on atoms (cartesian axes, Ry/au):\n' +
            '\n'.join(f'atom {index} type 1 force = 0 0 0' for index in range(1, 5)) +
            '\nJOB DONE.\n')
        manifest_time = Path(plan['manifest_file']).stat().st_mtime_ns
        os.utime(output, ns=(manifest_time - 1_000_000_000, manifest_time - 1_000_000_000))
        with self.assertRaisesRegex(ValueError, 'predates'):
            record_force_result(plan, job['id'])
        os.utime(output, ns=(manifest_time + 1_000_000_000, manifest_time + 1_000_000_000))
        record_force_result(plan, job['id'])
        self.assertTrue(Path(job['cache_file']).is_file())
        Path(job['input_file']).write_text(job['input_text'] + '\n! changed\n')
        with self.assertRaisesRegex(ValueError, 'input changed'):
            record_force_result(plan, job['id'])

    def test_scf_failure_preserves_completed_force_for_resume(self):
        plan = self.plan()
        with patch('nanoworks.qe_phonon.qe.run_pw_forces',
                   side_effect=[np.zeros((4, 3)), RuntimeError('SCF failed')]):
            with self.assertRaisesRegex(RuntimeError, 'SCF failed'):
                run_force_plan(plan)
        report = json.loads(Path(plan['struct'] + '-PHONON-QE-Result-Summary.json').read_text())
        self.assertEqual(report['status'], 'failed')
        self.assertEqual(report['completed_force_jobs'], 1)
        with patch('nanoworks.qe_phonon.qe.run_pw_forces', return_value=np.zeros((4, 3))) as run:
            with patch('nanoworks.qe_phonon.postprocess', return_value={}):
                run_force_plan(plan)
        self.assertEqual(run.call_count, len(plan['jobs']) - 1)

    def test_harmonic_forces_produce_band_dos_and_force_constants(self):
        plan = self.plan()
        phonon = make_phonon(plan['unitcell'], plan['supercell'], plan['displacement'])
        cells = [phonon.supercell, *phonon.supercells_with_displacements]
        positions = cells[0].positions
        residual = np.array([[.01, 0, 0], [-.01, 0, 0], [.01, 0, 0], [-.01, 0, 0]])
        for job, cell in zip(plan['jobs'], cells):
            _save_force(job, residual - 5.0 * (cell.positions - positions), binary_identity=self.binary)
        with patch('nanoworks.qe_phonon._plot', return_value=Path('phonon.png')):
            with patch('nanoworks.qe_phonon.pw_executable_identity', return_value=None) as identity:
                with patch('nanoworks.qe_phonon.qpoint_frequencies', wraps=qpoint_frequencies) as batch:
                    result = postprocess(plan)
        identity.assert_not_called()
        batch.assert_called_once()
        np.testing.assert_array_equal(batch.call_args.args[1], plan['band_path']['kpoints'])
        self.assertEqual(np.shape(result['frequencies']['frequencies_thz']),
                         (len(plan['band_path']['kpoints']), 6))
        self.assertEqual(result['frequencies']['nmodes'], 6)
        self.assertEqual(result['dos']['natoms'], 2)
        self.assertAlmostEqual(result['residual_force_max_ev_angstrom'], .01)
        self.assertGreater(np.max(result['frequencies']['frequencies_thz']), 0)
        self.assertTrue(Path(result['band_data_file']).is_file())
        self.assertTrue(Path(result['dos_data_file']).is_file())
        self.assertTrue(Path(result['mesh_data_file']).is_file())
        self.assertEqual(result['status'], 'complete')
        self.assertEqual(result['pw_executable_sha256'], self.binary['sha256'])
        self.assertGreater(np.linalg.norm(np.load(result['force_constants_file'])), 0)
