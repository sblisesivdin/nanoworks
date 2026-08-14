import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

from nanoworks.engine import load_engine_module
from ase import Atoms
from nanoworks.engine.qe import (
    QE_REFERENCE_VERSION,
    ev_to_rydberg,
    build_control_settings,
    build_system_settings,
    build_cell_parameters,
    build_atomic_positions,
    build_atomic_species,
    build_kpoint_settings,
    build_occupation_settings,
    build_electrons_settings,
    format_qe_value,
    render_namelist,
    render_scf_input,
    rydberg_to_ev,
    resolve_qe_executable,
    build_qe_command,
    run_qe_program,
    parse_pw_output,
)


class TestQEEngine(unittest.TestCase):

    def test_reference_version_is_qe_72(self):
        self.assertEqual(QE_REFERENCE_VERSION, (7, 2))

    def test_ev_to_rydberg(self):
        self.assertAlmostEqual(
            ev_to_rydberg(13.605693122994),
            1.0,
        )

    def test_control_settings_defaults(self):
        settings = build_control_settings()

        self.assertEqual(settings['calculation'], 'scf')
        self.assertEqual(settings['prefix'], 'nanoworks')
        self.assertNotIn('pseudo_dir', settings)
        self.assertNotIn('outdir', settings)

    def test_control_settings_optional_paths(self):
        settings = build_control_settings(
            calculation='SCF',
            prefix='GaAs',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe',
        )

        self.assertEqual(settings['calculation'], 'scf')
        self.assertEqual(settings['prefix'], 'GaAs')
        self.assertEqual(settings['pseudo_dir'], '/tmp/pseudos')
        self.assertEqual(settings['outdir'], '/tmp/qe')

    def test_system_settings_basic(self):
        settings = build_system_settings(
            cutoff_ev=340,
            nat=2,
            ntyp=2,
        )

        self.assertEqual(settings['ibrav'], 0)
        self.assertEqual(settings['nat'], 2)
        self.assertEqual(settings['ntyp'], 2)
        self.assertAlmostEqual(
            settings['ecutwfc'],
            340 / 13.605693122994,
        )

        self.assertNotIn('tot_charge', settings)
        self.assertNotIn('nbnd', settings)
        self.assertNotIn('nspin', settings)
        self.assertNotIn('ecutrho', settings)

    def test_system_settings_optional_values(self):
        settings = build_system_settings(
            cutoff_ev=400,
            nat=4,
            ntyp=2,
            total_charge=1.0,
            nbands=20,
            spinpol=True,
        )

        self.assertEqual(settings['tot_charge'], 1.0)
        self.assertEqual(settings['nbnd'], 20)
        self.assertEqual(settings['nspin'], 2)

    def test_qe_engine_can_be_loaded(self):
        engine = load_engine_module('qe')

        self.assertEqual(
            engine.QE_REFERENCE_VERSION,
            (7, 2),
        )
    
    def test_cell_parameters_are_built_in_angstrom(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
            cell=[
                [5.65, 0.0, 0.0],
                [0.0, 5.65, 0.0],
                [0.0, 0.0, 5.65],
            ],
            pbc=True,
        )

        card = build_cell_parameters(atoms)

        self.assertEqual(card['option'], 'angstrom')
        self.assertEqual(
            card['vectors'],
            [
                (5.65, 0.0, 0.0),
                (0.0, 5.65, 0.0),
                (0.0, 0.0, 5.65),
            ],
        )
    
    def test_atomic_positions_are_built_in_angstrom(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.4125, 1.4125, 1.4125],
            ],
            cell=[5.65, 5.65, 5.65],
            pbc=True,
        )

        card = build_atomic_positions(atoms)

        self.assertEqual(card['option'], 'angstrom')
        self.assertEqual(
            card['positions'],
            [
                ('Ga', 0.0, 0.0, 0.0),
                ('As', 1.4125, 1.4125, 1.4125),
            ],
        )
    
    def test_atomic_species_use_pseudopotential_mapping(self):
        atoms = Atoms(
            'GaAsGa',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
                [2.0, 2.0, 2.0],
            ],
        )

        species = build_atomic_species(
            atoms,
            {
                'Ga': 'Ga.upf',
                'As': 'As.upf',
            },
        )

        self.assertEqual(len(species), 2)

        self.assertEqual(species[0][0], 'Ga')
        self.assertEqual(species[0][2], 'Ga.upf')

        self.assertEqual(species[1][0], 'As')
        self.assertEqual(species[1][2], 'As.upf')

        self.assertGreater(species[0][1], 0.0)
        self.assertGreater(species[1][1], 0.0)
    
    def test_atomic_species_reject_missing_pseudopotential(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.0, 1.0, 1.0],
            ],
        )

        with self.assertRaisesRegex(
            ValueError,
            'Missing pseudopotential mapping for: As',
        ):
            build_atomic_species(
                atoms,
                {
                    'Ga': 'Ga.upf',
                },
            )
    
    def test_gamma_centered_kpoint_mesh(self):
        settings = build_kpoint_settings(
            (4, 4, 4),
            gamma=True,
        )

        self.assertEqual(settings['option'], 'automatic')
        self.assertEqual(settings['size'], (4, 4, 4))
        self.assertEqual(settings['shift'], (0, 0, 0))
    
    def test_monkhorst_pack_even_mesh_is_shifted(self):
        settings = build_kpoint_settings(
            (4, 6, 8),
            gamma=False,
        )

        self.assertEqual(settings['size'], (4, 6, 8))
        self.assertEqual(settings['shift'], (1, 1, 1))
    
    def test_monkhorst_pack_shift_depends_on_mesh_parity(self):
        settings = build_kpoint_settings(
            (4, 5, 6),
            gamma=False,
        )

        self.assertEqual(settings['shift'], (1, 0, 1))
    
    def test_kpoint_mesh_rejects_nonpositive_values(self):
        with self.assertRaises(ValueError):
            build_kpoint_settings((4, 0, 4))
    
    def test_fixed_occupation_settings(self):
        settings = build_occupation_settings('fixed')

        self.assertEqual(
            settings,
            {'occupations': 'fixed'},
        )
    
    def test_smearing_settings_convert_width_to_rydberg(self):
        settings = build_occupation_settings(
            occupations='smearing',
            smearing='fermi-dirac',
            width_ev=0.1,
        )

        self.assertEqual(
            settings['occupations'],
            'smearing',
        )
        self.assertEqual(
            settings['smearing'],
            'fermi-dirac',
        )
        self.assertAlmostEqual(
            settings['degauss'],
            0.1 / 13.605693122994,
        )
    
    def test_cold_smearing_alias(self):
        settings = build_occupation_settings(
            occupations='smearing',
            smearing='cold',
            width_ev=0.05,
        )

        self.assertEqual(
            settings['smearing'],
            'marzari-vanderbilt',
        )
    
    def test_tetrahedra_has_no_smearing_parameters(self):
        settings = build_occupation_settings(
            'tetrahedra'
        )

        self.assertEqual(
            settings,
            {'occupations': 'tetrahedra'},
        )
    
    def test_smearing_requires_width(self):
        with self.assertRaises(ValueError):
            build_occupation_settings(
                occupations='smearing',
                smearing='gaussian',
            )
    
    def test_electrons_settings_can_use_qe_defaults(self):
        settings = build_electrons_settings()

        self.assertEqual(settings, {})

    def test_electrons_settings_explicit_values(self):
        settings = build_electrons_settings(
            conv_thr=1.0e-8,
            mixing_beta=0.3,
            electron_maxstep=200,
            diagonalization='david',
        )

        self.assertEqual(settings['conv_thr'], 1.0e-8)
        self.assertEqual(settings['mixing_beta'], 0.3)
        self.assertEqual(settings['electron_maxstep'], 200)
        self.assertEqual(settings['diagonalization'], 'david')

    def test_qe_value_formatter(self):
        self.assertEqual(
            format_qe_value('scf'),
            "'scf'",
        )

        self.assertEqual(
            format_qe_value(True),
            '.true.',
        )

        self.assertEqual(
            format_qe_value(False),
            '.false.',
        )

        self.assertEqual(
            format_qe_value(4),
            '4',
        )

    def test_render_namelist(self):
        text = render_namelist(
            'control',
            {
                'calculation': 'scf',
                'prefix': 'GaAs',
            },
        )

        self.assertEqual(
            text,
            "\n".join([
                "&CONTROL",
                "  calculation = 'scf',",
                "  prefix = 'GaAs',",
                "/",
            ]),
        )

    def test_render_complete_scf_input(self):
        atoms = Atoms(
            'GaAs',
            positions=[
                [0.0, 0.0, 0.0],
                [1.4125, 1.4125, 1.4125],
            ],
            cell=[
                [2.825, 2.825, 0.0],
                [2.825, 0.0, 2.825],
                [0.0, 2.825, 2.825],
            ],
            pbc=True,
        )

        text = render_scf_input(
            atoms=atoms,
            pseudopotentials={
                'Ga': 'Ga.upf',
                'As': 'As.upf',
            },
            cutoff_ev=340,
            kpoint_size=(4, 4, 4),
            gamma=False,
            occupations='fixed',
            prefix='GaAs',
        )

        self.assertIn("&CONTROL", text)
        self.assertIn("&SYSTEM", text)
        self.assertIn("&ELECTRONS", text)

        self.assertIn(
            "calculation = 'scf'",
            text,
        )

        self.assertIn(
            "prefix = 'GaAs'",
            text,
        )

        self.assertIn(
            "ibrav = 0",
            text,
        )

        self.assertIn(
            "nat = 2",
            text,
        )

        self.assertIn(
            "ntyp = 2",
            text,
        )

        self.assertIn(
            "occupations = 'fixed'",
            text,
        )

        self.assertIn(
            "ATOMIC_SPECIES",
            text,
        )

        self.assertIn(
            "Ga.upf",
            text,
        )

        self.assertIn(
            "As.upf",
            text,
        )

        self.assertIn(
            "ATOMIC_POSITIONS angstrom",
            text,
        )

        self.assertIn(
            "K_POINTS automatic",
            text,
        )

        self.assertIn(
            "4 4 4 1 1 1",
            text,
        )

        self.assertIn(
            "CELL_PARAMETERS angstrom",
            text,
        )

    def test_render_scf_input_with_smearing(self):
        atoms = Atoms(
            'Al',
            positions=[[0.0, 0.0, 0.0]],
            cell=[4.05, 4.05, 4.05],
            pbc=True,
        )

        text = render_scf_input(
            atoms=atoms,
            pseudopotentials={
                'Al': 'Al.upf',
            },
            cutoff_ev=400,
            kpoint_size=(6, 6, 6),
            occupations='smearing',
            smearing='cold',
            width_ev=0.1,
        )

        self.assertIn(
            "occupations = 'smearing'",
            text,
        )

        self.assertIn(
            "smearing = 'marzari-vanderbilt'",
            text,
        )

        self.assertIn(
            "degauss = ",
            text,
        )
    
    def test_rydberg_to_ev(self):
        self.assertAlmostEqual(
            rydberg_to_ev(1.0),
            13.605693122994,
        )

    @patch('nanoworks.engine.qe.shutil.which')
    def test_resolve_qe_executable_from_path(
        self,
        which,
    ):
        which.return_value = '/opt/qe/bin/pw.x'

        resolved = resolve_qe_executable(
            'pw.x'
        )

        self.assertEqual(
            resolved,
            '/opt/qe/bin/pw.x',
        )

        which.assert_called_once_with(
            'pw.x'
        )

    @patch('nanoworks.engine.qe.resolve_qe_executable')
    def test_build_qe_command_with_mpi(
        self,
        resolve,
    ):
        resolve.return_value = '/opt/qe/bin/pw.x'

        command = build_qe_command(
            input_file='si.in',
            executable='pw.x',
            launcher=[
                'mpiexec',
                '-np',
                '8',
            ],
        )

        self.assertEqual(
            command,
            [
                'mpiexec',
                '-np',
                '8',
                '/opt/qe/bin/pw.x',
                '-i',
                'si.in',
            ],
        )

    @patch('nanoworks.engine.qe.resolve_qe_executable')
    def test_qe_launcher_rejects_shell_string(
        self,
        resolve,
    ):
        resolve.return_value = '/opt/qe/bin/pw.x'

        with self.assertRaises(TypeError):
            build_qe_command(
                input_file='si.in',
                launcher='mpiexec -np 8',
            )

    def test_parse_pw_output(self):
        output_text = """
         Program PWSCF v.7.2 starts

    !    total energy              =     -15.12345678 Ry

         JOB DONE.
    """

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'si.out'

            path.write_text(
                output_text
            )

            result = parse_pw_output(
                path
            )

            self.assertTrue(
                result['job_done']
            )

            self.assertAlmostEqual(
                result['total_energy_ry'],
                -15.12345678,
            )

            self.assertAlmostEqual(
                result['total_energy_ev'],
                -15.12345678
                * 13.605693122994,
            )

    @patch('nanoworks.engine.qe.subprocess.run')
    @patch('nanoworks.engine.qe.resolve_qe_executable')
    def test_run_qe_program(
        self,
        resolve,
        run,
    ):
        resolve.return_value = '/opt/qe/bin/pw.x'

        run.return_value.returncode = 0

        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)

            input_file = (
                tmp / 'si.in'
            )

            output_file = (
                tmp / 'si.out'
            )

            input_file.write_text(
                '&CONTROL\n/\n'
            )

            result = run_qe_program(
                input_file=input_file,
                output_file=output_file,
            )

            self.assertEqual(
                result['returncode'],
                0,
            )

            self.assertTrue(
                output_file.exists()
            )

            command = (
                run.call_args.args[0]
            )

            self.assertEqual(
                command[-2:],
                [
                    '-i',
                    str(input_file),
                ],
            )



if __name__ == '__main__':
    unittest.main()
