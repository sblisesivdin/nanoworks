import unittest
import tempfile
from pathlib import Path
from unittest.mock import patch

from nanoworks.engine import load_engine_module
from ase import Atoms
from ase.build import bulk
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
    render_pw_input,
    render_scf_input,
    render_nscf_input,
    rydberg_to_ev,
    build_qe_launcher,
    resolve_qe_executable,
    build_qe_command,
    run_qe_program,
    parse_pw_output,
    resolve_qe_kpoint_size,
    resolve_qe_occupation,
    validate_qe_version,
    validate_qe_xc,
    run_scf,
    run_nscf,
    has_qe_state,
    render_dos_input,
    run_dos,
    parse_dos_output,
    render_projwfc_input,
    run_projwfc,
    parse_projwfc_pdos_file,
    aggregate_projwfc_pdos,
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

    def test_render_dos_input_for_tetrahedra(self):
        text = render_dos_input(
            prefix='nanoworks',
            outdir='/tmp/qe-state',
            fildos='/tmp/gaas.dos',
            bz_sum='tetrahedra',
            emin=-10.0,
            emax=10.0,
            delta_e=0.02,
        )

        self.assertIn(
            '&DOS',
            text,
        )

        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )

        self.assertIn(
            "outdir = '/tmp/qe-state'",
            text,
        )

        self.assertIn(
            "bz_sum = 'tetrahedra'",
            text,
        )

        self.assertIn(
            'Emin = -10',
            text,
        )

        self.assertIn(
            'Emax = 10',
            text,
        )

        self.assertIn(
            'DeltaE = 0.02',
            text,
        )

        self.assertIn(
            "fildos = '/tmp/gaas.dos'",
            text,
        )

        self.assertNotIn(
            'degauss',
            text,
        )

    def test_render_dos_input_rejects_invalid_energy_range(self):
        with self.assertRaisesRegex(
            ValueError,
            'Emax must be greater than Emin',
        ):
            render_dos_input(
                emin=5.0,
                emax=-5.0,
            )

    def test_render_dos_input_rejects_invalid_bz_sum(self):
        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE DOS BZ summation method',
        ):
            render_dos_input(
                bz_sum='invalid',
            )

    def test_run_dos_requires_qe_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(
                tmpdir
            )

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE electronic state',
            ):
                run_dos(
                    input_file=tmpdir / 'dos.in',
                    output_file=tmpdir / 'dos.out',
                    state_dir=tmpdir / 'state',
                    dos_file=tmpdir / 'result.dos',
                )

    def test_render_nscf_input(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        text = render_nscf_input(
            atoms=atoms,
            pseudopotentials={
                'Si': 'Si.upf',
            },
            cutoff_ev=400.0,
            kpoint_size=(8, 8, 8),
            occupations='tetrahedra',
            prefix='nanoworks',
            pseudo_dir='/tmp/pseudos',
            outdir='/tmp/qe-state',
        )

        self.assertIn(
            "calculation = 'nscf'",
            text,
        )

        self.assertIn(
            "occupations = 'tetrahedra'",
            text,
        )

        self.assertIn(
            "8 8 8 1 1 1",
            text,
        )

        self.assertIn(
            "prefix = 'nanoworks'",
            text,
        )

        self.assertIn(
            "outdir = '/tmp/qe-state'",
            text,
        )

    def test_render_pw_input_rejects_unsupported_calculation(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with self.assertRaisesRegex(
            ValueError,
            'Unsupported QE pw.x calculation type',
        ):
            render_pw_input(
                calculation='relax',
                atoms=atoms,
                pseudopotentials={
                    'Si': 'Si.upf',
                },
                cutoff_ev=400.0,
                kpoint_size=(4, 4, 4),
            )
            
    def test_rydberg_to_ev(self):
        self.assertAlmostEqual(
            rydberg_to_ev(1.0),
            13.605693122994,
        )

    def test_has_qe_state_accepts_complete_saved_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_dir = Path(tmpdir)

            save_dir = (
                state_dir
                / 'nanoworks.save'
            )

            save_dir.mkdir()

            (
                save_dir
                / 'data-file-schema.xml'
            ).write_text(
                '<espresso/>',
                encoding='utf-8',
            )

            self.assertTrue(
                has_qe_state(
                    state_dir
                )
            )

    def test_has_qe_state_rejects_incomplete_saved_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_dir = Path(tmpdir)

            (
                state_dir
                / 'nanoworks.save'
            ).mkdir()

            self.assertFalse(
                has_qe_state(
                    state_dir
                )
            )
    
    def test_has_qe_state_rejects_missing_state_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            state_dir = (
                Path(tmpdir)
                / 'missing-state'
            )

            self.assertFalse(
                has_qe_state(
                    state_dir
                )
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
    
         the Fermi energy is     5.4321 ev

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
            
            self.assertEqual(
                result['qe_version'],
                (7, 2),
            )

            self.assertAlmostEqual(
                result['fermi_energy_ev'],
                5.4321,
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
    
    def test_parse_pw_output_with_patch_version(self):
        output_text = """
         Program PWSCF v.7.2.1 starts

         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmp:
            path = (
                Path(tmp)
                / 'pw.out'
            )

            path.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_output(
                path
            )

            self.assertEqual(
                result['qe_version'],
                (7, 2, 1),
            )
        
    def test_parse_pw_output_without_version(self):
        output_text = """
         JOB DONE.
        """

        with tempfile.TemporaryDirectory() as tmp:
            path = (
                Path(tmp)
                / 'pw.out'
            )

            path.write_text(
                output_text,
                encoding='utf-8',
            )

            result = parse_pw_output(
                path
            )

            self.assertIsNone(
                result['qe_version']
            )
            
    def test_parse_pw_output_without_fermi_energy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_file = (
                Path(tmpdir)
                / 'pw.out'
            )

            output_file.write_text(
                "\n".join([
                    "!    total energy = -10.00000000 Ry",
                    "JOB DONE.",
                ]),
                encoding='utf-8',
            )

            result = parse_pw_output(
                output_file
            )

            self.assertIsNone(
                result['fermi_energy_ev']
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

            environment = (
                run.call_args.kwargs['env']
            )

            self.assertEqual(
                environment['OMP_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['OPENBLAS_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['MKL_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['VECLIB_MAXIMUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['NUMEXPR_NUM_THREADS'],
                '1',
            )

            self.assertEqual(
                environment['OMP_DYNAMIC'],
                'FALSE',
            )
    
    def test_qe_launcher_is_none_for_single_core(self):
        launcher = build_qe_launcher(
            parallel_cores=1
        )

        self.assertIsNone(
            launcher
        )
    
    @patch('nanoworks.engine.qe.shutil.which')
    def test_qe_launcher_uses_mpiexec(
        self,
        which,
    ):
        which.side_effect = (
            lambda name: (
                '/usr/bin/mpiexec'
                if name == 'mpiexec'
                else None
            )
        )

        launcher = build_qe_launcher(
            parallel_cores=8
        )

        self.assertEqual(
            launcher,
            [
                '/usr/bin/mpiexec',
                '-np',
                '8',
            ],
        )

    @patch('nanoworks.engine.qe.shutil.which')
    def test_qe_launcher_uses_srun(
        self,
        which,
    ):
        available = {
            'mpiexec': None,
            'mpirun': None,
            'srun': '/usr/bin/srun',
        }

        which.side_effect = (
            lambda name: available[name]
        )

        launcher = build_qe_launcher(
            parallel_cores=16
        )

        self.assertEqual(
            launcher,
            [
                '/usr/bin/srun',
                '-n',
                '16',
            ],
        )

    def test_qe_launcher_rejects_nonpositive_core_count(self):
        with self.assertRaises(ValueError):
            build_qe_launcher(
                parallel_cores=0
            )

    @patch('nanoworks.engine.qe.shutil.which')
    def test_qe_launcher_requires_mpi_command(
        self,
        which,
    ):
        which.return_value = None

        with self.assertRaises(
            FileNotFoundError
        ):
            build_qe_launcher(
                parallel_cores=8
            )

    def test_qe_kpoint_density_resolves_mesh(self):
        atoms = Atoms(
            'Si2',
            positions=[
                [0.0, 0.0, 0.0],
                [1.35, 1.35, 1.35],
            ],
            cell=[
                [5.4, 0.0, 0.0],
                [0.0, 5.4, 0.0],
                [0.0, 0.0, 5.4],
            ],
            pbc=True,
        )

        mesh = resolve_qe_kpoint_size(
            atoms,
            density=2.5,
            size=(1, 1, 1),
        )

        self.assertEqual(
            len(mesh),
            3,
        )

        self.assertTrue(
            all(value > 0 for value in mesh)
        )

    def test_resolve_qe_occupation_supports_tetrahedra(self):
        result = resolve_qe_occupation(
            'tetrahedra'
        )

        self.assertEqual(
            result,
            {
                'occupations': 'tetrahedra',
                'smearing': None,
                'width_ev': None,
            },
        )

    def test_resolve_qe_occupation_supports_optimized_tetrahedra(self):
        result = resolve_qe_occupation(
            'tetrahedra-opt'
        )

        self.assertEqual(
            result['occupations'],
            'tetrahedra_opt',
        )

    def test_qe_occupation_translates_fermi_dirac(self):
        settings = resolve_qe_occupation(
            {
                'name': 'fermi-dirac',
                'width': 0.05,
            }
        )

        self.assertEqual(
            settings,
            {
                'occupations': 'smearing',
                'smearing': 'fermi-dirac',
                'width_ev': 0.05,
            },
        )

    def test_qe_xc_accepts_pbe(self):
        self.assertEqual(
            validate_qe_xc('PBE'),
            'pbe',
        )

    def test_qe_xc_accepts_pbe(self):
        self.assertEqual(
            validate_qe_xc('PBE'),
            'pbe',
        )
    
    def test_validate_qe_version_accepts_reference_version(self):
        result = validate_qe_version(
            (7, 2)
        )

        self.assertEqual(
            result,
            (7, 2),
        )
    
    def test_validate_qe_version_accepts_patch_version(self):
        result = validate_qe_version(
            (7, 2, 1)
        )

        self.assertEqual(
            result,
            (7, 2, 1),
        )

    def test_validate_qe_version_rejects_older_version(self):
        with self.assertRaisesRegex(
            ValueError,
            'requires Quantum ESPRESSO 7.2 or newer',
        ):
            validate_qe_version(
                (7, 1)
            )

    def test_validate_qe_version_rejects_missing_version(self):
        with self.assertRaisesRegex(
            ValueError,
            'version could not be detected',
        ):
            validate_qe_version(
                None
            )

    def test_run_nscf_requires_ground_state(self):
        atoms = bulk(
            'Si',
            'diamond',
            a=5.43,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            with self.assertRaisesRegex(
                FileNotFoundError,
                'valid QE ground-state result',
            ):
                run_nscf(
                    atoms=atoms,
                    input_file=tmpdir / 'nscf.in',
                    output_file=tmpdir / 'nscf.out',
                    state_dir=tmpdir / 'state',
                    pseudopotentials={
                        'Si': 'Si.upf',
                    },
                    pseudo_dir='/tmp/pseudos',
                    cutoff_ev=400.0,
                )

    def test_parse_dos_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# E (eV) DOS(E) Int DOS(E)\n"
                "-1.0000  0.1000  0.0100\n"
                " 0.0000  0.5000  0.2000\n"
                " 1.0000  0.2500  0.6000\n",
                encoding='utf-8',
            )

            result = parse_dos_output(
                dos_file
            )

        self.assertEqual(
            result['npoints'],
            3,
        )

        self.assertEqual(
            result['energies_ev'],
            [-1.0, 0.0, 1.0],
        )

        self.assertEqual(
            result['dos'],
            [0.1, 0.5, 0.25],
        )

        self.assertEqual(
            result['integrated_dos'],
            [0.01, 0.2, 0.6],
        )

    def test_parse_dos_output_supports_fortran_exponents(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# E DOS IntDOS\n"
                "1.000D+00 2.500D-01 6.000D-01\n",
                encoding='utf-8',
            )

            result = parse_dos_output(
                dos_file
            )

        self.assertEqual(
            result['energies_ev'],
            [1.0],
        )

        self.assertEqual(
            result['dos'],
            [0.25],
        )

    def test_parse_dos_output_rejects_missing_data(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dos_file = (
                Path(tmpdir)
                / 'nanoworks.dos'
            )

            dos_file.write_text(
                "# no DOS data\n",
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'No DOS data could be parsed',
            ):
                parse_dos_output(
                    dos_file
                )

def test_render_projwfc_input(self):
    text = render_projwfc_input(
        prefix='nanoworks',
        outdir='/tmp/qe-state',
        filpdos='/tmp/gaas-pdos',
        emin=1.0,
        emax=10.0,
        delta_e=0.02,
    )

    self.assertIn(
        '&PROJWFC',
        text,
    )

    self.assertIn(
        "prefix = 'nanoworks'",
        text,
    )

    self.assertIn(
        "outdir = '/tmp/qe-state'",
        text,
    )

    self.assertIn(
        "filpdos = '/tmp/gaas-pdos'",
        text,
    )

    self.assertIn(
        'Emin = 1',
        text,
    )

    self.assertIn(
        'Emax = 10',
        text,
    )

    self.assertIn(
        'DeltaE = 0.02',
        text,
    )

    self.assertIn(
        'lsym = .true.',
        text,
    )

    self.assertNotIn(
        'degauss',
        text,
    )
    
def test_render_projwfc_input_rejects_invalid_energy_range(self):
    with self.assertRaisesRegex(
        ValueError,
        'Emax must be greater than Emin',
    ):
        render_projwfc_input(
            emin=5.0,
            emax=-5.0,
        )

def test_run_projwfc_requires_qe_state(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        with self.assertRaisesRegex(
            FileNotFoundError,
            'valid QE electronic state',
        ):
            run_projwfc(
                input_file=tmpdir / 'pdos.in',
                output_file=tmpdir / 'pdos.out',
                state_dir=tmpdir / 'state',
                pdos_prefix=tmpdir / 'pdos',
            )

def test_parse_projwfc_pdos_p_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#2(As)_wfc#3(p)'
        )

        pdos_file.write_text(
            "# E (eV) ldos(E) pdos(E) pdos(E) pdos(E)\n"
            "2.048 0.186E-01 0.621E-02 0.622E-02 0.623E-02\n"
            "2.078 0.175E-01 0.583E-02 0.584E-02 0.585E-02\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertEqual(
        result['atom_index'],
        2,
    )

    self.assertEqual(
        result['symbol'],
        'As',
    )

    self.assertEqual(
        result['wfc_index'],
        3,
    )

    self.assertEqual(
        result['orbital'],
        'p',
    )

    self.assertEqual(
        result['components']['pz'],
        [0.00621, 0.00583],
    )

    self.assertEqual(
        result['components']['px'],
        [0.00622, 0.00584],
    )

    self.assertEqual(
        result['components']['py'],
        [0.00623, 0.00585],
    )

def test_parse_projwfc_pdos_d_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#1(Ga)_wfc#1(d)'
        )

        pdos_file.write_text(
            "# E ldos d1 d2 d3 d4 d5\n"
            "2.048 0.198E-02 "
            "0.654E-03 0.223E-03 0.224E-03 "
            "0.438E-03 0.439E-03\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertEqual(
        result['orbital'],
        'd',
    )

    self.assertEqual(
        result['components']['d3z2_r2'],
        [0.000654],
    )

    self.assertEqual(
        result['components']['dxz'],
        [0.000223],
    )

    self.assertEqual(
        result['components']['dyz'],
        [0.000224],
    )

    self.assertEqual(
        result['components']['dx2_y2'],
        [0.000438],
    )

    self.assertEqual(
        result['components']['dxy'],
        [0.000439],
    )

def test_parse_projwfc_pdos_s_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        pdos_file = (
            Path(tmpdir)
            / 'test.pdos_atm#1(Ga)_wfc#2(s)'
        )

        pdos_file.write_text(
            "# E ldos pdos\n"
            "2.048 0.842E-01 0.842E-01\n",
            encoding='utf-8',
        )

        result = parse_projwfc_pdos_file(
            pdos_file
        )

    self.assertEqual(
        result['components']['s'],
        [0.0842],
    )

def test_aggregate_projwfc_pdos(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        prefix = (
            tmpdir
            / 'nanoworks-pdos'
        )

        ga_s = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#2(s)'
        )

        ga_p = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#3(p)'
        )

        as_p = Path(
            str(prefix)
            + '.pdos_atm#2(As)_wfc#3(p)'
        )

        ga_s.write_text(
            "# E ldos s\n"
            "1.0 0.10 0.10\n"
            "2.0 0.20 0.20\n",
            encoding='utf-8',
        )

        ga_p.write_text(
            "# E ldos pz px py\n"
            "1.0 0.60 0.10 0.20 0.30\n"
            "2.0 0.90 0.20 0.30 0.40\n",
            encoding='utf-8',
        )

        as_p.write_text(
            "# E ldos pz px py\n"
            "1.0 0.30 0.05 0.10 0.15\n"
            "2.0 0.60 0.10 0.20 0.30\n",
            encoding='utf-8',
        )

        result = aggregate_projwfc_pdos(
            prefix
        )

    self.assertEqual(
        result['energies_ev'],
        [1.0, 2.0],
    )

    self.assertEqual(
        result['s_total'],
        [0.1, 0.2],
    )

    self.assertEqual(
        result['p_total'],
        [0.9, 1.5],
    )

    self.assertAlmostEqual(
        result['total'][0],
        1.0,
    )

    self.assertAlmostEqual(
        result['total'][1],
        1.7,
    )

    self.assertEqual(
        result['pz'],
        [0.15, 0.30],
    )

    self.assertEqual(
        result['px'],
        [0.30, 0.50],
    )

    self.assertEqual(
        result['py'],
        [0.45, 0.70],
    )

    self.assertEqual(
        result['d_total'],
        [0.0, 0.0],
    )

    self.assertEqual(
        result['f_total'],
        [0.0, 0.0],
    )

def test_aggregate_projwfc_pdos_rejects_mismatched_energy_grid(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(
            tmpdir
        )

        prefix = (
            tmpdir
            / 'nanoworks-pdos'
        )

        first = Path(
            str(prefix)
            + '.pdos_atm#1(Ga)_wfc#2(s)'
        )

        second = Path(
            str(prefix)
            + '.pdos_atm#2(As)_wfc#2(s)'
        )

        first.write_text(
            "# E ldos s\n"
            "1.0 0.10 0.10\n"
            "2.0 0.20 0.20\n",
            encoding='utf-8',
        )

        second.write_text(
            "# E ldos s\n"
            "1.0 0.10 0.10\n"
            "2.1 0.20 0.20\n",
            encoding='utf-8',
        )

        with self.assertRaisesRegex(
            ValueError,
            'same energy grid',
        ):
            aggregate_projwfc_pdos(
                prefix
            )

if __name__ == '__main__':
    unittest.main()
