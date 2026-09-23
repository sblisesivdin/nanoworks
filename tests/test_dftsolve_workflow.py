import tempfile
import unittest
import sys
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from ase import Atoms
from ase.io import write

with patch.object(
    sys,
    'argv',
    [sys.argv[0]],
):
    from nanoworks.dftsolve import (
        DFTConfig,
        GPAW_STAGE_GROUP_ENV,
        build_gpaw_process_command,
        check_dft_configuration,
        dftsolve as DFTSolver,
        format_dft_preflight_json,
        format_dft_preflight_report,
        load_slurm_profile,
        main,
        prepare_qe_dry_run,
        release_stage_resources,
        required_dft_executables,
        run_calculation_stages,
        run_gpaw_stage_processes,
        should_split_gpaw_optical,
        write_qe_slurm_script,
    )


class TestDFTSolveWorkflow(unittest.TestCase):

    def test_required_qe_executables_follow_selected_stages(self):
        config = DFTConfig(
            Engine='QE',
            Ground_calc=True,
            DOS_calc=True,
            Band_calc=True,
            Projected_band_plot=True,
            Density_calc=True,
            Phonon_calc=True,
            Optical_calc=True,
        )

        self.assertEqual(
            set(required_dft_executables(config)),
            {
                'bands.x',
                'dos.x',
                'epsilon.x',
                'matdyn.x',
                'ph.x',
                'pp.x',
                'projwfc.x',
                'pw.x',
                'q2r.x',
            },
        )

    def test_qe_preflight_accepts_complete_optical_workflow(self):
        config = DFTConfig(
            Engine='QE',
            Ground_calc=True,
            Optical_calc=True,
            bulk_configuration=Atoms(
                'Si2',
                cell=[5.4, 5.4, 5.4],
                pbc=True,
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.shutil.which',
                side_effect=lambda name: f'/usr/bin/{name}',
            ),
            patch(
                'nanoworks.dftsolve.get_qe_pseudo_dir',
                return_value=Path('/pseudos'),
            ),
            patch(
                'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                return_value={
                    'Si': 'Si.upf',
                },
            ),
        ):
            report = check_dft_configuration(
                config,
                struct='silicon',
                parallel_cores=2,
            )

        self.assertTrue(
            report['ok']
        )
        self.assertEqual(
            report['errors'],
            [],
        )
        self.assertIn(
            '[OK] executable:epsilon.x: /usr/bin/epsilon.x',
            format_dft_preflight_report(report),
        )

    def test_qe_preflight_reports_missing_executable(self):
        config = DFTConfig(
            Engine='QE',
            Ground_calc=True,
            Optical_calc=True,
            bulk_configuration=Atoms(
                'Si2',
                cell=[5.4, 5.4, 5.4],
                pbc=True,
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.shutil.which',
                side_effect=lambda name: (
                    None
                    if name == 'epsilon.x'
                    else f'/usr/bin/{name}'
                ),
            ),
            patch(
                'nanoworks.dftsolve.get_qe_pseudo_dir',
                return_value=Path('/pseudos'),
            ),
            patch(
                'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                return_value={
                    'Si': 'Si.upf',
                },
            ),
        ):
            report = check_dft_configuration(
                config,
                struct='silicon',
            )

        self.assertFalse(
            report['ok']
        )
        self.assertIn(
            'epsilon.x was not found in PATH.',
            [
                error['detail']
                for error in report['errors']
            ],
        )

        payload = json.loads(
            format_dft_preflight_json(report)
        )
        self.assertEqual(
            payload['schema_version'],
            1,
        )
        self.assertFalse(
            payload['ok']
        )
        self.assertEqual(
            payload['error_count'],
            len(payload['errors']),
        )
        self.assertIsInstance(
            payload['stages'],
            list,
        )

    def test_qe_preflight_rejects_unsupported_elastic_stage(self):
        config = DFTConfig(
            Engine='QE',
            Ground_calc=True,
            Elastic_calc=True,
            bulk_configuration=Atoms(
                'Si2',
                cell=[5.4, 5.4, 5.4],
                pbc=True,
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.shutil.which',
                return_value='/usr/bin/pw.x',
            ),
            patch(
                'nanoworks.dftsolve.get_qe_pseudo_dir',
                return_value=Path('/pseudos'),
            ),
            patch(
                'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                return_value={
                    'Si': 'Si.upf',
                },
            ),
        ):
            report = check_dft_configuration(
                config,
                struct='silicon',
            )

        self.assertFalse(
            report['ok']
        )
        self.assertIn(
            'elastic',
            [
                error['name']
                for error in report['errors']
            ],
        )

    def test_qe_preflight_requires_saved_state_when_ground_is_skipped(self):
        config = DFTConfig(
            Engine='QE',
            Ground_calc=False,
            Density_calc=True,
            bulk_configuration=Atoms(
                'Si2',
                cell=[5.4, 5.4, 5.4],
                pbc=True,
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.shutil.which',
                return_value='/usr/bin/pp.x',
            ),
            patch(
                'nanoworks.engine.qe.has_qe_state',
                return_value=False,
            ),
        ):
            report = check_dft_configuration(
                config,
                struct='silicon',
            )

        self.assertFalse(
            report['ok']
        )
        self.assertIn(
            'ground-state',
            [
                error['name']
                for error in report['errors']
            ],
        )

    def test_check_cli_does_not_create_output_or_run_calculations(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'preflight_input.py'
            geometry_file = tmpdir / 'silicon.cif'
            output_dir = tmpdir / 'check-results'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Ground_calc = True\n"
                "Outdirname = 'check-results'\n",
                encoding='utf-8',
            )
            write(
                geometry_file,
                Atoms(
                    'Si2',
                    scaled_positions=[
                        (0.0, 0.0, 0.0),
                        (0.25, 0.25, 0.25),
                    ],
                    cell=[5.4, 5.4, 5.4],
                    pbc=True,
                ),
            )
            report = {
                'ok': True,
                'engine': 'QE',
                'stages': ('ground',),
                'checks': [],
                'errors': [],
            }

            with (
                patch.object(
                    sys,
                    'argv',
                    [
                        'dftsolve',
                        '--check',
                        '-i',
                        str(input_file),
                        '-g',
                        str(geometry_file),
                    ],
                ),
                patch(
                    'nanoworks.dftsolve.check_dft_configuration',
                    return_value=report,
                ) as check,
                patch(
                    'nanoworks.dftsolve.format_dft_preflight_report',
                    return_value='Result: READY',
                ),
                patch(
                    'nanoworks.dftsolve.parprint',
                ) as output,
            ):
                return_code = main()

            self.assertEqual(
                return_code,
                0,
            )
            self.assertFalse(
                output_dir.exists()
            )
            check.assert_called_once()
            output.assert_any_call(
                'Result: READY'
            )

    def test_check_json_cli_uses_machine_readable_formatter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'preflight_input.py'
            geometry_file = tmpdir / 'silicon.cif'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Ground_calc = True\n",
                encoding='utf-8',
            )
            write(
                geometry_file,
                Atoms(
                    'Si2',
                    scaled_positions=[
                        (0.0, 0.0, 0.0),
                        (0.25, 0.25, 0.25),
                    ],
                    cell=[5.4, 5.4, 5.4],
                    pbc=True,
                ),
            )
            report = {
                'ok': False,
                'engine': 'QE',
                'stages': ('ground',),
                'checks': [],
                'errors': [{
                    'status': 'error',
                    'name': 'executable:pw.x',
                    'detail': 'pw.x was not found in PATH.',
                }],
            }

            with (
                patch.object(
                    sys,
                    'argv',
                    [
                        'dftsolve',
                        '--check',
                        '--json',
                        '-i',
                        str(input_file),
                        '-g',
                        str(geometry_file),
                    ],
                ),
                patch(
                    'nanoworks.dftsolve.check_dft_configuration',
                    return_value=report,
                ),
                patch(
                    'nanoworks.dftsolve.parprint',
                ) as output,
            ):
                return_code = main()

            self.assertEqual(
                return_code,
                2,
            )
            rendered = output.call_args.args[0]
            payload = json.loads(rendered)
            self.assertEqual(
                payload['error_count'],
                1,
            )
            self.assertEqual(
                payload['stages'],
                ['ground'],
            )

    def test_json_cli_requires_check_mode(self):
        with (
            patch.object(
                sys,
                'argv',
                [
                    'dftsolve',
                    '--json',
                ],
            ),
            patch(
                'nanoworks.dftsolve.parprint',
            ) as output,
        ):
            return_code = main()

        self.assertEqual(
            return_code,
            2,
        )
        output.assert_called_once_with(
            'ERROR: --json requires --check.'
        )

    def test_qe_dry_run_writes_inputs_plan_and_script_without_execution(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            struct = Path(tmpdir) / 'dry-run' / 'silicon'
            config = DFTConfig(
                Engine='QE',
                Ground_calc=True,
                DOS_calc=True,
                Band_calc=True,
                Band_path='GXG',
                Projected_band_plot=True,
                Density_calc=True,
                Phonon_calc=True,
                Phonon_path='GXG',
                Optical_calc=True,
                bulk_configuration=Atoms(
                    'Si2',
                    scaled_positions=[
                        (0.0, 0.0, 0.0),
                        (0.25, 0.25, 0.25),
                    ],
                    cell=[5.4, 5.4, 5.4],
                    pbc=True,
                ),
            )

            with (
                patch(
                    'nanoworks.dftsolve.get_qe_pseudo_dir',
                    return_value=Path('/pseudos'),
                ),
                patch(
                    'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                    return_value={
                        'Si': 'Si.upf',
                    },
                ),
                patch(
                    'nanoworks.dftsolve.shutil.which',
                    return_value=None,
                ),
                patch(
                    'nanoworks.engine.qe.subprocess.run',
                ) as execute,
            ):
                plan = prepare_qe_dry_run(
                    config,
                    struct=struct,
                    parallel_cores=4,
                )

            execute.assert_not_called()
            self.assertEqual(
                plan['schema_version'],
                1,
            )
            self.assertEqual(
                plan['parallel_cores'],
                4,
            )
            self.assertFalse(
                plan['spin_polarized']
            )
            job_ids = {
                job['id']
                for job in plan['jobs']
            }
            self.assertEqual(
                job_ids,
                {
                    'ground',
                    'dos-nscf',
                    'dos-total',
                    'dos-projected',
                    'band',
                    'band-projections',
                    'density-pseudo-total',
                    'phonon-grid',
                    'phonon-force-constants',
                    'phonon-band',
                    'phonon-dos',
                    'optical-nscf',
                    'optical-epsilon',
                },
            )

            for job in plan['jobs']:
                self.assertTrue(
                    Path(job['input_file']).is_file()
                )
                self.assertEqual(
                    job['command'][:3],
                    ['mpiexec', '-np', '4'],
                )

            plan_file = Path(plan['plan_file'])
            script_file = Path(plan['script_file'])
            self.assertTrue(
                plan_file.is_file()
            )
            self.assertTrue(
                script_file.is_file()
            )
            self.assertTrue(
                script_file.stat().st_mode & 0o111
            )
            syntax_check = subprocess.run(
                [
                    'bash',
                    '-n',
                    str(script_file),
                ],
                capture_output=True,
                check=False,
                text=True,
            )
            self.assertEqual(
                syntax_check.returncode,
                0,
                syntax_check.stderr,
            )
            stored_plan = json.loads(
                plan_file.read_text(encoding='utf-8')
            )
            self.assertEqual(
                stored_plan['jobs'],
                plan['jobs'],
            )
            self.assertIn(
                "calculation = 'scf'",
                Path(
                    plan['jobs'][0]['input_file']
                ).read_text(encoding='utf-8'),
            )

            slurm_script = write_qe_slurm_script(
                plan,
                wall_time='2-12:30:00',
                memory='64G',
                partition='compute',
                account='project123',
                qos='normal',
                modules=[
                    'gcc/13.2',
                    'quantum-espresso/7.3.1',
                ],
                job_name='Si combined workflow',
                profile_file='/profiles/truba.json',
            )
            slurm_text = slurm_script.read_text(
                encoding='utf-8'
            )
            self.assertIn(
                '#SBATCH --job-name=Si-combined-workflow',
                slurm_text,
            )
            self.assertIn(
                '#SBATCH --ntasks=4',
                slurm_text,
            )
            self.assertIn(
                '#SBATCH --time=2-12:30:00',
                slurm_text,
            )
            self.assertIn(
                '#SBATCH --mem=64G',
                slurm_text,
            )
            self.assertIn(
                '#SBATCH --partition=compute',
                slurm_text,
            )
            self.assertIn(
                '#SBATCH --account=project123',
                slurm_text,
            )
            self.assertIn(
                '#SBATCH --qos=normal',
                slurm_text,
            )
            self.assertIn(
                'module load gcc/13.2',
                slurm_text,
            )
            self.assertIn(
                'module load quantum-espresso/7.3.1',
                slurm_text,
            )
            self.assertEqual(
                slurm_text.count('srun -n 4'),
                len(plan['jobs']),
            )
            slurm_syntax = subprocess.run(
                [
                    'bash',
                    '-n',
                    str(slurm_script),
                ],
                capture_output=True,
                check=False,
                text=True,
            )
            self.assertEqual(
                slurm_syntax.returncode,
                0,
                slurm_syntax.stderr,
            )
            stored_slurm_plan = json.loads(
                plan_file.read_text(encoding='utf-8')
            )
            self.assertEqual(
                stored_slurm_plan['scheduler'],
                'slurm',
            )
            self.assertEqual(
                stored_slurm_plan['slurm']['ntasks'],
                4,
            )
            self.assertEqual(
                stored_slurm_plan['slurm']['profile_file'],
                '/profiles/truba.json',
            )

    def test_qe_spin_dry_run_writes_combined_spin_workflow(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            struct = Path(tmpdir) / 'spin-dry-run' / 'iron'
            config = DFTConfig(
                Engine='QE',
                Ground_calc=True,
                DOS_calc=True,
                Band_calc=True,
                Band_path='GX',
                Projected_band_plot=True,
                Density_calc=True,
                Spin_calc=True,
                Magmom_per_atom=2.0,
                bulk_configuration=Atoms(
                    'Fe',
                    cell=[2.87, 2.87, 2.87],
                    pbc=True,
                ),
            )

            with (
                patch(
                    'nanoworks.dftsolve.get_qe_pseudo_dir',
                    return_value=Path('/pseudos'),
                ),
                patch(
                    'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                    return_value={
                        'Fe': 'Fe.upf',
                    },
                ),
                patch(
                    'nanoworks.engine.qe.read_upf_z_valence',
                    return_value=8.0,
                ),
                patch(
                    'nanoworks.dftsolve.shutil.which',
                    return_value=None,
                ),
            ):
                plan = prepare_qe_dry_run(
                    config,
                    struct=struct,
                    parallel_cores=2,
                )

            self.assertTrue(
                plan['spin_polarized']
            )

            jobs = {
                job['id']: job
                for job in plan['jobs']
            }
            self.assertEqual(
                set(jobs),
                {
                    'ground',
                    'dos-nscf',
                    'dos-total',
                    'dos-projected',
                    'band',
                    'band-projections',
                    'density-pseudo-total',
                    'density-pseudo-up',
                    'density-pseudo-down',
                    'density-spin-density',
                },
            )

            for job_id in (
                'ground',
                'dos-nscf',
                'band',
            ):
                input_text = Path(
                    jobs[job_id]['input_file']
                ).read_text(
                    encoding='utf-8'
                )
                self.assertIn(
                    'nspin = 2',
                    input_text,
                )
                self.assertIn(
                    'starting_magnetization(1) = 0.25',
                    input_text,
                )

            density_components = {
                'density-pseudo-total': 0,
                'density-pseudo-up': 1,
                'density-pseudo-down': 2,
            }

            for job_id, spin_component in density_components.items():
                input_text = Path(
                    jobs[job_id]['input_file']
                ).read_text(
                    encoding='utf-8'
                )
                self.assertIn(
                    f'spin_component = {spin_component}',
                    input_text,
                )

            spin_density_text = Path(
                jobs[
                    'density-spin-density'
                ][
                    'input_file'
                ]
            ).read_text(
                encoding='utf-8'
            )
            self.assertIn(
                'plot_num = 6',
                spin_density_text,
            )
            self.assertNotIn(
                'spin_component',
                spin_density_text,
            )

            stored_plan = json.loads(
                Path(plan['plan_file']).read_text(
                    encoding='utf-8'
                )
            )
            self.assertTrue(
                stored_plan['spin_polarized']
            )

    def test_dry_run_cli_stops_before_calculation_stages(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'dry_run_input.py'
            geometry_file = tmpdir / 'silicon.cif'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Ground_calc = True\n",
                encoding='utf-8',
            )
            write(
                geometry_file,
                Atoms(
                    'Si2',
                    scaled_positions=[
                        (0.0, 0.0, 0.0),
                        (0.25, 0.25, 0.25),
                    ],
                    cell=[5.4, 5.4, 5.4],
                    pbc=True,
                ),
            )
            report = {
                'ok': True,
                'engine': 'QE',
                'stages': ('ground',),
                'checks': [],
                'errors': [],
            }
            plan = {
                'stages': ['ground'],
                'jobs': [{}],
                'plan_file': '/tmp/plan.json',
                'script_file': '/tmp/run.sh',
            }

            with (
                patch.object(
                    sys,
                    'argv',
                    [
                        'dftsolve',
                        '--dry-run',
                        '--scheduler',
                        'slurm',
                        '--slurm-time',
                        '12:00:00',
                        '--slurm-account',
                        'project123',
                        '--slurm-qos',
                        'normal',
                        '--slurm-module',
                        'gcc/13.2',
                        '--slurm-module',
                        'quantum-espresso/7.3.1',
                        '-i',
                        str(input_file),
                        '-g',
                        str(geometry_file),
                    ],
                ),
                patch(
                    'nanoworks.dftsolve.check_dft_configuration',
                    return_value=report,
                ) as check,
                patch(
                    'nanoworks.dftsolve.prepare_qe_dry_run',
                    return_value=plan,
                ) as prepare,
                patch(
                    'nanoworks.dftsolve.write_qe_slurm_script',
                ) as write_slurm,
                patch(
                    'nanoworks.dftsolve.run_calculation_stages',
                ) as execute,
                patch(
                    'nanoworks.dftsolve.parprint',
                ),
            ):
                return_code = main()

            self.assertEqual(
                return_code,
                0,
            )
            check.assert_called_once()
            prepare.assert_called_once()
            write_slurm.assert_called_once_with(
                plan,
                wall_time='12:00:00',
                memory=None,
                partition=None,
                account='project123',
                qos='normal',
                modules=[
                    'gcc/13.2',
                    'quantum-espresso/7.3.1',
                ],
                job_name=None,
                profile_file=None,
            )
            execute.assert_not_called()

    def test_slurm_scheduler_requires_dry_run(self):
        with (
            patch.object(
                sys,
                'argv',
                [
                    'dftsolve',
                    '--scheduler',
                    'slurm',
                ],
            ),
            patch(
                'nanoworks.dftsolve.parprint',
            ) as output,
        ):
            return_code = main()

        self.assertEqual(
            return_code,
            2,
        )
        output.assert_called_once_with(
            'ERROR: --scheduler requires --dry-run.'
        )

    def test_cluster_profile_implies_slurm_and_cli_overrides_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            input_file = tmpdir / 'dry_run_input.py'
            geometry_file = tmpdir / 'silicon.cif'
            profile_file = tmpdir / 'truba.json'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Ground_calc = True\n",
                encoding='utf-8',
            )
            write(
                geometry_file,
                Atoms(
                    'Si2',
                    scaled_positions=[
                        (0.0, 0.0, 0.0),
                        (0.25, 0.25, 0.25),
                    ],
                    cell=[5.4, 5.4, 5.4],
                    pbc=True,
                ),
            )
            profile_file.write_text(
                json.dumps({
                    'slurm': {
                        'time': '2-00:00:00',
                        'memory': '64G',
                        'partition': 'compute',
                        'account': 'project123',
                        'modules': ['qe/profile'],
                    },
                }),
                encoding='utf-8',
            )
            report = {
                'ok': True,
                'engine': 'QE',
                'stages': ('ground',),
                'checks': [],
                'errors': [],
            }
            plan = {
                'stages': ['ground'],
                'jobs': [{}],
                'plan_file': '/tmp/plan.json',
                'script_file': '/tmp/run.sh',
            }

            with (
                patch.object(
                    sys,
                    'argv',
                    [
                        'dftsolve',
                        '--dry-run',
                        '--cluster-profile',
                        str(profile_file),
                        '--slurm-time',
                        '06:00:00',
                        '--slurm-module',
                        'qe/override',
                        '-i',
                        str(input_file),
                        '-g',
                        str(geometry_file),
                    ],
                ),
                patch(
                    'nanoworks.dftsolve.check_dft_configuration',
                    return_value=report,
                ),
                patch(
                    'nanoworks.dftsolve.prepare_qe_dry_run',
                    return_value=plan,
                ),
                patch(
                    'nanoworks.dftsolve.write_qe_slurm_script',
                ) as write_slurm,
                patch(
                    'nanoworks.dftsolve.parprint',
                ),
            ):
                return_code = main()

            self.assertEqual(
                return_code,
                0,
            )
            write_slurm.assert_called_once_with(
                plan,
                wall_time='06:00:00',
                memory='64G',
                partition='compute',
                account='project123',
                qos=None,
                modules=['qe/override'],
                job_name=None,
                profile_file=str(profile_file.resolve()),
            )

    def test_slurm_writer_rejects_invalid_wall_time(self):
        plan = {
            'engine': 'QE',
            'dry_run': True,
            'parallel_cores': 4,
            'jobs': [{
                'executable': 'pw.x',
                'input_file': '/tmp/input.in',
                'output_file': '/tmp/output.out',
                'working_directory': None,
            }],
            'plan_file': '/tmp/plan.json',
            'script_file': '/tmp/run.sh',
        }

        with self.assertRaisesRegex(
            ValueError,
            'Slurm time',
        ):
            write_qe_slurm_script(
                plan,
                wall_time='12:90:00',
            )

    def test_load_slurm_profile_accepts_json_settings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            profile_file = Path(tmpdir) / 'truba.json'
            profile_file.write_text(
                json.dumps({
                    'slurm': {
                        'time': '2-00:00:00',
                        'memory': '64G',
                        'partition': 'compute',
                        'account': 'project123',
                        'qos': 'normal',
                        'modules': [
                            'gcc/13.2',
                            'quantum-espresso/7.3.1',
                        ],
                        'job_name': 'nanoworks-qe',
                    },
                }),
                encoding='utf-8',
            )

            profile = load_slurm_profile(
                profile_file
            )

        self.assertEqual(
            profile['time'],
            '2-00:00:00',
        )
        self.assertEqual(
            profile['modules'],
            [
                'gcc/13.2',
                'quantum-espresso/7.3.1',
            ],
        )
        self.assertEqual(
            profile['profile_file'],
            str(profile_file.resolve()),
        )

    def test_load_slurm_profile_resolves_named_user_profile(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            profile_file = (
                Path(tmpdir)
                / '.config'
                / 'nanoworks'
                / 'clusters'
                / 'truba.json'
            )
            profile_file.parent.mkdir(
                parents=True
            )
            profile_file.write_text(
                json.dumps({
                    'slurm': {
                        'account': 'project123',
                    },
                }),
                encoding='utf-8',
            )

            with patch(
                'nanoworks.dftsolve.Path.home',
                return_value=Path(tmpdir),
            ):
                profile = load_slurm_profile(
                    'truba'
                )

        self.assertEqual(
            profile['account'],
            'project123',
        )

    def test_load_slurm_profile_rejects_unknown_setting(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            profile_file = Path(tmpdir) / 'invalid.json'
            profile_file.write_text(
                json.dumps({
                    'slurm': {
                        'shell_command': 'unsafe',
                    },
                }),
                encoding='utf-8',
            )

            with self.assertRaisesRegex(
                ValueError,
                'shell_command',
            ):
                load_slurm_profile(
                    profile_file
                )

    def test_slurm_writer_rejects_unsafe_module_name(self):
        plan = {
            'engine': 'QE',
            'dry_run': True,
            'parallel_cores': 4,
            'jobs': [{
                'executable': 'pw.x',
                'input_file': '/tmp/input.in',
                'output_file': '/tmp/output.out',
                'working_directory': None,
            }],
            'plan_file': '/tmp/plan.json',
            'script_file': '/tmp/run.sh',
        }

        with self.assertRaisesRegex(
            ValueError,
            'module name',
        ):
            write_qe_slurm_script(
                plan,
                modules=['qe/7.3; touch unsafe'],
            )

    def test_qe_hybrid_dry_run_writes_stage_specific_scf_inputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            struct = Path(tmpdir) / 'hybrid-dry-run' / 'silicon'
            config = DFTConfig(
                Engine='QE',
                XC_calc='HSE06',
                Ground_calc=True,
                DOS_calc=True,
                Band_calc=True,
                Band_path='GXG',
                Band_npoints=5,
                Projected_band_plot=True,
                Ground_kpts_x=2,
                Ground_kpts_y=2,
                Ground_kpts_z=2,
                bulk_configuration=Atoms(
                    'Si2',
                    scaled_positions=[
                        (0.0, 0.0, 0.0),
                        (0.25, 0.25, 0.25),
                    ],
                    cell=[5.4, 5.4, 5.4],
                    pbc=True,
                ),
            )

            with (
                patch(
                    'nanoworks.dftsolve.get_qe_pseudo_dir',
                    return_value=Path('/pseudos'),
                ),
                patch(
                    'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                    return_value={
                        'Si': 'Si.upf',
                    },
                ),
                patch(
                    'nanoworks.engine.qe.subprocess.run',
                ) as execute,
            ):
                plan = prepare_qe_dry_run(
                    config,
                    struct=struct,
                )

            execute.assert_not_called()
            jobs = {
                job['id']: job
                for job in plan['jobs']
            }
            self.assertEqual(
                set(jobs),
                {
                    'ground',
                    'dos-hybrid-scf',
                    'dos-total',
                    'dos-projected',
                    'band-hybrid-scf',
                    'band-postprocess',
                    'band-projections',
                },
            )
            self.assertEqual(
                jobs['dos-total']['depends_on'],
                ['dos-hybrid-scf'],
            )
            self.assertEqual(
                jobs['band-postprocess']['depends_on'],
                ['band-hybrid-scf'],
            )
            self.assertEqual(
                jobs['band-projections']['depends_on'],
                ['band-postprocess'],
            )

            dos_scf = Path(
                jobs['dos-hybrid-scf']['input_file']
            ).read_text(encoding='utf-8')
            band_scf = Path(
                jobs['band-hybrid-scf']['input_file']
            ).read_text(encoding='utf-8')
            bands_postprocess = Path(
                jobs['band-postprocess']['input_file']
            ).read_text(encoding='utf-8')

            self.assertIn(
                "calculation = 'scf'",
                dos_scf,
            )
            self.assertNotIn(
                'ADDITIONAL_K_POINTS',
                dos_scf,
            )
            self.assertIn(
                'ADDITIONAL_K_POINTS crystal',
                band_scf,
            )
            self.assertIn(
                'nqx1 = 2',
                band_scf,
            )
            self.assertIn(
                '&BANDS',
                bands_postprocess,
            )
            self.assertTrue(
                jobs['band-hybrid-scf'][
                    'metadata'
                ][
                    'helper_count'
                ] > 0
            )
            self.assertEqual(
                jobs['band-postprocess'][
                    'metadata'
                ][
                    'band_indices'
                ],
                list(range(config.Band_npoints)),
            )
            hybrid_slurm = write_qe_slurm_script(
                plan,
                wall_time='08:00:00',
            ).read_text(encoding='utf-8')
            self.assertIn(
                'srun -n 1 bands.x -i',
                hybrid_slurm,
            )
            self.assertIn(
                'srun -n 1 projwfc.x -i',
                hybrid_slurm,
            )

    def test_opticalcalc_dispatches_to_gpaw(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Engine = 'GPAW'
        expected = object()
        solver._opticalcalc_gpaw = Mock(
            return_value=expected
        )

        result = solver.opticalcalc()

        self.assertIs(
            result,
            expected,
        )
        solver._opticalcalc_gpaw.assert_called_once_with()

    def test_opticalcalc_dispatches_to_qe(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Engine = 'QE'
        expected = object()
        solver._opticalcalc_qe = Mock(
            return_value=expected
        )

        result = solver.opticalcalc()

        self.assertIs(
            result,
            expected,
        )
        solver._opticalcalc_qe.assert_called_once_with()

    def test_qe_opticalcalc_runs_nscf_epsilon_and_writes_tables(self):
        optical_data = {
            'energies_ev': [0.0, 1.0],
            'directions': {
                direction: {
                    'epsilon_real': [1.0, 2.0],
                    'epsilon_imaginary': [0.0, 0.5],
                    'refractive_index': [1.0, 1.5],
                    'extinction_coefficient': [0.0, 0.2],
                    'absorption_cm_inverse': [0.0, 1000.0],
                    'reflectivity': [0.0, 0.05],
                }
                for direction in 'xyz'
            },
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            solver = object.__new__(
                DFTSolver
            )
            solver.Engine = 'QE'
            solver.Mode = 'PW'
            solver.SOC_calc = False
            solver.Opt_calc_type = 'RPA'
            solver.XC_calc = 'PBE'
            solver.struct = str(
                Path(tmpdir) / 'silicon'
            )
            solver.bulk_configuration = Atoms(
                'Si2',
                cell=[5.4, 5.4, 5.4],
                pbc=True,
            )
            solver.Gamma = False
            solver.Ground_gamma = None
            solver.Ground_kpts_density = None
            solver.Ground_kpts_x = 2
            solver.Ground_kpts_y = 2
            solver.Ground_kpts_z = 2
            solver.Opt_kpts_density = None
            solver.Opt_kpts_x = 4
            solver.Opt_kpts_y = 4
            solver.Opt_kpts_z = 4
            solver.Opt_gamma = None
            solver.Spin_calc = False
            solver.Cut_off_energy = 500.0
            solver.Total_charge = 0.0
            solver.Opt_num_of_bands = 16
            solver.Setup_params = None
            solver.XC_exx_fraction = None
            solver.XC_omega = None
            solver.Opt_FD_smearing = 0.05
            solver.Opt_eta = 0.1
            solver.Opt_BSE_min_en = 0.0
            solver.Opt_BSE_max_en = 10.0
            solver.Opt_BSE_num_of_data = 101
            solver.Opt_min_en = 0.0
            solver.Opt_max_en = 10.0
            solver.Opt_num_of_data = 101
            solver.Opt_shift_en = 0.2
            solver.parallel_cores = 2
            solver._generate_optical_figures = Mock()

            def write_tables(data, output_prefix):
                self.assertIs(
                    data,
                    optical_data,
                )
                output_files = {}

                for direction in 'xyz':
                    output_file = Path(
                        f"{output_prefix}-AllData_"
                        f"{direction}direction.dat"
                    )
                    output_file.write_text(
                        'header\n'
                        '0 1 0 1 0 0 0\n'
                        '1 2 0.5 1.5 0.2 1000 0.05\n',
                        encoding='utf-8',
                    )
                    output_files[direction] = output_file

                return output_files

            solver.engine = SimpleNamespace(
                validate_qe_xc=Mock(
                    return_value='pbe'
                ),
                has_qe_state=Mock(
                    return_value=True
                ),
                run_nscf=Mock(
                    return_value={'result': {}}
                ),
                run_epsilon=Mock(
                    return_value={
                        'optical_data': optical_data,
                    }
                ),
                write_epsilon_optical_data=Mock(
                    side_effect=write_tables
                ),
            )

            with (
                patch(
                    'nanoworks.dftsolve.get_qe_pseudo_dir',
                    return_value=Path('pseudos'),
                ),
                patch(
                    'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                    return_value={
                        'Si': 'Si.upf',
                    },
                ),
                patch(
                    'nanoworks.dftsolve.parprint',
                ),
            ):
                result = solver._opticalcalc_qe()

        nscf_call = solver.engine.run_nscf.call_args.kwargs
        self.assertEqual(
            nscf_call['kpoint_size'],
            (4, 4, 4),
        )
        self.assertTrue(
            nscf_call['nosym']
        )
        self.assertEqual(
            nscf_call['nbands'],
            16,
        )
        epsilon_call = solver.engine.run_epsilon.call_args.kwargs
        self.assertEqual(
            epsilon_call['wmax'],
            10.0,
        )
        self.assertEqual(
            epsilon_call['nw'],
            101,
        )
        self.assertEqual(
            epsilon_call['intersmear'],
            0.1,
        )
        self.assertEqual(
            result['epsilon']['optical_data'],
            optical_data,
        )
        self.assertEqual(
            solver._generate_optical_figures.call_count,
            3,
        )

    def test_qe_opticalcalc_rejects_bse(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Mode = 'PW'
        solver.SOC_calc = False
        solver.Opt_calc_type = 'BSE'

        with self.assertRaisesRegex(
            NotImplementedError,
            "Opt_calc_type = 'RPA' only",
        ):
            solver._opticalcalc_qe()

    def test_opticalcalc_rejects_unknown_engine(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Engine = 'UNKNOWN'

        with self.assertRaisesRegex(
            ValueError,
            'Unsupported optical engine: UNKNOWN',
        ):
            solver.opticalcalc()

    def test_run_calculation_stages_runs_ground_only_by_default(self):
        solver = Mock()
        config = SimpleNamespace()

        run_calculation_stages(
            solver,
            config,
        )

        solver.groundcalc.assert_called_once_with()
        solver.elasticcalc.assert_not_called()
        solver.doscalc.assert_not_called()
        solver.bandcalc.assert_not_called()
        solver.densitycalc.assert_not_called()
        solver.phononcalc.assert_not_called()
        solver.opticalcalc.assert_not_called()

    def test_run_calculation_stages_keeps_optical_last(self):
        calls = []
        solver = SimpleNamespace(
            groundcalc=lambda: calls.append('ground'),
            elasticcalc=lambda: calls.append('elastic'),
            doscalc=lambda: calls.append('dos'),
            bandcalc=lambda: calls.append('band'),
            densitycalc=lambda: calls.append('density'),
            phononcalc=lambda: calls.append('phonon'),
            opticalcalc=lambda: calls.append('optical'),
        )
        config = SimpleNamespace(
            Elastic_calc=True,
            DOS_calc=True,
            Band_calc=True,
            Density_calc=True,
            Phonon_calc=True,
            Optical_calc=True,
        )

        run_calculation_stages(
            solver,
            config,
        )

        self.assertEqual(
            calls,
            [
                'ground',
                'elastic',
                'dos',
                'band',
                'density',
                'phonon',
                'optical',
            ],
        )

    def test_run_calculation_stages_accepts_selected_stage_group(self):
        calls = []
        solver = SimpleNamespace(
            groundcalc=lambda: calls.append('ground'),
            doscalc=lambda: calls.append('dos'),
            opticalcalc=lambda: calls.append('optical'),
        )
        config = SimpleNamespace(
            DOS_calc=True,
            Optical_calc=True,
        )

        run_calculation_stages(
            solver,
            config,
            stages=('optical',),
        )

        self.assertEqual(
            calls,
            ['optical'],
        )

    def test_mixed_gpaw_optical_workflow_requires_process_split(self):
        self.assertTrue(
            should_split_gpaw_optical(
                SimpleNamespace(
                    Engine='GPAW',
                    Ground_calc=True,
                    Optical_calc=True,
                )
            )
        )
        self.assertFalse(
            should_split_gpaw_optical(
                SimpleNamespace(
                    Engine='GPAW',
                    Ground_calc=False,
                    Optical_calc=True,
                )
            )
        )
        self.assertFalse(
            should_split_gpaw_optical(
                SimpleNamespace(
                    Engine='QE',
                    Ground_calc=True,
                    Optical_calc=True,
                )
            )
        )

    def test_build_serial_gpaw_process_command_uses_fresh_python(self):
        command = build_gpaw_process_command(
            None,
            ['-i', 'input.py', '-g', 'geometry.cif'],
        )

        self.assertEqual(
            command[0],
            sys.executable,
        )
        self.assertTrue(
            command[1].endswith('nanoworks/dftsolve.py')
        )
        self.assertEqual(
            command[2:],
            ['-i', 'input.py', '-g', 'geometry.cif'],
        )

    def test_gpaw_stage_processes_are_sequential_and_stop_on_failure(self):
        completed = [
            SimpleNamespace(returncode=0),
            SimpleNamespace(returncode=7),
        ]

        with (
            patch(
                'nanoworks.dftsolve.build_gpaw_process_command',
                return_value=['python', 'dftsolve.py'],
            ),
            patch(
                'nanoworks.dftsolve.subprocess.run',
                side_effect=completed,
            ) as run,
            patch(
                'nanoworks.dftsolve.parprint',
            ) as output,
        ):
            return_code = run_gpaw_stage_processes(
                4,
                ['-i', 'input.py'],
            )

        self.assertEqual(
            return_code,
            7,
        )
        self.assertEqual(
            run.call_count,
            2,
        )
        self.assertEqual(
            run.call_args_list[0].kwargs['env'][GPAW_STAGE_GROUP_ENV],
            'electronic',
        )
        self.assertEqual(
            run.call_args_list[1].kwargs['env'][GPAW_STAGE_GROUP_ENV],
            'optical',
        )
        self.assertNotIn(
            'GPAW_MPI_BACKEND',
            run.call_args_list[0].kwargs['env'],
        )
        self.assertEqual(
            run.call_args_list[0].kwargs['env']['OMP_NUM_THREADS'],
            '1',
        )
        self.assertIn(
            'exit code 7',
            output.call_args.args[0],
        )

    def test_release_stage_resources_detaches_calculator_and_synchronizes(self):
        atoms = SimpleNamespace(
            calc=object(),
        )
        solver = SimpleNamespace(
            bulk_configuration=atoms,
        )

        with (
            patch(
                'nanoworks.dftsolve.gc.collect'
            ) as collect,
            patch(
                'nanoworks.dftsolve.world.barrier'
            ) as barrier,
        ):
            release_stage_resources(
                solver
            )

        self.assertIsNone(
            atoms.calc
        )
        collect.assert_called_once_with()
        barrier.assert_called_once_with()

    def test_run_calculation_stages_releases_resources_around_optical(self):
        atoms = SimpleNamespace(
            calc=object(),
        )
        calculator_seen_by_optical = []

        def run_optical():
            calculator_seen_by_optical.append(
                atoms.calc
            )
            atoms.calc = object()

        solver = SimpleNamespace(
            bulk_configuration=atoms,
            groundcalc=lambda: None,
            opticalcalc=run_optical,
        )
        config = SimpleNamespace(
            Optical_calc=True,
        )

        with (
            patch(
                'nanoworks.dftsolve.gc.collect'
            ) as collect,
            patch(
                'nanoworks.dftsolve.world.barrier'
            ) as barrier,
        ):
            run_calculation_stages(
                solver,
                config,
            )

        self.assertEqual(
            calculator_seen_by_optical,
            [None],
        )
        self.assertIsNone(
            atoms.calc
        )
        self.assertEqual(
            collect.call_count,
            2,
        )
        self.assertEqual(
            barrier.call_count,
            2,
        )

    def test_load_existing_final_structure(self):
        initial = Atoms(
            'Si',
            positions=[
                [0.0, 0.0, 0.0],
            ],
            cell=[4.0, 4.0, 4.0],
            pbc=True,
        )

        final = Atoms(
            'Si',
            positions=[
                [0.25, 0.25, 0.25],
            ],
            cell=[5.0, 5.0, 5.0],
            pbc=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            solver = object.__new__(
                DFTSolver
            )

            solver.struct = str(
                Path(tmpdir)
                / 'silicon'
            )
            solver.Engine = 'QE'
            solver.bulk_configuration = initial
            solver.config = SimpleNamespace(
                bulk_configuration=initial,
            )

            final_file = Path(
                solver.struct
                + '-GROUND-QE-Result-Final.cif'
            )

            write(
                final_file,
                final,
            )

            with patch(
                'nanoworks.dftsolve.parprint'
            ) as warning:
                loaded = (
                    solver
                    ._load_existing_final_structure()
                )

            self.assertTrue(
                loaded
            )
            self.assertAlmostEqual(
                solver.bulk_configuration.cell.lengths()[0],
                5.0,
            )
            self.assertIs(
                solver.config.bulk_configuration,
                solver.bulk_configuration,
            )
            self.assertIn(
                'WARNING:',
                warning.call_args.args[0],
            )
            self.assertIn(
                str(final_file),
                warning.call_args.args[0],
            )

    def test_load_existing_final_structure_when_missing(self):
        initial = Atoms(
            'Si',
            cell=[4.0, 4.0, 4.0],
            pbc=True,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            solver = object.__new__(
                DFTSolver
            )

            solver.struct = str(
                Path(tmpdir)
                / 'silicon'
            )
            solver.Engine = 'GPAW'
            solver.bulk_configuration = initial
            solver.config = SimpleNamespace(
                bulk_configuration=initial,
            )

            with patch(
                'nanoworks.dftsolve.parprint'
            ) as warning:
                loaded = (
                    solver
                    ._load_existing_final_structure()
                )

            self.assertFalse(
                loaded
            )
            self.assertIs(
                solver.bulk_configuration,
                initial,
            )
            warning.assert_not_called()

    def test_qe_engine_specific_defaults(self):
        config = DFTConfig(
            Engine='qe',
        )

        self.assertEqual(
            config.Engine,
            'QE',
        )
        self.assertEqual(
            config.XC_calc,
            'PBE',
        )
        self.assertEqual(
            config.DOS_occupation,
            'tetrahedra',
        )
        self.assertTrue(
            config.Fix_symmetry
        )
        self.assertIsNone(
            config.Mixer_type
        )
        self.assertIsNone(
            config.Phonon_PW_cutoff
        )
        self.assertIsNone(
            config.Phonon_kpts_x
        )
        self.assertIsNone(
            config.Phonon_kpts_y
        )
        self.assertIsNone(
            config.Phonon_kpts_z
        )
        self.assertEqual(
            config.Opt_calc_type,
            'RPA',
        )
        self.assertEqual(
            (
                config.Opt_min_en,
                config.Opt_max_en,
                config.Opt_num_of_data,
            ),
            (0.0, 20.0, 1001),
        )

    def test_gpaw_engine_specific_defaults(self):
        config = DFTConfig(
            Engine='GPAW',
            Mixer_type='custom-mixer',
        )

        self.assertEqual(
            config.XC_calc,
            'LDA',
        )
        self.assertIsNone(
            config.DOS_occupation
        )
        self.assertFalse(
            config.Fix_symmetry
        )
        self.assertEqual(
            config.Mixer_type,
            'custom-mixer',
        )
        self.assertEqual(
            config.Phonon_PW_cutoff,
            400,
        )
        self.assertEqual(
            config.Opt_calc_type,
            'BSE',
        )
        self.assertEqual(
            (
                config.Phonon_kpts_x,
                config.Phonon_kpts_y,
                config.Phonon_kpts_z,
            ),
            (3, 3, 3),
        )

    def test_explicit_values_override_engine_defaults(self):
        occupation = {
            'name': 'fermi-dirac',
            'width': 0.02,
        }

        config = DFTConfig(
            Engine='QE',
            XC_calc='LDA',
            DOS_occupation=occupation,
            Fix_symmetry=False,
            Phonon_PW_cutoff=500,
            Phonon_kpts_x=4,
            Phonon_kpts_y=5,
            Phonon_kpts_z=6,
            Opt_min_en=1.0,
            Opt_max_en=12.0,
            Opt_num_of_data=221,
        )

        self.assertEqual(
            config.XC_calc,
            'LDA',
        )
        self.assertIs(
            config.DOS_occupation,
            occupation,
        )
        self.assertFalse(
            config.Fix_symmetry
        )
        self.assertEqual(
            config.Phonon_PW_cutoff,
            500,
        )
        self.assertEqual(
            (
                config.Phonon_kpts_x,
                config.Phonon_kpts_y,
                config.Phonon_kpts_z,
            ),
            (4, 5, 6),
        )
        self.assertEqual(
            (
                config.Opt_min_en,
                config.Opt_max_en,
                config.Opt_num_of_data,
            ),
            (1.0, 12.0, 221),
        )

    def test_phononcalc_dispatches_to_gpaw(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Engine = 'GPAW'

        with patch.object(
            solver,
            '_phononcalc_gpaw',
            return_value='gpaw-phonons',
        ) as workflow:
            result = solver.phononcalc()

        workflow.assert_called_once_with()
        self.assertEqual(
            result,
            'gpaw-phonons',
        )

    def test_phononcalc_dispatches_to_qe(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Engine = 'QE'

        with patch.object(
            solver,
            '_phononcalc_qe',
            return_value='qe-phonons',
        ) as workflow:
            result = solver.phononcalc()

        workflow.assert_called_once_with()
        self.assertEqual(
            result,
            'qe-phonons',
        )

    def test_phononcalc_rejects_unknown_engine(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Engine = 'UNKNOWN'

        with self.assertRaisesRegex(
            ValueError,
            'Unsupported phonon engine',
        ):
            solver.phononcalc()

    def test_qe_phononcalc_runs_native_dfpt_chain(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            solver = object.__new__(
                DFTSolver
            )
            solver.struct = str(
                Path(tmpdir)
                / 'silicon'
            )
            solver.Engine = 'QE'
            solver.Mode = 'PW'
            solver.XC_calc = 'PBE'
            solver.Phonon_thermal_calc = True
            solver.Phonon_T_min = 0.0
            solver.Phonon_T_max = 600.0
            solver.Phonon_T_step = 20.0
            solver.Phonon_PW_cutoff = None
            solver.Phonon_kpts_x = None
            solver.Phonon_kpts_y = None
            solver.Phonon_kpts_z = None
            solver.Phonon_displacement = 1.0e-3
            solver.Phonon_supercell = (
                (2, 0, 0),
                (0, 3, 0),
                (0, 0, 4),
            )
            solver.Phonon_qpts_x = 12
            solver.Phonon_qpts_y = 13
            solver.Phonon_qpts_z = 14
            solver.Phonon_path = 'GX'
            solver.Phonon_npoints = 21
            solver.Phonon_acoustic_sum_rule = True
            solver.parallel_cores = 8
            solver.bulk_configuration = Atoms(
                'Si',
                cell=[5.4, 5.4, 5.4],
                pbc=True,
            )

            band_path = {
                'option': 'crystal',
                'kpoints': [
                    (0.0, 0.0, 0.0),
                    (0.5, 0.0, 0.0),
                ],
                'npoints': 2,
            }
            thermal_data = {
                'integrated_mode_weight': 5.9,
                'excluded_mode_weight': 0.1,
            }
            solver.engine = SimpleNamespace(
                THZ_PER_CM_MINUS_ONE=0.0299792458,
                validate_qe_xc=Mock(
                    return_value='pbe'
                ),
                has_qe_state=Mock(
                    return_value=True
                ),
                resolve_qe_phonon_qpoint_grid=Mock(
                    return_value=(2, 3, 4)
                ),
                build_band_path=Mock(
                    return_value=band_path
                ),
                run_ph=Mock(
                    return_value={'stage': 'ph'}
                ),
                run_q2r=Mock(
                    return_value={'stage': 'q2r'}
                ),
                run_matdyn_band=Mock(
                    return_value={
                        'stage': 'band',
                        'frequencies': {
                            'nqpoints': 2,
                        },
                    }
                ),
                run_matdyn_dos=Mock(
                    return_value={
                        'stage': 'dos',
                        'dos': {
                            'npoints': 2,
                        },
                    }
                ),
                write_matdyn_band_data=Mock(
                    side_effect=(
                        lambda output_file, **kwargs: output_file
                    )
                ),
                write_matdyn_dos_data=Mock(
                    side_effect=(
                        lambda output_file, **kwargs: output_file
                    )
                ),
                calculate_phonon_thermal_properties=Mock(
                    return_value=thermal_data
                ),
                write_phonon_thermal_properties=Mock(
                    side_effect=(
                        lambda output_file, **kwargs: output_file
                    )
                ),
            )
            solver._plot_qe_phonon_results = Mock(
                return_value=Path(
                    solver.struct
                    + '-PHONON-QE-Graph-Phonon.png'
                )
            )

            workflow = solver._phononcalc_qe()

        solver.engine.validate_qe_xc.assert_called_once_with(
            'PBE',
            pseudo_xc='pbe',
        )
        solver.engine.has_qe_state.assert_called_once_with(
            Path(
                solver.struct
                + '-GROUND-QE-Result-State'
            ),
            prefix='nanoworks',
        )
        solver.engine.run_ph.assert_called_once()
        solver.engine.run_q2r.assert_called_once()
        solver.engine.run_matdyn_band.assert_called_once()
        solver.engine.run_matdyn_dos.assert_called_once()
        solver.engine.write_matdyn_band_data.assert_called_once()
        solver.engine.write_matdyn_dos_data.assert_called_once()
        solver._plot_qe_phonon_results.assert_called_once()
        (
            solver.engine.calculate_phonon_thermal_properties
            .assert_called_once_with(
                {
                    'npoints': 2,
                },
                t_min=0.0,
                t_max=600.0,
                t_step=20.0,
            )
        )
        solver.engine.write_phonon_thermal_properties.assert_called_once()
        self.assertEqual(
            solver.engine.run_ph.call_args.kwargs[
                'qpoint_grid'
            ],
            (2, 3, 4),
        )
        self.assertEqual(
            solver.engine.run_matdyn_dos.call_args.kwargs[
                'qpoint_grid'
            ],
            (12, 13, 14),
        )
        self.assertIs(
            solver.engine.run_matdyn_band.call_args.kwargs[
                'band_path'
            ],
            band_path,
        )
        self.assertEqual(
            workflow['ph'],
            {'stage': 'ph'},
        )
        self.assertEqual(
            workflow['dos'],
            {
                'stage': 'dos',
                'dos': {
                    'npoints': 2,
                },
            },
        )
        self.assertEqual(
            workflow['band_data_file'],
            Path(
                solver.struct
                + '-PHONON-QE-Result-Band-THz.dat'
            ),
        )
        self.assertEqual(
            workflow['dos_data_file'],
            Path(
                solver.struct
                + '-PHONON-QE-Result-DOS-THz.dat'
            ),
        )
        self.assertEqual(
            workflow['graph_file'],
            Path(
                solver.struct
                + '-PHONON-QE-Graph-Phonon.png'
            ),
        )
        self.assertIs(
            workflow['thermal_data'],
            thermal_data,
        )
        self.assertEqual(
            workflow['thermal_data_file'],
            Path(
                solver.struct
                + '-PHONON-QE-Result-Thermal-Properties.csv'
            ),
        )

    def test_plot_qe_phonon_results(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            solver = object.__new__(
                DFTSolver
            )
            solver.engine = SimpleNamespace(
                THZ_PER_CM_MINUS_ONE=0.0299792458,
            )
            output_file = (
                Path(tmpdir)
                / 'phonon.png'
            )

            result = solver._plot_qe_phonon_results(
                output_file=output_file,
                band_path={
                    'distances': [0.0, 0.5, 1.0],
                    'special_distances': [0.0, 1.0],
                    'labels': ['G', 'X'],
                },
                frequencies={
                    'frequencies_thz': [
                        [-0.2, 1.0, 2.0],
                        [0.0, 1.5, 2.5],
                        [0.2, 2.0, 3.0],
                    ],
                },
                dos_data={
                    'frequencies_thz': [-0.2, 0.0, 3.0],
                    'dos': [0.0, 0.02, 0.0],
                },
            )

            self.assertEqual(
                result,
                output_file,
            )
            self.assertTrue(
                output_file.is_file()
            )
            self.assertGreater(
                output_file.stat().st_size,
                0,
            )

    def test_plot_qe_phonon_results_rejects_bad_band_shape(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.engine = SimpleNamespace(
            THZ_PER_CM_MINUS_ONE=0.0299792458,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(
                ValueError,
                'band data dimensions',
            ):
                solver._plot_qe_phonon_results(
                    output_file=(
                        Path(tmpdir)
                        / 'bad.png'
                    ),
                    band_path={
                        'distances': [0.0, 1.0],
                        'special_distances': [],
                        'labels': [],
                    },
                    frequencies={
                        'frequencies_thz': [[1.0, 2.0]],
                    },
                    dos_data={
                        'frequencies_thz': [0.0, 1.0],
                        'dos': [0.0, 1.0],
                    },
                )

    def test_qe_phononcalc_requires_ground_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            solver = object.__new__(
                DFTSolver
            )
            solver.struct = str(
                Path(tmpdir)
                / 'silicon'
            )
            solver.Mode = 'PW'
            solver.XC_calc = 'PBE'
            solver.Phonon_thermal_calc = False
            solver.engine = SimpleNamespace(
                validate_qe_xc=Mock(
                    return_value='pbe'
                ),
                has_qe_state=Mock(
                    return_value=False
                ),
            )

            with self.assertRaisesRegex(
                FileNotFoundError,
                'ground-state result',
            ):
                solver._phononcalc_qe()

    def test_qe_groundcalc_passes_hybrid_settings_to_scf(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Mode = 'PW'
        solver.Engine = 'QE'
        solver.struct = 'silicon'
        solver.XC_calc = 'HSE06'
        solver.XC_exx_fraction = 0.30
        solver.XC_omega = 0.12
        solver.Ground_calc = True
        solver.Geo_optim = False
        solver.Gamma = False
        solver.Ground_gamma = None
        solver.Spin_calc = False
        solver.bulk_configuration = Atoms(
            'Si',
            cell=[5.4, 5.4, 5.4],
            pbc=True,
        )
        solver.Cut_off_energy = 500.0
        solver.Ground_kpts_density = None
        solver.Ground_kpts_x = 2
        solver.Ground_kpts_y = 2
        solver.Ground_kpts_z = 2
        solver.Total_charge = 0.0
        solver.Ground_num_of_bands = None
        solver.Setup_params = None
        solver.Occupation = None
        solver.parallel_cores = 1
        solver.config = SimpleNamespace(
            vdW_calc='NONE',
        )
        solver.engine = SimpleNamespace(
            validate_qe_xc=Mock(
                return_value='hse06'
            ),
            run_scf=Mock(
                return_value={
                    'result': {
                        'total_energy_ev': -10.0,
                        'total_energy_ry': -0.75,
                    },
                }
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.get_qe_pseudo_dir',
                return_value=Path('pseudos'),
            ),
            patch(
                'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                return_value={
                    'Si': 'Si.upf',
                },
            ),
            patch(
                'nanoworks.dftsolve.write_cif',
            ),
            patch(
                'nanoworks.dftsolve.parprint',
            ),
        ):
            solver._groundcalc_qe()

        solver.engine.validate_qe_xc.assert_called_once_with(
            'HSE06',
            pseudo_xc='pbe',
            allow_hybrid=True,
        )
        call = solver.engine.run_scf.call_args.kwargs
        self.assertEqual(
            call['xc_calc'],
            'HSE06',
        )
        self.assertEqual(
            call['exx_fraction'],
            0.30,
        )
        self.assertEqual(
            call['omega'],
            0.12,
        )

    def test_qe_densitycalc_accepts_hybrid_ground_state(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.struct = 'silicon'
        solver.XC_calc = 'HSE06'
        solver.Spin_calc = False
        solver.parallel_cores = 1
        solver.engine = SimpleNamespace(
            validate_qe_xc=Mock(
                return_value='hse06'
            ),
            run_pp_density=Mock(
                side_effect=RuntimeError(
                    'stop after density validation'
                )
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.parprint',
            ),
            self.assertRaisesRegex(
                RuntimeError,
                'stop after density validation',
            ),
        ):
            solver._densitycalc_qe()

        solver.engine.validate_qe_xc.assert_called_once_with(
            'HSE06',
            pseudo_xc='pbe',
            allow_hybrid=True,
        )
        solver.engine.run_pp_density.assert_called_once()

    def test_qe_spin_densitycalc_dispatches_all_components(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.struct = 'iron'
        solver.XC_calc = 'PBE'
        solver.Spin_calc = True
        solver.parallel_cores = 4
        solver.engine = SimpleNamespace(
            validate_qe_xc=Mock(
                return_value='pbe'
            ),
            run_pp_density=Mock(
                side_effect=lambda **kwargs: kwargs
            ),
        )

        with patch(
            'nanoworks.dftsolve.parprint',
        ):
            outputs = solver._densitycalc_qe()

        solver.engine.validate_qe_xc.assert_called_once_with(
            'PBE',
            pseudo_xc='pbe',
            allow_hybrid=True,
        )
        self.assertEqual(
            set(outputs),
            {
                'Pseudo-Total',
                'Pseudo-Up',
                'Pseudo-Down',
                'Spin-Density',
            },
        )

        calls = [
            (
                call.kwargs['cube_file'].name,
                call.kwargs['plot_num'],
                call.kwargs['spin_component'],
            )
            for call in (
                solver.engine.run_pp_density.call_args_list
            )
        ]
        self.assertEqual(
            calls,
            [
                (
                    'iron-EDENSITY-QE-Result-Pseudo-Total.cube',
                    0,
                    0,
                ),
                (
                    'iron-EDENSITY-QE-Result-Pseudo-Up.cube',
                    0,
                    1,
                ),
                (
                    'iron-EDENSITY-QE-Result-Pseudo-Down.cube',
                    0,
                    2,
                ),
                (
                    'iron-EDENSITY-QE-Result-Spin-Density.cube',
                    6,
                    None,
                ),
            ],
        )

        for call in solver.engine.run_pp_density.call_args_list:
            self.assertEqual(
                call.kwargs['state_dir'],
                Path('iron-GROUND-QE-Result-State'),
            )
            self.assertEqual(
                call.kwargs['parallel_cores'],
                4,
            )

    def test_qe_doscalc_dispatches_hybrid_dos_workflow(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Mode = 'PW'
        solver.SOC_calc = False
        solver.Engine = 'QE'
        solver.struct = 'silicon'
        solver.XC_calc = 'HSE03'
        solver.XC_exx_fraction = 0.28
        solver.XC_omega = 0.15
        solver.Gamma = False
        solver.Ground_gamma = None
        solver.Ground_kpts_density = None
        solver.Ground_kpts_x = 2
        solver.Ground_kpts_y = 2
        solver.Ground_kpts_z = 2
        solver.DOS_kpts_density = None
        solver.DOS_kpts_x = 4
        solver.DOS_kpts_y = 4
        solver.DOS_kpts_z = 4
        solver.DOS_gamma = None
        solver.DOS_occupation = 'tetrahedra'
        solver.Occupation = None
        solver.Spin_calc = False
        solver.bulk_configuration = Atoms(
            'Si',
            cell=[5.4, 5.4, 5.4],
            pbc=True,
        )
        solver.Cut_off_energy = 500.0
        solver.Total_charge = 0.0
        solver.DOS_num_of_bands = 16
        solver.DOS_npoints = 161
        solver.Energy_min = -8.0
        solver.Energy_max = 8.0
        solver.Setup_params = None
        solver.parallel_cores = 1
        solver.engine = SimpleNamespace(
            validate_qe_xc=Mock(
                return_value='hse03'
            ),
            has_qe_state=Mock(
                return_value=True
            ),
            resolve_qe_occupation=Mock(
                return_value={
                    'occupations': 'tetrahedra',
                }
            ),
            run_nscf=Mock(
                side_effect=RuntimeError(
                    'stop after NSCF'
                )
            ),
            run_hybrid_dos=Mock(
                side_effect=RuntimeError(
                    'stop after hybrid DOS'
                )
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.get_qe_pseudo_dir',
                return_value=Path('pseudos'),
            ),
            patch(
                'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                return_value={
                    'Si': 'Si.upf',
                },
            ),
            patch(
                'nanoworks.dftsolve.parprint',
            ),
            self.assertRaisesRegex(
                RuntimeError,
                'stop after hybrid DOS',
            ),
        ):
            solver._doscalc_qe()

        solver.engine.validate_qe_xc.assert_called_once_with(
            'HSE03',
            pseudo_xc='pbe',
            allow_hybrid=True,
        )
        solver.engine.run_nscf.assert_not_called()
        solver.engine.run_hybrid_dos.assert_called_once()
        hybrid_call = solver.engine.run_hybrid_dos.call_args.kwargs
        self.assertEqual(
            hybrid_call['state_dir'],
            Path('silicon-DOS-QE-Result-State'),
        )
        self.assertEqual(
            hybrid_call['emin'],
            -8.0,
        )
        self.assertEqual(
            hybrid_call['emax'],
            8.0,
        )
        self.assertTrue(
            hybrid_call['relative_to_fermi']
        )

    def test_qe_bandcalc_dispatches_hybrid_projected_bands(self):
        solver = object.__new__(
            DFTSolver
        )
        solver.Mode = 'PW'
        solver.SOC_calc = False
        solver.struct = 'silicon'
        solver.XC_calc = 'PBE0'
        solver.XC_exx_fraction = 0.32
        solver.XC_omega = None
        solver.Projected_band_plot = True
        solver.Projections = []
        solver.Gamma = False
        solver.Ground_gamma = None
        solver.Ground_kpts_density = None
        solver.Ground_kpts_x = 2
        solver.Ground_kpts_y = 2
        solver.Ground_kpts_z = 2
        solver.Spin_calc = False
        solver.bulk_configuration = Atoms(
            'Si',
            cell=[5.4, 5.4, 5.4],
            pbc=True,
        )
        solver.Band_path = 'GX'
        solver.Band_npoints = 5
        solver.Cut_off_energy = 500.0
        solver.Total_charge = 0.0
        solver.Band_num_of_bands = 16
        solver.Setup_params = None
        solver.Occupation = None
        solver.parallel_cores = 1
        band_path = {
            'option': 'crystal',
            'kpoints': [
                (0.0, 0.0, 0.0),
                (0.5, 0.0, 0.0),
            ],
            'npoints': 2,
        }
        solver.engine = SimpleNamespace(
            validate_qe_xc=Mock(
                return_value='pbe0'
            ),
            has_qe_state=Mock(
                return_value=True
            ),
            build_band_path=Mock(
                return_value=band_path
            ),
            resolve_qe_kpoint_size=Mock(
                return_value=(2, 2, 2)
            ),
            run_hybrid_bands=Mock(
                side_effect=RuntimeError(
                    'stop after hybrid bands'
                )
            ),
        )

        with (
            patch(
                'nanoworks.dftsolve.get_qe_pseudo_dir',
                return_value=Path('pseudos'),
            ),
            patch(
                'nanoworks.dftsolve.resolve_qe_pseudopotentials',
                return_value={
                    'Si': 'Si.upf',
                },
            ),
            patch(
                'nanoworks.dftsolve.parprint',
            ),
            self.assertRaisesRegex(
                RuntimeError,
                'stop after hybrid bands',
            ),
        ):
            solver._bandcalc_qe()

        solver.engine.validate_qe_xc.assert_called_once_with(
            'PBE0',
            pseudo_xc='pbe',
            allow_hybrid=True,
        )
        solver.engine.resolve_qe_kpoint_size.assert_called_once_with(
            solver.bulk_configuration,
            density=None,
            size=(2, 2, 2),
        )
        solver.engine.run_hybrid_bands.assert_called_once()

        hybrid_call = (
            solver.engine.run_hybrid_bands.call_args.kwargs
        )
        self.assertEqual(
            hybrid_call['state_dir'],
            Path('silicon-BAND-QE-Result-State'),
        )
        self.assertEqual(
            hybrid_call['band_path'],
            band_path,
        )
        self.assertEqual(
            hybrid_call['qpoint_grid'],
            (2, 2, 2),
        )
        self.assertFalse(
            hybrid_call['gamma']
        )
        self.assertTrue(
            hybrid_call['projected_band']
        )
        self.assertEqual(
            hybrid_call['projections'],
            [],
        )
        self.assertEqual(
            hybrid_call['projection_prefix'],
            Path('silicon-BAND-QE-Result-Projections'),
        )

if __name__ == '__main__':
    unittest.main()
