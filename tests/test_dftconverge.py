import csv
import io
import json
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from unittest.mock import Mock, patch

from ase import Atoms

from nanoworks import dftconverge
from nanoworks.convergence import ConvergenceRunResult, StaticEnergyResult


class TestDFTConvergeCLI(unittest.TestCase):

    def test_check_prints_ordered_plan(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Convergence_tasks = ['lattice', 'cutoff', 'kpoints']\n"
                "Convergence_cutoffs = [300, 400, 500]\n"
                "Convergence_kpoints = [2.0, 3.0, 4.0]\n"
                "Convergence_lattice_scales = [0.98, 1.0, 1.02]\n",
                encoding='utf-8',
            )
            geometry_file.write_text('test geometry', encoding='utf-8')

            output = io.StringIO()
            with redirect_stdout(output):
                result = dftconverge.main([
                    '--check',
                    '-p',
                    '4',
                    '-i',
                    str(input_file),
                    '-g',
                    str(geometry_file),
                ])

        rendered = output.getvalue()
        self.assertEqual(result, 0)
        self.assertIn('Engine: QE', rendered)
        self.assertLess(
            rendered.index('1. cutoff'),
            rendered.index('2. kpoints'),
        )
        self.assertLess(
            rendered.index('2. kpoints'),
            rendered.index('3. lattice'),
        )
        self.assertIn('Parallel processes: 4', rendered)
        self.assertIn('no calculations executed', rendered)

    def test_check_rejects_missing_task_values_without_engines(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text(
                "Engine = 'GPAW'\n"
                "Convergence_tasks = ['cutoff']\n",
                encoding='utf-8',
            )
            geometry_file.write_text('not parsed by check', encoding='utf-8')

            errors = io.StringIO()
            with redirect_stderr(errors):
                with self.assertRaises(SystemExit) as raised:
                    dftconverge.main([
                        '--check',
                        '-i',
                        str(input_file),
                        '-g',
                        str(geometry_file),
                    ])

        self.assertEqual(raised.exception.code, 2)
        self.assertIn('Convergence_cutoffs is required', errors.getvalue())

    def test_check_rejects_mixed_kpoint_candidate_types(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Convergence_tasks = ['kpoints']\n"
                "Cut_off_energy = 500\n"
                "Convergence_kpoints = [2.0, (4, 4, 4), 5.0]\n",
                encoding='utf-8',
            )
            geometry_file.write_text('not parsed by check', encoding='utf-8')

            errors = io.StringIO()
            with redirect_stderr(errors):
                with self.assertRaises(SystemExit):
                    dftconverge.main([
                        '--check',
                        '-i',
                        str(input_file),
                        '-g',
                        str(geometry_file),
                    ])

        self.assertIn('cannot mix densities and meshes', errors.getvalue())

    def test_packaged_convergence_examples_pass_check(self):
        example_root = (
            Path(__file__).resolve().parents[1]
            / 'nanoworks'
            / 'examples'
            / 'Convergence'
        )

        for input_name, engine in (
            ('Si-GPAW-convergence.py', 'GPAW'),
            ('Si-QE-convergence.py', 'QE'),
        ):
            with self.subTest(engine=engine):
                output = io.StringIO()
                with redirect_stdout(output):
                    result = dftconverge.main([
                        '--check',
                        '-p',
                        '4',
                        '-i',
                        str(example_root / input_name),
                        '-g',
                        str(example_root / 'Si.cif'),
                    ])

                self.assertEqual(result, 0)
                self.assertIn('Engine: ' + engine, output.getvalue())
                self.assertIn('Result: VALID', output.getvalue())

    def test_gpaw_parallel_command_uses_gpaw_python(self):
        args = Mock(
            input='convergence.py',
            geometry='structure.cif',
        )

        def find_executable(name):
            return {
                'mpiexec': '/usr/bin/mpiexec',
                'gpaw': '/opt/gpaw/bin/gpaw',
            }.get(name)

        with patch.object(
            dftconverge.shutil,
            'which',
            side_effect=find_executable,
        ):
            command = dftconverge.build_gpaw_mpi_command(8, args)

        self.assertEqual(command[:7], [
            '/usr/bin/mpiexec',
            '-np',
            '8',
            '/opt/gpaw/bin/gpaw',
            'python',
            '--',
            str(Path(dftconverge.__file__).resolve()),
        ])
        self.assertEqual(
            command[-6:],
            ['-p', '8', '-i', 'convergence.py', '-g', 'structure.cif'],
        )

    def test_main_restarts_parallel_gpaw_before_execution(self):
        class Restarted(Exception):
            pass

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text(
                "Engine = 'GPAW'\n"
                "Convergence_tasks = ['cutoff']\n"
                "Convergence_cutoffs = [300, 400, 500]\n",
                encoding='utf-8',
            )
            geometry_file.write_text(
                'not read before restart',
                encoding='utf-8',
            )

            with patch.dict(
                dftconverge.os.environ,
                {dftconverge.GPAW_MPI_ENV: '0'},
            ), patch.object(
                dftconverge,
                'restart_gpaw_with_mpi',
                side_effect=Restarted,
            ) as restart, self.assertRaises(Restarted):
                dftconverge.main([
                    '-p',
                    '4',
                    '-i',
                    str(input_file),
                    '-g',
                    str(geometry_file),
                ])

        restart.assert_called_once()

    def test_cutoff_execution_maps_qe_configuration(self):
        class FakeBackend:
            name = 'QE'

            def __init__(self):
                self.calls = []

            def calculate_static_energy(self, atoms, **kwargs):
                call_index = len(self.calls)
                self.calls.append(kwargs)
                if call_index < 4:
                    energies = {
                        400.0: -20.0,
                        450.0: -20.02,
                        500.0: -20.021,
                        550.0: -20.0215,
                    }
                    energy = energies[kwargs['cutoff_ev']]
                elif call_index < 8:
                    density = kwargs['kpoint_settings']['density']
                    energies = {
                        2.0: -20.0,
                        3.0: -20.02,
                        4.0: -20.021,
                        5.0: -20.0215,
                    }
                    energy = energies[density]
                else:
                    scale = round(atoms.cell.lengths()[0] / 5.0, 2)
                    energies = {
                        0.98: -20.0,
                        1.00: -20.1,
                        1.02: -20.05,
                    }
                    energy = energies[scale]
                return StaticEnergyResult(
                    engine='QE',
                    total_energy_ev=energy,
                )

        backend = FakeBackend()
        backend_loader = Mock(return_value=backend)
        pseudo_dir_getter = Mock(return_value=Path('/pseudos'))
        pseudo_resolver = Mock(return_value={'Si': 'Si.upf'})

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text('', encoding='utf-8')
            geometry_file.write_text('', encoding='utf-8')
            plan = dftconverge.build_convergence_plan(
                config={
                    'Engine': 'QE',
                    'Convergence_tasks': [
                        'cutoff',
                        'kpoints',
                        'lattice',
                    ],
                },
                input_file=input_file,
                geometry_file=geometry_file,
                parallel_cores=6,
            )

            result = dftconverge.execute_convergence_plan(
                config={
                    'Engine': 'QE',
                    'Convergence_tasks': [
                        'cutoff',
                        'kpoints',
                        'lattice',
                    ],
                    'Convergence_cutoffs': [400, 450, 500, 550],
                    'Convergence_kpoints': [2.0, 3.0, 4.0, 5.0],
                    'Convergence_lattice_scales': [0.98, 1.0, 1.02],
                    'Ground_kpts_x': 6,
                    'Ground_kpts_y': 6,
                    'Ground_kpts_z': 2,
                    'XC_calc': 'PBE',
                },
                plan=plan,
                structure_reader=Mock(return_value=Atoms(
                    'Si2',
                    positions=[(0, 0, 0), (1, 1, 0)],
                    cell=(5, 5, 12),
                    pbc=(True, True, False),
                )),
                backend_loader=backend_loader,
                pseudo_dir_getter=pseudo_dir_getter,
                pseudo_resolver=pseudo_resolver,
            )
            artifacts = dftconverge.write_convergence_results(
                config={},
                plan=plan,
                result=result,
            )
            summary = json.loads(
                artifacts['summary'].read_text(encoding='utf-8')
            )
            with artifacts['table'].open(
                encoding='utf-8',
                newline='',
            ) as fd:
                rows = list(csv.DictReader(fd))
            optimized_structure_exists = (
                artifacts['optimized_structure'].is_file()
            )

        backend_loader.assert_called_once_with(
            'QE',
            pseudopotentials={'Si': 'Si.upf'},
            pseudo_dir=Path('/pseudos'),
            executable='pw.x',
        )
        self.assertEqual(len(backend.calls), 11)
        self.assertEqual(
            backend.calls[0]['kpoint_settings']['size'],
            (6, 6, 2),
        )
        self.assertEqual(backend.calls[0]['parallel_cores'], 6)
        self.assertEqual(result.cutoff.selection.value, 500.0)
        self.assertEqual(result.kpoints.selection.value, 4.0)
        self.assertEqual(result.lattice.selection.scale, 1.0)
        self.assertTrue(all(
            call['cutoff_ev'] == 500.0
            for call in backend.calls[4:]
        ))
        self.assertTrue(all(
            call['kpoint_settings']['density'] == 4.0
            for call in backend.calls[8:]
        ))
        self.assertEqual(summary['schema_version'], 1)
        self.assertEqual(summary['selected']['cutoff_ev'], 500.0)
        self.assertEqual(
            summary['selected']['kpoints'],
            {'density': 4.0},
        )
        self.assertEqual(summary['selected']['lattice_scale'], 1.0)
        self.assertEqual(len(rows), 11)
        self.assertEqual(rows[-1]['task'], 'lattice')
        self.assertTrue(optimized_structure_exists)

    def test_cutoff_execution_maps_gpaw_spin_configuration(self):
        atoms = Mock()
        atoms.__len__ = Mock(return_value=2)
        atoms.get_chemical_symbols.return_value = ['Fe', 'Fe']
        backend = Mock()
        backend.name = 'GPAW'
        backend.calculate_static_energy.side_effect = [
            StaticEnergyResult('GPAW', -20.0),
            StaticEnergyResult('GPAW', -20.001),
            StaticEnergyResult('GPAW', -20.0015),
        ]
        backend_loader = Mock(return_value=backend)

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text('', encoding='utf-8')
            geometry_file.write_text('', encoding='utf-8')
            plan = dftconverge.build_convergence_plan(
                config={
                    'Engine': 'GPAW',
                    'Convergence_tasks': ['cutoff'],
                },
                input_file=input_file,
                geometry_file=geometry_file,
            )

            dftconverge.execute_convergence_plan(
                config={
                    'Convergence_cutoffs': [300, 350, 400],
                    'Spin_calc': True,
                    'Magmom_per_atom': {'Fe': 2.5},
                },
                plan=plan,
                structure_reader=Mock(return_value=atoms),
                backend_loader=backend_loader,
            )

        backend_loader.assert_called_once_with('GPAW')
        settings = (
            backend.calculate_static_energy.call_args_list[0]
            .kwargs['settings']
        )
        self.assertEqual(settings['magnetic_moments'], [2.5, 2.5])
        self.assertEqual(settings['xc_calc'], 'LDA')

    def test_lattice_only_requires_a_fixed_cutoff(self):
        backend = Mock()
        backend_loader = Mock(return_value=backend)

        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text('', encoding='utf-8')
            geometry_file.write_text('', encoding='utf-8')
            plan = dftconverge.build_convergence_plan(
                config={
                    'Engine': 'GPAW',
                    'Convergence_tasks': ['lattice'],
                },
                input_file=input_file,
                geometry_file=geometry_file,
            )

            with self.assertRaisesRegex(
                ValueError,
                'requires Cut_off_energy',
            ):
                dftconverge.execute_convergence_plan(
                    config={
                        'Convergence_lattice_scales': [0.98, 1.0, 1.02],
                    },
                    plan=plan,
                    structure_reader=Mock(return_value=Atoms(
                        'Si',
                        cell=(5, 5, 5),
                        pbc=True,
                    )),
                    backend_loader=backend_loader,
                )

        backend.calculate_static_energy.assert_not_called()

    def test_main_executes_cutoff_workflow(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'convergence.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text(
                "Engine = 'QE'\n"
                "Convergence_tasks = ['cutoff']\n"
                "Convergence_cutoffs = [400, 450, 500]\n",
                encoding='utf-8',
            )
            geometry_file.write_text('', encoding='utf-8')
            fake_result = Mock()
            fake_result.points = (
                Mock(
                    cutoff_ev=400.0,
                    total_energy_ev=-10.0,
                    energy_ev_per_atom=-5.0,
                ),
            )
            fake_result.selection = Mock(value=400.0)
            fake_run_result = ConvergenceRunResult(cutoff=fake_result)

            output = io.StringIO()
            with patch.object(
                dftconverge,
                'execute_convergence_plan',
                return_value=fake_run_result,
            ) as execute, redirect_stdout(output):
                result = dftconverge.main([
                    '-i',
                    str(input_file),
                    '-g',
                    str(geometry_file),
                ])

        self.assertEqual(result, 0)
        execute.assert_called_once()
        self.assertIn('Selected cutoff: 400 eV', output.getvalue())


if __name__ == '__main__':
    unittest.main()
