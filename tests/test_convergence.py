import tempfile
import unittest
from pathlib import Path

from ase import Atoms

from nanoworks.convergence import (
    ORDERED_CONVERGENCE_TASKS,
    StaticEnergyResult,
    build_convergence_plan,
    normalize_convergence_tasks,
    run_cutoff_sweep,
    run_kpoint_sweep,
    run_lattice_sweep,
    select_converged_value,
)


class TestConvergenceCore(unittest.TestCase):

    def test_default_tasks_use_dependency_safe_order(self):
        self.assertEqual(
            normalize_convergence_tasks(),
            ORDERED_CONVERGENCE_TASKS,
        )

    def test_requested_tasks_are_reordered_and_deduplicated(self):
        self.assertEqual(
            normalize_convergence_tasks(
                ['lattice', 'cutoff', 'lattice']
            ),
            ('cutoff', 'lattice'),
        )

    def test_unknown_task_is_rejected(self):
        with self.assertRaisesRegex(
            ValueError,
            'Unsupported convergence task',
        ):
            normalize_convergence_tasks(['cutoff', 'phonon'])

    def test_plan_requires_an_explicit_supported_engine(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'input.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text('', encoding='utf-8')
            geometry_file.write_text('', encoding='utf-8')

            with self.assertRaisesRegex(ValueError, 'must define Engine'):
                build_convergence_plan(
                    config={},
                    input_file=input_file,
                    geometry_file=geometry_file,
                )

    def test_plan_normalizes_engine_and_tasks(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            input_file = root / 'input.py'
            geometry_file = root / 'structure.cif'
            input_file.write_text('', encoding='utf-8')
            geometry_file.write_text('', encoding='utf-8')

            plan = build_convergence_plan(
                config={
                    'Engine': ' qe ',
                    'Convergence_tasks': 'lattice, cutoff, kpoints',
                },
                input_file=input_file,
                geometry_file=geometry_file,
                parallel_cores=8,
            )

            self.assertEqual(plan.engine, 'QE')
            self.assertEqual(
                plan.tasks,
                ('cutoff', 'kpoints', 'lattice'),
            )
            self.assertEqual(plan.parallel_cores, 8)

    def test_first_stable_convergence_window_is_selected(self):
        selection = select_converged_value(
            values=[400, 450, 500, 550],
            total_energies_ev=[
                -100.0000,
                -100.0100,
                -100.0107,
                -100.0110,
            ],
            atom_count=1,
            tolerance_ev_per_atom=0.001,
            consecutive_points=2,
        )

        self.assertIsNotNone(selection)
        self.assertEqual(selection.value, 500)
        self.assertEqual(selection.index, 2)
        self.assertAlmostEqual(
            selection.delta_ev_per_atom,
            0.0007,
        )

    def test_energy_changes_are_normalized_per_atom(self):
        selection = select_converged_value(
            values=[400, 450, 500],
            total_energies_ev=[-20.0, -20.0015, -20.0020],
            atom_count=2,
            tolerance_ev_per_atom=0.001,
            consecutive_points=2,
        )

        self.assertIsNotNone(selection)
        self.assertEqual(selection.value, 450)

    def test_unconverged_series_returns_none(self):
        selection = select_converged_value(
            values=[300, 400, 500],
            total_energies_ev=[-10.0, -10.1, -10.2],
            atom_count=1,
            tolerance_ev_per_atom=0.001,
            consecutive_points=2,
        )

        self.assertIsNone(selection)

    def test_series_must_cover_requested_window(self):
        with self.assertRaisesRegex(ValueError, 'At least 3 values'):
            select_converged_value(
                values=[400, 450],
                total_energies_ev=[-10.0, -10.001],
                atom_count=1,
                consecutive_points=2,
            )

    def test_cutoff_sweep_uses_injected_backend(self):
        class FakeBackend:
            name = 'FAKE'

            def __init__(self):
                self.calls = []

            def calculate_static_energy(
                self,
                atoms,
                *,
                cutoff_ev,
                kpoint_settings,
                workdir,
                settings,
                parallel_cores,
            ):
                self.calls.append({
                    'cutoff_ev': cutoff_ev,
                    'kpoint_settings': dict(kpoint_settings),
                    'workdir': workdir,
                    'settings': dict(settings),
                    'parallel_cores': parallel_cores,
                })
                energies = {
                    400.0: -20.0000,
                    450.0: -20.0200,
                    500.0: -20.0214,
                    550.0: -20.0220,
                }
                return StaticEnergyResult(
                    engine=self.name,
                    total_energy_ev=energies[cutoff_ev],
                    metadata={'mock': True},
                )

        backend = FakeBackend()

        with tempfile.TemporaryDirectory() as temp_dir:
            result = run_cutoff_sweep(
                backend=backend,
                atoms=[object(), object()],
                cutoff_values=[400, 450, 500, 550],
                kpoint_settings={
                    'density': 3.0,
                    'gamma': True,
                },
                workdir=Path(temp_dir),
                settings={'xc_calc': 'PBE'},
                parallel_cores=4,
                tolerance_ev_per_atom=0.001,
                consecutive_points=2,
            )

        self.assertEqual(len(backend.calls), 4)
        self.assertEqual(backend.calls[0]['cutoff_ev'], 400.0)
        self.assertEqual(backend.calls[0]['parallel_cores'], 4)
        self.assertEqual(
            backend.calls[0]['kpoint_settings']['density'],
            3.0,
        )
        self.assertEqual(result.selection.value, 500.0)
        self.assertAlmostEqual(
            result.points[-1].energy_ev_per_atom,
            -10.011,
        )

    def test_cutoff_sweep_rejects_unsorted_values_before_running(self):
        class UnusedBackend:
            def calculate_static_energy(self, *args, **kwargs):
                raise AssertionError('backend must not be called')

        with self.assertRaisesRegex(ValueError, 'strictly increasing'):
            run_cutoff_sweep(
                backend=UnusedBackend(),
                atoms=[object()],
                cutoff_values=[400, 500, 450],
                kpoint_settings={'size': (4, 4, 4)},
                workdir=Path('unused'),
            )

    def test_kpoint_density_sweep_uses_common_backend(self):
        class FakeBackend:
            name = 'FAKE'

            def __init__(self):
                self.calls = []

            def calculate_static_energy(self, atoms, **kwargs):
                self.calls.append(kwargs)
                energies = {
                    2.0: -20.0,
                    3.0: -20.02,
                    4.0: -20.021,
                    5.0: -20.0215,
                }
                density = kwargs['kpoint_settings']['density']
                return StaticEnergyResult(
                    engine=self.name,
                    total_energy_ev=energies[density],
                )

        backend = FakeBackend()
        with tempfile.TemporaryDirectory() as temp_dir:
            result = run_kpoint_sweep(
                backend=backend,
                atoms=[object(), object()],
                kpoint_values=[2.0, 3.0, 4.0, 5.0],
                cutoff_ev=500,
                workdir=Path(temp_dir),
                parallel_cores=4,
                gamma=True,
            )

        self.assertEqual(len(backend.calls), 4)
        self.assertEqual(backend.calls[0]['cutoff_ev'], 500.0)
        self.assertTrue(
            backend.calls[0]['kpoint_settings']['gamma']
        )
        self.assertNotIn(
            'size',
            backend.calls[0]['kpoint_settings'],
        )
        self.assertEqual(result.selection.value, 4.0)

    def test_kpoint_mesh_sweep_preserves_meshes(self):
        class FakeBackend:
            name = 'FAKE'

            def calculate_static_energy(self, atoms, **kwargs):
                mesh = kwargs['kpoint_settings']['size']
                energies = {
                    (2, 2, 1): -10.0,
                    (4, 4, 1): -10.001,
                    (6, 6, 1): -10.0015,
                }
                return StaticEnergyResult(
                    engine=self.name,
                    total_energy_ev=energies[mesh],
                )

        with tempfile.TemporaryDirectory() as temp_dir:
            result = run_kpoint_sweep(
                backend=FakeBackend(),
                atoms=[object()],
                kpoint_values=[(2, 2, 1), (4, 4, 1), (6, 6, 1)],
                cutoff_ev=450,
                workdir=Path(temp_dir),
            )

        self.assertEqual(result.selection.value, (4, 4, 1))
        self.assertEqual(
            result.points[-1].kpoint_settings['size'],
            (6, 6, 1),
        )
        self.assertNotIn(
            'density',
            result.points[-1].kpoint_settings,
        )

    def test_kpoint_sweep_rejects_unsorted_meshes(self):
        class UnusedBackend:
            def calculate_static_energy(self, *args, **kwargs):
                raise AssertionError('backend must not be called')

        with self.assertRaisesRegex(ValueError, 'strictly increasing'):
            run_kpoint_sweep(
                backend=UnusedBackend(),
                atoms=[object()],
                kpoint_values=[(4, 4, 1), (2, 2, 1), (6, 6, 1)],
                cutoff_ev=450,
                workdir=Path('unused'),
            )

    def test_lattice_sweep_scales_only_selected_periodic_axes(self):
        class FakeBackend:
            name = 'FAKE'

            def __init__(self):
                self.cells = []

            def calculate_static_energy(self, atoms, **kwargs):
                self.cells.append(atoms.cell.copy())
                scale = round(atoms.cell.lengths()[0] / 4.0, 2)
                energies = {
                    0.98: -10.0,
                    1.00: -10.1,
                    1.02: -10.05,
                }
                return StaticEnergyResult(
                    engine=self.name,
                    total_energy_ev=energies[scale],
                )

        atoms = Atoms(
            'Si2',
            positions=[(0, 0, 0), (1, 1, 0)],
            cell=(4, 4, 10),
            pbc=(True, True, False),
        )
        backend = FakeBackend()

        with tempfile.TemporaryDirectory() as temp_dir:
            result = run_lattice_sweep(
                backend=backend,
                atoms=atoms,
                lattice_scales=[0.98, 1.0, 1.02],
                cutoff_ev=500,
                kpoint_settings={'size': (6, 6, 1)},
                workdir=Path(temp_dir),
            )

        self.assertEqual(result.selection.scale, 1.0)
        self.assertAlmostEqual(backend.cells[0].lengths()[0], 3.92)
        self.assertAlmostEqual(backend.cells[0].lengths()[1], 3.92)
        self.assertAlmostEqual(backend.cells[0].lengths()[2], 10.0)
        self.assertAlmostEqual(
            result.optimized_atoms.cell.lengths()[0],
            4.0,
        )
        self.assertAlmostEqual(atoms.cell.lengths()[0], 4.0)

    def test_lattice_sweep_rejects_unbracketed_minimum(self):
        class FakeBackend:
            name = 'FAKE'

            def calculate_static_energy(self, atoms, **kwargs):
                scale = round(atoms.cell.lengths()[0] / 4.0, 2)
                return StaticEnergyResult(
                    engine=self.name,
                    total_energy_ev=-scale,
                )

        atoms = Atoms(
            'Si',
            cell=(4, 4, 4),
            pbc=True,
        )
        with tempfile.TemporaryDirectory() as temp_dir:
            result = run_lattice_sweep(
                backend=FakeBackend(),
                atoms=atoms,
                lattice_scales=[0.98, 1.0, 1.02],
                cutoff_ev=450,
                kpoint_settings={'density': 3.0},
                workdir=Path(temp_dir),
            )

        self.assertIsNone(result.selection)
        self.assertIsNone(result.optimized_atoms)


if __name__ == '__main__':
    unittest.main()
