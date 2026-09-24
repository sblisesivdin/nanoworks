import tempfile
import unittest
from pathlib import Path

from nanoworks.convergence import (
    ORDERED_CONVERGENCE_TASKS,
    StaticEnergyResult,
    build_convergence_plan,
    normalize_convergence_tasks,
    run_cutoff_sweep,
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


if __name__ == '__main__':
    unittest.main()
