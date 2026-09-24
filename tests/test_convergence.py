import tempfile
import unittest
from pathlib import Path

from nanoworks.convergence import (
    ORDERED_CONVERGENCE_TASKS,
    build_convergence_plan,
    normalize_convergence_tasks,
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


if __name__ == '__main__':
    unittest.main()
