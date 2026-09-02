import tempfile
import unittest
import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from ase import Atoms
from ase.io import write

with patch.object(
    sys,
    'argv',
    [sys.argv[0]],
):
    from nanoworks.dftsolve import (
        DFTConfig,
        dftsolve as DFTSolver,
    )


class TestDFTSolveWorkflow(unittest.TestCase):

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

if __name__ == '__main__':
    unittest.main()
