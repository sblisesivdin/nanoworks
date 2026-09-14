import tempfile
import unittest
import sys
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

if __name__ == '__main__':
    unittest.main()
