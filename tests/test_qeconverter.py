import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from ase.units import Bohr
from nanoworks.qeconverter import (
    QEInputSettings,
    build_config_lines,
    determine_system_name,
    parse_qe_input,
    _parse_qe_bool,
    RY_TO_EV,
    _build_workflow_lines,
    _build_kpoint_lines,
)


class TestQEConverter(unittest.TestCase):

    def test_parse_qe_input_reads_basic_scf_settings(self):
        input_text = """
&CONTROL
  calculation = 'scf',
/
&SYSTEM
  ibrav = 0,
  nat = 2,
  ntyp = 1,
  ecutwfc = 30.0,
  tot_charge = -1.0,
  nbnd = 24,
  occupations = 'smearing',
  smearing = 'mp',
  degauss = 0.01,
  nspin = 2,
  starting_magnetization(1) = 0.5,
/
&ELECTRONS
  conv_thr = 1.0d-8,
/

K_POINTS automatic
4 5 6 0 0 0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = (
                Path(tmpdir)
                / 'silicon.scf.in'
            )

            input_file.write_text(
                input_text,
                encoding='utf-8',
            )

            settings = parse_qe_input(
                input_file
            )

        self.assertEqual(
            settings.calculation,
            'scf',
        )
        self.assertAlmostEqual(
            settings.ecutwfc,
            30.0,
        )
        self.assertEqual(
            settings.occupations,
            'smearing',
        )
        self.assertEqual(
            settings.smearing,
            'mp',
        )
        self.assertAlmostEqual(
            settings.degauss,
            0.01,
        )
        self.assertEqual(
            settings.nspin,
            2,
        )
        self.assertAlmostEqual(
            settings.starting_magnetization['1'],
            0.5,
        )
        self.assertAlmostEqual(
            settings.conv_thr,
            1.0e-8,
        )
        self.assertEqual(
            settings.k_mesh,
            [4, 5, 6],
        )
        self.assertEqual(
            settings.k_shift,
            [0, 0, 0],
        )
        self.assertAlmostEqual(
            settings.total_charge,
            -1.0,
        )
        self.assertEqual(
            settings.nbands,
            24,
        )

    def test_build_config_lines_emits_basic_settings(self):
        settings = QEInputSettings(
            calculation='scf',
            ecutwfc=30.0,
            occupations='smearing',
            smearing='mp',
            degauss=0.01,
            nspin=2,
            total_charge=-1.0,
            nbands=24,
            k_mesh=[
                4,
                5,
                6,
            ],
            k_shift=[
                0,
                0,
                0,
            ],
        )

        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        lines = build_config_lines(
            name='Silicon',
            geom_filename='Silicon.cif',
            settings=settings,
            args=args,
        )

        text = '\n'.join(
            lines
        )

        self.assertIn(
            "Outdirname = 'Silicon-results'",
            text,
        )
        self.assertIn(
            "Mode = 'PW'",
            text,
        )
        self.assertIn(
            "Ground_calc = True",
            text,
        )
        self.assertIn(
            "Geo_optim = False",
            text,
        )
        self.assertIn(
            "Cut_off_energy = 408.2",
            text,
        )
        self.assertIn(
            "Ground_kpts_x = 4",
            text,
        )
        self.assertIn(
            "Ground_kpts_y = 5",
            text,
        )
        self.assertIn(
            "Ground_kpts_z = 6",
            text,
        )
        self.assertIn(
            "Ground_gamma = True",
            text,
        )
        self.assertIn(
            "XC_calc = 'PBE'",
            text,
        )
        self.assertIn(
            (
                "Occupation = {"
                "'name': 'methfessel-paxton', "
                "'width': 0.13605693009}"
            ),
            text,
        )
        self.assertIn(
            "Spin_calc = True",
            text,
        )
        self.assertIn(
            (
                "# Geometry file to use with "
                "dftsolve: Silicon.cif"
            ),
            text,
        )
        
        self.assertIn(
            "Engine = 'QE'",
            text,
        )

        self.assertIn(
            "Total_charge = -1",
            text,
        )
        self.assertIn(
            "Ground_num_of_bands = 24",
            text,
        )

    def test_determine_system_name_sanitizes_name(self):
        name = determine_system_name(
            Path(
                'sample.scf.in'
            ),
            'Fe spin calculation',
        )

        self.assertEqual(
            name,
            'Fe_spin_calculation',
        )

    def test_determine_system_name_uses_input_stem(self):
        name = determine_system_name(
            Path(
                'gaas.relax.in'
            ),
            None,
        )

        self.assertEqual(
            name,
            'gaas_relax',
        )

    def test_parse_qe_input_handles_multiple_assignments(self):
        input_text = """
&CONTROL
  calculation = 'bands',
/
&SYSTEM
  ibrav = 0, nat = 2, ntyp = 1,
  ecutwfc = 4.0D1, occupations = 'fixed', nspin = 1,
/
K_POINTS {automatic}
6 6 4 1 1 0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = (
                Path(tmpdir)
                / 'bands.in'
            )

            input_file.write_text(
                input_text,
                encoding='utf-8',
            )

            settings = parse_qe_input(
                input_file
            )

        self.assertEqual(
            settings.calculation,
            'bands',
        )
        self.assertAlmostEqual(
            settings.ecutwfc,
            40.0,
        )
        self.assertEqual(
            settings.occupations,
            'fixed',
        )
        self.assertEqual(
            settings.nspin,
            1,
        )
        self.assertEqual(
            settings.k_mesh,
            [
                6,
                6,
                4,
            ],
        )
        self.assertEqual(
            settings.k_shift,
            [
                1,
                1,
                0,
            ],
        )

    def test_parse_qe_input_reads_gamma_card_without_data_row(self):
        input_text = """
&CONTROL
  calculation = 'scf',
/
&SYSTEM
  ecutwfc = 30.0,
/
K_POINTS gamma
CELL_PARAMETERS angstrom
5.0 0.0 0.0
0.0 5.0 0.0
0.0 0.0 5.0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = (
                Path(tmpdir)
                / 'gamma.in'
            )

            input_file.write_text(
                input_text,
                encoding='utf-8',
            )

            settings = parse_qe_input(
                input_file
            )

        self.assertEqual(
            settings.k_mesh,
            [
                1,
                1,
                1,
            ],
        )
        self.assertEqual(
            settings.k_shift,
            [
                0,
                0,
                0,
            ],
        )
        self.assertAlmostEqual(
            settings.ecutwfc,
            30.0,
        )

    def test_build_config_lines_preserves_fixed_occupation(self):
        settings = QEInputSettings(
            calculation='scf',
            occupations='fixed',
        )

        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        text = '\n'.join(
            build_config_lines(
                name='Silicon',
                geom_filename='Silicon.cif',
                settings=settings,
                args=args,
            )
        )

        self.assertIn(
            "Occupation = 'fixed'",
            text,
        )

    def test_build_config_lines_preserves_tetrahedra(self):
        settings = QEInputSettings(
            calculation='scf',
            occupations='tetrahedra_opt',
        )

        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        text = '\n'.join(
            build_config_lines(
                name='Silicon',
                geom_filename='Silicon.cif',
                settings=settings,
                args=args,
            )
        )

        self.assertIn(
            "Occupation = 'tetrahedra_opt'",
            text,
        )

    def test_parse_qe_input_reads_relaxation_settings(self):
        input_text = """
&CONTROL
  calculation = 'vc-relax',
  forc_conv_thr = 1.0D-3,
/
&SYSTEM
  nosym = .false.,
/
&IONS
  ion_dynamics = 'bfgs',
  trust_radius_max = 2.0D-1,
/
&CELL
  cell_dynamics = 'bfgs',
  cell_dofree = 'all',
  press = 5.0,
/
K_POINTS automatic
4 4 4 0 0 0
"""

        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = (
                Path(tmpdir)
                / 'vc-relax.in'
            )

            input_file.write_text(
                input_text,
                encoding='utf-8',
            )

            settings = parse_qe_input(
                input_file
            )

        self.assertEqual(
            settings.calculation,
            'vc-relax',
        )
        self.assertAlmostEqual(
            settings.forc_conv_thr,
            1.0e-3,
        )
        self.assertEqual(
            settings.ion_dynamics,
            'bfgs',
        )
        self.assertAlmostEqual(
            settings.trust_radius_max,
            0.2,
        )
        self.assertEqual(
            settings.cell_dynamics,
            'bfgs',
        )
        self.assertEqual(
            settings.cell_dofree,
            'all',
        )
        self.assertAlmostEqual(
            settings.pressure_kbar,
            5.0,
        )
        self.assertFalse(
            settings.nosym
        )

    def test_parse_qe_bool_supports_fortran_values(self):
        self.assertTrue(
            _parse_qe_bool(
                '.true.'
            )
        )
        self.assertFalse(
            _parse_qe_bool(
                '.false.'
            )
        )

        with self.assertRaisesRegex(
            ValueError,
            'Unable to parse QE logical value',
        ):
            _parse_qe_bool(
                'maybe'
            )

    def test_build_config_lines_converts_vc_relax_settings(self):
        settings = QEInputSettings(
            calculation='vc-relax',
            forc_conv_thr=1.0e-3,
            ion_dynamics='bfgs',
            trust_radius_max=0.2,
            cell_dynamics='bfgs',
            cell_dofree='all',
            pressure_kbar=5.0,
            nosym=False,
        )

        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        text = '\n'.join(
            build_config_lines(
                name='GaAs',
                geom_filename='GaAs.cif',
                settings=settings,
                args=args,
            )
        )

        self.assertIn(
            "Geo_optim = True",
            text,
        )
        self.assertIn(
            "Optimizer = 'LBFGS'",
            text,
        )
        self.assertIn(
            (
                "Relax_cell = [True, True, True, "
                "True, True, True]"
            ),
            text,
        )
        self.assertIn(
            "Fix_symmetry = True",
            text,
        )
        self.assertIn(
            "Hydrostatic_pressure = 0.5",
            text,
        )

        force_line = next(
            line
            for line in text.splitlines()
            if line.startswith(
                'Max_F_tolerance ='
            )
        )

        max_force = float(
            force_line.split(
                '=',
                1,
            )[1]
        )

        self.assertAlmostEqual(
            max_force,
            1.0e-3 * RY_TO_EV / Bohr,
        )

        step_line = next(
            line
            for line in text.splitlines()
            if line.startswith(
                'Max_step ='
            )
        )

        max_step = float(
            step_line.split(
                '=',
                1,
            )[1]
        )

        self.assertAlmostEqual(
            max_step,
            0.2 * Bohr,
        )

    def test_build_config_lines_marks_volume_approximation(self):
        settings = QEInputSettings(
            calculation='vc-relax',
            cell_dofree='volume',
            nosym=True,
        )

        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        text = '\n'.join(
            build_config_lines(
                name='GaAs',
                geom_filename='GaAs.cif',
                settings=settings,
                args=args,
            )
        )

        self.assertIn(
            (
                "# NOTICE: QE cell_dofree = 'volume' "
                "cannot be represented exactly"
            ),
            text,
        )
        self.assertIn(
            (
                "Relax_cell = [True, True, True, "
                "True, True, True]"
            ),
            text,
        )
        self.assertIn(
            "Fix_symmetry = True",
            text,
        )

    def test_build_config_lines_falls_back_for_unknown_cell_dofree(
        self,
    ):
        settings = QEInputSettings(
            calculation='vc-relax',
            cell_dofree='custom-mode',
        )

        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        text = '\n'.join(
            build_config_lines(
                name='GaAs',
                geom_filename='GaAs.cif',
                settings=settings,
                args=args,
            )
        )

        self.assertIn(
            (
                "# NOTICE: Unknown QE cell_dofree = "
                "'custom-mode'"
            ),
            text,
        )
        self.assertIn(
            (
                "Relax_cell = [True, True, True, "
                "True, True, True]"
            ),
            text,
        )

    def test_build_workflow_lines_maps_nscf_to_dos(self):
        settings = QEInputSettings(
            calculation='nscf',
        )

        text = '\n'.join(
            _build_workflow_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_calc = True",
            text,
        )
        self.assertIn(
            "DOS_calc = True",
            text,
        )
        self.assertIn(
            "Band_calc = False",
            text,
        )
        self.assertIn(
            (
                "# NOTICE: QE calculation = 'nscf' "
                "is interpreted as a DOS workflow."
            ),
            text,
        )

    def test_build_workflow_lines_maps_bands(self):
        settings = QEInputSettings(
            calculation='bands',
        )

        text = '\n'.join(
            _build_workflow_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_calc = True",
            text,
        )
        self.assertIn(
            "DOS_calc = False",
            text,
        )
        self.assertIn(
            "Band_calc = True",
            text,
        )
        self.assertIn(
            "review Band_path",
            text,
        )

    def test_build_workflow_lines_falls_back_for_unknown_calculation(
        self,
    ):
        settings = QEInputSettings(
            calculation='md',
        )

        text = '\n'.join(
            _build_workflow_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_calc = True",
            text,
        )
        self.assertIn(
            "Geo_optim = False",
            text,
        )
        self.assertIn(
            "DOS_calc = False",
            text,
        )
        self.assertIn(
            "Band_calc = False",
            text,
        )
        self.assertIn(
            (
                "# NOTICE: QE calculation = 'md' "
                "has no direct Nanoworks workflow mapping"
            ),
            text,
        )

    def test_build_config_lines_assigns_nbands_to_workflow_stage(
        self,
    ):
        args = SimpleNamespace(
            outdirname=None,
            xc=None,
        )

        cases = [
            (
                'scf',
                'Ground_num_of_bands = 32',
            ),
            (
                'relax',
                'Ground_num_of_bands = 32',
            ),
            (
                'nscf',
                'DOS_num_of_bands = 32',
            ),
            (
                'bands',
                'Band_num_of_bands = 32',
            ),
        ]

        for calculation, expected in cases:
            with self.subTest(
                calculation=calculation
            ):
                settings = QEInputSettings(
                    calculation=calculation,
                    nbands=32,
                )

                text = '\n'.join(
                    build_config_lines(
                        name='Silicon',
                        geom_filename='Silicon.cif',
                        settings=settings,
                        args=args,
                    )
                )

                self.assertIn(
                    expected,
                    text,
                )

    def test_build_kpoint_lines_preserves_parity_shift(self):
        settings = QEInputSettings(
            k_mesh=[
                4,
                5,
                6,
            ],
            k_shift=[
                1,
                0,
                1,
            ],
        )

        text = '\n'.join(
            _build_kpoint_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_kpts_x = 4",
            text,
        )
        self.assertIn(
            "Ground_kpts_y = 5",
            text,
        )
        self.assertIn(
            "Ground_kpts_z = 6",
            text,
        )
        self.assertIn(
            "Ground_gamma = False",
            text,
        )
        self.assertNotIn(
            '# NOTICE:',
            text,
        )

    def test_build_kpoint_lines_preserves_zero_shift(self):
        settings = QEInputSettings(
            k_mesh=[
                4,
                4,
                4,
            ],
            k_shift=[
                0,
                0,
                0,
            ],
        )

        text = '\n'.join(
            _build_kpoint_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_gamma = True",
            text,
        )
        self.assertNotIn(
            '# NOTICE:',
            text,
        )

    def test_build_kpoint_lines_approximates_mixed_shift(self):
        settings = QEInputSettings(
            k_mesh=[
                6,
                6,
                4,
            ],
            k_shift=[
                1,
                1,
                0,
            ],
        )

        text = '\n'.join(
            _build_kpoint_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_gamma = False",
            text,
        )
        self.assertIn(
            (
                "# NOTICE: QE k-point shift [1, 1, 0] "
                "cannot be represented exactly"
            ),
            text,
        )
        self.assertIn(
            "nearest available shift [1, 1, 1]",
            text,
        )

    def test_build_kpoint_lines_marks_unrepresentable_odd_shift(
        self,
    ):
        settings = QEInputSettings(
            k_mesh=[
                5,
                5,
                5,
            ],
            k_shift=[
                1,
                1,
                1,
            ],
        )

        text = '\n'.join(
            _build_kpoint_lines(
                settings
            )
        )

        self.assertIn(
            "Ground_gamma = True",
            text,
        )
        self.assertIn(
            "nearest available shift [0, 0, 0]",
            text,
        )

if __name__ == '__main__':
    unittest.main()
