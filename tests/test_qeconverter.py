import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from nanoworks.qeconverter import (
    QEInputSettings,
    build_config_lines,
    determine_system_name,
    parse_qe_input,
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
            "Gamma = True",
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

if __name__ == '__main__':
    unittest.main()
