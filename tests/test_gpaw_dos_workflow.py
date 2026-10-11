# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression checks for GPAW result exports using in-memory calculators."""

import importlib
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np


class TestGPAWDOSWorkflow(unittest.TestCase):
    def setUp(self):
        with patch.object(sys, 'argv', [sys.argv[0]]):
            self.module = importlib.import_module('nanoworks.dftsolve')
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.calc = Mock()
        self.calc.get_fermi_level.return_value = 2.0
        self.engine = Mock()
        self.engine.is_hybrid.return_value = False
        self.engine.prepare_dos_calc.return_value = self.calc
        self.rawdos = Mock()
        self.rawdos.get_energies.return_value = np.array([3., 4., 5.])
        self.rawdos.raw_pdos.side_effect = lambda energies, **kw: np.full(
            len(energies), .1 if kw['spin'] in (0, None) else .2)
        self.rawdos.raw_dos.side_effect = lambda energies, **kw: np.full(
            len(energies), {None: 30., 0: 10., 1: 20.}[kw.get('spin')])
        self.factory = Mock(return_value=self.rawdos)
        self.gpaw_dos = types.ModuleType('gpaw.dos')
        self.gpaw_dos.DOSCalculator = types.SimpleNamespace(from_calculator=self.factory)
        self.atoms = Mock()
        self.atoms.get_chemical_symbols.return_value = ['Si', 'Si']
        self.atoms.get_global_number_of_atoms.return_value = 2
        self.solver = types.SimpleNamespace(
            engine=self.engine, XC_calc='PBE', Gamma=False, Ground_gamma=None,
            DOS_kpts_density=3., DOS_kpts_x=0, DOS_kpts_y=0, DOS_kpts_z=0,
            DOS_gamma=False, Ground_kpts_density=2., Ground_kpts_x=0,
            Ground_kpts_y=0, Ground_kpts_z=0, DOS_integration='smearing',
            DOS_width=.001, DOS_num_of_bands=4, DOS_npoints=3, Occupation={},
            struct=str(Path(self.directory.name) / 'Si'), Engine='GPAW',
            scf_settings={'convergence': {}}, bulk_configuration=self.atoms,
            Spin_calc=False, SOC_calc=False, Energy_min=-1., Energy_max=1.,
            hybrid_fermi_level=Mock(return_value=2.), _t=lambda key: key)

    def run_dos(self):
        with patch.dict(sys.modules, {'gpaw': types.ModuleType('gpaw'), 'gpaw.dos': self.gpaw_dos}), \
                patch.object(self.module, 'world', types.SimpleNamespace(rank=1)), \
                patch.object(self.module, 'paropen', open), \
                patch.object(self.module, 'parprint'), \
                patch.object(self.module, 'resolve_stage_kpoint_settings', return_value=(3., None, False)):
            self.module.dftsolve._doscalc_gpaw(self.solver)

    def table(self, suffix):
        return np.loadtxt(self.solver.struct + '-DOS-GPAW-Result-' + suffix + '.dat', skiprows=1)

    def test_pdos_uses_prepared_calculator_and_exports_only_tabular_rows(self):
        for spin, hybrid in [(False, False), (True, False), (False, True)]:
            with self.subTest(spin=spin, hybrid=hybrid):
                self.solver.Spin_calc = spin
                self.engine.is_hybrid.return_value = hybrid
                self.factory.reset_mock()
                self.run_dos()
                self.factory.assert_called_once_with(
                    self.calc, soc=False, theta=0.0, phi=0.0, shift_fermi_level=False)
                suffixes = ['-Down', '-Up'] if spin else ['']
                for suffix in suffixes:
                    np.testing.assert_array_equal(self.table('PDOS' + suffix)[:, 0], [1., 2., 3.])
                    rows = Path(self.solver.struct + '-DOS-GPAW-Result-Raw-PDOS-Each-Atom'
                                + suffix + '.dat').read_text().splitlines()
                    self.assertEqual(len(rows), 7)
                    self.assertTrue(all(len(row.split()) == len(rows[0].split()) for row in rows))
                    self.assertEqual([row.split()[0] for row in rows[1:]], ['1'] * 3 + ['2'] * 3)

    def test_total_dos_is_not_the_sum_of_atomic_projections(self):
        for spin in (False, True):
            with self.subTest(spin=spin):
                self.solver.Spin_calc = spin
                self.run_dos()
                expected = [('DOS-Down', 'PDOS-Down', 20.), ('DOS-Up', 'PDOS-Up', 10.)] if spin else [('DOS', 'PDOS', 30.)]
                for dos_name, pdos_name, value in expected:
                    np.testing.assert_array_equal(self.table(dos_name)[:, 1], np.full(3, value))
                    self.assertTrue(np.all(self.table(pdos_name)[:, -1] < value))
                for call in self.rawdos.raw_dos.call_args_list:
                    self.assertEqual(call.kwargs['width'], .001)

    def test_tetrahedron_passes_zero_width_to_both_dos_and_pdos(self):
        self.solver.DOS_integration = 'tetrahedron'
        self.run_dos()
        for mock in (self.rawdos.raw_dos, self.rawdos.raw_pdos):
            self.assertTrue(mock.called)
            self.assertTrue(all(call.kwargs['width'] == 0. for call in mock.call_args_list))

    def test_soc_uses_gpaw_weighted_dos_and_exact_requested_width(self):
        self.solver.SOC_calc = True
        self.run_dos()
        self.factory.assert_called_once_with(
            self.calc, soc=True, theta=0.0, phi=0.0, shift_fermi_level=False)
        np.testing.assert_array_equal(self.table('DOS-SOC')[:, 0], [-1., 0., 1.])
        np.testing.assert_array_equal(self.table('DOS-SOC')[:, 1], [30., 30., 30.])
        call = self.rawdos.raw_dos.call_args
        np.testing.assert_array_equal(call.args[0], [1., 2., 3.])
        self.assertEqual(call.kwargs, {'spin': None, 'width': .001})
        self.rawdos.raw_pdos.assert_not_called()

    def test_soc_invalid_reference_or_spectrum_preserves_previous_table(self):
        self.solver.SOC_calc = True
        output = Path(self.solver.struct + '-DOS-GPAW-Result-DOS-SOC.dat')
        for reference, values in [(np.nan, [1., 2., 3.]), (np.inf, [1., 2., 3.]),
                                  (None, [1., 2., 3.]), (2., [1., np.nan, 3.]),
                                  (2., [1., 2.])]:
            with self.subTest(reference=reference, values=values):
                output.write_bytes(b'previous spectrum')
                self.calc.get_fermi_level.return_value = reference
                self.rawdos.raw_dos.side_effect = None
                self.rawdos.raw_dos.return_value = np.array(values)
                with self.assertRaises(ValueError):
                    self.run_dos()
                self.assertEqual(output.read_bytes(), b'previous spectrum')

    def test_spin_pdos_labels_match_gpaw_channels(self):
        self.solver.Spin_calc = True
        self.run_dos()
        np.testing.assert_allclose(self.table('PDOS-Up')[:, 1], .2)
        np.testing.assert_allclose(self.table('PDOS-Down')[:, 1], .4)

    def test_band_spin_files_match_the_same_channels(self):
        from nanoworks.result_plotting import load_band_table
        self.solver.Spin_calc = True
        self.solver.Band_path = 'GX'
        self.solver.Band_npoints = 99
        self.solver.Band_num_of_bands = 2
        self.solver.Projected_band_plot = False
        bs = types.SimpleNamespace(path=types.SimpleNamespace(
            kpts=np.zeros((3, 3)), get_linear_kpoint_axis=lambda: ([0., .5, 1.], [0., 1.], ['G', 'X'])),
            write=Mock())
        self.calc.band_structure.return_value = bs
        self.calc.get_number_of_bands.return_value = 2
        self.calc.get_eigenvalues.side_effect = lambda k, spin: np.array([k + spin * 10., k + spin * 10. + 1.])
        self.engine.prepare_band_calc.return_value = self.calc
        saved = types.SimpleNamespace(write=Mock(), plot=Mock())
        factory = Mock(return_value=saved)
        ase_bands = types.ModuleType('ase.spectrum.band_structure')
        ase_bands.BandStructure = factory
        with patch.dict(sys.modules, {'ase.spectrum.band_structure': ase_bands}), \
                patch.object(self.module, 'world', types.SimpleNamespace(rank=1)), \
                patch.object(self.module, 'paropen', open), patch.object(self.module, 'parprint'), \
                patch.object(self.module, 'run_optical_exports', side_effect=lambda callback: callback()):
            self.module.dftsolve._bandcalc_gpaw(self.solver)
        self.assertEqual(factory.call_args.kwargs['reference'], 0.)
        self.assertEqual(saved.plot.call_args.kwargs['emin'], -1.)
        self.assertEqual(saved.plot.call_args.kwargs['emax'], 1.)
        np.testing.assert_array_equal(factory.call_args.kwargs['energies'][0],
                                      [[-2., -1.], [-1., 0.], [0., 1.]])
        for label, spin in [('Up', 0), ('Down', 1)]:
            expected = np.array([[k + spin * 10. - 2., k + spin * 10. - 1.] for k in range(3)])
            for suffix in ['', '-XYYY']:
                _, values = load_band_table(self.solver.struct + '-BAND-GPAW-Result-Band-' + label + suffix + '.dat')
                np.testing.assert_array_equal(values, expected)

    def test_density_spin_labels_and_polarization_have_correct_sign(self):
        self.solver.Spin_calc = True
        self.solver.Refine_grid = 2
        self.engine.load_gpaw_calc.return_value = self.calc
        self.calc.get_pseudo_density.side_effect = lambda spin=None: np.full((2, 2, 2), {None: 4., 0: 3., 1: 1.}[spin])
        self.calc.get_all_electron_density.side_effect = lambda spin=None, **kw: np.full((2, 2, 2), {None: 4., 0: 3., 1: 1.}[spin])
        with patch.object(self.module, 'paropen', open), patch.object(self.module, 'parprint'), \
                patch.object(self.module, 'write') as write:
            self.module.dftsolve._densitycalc_gpaw(self.solver)
        for call in write.call_args_list:
            name = str(call.args[0])
            expected = 3. if '-Up.cube' in name else 1. if '-Down.cube' in name else .5 if '-Zeta.cube' in name else 4.
            scale = 1. if '-Zeta.cube' in name else self.module.Bohr ** 3
            np.testing.assert_allclose(call.kwargs['data'], expected * scale)

    def test_density_polarization_is_finite_in_vacuum(self):
        self.solver.Spin_calc = True
        self.solver.Refine_grid = 2
        self.engine.load_gpaw_calc.return_value = self.calc
        self.calc.get_pseudo_density.side_effect = lambda **kw: np.zeros((2, 2, 2))
        self.calc.get_all_electron_density.side_effect = lambda **kw: np.zeros((2, 2, 2))
        with patch.object(self.module, 'paropen', open), patch.object(self.module, 'parprint'), \
                patch.object(self.module, 'write') as write, np.errstate(divide='raise', invalid='raise'):
            self.module.dftsolve._densitycalc_gpaw(self.solver)
        for call in write.call_args_list:
            self.assertTrue(np.isfinite(call.kwargs['data']).all())
            np.testing.assert_array_equal(call.kwargs['data'], np.zeros((2, 2, 2)))

    def test_hybrid_reference_accepts_zero_and_rejects_unavailable_values(self):
        self.solver.Ground_fermi_level = 0.
        self.assertEqual(self.module.dftsolve.hybrid_fermi_level(self.solver), 0.)
        self.engine.load_gpaw_calc.assert_not_called()
        self.solver.Ground_fermi_level = np.inf
        ground = Mock()
        ground.get_fermi_level.return_value = np.nan
        self.engine.load_gpaw_calc.return_value = ground
        self.calc.get_fermi_level.return_value = -3.
        self.assertEqual(self.module.dftsolve.hybrid_fermi_level(self.solver, self.calc), -3.)
        self.calc.get_fermi_level.return_value = np.inf
        with self.assertRaisesRegex(ValueError, 'finite'):
            self.module.dftsolve.hybrid_fermi_level(self.solver, self.calc)

    def test_band_invalid_reference_stops_before_result_exports(self):
        self.engine.prepare_band_calc.return_value = self.calc
        self.solver.Band_path = 'GX'
        self.solver.Band_npoints = 3
        self.solver.Band_num_of_bands = 2
        for reference in [None, np.nan, np.inf]:
            self.calc.get_fermi_level.return_value = reference
            with patch.object(self.module, 'parprint'), self.assertRaisesRegex(ValueError, 'finite'):
                self.module.dftsolve._bandcalc_gpaw(self.solver)
        self.assertFalse(list(Path(self.directory.name).glob('*BAND*')))

    def test_magnetic_soc_bands_export_soc_eigenstates(self):
        from nanoworks.result_plotting import load_band_table
        self.solver.Spin_calc = self.solver.SOC_calc = True
        self.solver.Band_path = 'GX'
        self.solver.Band_npoints = 3
        self.solver.Band_num_of_bands = 2
        self.solver.Projected_band_plot = False
        bs = types.SimpleNamespace(path=types.SimpleNamespace(kpts=np.zeros((3, 3)),
            get_linear_kpoint_axis=lambda: ([0., .5, 1.], [0., 1.], ['G', 'X'])))
        self.calc.band_structure.return_value = bs
        self.calc.get_number_of_bands.return_value = 2
        self.engine.prepare_band_calc.return_value = self.calc
        values = np.arange(12.).reshape(3, 4)
        soc_module = types.ModuleType('gpaw.spinorbit')
        soc_module.soc_eigenstates = Mock(return_value=types.SimpleNamespace(eigenvalues=lambda: values))
        ase_module = types.ModuleType('ase.spectrum.band_structure')
        ase_module.BandStructure = Mock(return_value=types.SimpleNamespace(write=Mock(), plot=Mock()))
        with patch.dict(sys.modules, {'gpaw': types.ModuleType('gpaw'), 'gpaw.spinorbit': soc_module,
                                     'ase.spectrum.band_structure': ase_module}), \
                patch.object(self.module, 'paropen', open), patch.object(self.module, 'parprint'), \
                patch.object(self.module, 'run_optical_exports', side_effect=lambda callback: callback()):
            self.module.dftsolve._bandcalc_gpaw(self.solver)
        _, actual = load_band_table(self.solver.struct + '-BAND-GPAW-Result-Band-XYYY.dat')
        np.testing.assert_array_equal(actual, values - 2.)
        self.calc.get_eigenvalues.assert_not_called()
        self.assertFalse(list(Path(self.directory.name).glob('*Band-Up*')))

    def test_soc_projected_bands_stop_before_engine_preparation(self):
        self.solver.SOC_calc = self.solver.Projected_band_plot = True
        with self.assertRaisesRegex(NotImplementedError, 'SOC projected'):
            self.module.dftsolve._bandcalc_gpaw(self.solver)
        self.engine.prepare_band_calc.assert_not_called()
        with self.assertRaisesRegex(NotImplementedError, 'SOC projected'):
            self.module.dftsolve._draw_projectedband(self.solver, self.calc)
