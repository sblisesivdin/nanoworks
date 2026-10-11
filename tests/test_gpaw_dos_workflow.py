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
            hybrid_fermi_level=Mock(return_value=2.))

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
