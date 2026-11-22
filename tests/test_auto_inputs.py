# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import ast
import tempfile
import unittest
import sys
from pathlib import Path

from ase import Atoms
from ase.build import bulk
from ase.calculators.calculator import kptdensity2monkhorstpack
from ase.io import write

from nanoworks.dftsolve import build_auto_input_lines, struct_from_auto, struct_from_file


class TestAutoInputs(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def generate(self, atoms, name='sample', overrides=None):
        geometry = self.root / f'{name}.cif'
        write(geometry, atoms)
        struct, config = struct_from_auto(geometry, overrides=overrides,
                                         report_structure=False)
        path = Path(struct + '-CONFIG-NANOWORKS-Input-Auto.py')
        return geometry, path, config

    def test_auto_mesh_uses_reciprocal_density_and_shrinks_for_supercells(self):
        primitive = bulk('Si', 'diamond', a=5.43, cubic=True)
        meshes = []
        for index, atoms in enumerate((primitive, primitive.repeat((2, 1, 1)))):
            _, _, config = self.generate(atoms, f'cell{index}')
            mesh = (config.Ground_kpts_x, config.Ground_kpts_y, config.Ground_kpts_z)
            expected = kptdensity2monkhorstpack(atoms, 3., even=None)
            self.assertEqual(mesh, tuple(expected))
            meshes.append(mesh)
        self.assertLess(meshes[1][0], meshes[0][0])
        self.assertEqual(meshes[1][1:], meshes[0][1:])

    def test_skewed_cell_mesh_and_vacuum_direction(self):
        atoms = bulk('Si', 'diamond', a=5.43)
        _, _, config = self.generate(atoms, 'skewed')
        self.assertEqual((config.Ground_kpts_x, config.Ground_kpts_y, config.Ground_kpts_z),
                         tuple(kptdensity2monkhorstpack(atoms, 3., even=None)))
        _, _, slab = self.generate(Atoms('Si', cell=[5, 5, 20], pbc=True), 'slab')
        self.assertEqual(slab.Ground_kpts_z, 1)
        self.assertEqual(slab.Electrostatic_boundary, 'periodic')

    def test_minimal_generated_input_round_trips_with_the_geometry(self):
        geometry, path, config = self.generate(bulk('Si', 'diamond', a=5.43),
                                              overrides={'Engine': 'QE'})
        text = path.read_text()
        ast.parse(text)
        for keyword in ('Pseudo_family', 'Pseudo_xc', 'EXX_cutoff', 'XC_omega',
                        'Density_cutoff_ratio', 'Geometry_optimizer', 'Mode'):
            self.assertNotIn(keyword + ' =', text)
        self.assertNotIn('import ', text)
        _, loaded = struct_from_file(path, geometry, create_output=False,
                                     report_structure=False)
        self.assertEqual(build_auto_input_lines(config), build_auto_input_lines(loaded))
        self.assertTrue(loaded.Ground_calc and loaded.DOS_calc and loaded.Band_calc)

    def test_explicit_default_overrides_are_retained_once_and_arrays_are_literals(self):
        overrides = {'Engine': 'QE', 'DOS_calc': False, 'Gamma': True,
                     'Density_cutoff_ratio': None, 'Ground_kpts_x': 7,
                     'Phonon_supercell': [[3, 0, 0], [0, 2, 0], [0, 0, 1]]}
        geometry, path, config = self.generate(bulk('Si', 'diamond', a=5.43),
                                              overrides=overrides)
        text = path.read_text()
        for name in overrides:
            self.assertEqual(text.count(name + ' ='), 1)
        _, loaded = struct_from_file(path, geometry, create_output=False,
                                     report_structure=False)
        self.assertEqual(build_auto_input_lines(config, overrides),
                         build_auto_input_lines(loaded, overrides))
        self.assertFalse(loaded.DOS_calc)
        self.assertEqual(loaded.Ground_kpts_x, 7)

    def test_loading_same_input_after_edit_and_same_name_in_another_folder(self):
        for dirname, cutoff in (('first', 400), ('second', 600), ('first', 500)):
            directory = self.root / dirname
            directory.mkdir(exist_ok=True)
            path = directory / 'input.py'
            path.write_text(f"Engine = 'QE'\nWavefunction_cutoff = {cutoff}\n", encoding='utf-8')
            _, config = struct_from_file(path, None, create_output=False, report_structure=False)
            self.assertEqual(config.Wavefunction_cutoff, cutoff)

    def test_input_local_imports_work_and_search_path_is_restored_on_errors(self):
        (self.root / '_auto_input_helper.py').write_text('cutoff = 525\n', encoding='utf-8')
        path = self.root / 'local.py'
        path.write_text("from _auto_input_helper import cutoff\n"
                        "Engine = 'QE'\nWavefunction_cutoff = cutoff\n", encoding='utf-8')
        before = sys.path[:]
        _, config = struct_from_file(path, None, create_output=False, report_structure=False)
        self.assertEqual(config.Wavefunction_cutoff, 525)
        self.assertEqual(sys.path, before)
        path.write_text("raise RuntimeError('input failed')\n", encoding='utf-8')
        with self.assertRaisesRegex(RuntimeError, 'input failed'):
            struct_from_file(path, None, create_output=False, report_structure=False)
        self.assertEqual(sys.path, before)
