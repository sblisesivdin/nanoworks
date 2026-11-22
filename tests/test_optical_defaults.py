# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Optical defaults tested with small UPF fixtures and no DFT execution."""

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

from ase import Atoms

from nanoworks.optical import (
    resolve_optical_occupation, optical_engine_occupation, resolve_qe_optical_bands,
)
from nanoworks.dftsolve import (
    DFTConfig, dftsolve as DFTSolver, check_dft_configuration, prepare_qe_dry_run, get_band_path,
)


class TestOpticalDefaults(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for symbol, valence in (('Si', 4), ('O', 6)):
            (self.root / f'{symbol}.upf').write_text(
                f'<UPF><PP_HEADER pseudo_type="NC" z_valence="{valence}"/></UPF>',
                encoding='utf-8',
            )

    def config(self, **changes):
        values = dict(Engine='QE', Optical_calc=True,
            bulk_configuration=Atoms('Si2', cell=[5.4, 5.4, 5.4], pbc=True),
            Pseudo_dir=str(self.root), Pseudopotentials={'Si': 'Si.upf'})
        values.update(changes)
        return DFTConfig(**values)

    def test_optical_occupation_inherits_ground_scheme_and_width(self):
        for engine in ('QE', 'GPAW'):
            for scheme in ('fixed', 'fermi-dirac', 'marzari-vanderbilt', 'methfessel-paxton'):
                with self.subTest(engine=engine, scheme=scheme):
                    config = self.config(Engine=engine, Opt_calc_type='RPA',
                                         Occupation_scheme=scheme, Smearing_width=.12)
                    resolved = resolve_optical_occupation(config)
                    self.assertEqual(resolved['scheme'], scheme)
                    self.assertEqual(resolved['width'], None if scheme == 'fixed' else .12)
                    self.assertIsNone(config.Opt_FD_smearing)
                    native = optical_engine_occupation(config)
                    if scheme == 'fixed':
                        self.assertEqual(native, 'fixed' if engine == 'QE' else {'name': 'fixed-uniform'})
                    else:
                        self.assertEqual(native, {'name': scheme, 'width': .12})

    def test_explicit_optical_fd_width_overrides_ground_occupation(self):
        config = self.config(Occupation_scheme='fixed', Opt_FD_smearing=.02)
        self.assertEqual(optical_engine_occupation(config), {'name': 'fermi-dirac', 'width': .02})
        # Preserve GPAW's existing explicit zero-width option.
        config = self.config(Engine='GPAW', Opt_calc_type='RPA', Opt_FD_smearing=0.)
        self.assertEqual(optical_engine_occupation(config), {'name': 'fermi-dirac', 'width': 0.})

    def test_auto_band_count_scales_with_composition_and_charge(self):
        for formula, charge, expected_electrons in (('Si2', 0., 8.), ('Si8', 0., 32.),
                                                   ('SiO2', 1., 15.), ('SiO2', -.5, 16.5)):
            with self.subTest(formula=formula, charge=charge):
                atoms = Atoms(formula)
                bands = resolve_qe_optical_bands(atoms,
                    {'Si': 'Si.upf', 'O': str(self.root / 'O.upf')}, self.root,
                    total_charge=charge)
                self.assertEqual(bands['electrons'], expected_electrons)
                self.assertGreater(bands['nbands'], bands['occupied_minimum'])
                if formula == 'Si8':
                    self.assertGreater(bands['nbands'], 8)
        polarized = resolve_qe_optical_bands(Atoms('Si2'), {'Si': 'Si.upf'}, self.root, spinpol=True)
        self.assertGreater(polarized['nbands'], 8)

    def test_explicit_band_count_is_preserved_or_rejected_before_qe(self):
        for count in (1, 4, True, 8.5):
            with self.subTest(count=count), self.assertRaisesRegex(ValueError, 'Opt_num_of_bands'):
                resolve_qe_optical_bands(Atoms('Si2'), {'Si': 'Si.upf'}, self.root, nbands=count)
        bands = resolve_qe_optical_bands(Atoms('Si2'), {'Si': 'Si.upf'}, self.root, nbands=17)
        self.assertEqual(bands['nbands'], 17)
        with self.assertRaisesRegex(ValueError, 'positive electron count'):
            resolve_qe_optical_bands(Atoms('Si2'), {'Si': 'Si.upf'}, self.root, total_charge=8.)

    def test_preflight_and_dry_run_share_automatic_bands_and_fixed_occupation(self):
        config = self.config(Occupation_scheme='fixed',
                             bulk_configuration=Atoms('Si8', cell=[5.4, 5.4, 5.4], pbc=True))
        self.assertIsNone(config.Opt_num_of_bands)
        report = check_dft_configuration(config, str(self.root / 'sample'),
            check_executables=False, check_saved_state=False)
        self.assertTrue(report['ok'])
        self.assertIn('20 bands', next(item['detail'] for item in report['checks']
                                      if item['name'] == 'optical-bands'))
        with patch('nanoworks.engine.qe.run_qe_program') as run:
            plan = prepare_qe_dry_run(config, str(self.root / 'sample'))
        run.assert_not_called()
        job = next(job for job in plan['jobs'] if job['id'] == 'optical-nscf')
        text = Path(job['input_file']).read_text(encoding='utf-8')
        self.assertIn('nbnd = 20', text)
        self.assertIn("occupations = 'fixed'", text)
        self.assertNotIn('degauss', text)
        self.assertEqual(job['metadata']['bands']['nbands'], 20)
        self.assertIsNone(config.Opt_num_of_bands)
        self.assertEqual(config.DOS_integration, 'tetrahedron')

    def test_preflight_rejects_occupied_only_basis_and_missing_valence_metadata(self):
        for config in (self.config(Opt_num_of_bands=4),
                       self.config(Pseudopotentials={'Si': 'missing.upf'})):
            with self.subTest(config=config):
                report = check_dft_configuration(config, str(self.root / 'sample'),
                    check_executables=False, check_saved_state=False)
                self.assertFalse(report['ok'])
                self.assertTrue(any(error['name'] == 'optical-bands' for error in report['errors']))

    def test_optical_execution_uses_the_same_automatic_choices_as_dry_run(self):
        config = self.config(Occupation_scheme='cold', Smearing_width=.12,
                            bulk_configuration=Atoms('Si8', cell=[5.4, 5.4, 5.4], pbc=True))
        solver = object.__new__(DFTSolver)
        solver.__dict__.update(vars(config))
        solver.config = config
        solver.struct = str(self.root / 'runtime')
        solver.parallel_cores = 1
        solver._write_optical_summary = Mock(return_value='summary.json')
        solver.engine = SimpleNamespace(
            validate_qe_xc=Mock(), has_qe_state=Mock(return_value=True),
            run_nscf=Mock(return_value={}),
            run_epsilon=Mock(return_value={'optical_data': {}}),
            write_epsilon_optical_data=Mock(return_value={}),
        )
        solver._opticalcalc_qe()
        arguments = solver.engine.run_nscf.call_args.kwargs
        self.assertEqual(arguments['nbands'], 20)
        self.assertEqual(arguments['occupation'], {'name': 'marzari-vanderbilt', 'width': .12})
        self.assertEqual(solver._resolved_optical_bands['nbands'], 20)
        self.assertIsNone(config.Opt_num_of_bands)

    def test_gpaw_fixed_density_receives_the_inherited_occupation(self):
        # Import the real adapter with lightweight GPAW stubs, without installing
        # or executing a numerical GPAW engine.
        gpaw = SimpleNamespace(GPAW=Mock(), PW=Mock(), MixerSum=Mock(), FermiDirac=Mock())
        with patch.dict('sys.modules', {'gpaw': gpaw,
                                      'gpaw.eigensolvers': SimpleNamespace(Davidson=Mock())}):
            path = Path(__file__).parents[1] / 'nanoworks/engine/gpaw.py'
            spec = importlib.util.spec_from_file_location('_optical_gpaw_adapter', path)
            adapter = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(adapter)
        for occupation in ({'name': 'fixed-uniform'}, {'name': 'marzari-vanderbilt', 'width': .12}):
            with self.subTest(occupation=occupation), patch.object(adapter, 'load_gpaw_calc') as load:
                adapter.prepare_optical_calc('sample.gpw', False, 'opt.txt', 16,
                    None, None, (2, 2, 2), True, occupation=occupation)
                self.assertEqual(load.return_value.fixed_density.call_args.kwargs['occupations'], occupation)
        gpaw.FermiDirac.assert_not_called()

    def test_gpaw_phonopy_path_matches_the_final_cell(self):
        atoms = Atoms('Si', cell=[4, 4, 4], pbc=True)
        qpoints, labels, connections = get_band_path(atoms, None, 5)
        self.assertTrue(len(qpoints))
        self.assertIn('X', labels)
        self.assertEqual(len(connections), len(qpoints))

