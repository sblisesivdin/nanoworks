# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
"""JTH installation and workflow checks; no network or DFT engine execution."""

import contextlib
import hashlib
import io
import json
from pathlib import Path
import shutil
import tarfile
import tempfile
import unittest
from unittest.mock import Mock, patch

from ase import Atoms

from nanoworks import cli, pseudos
from nanoworks.dftsolve import (DFTConfig, check_dft_configuration,
                               prepare_qe_dry_run, resolve_qe_pseudo_configuration,
                               dftsolve as DFTSolver)
from nanoworks.dftconverge import _build_static_energy_settings, _build_backend
from nanoworks.engine.qe import resolve_qe_xc_settings, validate_qe_xc
from nanoworks.occupations import resolve_engine_occupation


class TestJTHPseudos(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.specs = {xc: dict(spec, count=2) for xc, spec in pseudos.JTH_SETS.items()}
        self.archives = {}
        for xc in self.specs:
            self.make_archive(xc)
        self.home = patch('nanoworks.pseudos.Path.home', return_value=self.root).start()
        self.addCleanup(patch.stopall)
        patch.dict(pseudos.JTH_SETS, self.specs).start()
        self.download = patch('nanoworks.pseudos.urlretrieve', side_effect=self.download_archive).start()
        self.atoms = Atoms('Si', cell=[5, 5, 5], pbc=True)

    def make_archive(self, xc, fault=None):
        path = self.root / (xc + '.tar.gz')
        self.archives[xc] = path
        spec = self.specs[xc]
        with tarfile.open(path, 'w:gz') as archive:
            for symbol, valence in (('Si', 4), ('Ni', 10)):
                if fault == 'missing' and symbol == 'Ni':
                    continue
                name = f'{symbol}.{spec["filename_xc"]}-JTH.UPF'
                if fault == 'traversal':
                    name = '../' + name
                else:
                    name = 'JTH-table/' + name
                typ = 'NC' if fault == 'type' else 'PAW'
                data = (f'<UPF><PP_HEADER element="{symbol}" pseudo_type="{typ}" '
                        f'z_valence="{valence}" has_so="F"/></UPF>').encode()
                info = tarfile.TarInfo(name)
                info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
            # This optional variant must not become the default for Si.
            data = b'<UPF><PP_HEADER element="Si" pseudo_type="PAW" z_valence="12" has_so="F"/></UPF>'
            info = tarfile.TarInfo(f'JTH-table/Si.{spec["filename_xc"]}-JTH_sp.UPF')
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        spec['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()

    def download_archive(self, url, destination):
        xc = next(xc for xc, spec in self.specs.items() if spec['url'] == url)
        shutil.copyfile(self.archives[xc], destination)

    def install(self):
        return pseudos.install_qe_pseudopotentials(family='jth')

    def test_installs_both_tables_without_replacing_the_nc_default(self):
        result = self.install()
        self.assertEqual(set(result), {'pbe', 'lda'})
        for xc in result:
            self.assertEqual(result[xc]['count'], 2)
            manifest = pseudos.load_qe_pseudo_manifest(family='jth', xc=xc)
            self.assertEqual(manifest['family'], 'jth')
            self.assertEqual(manifest['source_sha256'], self.specs[xc]['sha256'])
            self.assertFalse(any('_sp' in f for f in manifest['files'].values()))
            config = DFTConfig(Engine='QE', Pseudo_family='jth', XC_calc=xc.upper())
            directory, files = resolve_qe_pseudo_configuration(config, self.atoms)
            self.assertEqual(config.Pseudo_xc, xc)
            self.assertEqual(directory, result[xc]['directory'])
            self.assertEqual(set(files), {'Si'})
        default = DFTConfig(Engine='QE')
        self.assertEqual(default.Pseudo_family, 'pseudodojo')
        self.assertEqual(default.Pseudo_xc, 'pbe')
        # Having JTH installed must not trigger an automatic fallback to PAW.
        with self.assertRaisesRegex(FileNotFoundError, 'pseudodojo'):
            resolve_qe_pseudo_configuration(default, self.atoms)

    def test_reuses_complete_tables_and_repairs_only_the_broken_table(self):
        result = self.install()
        self.download.reset_mock()
        self.assertTrue(all(r['skipped'] for r in self.install().values()))
        self.download.assert_not_called()
        (result['lda']['directory'] / 'Si.LDA_PW-JTH.UPF').unlink()
        repaired = self.install()
        self.assertTrue(repaired['pbe']['skipped'])
        self.assertFalse(repaired['lda']['skipped'])
        self.download.assert_called_once()

    def test_bad_checksum_does_not_install_files(self):
        self.specs['pbe']['sha256'] = '0' * 64
        with self.assertRaisesRegex(RuntimeError, 'SHA256 mismatch'):
            self.install()
        self.assertFalse(list((self.root / '.nanoworks').rglob('*.UPF')))
        self.assertFalse(list((self.root / '.nanoworks').rglob('manifest.json')))

    def test_wrong_element_registry_is_rejected_and_repaired(self):
        result = self.install()
        path = result['pbe']['manifest']
        manifest = json.loads(path.read_text())
        manifest['files']['Si'] = manifest['files']['Ni']
        path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, '--pseudo-family jth'):
            pseudos.resolve_qe_pseudopotentials(self.atoms, family='jth')
        repaired = self.install()
        self.assertFalse(repaired['pbe']['skipped'])
        self.assertTrue(repaired['lda']['skipped'])
        self.assertEqual(pseudos.resolve_qe_pseudopotentials(self.atoms, family='jth'),
                         {'Si': 'Si.GGA_PBE-JTH.UPF'})

    def test_incomplete_wrong_type_and_unsafe_archives_are_rejected(self):
        for fault, message in (('missing', 'expected 2'), ('type', 'Invalid scalar'),
                               ('traversal', 'Unsafe JTH archive')):
            with self.subTest(fault=fault):
                self.make_archive('pbe', fault=fault)
                with self.assertRaisesRegex(RuntimeError, message):
                    self.install()
                self.assertFalse(list((self.root / '.nanoworks').rglob('*.UPF')))

    def test_download_error_has_an_actionable_library_name(self):
        self.download.side_effect = pseudos.URLError('offline')
        with self.assertRaisesRegex(RuntimeError, 'download JTH PBE'):
            self.install()

    def test_missing_jth_reports_the_jth_install_command(self):
        with self.assertRaisesRegex(FileNotFoundError, '--pseudo-family jth'):
            pseudos.resolve_qe_pseudopotentials(self.atoms, family='jth')

    def test_xc_is_inferred_and_conflicting_references_are_rejected(self):
        for name, expected in (('PBE', 'pbe'), ('HSE06', 'pbe'), ('LDA', 'lda'),
                               ('LDA-PW', 'lda'), ('PW', 'lda')):
            with self.subTest(xc=name):
                self.assertEqual(pseudos.resolve_qe_pseudo_options(
                    family='jth', xc_calc=name)['xc'], expected)
        for method, reference in (('PBE', 'lda'), ('LDA', 'pbe')):
            with self.assertRaises(ValueError):
                pseudos.resolve_qe_pseudo_options(family='jth', xc_calc=method, xc=reference)
        with self.assertRaisesRegex(ValueError, "Pseudo_family='jth'"):
            DFTConfig(Engine='QE', XC_calc='LDA')
        self.assertEqual(resolve_qe_xc_settings('LDA', pseudo_xc='lda')['input_dft'], 'PW')

    def test_optics_soc_and_reduced_exchange_cutoff_stop_before_execution(self):
        with patch('nanoworks.engine.qe.run_qe_program') as run:
            for flags, error in (({'SOC_calc': True}, 'scalar-relativistic'),
                                 ({'Optical_calc': True}, 'epsilon.x'),
                                 ({'Pseudo_relativistic': 'full'}, 'scalar-relativistic'),
                                 ({'EXX_cutoff': 800}, 'EXX_cutoff')):
                with self.subTest(flags=flags), self.assertRaisesRegex((ValueError, NotImplementedError), error):
                    DFTConfig(Engine='QE', Pseudo_family='jth', **flags)
            run.assert_not_called()
        config = DFTConfig(Engine='QE', Pseudo_family='jth', Ground_calc=True,
                           bulk_configuration=self.atoms)
        config.Optical_calc = True
        report = check_dft_configuration(config, str(self.root / 'sample'),
                                         check_executables=False, check_saved_state=False)
        self.assertFalse(report['ok'])
        self.assertTrue(any('epsilon.x' in item['detail'] for item in report['checks']))

    def test_lda_preflight_and_dry_run_use_pw92_without_running_qe(self):
        self.install()
        config = DFTConfig(Engine='QE', Pseudo_family='jth', XC_calc='LDA',
                           Ground_calc=True, bulk_configuration=self.atoms,
                           Density_cutoff_ratio=8)
        with patch('nanoworks.engine.qe.run_qe_program') as run:
            report = check_dft_configuration(config, str(self.root / 'sample'),
                                             check_executables=False, check_saved_state=False)
            self.assertTrue(report['ok'], report['checks'])
            plan = prepare_qe_dry_run(config, self.root / 'sample')
            run.assert_not_called()
        deck = Path(plan['jobs'][0]['input_file']).read_text()
        self.assertIn("input_dft = 'PW'", deck)
        self.assertIn('Si.LDA_PW-JTH.UPF', deck)

    def test_convergence_keeps_the_jth_lda_reference(self):
        config = {'Pseudo_family': 'jth', 'XC_calc': 'LDA'}
        settings = _build_static_energy_settings(config, self.atoms, 'QE')
        self.assertEqual(settings['pseudo_xc'], 'lda')
        self.install()
        loader = Mock()
        plan = type('Plan', (), {'engine': 'QE'})()
        _build_backend(config, plan, self.atoms, backend_loader=loader)
        self.assertEqual(loader.call_args.kwargs['pseudopotentials'], {'Si': 'Si.LDA_PW-JTH.UPF'})
        config['EXX_cutoff'] = 800
        with self.assertRaisesRegex(NotImplementedError, 'EXX_cutoff'):
            _build_backend(config, plan, self.atoms, backend_loader=loader)

    def test_runtime_ground_and_relaxation_forward_the_lda_reference(self):
        self.install()
        for relax in (False, True):
            with self.subTest(relax=relax):
                config = DFTConfig(Engine='QE', Pseudo_family='jth', XC_calc='LDA',
                                   Ground_calc=True, Geo_optim=relax, bulk_configuration=self.atoms)
                solver = object.__new__(DFTSolver)
                solver.__dict__.update(vars(config))
                solver.struct = str(self.root / 'sample')
                solver.parallel_cores = 1
                solver.Occupation = resolve_engine_occupation('QE')
                solver.engine = Mock(validate_qe_xc=validate_qe_xc)
                runner = solver.engine.run_relax if relax else solver.engine.run_scf
                runner.side_effect = RuntimeError('stop before the external calculation')
                with contextlib.redirect_stdout(io.StringIO()), self.assertRaisesRegex(RuntimeError, 'stop before'):
                    solver._groundcalc_qe()
                self.assertEqual(runner.call_args.kwargs['pseudo_xc'], 'lda')
                self.assertEqual(runner.call_args.kwargs['pseudopotentials'], {'Si': 'Si.LDA_PW-JTH.UPF'})

    def test_lda_dry_run_keeps_one_reference_across_electronic_and_response_stages(self):
        self.install()
        config = DFTConfig(Engine='QE', Pseudo_family='jth', XC_calc='LDA',
                           Ground_calc=True, Geo_optim=True, DOS_calc=True, Band_calc=True,
                           Elastic_calc=True, Phonon_calc=True, bulk_configuration=self.atoms,
                           Density_cutoff_ratio=8)
        with patch('nanoworks.engine.qe.run_qe_program') as run:
            plan = prepare_qe_dry_run(config, self.root / 'combined')
            run.assert_not_called()
        decks = [Path(job['input_file']).read_text() for job in plan['jobs']
                 if job['executable'] in ('pw.x', 'thermo_pw.x')]
        self.assertEqual(len(decks), 4)
        for text in decks:
            self.assertIn("input_dft = 'PW'", text)
            self.assertIn('Si.LDA_PW-JTH.UPF', text)

    def test_cli_installs_selected_family_and_reports_paw(self):
        output = io.StringIO()
        with patch('sys.argv', ['nanoworks', '--install-qe-pseudos', '--pseudo-family', 'jth']), \
                contextlib.redirect_stdout(output), self.assertRaises(SystemExit) as exit_status:
            cli.main()
        self.assertEqual(exit_status.exception.code, 0)
        self.assertIn('PAW UPF', output.getvalue())
        self.assertIn('still uses PseudoDojo NC', output.getvalue())
        with patch('sys.argv', ['nanoworks', '--install-qe-pseudos']), \
                patch('nanoworks.cli.install_qe_pseudopotentials', return_value={}) as install, \
                patch('nanoworks.cli.format_qe_pseudo_installation', return_value='ready'), \
                contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit):
            cli.main()
        install.assert_called_once_with(family='pseudodojo')
