# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Cutoff policy checks using UPF headers, without running QE."""

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ase import Atoms

from nanoworks.cutoffs import EV_PER_RYDBERG, resolve_qe_cutoff_settings
from nanoworks.engine.qe import render_scf_input, render_nscf_input
from nanoworks.pseudos import read_upf_cutoff_metadata


class TestQECutoffPolicy(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def upf(self, name='Si.upf', attributes='pseudo_type="NC"'):
        path = self.root / name
        path.write_text('<UPF><PP_HEADER ' + attributes + '/></UPF>', encoding='utf-8')
        return path

    def settings(self, cutoff_ry=50, ratio=None, pseudos=None):
        return resolve_qe_cutoff_settings(
            cutoff_ry * EV_PER_RYDBERG,
            pseudos or {'Si': 'Si.upf'}, self.root, ratio,
        )

    def test_header_hints_use_rydbergs_and_fortran_exponents(self):
        path = self.upf(attributes="pseudo_type='NC' wfc_cutoff='4.0D1' rho_cutoff='1.6d2'")
        metadata = read_upf_cutoff_metadata(path)
        self.assertEqual(metadata, {'pseudo_type': 'NC', 'wfc_cutoff': 40, 'rho_cutoff': 160})

    def test_unusable_hints_and_text_outside_header_are_ignored(self):
        for hint in ('0', '-3', 'nan', 'inf', 'bad'):
            with self.subTest(hint=hint):
                path = self.upf(attributes=f'wfc_cutoff="{hint}" rho_cutoff="{hint}"')
                self.assertIsNone(read_upf_cutoff_metadata(path)['rho_cutoff'])
        path.write_text('<PP_INFO wfc_cutoff="40" rho_cutoff="500"/>', encoding='utf-8')
        self.assertIsNone(read_upf_cutoff_metadata(path)['rho_cutoff'])

    def test_norm_conserving_without_hints_retains_fourfold_baseline(self):
        self.upf()
        settings = self.settings()
        self.assertEqual(settings['density_ratio'], 4)
        self.assertEqual(settings['notices'], [])

    def test_mixture_uses_largest_recommended_ratio_and_absolute_floor(self):
        self.upf(attributes='pseudo_type="NC" wfc_cutoff="40" rho_cutoff="160"')
        oxygen = self.upf('O.upf', 'pseudo_type="US" wfc_cutoff="60" rho_cutoff="720"')
        pseudos = {'Si': 'Si.upf', 'O': str(oxygen)}
        low = self.settings(cutoff_ry=30, pseudos=pseudos)
        high = self.settings(cutoff_ry=100, pseudos=pseudos)
        self.assertAlmostEqual(low['density_ev'] / EV_PER_RYDBERG, 720)
        self.assertEqual(high['density_ratio'], 12)
        self.assertAlmostEqual(low['wavefunction_ev'] / EV_PER_RYDBERG, 30)
        self.assertEqual(len(low['notices']), 2)
        self.assertEqual(high['notices'], [])

    def test_density_hint_without_wavefunction_hint_is_an_absolute_floor(self):
        self.upf(attributes='pseudo_type="PAW" rho_cutoff="600"')
        self.assertAlmostEqual(self.settings()['density_ratio'], 12)

    def test_paw_and_ultrasoft_flags_override_nc_label(self):
        for flag in ('is_paw="T"', 'is_ultrasoft=".true."'):
            with self.subTest(flag=flag):
                self.upf(attributes='pseudo_type="NC" ' + flag)
                settings = self.settings()
                self.assertTrue(settings['notices'])
                self.assertIn('no usable UPF density-cutoff', settings['notices'][0])

    def test_unknown_or_missing_metadata_reports_fallback(self):
        for content in (None, '<UPF><PP_HEADER/></UPF>', '<PP_HEADER>old UPF</PP_HEADER>'):
            with self.subTest(content=content):
                if content is not None:
                    (self.root / 'Si.upf').write_text(content, encoding='utf-8')
                settings = self.settings()
                self.assertEqual(settings['density_ratio'], 4)
                self.assertTrue(settings['notices'])

    def test_explicit_ratio_is_preserved_without_reading_hints(self):
        self.upf(attributes='pseudo_type="US" wfc_cutoff="60" rho_cutoff="720"')
        with patch('nanoworks.pseudos.read_upf_cutoff_metadata') as read:
            settings = self.settings(ratio=8)
        read.assert_not_called()
        self.assertEqual(settings['density_ratio'], 8)
        self.assertEqual(settings['source'], 'explicit Density_cutoff_ratio')

    def test_renderers_share_auto_policy_and_record_actual_cutoffs(self):
        self.upf(attributes='pseudo_type="US" wfc_cutoff="40" rho_cutoff="400"')
        for render in (render_scf_input, render_nscf_input):
            with self.subTest(render=render.__name__):
                text = render(
                    atoms=Atoms('Si', cell=[5, 5, 5], pbc=True),
                    pseudopotentials={'Si': 'Si.upf'}, pseudo_dir=self.root,
                    cutoff_ev=50 * EV_PER_RYDBERG, kpoint_size=(2, 2, 2),
                )
                self.assertIn('ecutrho = 500', text)
                self.assertIn('ratio=10 (automatic UPF policy)', text)

    def test_preflight_and_dry_run_agree_and_preserve_config(self):
        from nanoworks.dftsolve import DFTConfig, check_dft_configuration, prepare_qe_dry_run

        self.upf(attributes='pseudo_type="US" wfc_cutoff="40" rho_cutoff="400"')
        config = DFTConfig(
            Engine='QE', Ground_calc=True,
            bulk_configuration=Atoms('Si', cell=[5, 5, 5], pbc=True),
            Pseudo_dir=str(self.root), Pseudopotentials={'Si': 'Si.upf'},
            Wavefunction_cutoff=50 * EV_PER_RYDBERG,
        )
        report = check_dft_configuration(
            config, str(self.root / 'sample'),
            check_executables=False, check_saved_state=False,
        )
        self.assertTrue(report['ok'])
        density = next(check for check in report['checks'] if check['name'] == 'density-cutoff')
        self.assertIn('ratio 10', density['detail'])
        with patch('nanoworks.engine.qe.run_qe_program') as run:
            plan = prepare_qe_dry_run(config, str(self.root / 'sample'))
        run.assert_not_called()
        text = Path(plan['jobs'][0]['input_file']).read_text(encoding='utf-8')
        self.assertIn('ecutrho = 500', text)
        self.assertIsNone(config.Density_cutoff_ratio)
        self.assertEqual(config.Wavefunction_cutoff, 50 * EV_PER_RYDBERG)

    def test_convergence_runner_resolves_each_cutoff_and_preserves_overrides(self):
        from nanoworks.dftconverge import _build_static_energy_settings
        from nanoworks.convergence_backends.qe import QEStaticEnergyBackend

        self.upf(attributes='pseudo_type="US" wfc_cutoff="40" rho_cutoff="400"')
        atoms = Atoms('Si', cell=[5, 5, 5], pbc=True)
        backend = QEStaticEnergyBackend({'Si': 'Si.upf'}, self.root)
        for cutoff, ratio, expected_rho in ((50, None, 500), (80, None, 800), (50, 6, 300)):
            with self.subTest(cutoff=cutoff, ratio=ratio):
                config = {} if ratio is None else {'Density_cutoff_ratio': ratio}
                settings = _build_static_energy_settings(config, atoms, 'QE')
                with patch('nanoworks.engine.qe.run_qe_program') as run, patch(
                    'nanoworks.engine.qe.parse_pw_output', return_value={
                        'total_energy_ev': -42., 'qe_version': (7, 4, 1), 'job_done': True,
                    },
                ):
                    result = backend.calculate_static_energy(
                        atoms, cutoff_ev=cutoff * EV_PER_RYDBERG,
                        kpoint_settings={'size': (2, 2, 2)},
                        workdir=self.root / f'point-{cutoff}-{ratio}',
                        settings=settings, parallel_cores=1,
                    )
                run.assert_called_once()
                text = Path(run.call_args.kwargs['input_file']).read_text(encoding='utf-8')
                self.assertIn(f'ecutrho = {expected_rho}', text)
                self.assertEqual(result.total_energy_ev, -42.)
                self.assertEqual(settings['density_cutoff_ratio'], ratio)

