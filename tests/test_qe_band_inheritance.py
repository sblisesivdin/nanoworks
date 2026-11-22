# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from ase import Atoms

from nanoworks.dftsolve import DFTConfig, prepare_qe_dry_run
from nanoworks.engine.qe import (
    resolve_qe_stage_nbands, run_nscf, run_bands, build_band_path,
)


class TestQEBandInheritance(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.struct = self.root / 'sample'
        self.state = Path(str(self.struct) + '-GROUND-QE-Result-State')
        self.schema = self.state / 'nanoworks.save/data-file-schema.xml'
        self.schema.parent.mkdir(parents=True)
        (self.root / 'Si.upf').write_text(
            '<UPF><PP_HEADER pseudo_type="NC" z_valence="4"/></UPF>', encoding='utf-8')
        self.atoms = Atoms('Si2', cell=[5.4, 5.4, 5.4], pbc=True)

    def write_state(self, fields='<nbnd>32</nbnd>'):
        self.schema.write_text('<q:espresso xmlns:q="urn:qe"><q:input><nbnd>99</nbnd>'
            '</q:input><q:output><q:band_structure>' + fields
            + '</q:band_structure></q:output></q:espresso>', encoding='utf-8')

    def config(self, **changes):
        values = dict(Engine='QE', DOS_calc=True, Band_calc=True, Band_path='GX', Band_npoints=5,
            bulk_configuration=self.atoms, Pseudo_dir=str(self.root),
            Pseudopotentials={'Si': 'Si.upf'}, Ground_kpts_x=1, Ground_kpts_y=1, Ground_kpts_z=1)
        values.update(changes)
        return DFTConfig(**values)

    def test_only_actual_output_counts_are_inherited_and_overrides_win(self):
        self.write_state()
        self.assertEqual(resolve_qe_stage_nbands(None, self.state, ground_nbands=12), 32)
        self.assertEqual(resolve_qe_stage_nbands(40, self.state, ground_nbands=12), 40)
        self.write_state('<nbnd_up>24</nbnd_up><nbnd_dw>36</nbnd_dw>')
        self.assertEqual(resolve_qe_stage_nbands(None, self.state), 36)

    def test_missing_counts_keep_planned_or_native_defaults(self):
        self.write_state('')
        self.assertEqual(resolve_qe_stage_nbands(None, self.state, ground_nbands=22), 22)
        self.assertIsNone(resolve_qe_stage_nbands(None, self.state))
        self.assertEqual(resolve_qe_stage_nbands(None, self.root / 'missing', ground_nbands=22), 22)

    def test_bad_saved_metadata_is_rejected_without_overriding_explicit_counts(self):
        for text in ('0', '-2', '4.5', 'bad'):
            self.write_state(f'<nbnd>{text}</nbnd>')
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, 'positive integer'):
                resolve_qe_stage_nbands(None, self.state)
        self.schema.write_text('not XML', encoding='utf-8')
        self.assertEqual(resolve_qe_stage_nbands(48, self.state), 48)
        with self.assertRaisesRegex(ValueError, 'not valid XML'):
            resolve_qe_stage_nbands(None, self.state)

    def test_dry_run_uses_new_ground_count_or_saved_count_for_both_xcs(self):
        self.write_state()
        for xc in ('PBE', 'HSE06'):
            for ground in (True, False):
                with self.subTest(xc=xc, ground=ground):
                    config = self.config(XC_calc=xc, Ground_calc=ground, Ground_num_of_bands=24)
                    plan = prepare_qe_dry_run(config, self.struct)
                    expected = 24 if ground else 32
                    electronic = [job for job in plan['jobs'] if job['id'] in
                        ('dos-nscf', 'dos-hybrid-scf', 'band', 'band-hybrid-scf')]
                    self.assertEqual(len(electronic), 2)
                    for job in electronic:
                        self.assertIn(f'nbnd = {expected}', Path(job['input_file']).read_text())

    def test_explicit_stage_counts_override_shared_ground_setting(self):
        self.write_state()
        plan = prepare_qe_dry_run(self.config(Ground_calc=True, Ground_num_of_bands=24,
            DOS_num_of_bands=40, Band_num_of_bands=44), self.struct)
        for job in plan['jobs']:
            if job['id'] in ('dos-nscf', 'band'):
                expected = 40 if job['id'] == 'dos-nscf' else 44
                self.assertIn(f'nbnd = {expected}', Path(job['input_file']).read_text())

    def test_nscf_and_band_runners_write_the_inherited_allocation(self):
        self.write_state()
        path = build_band_path(self.atoms, 'GX', 5)
        for runner in (run_nscf, run_bands):
            for explicit in (None, 48):
                with self.subTest(runner=runner.__name__, explicit=explicit):
                    kwargs = dict(atoms=self.atoms, input_file=self.root / 'stage.in',
                        output_file=self.root / 'stage.out', state_dir=self.state,
                        pseudopotentials={'Si': 'Si.upf'}, pseudo_dir=self.root,
                        cutoff_ev=500, density_cutoff_ratio=4, nbands=explicit)
                    if runner is run_bands:
                        kwargs['band_path'] = path
                    with patch('nanoworks.engine.qe.run_qe_program') as run, patch(
                        'nanoworks.engine.qe.parse_pw_output', return_value={
                            'qe_version': (7, 4, 1), 'job_done': True}), patch(
                        'nanoworks.engine.qe.parse_pw_bands_output',
                        return_value={'nkpoints': path['npoints']}):
                        runner(**kwargs)
                    run.assert_called_once()
                    expected = 32 if explicit is None else explicit
                    self.assertIn(f'nbnd = {expected}', (self.root / 'stage.in').read_text())
