# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from ase.build import bulk

from nanoworks.engine.qe import (
    EV_PER_RYDBERG, Bohr, run_relax, validate_relaxed_final_forces,
    validate_relaxed_final_stress,
)

DATA = Path(__file__).parent / 'data'


def test_recorded_final_rechecks_distinguish_valid_and_failed_geometry():
    reference = (DATA / 'qe_zno_reference_final_recheck.out').read_text()
    bad = (DATA / 'qe_zno_hubbard_final_recheck.out').read_text()
    assert 'bfgs converged' in bad and 'JOB DONE.' in bad
    assert validate_relaxed_final_forces(reference, 4, .05) < .05
    with pytest.raises(RuntimeError, match='despite the earlier BFGS'):
        validate_relaxed_final_forces(bad, 4, .05)


@pytest.mark.parametrize('case', ['reference', 'hubbard'])
def test_run_relax_checks_recorded_final_forces_before_returning_structure(tmp_path, case):
    text = (DATA / f'qe_zno_{case}_final_recheck.out').read_text()

    def fake_program(**kwargs):
        Path(kwargs['output_file']).write_text(text)
        return {'returncode': 0}

    with patch('nanoworks.engine.qe.run_qe_program', side_effect=fake_program):
        kwargs = dict(atoms=bulk('ZnO', 'wurtzite', a=3.25, c=5.2),
                      input_file=tmp_path/'relax.in', output_file=tmp_path/'relax.out',
                      state_dir=tmp_path/'state', pseudopotentials={'Zn':'Zn.UPF','O':'O.UPF'},
                      pseudo_dir=tmp_path, cutoff_ev=1400, optimizer='LBFGS',
                      max_force=.05, max_step=.1, relax_cell=[True]*6)
        if case == 'hubbard':
            with pytest.raises(RuntimeError, match='final geometry force component'):
                run_relax(**kwargs)
        else:
            result = run_relax(**kwargs)
            assert result['geometry_converged']
            assert result['final_max_force_component_ev_angstrom'] < .05
            assert result['atoms'].get_volume() == pytest.approx(49.25725, abs=1e-4)


def test_qe_force_component_convention_not_vector_norm():
    value = .04 * Bohr / EV_PER_RYDBERG
    text = f'''convergence has been achieved
Forces acting on atoms (cartesian axes, Ry/au):
atom 1 type 1 force = {value:.12e} {value:.12e} 0.0
'''
    assert np.sqrt(2)*.04 > .05
    assert validate_relaxed_final_forces(text, 1, .05) == pytest.approx(.04)


@pytest.mark.parametrize('ending', [
    'convergence NOT achieved',
    'convergence has been achieved',
    'Forces acting on atoms (cartesian axes, Ry/au):\n',
    'Forces acting on atoms (cartesian axes, eV/Angstrom):\natom 1 type 1 force = 0 0 0',
])
def test_late_invalid_recheck_cannot_fall_back_to_earlier_good_forces(ending):
    text = '''convergence has been achieved
Forces acting on atoms (cartesian axes, Ry/au):
atom 1 type 1 force = 0 0 0
bfgs converged in 5 scf cycles and 4 bfgs steps
''' + ending
    with pytest.raises(RuntimeError):
        validate_relaxed_final_forces(text, 1, .05)


def test_qe_vc_relax_stress_recheck_accepts_recorded_valid_zno():
    reference = (DATA / 'qe_zno_reference_final_recheck.out').read_text()
    # Only the final SCF stress is accepted, not an earlier BFGS tensor.
    assert validate_relaxed_final_stress(reference) == pytest.approx(.34)


def test_qe_mos2_d3_final_stress_and_forces_reject_earlier_bfgs_convergence():
    text = (DATA / 'qe_mos2_d3_final_recheck.out').read_text()
    assert 'bfgs converged' in text and 'JOB DONE.' in text
    with pytest.raises(RuntimeError, match='despite the earlier BFGS'):
        validate_relaxed_final_forces(text, 6, .05)
    with pytest.raises(RuntimeError, match='283.1 kbar'):
        validate_relaxed_final_stress(text)


def test_qe_vc_relax_stress_recheck_detects_false_positive_force_convergence():
    # A 0.08 eV/Ang force target would accept the archived MoS2 forces,
    # but the huge final cell stress must still reject the relaxed cell.
    text = (DATA / 'qe_mos2_d3_final_recheck.out').read_text()
    assert validate_relaxed_final_forces(text, 6, .08) < .08
    with pytest.raises(RuntimeError, match='exceeds 0.5 kbar'):
        validate_relaxed_final_stress(text)


@pytest.mark.parametrize('ending', [
    'convergence NOT achieved',
    'convergence has been achieved',
    'convergence has been achieved\n'
    'total stress (Ry/bohr**3) (kbar) P=nan',
])
def test_qe_vc_relax_stress_recheck_requires_completed_last_block(ending):
    text = ('convergence has been achieved\n'
            'total stress (Ry/bohr**3) (kbar) P=0.01\n'
            '0 0 0  0 0 0\n'
            '0 0 0  0 0 0\n'
            '0 0 0  0 0 0\n' + ending)
    with pytest.raises(RuntimeError):
        validate_relaxed_final_stress(text)


def test_qe_vc_relax_stress_recheck_supports_hydrostatic_target():
    text = ('convergence has been achieved\n'
            'total stress (Ry/bohr**3) (kbar) P=10\n'
            '0 0 0  -10 0 0\n'
            '0 0 0  0 -10 0\n'
            '0 0 0  0 0 -10\n')
    assert validate_relaxed_final_stress(text, target_pressure_kbar=10) == 0.0
    with pytest.raises(RuntimeError, match='exceeds'):
        validate_relaxed_final_stress(text, target_pressure_kbar=0)
