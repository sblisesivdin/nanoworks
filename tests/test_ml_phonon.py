# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Analytic force-constant and failure-path checks without ML model downloads."""

import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import _amu, _e

from nanoworks.ml_phonon import EV_TO_THZ, run_phonon, validate_phonon
from nanoworks.mlsolve import MLConfig


class PairSpring(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, spring=2.0, **kwargs):
        super().__init__(**kwargs)
        self.spring = spring

    def calculate(self, atoms=None, properties=('energy', 'forces'), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        displacement = atoms.positions[1] - atoms.positions[0] - [1, 0, 0]
        force = self.spring * displacement
        self.results = {'energy': self.spring * displacement @ displacement / 2,
                        'forces': np.array([force, -force])}


def setup(tmp_path, spring=2.0):
    atoms = Atoms('H2', positions=[[0, 0, 0], [1, 0, 0]], cell=[4, 4, 4], pbc=True)
    atoms.set_masses([1, 1])
    atoms.calc = PairSpring(spring)
    config = MLConfig(task='phonon', phonon_supercell=(1, 1, 1), phonon_mesh=(2, 2, 2),
                      phonon_npoints=10, phonon_dos_bins=20, phonon_acoustic=False,
                      out_file=str(tmp_path / 'optimized.cif'))
    return atoms, config


def report(tmp_path):
    return json.loads((tmp_path / 'pair-ML-PHONON-Result.json').read_text())


def test_analytic_optical_modes_units_dos_and_geometry(tmp_path):
    atoms, config = setup(tmp_path)
    positions = atoms.positions.copy()
    assert run_phonon(atoms, config, 'pair') == 0
    result = report(tmp_path)
    # Two equal masses coupled by k: optical angular frequency sqrt(2k/m).
    expected = np.sqrt(4 * _e / _amu) * 1e10 / (2 * np.pi * 1e12)
    np.testing.assert_allclose(result['gamma_frequencies_THz'][:3], 0, atol=1e-6)
    np.testing.assert_allclose(result['gamma_frequencies_THz'][3:], expected, rtol=1e-8)
    assert EV_TO_THZ == pytest.approx(241.7989, rel=1e-5)
    assert result['dos_integrated_modes'] == pytest.approx(6)
    assert result['has_sampled_imaginary_modes'] is False
    assert result['status'] == 'success'
    np.testing.assert_array_equal(atoms.positions, positions)
    assert (tmp_path / 'pair-ML-PHONON-Graph.png').is_file()
    assert list(Path(result['force_cache']).glob('cache.*.json'))


def test_unstable_modes_preserve_negative_sign(tmp_path):
    atoms, config = setup(tmp_path, spring=-2)
    assert run_phonon(atoms, config, 'pair') == 0
    result = report(tmp_path)
    assert result['has_sampled_imaginary_modes'] is True
    assert result['minimum_frequency_THz'] < -1
    assert result['imaginary_modes_on_mesh'] == 3 * 8


def test_unrelaxed_reference_removes_old_results(tmp_path):
    atoms, config = setup(tmp_path)
    assert run_phonon(atoms, config, 'pair') == 0
    old_cache = report(tmp_path)['force_cache']
    atoms.positions[1, 0] += 0.2
    assert run_phonon(atoms, config, 'pair') == 3
    assert report(tmp_path)['status'] == 'not_converged'
    assert report(tmp_path)['force_cache'] != old_cache
    assert not (tmp_path / 'pair-ML-PHONON-Graph.png').exists()
    assert not (tmp_path / 'pair-ML-PHONON-Bands.dat').exists()


def test_nonfinite_displaced_forces_retain_cache(tmp_path):
    class InvalidSpring(PairSpring):
        def calculate(self, *args, **kwargs):
            super().calculate(*args, **kwargs)
            if np.linalg.norm(self.results['forces']) > 0:
                self.results['forces'][:] = np.nan
    atoms, config = setup(tmp_path)
    atoms.calc = InvalidSpring()
    assert run_phonon(atoms, config, 'pair') == 1
    assert 'Nonfinite' in report(tmp_path)['message']
    assert not (tmp_path / 'pair-ML-PHONON-ForceConstants.npz').exists()


@pytest.mark.parametrize('name,value', [('phonon_supercell', (2, True, 2)),
    ('phonon_mesh', (2, 0, 2)), ('phonon_delta', 0), ('phonon_npoints', 1),
    ('phonon_dos_bins', 1), ('phonon_acoustic', 'yes'),
    ('phonon_imaginary_tolerance', -1), ('phonon_path', '')])
def test_invalid_settings(tmp_path, name, value):
    atoms, config = setup(tmp_path)
    setattr(config, name, value)
    with pytest.raises(ValueError):
        validate_phonon(atoms, config)


def test_nonperiodic_and_constrained_inputs_rejected(tmp_path):
    from ase.constraints import FixAtoms
    atoms, config = setup(tmp_path)
    atoms.pbc[2] = False
    with pytest.raises(ValueError, match='3D bulk'):
        validate_phonon(atoms, config)
    atoms.pbc = True
    atoms.set_constraint(FixAtoms(indices=[0]))
    with pytest.raises(ValueError, match='unconstrained'):
        validate_phonon(atoms, config)


def test_cli_validation_before_model_loading(tmp_path, monkeypatch):
    from nanoworks import mlsolve
    atoms, config = setup(tmp_path)
    config.phonon_delta = -1
    monkeypatch.setattr(mlsolve, 'config_from_file', lambda **kwargs: ('pair', config))
    config.bulk_configuration = atoms
    monkeypatch.setattr(mlsolve, 'get_ml_calculator', lambda *args, **kwargs: pytest.fail('Model loaded'))
    monkeypatch.setattr('sys.argv', ['mlsolve', '-g', 'pair.cif', '-i', 'input.py'])
    with pytest.raises(SystemExit) as exc:
        mlsolve.main()
    assert exc.value.code == 2
