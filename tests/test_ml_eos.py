# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""EOS validation, units, persistence, and isolation using an analytic potential."""

import json
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.eos import birchmurnaghan
from ase.io import read
from ase.units import GPa

from nanoworks import mlsolve
from nanoworks.ml_eos import fit_eos, run_eos, validate_eos


class AnalyticEOS(Calculator):
    implemented_properties = ['energy']

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results = {'energy': birchmurnaghan(atoms.get_volume(), -3.0, 0.7, 4.1, 10.0)}


def setup_eos(tmp_path, **options):
    atoms = Atoms('Cu', scaled_positions=[[0.1, 0.2, 0.3]], cell=np.eye(3) * 10 ** (1 / 3), pbc=True)
    atoms.calc = AnalyticEOS()
    config = mlsolve.MLConfig(task='eos', out_file=str(tmp_path / 'optimized.cif'),
                             eos_relax_atoms=False, **options)
    return atoms, config


def test_eos_recovers_known_bulk_properties_and_preserves_input(tmp_path):
    atoms, config = setup_eos(tmp_path)
    cell, positions = atoms.cell.copy(), atoms.positions.copy()

    assert run_eos(atoms, config, 'Cu') == 0
    report = json.loads((tmp_path / 'Cu-ML-EOS-Fit.json').read_text())
    assert report['status'] == 'success'
    assert report['result']['volume_A3'] == pytest.approx(10.0, rel=1e-5)
    assert report['result']['bulk_modulus_GPa'] == pytest.approx(0.7 / GPa, rel=1e-5)
    assert report['result']['bulk_modulus_pressure_derivative'] == pytest.approx(4.1, rel=1e-4)
    assert report['result']['fit_rmse_eV'] < 1e-8
    data = np.loadtxt(tmp_path / 'Cu-ML-EOS-Result.dat')
    assert data.shape == (11, 5)
    assert data[:, 1] == pytest.approx(data[:, 0] * 10)
    assert (tmp_path / 'Cu-ML-EOS-Graph.png').stat().st_size > 0
    structures = read(tmp_path / 'Cu-ML-EOS-Structures.traj', index=':')
    assert len(structures) == 11
    assert structures[0].get_potential_energy() == pytest.approx(data[0, 2])
    np.testing.assert_allclose(atoms.cell, cell)
    np.testing.assert_allclose(atoms.positions, positions)


@pytest.mark.parametrize('options', [
    {'eos_scale': (1.1, 0.9)}, {'eos_scale': (0, 1)}, {'eos_scale': (1, float('nan'))},
    {'eos_points': 4}, {'eos_points': 5.5}, {'eos_fit': 'typo'},
])
def test_invalid_eos_settings(tmp_path, options):
    atoms, config = setup_eos(tmp_path, **options)
    with pytest.raises(ValueError):
        validate_eos(atoms, config)


def test_eos_rejects_non_bulk(tmp_path):
    atoms, config = setup_eos(tmp_path)
    atoms.pbc = (True, True, False)
    with pytest.raises(ValueError, match='3D bulk'):
        validate_eos(atoms, config)


def test_unconverged_point_preserves_samples_without_fit(tmp_path):
    atoms, config = setup_eos(tmp_path, eos_points=5)
    config.eos_relax_atoms = True
    dynamics = Mock()
    dynamics.run.side_effect = [True, True, False, True, True]
    optimizer = Mock(return_value=dynamics)

    assert run_eos(atoms, config, 'Cu', optimizer) == 3
    report = json.loads((tmp_path / 'Cu-ML-EOS-Fit.json').read_text())
    assert report['status'] == 'not_converged'
    assert len(report['samples']) == 5
    assert 'result' not in report
    assert not (tmp_path / 'Cu-ML-EOS-Graph.png').exists()
    for call in optimizer.call_args_list:
        assert isinstance(call.args[0], Atoms)


def test_calculator_failure_keeps_completed_samples(tmp_path):
    atoms, config = setup_eos(tmp_path)
    original = atoms.calc.calculate
    calls = 0

    def fail_at_third(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise RuntimeError('model failed')
        original(*args, **kwargs)

    atoms.calc.calculate = fail_at_third
    assert run_eos(atoms, config, 'Cu') == 1
    report = json.loads((tmp_path / 'Cu-ML-EOS-Fit.json').read_text())
    assert report['status'] == 'failed'
    assert len(report['samples']) == 2
    assert 'result' not in report
    assert np.loadtxt(tmp_path / 'Cu-ML-EOS-Result.dat').shape == (2, 5)


def test_cli_rejects_invalid_eos_before_loading_model(tmp_path, monkeypatch):
    atoms, config = setup_eos(tmp_path, eos_points=2)
    config.bulk_configuration = atoms
    monkeypatch.setattr('sys.argv', ['mlsolve', '-g', 'Cu.cif', '-i', 'input.py'])
    monkeypatch.setattr(mlsolve, 'config_from_file', lambda **kwargs: ('Cu', config))
    factory = Mock()
    monkeypatch.setattr(mlsolve, 'get_ml_calculator', factory)
    with pytest.raises(SystemExit) as error:
        mlsolve.main()
    assert error.value.code == 2
    factory.assert_not_called()


def test_fit_outside_scan_is_rejected(monkeypatch):
    eos = Mock()
    eos.fit.return_value = (15.0, -3.0, 0.7)
    eos.eos_parameters = [-3.0, 0.7, 4.1, 15.0]
    monkeypatch.setattr('nanoworks.ml_eos.EquationOfState', Mock(return_value=eos))
    with pytest.raises(ValueError, match='outside'):
        fit_eos(np.linspace(9, 11, 5), np.zeros(5), 'birchmurnaghan')


def test_unconverged_rerun_does_not_keep_previous_fit_or_graph(tmp_path):
    atoms, config = setup_eos(tmp_path, eos_points=5)
    assert run_eos(atoms, config, 'Cu') == 0
    config.eos_relax_atoms = True
    optimizer = Mock(return_value=Mock(run=Mock(return_value=False)))

    assert run_eos(atoms, config, 'Cu', optimizer) == 3
    report = json.loads((tmp_path / 'Cu-ML-EOS-Fit.json').read_text())
    assert report['status'] == 'not_converged'
    assert 'result' not in report
    assert not (tmp_path / 'Cu-ML-EOS-Graph.png').exists()
