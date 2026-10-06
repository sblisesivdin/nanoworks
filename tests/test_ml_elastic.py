# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Known stress responses verify shear convention, slab units, and workflow status."""

import json
from unittest.mock import Mock

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.units import GPa

from nanoworks.ml_elastic import run_elastic, validate_elastic
from nanoworks.mlsolve import MLConfig


class LinearStress(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, reference_cell, stiffness_gpa, prestress=None):
        super().__init__()
        self.reference_cell = np.array(reference_cell)
        self.stiffness = np.array(stiffness_gpa)
        self.prestress = np.zeros(6) if prestress is None else np.array(prestress)

    def calculate(self, atoms=None, properties=('energy',), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        deformation = atoms.cell.T @ np.linalg.inv(self.reference_cell.T)
        strain = (deformation + deformation.T) / 2 - np.eye(3)
        voigt = np.array([strain[0, 0], strain[1, 1], strain[2, 2],
                          2 * strain[1, 2], 2 * strain[0, 2], 2 * strain[0, 1]])
        self.results = {'energy': float(voigt @ self.stiffness @ voigt) / 2,
                        'forces': np.zeros((len(atoms), 3)),
                        'stress': (self.stiffness @ voigt + self.prestress) * GPa}


def setup_elastic(tmp_path, cell=None, tensor=None, **options):
    cell = np.diag([3., 4., 5.]) if cell is None else np.array(cell)
    if tensor is None:
        tensor = np.zeros((6, 6))
        tensor[:3, :3] = 80
        np.fill_diagonal(tensor[:3, :3], 200)
        tensor[3:, 3:] = np.eye(3) * 60
    atoms = Atoms('Cu', scaled_positions=[[0, 0, 0]], cell=cell, pbc=True)
    atoms.calc = LinearStress(cell, tensor)
    config = MLConfig(task='elastic', elastic_relax_internal=False,
                      out_file=str(tmp_path / 'optimized.cif'), **options)
    return atoms, config, tensor


def read_result(tmp_path):
    return json.loads((tmp_path / 'Cu-ML-ELASTIC-Result.json').read_text())


def test_3d_recovers_tensor_and_moduli_without_mutating_input(tmp_path):
    atoms, config, tensor = setup_elastic(tmp_path, elastic_dimensionality='3D')
    cell, positions = atoms.cell.copy(), atoms.positions.copy()
    assert run_elastic(atoms, config, 'Cu') == 0
    report = read_result(tmp_path)
    assert report['status'] == 'success'
    assert report['units'] == 'GPa'
    np.testing.assert_allclose(report['result']['raw_stiffness'], tensor, atol=1e-9)
    hill = report['result']['properties']['hill']
    assert hill['bulk_modulus_gpa'] == pytest.approx(120)
    assert hill['shear_modulus_gpa'] == pytest.approx(60)
    assert hill['young_modulus_gpa'] == pytest.approx(9 * 120 * 60 / 420)
    assert hill['poisson_ratio'] == pytest.approx(240 / 840)
    assert report['result']['stability']['mechanically_stable'] is True
    assert len(report['samples']) == 25
    np.testing.assert_allclose(atoms.cell, cell)
    np.testing.assert_allclose(atoms.positions, positions)


@pytest.mark.parametrize('normal', ['x', 'y', 'z'])
def test_2d_stiffness_is_independent_of_vacuum_height(tmp_path, normal):
    normal_index = ('x', 'y', 'z').index(normal)
    indices = {'x': (1, 2, 3), 'y': (0, 2, 4), 'z': (0, 1, 5)}[normal]
    target = np.array([[160., 50., 0.], [50., 160., 0.], [0., 0., 55.]])
    for height in (15., 30.):
        lengths = [3., 4., 5.]
        lengths[normal_index] = height
        tensor = np.zeros((6, 6))
        tensor[np.ix_(indices, indices)] = target / (height * 0.1)
        atoms, config, _ = setup_elastic(tmp_path, cell=np.diag(lengths), tensor=tensor,
                                          elastic_dimensionality='2D', elastic_normal_axis=normal)
        assert run_elastic(atoms, config, 'Cu') == 0
        report = read_result(tmp_path)
        assert report['units'] == 'N/m'
        np.testing.assert_allclose(report['result']['symmetric_stiffness'], target, atol=1e-9)
        assert 'full_tensor_gpa' not in report['result']
        assert report['result']['properties']['poisson_ratio_first_second'] == pytest.approx(50 / 160)


def test_reference_stress_does_not_certify_zero_stress_stability(tmp_path):
    atoms, config, _ = setup_elastic(tmp_path, elastic_dimensionality='3D')
    atoms.calc.prestress = np.array([2., 0., 0., 0., 0., 0.])
    assert run_elastic(atoms, config, 'Cu') == 0
    stability = read_result(tmp_path)['result']['stability']
    assert stability['positive_definite'] is True
    assert stability['zero_prestress_condition_met'] is False
    assert stability['mechanically_stable'] is None


def test_unconverged_reference_keeps_samples_without_tensor(tmp_path):
    atoms, config, _ = setup_elastic(tmp_path)
    config.elastic_relax_internal = True
    dynamics = Mock(run=Mock(return_value=False))
    assert run_elastic(atoms, config, 'Cu', Mock(return_value=dynamics)) == 3
    report = read_result(tmp_path)
    assert report['status'] == 'not_converged'
    assert len(report['samples']) == 1
    assert 'result' not in report
    assert not (tmp_path / 'Cu-ML-ELASTIC-Tensor.dat').exists()


def test_missing_stress_is_a_calculation_error(tmp_path):
    atoms, config, _ = setup_elastic(tmp_path)
    atoms.calc.implemented_properties = ['energy', 'forces']
    assert run_elastic(atoms, config, 'Cu') == 1
    assert read_result(tmp_path)['status'] == 'failed'


def test_singular_tensor_is_saved_without_invalid_moduli(tmp_path):
    atoms, config, _ = setup_elastic(tmp_path, tensor=np.zeros((6, 6)),
                                      elastic_dimensionality='3D')
    assert run_elastic(atoms, config, 'Cu') == 0
    result = read_result(tmp_path)['result']
    assert 'properties' not in result
    assert result['stability']['mechanically_stable'] is False
    assert result['stability']['condition_number'] is None
    assert (tmp_path / 'Cu-ML-ELASTIC-Tensor.dat').exists()


def test_cli_rejects_bad_strain_before_model_loading(tmp_path, monkeypatch):
    from nanoworks import mlsolve
    atoms, config, _ = setup_elastic(tmp_path, elastic_strain=0)
    config.bulk_configuration = atoms
    monkeypatch.setattr('sys.argv', ['mlsolve', '-g', 'Cu.cif', '-i', 'input.py'])
    monkeypatch.setattr(mlsolve, 'config_from_file', lambda **kwargs: ('Cu', config))
    factory = Mock()
    monkeypatch.setattr(mlsolve, 'get_ml_calculator', factory)
    with pytest.raises(SystemExit) as error:
        mlsolve.main()
    assert error.value.code == 2
    factory.assert_not_called()


@pytest.mark.parametrize('options', [
    {'elastic_strain': 0}, {'elastic_strain': float('nan')},
    {'elastic_points': 4}, {'elastic_points': 3.5}, {'elastic_dimensionality': '1D'},
])
def test_bad_elastic_inputs_are_rejected(tmp_path, options):
    atoms, config, _ = setup_elastic(tmp_path, **options)
    with pytest.raises(ValueError):
        validate_elastic(atoms, config)


def test_tilted_slab_normal_is_rejected(tmp_path):
    atoms, config, _ = setup_elastic(tmp_path, cell=[[3, 0, 0], [0, 4, 0], [1, 0, 20]],
                                     elastic_dimensionality='2D')
    with pytest.raises(ValueError, match='reorient'):
        validate_elastic(atoms, config)
