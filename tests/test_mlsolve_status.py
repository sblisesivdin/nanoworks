# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Check that ML workflow failures cannot be reported as successful runs."""

from unittest.mock import Mock

import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.filters import FrechetCellFilter

from nanoworks import mlsolve


def prepare_run(monkeypatch, tmp_path, **options):
    atoms = bulk('Cu')
    config = mlsolve.MLConfig(
        bulk_configuration=atoms, cell_relax=False,
        out_file=str(tmp_path / 'optimized.cif'), **options,
    )
    monkeypatch.setattr('sys.argv', ['mlsolve', '-g', 'Cu.cif', '-i', 'input.py'])
    monkeypatch.setattr(mlsolve, 'config_from_file', lambda **kwargs: ('Cu', config))
    calculator_factory = Mock(return_value=EMT())
    monkeypatch.setattr(mlsolve, 'get_ml_calculator', calculator_factory)
    writer = Mock()
    monkeypatch.setattr(mlsolve, 'write', writer)
    return config, calculator_factory, writer


@pytest.mark.parametrize('converged, status', [(True, 0), (False, 3)])
def test_optimization_status(monkeypatch, tmp_path, capsys, converged, status):
    config, _, writer = prepare_run(monkeypatch, tmp_path)
    dynamics = Mock()
    dynamics.run.return_value = converged
    monkeypatch.setattr(mlsolve, 'BFGS', Mock(return_value=dynamics))

    assert mlsolve.main() == status
    writer.assert_called_once_with(config.out_file, config.bulk_configuration)
    output = capsys.readouterr().out
    if not converged:
        assert 'did not converge' in output


def test_optimization_failure_saves_crash_in_output_directory(monkeypatch, tmp_path):
    config, _, writer = prepare_run(monkeypatch, tmp_path)
    dynamics = Mock()
    dynamics.run.side_effect = RuntimeError('calculator failed')
    monkeypatch.setattr(mlsolve, 'BFGS', Mock(return_value=dynamics))

    assert mlsolve.main() == 1
    writer.assert_called_once_with(str(tmp_path / 'crash_dump.cif'), config.bulk_configuration)


def test_cell_relaxation_passes_filter_to_optimizer_and_saves_atoms(monkeypatch, tmp_path):
    config, _, writer = prepare_run(monkeypatch, tmp_path)
    config.cell_relax = True
    dynamics = Mock()
    dynamics.run.return_value = True
    optimizer = Mock(return_value=dynamics)
    monkeypatch.setattr(mlsolve, 'BFGS', optimizer)

    assert mlsolve.main() == 0
    target = optimizer.call_args.args[0]
    assert isinstance(target, FrechetCellFilter)
    assert target.atoms is config.bulk_configuration
    # Exercise the combined atomic/cell force interface using EMT stress.
    assert target.get_forces().shape == (len(config.bulk_configuration) + 3, 3)
    writer.assert_called_once_with(config.out_file, config.bulk_configuration)


def test_static_failure_returns_error(monkeypatch, tmp_path):
    config, _, _ = prepare_run(monkeypatch, tmp_path, task='static')
    monkeypatch.setattr(
        config.bulk_configuration, 'get_potential_energy',
        Mock(side_effect=RuntimeError('calculator failed')),
    )

    assert mlsolve.main() == 1


@pytest.mark.parametrize('options', [{'task': 'md'}, {'optimizer': 'typo'}])
def test_invalid_selection_fails_before_model_loading(monkeypatch, tmp_path, options):
    _, calculator_factory, _ = prepare_run(monkeypatch, tmp_path, **options)

    with pytest.raises(SystemExit) as error:
        mlsolve.main()
    assert error.value.code == 2
    calculator_factory.assert_not_called()
