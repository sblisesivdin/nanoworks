"""Input paths and current contents must determine the ML configuration."""

import sys

import pytest

from nanoworks.mlsolve import config_from_file


def test_model_options_from_input_file(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / "ml_model_options.py"
    input_file.write_text(
        "variant = 'small'\n"
        "dtype = 'float32'\n"
        "organic = True\n"
        "dispersion = True\n"
        "model_path = 'custom.pt'\n"
        "model_name = '7net-custom'\n"
    )

    _, config = config_from_file(input_file, geometryfile=None)

    assert config.variant == 'small'
    assert config.dtype == 'float32'
    assert config.organic is True
    assert config.dispersion is True
    assert config.model_path == 'custom.pt'
    assert config.model_name == '7net-custom'


def test_same_named_inputs_in_different_directories(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for dirname, model in [('first', 'mace'), ('second', 'chgnet')]:
        directory = tmp_path / dirname
        directory.mkdir()
        input_file = directory / 'input.py'
        input_file.write_text(f"model = '{model}'\n")

        _, config = config_from_file(input_file, geometryfile=None)
        assert config.model == model


def test_input_is_reloaded_after_edit(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / 'input.py'
    input_file.write_text("fmax = 0.05\n")
    _, first = config_from_file(input_file, geometryfile=None)
    input_file.write_text("fmax = 0.01\n")
    _, second = config_from_file(input_file, geometryfile=None)

    assert first.fmax == 0.05
    assert second.fmax == 0.01


def test_input_filename_and_sibling_import(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    directory = tmp_path / 'inputs'
    directory.mkdir()
    (directory / 'ml_test_settings.py').write_text("model = 'sevennet'\n")
    input_file = directory / 'ml-input.v1.py'
    input_file.write_text("from ml_test_settings import model\n")
    original_path = sys.path[:]
    try:
        _, config = config_from_file(input_file, geometryfile=None)
        assert config.model == 'sevennet'
        assert sys.path == original_path
    finally:
        sys.modules.pop('ml_test_settings', None)


def test_import_path_restored_on_input_error(tmp_path):
    input_file = tmp_path / 'broken.py'
    input_file.write_text("raise RuntimeError('broken input')\n")
    original_path = sys.path[:]

    with pytest.raises(RuntimeError, match='broken input'):
        config_from_file(input_file, geometryfile=None)
    assert sys.path == original_path
