"""Input-file model options must reach the ML calculator configuration."""

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
