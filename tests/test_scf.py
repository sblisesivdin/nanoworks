# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Tests for engine-neutral SCF controls."""

import pytest

from nanoworks.scf import (
    qe_conv_thr_to_accuracy,
    resolve_gpaw_scf_settings,
    resolve_qe_scf_settings,
    validate_scf_settings,
)


def test_qe_adapter_translates_portable_settings():
    assert resolve_qe_scf_settings(
        accuracy='tight',
        max_steps=220,
        mixing=0.25,
        solver='robust',
    ) == {
        'conv_thr': 1.0e-8,
        'mixing_beta': 0.25,
        'electron_maxstep': 220,
        'diagonalization': 'cg',
    }


def test_gpaw_adapter_translates_portable_settings():
    settings = resolve_gpaw_scf_settings(
        accuracy='tight',
        max_steps=220,
        mixing=0.25,
        solver='fast',
    )

    assert settings['convergence']['energy'] == 1.0e-6
    assert settings['mixing'] == 0.25
    assert settings['maxiter'] == 220
    assert settings['eigensolver'] == 'rmm-diis'


@pytest.mark.parametrize(
    ('conv_thr', 'accuracy'),
    [
        (1.0e-4, 'loose'),
        (1.0e-6, 'normal'),
        (1.0e-8, 'tight'),
        (1.0e-10, 'very-tight'),
    ],
)
def test_qe_threshold_maps_to_portable_profile(conv_thr, accuracy):
    assert qe_conv_thr_to_accuracy(conv_thr) == accuracy


@pytest.mark.parametrize('value', [True, False, 0, -1, 2.9, 10.0, '10', float('nan'), float('inf')])
def test_scf_steps_must_be_positive_integers_without_truncation(value):
    with pytest.raises(ValueError, match='SCF_max_steps'):
        validate_scf_settings(max_steps=value)


@pytest.mark.parametrize('value', [True, False, None, 0, -0.1, 1.1, float('nan'), float('inf'), 'invalid'])
def test_scf_mixing_requires_finite_non_boolean_fraction(value):
    if value is None:
        assert validate_scf_settings(mixing=None)['mixing'] is None
    else:
        with pytest.raises(ValueError, match='SCF_mixing'):
            validate_scf_settings(mixing=value)


@pytest.mark.parametrize('value', [0, -1e-8, True, False, float('nan'), float('inf'), '-inf', 'invalid', None])
def test_invalid_qe_convergence_threshold_is_not_mapped_to_loose(value):
    with pytest.raises(ValueError, match='conv_thr'):
        qe_conv_thr_to_accuracy(value)


def test_valid_integer_scf_steps_and_qe_thresholds_still_map():
    assert validate_scf_settings(max_steps=150, mixing=0.3)['max_steps'] == 150
    assert validate_scf_settings(max_steps=150, mixing=0.3)['mixing'] == 0.3
    assert qe_conv_thr_to_accuracy(1e-6) == 'normal'


@pytest.mark.parametrize(
    ('keyword', 'value'),
    [
        ('accuracy', 'unknown'),
        ('max_steps', 0),
        ('mixing', 1.1),
        ('solver', 'david'),
    ],
)
def test_invalid_portable_scf_setting_is_rejected(keyword, value):
    kwargs = {keyword: value}

    with pytest.raises(ValueError):
        validate_scf_settings(**kwargs)
