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
