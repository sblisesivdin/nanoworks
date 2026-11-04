# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from types import SimpleNamespace

import numpy as np
import pytest

from nanoworks.dos import ensure_gpaw_soc_dos_weights


@pytest.mark.parametrize('local_points', [0, 1, 2, 4, 9])
def test_full_bz_states_survive_rank_local_weight_lengths(local_points):
    # Every point has different eigenvalues: scaling a partial sum would
    # still leave the wrong spectrum even if its integral were repaired.
    energies = np.linspace(-8, 8, 16001)
    eigenvalues = np.arange(18, dtype=float).reshape(1, 9, 2) / 3 - 3
    rawdos = SimpleNamespace(eig_skn=eigenvalues,
                             weight_k=np.full(local_points, 1 / 9))
    changed = ensure_gpaw_soc_dos_weights(rawdos, 9)
    assert changed == (local_points != 9)
    width = .08
    actual = sum(w * np.exp(-((energies[:, None] - eig) / width)**2).sum(axis=1)
                 / (np.sqrt(np.pi) * width)
                 for eig, w in zip(rawdos.eig_skn[0], rawdos.weight_k))
    expected = np.exp(-((energies[:, None] - eigenvalues.ravel()) / width)**2).sum(axis=1)
    expected /= 9 * np.sqrt(np.pi) * width
    np.testing.assert_allclose(actual, expected, atol=1e-14)
    assert np.trapezoid(actual, energies) == pytest.approx(2.0, abs=1e-12)
    # Repeated calls and an already-correct serial/full array are unchanged.
    weights = rawdos.weight_k
    assert not ensure_gpaw_soc_dos_weights(rawdos, 9)
    assert rawdos.weight_k is weights


@pytest.mark.parametrize('weights', [[np.nan], [-1 / 3], [1 / 2, 1 / 2], [1 / 3] * 4])
def test_unknown_weight_layouts_are_not_silently_rescaled(weights):
    rawdos = SimpleNamespace(eig_skn=np.zeros((1, 3, 2)), weight_k=np.array(weights))
    with pytest.raises(ValueError, match='unexpected full-BZ'):
        ensure_gpaw_soc_dos_weights(rawdos, 3)


@pytest.mark.parametrize('eigenvalues', [np.zeros((1, 2, 4)), np.zeros((2, 3, 4)),
                                       np.full((1, 3, 4), np.nan)])
def test_missing_or_nonfinite_global_eigenvalues_are_rejected(eigenvalues):
    rawdos = SimpleNamespace(eig_skn=eigenvalues, weight_k=np.full(3, 1 / 3))
    with pytest.raises(ValueError, match='eigenvalues for the full BZ'):
        ensure_gpaw_soc_dos_weights(rawdos, 3)
