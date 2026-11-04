# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Cartesian-frame elasticity regressions without a DFT engine."""
import io
from pathlib import Path
import zipfile

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read
from ase.stress import voigt_6_to_full_3x3_stress
from ase.units import GPa

from nanoworks.elasticity import fit_cartesian_elastic_tensor


PAIRS = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def strain_tensor(e):
    result = voigt_6_to_full_3x3_stress(e).copy()
    for i, j in PAIRS[3:]:
        result[i, j] /= 2
        result[j, i] /= 2
    return result


def rotate_stiffness(c, q):
    """Independent fourth-rank rotation; no engineering-shear shortcuts."""
    indices = {(i, j): a for a, (i, j) in enumerate(PAIRS)}
    indices.update({(j, i): a for a, (i, j) in enumerate(PAIRS)})
    t = np.empty((3, 3, 3, 3))
    for i, j, k, l in np.ndindex(3, 3, 3, 3):
        t[i, j, k, l] = c[indices[i, j], indices[k, l]]
    rotated = np.einsum('ia,jb,kc,ld,abcd->ijkl', q, q, q, q, t)
    return np.array([[rotated[i, j, k, l] for k, l in PAIRS] for i, j in PAIRS])


def samples(response, reference_stress=None):
    ref = Atoms('TiC', scaled_positions=[[0, 0, 0], [.5, .5, .5]],
                cell=[[3, 0, 0], [1.5, 2.6, 0], [1.5, .9, 2.5]], pbc=True)
    stress0 = np.zeros(6) if reference_stress is None else reference_stress
    ref.calc = SinglePointCalculator(ref, energy=0, stress=stress0)
    frames = []
    for mode in range(6):
        # Unequal/fallback amplitudes and a nonzero intercept must work.
        for value in [-.008, -.002, 0, .001, .006]:
            e = np.zeros(6)
            e[mode] = value
            atom = ref.copy()
            atom.set_cell(ref.cell @ (np.eye(3) + strain_tensor(e)), scale_atoms=True)
            atom.calc = SinglePointCalculator(atom, energy=0, stress=stress0 + response @ e)
            frames.append(atom)
    return ref, frames


@pytest.mark.parametrize('crystal', ['cubic', 'hexagonal', 'triclinic'])
def test_rotated_crystals_retain_full_cartesian_tensor(crystal):
    if crystal == 'cubic':
        c = np.diag([520., 520., 520., 170., 170., 170.])
        c[:3, :3] += 120 * (np.ones((3, 3)) - np.eye(3))
    elif crystal == 'hexagonal':
        c = np.diag([200., 200., 300., 55., 55., 70.])
        c[0, 1] = c[1, 0] = 60
        c[0, 2] = c[2, 0] = c[1, 2] = c[2, 1] = 40
    else:
        a = np.random.default_rng(19).normal(size=(6, 6))
        c = 10 * a @ a.T + 50 * np.eye(6)
    q, _ = np.linalg.qr(np.random.default_rng(23).normal(size=(3, 3)))
    expected = rotate_stiffness(c, q) * GPa
    assert np.max(abs(expected[:3, 3:])) > .01 * GPa
    ref, frames = samples(expected)
    fitted, info = fit_cartesian_elastic_tensor(ref, frames)
    np.testing.assert_allclose(fitted, expected, atol=1e-11, rtol=0)
    assert info['rank'] == 6
    assert info['rms_residual'] < 1e-13


def test_rotating_entire_trajectory_rotates_recovered_response():
    c = np.diag([500., 510., 520., 100., 110., 120.]) * GPa
    ref, frames = samples(c)
    angle = .71
    q = np.array([[np.cos(angle), -np.sin(angle), 0],
                  [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    rotated = []
    for atom in [ref, *frames]:
        stress = q @ voigt_6_to_full_3x3_stress(atom.get_stress()) @ q.T
        new = atom.copy()
        new.set_cell(atom.cell @ q.T, scale_atoms=True)
        new.calc = SinglePointCalculator(new, energy=0, stress=stress)
        rotated.append(new)
    fitted, _ = fit_cartesian_elastic_tensor(rotated[0], rotated[1:])
    np.testing.assert_allclose(fitted, rotate_stiffness(c, q), atol=1e-11, rtol=0)


def test_reference_subtraction_symmetrization_and_pressure_correction():
    response = np.diag([4., 5., 6., 1., 2., 3.])
    response[0, 3], response[3, 0] = .3, .1
    stress0 = np.array([-.02, -.02, -.02, 0, 0, 0])
    ref, frames = samples(response, stress0)
    tensor, info = fit_cartesian_elastic_tensor(ref, frames)
    expected = (response + response.T) / 2
    expected += .02 * np.diag([1, 1, 1, 1, 1, 1])
    for i, j in [(0, 1), (0, 2), (1, 2)]:
        expected[i, j] -= .02
        expected[j, i] -= .02
    np.testing.assert_allclose(tensor, expected, atol=1e-12, rtol=0)
    np.testing.assert_allclose(info['response'], response, atol=1e-12, rtol=0)
    assert info['maximum_asymmetry'] == pytest.approx(.2)
    assert info['pressure'] == pytest.approx(.02)
    assert info['symmetric_rms_residual'] > info['rms_residual']


def test_reference_deviatoric_stress_is_subtracted_componentwise():
    c = np.eye(6)
    stress0 = np.array([.01, -.01, 0, .002, .003, .004])
    ref, frames = samples(c, stress0)
    fitted, _ = fit_cartesian_elastic_tensor(ref, frames)
    np.testing.assert_allclose(fitted, c, atol=1e-12)


def test_missing_shear_mode_is_rejected():
    ref, frames = samples(np.eye(6))
    with pytest.raises(RuntimeError, match='strain rank 5 of 6'):
        fit_cartesian_elastic_tensor(ref, frames[:-5])


@pytest.mark.parametrize('invalid', [np.nan, np.inf])
def test_nonfinite_stress_is_rejected(invalid):
    ref, frames = samples(np.eye(6))
    frames[0].calc = SinglePointCalculator(frames[0], energy=0, stress=[invalid] * 6)
    with pytest.raises(ValueError, match='six finite stresses'):
        fit_cartesian_elastic_tensor(ref, frames)


def test_empty_and_singular_data_are_rejected():
    ref, frames = samples(np.eye(6))
    with pytest.raises(ValueError, match='requires deformation'):
        fit_cartesian_elastic_tensor(ref, [])
    ref.cell[2] = ref.cell[1]
    with pytest.raises(ValueError, match='full rank'):
        fit_cartesian_elastic_tensor(ref, frames)


def test_archived_tic_retains_normal_shear_coupling():
    archive = Path(__file__).resolve().parents[1] / 'docs/source/_static/examples/tic-gpaw-elastic-audit-reference.zip'
    with zipfile.ZipFile(archive) as z:
        name = next(n for n in z.namelist() if n.endswith('deformations.traj'))
        frames = read(io.BytesIO(z.read(name)), format='traj', index=':')
    tensor, info = fit_cartesian_elastic_tensor(frames[2], frames)
    c = tensor / GPa
    np.testing.assert_allclose(np.diag(c),
                               [490.18008, 490.00502, 480.56197, 189.17956, 189.16439, 179.63834], atol=1e-5)
    assert c[0, 3] == pytest.approx(-13.31059, abs=1e-5)
    assert c[3, 0] == c[0, 3]
    assert info['rms_residual'] / GPa == pytest.approx(.04714523281)
    assert info['maximum_asymmetry'] / GPa == pytest.approx(.13574856)
