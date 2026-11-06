# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Numerical checks on the recorded QE 7.4.1 Al DFPT phonon exports."""

import hashlib
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).resolve().parents[1] / "docs/source/_static/examples"

HASHES = {
    "al-qe-dfpt-band-thz.dat":
        "88fe89b9e1f0086287f05509c5c44303d53e0be0e7121323d7f37d869dee8212",
    "al-qe-dfpt-dos-thz.dat":
        "de16294add000cb196dc718991758191c653e0d5743c9d472c3efe8d0eeca5f5",
    "al-qe-dfpt-thermal.dat":
        "18766876841480284dabe433c0d50673472bc72c06ab6450b5cc3c07766b0e86",
    "al-qe-dfpt-final.cif":
        "e744c6f4f1ebc21b9e45d1140bd8d9a39cdb47a32497730e6d5210722a73b17f",
}


@pytest.mark.parametrize("name,expected", HASHES.items())
def test_al_qe_dfpt_source_files_preserved(name, expected):
    assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == expected


def test_al_qe_dfpt_three_acoustic_branches_and_path():
    data = np.loadtxt(DATA / "al-qe-dfpt-band-thz.dat", skiprows=1)
    assert data.shape == (61, 7)
    np.testing.assert_allclose(data[[0, 45], 4:], 0, atol=1e-9)
    assert np.isfinite(data).all()
    assert data[:, 4:].min() >= -1e-8
    assert data[:, 4:].max() == pytest.approx(11.6353500148, abs=1e-8)
    assert int(np.argmax(data[:, 4:].max(axis=1))) == 60
    np.testing.assert_allclose(data[15, 4:], [
        5.6560853736, 5.6568048755, 11.1692386969,
    ], atol=1e-8)


def test_al_qe_dfpt_dos_normalization_and_one_atom_projection():
    data = np.loadtxt(DATA / "al-qe-dfpt-dos-thz.dat", skiprows=1)
    assert data.shape == (389, 3)
    assert np.isfinite(data).all()
    assert np.min(data[:, 1:]) >= 0
    assert np.trapezoid(data[:, 1], data[:, 0]) == pytest.approx(
        2.99900626, abs=1e-7
    )
    np.testing.assert_allclose(data[:, 1], data[:, 2], atol=2.2e-5)


def test_al_qe_dfpt_recorded_vibrational_thermodynamics():
    data = np.loadtxt(DATA / "al-qe-dfpt-thermal.dat", skiprows=1)
    assert data.shape == (101, 5)
    np.testing.assert_allclose(data[:, 0], np.arange(0, 1001, 10))
    assert data[0, 1] == pytest.approx(4.0430071172)
    assert data[0, 3] == 0
    assert data[30, 4] == pytest.approx(22.41365185, abs=1e-7)
    assert data[-1, 4] == pytest.approx(24.6879071, abs=1e-7)
    assert data[-1, 4] < 3 * 8.314462618
