# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Verify the numerical QE Si DFPT tutorial exports without running pw.x."""

import hashlib
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).resolve().parents[1] / "docs/source/_static/examples"
HASHES = {
    "si-qe-dfpt-band-thz.dat": "fb2912a0e6d57ccc61d31780c106997b3f40f914d081ac15bdfbfc9651c697d5",
    "si-qe-dfpt-dos-total.dat": "213314f7fbc4101c1802074afa332852fa1bbe897ea7ad23db417473ee57d089",
    "si-qe-dfpt-final.cif": "372a4d3ca5c522f5c370ce3eb95ac5e25dca49a54cd206ac075e91525255765f",
}


@pytest.mark.parametrize("name,sha256", HASHES.items())
def test_archived_si_dfpt_data_hashes(name, sha256):
    assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == sha256


def test_six_si_branches_and_gamma_acoustic_sum_rule():
    b = np.loadtxt(DATA / "si-qe-dfpt-band-thz.dat", skiprows=1)
    assert b.shape == (61, 10)
    assert np.isfinite(b).all()
    np.testing.assert_allclose(b[[0, 60], 4:7], 0, atol=1e-9)
    np.testing.assert_allclose(b[[0, 60], 7:], np.tile([15.08844049, 15.08844049, 15.0884464859], (2, 1)), atol=1e-8)
    assert b[:, 4:].min() >= -1e-9
    assert b[:, 4:].max() == pytest.approx(15.0884464859, abs=1e-8)


def test_recorded_si_path_endpoint_indices():
    b = np.loadtxt(DATA / "si-qe-dfpt-band-thz.dat", skiprows=1)
    for index, q in [(0, [0, 0, 0]), (21, [.5, -.288675, .408248]),
                     (28, [.625, -.072169, .408248]), (41, [.5, .288675, .204124]),
                     (60, [0, 0, 0])]:
        np.testing.assert_allclose(b[index, 1:4], q, atol=1e-6)


def test_dos_integral_close_to_six_modes_per_primitive_cell():
    dos = np.loadtxt(DATA / "si-qe-dfpt-dos-total.dat")
    assert dos.shape == (505, 2)
    assert np.isfinite(dos).all()
    assert np.all(np.diff(dos[:, 0]) > 0)
    assert np.min(dos[:, 1]) >= 0
    integral = np.sum(np.diff(dos[:, 0]) * (dos[:-1, 1] + dos[1:, 1]) / 2)
    assert integral == pytest.approx(5.992795378535916, abs=1e-8)
