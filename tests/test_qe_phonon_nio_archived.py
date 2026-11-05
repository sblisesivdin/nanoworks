# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Independent numerical checks for the archived NiO QE Hubbard tutorial."""

import hashlib
from pathlib import Path

import numpy as np
import pytest


DATA = Path(__file__).resolve().parents[1] / "docs/source/_static/examples"
REFERENCES = {
    "nio-qe-hubbard-band-thz.dat":
        "979abe1ddbfb06a35bbc8de75ba0673b5a7664a1dbf0664b571ea4522c321531",
    "nio-qe-hubbard-dos-thz.dat":
        "5dbc67e77c80736bd28f8368d62c9cb5e19c8828fc20b661b06a1c199d0c0517",
    "nio-qe-hubbard-thermal.dat":
        "11797716ac445ff504eb121775d2462b9ed33dedd081ae2435cec474def38d04",
}


@pytest.mark.parametrize("filename,expected", REFERENCES.items())
def test_archived_nio_tables_preserve_original_sha256(filename, expected):
    assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == expected


def test_archived_nio_phonon_band_count_and_segment_boundaries():
    band = np.loadtxt(DATA / "nio-qe-hubbard-band-thz.dat", skiprows=1)
    assert band.shape == (40, 16)  # distance, fractional q, 12 frequencies
    assert np.isfinite(band).all()
    assert np.max(np.abs(band[0, 4:7])) < 1e-5  # three Gamma acoustic modes
    assert band[:, 4:].min() > -1e-5
    assert band[:, 4:].max() == pytest.approx(16.2245502158)
    assert np.count_nonzero(np.isclose(np.diff(band[:, 0]), 0, atol=1e-9)) == 4


def test_archived_nio_dos_integral_and_atom_projections():
    dos = np.loadtxt(DATA / "nio-qe-hubbard-dos-thz.dat", skiprows=1)
    assert dos.shape == (201, 6)
    assert np.isfinite(dos).all()
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    assert integrate(dos[:, 1], dos[:, 0]) == pytest.approx(12.0, abs=1e-4)
    np.testing.assert_allclose(
        integrate(dos[:, 2:], dos[:, 0], axis=0), [3, 3, 3, 3],
        rtol=0, atol=1e-4,
    )
    np.testing.assert_allclose(dos[:, 2:].sum(axis=1), dos[:, 1],
                               rtol=0, atol=1e-8)


def test_archived_nio_harmonic_thermal_table():
    thermal = np.loadtxt(DATA / "nio-qe-hubbard-thermal.dat", skiprows=1)
    assert thermal.shape == (101, 5)
    np.testing.assert_allclose(thermal[:, 0], np.arange(0, 1001, 10))
    assert thermal[0, 1] == pytest.approx(23.9574185065)
    assert thermal[30, 1] == pytest.approx(15.8199302, abs=1e-5)
    assert thermal[30, 3] == pytest.approx(68.97111529, abs=1e-5)
    assert thermal[30, 4] == pytest.approx(80.13272975, abs=1e-5)
