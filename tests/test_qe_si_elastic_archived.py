# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Regression audit of the recorded QE/thermo_pw silicon elasticity smoke run."""

import hashlib
from pathlib import Path

import numpy as np
import pytest

from nanoworks.elasticity import analyze_elastic_stability


DATA = Path(__file__).resolve().parents[1] / "docs/source/_static/examples"
CANONICAL = DATA / "si-qe-elastic-canonical.dat"
TENSOR = DATA / "si-qe-elastic-tensor-gpa.dat"
ARCHIVE_SHA256 = (
    "b13d4a306bd2728e01f470cefe1e6069c7496a955c9e44e58dd60f5c1d56e646"
)


def _archived_tensor():
    # thermo_pw canonical .g1 contains 36 stiffness and 36 compliance
    # values in kbar and 1/Mbar, respectively; only the first 36 are Cij.
    values = np.fromstring(CANONICAL.read_text(encoding="utf-8"), sep=" ")
    assert values.shape == (72,)
    return values[:36].reshape(6, 6) / 10  # kbar -> GPa


def test_silicon_thermo_pw_raw_data_and_export_are_identical():
    assert hashlib.sha256(CANONICAL.read_bytes()).hexdigest() == ARCHIVE_SHA256
    canonical = _archived_tensor()
    recorded = np.loadtxt(TENSOR)
    assert recorded.shape == (6, 6)
    np.testing.assert_allclose(recorded, canonical, rtol=0, atol=6e-9)


def test_silicon_smoke_tensor_fails_major_symmetry_quality_check():
    tensor = _archived_tensor()
    diagnostic = analyze_elastic_stability(tensor)
    assert diagnostic["maximum_tensor_asymmetry_gpa"] == pytest.approx(
        1.46196117, abs=1e-7
    )
    assert diagnostic["tensor_symmetry_tolerance_gpa"] == pytest.approx(
        0.002003388413, abs=1e-9
    )
    assert not diagnostic["tensor_symmetric"]
    # Positive eigenvalues of (C+C.T)/2 are not a reliable material stability
    # claim while the raw stress-strain tensor lacks major symmetry.
    assert diagnostic["mechanically_stable"]
    assert diagnostic["interpretation_status"] == "inconclusive"
    assert diagnostic["minimum_eigenvalue"] == pytest.approx(
        100.1403234, abs=1e-6
    )


def test_silicon_smoke_cubic_averages_are_diagnostic_only():
    tensor = _archived_tensor()
    cubic_C11 = np.mean(np.diag(tensor)[:3])
    cubic_C12 = tensor[:3, :3][~np.eye(3, dtype=bool)].mean()
    cubic_C44 = np.diag(tensor)[3:].mean()
    assert cubic_C11 == pytest.approx(199.7021091667, abs=1e-7)
    assert cubic_C12 == pytest.approx(78.9518324467, abs=1e-7)
    assert cubic_C44 == pytest.approx(100.1786682667, abs=1e-7)


def test_silicon_archived_summary_reports_asymmetry_and_hill_moduli():
    summary = (DATA / "si-qe-elastic-allresults.txt").read_text(encoding="utf-8")
    assert "Maximum |Cij-Cji|: 1.461961 GPa" in summary
    assert "Tensor symmetric within tolerance: no" in summary
    assert "Mechanically stable: yes" in summary  # historical wording only
    assert "Bulk modulus: 119.435569 GPa" in summary
    assert "Shear modulus: 81.737519 GPa" in summary
    assert "Young modulus: 199.635001 GPa" in summary
