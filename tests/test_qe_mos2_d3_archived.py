# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Reproduce numerical diagnostics from the rejected QE MoS2-D3 run."""

import hashlib
from pathlib import Path
import re

import numpy as np
import pytest


DATA = Path(__file__).resolve().parents[1] / "docs/source/_static/examples"
HASHES = {
    "mos2-qe-d3-bands.dat":
        "f7016a0bbbe8cf9afd529c7b603119c9f7070ddfe0cf2903f2db763aa9712d43",
    "mos2-qe-d3-kpath.dat":
        "968badee9a6070273ca6a6daeffc6ff92a9dcf1210115a58c09a0dd6bc81cf7e",
    "mos2-qe-d3-dos.dat":
        "841b29265c5d77e263b065adbd6653dc35a5c13f929b59482769ed85320a84ae",
    "mos2-qe-d3-pdos.dat":
        "1dee642a04ed4f9d1f2cf989a01199e8ff34b75e08ccd0779281d447f288a13c",
    "mos2-qe-d3-unverified-final.cif":
        "5326259ad6b575f54bd54c06e822c0f6eab43461a921e608efd6db00a1fe5041",
}


@pytest.mark.parametrize("filename,expected", HASHES.items())
def test_archive_files_are_unmodified(filename, expected):
    assert hashlib.sha256((DATA / filename).read_bytes()).hexdigest() == expected


def _cell_length(cif, name):
    text = Path(cif).read_text(encoding="utf-8")
    match = re.search(r"^_cell_length_" + name + r"\s+([0-9.]+)",
                      text, flags=re.MULTILINE)
    assert match is not None
    return float(match.group(1))


def test_archived_d3_cell_is_not_a_converged_reference():
    final = DATA / "mos2-qe-d3-unverified-final.cif"
    a, c = _cell_length(final, "a"), _cell_length(final, "c")
    assert a == pytest.approx(2.903446769, abs=1e-9)
    assert c == pytest.approx(10.694214657, abs=1e-9)
    assert 100 * (a / 3.16 - 1) == pytest.approx(-8.118773, abs=1e-5)
    assert 100 * (c / 12.29 - 1) == pytest.approx(-12.984421, abs=1e-5)


def test_recorded_d3_band_path_overlap_is_not_a_valid_gap():
    band = np.loadtxt(DATA / "mos2-qe-d3-bands.dat", skiprows=1)
    assert band.shape == (61, 32)
    # QE run used 52 electrons, hence 26 occupied non-spin-polarized bands.
    valence = band[:, 26]
    conduction = band[:, 27]
    assert int(valence.argmax()) == 0  # Gamma
    assert int(conduction.argmin()) == 33  # K
    assert valence.max() == pytest.approx(0.4022, abs=1e-8)
    assert conduction.min() == pytest.approx(-0.1715, abs=1e-8)
    assert conduction.min() - valence.max() == pytest.approx(-0.5737, abs=1e-8)
    assert conduction[33] - valence[33] == pytest.approx(0.5123, abs=1e-8)


def test_recorded_d3_dos_and_pdos_sampling_ranges_differ():
    dos = np.loadtxt(DATA / "mos2-qe-d3-dos.dat", skiprows=1)
    pdos = np.loadtxt(DATA / "mos2-qe-d3-pdos.dat", skiprows=1)
    assert dos.shape == (501, 2)
    assert pdos.shape == (424, 14)
    np.testing.assert_allclose(dos[:len(pdos), 0], pdos[:, 0], atol=1e-12)
    assert dos[250, 0] == pytest.approx(0.0003, abs=1e-9)
    assert dos[250, 1] == pytest.approx(0.1256, abs=1e-8)


def test_archived_d3_kpath_contains_61_samples():
    path = (DATA / "mos2-qe-d3-kpath.dat").read_text(encoding="utf-8").splitlines()
    assert len(path) == 62  # header plus 61 records
    labeled = [(int(cols[0]), cols[-1]) for line in path[1:]
               if (cols := line.split())[-1] != "-"]
    assert labeled == [(0, "G"), (21, "M"), (33, "K"), (60, "G")]
