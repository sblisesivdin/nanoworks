# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Headless regression checks for archived QE WSe2 band tutorial data."""

import hashlib
from pathlib import Path

import numpy as np
import pytest


DATA = Path(__file__).resolve().parents[1] / "docs/source/_static/examples"
FILES = {
    "scalar": (
        "wse2-qe-bands-scalar.dat", 28, 23,
        "ec1bb048b3b92dafca3582d4b87f719dcd54154bd3e1c676fa7fac0d31db2343",
    ),
    "soc": (
        "wse2-qe-bands-soc.dat", 56, 46,
        "665fccab207e5c4d0f6860e4b6194dedbee9f12d0ebec4e79fb416c41958b3ec",
    ),
}


def _read_bands(path, nbands):
    rows = np.loadtxt(path, skiprows=1)
    assert rows.shape == (100 * nbands, 2)
    bands = np.stack([
        np.sort(rows[rows[:, 0] == i, 1]) for i in range(100)
    ])
    assert bands.shape == (100, nbands)
    assert np.isfinite(bands).all()
    return bands


@pytest.mark.parametrize("mode", ["scalar", "soc"])
def test_recorded_band_tables_match_archive_and_path_gaps(mode):
    filename, nbands, filled, expected_sha = FILES[mode]
    source = DATA / filename
    assert hashlib.sha256(source.read_bytes()).hexdigest() == expected_sha
    bands = _read_bands(source, nbands)
    valence = bands[:, filled - 1]
    conduction = bands[:, filled]
    assert int(np.argmax(valence)) == 0  # Gamma
    assert int(np.argmin(conduction)) == 56  # K
    expected_indirect = {"scalar": 1.2476, "soc": 1.1911}[mode]
    expected_direct_k = {"scalar": 1.5916, "soc": 1.3107}[mode]
    assert conduction.min() - valence.max() == pytest.approx(expected_indirect, abs=1e-8)
    assert conduction[56] - valence[56] == pytest.approx(expected_direct_k, abs=1e-8)


def test_recorded_soc_k_splittings_use_46_occupied_spinors():
    bands = _read_bands(DATA / FILES["soc"][0], 56)
    assert bands[56, 45] - bands[56, 44] == pytest.approx(0.5034, abs=1e-8)
    assert bands[56, 47] - bands[56, 46] == pytest.approx(0.0070, abs=1e-8)


def test_recorded_path_is_identical_and_special_points_are_sampled():
    source = DATA / "wse2-qe-kpath.dat"
    assert hashlib.sha256(source.read_bytes()).hexdigest() == (
        "122242ba20d076a076fbac2bfff8cc97cab15822783efdcb096e55faabdc7937"
    )
    rows = [line.split() for line in source.read_text().splitlines()[1:]]
    assert len(rows) == 100
    assert [int(row[0]) for row in rows] == list(range(100))
    assert [(int(row[0]), row[-1]) for row in rows if row[-1] != "-"] == [
        (0, "G"), (36, "M"), (56, "K"), (99, "G"),
    ]
