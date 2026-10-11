# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Portable band-path coordinates and labels in an ordinary DAT table."""

from pathlib import Path
from urllib.parse import unquote

import numpy as np

from nanoworks.tableio import DatWriter


BAND_PATH_COLUMNS = ['KPointIndex(1)', 'Distance(1/A)', 'BreakBefore(1)', 'Label(1)']


def write_band_energy_tables(prefix, energies, spin_labels=None):
    """Stage long and wide DAT tables directly from one referenced energy array.

    Call on the output rank only. All channels are validated and staged before
    replacing any existing table; no process needs to read another rank's output.
    Explicit labels preserve the caller's established spin convention.
    """
    energies = np.asarray(energies, dtype=float)
    if (energies.ndim != 3 or energies.shape[0] not in (1, 2)
            or energies.shape[1] < 2 or energies.shape[2] < 1
            or not np.isfinite(energies).all()):
        raise ValueError('Band energies require finite (1 or 2 spins, kpoints, bands) data.')
    labels = list(spin_labels) if spin_labels is not None else ['']
    if (len(labels) != energies.shape[0]
            or (energies.shape[0] == 1 and labels != [''])
            or (energies.shape[0] == 2 and set(labels) != {'Up', 'Down'})):
        raise ValueError('Band spin labels must match the energy channels.')
    staged = []
    try:
        for label, channel in zip(labels, energies):
            suffix = '-' + label if label else ''
            for wide in (False, True):
                filename = Path(str(prefix) + '-Result-Band' + suffix
                                + ('-XYYY' if wide else '') + '.dat')
                filename.parent.mkdir(parents=True, exist_ok=True)
                temporary = filename.with_name(filename.name + '.tmp')
                staged.append((temporary, filename))
                with temporary.open('w', encoding='utf-8') as stream:
                    writer = DatWriter(stream)
                    if wide:
                        writer.writerow(['KPointIndex(1)'] + [
                            f'Band_{band + 1}(eV)' for band in range(channel.shape[1])])
                        for index, row in enumerate(channel):
                            writer.writerow([index] + [format(value, '.17g') for value in row])
                    else:
                        writer.writerow(['KPointIndex(1)', 'Energy(eV)'])
                        for band in range(channel.shape[1]):
                            for index, value in enumerate(channel[:, band]):
                                writer.writerow([index, format(value, '.17g')])
                            stream.write('\n')
        for temporary, filename in staged:
            temporary.replace(filename)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)
    return [str(filename) for _, filename in staged]


def _validate_axis(indices, distances, breaks):
    if (indices.ndim != 1 or distances.shape != indices.shape or breaks.shape != indices.shape
            or len(indices) < 2 or not np.isfinite(indices).all()
            or not np.isfinite(distances).all() or not np.isfinite(breaks).all()
            or not np.array_equal(indices, np.arange(len(indices)))
            or np.any(distances < 0) or np.any(np.diff(distances) < 0)
            or distances[-1] <= distances[0] or np.any(~np.isin(breaks, [0, 1]))
            or breaks[0] != 0):
        raise ValueError('Band paths require consecutive indices, nondecreasing finite distances and 0/1 break flags.')


def write_band_path_table(filename, distances, special_distances, labels):
    """Export one shared axis for all spin and band-energy tables."""
    distances = np.asarray(distances, dtype=float)
    special = np.asarray(special_distances, dtype=float)
    if (distances.ndim != 1 or special.ndim != 1 or len(special) != len(labels)
            or not np.isfinite(special).all()):
        raise ValueError('Band-path tick distances must match their labels.')
    indices = np.arange(len(distances))
    breaks = np.r_[0, np.diff(distances) == 0].astype(int)
    _validate_axis(indices, distances, breaks)
    row_labels = [[] for _ in indices]
    for position, label in zip(special, labels):
        if not isinstance(label, str) or not label.strip() or label == '-':
            raise ValueError('Band-path special-point labels must be nonempty strings.')
        matches = np.flatnonzero(np.isclose(distances, position, rtol=1e-8, atol=1e-10))
        if not len(matches):
            raise ValueError('Band-path special points must occur on the sampled distance axis.')
        if label not in row_labels[matches[0]]:
            row_labels[matches[0]].append(label)
    filename = Path(filename)
    filename.parent.mkdir(parents=True, exist_ok=True)
    temporary = filename.with_name(filename.name + '.tmp')
    try:
        with temporary.open('w', encoding='utf-8') as stream:
            writer = DatWriter(stream)
            writer.writerow(BAND_PATH_COLUMNS)
            for index, distance, before, names in zip(indices, distances, breaks, row_labels):
                writer.writerow([index, format(distance, '.17g'), before, ','.join(names) or '-'])
        temporary.replace(filename)
    finally:
        temporary.unlink(missing_ok=True)
    return str(filename)


def load_band_path_table(filename, expected_indices):
    """Require an axis matching every saved k-point before using its labels."""
    lines = Path(filename).read_text(encoding='utf-8-sig').splitlines()
    if not lines or lines[0].split() != BAND_PATH_COLUMNS:
        raise ValueError('Unsupported band-path DAT header.')
    rows = [line.split() for line in lines[1:] if line.strip()]
    if not rows or any(len(row) != len(BAND_PATH_COLUMNS) for row in rows):
        raise ValueError('Band-path rows must match the four header columns.')
    values = np.asarray([row[:3] for row in rows], dtype=float)
    indices, distances, breaks = values.T
    _validate_axis(indices, distances, breaks)
    if not np.array_equal(indices, expected_indices):
        raise ValueError('Band-path indices do not match the saved band-energy table.')
    ticks = {}
    for distance, row in zip(distances, rows):
        label = unquote(row[3])
        if label != '-':
            ticks.setdefault(float(distance), []).append(label)
    return {'distances': distances, 'breaks': breaks.astype(bool),
            'ticks': list(ticks), 'labels': [','.join(names) for names in ticks.values()]}
