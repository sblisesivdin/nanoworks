# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Plot saved DFT results without importing an electronic-structure engine."""

from pathlib import Path
import re
from io import StringIO

import numpy as np


def load_optical_table(filename):
    """Read a seven-column optical DAT/CSV using its explicit column header."""
    from nanoworks.optical import OPTICAL_TABLE_HEADER, LEGACY_OPTICAL_TABLE_HEADER, validate_optical_table
    filename = Path(filename)
    delimiter = ',' if filename.suffix.lower() == '.csv' else None
    with filename.open(encoding='utf-8-sig') as stream:
        header = stream.readline().strip().removeprefix('#').strip()
        columns = [column.strip() for column in header.split(delimiter)]
        if columns not in (OPTICAL_TABLE_HEADER.split(), LEGACY_OPTICAL_TABLE_HEADER.split()):
            raise ValueError('Expected a seven-column optical table with Energy(eV), '
                             'dielectric, refractive/extinction, Absorption(1/cm) and Reflectivity headers.')
        table = np.loadtxt(stream, delimiter=delimiter, ndmin=2)
    return validate_optical_table(table)


def plot_optical_table(filename):
    """Plot a standalone optical table without a summary or calculator state."""
    from nanoworks.optical import generate_optical_figures
    filename = Path(filename).expanduser().resolve()
    table = load_optical_table(filename)
    prefix = filename.with_suffix('')
    figures = generate_optical_figures(table, str(prefix), 'Optical')
    return {'figure_files': {'Optical': figures}, 'output_directory': str(filename.parent)}


def load_dos_table(filename):
    """Read common GPAW/QE DOS DAT or legacy CSV, retaining its saved energy reference."""
    filename = Path(filename)
    lines = filename.read_text(encoding='utf-8-sig').splitlines()
    delimiter = ',' if filename.suffix.lower() == '.csv' else None
    if lines and lines[0].strip().removeprefix('#').strip().startswith('Energy('):
        lines = lines[1:]
    table = np.loadtxt(StringIO('\n'.join(lines)), delimiter=delimiter, ndmin=2)
    if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2
            or not np.isfinite(table).all() or np.any(np.diff(table[:, 0]) <= 0)):
        raise ValueError('DOS results require finite two-column DOS data on an increasing energy grid.')
    return table


def plot_dos_result(filename):
    """Validate both spin channels before replacing the existing DOS figure."""
    filename = Path(filename).expanduser().resolve()
    match = re.fullmatch(r'(.+-DOS-(?:GPAW|QE))-Result-DOS(?:-(Up|Down))?\.(?:dat|csv)', filename.name)
    if match is None:
        raise ValueError('Expected a Nanoworks *-DOS-GPAW/QE-Result-DOS[-Up/Down].dat file (legacy CSV is accepted).')
    prefix, spin = match.groups()
    extension = filename.suffix
    if spin:
        up = load_dos_table(filename.with_name(prefix + '-Result-DOS-Up' + extension))
        down = load_dos_table(filename.with_name(prefix + '-Result-DOS-Down' + extension))
        if up.shape != down.shape or not np.array_equal(up[:, 0], down[:, 0]):
            raise ValueError('Spin-resolved DOS tables must have identical energy grids.')
        channels = [(up, 1., 'Spin Up', 'blue'), (down, -1., 'Spin Down', 'red')]
    else:
        channels = [(load_dos_table(filename), 1., 'DOS', 'blue')]

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = filename.with_name(prefix + '-Graph-DOS.png')
    fig, ax = plt.subplots(figsize=(8, 6))
    try:
        for table, sign, label, color in channels:
            ax.plot(table[:, 0], sign * table[:, 1], color=color, linewidth=1.5, label=label)
            ax.fill_between(table[:, 0], 0, sign * table[:, 1], color=color, alpha=.2)
        ax.set_xlabel('Energy relative to Fermi level (eV)')
        ax.set_ylabel('Density of states (states/eV)')
        ax.axvline(0, color='k', linestyle='--', linewidth=1)
        ax.axhline(0, color='k', linewidth=.8)
        ax.set_xlim(channels[0][0][0, 0], channels[0][0][-1, 0])
        if spin:
            ax.legend()
        fig.tight_layout()
        fig.savefig(output, dpi=300)
    finally:
        plt.close(fig)
    return {'figure_files': {'DOS': [str(output)]}, 'output_directory': str(filename.parent)}


def plot_saved_result(filename):
    """Dispatch the common --plot command using supported saved result formats."""
    filename = Path(filename).expanduser().resolve()
    if filename.suffix.lower() == '.json':
        from nanoworks.optical import plot_optical_summary
        return plot_optical_summary(filename)
    if filename.suffix.lower() in ('.dat', '.csv') and re.fullmatch(
            r'.+-DOS-(?:GPAW|QE)-Result-DOS(?:-(Up|Down))?\.(?:dat|csv)', filename.name):
        return plot_dos_result(filename)
    if filename.suffix.lower() in ('.dat', '.csv'):
        return plot_optical_table(filename)
    raise ValueError('Supported saved results are optical DAT/CSV tables and Nanoworks DOS DAT files (legacy CSV is supported); '
                     'optical summary JSON files are optional.')
