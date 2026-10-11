# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Plot saved DFT results without importing an electronic-structure engine."""

from pathlib import Path
import re

import numpy as np


def load_dos_table(filename):
    """Read the common GPAW/QE DOS CSV, retaining its saved energy reference."""
    table = np.loadtxt(filename, delimiter=',', ndmin=2)
    if (table.ndim != 2 or table.shape[1] != 2 or len(table) < 2
            or not np.isfinite(table).all() or np.any(np.diff(table[:, 0]) <= 0)):
        raise ValueError('DOS results require finite two-column CSV data on an increasing energy grid.')
    return table


def plot_dos_result(filename):
    """Validate both spin channels before replacing the existing DOS figure."""
    filename = Path(filename).expanduser().resolve()
    match = re.fullmatch(r'(.+-DOS-(?:GPAW|QE))-Result-DOS(?:-(Up|Down))?\.csv', filename.name)
    if match is None:
        raise ValueError('Expected a Nanoworks *-DOS-GPAW/QE-Result-DOS[-Up/Down].csv file.')
    prefix, spin = match.groups()
    if spin:
        up = load_dos_table(filename.with_name(prefix + '-Result-DOS-Up.csv'))
        down = load_dos_table(filename.with_name(prefix + '-Result-DOS-Down.csv'))
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
    if filename.suffix.lower() == '.csv':
        return plot_dos_result(filename)
    raise ValueError('Supported saved results are optical summary JSON files and Nanoworks DOS CSV files.')
