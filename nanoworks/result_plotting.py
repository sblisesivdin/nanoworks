# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Plot saved DFT results without importing an electronic-structure engine."""

from pathlib import Path
import re
from io import StringIO

import numpy as np


def _save_figure(figure, filename):
    """Replace a saved PNG only after Matplotlib finishes writing it."""
    filename = Path(filename)
    temporary = filename.with_name(filename.name + '.tmp')
    try:
        figure.savefig(temporary, format='png', dpi=300)
        temporary.replace(filename)
    finally:
        temporary.unlink(missing_ok=True)


def load_band_table(filename):
    """Read long band blocks or wide XYYY columns without a calculator state."""
    filename = Path(filename)
    lines = filename.read_text(encoding='utf-8-sig').splitlines()
    wide = filename.stem.endswith('-XYYY')
    if lines and lines[0].strip().removeprefix('#').strip().startswith('KPointIndex('):
        header = lines.pop(0).strip().removeprefix('#').strip().split()
        if wide:
            if (len(header) < 2 or header[0] != 'KPointIndex(1)'
                    or header[1:] != [f'Band_{i + 1}(eV)' for i in range(len(header) - 1)]):
                raise ValueError('Invalid wide band table header.')
        elif header != ['KPointIndex(1)', 'Energy(eV)']:
            raise ValueError('Invalid long band table header.')
    else:
        header = None  # Older Nanoworks band tables had no header.
    if wide:
        table = np.loadtxt(StringIO('\n'.join(lines)), ndmin=2)
        if table.shape[1] < 2 or (header is not None and table.shape[1] != len(header)):
            raise ValueError('Wide band tables must match their numbered energy columns.')
        indices, energies = table[:, 0], table[:, 1:]
    else:
        blocks, current = [], []
        for line in lines + ['']:
            if line.strip():
                current.append(line)
            elif current:
                block = np.loadtxt(StringIO('\n'.join(current)), ndmin=2)
                if block.shape[1] != 2:
                    raise ValueError('Long band tables require two columns per band block.')
                blocks.append(block)
                current = []
        if not blocks or any(not np.array_equal(block[:, 0], blocks[0][:, 0]) for block in blocks):
            raise ValueError('Band blocks must have identical k-point grids.')
        indices = blocks[0][:, 0]
        energies = np.column_stack([block[:, 1] for block in blocks])
    if (len(indices) < 2 or not np.isfinite(indices).all() or not np.isfinite(energies).all()
            or np.any(indices < 0) or np.any(indices != np.floor(indices))
            or np.any(np.diff(indices) <= 0)):
        raise ValueError('Band tables require finite energies on increasing nonnegative integer k-point indices.')
    return indices, energies


def plot_band_result(filename):
    """Draw saved band energies on their recorded k-point index axis."""
    filename = Path(filename).expanduser().resolve()
    match = re.fullmatch(r'(.+-BAND-(?:GPAW|QE))-Result-Band(?:-(Up|Down))?(-XYYY)?\.dat', filename.name)
    if match is None:
        raise ValueError('Expected a Nanoworks *-BAND-GPAW/QE-Result-Band[-Up/Down][-XYYY].dat file.')
    prefix, spin, wide = match.groups()
    if spin:
        tables = [load_band_table(filename.with_name(prefix + '-Result-Band-' + channel
                                                   + (wide or '') + '.dat')) for channel in ('Up', 'Down')]
        if (not np.array_equal(tables[0][0], tables[1][0])
                or tables[0][1].shape != tables[1][1].shape):
            raise ValueError('Spin-resolved band tables must have identical k-point grids and band counts.')
    else:
        tables = [load_band_table(filename)]
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = filename.with_name(prefix + '-Graph-Band-Index.png')
    fig, ax = plt.subplots(figsize=(8, 6))
    try:
        for (indices, energies), color, label in zip(tables, ('blue', 'red'), ('Spin Up', 'Spin Down')):
            for band in range(energies.shape[1]):
                ax.plot(indices, energies[:, band], color=color, linewidth=1,
                        label=label if spin and band == 0 else None)
        ax.set_xlabel('K-point index')
        ax.set_ylabel('Energy relative to saved reference (eV)')
        ax.axhline(0, color='k', linestyle='--', linewidth=.8)
        ax.set_xlim(tables[0][0][0], tables[0][0][-1])
        if spin:
            ax.legend()
        fig.tight_layout()
        _save_figure(fig, output)
    finally:
        plt.close(fig)
    return {'figure_files': {'Band': [str(output)]}, 'output_directory': str(filename.parent)}


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
        _save_figure(fig, output)
    finally:
        plt.close(fig)
    return {'figure_files': {'DOS': [str(output)]}, 'output_directory': str(filename.parent)}


def plot_saved_result(filename):
    """Dispatch the common --plot command using supported saved result formats."""
    filename = Path(filename).expanduser().resolve()
    if filename.suffix.lower() == '.json':
        from nanoworks.optical import plot_optical_summary
        return plot_optical_summary(filename)
    if re.fullmatch(r'.+-BAND-(?:GPAW|QE)-Result-Band(?:-(Up|Down))?(-XYYY)?\.dat', filename.name):
        return plot_band_result(filename)
    if filename.suffix.lower() in ('.dat', '.csv') and re.fullmatch(
            r'.+-DOS-(?:GPAW|QE)-Result-DOS(?:-(Up|Down))?\.(?:dat|csv)', filename.name):
        return plot_dos_result(filename)
    if filename.suffix.lower() in ('.dat', '.csv'):
        return plot_optical_table(filename)
    raise ValueError('Supported saved results are optical, DOS and band DAT files (legacy CSV is supported); '
                     'optical summary JSON files are optional.')
