# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Whitespace-separated result tables with one token per named column."""

import csv
import re
from urllib.parse import quote


_UNITS = {
    'Energy': 'Energy(eV)', 'energy_ev': 'Energy(eV)', 'energy_eV': 'Energy(eV)',
    'total_energy_ev': 'TotalEnergy(eV)', 'energy_ev_per_atom': 'EnergyPerAtom(eV/atom)',
    'cell_volume_angstrom3': 'CellVolume(A^3)', 'DOS': 'DOS(states/eV)',
    'distance': 'Distance(1/A)', 'weight': 'Weight(1)', 'energy_eV_per_atom': 'EnergyPerAtom(eV/atom)',
    'parameter_value': 'ParameterValue(eV_or_points/A_or_ratio_or_mesh)',
    **{label: label + '(states/eV)' for label in
       ('s-total','p-total','pz','px','py','d-total','d3z2_r2','dzx','dxz','dyz','dx2_y2','dxy','f-total','TOTAL')},
    'engineering_strain': 'EngineeringStrain(1)', 'g(r)': 'g_r(1)',
}


def table_header(columns):
    """Return a single header row; spaces in labels never become extra columns."""
    if isinstance(columns, str):
        columns = next(csv.reader([columns])) if ',' in columns else columns.split()
    tokens = []
    for column in columns:
        column = str(column).strip().removeprefix('#').strip()
        column = _UNITS.get(column, column)
        if column.endswith('_GPa'):
            column = column[:-4] + '(GPa)'
        column = re.sub(r'\s+', '_', column)
        if not column:
            raise ValueError('Table columns must have nonempty names.')
        if '(' not in column and '[' not in column:
            column += '(1)'
        tokens.append(column)
    if len(set(tokens)) != len(tokens):
        raise ValueError('Table column names must be unique.')
    return ' '.join(tokens)


def table_row(values):
    """Keep missing entries explicit and escape whitespace in text fields."""
    return ' '.join('nan' if value is None or str(value) == '' else
                    quote(str(value), safe='+-./:;_=()[]') for value in values)


def csv_line_to_dat(line):
    """Adapt existing formatted MD rows without changing their numeric precision."""
    return table_row(next(csv.reader([line.rstrip('\n')]))) + '\n'


class DatWriter:
    """Small row-writer interface for result exporters previously using CSV."""

    def __init__(self, stream):
        self.stream = stream
        self.ncolumns = None

    def writerow(self, row):
        row = list(row)
        if self.ncolumns is None:
            self.ncolumns = len(row)
            self.stream.write(table_header(row) + '\n')
        else:
            if len(row) != self.ncolumns:
                raise ValueError('Result table row does not match its column header.')
            self.stream.write(table_row(row) + '\n')
