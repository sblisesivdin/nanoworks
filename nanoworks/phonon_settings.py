# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Portable phonon parameter validation before electronic calculations."""

import math
import operator

import numpy as np


def validate_atomic_masses(masses, natoms):
    """Require one finite positive atomic mass per phonon unit-cell site."""
    try:
        values = np.asarray(masses, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('Phonon atomic masses must be finite positive values for every atom.') from exc
    if natoms <= 0 or values.shape != (natoms,) or not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('Phonon atomic masses must be finite positive values for every atom.')
    return values.copy()


def positive_integer(value, name, minimum=1):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f'{name} must be an integer >= {minimum}.')
    try:
        result = operator.index(value)
    except TypeError as exc:
        raise ValueError(f'{name} must be an integer >= {minimum}.') from exc
    if result < minimum:
        raise ValueError(f'{name} must be an integer >= {minimum}.')
    return result


def finite_number(value, name, minimum=0., strict=False):
    try:
        if isinstance(value, (bool, np.bool_, str, bytes)):
            raise ValueError
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f'{name} must be a finite number.') from exc
    if not math.isfinite(result) or result < minimum or (strict and result == minimum):
        relation = '>' if strict else '>='
        raise ValueError(f'{name} must be finite and {relation} {minimum:g}.')
    return result


def normalize_supercell(value):
    try:
        matrix = np.asarray(value, dtype=object)
    except (TypeError, ValueError) as exc:
        raise ValueError('Phonon_supercell must contain three integers or a 3x3 integer matrix.') from exc
    if matrix.shape == (3,):
        matrix = np.diag([positive_integer(item, 'Phonon_supercell') for item in matrix])
    if matrix.shape != (3, 3):
        raise ValueError('Phonon_supercell must contain three integers or a 3x3 integer matrix.')
    rows = []
    for row in matrix:
        try:
            if any(isinstance(item, (bool, np.bool_)) for item in row):
                raise TypeError
            rows.append([operator.index(item) for item in row])
        except TypeError as exc:
            raise ValueError('Phonon_supercell matrix values must be integers.') from exc
    a, b, c = rows
    determinant = (a[0] * (b[1] * c[2] - b[2] * c[1])
                   - a[1] * (b[0] * c[2] - b[2] * c[0])
                   + a[2] * (b[0] * c[1] - b[1] * c[0]))
    if determinant <= 0:
        raise ValueError('Phonon_supercell must have a positive determinant.')
    try:
        normalized = np.asarray(rows, dtype=np.int64)
    except OverflowError as exc:
        raise ValueError('Phonon_supercell entries exceed the supported integer range.') from exc
    return normalized, determinant


def validate_phonon_settings(config, engine=None, validate_structure=True):
    engine = engine or config.Engine
    finite_displacement = engine == 'GPAW' or bool(config.Hubbard_U)
    matrix, multiplier = normalize_supercell(config.Phonon_supercell)
    if not finite_displacement:
        if not np.array_equal(matrix, np.diag(np.diag(matrix))) or np.any(np.diag(matrix) <= 0):
            raise ValueError('Native QE DFPT requires a positive diagonal Phonon_supercell.')
    result = {'Phonon_supercell': matrix}
    for name in ('Phonon_qpts_x', 'Phonon_qpts_y', 'Phonon_qpts_z', 'Phonon_npoints'):
        result[name] = positive_integer(getattr(config, name), name, 2 if name == 'Phonon_npoints' else 1)
    for name in ('Phonon_kpts_x', 'Phonon_kpts_y', 'Phonon_kpts_z'):
        value = getattr(config, name)
        result[name] = None if value is None else positive_integer(value, name)
    if config.Phonon_PW_cutoff is not None:
        result['Phonon_PW_cutoff'] = finite_number(config.Phonon_PW_cutoff, 'Phonon_PW_cutoff', strict=True)
    if finite_displacement:
        result['Phonon_displacement'] = finite_number(config.Phonon_displacement, 'Phonon_displacement', strict=True)
    if config.Phonon_thermal_calc:
        for name in ('Phonon_T_min', 'Phonon_T_max', 'Phonon_T_step'):
            result[name] = finite_number(getattr(config, name), name, strict=name == 'Phonon_T_step')
        if result['Phonon_T_max'] < result['Phonon_T_min']:
            raise ValueError('Phonon_T_max must be >= Phonon_T_min.')
    details = {'method': 'finite-displacement' if finite_displacement else 'DFPT',
        'supercell_matrix': matrix.tolist(), 'dos_mesh': [result[f'Phonon_qpts_{axis}'] for axis in 'xyz'],
        'thermal_enabled': bool(config.Phonon_thermal_calc)}
    if finite_displacement:
        details['supercell_multiplier'] = multiplier
        details['displacement_angstrom'] = result['Phonon_displacement']
        atoms = getattr(config, 'bulk_configuration', None)
        if atoms is not None:
            details['supercell_atoms'] = len(atoms) * multiplier
            if validate_structure and hasattr(atoms, 'get_masses'):
                masses = validate_atomic_masses(atoms.get_masses(), len(atoms))
                details['atomic_masses_amu'] = masses.tolist()
                details['mass_source'] = ('ASE explicit masses' if atoms.has('masses')
                                          else 'ASE elemental defaults')
    else:
        details['dfpt_qpoint_grid'] = np.diag(matrix).tolist()
    if config.Phonon_thermal_calc:
        details['temperature_range_kelvin'] = [result[name] for name in
            ('Phonon_T_min', 'Phonon_T_max', 'Phonon_T_step')]
    return result, details
