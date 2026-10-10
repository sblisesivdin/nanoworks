# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Backend-independent validation of optical inputs and numerical results."""

import math
from numbers import Integral, Real

import numpy as np


def validate_optical_settings(config):
    """Validate only parameters used by the selected optical method."""
    method = str(config.Opt_calc_type).strip().upper()
    if method not in ('RPA', 'BSE'):
        raise ValueError('Opt_calc_type must be RPA or BSE.')
    result = {'Opt_calc_type': method}

    def number(name, minimum=None, strict=False):
        value = getattr(config, name)
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
            raise ValueError(f'{name} must be a finite real number.')
        value = float(value)
        if (not math.isfinite(value) or (minimum is not None and
                (value < minimum or (strict and value == minimum)))):
            raise ValueError(f'{name} is outside its finite supported range.')
        result[name] = value

    def count(name, minimum):
        value = getattr(config, name)
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
            raise ValueError(f'{name} must be an integer >= {minimum}.')
        result[name] = int(value)

    count('Opt_num_of_bands', 1)
    number('Opt_eta', 0., strict=config.Engine == 'QE')
    number('Opt_FD_smearing', 0., strict=config.Engine == 'QE')
    if config.Engine == 'QE' or method == 'BSE':
        number('Opt_min_en', 0.)
        number('Opt_max_en', 0., strict=True)
        count('Opt_num_of_data', 2)
        if result['Opt_max_en'] <= result['Opt_min_en']:
            raise ValueError('Opt_max_en must be greater than Opt_min_en.')
    if config.Engine == 'QE':
        number('Opt_shift_en')
    else:
        number('Opt_cut_of_energy', 0., strict=True)
        if method == 'RPA':
            number('Opt_domega0', 0., strict=True)
            number('Opt_omega2', 0., strict=True)
    return result


def validate_optical_table(data):
    """Require finite seven-column spectra on an increasing nonnegative grid."""
    try:
        table = np.asarray(data, dtype=float)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('Optical tables require finite N x 7 arrays with at least two energy points.') from exc
    if (table.ndim != 2 or table.shape[1] != 7 or len(table) < 2
            or not np.isfinite(table).all() or np.any(table[:, 0] < 0)
            or np.any(np.diff(table[:, 0]) <= 0)):
        raise ValueError('Optical tables require finite N x 7 arrays on an increasing nonnegative energy grid.')
    return table.copy()


def run_optical_exports(callback):
    """Export on MPI root and propagate output failures to every rank."""
    from ase.parallel import broadcast, world
    error = result = None
    if world.rank == 0:
        try:
            result = callback()
        except (Exception, KeyboardInterrupt) as exc:
            error = f'{type(exc).__name__}: {exc}'
    error, result = broadcast((error, result), root=0, comm=world)
    if error is not None:
        raise RuntimeError('Optical export failed: ' + error)
    return result


REDUCED_PLANCK_EV_SECONDS = 6.582119569e-16
SPEED_OF_LIGHT_CM_PER_SECOND = 2.99792458e10


def derive_optical_table(energies, epsilon_real, epsilon_imaginary):
    """Derive seven-column optical spectra using the existing positive-k convention."""
    try:
        energy, real, imaginary = [np.asarray(values, dtype=float)
                                   for values in (energies, epsilon_real, epsilon_imaginary)]
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError('Dielectric spectra require matching finite one-dimensional arrays.') from exc
    if (energy.ndim != 1 or len(energy) < 2 or real.shape != energy.shape
            or imaginary.shape != energy.shape
            or not all(np.isfinite(values).all() for values in (energy, real, imaginary))):
        raise ValueError('Dielectric spectra require matching finite one-dimensional arrays.')
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            # Complex square roots avoid cancellation in |epsilon| - epsilon_real
            # for weak absorption. Retain the historical nonnegative extinction.
            index = np.sqrt(real + 1j * imaginary)
            refractive, extinction = index.real, np.abs(index.imag)
            reflection = np.empty_like(energy)
            regular = np.maximum(refractive, extinction) < 1e150
            n, k = refractive[regular], extinction[regular]
            reflection[regular] = ((n - 1.) ** 2 + k ** 2) / ((n + 1.) ** 2 + k ** 2)
            n, k = refractive[~regular], extinction[~regular]
            denominator = np.hypot(n + 1., k)
            reflection[~regular] = ((n - 1.) / denominator) ** 2 + (k / denominator) ** 2
            absorption = (2. * energy / (REDUCED_PLANCK_EV_SECONDS
                          * SPEED_OF_LIGHT_CM_PER_SECOND)) * extinction
            table = np.column_stack([energy, real, imaginary, refractive, extinction, absorption, reflection])
    except FloatingPointError as exc:
        raise ValueError('Derived optical properties exceed the supported finite numerical range.') from exc
    return validate_optical_table(table)


OPTICAL_TABLE_HEADER = ('Energy(eV) Eps_real Eps_img Refractive_Index '
                        'Extinction_Index Absorption(1/cm) Reflectivity')


def write_optical_table(path, data):
    """Publish a validated optical table without truncating a previous result."""
    from pathlib import Path
    table = validate_optical_table(data)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + '.tmp')
    try:
        with temporary.open('w', encoding='utf-8') as stream:
            np.savetxt(stream, table, fmt='%.12g', header=OPTICAL_TABLE_HEADER, comments='')
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)
    return path
