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
