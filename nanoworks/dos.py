# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Engine-neutral density-of-states integration settings."""

import math

import numpy as np


def ensure_gpaw_soc_dos_weights(rawdos, nbzkpts):
    """Match globally collected SOC eigenvalues to full-BZ uniform weights.

    GPAW 26.7 BZWaveFunctions.eigenvalues() collects every BZ point, but
    weights() uses the number of rank-local wavefunctions. Gaussian DOS
    then silently truncates its sum through zip(eigenvalues, weights).
    Repair only that known layout; reject other inconsistent data rather
    than renormalizing a partial spectrum. SOC spinors count once each.
    """
    eigenvalues = np.asarray(rawdos.eig_skn)
    weights = np.asarray(rawdos.weight_k, dtype=float)
    if (nbzkpts < 1 or eigenvalues.ndim != 3
            or eigenvalues.shape[0] != 1 or eigenvalues.shape[1] != nbzkpts
            or not np.isfinite(eigenvalues).all()):
        raise ValueError('GPAW SOC DOS requires finite eigenvalues for the full BZ.')
    if (weights.ndim != 1 or len(weights) > nbzkpts
            or not np.isfinite(weights).all()
            or not np.allclose(weights, 1.0 / nbzkpts, rtol=1e-10, atol=1e-14)):
        raise ValueError('GPAW SOC DOS has unexpected full-BZ k-point weights.')
    if len(weights) == nbzkpts:
        return False
    rawdos.weight_k = np.full(nbzkpts, 1.0 / nbzkpts)
    return True


DOS_INTEGRATION_METHODS = (
    'smearing',
    'tetrahedron',
)

_ALIASES = {
    'gaussian': 'smearing',
    'smearing': 'smearing',
    'tetrahedron': 'tetrahedron',
    'tetrahedra': 'tetrahedron',
}


def validate_dos_settings(integration='tetrahedron', width=0.1):
    """Normalize portable DOS integration settings."""
    key = str(integration).strip().lower()
    try:
        integration = _ALIASES[key]
    except KeyError as exc:
        raise ValueError(
            'DOS_integration must be one of: '
            + ', '.join(DOS_INTEGRATION_METHODS)
            + '.'
        ) from exc

    if integration == 'tetrahedron':
        return {
            'integration': integration,
            'width_ev': 0.0,
        }

    if isinstance(width, bool):
        raise ValueError(
            'DOS_width must be finite and greater than zero for smearing.'
        )
    try:
        width = float(width)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            'DOS_width must be finite and greater than zero for smearing.'
        ) from exc
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError(
            'DOS_width must be finite and greater than zero for smearing.'
        )

    return {
        'integration': integration,
        'width_ev': width,
    }


def resolve_dos_settings(engine, integration, width, ground_occupation):
    """Translate portable DOS integration intent for one DFT engine."""
    resolved = validate_dos_settings(integration, width)
    engine = str(engine).strip().upper()

    settings = {
        **resolved,
        'electronic_occupation': ground_occupation,
        'bz_sum': None,
        'degauss_ev': None,
        'ngauss': None,
    }

    if engine == 'GPAW':
        return settings

    if engine == 'QE':
        if resolved['integration'] == 'tetrahedron':
            settings['electronic_occupation'] = 'tetrahedra'
            settings['bz_sum'] = 'tetrahedra'
        else:
            settings['bz_sum'] = 'smearing'
            settings['degauss_ev'] = resolved['width_ev']
            settings['ngauss'] = 0
        return settings

    raise ValueError(f'Unsupported DFT engine: {engine}')
