# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Engine-neutral electronic occupation settings and backend adapters."""

import math


OCCUPATION_SCHEMES = (
    'fixed',
    'fermi-dirac',
    'methfessel-paxton',
    'marzari-vanderbilt',
)

_ALIASES = {
    'fixed': 'fixed',
    'fixed-uniform': 'fixed',
    'fd': 'fermi-dirac',
    'fermi-dirac': 'fermi-dirac',
    'fermi_dirac': 'fermi-dirac',
    'mp': 'methfessel-paxton',
    'methfessel-paxton': 'methfessel-paxton',
    'methfessel_paxton': 'methfessel-paxton',
    'cold': 'marzari-vanderbilt',
    'marzari-vanderbilt': 'marzari-vanderbilt',
    'marzari_vanderbilt': 'marzari-vanderbilt',
}


def validate_occupation_settings(scheme='fermi-dirac', width=0.05):
    """Normalize the portable occupation scheme and smearing width."""
    key = str(scheme).strip().lower()
    try:
        scheme = _ALIASES[key]
    except KeyError:
        raise ValueError(
            'Occupation_scheme must be one of: '
            + ', '.join(OCCUPATION_SCHEMES)
            + '.'
        )

    if scheme == 'fixed':
        return {
            'scheme': scheme,
            'width': None,
        }

    if isinstance(width, bool):
        raise ValueError(
            'Smearing_width must be finite and greater than zero.'
        )
    try:
        width = float(width)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            'Smearing_width must be finite and greater than zero.'
        ) from exc
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError(
            'Smearing_width must be finite and greater than zero.'
        )

    return {
        'scheme': scheme,
        'width': width,
    }


def resolve_gpaw_occupation(**settings):
    """Translate portable occupation intent to GPAW settings."""
    resolved = validate_occupation_settings(**settings)
    if resolved['scheme'] == 'fixed':
        return {'name': 'fixed-uniform'}
    return {
        'name': resolved['scheme'],
        'width': resolved['width'],
    }


def resolve_qe_occupation_input(**settings):
    """Translate portable occupation intent to QE-adapter input."""
    resolved = validate_occupation_settings(**settings)
    if resolved['scheme'] == 'fixed':
        return 'fixed'
    return {
        'name': resolved['scheme'],
        'width': resolved['width'],
    }


def resolve_engine_occupation(engine, **settings):
    """Resolve portable occupation intent for the selected DFT engine."""
    engine = str(engine).strip().upper()
    if engine == 'GPAW':
        return resolve_gpaw_occupation(**settings)
    if engine == 'QE':
        return resolve_qe_occupation_input(**settings)
    raise ValueError(f'Unsupported DFT engine: {engine}')
