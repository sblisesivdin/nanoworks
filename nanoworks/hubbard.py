# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Engine-neutral Hubbard-U setting helpers."""

import math
import re


HUBBARD_MANIFOLD = re.compile(
    r'([A-Z][a-z]?)-(\d+)([spdfg])',
)


def normalize_hubbard_u(hubbard_u):
    """Validate and normalize a ``manifold -> U (eV)`` mapping."""
    if hubbard_u is None:
        return {}

    if not isinstance(hubbard_u, dict):
        raise TypeError(
            'Hubbard_U must be a manifold-to-energy dictionary.'
        )

    normalized = {}

    for raw_manifold, raw_value in hubbard_u.items():
        if not isinstance(raw_manifold, str):
            raise TypeError(
                'Hubbard_U manifold names must be strings.'
            )

        manifold = raw_manifold.strip()
        match = HUBBARD_MANIFOLD.fullmatch(manifold)

        if match is None:
            raise ValueError(
                'Hubbard_U keys must use an explicit '
                "'Element-nl' manifold such as 'O-2p' or 'Zn-3d': "
                f'{raw_manifold!r}'
            )

        if isinstance(raw_value, bool):
            raise TypeError(
                'Hubbard_U energies must be numeric values in eV.'
            )

        try:
            value_ev = float(raw_value)
        except (TypeError, ValueError) as exc:
            raise TypeError(
                'Hubbard_U energies must be numeric values in eV.'
            ) from exc

        if not math.isfinite(value_ev) or value_ev <= 0.0:
            raise ValueError(
                'Hubbard_U energies must be finite and greater than zero.'
            )

        symbol, principal_number, orbital = match.groups()
        canonical = f'{symbol}-{int(principal_number)}{orbital}'

        if canonical in normalized:
            raise ValueError(
                f'Duplicate Hubbard_U manifold: {canonical}'
            )

        normalized[canonical] = value_ev

    return normalized


def resolve_gpaw_hubbard(hubbard_u):
    """Translate portable Hubbard manifolds to GPAW setup strings."""
    normalized = normalize_hubbard_u(hubbard_u)
    setups = {}

    for manifold, value_ev in normalized.items():
        symbol, orbital = manifold.split('-', 1)

        if symbol in setups:
            raise ValueError(
                'GPAW supports one Hubbard_U correction per element: '
                f'{symbol}'
            )

        setups[symbol] = f':{orbital[-1]},{value_ev}'

    return setups
