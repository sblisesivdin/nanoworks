"""Engine-neutral density-of-states integration settings."""

import math


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


def validate_dos_settings(integration='smearing', width=0.1):
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
