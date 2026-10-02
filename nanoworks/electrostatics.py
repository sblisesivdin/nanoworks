"""Engine-neutral electrostatic boundary settings."""


ELECTROSTATIC_BOUNDARIES = (
    'periodic',
    'isolated-2d',
)

ELECTROSTATIC_NORMAL_AXES = (
    'x',
    'y',
    'z',
)

_BOUNDARY_ALIASES = {
    'periodic': 'periodic',
    'pbc': 'periodic',
    '3d': 'periodic',
    'isolated-2d': 'isolated-2d',
    'isolated_2d': 'isolated-2d',
    '2d': 'isolated-2d',
}


def validate_electrostatic_settings(
    boundary='periodic',
    normal_axis='z',
    dipole_correction=False,
):
    """Normalize portable electrostatic boundary and dipole intent.

    ``isolated-2d`` describes Coulomb isolation normal to a periodic
    two-dimensional plane.  Dipole correction is kept as a separate
    setting because it is not physically equivalent to Coulomb truncation.
    """
    boundary_key = str(boundary).strip().lower()
    try:
        boundary = _BOUNDARY_ALIASES[boundary_key]
    except KeyError:
        raise ValueError(
            'Electrostatic_boundary must be one of: '
            + ', '.join(ELECTROSTATIC_BOUNDARIES)
            + '.'
        )

    normal_axis = str(normal_axis).strip().lower()
    if normal_axis not in ELECTROSTATIC_NORMAL_AXES:
        raise ValueError(
            'Electrostatic_normal_axis must be one of: '
            + ', '.join(ELECTROSTATIC_NORMAL_AXES)
            + '.'
        )

    if not isinstance(dipole_correction, bool):
        raise TypeError('Dipole_correction must be a boolean value.')

    periodic_axes = [True, True, True]
    if boundary == 'isolated-2d':
        periodic_axes[ELECTROSTATIC_NORMAL_AXES.index(normal_axis)] = False

    return {
        'boundary': boundary,
        'normal_axis': normal_axis,
        'dipole_correction': dipole_correction,
        'periodic_axes': tuple(periodic_axes),
    }
