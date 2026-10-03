"""Engine-neutral helpers for elastic-property reporting."""

import math

import numpy as np


_AXIS_INDEX = {
    'x': 0,
    'y': 1,
    'z': 2,
}

_IN_PLANE_VOIGT = {
    'x': (1, 2, 3),
    'y': (0, 2, 4),
    'z': (0, 1, 5),
}


def normalize_elastic_dimensionality(value):
    """Return the supported ``auto``, ``2D``, or ``3D`` spelling."""
    normalized = str(value).strip().lower()

    if normalized == 'auto':
        return 'auto'
    if normalized == '2d':
        return '2D'
    if normalized == '3d':
        return '3D'

    raise ValueError(
        "Elastic_dimensionality must be 'auto', '2D', or '3D'."
    )


def normalize_elastic_normal_axis(value):
    """Validate and normalize the non-periodic axis used for 2D results."""
    axis = str(value).strip().lower()

    if axis not in _AXIS_INDEX:
        raise ValueError(
            "Elastic_normal_axis must be 'x', 'y', or 'z'."
        )

    return axis


def estimate_vacuum_gap(atoms, normal_axis='z'):
    """Estimate the largest periodic nuclear gap along one cell vector."""
    axis = normalize_elastic_normal_axis(normal_axis)
    axis_index = _AXIS_INDEX[axis]
    cell = np.asarray(atoms.cell, dtype=float)
    cell_length = float(np.linalg.norm(cell[axis_index]))
    other_lengths = [
        float(np.linalg.norm(cell[index]))
        for index in range(3)
        if index != axis_index
    ]

    if not math.isfinite(cell_length) or cell_length <= 0.0:
        raise ValueError(
            f'Elastic {axis}-axis cell length must be positive and finite.'
        )

    if any(
        not math.isfinite(length) or length <= 0.0
        for length in other_lengths
    ):
        raise ValueError(
            'Elastic cell-vector lengths must be positive and finite.'
        )

    scaled = np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float)

    if scaled.ndim != 2 or scaled.shape[0] == 0 or scaled.shape[1] != 3:
        raise ValueError('Elastic dimensionality requires at least one atom.')

    coordinates = np.sort(np.mod(scaled[:, axis_index], 1.0))
    periodic_gaps = np.diff(
        np.concatenate((coordinates, coordinates[:1] + 1.0))
    )
    vacuum_fraction = float(np.max(periodic_gaps))

    return {
        'normal_axis': axis,
        'cell_length_angstrom': cell_length,
        'vacuum_gap_angstrom': vacuum_fraction * cell_length,
        'vacuum_fraction': vacuum_fraction,
        'cell_aspect_ratio': cell_length / max(other_lengths),
    }


def resolve_elastic_dimensionality(
    atoms,
    dimensionality='auto',
    normal_axis='z',
    minimum_vacuum_angstrom=5.0,
    minimum_vacuum_fraction=0.30,
    minimum_cell_aspect_ratio=1.50,
):
    """Resolve whether elastic results should be reported as 2D or 3D."""
    requested = normalize_elastic_dimensionality(dimensionality)
    geometry = estimate_vacuum_gap(atoms, normal_axis=normal_axis)

    if requested == 'auto':
        is_2d = (
            geometry['vacuum_gap_angstrom'] >= minimum_vacuum_angstrom
            and geometry['vacuum_fraction'] >= minimum_vacuum_fraction
            and geometry['cell_aspect_ratio'] >= minimum_cell_aspect_ratio
        )
        resolved = '2D' if is_2d else '3D'
    else:
        resolved = requested

    return {
        'requested': requested,
        'resolved': resolved,
        'automatic': requested == 'auto',
        **geometry,
    }


def calculate_2d_elastic_properties(
    elastic_tensor_gpa,
    normal_cell_length_angstrom,
    normal_axis='z',
):
    """Convert supercell-normalized stiffness to intrinsic 2D properties."""
    axis = normalize_elastic_normal_axis(normal_axis)
    tensor = np.asarray(elastic_tensor_gpa, dtype=float)

    if tensor.shape != (6, 6) or not np.all(np.isfinite(tensor)):
        raise ValueError(
            'The elastic tensor must be a finite 6x6 matrix.'
        )

    conversion_length = float(normal_cell_length_angstrom)

    if not math.isfinite(conversion_length) or conversion_length <= 0.0:
        raise ValueError(
            'The 2D elastic conversion length must be positive and finite.'
        )

    indices = _IN_PLANE_VOIGT[axis]
    stiffness = (
        tensor[np.ix_(indices, indices)]
        * conversion_length
        * 0.1
    )

    try:
        compliance = np.linalg.inv(stiffness)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            'The in-plane elastic stiffness matrix is singular.'
        ) from exc

    if (
        not np.all(np.isfinite(compliance))
        or any(
            abs(compliance[index, index]) < 1.0e-15
            for index in range(3)
        )
    ):
        raise ValueError(
            'The in-plane elastic compliance matrix is not finite.'
        )

    young_first = 1.0 / compliance[0, 0]
    young_second = 1.0 / compliance[1, 1]
    shear = 1.0 / compliance[2, 2]

    return {
        'normal_axis': axis,
        'normal_cell_length_angstrom': conversion_length,
        'voigt_indices': indices,
        'stiffness_n_per_m': stiffness,
        'young_modulus_first_n_per_m': young_first,
        'young_modulus_second_n_per_m': young_second,
        'shear_modulus_n_per_m': shear,
        'poisson_ratio_first_second': -compliance[0, 1] / compliance[0, 0],
        'poisson_ratio_second_first': -compliance[0, 1] / compliance[1, 1],
    }
