"""Engine-neutral hybrid-functional controls."""

import math


def validate_exx_kpoint_density(value):
    """Validate an optional exact-exchange q-point density."""
    if value is None:
        return None

    if isinstance(value, bool):
        raise TypeError(
            'EXX_kpoint_density must be a numeric value.'
        )

    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            'EXX_kpoint_density must be a numeric value.'
        ) from exc

    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(
            'EXX_kpoint_density must be finite and greater than zero.'
        )

    return value


def validate_exx_cutoff(value, wavefunction_cutoff):
    """Validate an optional exact-exchange cutoff in eV."""
    if value is None:
        return None

    if isinstance(value, bool):
        raise TypeError('EXX_cutoff must be a numeric value in eV.')

    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            'EXX_cutoff must be a numeric value in eV.'
        ) from exc

    if not math.isfinite(value):
        raise ValueError('EXX_cutoff must be finite.')

    wavefunction_cutoff = float(wavefunction_cutoff)
    if value <= wavefunction_cutoff:
        raise ValueError(
            'EXX_cutoff must be greater than Wavefunction_cutoff.'
        )

    return value
