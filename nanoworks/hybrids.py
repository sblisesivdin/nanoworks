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
