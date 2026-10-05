# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Engine-neutral plane-wave cutoff helpers."""

import math


def validate_cutoff_settings(
    wavefunction_cutoff,
    density_cutoff_ratio=4.0,
):
    """Validate plane-wave and charge-density cutoff settings."""
    if isinstance(wavefunction_cutoff, bool):
        raise TypeError(
            'Wavefunction_cutoff must be a numeric value in eV.'
        )

    try:
        wavefunction_cutoff = float(wavefunction_cutoff)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            'Wavefunction_cutoff must be a numeric value in eV.'
        ) from exc

    if (
        not math.isfinite(wavefunction_cutoff)
        or wavefunction_cutoff <= 0.0
    ):
        raise ValueError(
            'Wavefunction_cutoff must be finite and greater than zero.'
        )

    if isinstance(density_cutoff_ratio, bool):
        raise TypeError(
            'Density_cutoff_ratio must be a numeric value.'
        )

    try:
        density_cutoff_ratio = float(density_cutoff_ratio)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            'Density_cutoff_ratio must be a numeric value.'
        ) from exc

    if (
        not math.isfinite(density_cutoff_ratio)
        or density_cutoff_ratio < 1.0
    ):
        raise ValueError(
            'Density_cutoff_ratio must be finite and at least 1.0.'
        )

    return {
        'wavefunction_ev': wavefunction_cutoff,
        'density_ratio': density_cutoff_ratio,
        'density_ev': wavefunction_cutoff * density_cutoff_ratio,
    }
