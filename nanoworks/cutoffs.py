# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Engine-neutral plane-wave cutoff helpers."""

import math
from pathlib import Path

EV_PER_RYDBERG = 13.605693122994


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

    # None defers the density choice until the QE UPFs are available.
    if density_cutoff_ratio is None:
        return {
            'wavefunction_ev': wavefunction_cutoff,
            'density_ratio': None,
            'density_ev': None,
        }

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


def resolve_qe_cutoff_settings(
    wavefunction_cutoff, pseudopotentials, pseudo_dir=None,
    density_cutoff_ratio=None,
):
    """Choose a QE density cutoff from UPF hints unless explicitly set.

    Norm-conserving UPFs retain the fourfold baseline. Positive header
    recommendations supply both an absolute density floor and, when both
    cutoffs exist, a ratio maintained during wavefunction-cutoff sweeps.
    No recommendation is a convergence guarantee. Missing information is
    reported through notices, without changing the wavefunction cutoff.
    """
    from nanoworks.pseudos import read_upf_cutoff_metadata

    settings = validate_cutoff_settings(wavefunction_cutoff, density_cutoff_ratio)
    settings['notices'] = []
    if density_cutoff_ratio is not None:
        settings['source'] = 'explicit Density_cutoff_ratio'
        return settings

    ratio = 4.0
    wavefunction_ry = settings['wavefunction_ev'] / EV_PER_RYDBERG
    for symbol, filename in sorted(pseudopotentials.items()):
        path = Path(filename).expanduser()
        if not path.is_absolute() and pseudo_dir is not None:
            path = Path(pseudo_dir).expanduser() / path
        try:
            metadata = read_upf_cutoff_metadata(path)
        except OSError:
            metadata = {'wfc_cutoff': None, 'rho_cutoff': None, 'pseudo_type': ''}
        wfc = metadata['wfc_cutoff']
        rho = metadata['rho_cutoff']
        if wfc is not None and wavefunction_ry < wfc:
            settings['notices'].append(
                f'{symbol}: Wavefunction_cutoff={settings["wavefunction_ev"]:g} eV '
                f'is below the UPF recommendation ({wfc * EV_PER_RYDBERG:g} eV); '
                'the supplied value is retained. Check cutoff convergence.'
            )
        if rho is not None:
            ratio = max(ratio, rho / wavefunction_ry)
            if wfc is not None:
                ratio = max(ratio, rho / wfc)
        elif metadata['pseudo_type'] not in ('NC', 'SL', 'NORM-CONSERVING'):
            settings['notices'].append(
                f'{symbol}: no usable UPF density-cutoff recommendation in {path} '
                f'(type {metadata["pseudo_type"] or "unknown"}); using the fourfold '
                'baseline unless another species requires more. Set '
                'Density_cutoff_ratio explicitly and check convergence for '
                'ultrasoft/PAW pseudopotentials.'
            )
    settings.update(validate_cutoff_settings(wavefunction_cutoff, ratio))
    settings['source'] = 'automatic UPF policy'
    return settings
