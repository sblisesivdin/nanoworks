# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""Engine-neutral hybrid-functional controls."""

import math


HYBRID_ALIASES = {
    'HSE': 'HSE06',
    'HSE06': 'HSE06',
    'HSE03': 'HSE03',
    'PBE0': 'PBE0',
    'B3LYP': 'B3LYP',
    'EXX': 'EXX',
}

HYBRID_DEFAULTS = {
    'HSE06': {
        'exx_fraction': 0.25,
        'omega': 0.106,
        'screened': True,
    },
    'HSE03': {
        'exx_fraction': 0.25,
        'omega': 0.15,
        'screened': True,
    },
    'PBE0': {
        'exx_fraction': 0.25,
        'omega': None,
        'screened': False,
    },
    'B3LYP': {
        'exx_fraction': None,
        'omega': None,
        'screened': False,
    },
    'EXX': {
        'exx_fraction': None,
        'omega': None,
        'screened': False,
    },
}

HYBRID_CAPABILITIES = {
    'GPAW': frozenset(HYBRID_DEFAULTS),
    'QE': frozenset({'HSE06', 'HSE03', 'PBE0'}),
}

HYBRID_STAGE_CAPABILITIES = {
    'GPAW': frozenset({
        'ground',
        'geometry',
        'elastic',
        'dos',
        'band',
        'density',
        'optical',
    }),
    'QE': frozenset({
        'ground',
        'dos',
        'band',
        'density',
    }),
}

HYBRID_STAGES = frozenset({
    'ground',
    'geometry',
    'cell-relaxation',
    'elastic',
    'dos',
    'band',
    'density',
    'phonon',
    'optical',
})


def _extract_xc_name(xc_calc):
    """Extract an XC name from a string or calculator dictionary."""
    if isinstance(xc_calc, dict):
        xc_calc = xc_calc.get('name', xc_calc.get('xc'))
    if xc_calc is None:
        return ''
    return str(xc_calc).strip()


def normalize_hybrid_name(xc_calc):
    """Return the canonical hybrid name, or ``None`` when non-hybrid."""
    name = (
        _extract_xc_name(xc_calc)
        .upper()
        .replace('_', '')
        .replace('-', '')
    )
    return HYBRID_ALIASES.get(name)


def is_hybrid_functional(xc_calc):
    """Return whether an XC specification selects a known hybrid."""
    return normalize_hybrid_name(xc_calc) is not None


def _validate_optional_real(value, keyword, *, maximum=None):
    """Validate one optional finite positive hybrid control."""
    if value is None:
        return None
    if isinstance(value, bool):
        raise TypeError(f'{keyword} must be a numeric value.')
    try:
        value = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f'{keyword} must be a numeric value.') from exc
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f'{keyword} must be finite and greater than zero.')
    if maximum is not None and value > maximum:
        raise ValueError(f'{keyword} must be less than or equal to {maximum}.')
    return value


def resolve_hybrid_settings(
    xc_calc,
    exx_fraction=None,
    omega=None,
    engine=None,
):
    """Normalize portable hybrid intent and validate backend capability."""
    name = normalize_hybrid_name(xc_calc)
    if name is None:
        if exx_fraction is not None or omega is not None:
            raise ValueError(
                'XC_exx_fraction and XC_omega require a hybrid '
                'XC_calc.'
            )
        return None

    if engine is not None:
        engine = str(engine).strip().upper()
        try:
            supported = HYBRID_CAPABILITIES[engine]
        except KeyError:
            raise ValueError(f'Unsupported DFT engine: {engine}')
        if name not in supported:
            supported_names = ', '.join(sorted(supported))
            raise NotImplementedError(
                f'{engine} does not support hybrid functional {name}. '
                f'Supported hybrids: {supported_names}.'
            )

    exx_fraction = _validate_optional_real(
        exx_fraction,
        'XC_exx_fraction',
        maximum=1.0,
    )
    omega = _validate_optional_real(omega, 'XC_omega')
    defaults = HYBRID_DEFAULTS[name]
    if omega is not None and not defaults['screened']:
        raise ValueError(
            'XC_omega is only valid for screened HSE functionals.'
        )

    return {
        'name': name,
        'screened': defaults['screened'],
        'default_exx_fraction': defaults['exx_fraction'],
        'default_omega': defaults['omega'],
        'exx_fraction': exx_fraction,
        'omega': omega,
    }


def get_unsupported_hybrid_stages(xc_calc, engine, stages):
    """Return requested stages unavailable for a backend hybrid workflow."""
    settings = resolve_hybrid_settings(xc_calc, engine=engine)
    if settings is None:
        return ()

    engine = str(engine).strip().upper()
    normalized_stages = []
    for stage in stages:
        normalized = str(stage).strip().lower().replace('_', '-')
        if normalized not in HYBRID_STAGES:
            raise ValueError(f'Unknown hybrid calculation stage: {stage}')
        if normalized not in normalized_stages:
            normalized_stages.append(normalized)

    supported = HYBRID_STAGE_CAPABILITIES[engine]
    return tuple(
        stage
        for stage in normalized_stages
        if stage not in supported
    )


def validate_hybrid_stage_support(xc_calc, engine, stages):
    """Validate stage support and return unsupported-stage errors early."""
    unsupported = get_unsupported_hybrid_stages(
        xc_calc,
        engine,
        stages,
    )
    if unsupported:
        engine = str(engine).strip().upper()
        name = normalize_hybrid_name(xc_calc)
        stage_names = ', '.join(unsupported)
        raise NotImplementedError(
            f'{engine} hybrid functional {name} does not support '
            f'the following calculation stage(s): {stage_names}.'
        )
    return normalize_hybrid_name(xc_calc)


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
