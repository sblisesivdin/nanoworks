"""GPAW computation engine helpers."""

from gpaw import GPAW, PW, MixerSum
from gpaw.eigensolvers import Davidson

HYBRID_XC = ('HSE06', 'HSE03', 'B3LYP', 'PBE0', 'EXX')


def _get_xc_name(xc):
    """Extract the functional name from a GPAW XC specification."""
    if isinstance(xc, dict):
        xc = xc.get('name', xc.get('xc'))

    if xc is None:
        return ''

    return str(xc).strip().upper()


def is_hybrid(xc_calc):
    """Return whether the requested XC functional uses the hybrid workflow."""
    return _get_xc_name(xc_calc) in HYBRID_XC

def build_hybrid_xc(
    xc_calc,
    exx_fraction=None,
    omega=None,
    backend='pw',
):
    """Build the GPAW dictionary specification for a hybrid functional."""
    xc = {
        'name': str(xc_calc).upper(),
        'backend': backend,
    }

    if exx_fraction is not None:
        xc['fraction'] = exx_fraction

    if omega is not None:
        xc['omega'] = omega

    return xc

def create_gpaw_calc(*args, **kwargs):
    """Create a GPAW calculator using legacy GPAW for hybrid functionals."""
    if is_hybrid(kwargs.get('xc')) and 'legacy_gpaw' not in kwargs:
        kwargs['legacy_gpaw'] = True

    return GPAW(*args, **kwargs)

def load_gpaw_calc(filename, hybrid=False, **kwargs):
    """Load a GPAW calculator from a state file."""
    if hybrid and 'legacy_gpaw' not in kwargs:
        kwargs['legacy_gpaw'] = True

    return GPAW(filename, **kwargs)

def resolve_xc_and_setups(xc_input, user_setups=None):
    """Resolve GPAW XC and setup specifications."""

    if user_setups is None:
        setups = {}
    else:
        setups = dict(user_setups)

    is_libxc = False

    if isinstance(xc_input, dict):
        xc_str = str(
            xc_input.get('name', xc_input.get('xc', ''))
        ).strip()
    else:
        xc_str = str(xc_input).strip()

    if xc_str.lower().startswith('libxc:'):
        xc_str = xc_str[6:].strip()
        is_libxc = True

    elif any(
        xc_str.startswith(prefix)
        for prefix in ('MGGA_', 'GGA_', 'HYB_', 'LDA_')
    ):
        is_libxc = True

    elif '+' in xc_str:
        is_libxc = True

    return xc_str, setups, is_libxc

def build_kpoint_spec(density, size, gamma):
    """Build a GPAW k-point specification."""
    if density is not None:
        return {
            'density': density,
            'gamma': gamma,
        }

    return {
        'size': tuple(size),
        'gamma': gamma,
    }
    
def build_grid_spec(spacing, size):
    """Build a GPAW real-space grid specification."""
    if spacing is not None:
        return {
            'h': spacing,
        }

    return {
        'gpts': tuple(size),
    }

def build_ground_common_kwargs(
    mixer,
    charge,
    spinpol,
    txt,
    convergence,
    occupations,
    nbands='200%',
):
    """Build calculator arguments shared by GPAW ground-state modes."""
    return {
        'nbands': nbands,
        'mixer': mixer,
        'charge': charge,
        'spinpol': spinpol,
        'txt': txt,
        'convergence': convergence,
        'occupations': occupations,
    }

def create_regular_pw_ground_calc(
    cutoff,
    xc,
    setups,
    parallel,
    mixer,
    charge,
    spinpol,
    txt,
    convergence,
    occupations,
    kpoint_density,
    kpoint_size,
    gamma,
):
    """Create a regular GPAW plane-wave ground-state calculator."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
    )

    kwargs.update({
        'mode': PW(
            ecut=cutoff,
            force_complex_dtype=True,
        ),
        'xc': xc,
        'setups': setups,
        'parallel': parallel,
        'kpts': build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    })

    return create_gpaw_calc(**kwargs)

def create_hybrid_pw_ground_calc(
    cutoff,
    xc_calc,
    exx_fraction,
    omega,
    backend,
    mixer,
    charge,
    spinpol,
    txt,
    convergence,
    occupations,
    kpoint_density,
    kpoint_size,
    gamma,
):
    """Create a hybrid GPAW plane-wave ground-state calculator."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
    )

    kwargs.update({
        'mode': PW(
            ecut=cutoff,
            force_complex_dtype=True,
        ),
        'xc': build_hybrid_xc(
            xc_calc,
            exx_fraction,
            omega,
            backend,
        ),
        'parallel': {
            'band': 1,
            'kpt': 1,
        },
        'eigensolver': Davidson(niter=1),
        'kpts': build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    })

    return create_gpaw_calc(**kwargs)

def create_lcao_ground_calc(
    setups,
    parallel,
    mixer,
    charge,
    spinpol,
    txt,
    convergence,
    occupations,
    kpoint_density,
    kpoint_size,
    gamma,
    grid_spacing,
    grid_size,
    basis='dzp',
):
    """Create a GPAW LCAO ground-state calculator."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
    )

    kwargs.update({
        'mode': 'lcao',
        'basis': basis,
        'setups': setups,
        'parallel': parallel,
        'kpts': build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    })

    kwargs.update(
        build_grid_spec(
            spacing=grid_spacing,
            size=grid_size,
        )
    )

    return create_gpaw_calc(**kwargs)

def create_elastic_calc(
    cutoff,
    xc,
    setups,
    parallel,
    spinpol,
    kpoint_size,
    gamma,
    mixer,
    txt,
    charge,
    convergence,
    occupations,
    hybrid=False,
):
    """Create a GPAW calculator for elastic deformations."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
    )

    kwargs.update({
        'mode': PW(
            ecut=cutoff,
            force_complex_dtype=True,
        ),
        'xc': xc,
        'setups': setups,
        'parallel': parallel,
        'kpts': {
            'size': tuple(kpoint_size),
            'gamma': gamma,
        },
    })

    if hybrid:
        kwargs['eigensolver'] = Davidson(niter=1)

    return create_gpaw_calc(**kwargs)

def create_phonon_calc(
    cutoff,
    kpoint_size,
    txt,
):
    """Create a GPAW calculator for finite-displacement phonons."""
    return create_gpaw_calc(
        mode=PW(cutoff),
        kpts={
            'size': tuple(kpoint_size),
        },
        txt=txt,
    )

def resolve_elastic_settings(
    xc_calc,
    setups,
    world_size,
    exx_fraction=None,
    omega=None,
    backend='pw',
):
    """Resolve XC, setups, parallel settings, and hybrid state for elasticity."""
    actual_xc, resolved_setups, _ = resolve_xc_and_setups(
        xc_calc,
        setups,
    )

    hybrid = is_hybrid(xc_calc)

    if hybrid:
        elastic_xc = build_hybrid_xc(
            xc_calc,
            exx_fraction,
            omega,
            backend,
        )
        parallel = {
            'band': 1,
            'kpt': 1,
        }
    else:
        elastic_xc = actual_xc
        parallel = {
            'domain': world_size,
        }

    return elastic_xc, resolved_setups, parallel, hybrid

def create_default_mixer():
    """Create the default GPAW density mixer used by Nanoworks."""
    return MixerSum(
        beta=0.1,
        nmaxold=3,
        weight=50,
    )
