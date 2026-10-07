# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

"""GPAW computation engine helpers."""

from gpaw import GPAW, PW, MixerSum, FermiDirac
from gpaw.eigensolvers import Davidson
from numbers import Integral
from nanoworks.hubbard import resolve_gpaw_hubbard
from nanoworks.hybrids import (
    is_hybrid_functional,
    resolve_hybrid_settings,
)


def is_hybrid(xc_calc):
    """Return whether the requested XC functional uses the hybrid workflow."""
    return is_hybrid_functional(xc_calc)

def build_hybrid_xc(
    xc_calc,
    exx_fraction=None,
    omega=None,
    backend='pw',
):
    """Build the GPAW dictionary specification for a hybrid functional."""
    settings = resolve_hybrid_settings(
        xc_calc,
        exx_fraction=exx_fraction,
        omega=omega,
        engine='GPAW',
    )
    if settings is None:
        raise ValueError('build_hybrid_xc requires a hybrid XC functional.')

    xc = {
        'name': settings['name'],
        'backend': backend,
    }

    if settings['exx_fraction'] is not None:
        xc['fraction'] = settings['exx_fraction']

    if settings['omega'] is not None:
        xc['omega'] = settings['omega']

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

def resolve_xc_and_setups(xc_input, hubbard_u=None):
    """Resolve GPAW XC and setup specifications."""
    setups = resolve_gpaw_hubbard(hubbard_u)

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
    maxiter=None,
    eigensolver=None,
):
    """Build calculator arguments shared by GPAW ground-state modes."""
    kwargs = {
        'nbands': nbands,
        'mixer': mixer,
        'charge': charge,
        'spinpol': spinpol,
        'txt': txt,
        'convergence': convergence,
        'occupations': occupations,
    }

    if maxiter is not None:
        kwargs['maxiter'] = maxiter

    if eigensolver is not None:
        kwargs['eigensolver'] = eigensolver

    return kwargs

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
    nbands=None,
    maxiter=None,
    eigensolver=None,
):
    """Create a regular GPAW plane-wave ground-state calculator."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
        nbands='200%' if nbands is None else nbands,
        maxiter=maxiter,
        eigensolver=eigensolver,
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
    nbands=None,
    maxiter=None,
    eigensolver=None,
):
    """Create a hybrid GPAW plane-wave ground-state calculator."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
        nbands='200%' if nbands is None else nbands,
        maxiter=maxiter,
        eigensolver=eigensolver,
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
        'kpts': build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    })

    if eigensolver is None:
        kwargs['eigensolver'] = Davidson(niter=1)

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
    nbands=None,
    maxiter=None,
    eigensolver=None,
):
    """Create a GPAW LCAO ground-state calculator."""
    kwargs = build_ground_common_kwargs(
        mixer=mixer,
        charge=charge,
        spinpol=spinpol,
        txt=txt,
        convergence=convergence,
        occupations=occupations,
        nbands='200%' if nbands is None else nbands,
        maxiter=maxiter,
        eigensolver=eigensolver,
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
    kpoint_density,
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
        'kpts': build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    })

    if hybrid:
        kwargs['eigensolver'] = Davidson(niter=1)

    return create_gpaw_calc(**kwargs)

def create_phonon_calc(
    cutoff,
    kpoint_size,
    txt,
    ground_calc=None,
    supercell_multiplier=1,
):
    """Create forces on the ground-state electronic energy surface.

    ``new`` retains XC, Hubbard setups, spin, charge, occupations,
    convergence, and mixer settings without reusing the unit-cell density.
    Extensive quantities (charge and absolute band count) scale with the
    number of repeated unit cells.
    """
    if ground_calc is not None:
        kwargs = {
            'mode': PW(cutoff),
            'kpts': {'size': tuple(kpoint_size)},
            'txt': txt,
        }
        # Absolute band counts belong to the unit cell. Relative counts
        # such as '200%' already adjust to the displaced supercell.
        nbands = ground_calc.parameters.get('nbands')
        if isinstance(nbands, Integral):
            kwargs['nbands'] = int(nbands) * int(supercell_multiplier)
        charge = ground_calc.parameters.get('charge', 0.0)
        if charge:
            kwargs['charge'] = charge * int(supercell_multiplier)
        return ground_calc.new(**kwargs)
    return create_gpaw_calc(
        mode=PW(cutoff),
        kpts={
            'size': tuple(kpoint_size),
        },
        txt=txt,
    )

def resolve_elastic_settings(
    xc_calc,
    hubbard_u,
    world_size,
    exx_fraction=None,
    omega=None,
    backend='pw',
):
    """Resolve XC, setups, parallel settings, and hybrid state for elasticity."""
    actual_xc, resolved_setups, _ = resolve_xc_and_setups(
        xc_calc,
        hubbard_u,
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
    return create_mixer(0.1)


def create_mixer(beta):
    """Create the GPAW density mixer for a portable mixing value."""
    return MixerSum(
        beta=float(beta),
        nmaxold=3,
        weight=50,
    )

def prepare_optical_calc(
    filename,
    hybrid,
    txt,
    nbands,
    smearing,
    kpoint_density,
    kpoint_size,
    gamma,
):
    """Prepare a GPAW calculator for optical-response calculations."""
    parallel = {
        'domain': 1,
    }

    if hybrid:
        return load_gpaw_calc(
            filename,
            hybrid=True,
            txt=txt,
            parallel=parallel,
        )

    return load_gpaw_calc(
        filename,
    ).fixed_density(
        txt=txt,
        nbands=nbands,
        parallel=parallel,
        occupations=FermiDirac(smearing),
        kpts=build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    )

def prepare_dos_calc(
    filename,
    hybrid,
    txt,
    convergence,
    occupations,
    kpoint_density,
    kpoint_size,
    gamma,
    nbands,
):
    """Prepare a GPAW calculator for DOS calculations."""
    if hybrid:
        return load_gpaw_calc(
            filename,
            hybrid=True,
        )

    calc = load_gpaw_calc(filename)

    # Safe parameter cleanup compatible with GPAW 26.7.0+
    try:
        if hasattr(calc.parameters, 'pop'):
            calc.parameters.pop('extensions', None)
        elif hasattr(calc.parameters, 'extensions'):
            delattr(calc.parameters, 'extensions')
    except Exception:
        pass

    kwargs = {
        'txt': txt,
        'convergence': convergence,
        'occupations': occupations,
        'kpts': build_kpoint_spec(
            density=kpoint_density,
            size=kpoint_size,
            gamma=gamma,
        ),
    }

    if nbands is not None:
        kwargs['nbands'] = nbands

    return calc.fixed_density(**kwargs)

def prepare_band_calc(
    filename,
    hybrid,
    path,
    npoints,
    txt,
    occupations,
    convergence,
    nbands,
):
    """Prepare a GPAW calculator for band-structure calculations."""
    kpts = {
        'path': path,
        'npoints': npoints,
    }

    if hybrid:
        return load_gpaw_calc(
            filename,
            hybrid=True,
            symmetry='off',
            kpts=kpts,
            parallel={
                'band': 1,
                'kpt': 1,
            },
            occupations=occupations,
            txt=txt,
            convergence=convergence,
        )

    kwargs = {
        'kpts': kpts,
        'txt': txt,
        'symmetry': 'off',
        'occupations': occupations,
        'convergence': convergence,
    }

    if nbands is not None:
        kwargs['nbands'] = nbands

    return load_gpaw_calc(
        filename,
    ).fixed_density(**kwargs)
