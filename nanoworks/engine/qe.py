"""Quantum ESPRESSO computation engine helpers."""

import re
import os
import shutil
import subprocess
from pathlib import Path

from ase.data import atomic_masses, atomic_numbers
from ase.calculators.calculator import kptdensity2monkhorstpack

QE_REFERENCE_VERSION = (7, 2)

# CODATA-compatible conversion used by ASE and QE-related workflows.
EV_PER_RYDBERG = 13.605693122994


def ev_to_rydberg(value):
    """Convert an energy value from electron-volts to Rydberg."""
    return float(value) / EV_PER_RYDBERG

def rydberg_to_ev(value):
    """Convert an energy value from Rydberg to electron-volts."""
    return float(value) * EV_PER_RYDBERG


def build_control_settings(
    calculation='scf',
    prefix='nanoworks',
    pseudo_dir=None,
    outdir=None,
):
    """Build the QE &CONTROL namelist settings."""
    calculation = str(
        calculation
    ).lower()

    settings = {
        'calculation': calculation,
        'prefix': str(prefix),
    }

    if calculation == 'bands':
        settings['verbosity'] = 'high'

    if pseudo_dir is not None:
        settings['pseudo_dir'] = str(pseudo_dir)

    if outdir is not None:
        settings['outdir'] = str(outdir)

    return settings


def build_system_settings(
    cutoff_ev,
    nat,
    ntyp,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
):
    """Build the basic QE &SYSTEM namelist settings."""
    settings = {
        'ibrav': 0,
        'nat': int(nat),
        'ntyp': int(ntyp),
        'ecutwfc': ev_to_rydberg(cutoff_ev),
    }

    if total_charge != 0.0:
        settings['tot_charge'] = float(total_charge)

    if nbands is not None:
        settings['nbnd'] = int(nbands)

    if spinpol:
        settings['nspin'] = 2

    return settings

def build_cell_parameters(atoms):
    """Build QE CELL_PARAMETERS data from an ASE Atoms object."""
    vectors = [
        tuple(float(value) for value in vector)
        for vector in atoms.cell.array
    ]

    return {
        'option': 'angstrom',
        'vectors': vectors,
    }


def build_atomic_positions(atoms):
    """Build QE ATOMIC_POSITIONS data from an ASE Atoms object."""
    positions = []

    for symbol, position in zip(
        atoms.get_chemical_symbols(),
        atoms.get_positions(),
    ):
        positions.append(
            (
                symbol,
                float(position[0]),
                float(position[1]),
                float(position[2]),
            )
        )

    return {
        'option': 'angstrom',
        'positions': positions,
    }


def build_atomic_species(atoms, pseudopotentials):
    """Build QE ATOMIC_SPECIES data from an ASE Atoms object."""
    symbols = list(dict.fromkeys(atoms.get_chemical_symbols()))

    missing = [
        symbol
        for symbol in symbols
        if symbol not in pseudopotentials
    ]

    if missing:
        raise ValueError(
            "Missing pseudopotential mapping for: "
            + ", ".join(missing)
        )

    species = []

    for symbol in symbols:
        atomic_number = atomic_numbers[symbol]
        mass = float(atomic_masses[atomic_number])

        species.append(
            (
                symbol,
                mass,
                str(pseudopotentials[symbol]),
            )
        )

    return species

def resolve_qe_kpoint_size(
    atoms,
    density=None,
    size=(5, 5, 5),
):
    """Resolve Nanoworks k-point settings to an explicit QE mesh."""
    if density is not None:
        density = float(density)

        if density <= 0.0:
            raise ValueError(
                "QE k-point density must be greater than zero."
            )

        mesh = kptdensity2monkhorstpack(
            atoms,
            kptdensity=density,
            even=None,
        )

        return tuple(
            int(value)
            for value in mesh
        )

    mesh = tuple(
        int(value)
        for value in size
    )

    if len(mesh) != 3:
        raise ValueError(
            "QE k-point mesh must contain exactly 3 values."
        )

    if any(value <= 0 for value in mesh):
        raise ValueError(
            "QE k-point mesh values must be positive integers."
        )

    return mesh

def build_kpoint_settings(size, gamma=False):
    """Build a QE automatic K_POINTS mesh.

    Nanoworks ``gamma=True`` means a Gamma-centered mesh, not a
    Gamma-only calculation.
    """
    mesh = tuple(int(value) for value in size)

    if len(mesh) != 3:
        raise ValueError("QE k-point mesh must contain exactly 3 values.")

    if any(value <= 0 for value in mesh):
        raise ValueError("QE k-point mesh values must be positive integers.")

    if gamma:
        shifts = (0, 0, 0)
    else:
        shifts = tuple(
            1 if value % 2 == 0 else 0
            for value in mesh
        )

    return {
        'option': 'automatic',
        'size': mesh,
        'shift': shifts,
    }

def build_band_path(atoms, path, npoints):
    """Build an explicit QE band path from an ASE Atoms object."""
    npoints = int(npoints)

    if npoints < 2:
        raise ValueError(
            "QE band path must contain at least 2 k-points."
        )

    if path is not None:
        path = str(path).strip()

        if not path:
            raise ValueError(
                "QE band path must not be empty."
            )

    band_path = atoms.cell.bandpath(
        path=path,
        npoints=npoints,
    )

    distances, special_distances, labels = (
        band_path.get_linear_kpoint_axis()
    )

    kpoints = [
        tuple(float(value) for value in kpoint)
        for kpoint in band_path.kpts
    ]

    return {
        'option': 'crystal',
        'path': band_path.path,
        'kpoints': kpoints,
        'distances': [
            float(value)
            for value in distances
        ],
        'special_distances': [
            float(value)
            for value in special_distances
        ],
        'labels': list(labels),
        'npoints': len(kpoints),
    }

def render_band_kpoints(settings):
    """Render an explicit QE K_POINTS crystal card."""
    option = str(
        settings.get('option', '')
    ).strip().lower()

    if option != 'crystal':
        raise ValueError(
            "QE explicit band k-points must use crystal coordinates."
        )

    kpoints = list(
        settings.get('kpoints', [])
    )

    if not kpoints:
        raise ValueError(
            "QE band path does not contain any k-points."
        )

    declared_npoints = settings.get('npoints')

    if (
        declared_npoints is not None
        and int(declared_npoints) != len(kpoints)
    ):
        raise ValueError(
            "QE band path point count does not match its metadata."
        )

    lines = [
        'K_POINTS crystal',
        str(len(kpoints)),
    ]

    for kpoint in kpoints:
        if len(kpoint) != 3:
            raise ValueError(
                "Each QE band k-point must contain exactly 3 coordinates."
            )

        x, y, z = (
            float(value)
            for value in kpoint
        )

        lines.append(
            f"{x:.12f} {y:.12f} {z:.12f} 1.0"
        )

    return "\n".join(lines)

def validate_qe_version(
    version,
    minimum=QE_REFERENCE_VERSION,
):
    """Validate the Quantum ESPRESSO version used by a calculation."""
    if version is None:
        raise ValueError(
            "Quantum ESPRESSO version could not be detected "
            "from the pw.x output."
        )

    version = tuple(
        int(value)
        for value in version
    )

    minimum = tuple(
        int(value)
        for value in minimum
    )

    version_major_minor = version[:2]

    if version_major_minor < minimum:
        detected = '.'.join(
            str(value)
            for value in version
        )

        required = '.'.join(
            str(value)
            for value in minimum
        )

        raise ValueError(
            "Unsupported Quantum ESPRESSO version "
            f"{detected}. Nanoworks currently requires "
            f"Quantum ESPRESSO {required} or newer."
        )

    return version

def validate_qe_xc(
    xc_calc,
    pseudo_xc='pbe',
):
    """Validate XC compatibility with the installed QE pseudo set."""
    xc = str(
        xc_calc
    ).strip().lower()

    pseudo_xc = str(
        pseudo_xc
    ).strip().lower()

    aliases = {
        'pbe': 'pbe',
    }

    try:
        normalized = aliases[xc]
    except KeyError:
        raise ValueError(
            "The current Nanoworks QE pseudopotential library "
            f"supports PBE calculations only. Requested XC: {xc_calc}"
        )

    if normalized != pseudo_xc:
        raise ValueError(
            f"QE XC '{xc_calc}' is incompatible with "
            f"the installed '{pseudo_xc}' pseudopotentials."
        )

    return normalized

def resolve_qe_occupation(occupation):
    """Translate Nanoworks/GPAW-style occupation settings to QE settings."""
    if occupation is None:
        return {
            'occupations': 'fixed',
            'smearing': None,
            'width_ev': None,
        }

    if isinstance(occupation, str):
        name = occupation.strip().lower()
        width = None

    elif isinstance(occupation, dict):
        name = str(
            occupation.get(
                'name',
                'fixed',
            )
        ).strip().lower()

        width = occupation.get(
            'width'
        )

    else:
        raise TypeError(
            "QE occupation settings must be a string, "
            "dictionary, or None."
        )

    fixed_names = {
        'fixed',
    }

    if name in fixed_names:
        return {
            'occupations': 'fixed',
            'smearing': None,
            'width_ev': None,
        }
    
    tetrahedra_names = {
        'tetrahedra': 'tetrahedra',
        'tetrahedra_lin': 'tetrahedra_lin',
        'tetrahedra-lin': 'tetrahedra_lin',
        'tetrahedra_opt': 'tetrahedra_opt',
        'tetrahedra-opt': 'tetrahedra_opt',
    }

    if name in tetrahedra_names:
        return {
            'occupations': tetrahedra_names[name],
            'smearing': None,
            'width_ev': None,
        }

    smearing_aliases = {
        'fermi-dirac': 'fermi-dirac',
        'fermi_dirac': 'fermi-dirac',
        'fd': 'fermi-dirac',

        'gaussian': 'gaussian',
        'gauss': 'gaussian',

        'methfessel-paxton': 'methfessel-paxton',
        'methfessel_paxton': 'methfessel-paxton',
        'mp': 'methfessel-paxton',

        'marzari-vanderbilt': 'marzari-vanderbilt',
        'marzari_vanderbilt': 'marzari-vanderbilt',
        'cold': 'marzari-vanderbilt',
    }

    try:
        smearing = smearing_aliases[name]
    except KeyError:
        raise ValueError(
            f"Unsupported Nanoworks occupation scheme for QE: {name}"
        )

    if width is None:
        raise ValueError(
            f"QE smearing occupation '{name}' requires a width."
        )

    return {
        'occupations': 'smearing',
        'smearing': smearing,
        'width_ev': float(width),
    }

def build_occupation_settings(
    occupations='fixed',
    smearing=None,
    width_ev=None,
):
    """Build QE occupation-related &SYSTEM settings."""
    occupation = str(occupations).strip().lower()

    allowed_occupations = {
        'fixed',
        'smearing',
        'tetrahedra',
        'tetrahedra_lin',
        'tetrahedra_opt',
    }

    if occupation not in allowed_occupations:
        raise ValueError(
            f"Unsupported QE occupation scheme: {occupations}"
        )

    settings = {
        'occupations': occupation,
    }

    if occupation != 'smearing':
        if smearing is not None or width_ev is not None:
            raise ValueError(
                "Smearing type and width can only be used with "
                "occupations='smearing'."
            )

        return settings

    if smearing is None:
        raise ValueError(
            "QE smearing calculations require a smearing type."
        )

    if width_ev is None:
        raise ValueError(
            "QE smearing calculations require a smearing width."
        )

    aliases = {
        'gaussian': 'gaussian',
        'gauss': 'gaussian',

        'methfessel-paxton': 'methfessel-paxton',
        'm-p': 'methfessel-paxton',
        'mp': 'methfessel-paxton',

        'marzari-vanderbilt': 'marzari-vanderbilt',
        'cold': 'marzari-vanderbilt',
        'm-v': 'marzari-vanderbilt',
        'mv': 'marzari-vanderbilt',

        'fermi-dirac': 'fermi-dirac',
        'f-d': 'fermi-dirac',
        'fd': 'fermi-dirac',
    }

    smearing_key = str(smearing).strip().lower()

    try:
        qe_smearing = aliases[smearing_key]
    except KeyError:
        raise ValueError(
            f"Unsupported QE smearing type: {smearing}"
        )

    width_ev = float(width_ev)

    if width_ev <= 0.0:
        raise ValueError(
            "QE smearing width must be greater than zero."
        )

    settings['smearing'] = qe_smearing
    settings['degauss'] = ev_to_rydberg(width_ev)

    return settings

def build_electrons_settings(
    conv_thr=None,
    mixing_beta=None,
    electron_maxstep=None,
    diagonalization=None,
):
    """Build the QE &ELECTRONS namelist settings."""
    settings = {}

    if conv_thr is not None:
        conv_thr = float(conv_thr)

        if conv_thr <= 0.0:
            raise ValueError(
                "QE electronic convergence threshold must be greater than zero."
            )

        settings['conv_thr'] = conv_thr

    if mixing_beta is not None:
        mixing_beta = float(mixing_beta)

        if not 0.0 < mixing_beta <= 1.0:
            raise ValueError(
                "QE mixing_beta must be greater than zero and at most one."
            )

        settings['mixing_beta'] = mixing_beta

    if electron_maxstep is not None:
        electron_maxstep = int(electron_maxstep)

        if electron_maxstep <= 0:
            raise ValueError(
                "QE electron_maxstep must be a positive integer."
            )

        settings['electron_maxstep'] = electron_maxstep

    if diagonalization is not None:
        diagonalization = str(diagonalization).strip().lower()

        allowed = {
            'david',
            'cg',
        }

        if diagonalization not in allowed:
            raise ValueError(
                f"Unsupported QE diagonalization method: {diagonalization}"
            )

        settings['diagonalization'] = diagonalization

    return settings

def format_qe_value(value):
    """Format a Python value for a QE namelist."""
    if isinstance(value, bool):
        return '.true.' if value else '.false.'

    if isinstance(value, str):
        escaped = value.replace("'", "''")
        return f"'{escaped}'"

    if isinstance(value, int):
        return str(value)

    if isinstance(value, float):
        return f"{value:.12g}"

    raise TypeError(
        f"Unsupported QE namelist value type: {type(value).__name__}"
    )

def render_namelist(name, settings):
    """Render one QE namelist."""
    lines = [f"&{str(name).upper()}"]

    for key, value in settings.items():
        lines.append(
            f"  {key} = {format_qe_value(value)},"
        )

    lines.append("/")

    return "\n".join(lines)

def render_pw_input(
    calculation,
    atoms,
    pseudopotentials,
    cutoff_ev,
    kpoint_size,
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupations='fixed',
    smearing=None,
    width_ev=None,
    prefix='nanoworks',
    pseudo_dir=None,
    outdir=None,
    conv_thr=None,
    mixing_beta=None,
    electron_maxstep=None,
    diagonalization=None,
    band_path=None,
):
    """Render a complete QE pw.x input."""

    calculation = str(
        calculation
    ).strip().lower()

    allowed_calculations = {
        'scf',
        'nscf',
        'bands',
    }

    if calculation not in allowed_calculations:
        raise ValueError(
            f"Unsupported QE pw.x calculation type: {calculation}"
        )
    
    species = build_atomic_species(
        atoms,
        pseudopotentials,
    )

    control = build_control_settings(
        calculation=calculation,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=outdir,
    )

    system = build_system_settings(
        cutoff_ev=cutoff_ev,
        nat=len(atoms),
        ntyp=len(species),
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
    )

    system.update(
        build_occupation_settings(
            occupations=occupations,
            smearing=smearing,
            width_ev=width_ev,
        )
    )

    electrons = build_electrons_settings(
        conv_thr=conv_thr,
        mixing_beta=mixing_beta,
        electron_maxstep=electron_maxstep,
        diagonalization=diagonalization,
    )

    positions = build_atomic_positions(atoms)
    cell = build_cell_parameters(atoms)

    if calculation == 'bands':
        if band_path is None:
            raise ValueError(
                "QE bands calculation requires an explicit band path."
            )

        kpoint_card = render_band_kpoints(
            band_path
        )

    else:
        if band_path is not None:
            raise ValueError(
                "QE band paths can only be used with "
                "calculation='bands'."
            )

        kpoints = build_kpoint_settings(
            kpoint_size,
            gamma=gamma,
        )

        nk1, nk2, nk3 = kpoints['size']
        sk1, sk2, sk3 = kpoints['shift']

        kpoint_card = "\n".join([
            f"K_POINTS {kpoints['option']}",
            f"{nk1} {nk2} {nk3} {sk1} {sk2} {sk3}",
        ])

    lines = [
        render_namelist('CONTROL', control),
        render_namelist('SYSTEM', system),
        render_namelist('ELECTRONS', electrons),
        '',
        'ATOMIC_SPECIES',
    ]

    for symbol, mass, pseudo in species:
        lines.append(
            f"{symbol} {mass:.8f} {pseudo}"
        )

    lines.extend([
        '',
        f"ATOMIC_POSITIONS {positions['option']}",
    ])

    for symbol, x, y, z in positions['positions']:
        lines.append(
            f"{symbol} {x:.12f} {y:.12f} {z:.12f}"
        )

    lines.extend([
        '',
        kpoint_card,
    ])

    lines.extend([
        '',
        f"CELL_PARAMETERS {cell['option']}",
    ])

    for x, y, z in cell['vectors']:
        lines.append(
            f"{x:.12f} {y:.12f} {z:.12f}"
        )

    return "\n".join(lines) + "\n"

def render_scf_input(
    atoms,
    pseudopotentials,
    cutoff_ev,
    kpoint_size,
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupations='fixed',
    smearing=None,
    width_ev=None,
    prefix='nanoworks',
    pseudo_dir=None,
    outdir=None,
    conv_thr=None,
    mixing_beta=None,
    electron_maxstep=None,
    diagonalization=None,
):
    """Render a complete QE pw.x SCF input."""
    return render_pw_input(
        calculation='scf',
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=kpoint_size,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        occupations=occupations,
        smearing=smearing,
        width_ev=width_ev,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=outdir,
        conv_thr=conv_thr,
        mixing_beta=mixing_beta,
        electron_maxstep=electron_maxstep,
        diagonalization=diagonalization,
    )

def render_nscf_input(
    atoms,
    pseudopotentials,
    cutoff_ev,
    kpoint_size,
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupations='fixed',
    smearing=None,
    width_ev=None,
    prefix='nanoworks',
    pseudo_dir=None,
    outdir=None,
    conv_thr=None,
    mixing_beta=None,
    electron_maxstep=None,
    diagonalization=None,
):
    """Render a complete QE pw.x NSCF input."""
    return render_pw_input(
        calculation='nscf',
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=kpoint_size,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        occupations=occupations,
        smearing=smearing,
        width_ev=width_ev,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=outdir,
        conv_thr=conv_thr,
        mixing_beta=mixing_beta,
        electron_maxstep=electron_maxstep,
        diagonalization=diagonalization,
    )

def render_bands_input(
    atoms,
    pseudopotentials,
    cutoff_ev,
    band_path,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupations='fixed',
    smearing=None,
    width_ev=None,
    prefix='nanoworks',
    pseudo_dir=None,
    outdir=None,
    conv_thr=None,
    mixing_beta=None,
    electron_maxstep=None,
    diagonalization=None,
):
    """Render a complete QE pw.x bands input."""
    return render_pw_input(
        calculation='bands',
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=None,
        gamma=False,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        occupations=occupations,
        smearing=smearing,
        width_ev=width_ev,
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=outdir,
        conv_thr=conv_thr,
        mixing_beta=mixing_beta,
        electron_maxstep=electron_maxstep,
        diagonalization=diagonalization,
        band_path=band_path,
    )

def render_dos_input(
    prefix='nanoworks',
    outdir=None,
    fildos='nanoworks.dos',
    bz_sum=None,
    emin=None,
    emax=None,
    delta_e=None,
    degauss=None,
    ngauss=None,
):
    """Render a complete Quantum ESPRESSO dos.x input."""
    settings = {
        'prefix': str(prefix),
    }

    if outdir is not None:
        settings['outdir'] = str(outdir)

    if bz_sum is not None:
        bz_sum = str(
            bz_sum
        ).strip().lower()

        allowed_bz_sum = {
            'smearing',
            'tetrahedra',
            'tetrahedra_lin',
            'tetrahedra_opt',
        }

        if bz_sum not in allowed_bz_sum:
            raise ValueError(
                f"Unsupported QE DOS BZ summation method: {bz_sum}"
            )

        settings['bz_sum'] = bz_sum

    if emin is not None:
        settings['Emin'] = float(
            emin
        )

    if emax is not None:
        settings['Emax'] = float(
            emax
        )

    if (
        emin is not None
        and emax is not None
        and float(emax) <= float(emin)
    ):
        raise ValueError(
            "QE DOS Emax must be greater than Emin."
        )

    if delta_e is not None:
        delta_e = float(
            delta_e
        )

        if delta_e <= 0.0:
            raise ValueError(
                "QE DOS energy step must be greater than zero."
            )

        settings['DeltaE'] = delta_e

    if degauss is not None:
        degauss = float(
            degauss
        )

        if degauss <= 0.0:
            raise ValueError(
                "QE DOS degauss must be greater than zero."
            )

        settings['degauss'] = degauss

    if ngauss is not None:
        settings['ngauss'] = int(
            ngauss
        )

    settings['fildos'] = str(
        fildos
    )

    return (
        render_namelist(
            'DOS',
            settings,
        )
        + '\n'
    )

def render_projwfc_input(
    prefix='nanoworks',
    outdir=None,
    filpdos='nanoworks',
    emin=None,
    emax=None,
    delta_e=None,
    degauss=None,
    ngauss=None,
    lsym=True,
    diag_basis=False,
):
    """Render a complete Quantum ESPRESSO projwfc.x input."""
    settings = {
        'prefix': str(prefix),
    }

    if outdir is not None:
        settings['outdir'] = str(outdir)

    if emin is not None:
        settings['Emin'] = float(
            emin
        )

    if emax is not None:
        settings['Emax'] = float(
            emax
        )

    if (
        emin is not None
        and emax is not None
        and float(emax) <= float(emin)
    ):
        raise ValueError(
            "QE PDOS Emax must be greater than Emin."
        )

    if delta_e is not None:
        delta_e = float(
            delta_e
        )

        if delta_e <= 0.0:
            raise ValueError(
                "QE PDOS energy step must be greater than zero."
            )

        settings['DeltaE'] = delta_e

    if degauss is not None:
        degauss = float(
            degauss
        )

        if degauss <= 0.0:
            raise ValueError(
                "QE PDOS degauss must be greater than zero."
            )

        settings['degauss'] = degauss

    if ngauss is not None:
        settings['ngauss'] = int(
            ngauss
        )

    settings['lsym'] = bool(
        lsym
    )

    settings['diag_basis'] = bool(
        diag_basis
    )

    settings['filpdos'] = str(
        filpdos
    )

    return (
        render_namelist(
            'PROJWFC',
            settings,
        )
        + '\n'
    )

def build_qe_launcher(
    parallel_cores=1,
):
    """Build the MPI launcher for a Quantum ESPRESSO calculation."""
    parallel_cores = int(
        parallel_cores
    )

    if parallel_cores <= 0:
        raise ValueError(
            "QE parallel core count must be a positive integer."
        )

    if parallel_cores == 1:
        return None

    mpi_exe = (
        shutil.which('mpiexec')
        or shutil.which('mpirun')
        or shutil.which('srun')
    )

    if mpi_exe is None:
        raise FileNotFoundError(
            "mpiexec, mpirun, or srun was not found "
            "for QE parallel execution."
        )

    if 'srun' in Path(mpi_exe).name:
        flag = '-n'
    else:
        flag = '-np'

    return [
        mpi_exe,
        flag,
        str(parallel_cores),
    ]

def has_qe_state(
    state_dir,
    prefix='nanoworks',
):
    """Return True when a usable QE saved state is present."""
    state_dir = Path(
        state_dir
    )

    save_dir = (
        state_dir
        / f'{prefix}.save'
    )

    if not save_dir.is_dir():
        return False

    schema_file = (
        save_dir
        / 'data-file-schema.xml'
    )

    if not schema_file.is_file():
        return False

    return True

def resolve_qe_executable(
    executable='pw.x',
):
    """Resolve a Quantum ESPRESSO executable."""
    executable = str(executable)

    path = Path(executable).expanduser()

    if path.parent != Path('.'):
        if not path.exists():
            raise FileNotFoundError(
                f"Quantum ESPRESSO executable was not found: {path}"
            )

        if not path.is_file():
            raise FileNotFoundError(
                f"Quantum ESPRESSO executable is not a file: {path}"
            )

        return str(path.resolve())

    resolved = shutil.which(executable)

    if resolved is None:
        raise FileNotFoundError(
            f"Quantum ESPRESSO executable '{executable}' "
            "was not found in PATH."
        )

    return resolved

def build_qe_command(
    input_file,
    executable='pw.x',
    launcher=None,
):
    """Build a Quantum ESPRESSO execution command."""
    executable = resolve_qe_executable(
        executable
    )

    command = []

    if launcher is not None:
        if isinstance(launcher, str):
            raise TypeError(
                "QE launcher must be a sequence of command arguments, "
                "not a shell command string."
            )

        command.extend(
            str(value)
            for value in launcher
        )

    command.extend([
        executable,
        '-i',
        str(input_file),
    ])

    return command

def run_qe_program(
    input_file,
    output_file,
    executable='pw.x',
    launcher=None,
    cwd=None,
):
    """Run a Quantum ESPRESSO program and write its output to a file."""
    input_file = Path(
        input_file
    ).expanduser()

    output_file = Path(
        output_file
    ).expanduser()

    if not input_file.exists():
        raise FileNotFoundError(
            f"QE input file was not found: {input_file}"
        )

    if cwd is not None:
        cwd = Path(
            cwd
        ).expanduser()

    command = build_qe_command(
        input_file=input_file,
        executable=executable,
        launcher=launcher,
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    child_env = os.environ.copy()

    # Prevent MPI ranks from spawning additional BLAS/OpenMP threads.
    child_env['OMP_NUM_THREADS'] = '1'
    child_env['OPENBLAS_NUM_THREADS'] = '1'
    child_env['MKL_NUM_THREADS'] = '1'
    child_env['VECLIB_MAXIMUM_THREADS'] = '1'
    child_env['NUMEXPR_NUM_THREADS'] = '1'
    child_env['OMP_DYNAMIC'] = 'FALSE'

    with output_file.open(
        'w',
        encoding='utf-8',
    ) as fd:
        result = subprocess.run(
            command,
            cwd=cwd,
            stdout=fd,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
            env=child_env,
        )

    if result.returncode != 0:
        raise RuntimeError(
            "Quantum ESPRESSO calculation failed with "
            f"return code {result.returncode}. "
            f"See '{output_file}'."
        )

    return {
        'command': command,
        'returncode': result.returncode,
        'output_file': output_file,
    }

def parse_pw_output(output):
    """Parse basic results from pw.x output."""
    output = Path(
        output
    )

    if not output.exists():
        raise FileNotFoundError(
            f"QE output file was not found: {output}"
        )

    text = output.read_text(
        encoding='utf-8',
        errors='replace',
    )
    
    version_match = re.search(
        r'Program\s+PWSCF\s+v\.'
        r'(\d+)\.(\d+)(?:\.(\d+))?',
        text,
        flags=re.IGNORECASE,
    )

    qe_version = None

    if version_match:
        qe_version = tuple(
            int(value)
            for value in version_match.groups()
            if value is not None
        )

    energy_matches = re.findall(
        r'!\s+total energy\s*=\s*'
        r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)'
        r'\s+Ry',
        text,
        flags=re.IGNORECASE,
    )

    total_energy_ry = None
    total_energy_ev = None

    if energy_matches:
        value = (
            energy_matches[-1]
            .replace('D', 'E')
            .replace('d', 'e')
        )

        total_energy_ry = float(
            value
        )

        total_energy_ev = rydberg_to_ev(
            total_energy_ry
        )
    
    fermi_matches = re.findall(
        r'the\s+Fermi\s+energy\s+is\s+'
        r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)'
        r'\s+ev',
        text,
        flags=re.IGNORECASE,
    )

    fermi_energy_ev = None

    if fermi_matches:
        value = (
            fermi_matches[-1]
            .replace('D', 'E')
            .replace('d', 'e')
        )

        fermi_energy_ev = float(
            value
        )

    number_pattern = (
        r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)'
        r'(?:[EeDd][-+]?\d+)?'
    )

    band_edge_matches = re.findall(
        r'highest\s+occupied\s*,\s*'
        r'lowest\s+unoccupied\s+level\s*'
        r'\(\s*ev\s*\)\s*:\s*'
        rf'({number_pattern})\s+'
        rf'({number_pattern})',
        text,
        flags=re.IGNORECASE,
    )

    highest_occupied_matches = re.findall(
        r'highest\s+occupied\s+level\s*'
        r'\(\s*ev\s*\)\s*:\s*'
        rf'({number_pattern})',
        text,
        flags=re.IGNORECASE,
    )

    highest_occupied_ev = None
    lowest_unoccupied_ev = None

    if band_edge_matches:
        highest, lowest = (
            band_edge_matches[-1]
        )

        highest_occupied_ev = float(
            highest
            .replace('D', 'E')
            .replace('d', 'e')
        )

        lowest_unoccupied_ev = float(
            lowest
            .replace('D', 'E')
            .replace('d', 'e')
        )

    elif highest_occupied_matches:
        highest = (
            highest_occupied_matches[-1]
        )

        highest_occupied_ev = float(
            highest
            .replace('D', 'E')
            .replace('d', 'e')
        )

    job_done = (
        'JOB DONE.' in text
    )

    return {
        'job_done': job_done,
        'qe_version': qe_version,
        'total_energy_ry': total_energy_ry,
        'total_energy_ev': total_energy_ev,
        'fermi_energy_ev': fermi_energy_ev,
        'highest_occupied_ev': highest_occupied_ev,
        'lowest_unoccupied_ev': lowest_unoccupied_ev,
    }

def resolve_qe_band_reference(result):
    """Resolve the energy reference used for QE band outputs."""
    fermi_energy = result.get(
        'fermi_energy_ev'
    )

    if fermi_energy is not None:
        return {
            'energy_ev': float(
                fermi_energy
            ),
            'source': 'fermi',
        }

    highest_occupied = result.get(
        'highest_occupied_ev'
    )

    lowest_unoccupied = result.get(
        'lowest_unoccupied_ev'
    )

    if (
        highest_occupied is not None
        and lowest_unoccupied is not None
    ):
        highest_occupied = float(
            highest_occupied
        )

        lowest_unoccupied = float(
            lowest_unoccupied
        )

        if lowest_unoccupied < highest_occupied:
            raise ValueError(
                "QE lowest unoccupied level is below "
                "the highest occupied level."
            )

        return {
            'energy_ev': (
                highest_occupied
                + lowest_unoccupied
            ) / 2.0,
            'source': 'midgap',
        }

    if highest_occupied is not None:
        return {
            'energy_ev': float(
                highest_occupied
            ),
            'source': 'highest_occupied',
        }

    raise ValueError(
        "QE band reference energy could not be resolved."
    )

def parse_pw_bands_output(output):
    """Parse non-spin band eigenvalues from a QE pw.x output."""
    output = Path(
        output
    )

    if not output.exists():
        raise FileNotFoundError(
            f"QE bands output file was not found: {output}"
        )

    text = output.read_text(
        encoding='utf-8',
        errors='replace',
    )

    if re.search(
        r'\bSPIN\s+(?:UP|DOWN)\b',
        text,
        flags=re.IGNORECASE,
    ):
        raise NotImplementedError(
            "Spin-polarized QE band output parsing "
            "is not supported yet."
        )

    number_pattern = (
        r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)'
        r'(?:[EeDd][-+]?\d+)?'
    )

    kpoint_pattern = re.compile(
        r'^\s*k\s*=\s*'
        rf'({number_pattern})\s*'
        rf'({number_pattern})\s*'
        rf'({number_pattern})(?=\s|\()',
        flags=re.IGNORECASE,
    )

    bands_pattern = re.compile(
        r'bands\s*\(\s*ev\s*\)\s*:',
        flags=re.IGNORECASE,
    )

    number_line_pattern = re.compile(
        rf'^(?:\s*{number_pattern})+\s*$'
    )

    kpoints = []
    eigenvalues = []

    current_kpoint = None
    current_eigenvalues = None

    for line in text.splitlines():
        kpoint_match = kpoint_pattern.match(
            line
        )

        if kpoint_match:
            if (
                current_kpoint is not None
                and current_eigenvalues
            ):
                kpoints.append(
                    current_kpoint
                )
                eigenvalues.append(
                    current_eigenvalues
                )

            current_kpoint = tuple(
                float(
                    value
                    .replace('D', 'E')
                    .replace('d', 'e')
                )
                for value in kpoint_match.groups()
            )

            current_eigenvalues = None


        bands_match = bands_pattern.search(
            line
        )

        if (
            bands_match
            and current_kpoint is not None
        ):
            current_eigenvalues = []

            remainder = line[
                bands_match.end():
            ].strip()

            if remainder:
                current_eigenvalues.extend(
                    float(
                        value
                        .replace('D', 'E')
                        .replace('d', 'e')
                    )
                    for value in remainder.split()
                )

            continue

        if current_eigenvalues is None:
            continue

        stripped = line.strip()

        if not stripped:
            if current_eigenvalues:
                kpoints.append(
                    current_kpoint
                )
                eigenvalues.append(
                    current_eigenvalues
                )

                current_kpoint = None
                current_eigenvalues = None

            continue

        if number_line_pattern.fullmatch(
            line
        ):
            current_eigenvalues.extend(
                float(
                    value
                    .replace('D', 'E')
                    .replace('d', 'e')
                )
                for value in stripped.split()
            )
            continue

        if current_eigenvalues:
            kpoints.append(
                current_kpoint
            )
            eigenvalues.append(
                current_eigenvalues
            )

        current_kpoint = None
        current_eigenvalues = None

    if (
        current_kpoint is not None
        and current_eigenvalues
    ):
        kpoints.append(
            current_kpoint
        )
        eigenvalues.append(
            current_eigenvalues
        )

    if not eigenvalues:
        raise ValueError(
            "No QE band eigenvalues were found in "
            f"the output file: {output}"
        )

    band_counts = {
        len(values)
        for values in eigenvalues
    }

    if len(band_counts) != 1:
        raise ValueError(
            "QE band output contains inconsistent "
            "numbers of bands between k-points."
        )

    return {
        'spin_polarized': False,
        'nspins': 1,
        'nkpoints': len(kpoints),
        'nbands': len(eigenvalues[0]),
        'kpoints': kpoints,
        'eigenvalues_ev': [
            eigenvalues
        ],
    }

def prepare_qe_band_data(
    bands,
    band_path,
    reference_energy,
):
    """Prepare non-spin QE band data for output and plotting."""
    if (
        bands.get('spin_polarized')
        or int(bands.get('nspins', 0)) != 1
    ):
        raise NotImplementedError(
            "Spin-polarized QE band data preparation "
            "is not supported yet."
        )

    nkpoints = int(
        bands['nkpoints']
    )

    nbands = int(
        bands['nbands']
    )

    eigenvalues_by_spin = bands[
        'eigenvalues_ev'
    ]

    if len(eigenvalues_by_spin) != 1:
        raise ValueError(
            "Non-spin QE band data must contain "
            "exactly one spin channel."
        )

    eigenvalues = eigenvalues_by_spin[
        0
    ]

    if len(eigenvalues) != nkpoints:
        raise ValueError(
            "QE band eigenvalue count does not match "
            "the reported number of k-points."
        )

    if any(
        len(values) != nbands
        for values in eigenvalues
    ):
        raise ValueError(
            "QE band eigenvalue rows do not match "
            "the reported number of bands."
        )

    distances = [
        float(value)
        for value in band_path[
            'distances'
        ]
    ]

    if len(distances) != nkpoints:
        raise ValueError(
            "QE band-path distances do not match "
            "the number of parsed k-points."
        )

    special_distances = [
        float(value)
        for value in band_path[
            'special_distances'
        ]
    ]

    labels = list(
        band_path[
            'labels'
        ]
    )

    if len(special_distances) != len(labels):
        raise ValueError(
            "QE band-path labels and special-point "
            "distances do not match."
        )

    reference_energy = float(
        reference_energy
    )

    shifted_eigenvalues = [
        [
            float(value) - reference_energy
            for value in values
        ]
        for values in eigenvalues
    ]

    return {
        'nkpoints': nkpoints,
        'nbands': nbands,
        'distances': distances,
        'special_distances': special_distances,
        'labels': labels,
        'reference_energy_ev': reference_energy,
        'eigenvalues_ev': shifted_eigenvalues,
    }

def parse_dos_output(dos_file):
    """Parse a non-spin Quantum ESPRESSO dos.x data file."""
    dos_file = Path(
        dos_file
    )

    if not dos_file.is_file():
        raise FileNotFoundError(
            f"QE DOS data file was not found: {dos_file}"
        )

    energies = []
    dos_values = []
    integrated_dos = []

    with dos_file.open(
        'r',
        encoding='utf-8',
        errors='replace',
    ) as fd:
        for line in fd:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith('#'):
                continue

            parts = stripped.split()

            if len(parts) < 3:
                continue

            try:
                energy = float(
                    parts[0]
                    .replace('D', 'E')
                    .replace('d', 'e')
                )

                dos_value = float(
                    parts[1]
                    .replace('D', 'E')
                    .replace('d', 'e')
                )

                integrated = float(
                    parts[2]
                    .replace('D', 'E')
                    .replace('d', 'e')
                )
            except ValueError:
                continue

            energies.append(
                energy
            )

            dos_values.append(
                dos_value
            )

            integrated_dos.append(
                integrated
            )

    if not energies:
        raise ValueError(
            f"No DOS data could be parsed from '{dos_file}'."
        )

    return {
        'energies_ev': energies,
        'dos': dos_values,
        'integrated_dos': integrated_dos,
        'npoints': len(energies),
    }

def parse_projwfc_pdos_file(pdos_file):
    """Parse one non-spin Quantum ESPRESSO projwfc.x PDOS file."""
    pdos_file = Path(
        pdos_file
    )

    if not pdos_file.is_file():
        raise FileNotFoundError(
            f"QE PDOS data file was not found: {pdos_file}"
        )

    name = pdos_file.name

    match = re.search(
        r'\.pdos_atm#(\d+)\(([^)]+)\)_wfc#(\d+)\(([spdf])\)$',
        name,
    )

    if match is None:
        raise ValueError(
            f"Unsupported QE PDOS filename: {name}"
        )

    atom_index = int(
        match.group(1)
    )

    symbol = match.group(2)

    wfc_index = int(
        match.group(3)
    )

    orbital = match.group(4)

    component_names = {
        's': [
            's',
        ],
        'p': [
            'pz',
            'px',
            'py',
        ],
        'd': [
            'd3z2_r2',
            'dxz',
            'dyz',
            'dx2_y2',
            'dxy',
        ],
    }

    if orbital not in component_names:
        raise NotImplementedError(
            "QE PDOS parsing currently supports "
            "s, p, and d orbitals only."
        )

    energies = []
    ldos = []

    components = {
        name: []
        for name in component_names[orbital]
    }

    expected_columns = (
        2
        + len(component_names[orbital])
    )

    with pdos_file.open(
        'r',
        encoding='utf-8',
        errors='replace',
    ) as fd:
        for line in fd:
            stripped = line.strip()

            if not stripped:
                continue

            if stripped.startswith('#'):
                continue

            parts = stripped.split()

            if len(parts) < expected_columns:
                continue

            try:
                values = [
                    float(
                        value
                        .replace('D', 'E')
                        .replace('d', 'e')
                    )
                    for value in parts[
                        :expected_columns
                    ]
                ]
            except ValueError:
                continue

            energies.append(
                values[0]
            )

            ldos.append(
                values[1]
            )

            for index, component_name in enumerate(
                component_names[orbital]
            ):
                components[
                    component_name
                ].append(
                    values[index + 2]
                )

    if not energies:
        raise ValueError(
            f"No PDOS data could be parsed from '{pdos_file}'."
        )

    return {
        'atom_index': atom_index,
        'symbol': symbol,
        'wfc_index': wfc_index,
        'orbital': orbital,
        'energies_ev': energies,
        'ldos': ldos,
        'components': components,
        'npoints': len(energies),
    }

def aggregate_projwfc_pdos(pdos_prefix):
    """Aggregate non-spin QE projwfc.x atomic PDOS files."""
    pdos_prefix = Path(
        pdos_prefix
    )

    pattern = (
        pdos_prefix.name
        + '.pdos_atm#*'
    )

    pdos_files = sorted(
        pdos_prefix.parent.glob(
            pattern
        )
    )

    if not pdos_files:
        raise FileNotFoundError(
            "No Quantum ESPRESSO orbital PDOS files were found "
            f"for prefix: {pdos_prefix}"
        )

    energies = None

    totals = {
        's': None,
        'p': None,
        'd': None,
        'f': None,
    }

    components = {
        'pz': None,
        'px': None,
        'py': None,
        'd3z2_r2': None,
        'dxz': None,
        'dyz': None,
        'dx2_y2': None,
        'dxy': None,
    }

    parsed_files = []

    for pdos_file in pdos_files:
        parsed = parse_projwfc_pdos_file(
            pdos_file
        )

        parsed_files.append(
            parsed
        )

        file_energies = parsed[
            'energies_ev'
        ]

        if energies is None:
            energies = list(
                file_energies
            )

            npoints = len(
                energies
            )

            for orbital in totals:
                totals[orbital] = [
                    0.0
                ] * npoints

            for component in components:
                components[component] = [
                    0.0
                ] * npoints

        else:
            if len(file_energies) != len(energies):
                raise ValueError(
                    "QE PDOS files do not use the same "
                    "number of energy points."
                )

            for reference, value in zip(
                energies,
                file_energies,
            ):
                if abs(reference - value) > 1.0e-8:
                    raise ValueError(
                        "QE PDOS files do not use the same "
                        "energy grid."
                    )

        orbital = parsed[
            'orbital'
        ]

        if orbital == 'f':
            raise NotImplementedError(
                "QE f-orbital PDOS aggregation is not "
                "implemented yet."
            )

        for index, value in enumerate(
            parsed['ldos']
        ):
            totals[
                orbital
            ][index] += value

        for component_name, values in (
            parsed['components'].items()
        ):
            if orbital == 's' and component_name == 's':
                continue

            if component_name not in components:
                raise ValueError(
                    "Unexpected QE PDOS component: "
                    f"{component_name}"
                )

            for index, value in enumerate(
                values
            ):
                components[
                    component_name
                ][index] += value

    total_projected = []

    for index in range(
        len(energies)
    ):
        total_projected.append(
            totals['s'][index]
            + totals['p'][index]
            + totals['d'][index]
            + totals['f'][index]
        )

    return {
        'energies_ev': energies,
        's_total': totals['s'],
        'p_total': totals['p'],
        'pz': components['pz'],
        'px': components['px'],
        'py': components['py'],
        'd_total': totals['d'],
        'd3z2_r2': components['d3z2_r2'],
        'dxz': components['dxz'],
        'dyz': components['dyz'],
        'dx2_y2': components['dx2_y2'],
        'dxy': components['dxy'],
        'f_total': totals['f'],
        'total': total_projected,
        'files': pdos_files,
        'parsed_files': parsed_files,
        'npoints': len(energies),
    }

def run_scf(
    atoms,
    input_file,
    output_file,
    state_dir,
    pseudopotentials,
    pseudo_dir,
    cutoff_ev,
    kpoint_density=None,
    kpoint_size=(5, 5, 5),
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupation=None,
    parallel_cores=1,
    executable='pw.x',
    prefix='nanoworks',
):
    """Render, execute, and parse one QE pw.x SCF calculation."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    state_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    mesh = resolve_qe_kpoint_size(
        atoms,
        density=kpoint_density,
        size=kpoint_size,
    )

    occupation_settings = (
        resolve_qe_occupation(
            occupation
        )
    )

    input_text = render_scf_input(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=mesh,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        occupations=occupation_settings['occupations'],
        smearing=occupation_settings['smearing'],
        width_ev=occupation_settings['width_ev'],
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=state_dir,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_file.write_text(
        input_text,
        encoding='utf-8',
    )

    launcher = build_qe_launcher(
        parallel_cores=parallel_cores
    )

    execution = run_qe_program(
        input_file=input_file,
        output_file=output_file,
        executable=executable,
        launcher=launcher,
    )

    result = parse_pw_output(
        output_file
    )

    try:
        validate_qe_version(
            result['qe_version']
        )
    except ValueError as exc:
        raise RuntimeError(
            f"{exc} See '{output_file}'."
        ) from exc

    if not result['job_done']:
        raise RuntimeError(
            "Quantum ESPRESSO finished without a "
            "'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'kpoint_size': mesh,
        'execution': execution,
        'result': result,
    }

def run_nscf(
    atoms,
    input_file,
    output_file,
    state_dir,
    pseudopotentials,
    pseudo_dir,
    cutoff_ev,
    kpoint_density=None,
    kpoint_size=(5, 5, 5),
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupation=None,
    parallel_cores=1,
    executable='pw.x',
    prefix='nanoworks',
):
    """Render, execute, and parse one QE pw.x NSCF calculation."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE ground-state result is required "
            f"for the NSCF calculation: {state_dir}"
        )

    mesh = resolve_qe_kpoint_size(
        atoms,
        density=kpoint_density,
        size=kpoint_size,
    )

    occupation_settings = (
        resolve_qe_occupation(
            occupation
        )
    )

    input_text = render_nscf_input(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=mesh,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        occupations=occupation_settings['occupations'],
        smearing=occupation_settings['smearing'],
        width_ev=occupation_settings['width_ev'],
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=state_dir,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_file.write_text(
        input_text,
        encoding='utf-8',
    )

    launcher = build_qe_launcher(
        parallel_cores=parallel_cores
    )

    execution = run_qe_program(
        input_file=input_file,
        output_file=output_file,
        executable=executable,
        launcher=launcher,
    )

    result = parse_pw_output(
        output_file
    )

    try:
        validate_qe_version(
            result['qe_version']
        )
    except ValueError as exc:
        raise RuntimeError(
            f"{exc} See '{output_file}'."
        ) from exc

    if not result['job_done']:
        raise RuntimeError(
            "Quantum ESPRESSO NSCF calculation finished without a "
            "'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'kpoint_size': mesh,
        'execution': execution,
        'result': result,
    }

def run_bands(
    atoms,
    input_file,
    output_file,
    state_dir,
    pseudopotentials,
    pseudo_dir,
    cutoff_ev,
    band_path,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    occupation=None,
    parallel_cores=1,
    executable='pw.x',
    prefix='nanoworks',
):
    """Render, execute, and parse one QE pw.x bands calculation."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE ground-state result is required "
            f"for the bands calculation: {state_dir}"
        )

    occupation_settings = (
        resolve_qe_occupation(
            occupation
        )
    )

    input_text = render_bands_input(
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        band_path=band_path,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        occupations=occupation_settings['occupations'],
        smearing=occupation_settings['smearing'],
        width_ev=occupation_settings['width_ev'],
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=state_dir,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_file.write_text(
        input_text,
        encoding='utf-8',
    )

    launcher = build_qe_launcher(
        parallel_cores=parallel_cores
    )

    execution = run_qe_program(
        input_file=input_file,
        output_file=output_file,
        executable=executable,
        launcher=launcher,
    )

    result = parse_pw_output(
        output_file
    )

    try:
        validate_qe_version(
            result['qe_version']
        )
    except ValueError as exc:
        raise RuntimeError(
            f"{exc} See '{output_file}'."
        ) from exc

    if not result['job_done']:
        raise RuntimeError(
            "Quantum ESPRESSO bands calculation finished without a "
            "'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    bands = parse_pw_bands_output(
        output_file
    )

    requested_npoints = int(
        band_path['npoints']
    )

    if bands['nkpoints'] != requested_npoints:
        raise RuntimeError(
            "Quantum ESPRESSO bands output contains "
            f"{bands['nkpoints']} k-points, but "
            f"{requested_npoints} were requested."
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'band_path': band_path,
        'execution': execution,
        'result': result,
        'bands': bands,
    }

def run_dos(
    input_file,
    output_file,
    state_dir,
    dos_file,
    emin=None,
    emax=None,
    delta_e=None,
    bz_sum=None,
    degauss=None,
    ngauss=None,
    parallel_cores=1,
    executable='dos.x',
    prefix='nanoworks',
):
    """Render and execute one Quantum ESPRESSO dos.x calculation."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    dos_file = Path(
        dos_file
    )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE electronic state is required "
            f"for the DOS calculation: {state_dir}"
        )

    input_text = render_dos_input(
        prefix=prefix,
        outdir=state_dir,
        fildos=dos_file,
        bz_sum=bz_sum,
        emin=emin,
        emax=emax,
        delta_e=delta_e,
        degauss=degauss,
        ngauss=ngauss,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_file.write_text(
        input_text,
        encoding='utf-8',
    )

    launcher = build_qe_launcher(
        parallel_cores=parallel_cores
    )

    execution = run_qe_program(
        input_file=input_file,
        output_file=output_file,
        executable=executable,
        launcher=launcher,
    )

    text = output_file.read_text(
        encoding='utf-8',
        errors='replace',
    )

    if 'JOB DONE.' not in text:
        raise RuntimeError(
            "Quantum ESPRESSO DOS calculation finished without a "
            "'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    if not dos_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO dos.x finished but the DOS data file "
            f"was not created: {dos_file}"
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'dos_file': dos_file,
        'execution': execution,
    }

def run_projwfc(
    input_file,
    output_file,
    state_dir,
    pdos_prefix,
    emin=None,
    emax=None,
    delta_e=None,
    degauss=None,
    ngauss=None,
    parallel_cores=1,
    executable='projwfc.x',
    prefix='nanoworks',
):
    """Render and execute one Quantum ESPRESSO projwfc.x calculation."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    pdos_prefix = Path(
        pdos_prefix
    )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE electronic state is required "
            f"for the PDOS calculation: {state_dir}"
        )

    input_text = render_projwfc_input(
        prefix=prefix,
        outdir=state_dir,
        filpdos=pdos_prefix,
        emin=emin,
        emax=emax,
        delta_e=delta_e,
        degauss=degauss,
        ngauss=ngauss,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pdos_prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_file.write_text(
        input_text,
        encoding='utf-8',
    )

    launcher = build_qe_launcher(
        parallel_cores=parallel_cores
    )

    execution = run_qe_program(
        input_file=input_file,
        output_file=output_file,
        executable=executable,
        launcher=launcher,
    )

    text = output_file.read_text(
        encoding='utf-8',
        errors='replace',
    )

    if 'JOB DONE.' not in text:
        raise RuntimeError(
            "Quantum ESPRESSO PDOS calculation finished without a "
            "'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    pdos_tot_file = Path(
        str(pdos_prefix)
        + '.pdos_tot'
    )

    if not pdos_tot_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO projwfc.x finished but the "
            f"PDOS summary file was not created: {pdos_tot_file}"
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'pdos_prefix': pdos_prefix,
        'pdos_tot_file': pdos_tot_file,
        'execution': execution,
    }
