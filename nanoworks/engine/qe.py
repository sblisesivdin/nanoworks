"""Quantum ESPRESSO computation engine helpers."""

from ase.data import atomic_masses, atomic_numbers

QE_REFERENCE_VERSION = (7, 2)

# CODATA-compatible conversion used by ASE and QE-related workflows.
EV_PER_RYDBERG = 13.605693122994


def ev_to_rydberg(value):
    """Convert an energy value from electron-volts to Rydberg."""
    return float(value) / EV_PER_RYDBERG


def build_control_settings(
    calculation='scf',
    prefix='nanoworks',
    pseudo_dir=None,
    outdir=None,
):
    """Build the QE &CONTROL namelist settings."""
    settings = {
        'calculation': str(calculation).lower(),
        'prefix': str(prefix),
    }

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
    species = build_atomic_species(
        atoms,
        pseudopotentials,
    )

    control = build_control_settings(
        calculation='scf',
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
    kpoints = build_kpoint_settings(
        kpoint_size,
        gamma=gamma,
    )

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
        f"K_POINTS {kpoints['option']}",
    ])

    nk1, nk2, nk3 = kpoints['size']
    sk1, sk2, sk3 = kpoints['shift']

    lines.append(
        f"{nk1} {nk2} {nk3} {sk1} {sk2} {sk3}"
    )

    lines.extend([
        '',
        f"CELL_PARAMETERS {cell['option']}",
    ])

    for x, y, z in cell['vectors']:
        lines.append(
            f"{x:.12f} {y:.12f} {z:.12f}"
        )

    return "\n".join(lines) + "\n"
