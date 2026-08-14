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
