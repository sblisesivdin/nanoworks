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

    job_done = (
        'JOB DONE.' in text
    )

    return {
        'job_done': job_done,
        'qe_version': qe_version,
        'total_energy_ry': total_energy_ry,
        'total_energy_ev': total_energy_ev,
        'fermi_energy_ev': fermi_energy_ev,
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
