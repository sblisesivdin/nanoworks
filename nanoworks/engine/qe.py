"""Quantum ESPRESSO computation engine helpers."""

import math
import operator
import re
import os
import shutil
import subprocess
from pathlib import Path
import numpy as np
from ase.units import Bohr
from ase.data import atomic_masses, atomic_numbers
from ase.calculators.calculator import kptdensity2monkhorstpack
from nanoworks.pseudos import (
    read_upf_atomic_manifolds,
    read_upf_z_valence,
)

QE_REFERENCE_VERSION = (7, 2)

# CODATA-compatible conversion used by ASE and QE-related workflows.
EV_PER_RYDBERG = 13.605693122994
THZ_PER_CM_MINUS_ONE = 0.0299792458
BOLTZMANN_EV_PER_K = 8.617333262145e-5
EV_PER_THZ = 4.135667696e-3
KJ_PER_MOL_PER_EV = 96.48533212331002


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
    xc_calc=None,
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
    exx_qpoint_grid=None,
):
    """Build the basic QE &SYSTEM namelist settings."""
    settings = {
        'ibrav': 0,
        'nat': int(nat),
        'ntyp': int(ntyp),
        'ecutwfc': ev_to_rydberg(cutoff_ev),
    }

    xc_settings = None

    if xc_calc is not None:
        xc_settings = resolve_qe_xc_settings(
            xc_calc=xc_calc,
            pseudo_xc=pseudo_xc,
            exx_fraction=exx_fraction,
            omega=omega,
        )

        settings['input_dft'] = xc_settings[
            'input_dft'
        ]

        if xc_settings['hybrid']:
            settings['exx_fraction'] = xc_settings[
                'exx_fraction'
            ]

            if (
                xc_settings['screening_parameter']
                is not None
            ):
                settings['screening_parameter'] = (
                    xc_settings[
                        'screening_parameter'
                    ]
                )

    if exx_qpoint_grid is not None:
        if (
            xc_settings is None
            or not xc_settings['hybrid']
        ):
            raise ValueError(
                "QE EXX q-point grids require a hybrid functional."
            )

        nq1, nq2, nq3 = (
            _normalize_qe_exx_qpoint_grid(
                exx_qpoint_grid
            )
        )

        settings.update({
            'nqx1': nq1,
            'nqx2': nq2,
            'nqx3': nq3,
        })

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

def build_qe_magnetic_species(
    atoms,
    pseudopotentials,
    pseudo_dir,
    magnetic_moments,
):
    """Build QE species and positions for collinear magnetism."""
    chemical_symbols = (
        atoms.get_chemical_symbols()
    )

    moments = [
        float(value)
        for value in magnetic_moments
    ]

    if len(moments) != len(atoms):
        raise ValueError(
            "QE magnetic moment count must match "
            "the number of atoms."
        )

    if not any(
        abs(moment) > 1.0e-12
        for moment in moments
    ):
        raise ValueError(
            "QE spin-polarized calculations require at least "
            "one non-zero initial magnetic moment."
        )

    missing = sorted({
        symbol
        for symbol in chemical_symbols
        if symbol not in pseudopotentials
    })

    if missing:
        raise ValueError(
            "Missing pseudopotential mapping for: "
            + ", ".join(missing)
        )

    moment_groups = {}
    atom_group_keys = []

    for symbol, moment in zip(
        chemical_symbols,
        moments,
    ):
        groups = moment_groups.setdefault(
            symbol,
            [],
        )

        group_index = None

        for index, existing_moment in enumerate(
            groups
        ):
            if abs(
                existing_moment
                - moment
            ) <= 1.0e-12:
                group_index = index
                break

        if group_index is None:
            groups.append(
                moment
            )

            group_index = (
                len(groups) - 1
            )

        atom_group_keys.append(
            (
                symbol,
                group_index,
            )
        )

    labels = {}

    for symbol, groups in moment_groups.items():
        if len(groups) == 1:
            labels[
                (
                    symbol,
                    0,
                )
            ] = symbol

            continue

        if len(groups) > 9:
            raise ValueError(
                "QE magnetic species labeling supports at most "
                f"nine distinct moments for element {symbol}."
            )

        for group_index in range(
            len(groups)
        ):
            label = (
                f"{symbol}"
                f"{group_index + 1}"
            )

            if len(label) > 3:
                raise ValueError(
                    "QE atomic species labels must not exceed "
                    f"three characters: {label}"
                )

            labels[
                (
                    symbol,
                    group_index,
                )
            ] = label

    z_valence_by_symbol = {}

    for symbol in moment_groups:
        pseudo_file = Path(
            pseudopotentials[symbol]
        )

        if not pseudo_file.is_absolute():
            pseudo_file = (
                Path(pseudo_dir)
                / pseudo_file
            )

        z_valence_by_symbol[
            symbol
        ] = read_upf_z_valence(
            pseudo_file
        )

    species = []
    starting_magnetizations = {}
    seen_group_keys = set()

    for group_key in atom_group_keys:
        if group_key in seen_group_keys:
            continue

        seen_group_keys.add(
            group_key
        )

        symbol, group_index = group_key

        moment = moment_groups[
            symbol
        ][group_index]

        fraction = (
            moment
            / z_valence_by_symbol[symbol]
        )

        if abs(fraction) > 1.0:
            raise ValueError(
                "QE 7.2 starting magnetization must be "
                "between -1 and 1. The requested moment is "
                f"too large for species {symbol}."
            )

        atomic_number = atomic_numbers[
            symbol
        ]

        species.append(
            (
                labels[group_key],
                float(
                    atomic_masses[
                        atomic_number
                    ]
                ),
                str(
                    pseudopotentials[
                        symbol
                    ]
                ),
            )
        )

        species_index = len(
            species
        )

        starting_magnetizations[
            f'starting_magnetization({species_index})'
        ] = fraction

    if len(species) > 10:
        raise ValueError(
            "QE 7.2 supports at most ten atomic species."
        )

    positions = []

    for group_key, position in zip(
        atom_group_keys,
        atoms.get_positions(),
    ):
        positions.append(
            (
                labels[group_key],
                float(position[0]),
                float(position[1]),
                float(position[2]),
            )
        )

    return {
        'species': species,
        'positions': {
            'option': 'angstrom',
            'positions': positions,
        },
        'starting_magnetizations': (
            starting_magnetizations
        ),
        'species_labels': [
            labels[group_key]
            for group_key in atom_group_keys
        ],
        'magnetic_moments': moments,
        'ntyp': len(species),
    }

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

def _normalize_qe_exx_qpoint_grid(qpoint_grid):
    """Validate and normalize a QE EXX q-point grid."""
    try:
        raw_grid = tuple(qpoint_grid)
    except TypeError as error:
        raise ValueError(
            "QE EXX q-point grid must contain exactly 3 values."
        ) from error

    if len(raw_grid) != 3:
        raise ValueError(
            "QE EXX q-point grid must contain exactly 3 values."
        )

    try:
        grid = tuple(
            int(value)
            for value in raw_grid
        )
    except (TypeError, ValueError) as error:
        raise ValueError(
            "QE EXX q-point grid values must be positive integers."
        ) from error

    if any(
        float(value) != parsed
        for value, parsed in zip(raw_grid, grid)
    ):
        raise ValueError(
            "QE EXX q-point grid values must be positive integers."
        )

    if any(value <= 0 for value in grid):
        raise ValueError(
            "QE EXX q-point grid values must be positive integers."
        )

    return grid


def build_qe_exx_additional_kpoints(
    band_path,
    qpoint_grid,
):
    """Build zero-weight band and k+q points for a hybrid SCF run."""
    band_kpoints = list(
        band_path.get('kpoints', [])
    )

    if not band_kpoints:
        raise ValueError(
            "QE hybrid band path does not contain any k-points."
        )

    declared_npoints = band_path.get(
        'npoints'
    )

    if (
        declared_npoints is not None
        and int(declared_npoints) != len(band_kpoints)
    ):
        raise ValueError(
            "QE hybrid band path point count does not match "
            "its metadata."
        )

    grid = _normalize_qe_exx_qpoint_grid(
        qpoint_grid
    )

    normalized_band_kpoints = []

    for kpoint in band_kpoints:
        if len(kpoint) != 3:
            raise ValueError(
                "Each QE hybrid band k-point must contain exactly "
                "3 coordinates."
            )

        coordinates = tuple(
            float(value)
            for value in kpoint
        )

        if not all(
            math.isfinite(value)
            for value in coordinates
        ):
            raise ValueError(
                "QE hybrid band k-point coordinates must be finite."
            )

        normalized_band_kpoints.append(
            coordinates
        )

    additional_kpoints = list(
        normalized_band_kpoints
    )
    seen = {
        tuple(round(value, 12) for value in kpoint)
        for kpoint in normalized_band_kpoints
    }

    nq1, nq2, nq3 = grid

    for kpoint in normalized_band_kpoints:
        for iq1 in range(nq1):
            for iq2 in range(nq2):
                for iq3 in range(nq3):
                    if iq1 == iq2 == iq3 == 0:
                        continue

                    shifted = tuple(
                        value % 1.0
                        for value in (
                            kpoint[0] + iq1 / nq1,
                            kpoint[1] + iq2 / nq2,
                            kpoint[2] + iq3 / nq3,
                        )
                    )
                    key = tuple(
                        round(value, 12)
                        for value in shifted
                    )

                    if key in seen:
                        continue

                    seen.add(
                        key
                    )
                    additional_kpoints.append(
                        shifted
                    )

    band_point_count = len(
        normalized_band_kpoints
    )

    return {
        'option': 'crystal',
        'kpoints': additional_kpoints,
        'npoints': len(additional_kpoints),
        'band_indices': list(
            range(band_point_count)
        ),
        'helper_indices': list(
            range(
                band_point_count,
                len(additional_kpoints),
            )
        ),
        'qpoint_grid': grid,
    }

def render_qe_exx_additional_kpoints(settings):
    """Render a QE ADDITIONAL_K_POINTS card for a hybrid SCF run."""
    option = str(
        settings.get('option', '')
    ).strip().lower()

    if option != 'crystal':
        raise ValueError(
            "QE EXX additional k-points must use crystal coordinates."
        )

    kpoints = list(
        settings.get('kpoints', [])
    )

    if not kpoints:
        raise ValueError(
            "QE EXX additional k-point list must not be empty."
        )

    declared_npoints = settings.get(
        'npoints'
    )

    if (
        declared_npoints is not None
        and int(declared_npoints) != len(kpoints)
    ):
        raise ValueError(
            "QE EXX additional k-point count does not match "
            "its metadata."
        )

    lines = [
        'ADDITIONAL_K_POINTS crystal',
        str(len(kpoints)),
    ]

    for kpoint in kpoints:
        if len(kpoint) != 3:
            raise ValueError(
                "Each QE EXX additional k-point must contain exactly "
                "3 coordinates."
            )

        x, y, z = (
            float(value)
            for value in kpoint
        )

        if not all(
            math.isfinite(value)
            for value in (x, y, z)
        ):
            raise ValueError(
                "QE EXX additional k-point coordinates must be finite."
            )

        lines.append(
            f"{x:.12f} {y:.12f} {z:.12f} 0.0"
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

def resolve_qe_xc_settings(
    xc_calc,
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
):
    """Resolve Nanoworks XC names to QE input_dft hybrid settings."""
    xc = str(
        xc_calc
    ).strip().lower().replace(
        '_',
        '-',
    )
    pseudo_xc = str(
        pseudo_xc
    ).strip().lower()

    aliases = {
        'pbe': {
            'name': 'PBE',
            'input_dft': 'PBE',
            'hybrid': False,
            'exx_fraction': None,
            'screening_parameter': None,
        },
        'hse': {
            'name': 'HSE06',
            'input_dft': 'HSE',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': 0.106,
        },
        'hse06': {
            'name': 'HSE06',
            'input_dft': 'HSE',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': 0.106,
        },
        'hse-06': {
            'name': 'HSE06',
            'input_dft': 'HSE',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': 0.106,
        },
        'hse03': {
            'name': 'HSE03',
            'input_dft': 'HSE',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': 0.15,
        },
        'hse-03': {
            'name': 'HSE03',
            'input_dft': 'HSE',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': 0.15,
        },
        'pbe0': {
            'name': 'PBE0',
            'input_dft': 'PBE0',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': None,
        },
        'pbe-0': {
            'name': 'PBE0',
            'input_dft': 'PBE0',
            'hybrid': True,
            'exx_fraction': 0.25,
            'screening_parameter': None,
        },
    }

    try:
        settings = dict(
            aliases[xc]
        )
    except KeyError:
        raise ValueError(
            "The current Nanoworks QE backend supports PBE, "
            "HSE06, HSE03, and PBE0. "
            f"Requested XC: {xc_calc}"
        )

    if pseudo_xc != 'pbe':
        raise ValueError(
            f"QE XC '{xc_calc}' requires the installed "
            f"PBE pseudopotentials, not '{pseudo_xc}'."
        )

    if not settings['hybrid']:
        if exx_fraction is not None or omega is not None:
            raise ValueError(
                "XC_exx_fraction and XC_omega can only be used "
                "with a QE hybrid functional."
            )

        settings['pseudo_xc'] = pseudo_xc
        return settings

    if exx_fraction is not None:
        if isinstance(exx_fraction, bool):
            raise TypeError(
                "QE exact-exchange fraction must be a real number."
            )

        try:
            exx_fraction = float(
                exx_fraction
            )
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "QE exact-exchange fraction must be a real number."
            ) from exc

        if (
            not math.isfinite(exx_fraction)
            or exx_fraction <= 0.0
            or exx_fraction > 1.0
        ):
            raise ValueError(
                "QE exact-exchange fraction must satisfy "
                "0 < XC_exx_fraction <= 1."
            )

        settings['exx_fraction'] = exx_fraction

    if omega is not None:
        if settings['name'] == 'PBE0':
            raise ValueError(
                "XC_omega is only valid for screened HSE functionals."
            )

        if isinstance(omega, bool):
            raise TypeError(
                "QE HSE screening parameter must be a real number."
            )

        try:
            omega = float(
                omega
            )
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "QE HSE screening parameter must be a real number."
            ) from exc

        if not math.isfinite(omega) or omega <= 0.0:
            raise ValueError(
                "QE HSE screening parameter must be positive."
            )

        settings['screening_parameter'] = omega

    settings['pseudo_xc'] = pseudo_xc
    return settings


def validate_qe_xc(
    xc_calc,
    pseudo_xc='pbe',
    allow_hybrid=False,
):
    """Validate XC compatibility with the installed QE pseudo set."""
    settings = resolve_qe_xc_settings(
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
    )

    if settings['hybrid'] and not allow_hybrid:
        raise ValueError(
            f"QE hybrid XC '{settings['name']}' is not enabled "
            "for this calculation stage yet."
        )

    return settings['name'].lower()

def resolve_qe_hubbard(
    setup_params,
    pseudopotentials,
    pseudo_dir,
):
    """Translate Nanoworks Setup_params to QE Hubbard settings."""
    if not setup_params:
        return None

    if not isinstance(
        setup_params,
        dict,
    ):
        raise TypeError(
            "QE Setup_params must be a dictionary."
        )

    parameters = []
    notices = []

    for symbol, specification in setup_params.items():
        if symbol not in pseudopotentials:
            raise ValueError(
                "QE Hubbard parameters refer to an element "
                f"without a pseudopotential: {symbol}"
            )

        if not isinstance(
            specification,
            str,
        ):
            raise TypeError(
                "QE Hubbard setup specifications must be strings."
            )

        fields = [
            field.strip()
            for field in specification.lstrip(':').split(',')
        ]

        if len(fields) not in {
            2,
            3,
        }:
            raise ValueError(
                "QE Hubbard setup must use "
                "':orbital,U' or ':orbital,U,flag': "
                f"{symbol}={specification!r}"
            )

        orbital = fields[0].lower()

        orbital_match = re.fullmatch(
            r'(\d+)?([spdfg])',
            orbital,
        )

        if orbital_match is None:
            raise ValueError(
                "Unsupported QE Hubbard orbital specification: "
                f"{symbol}={orbital!r}"
            )

        try:
            value_ev = float(
                fields[1]
            )
        except ValueError:
            raise ValueError(
                "QE Hubbard U must be a numeric value: "
                f"{symbol}={fields[1]!r}"
            )

        if value_ev <= 0.0:
            raise ValueError(
                "QE Hubbard U must be greater than zero: "
                f"{symbol}={value_ev}"
            )

        pseudo_file = Path(
            pseudopotentials[
                symbol
            ]
        )

        if not pseudo_file.is_absolute():
            if pseudo_dir is None:
                raise ValueError(
                    "QE Hubbard manifold resolution requires "
                    "a pseudopotential directory."
                )

            pseudo_file = (
                Path(pseudo_dir)
                / pseudo_file
            )

        available_manifolds = (
            read_upf_atomic_manifolds(
                pseudo_file
            )
        )

        principal_number = (
            orbital_match.group(1)
        )

        angular_orbital = (
            orbital_match.group(2)
        )

        if principal_number is not None:
            manifold = (
                principal_number
                + angular_orbital
            )

            if manifold not in available_manifolds:
                raise ValueError(
                    f"QE Hubbard manifold {symbol}-{manifold} "
                    "was not found in pseudopotential "
                    f"'{pseudo_file}'. Available manifolds: "
                    + ", ".join(
                        available_manifolds
                    )
                )

        else:
            candidates = [
                manifold
                for manifold in available_manifolds
                if manifold.endswith(
                    angular_orbital
                )
            ]

            if not candidates:
                raise ValueError(
                    f"QE Hubbard orbital {symbol}-{angular_orbital} "
                    "was not found in pseudopotential "
                    f"'{pseudo_file}'. Available manifolds: "
                    + ", ".join(
                        available_manifolds
                    )
                )

            if len(candidates) > 1:
                raise ValueError(
                    f"QE Hubbard orbital {symbol}-{angular_orbital} "
                    "matches multiple pseudopotential manifolds: "
                    + ", ".join(
                        candidates
                    )
                    + ". Specify the principal quantum number "
                    "explicitly."
                )

            manifold = candidates[0]

        if len(fields) == 3:
            notices.append(
                "The GPAW-specific Hubbard normalization flag "
                f"for {symbol}-{manifold} cannot be represented "
                "exactly in QE; ortho-atomic projectors are used."
            )

        parameters.append({
            'symbol': symbol,
            'manifold': manifold,
            'value_ev': value_ev,
        })

    return {
        'projector': 'ortho-atomic',
        'parameters': parameters,
        'notices': notices,
    }


def render_qe_hubbard_card(
    settings,
    species_by_element=None,
):
    """Render a QE 7.2 HUBBARD card."""
    if settings is None:
        return ''

    species_by_element = (
        species_by_element
        or {}
    )

    lines = [
        f"! NOTICE: {notice}"
        for notice in settings[
            'notices'
        ]
    ]

    lines.append(
        "HUBBARD "
        f"({settings['projector']})"
    )

    for parameter in settings[
        'parameters'
    ]:
        symbol = parameter[
            'symbol'
        ]

        labels = species_by_element.get(
            symbol,
            [
                symbol,
            ],
        )

        for label in labels:
            lines.append(
                "U "
                f"{label}-{parameter['manifold']} "
                f"{parameter['value_ev']:.12g}"
            )

    return "\n".join(
        lines
    )

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

def _qe_cell_axes_are_orthogonal(
    atoms,
    tolerance=1.0e-6,
):
    """Return whether all three cell angles are approximately 90 degrees."""
    angles = atoms.cell.angles()

    return all(
        abs(
            float(angle)
            - 90.0
        ) <= tolerance
        for angle in angles
    )

def resolve_qe_cell_dofree(
    relax_cell,
):
    """Map a Nanoworks strain mask to QE cell_dofree."""
    mask = tuple(
        relax_cell
    )

    if len(mask) != 6:
        raise ValueError(
            "QE Relax_cell must contain exactly six "
            "boolean components."
        )

    if not all(
        isinstance(value, bool)
        for value in mask
    ):
        raise TypeError(
            "QE Relax_cell components must be boolean values."
        )

    mappings = {
        (
            False,
            False,
            False,
            False,
            False,
            False,
        ): None,
        (
            True,
            False,
            False,
            False,
            False,
            False,
        ): 'x',
        (
            False,
            True,
            False,
            False,
            False,
            False,
        ): 'y',
        (
            False,
            False,
            True,
            False,
            False,
            False,
        ): 'z',
        (
            True,
            True,
            False,
            False,
            False,
            False,
        ): 'xy',
        (
            True,
            False,
            True,
            False,
            False,
            False,
        ): 'xz',
        (
            False,
            True,
            True,
            False,
            False,
            False,
        ): 'yz',
        (
            True,
            True,
            True,
            False,
            False,
            False,
        ): 'xyz',
        (
            True,
            True,
            False,
            False,
            False,
            True,
        ): '2Dxy',
        (
            True,
            True,
            True,
            True,
            True,
            True,
        ): 'all',
    }

    try:
        return mappings[
            mask
        ]
    except KeyError:
        raise NotImplementedError(
            "The requested Relax_cell mask cannot be "
            "represented safely by QE cell_dofree."
        )


def resolve_qe_relaxation_settings(
    optimizer,
    max_force,
    max_step,
    relax_cell,
    hydrostatic_pressure=0.0,
    fix_symmetry=False,
    atoms=None,
):
    """Resolve Nanoworks geometry settings to QE namelists."""
    optimizer_key = (
        str(optimizer)
        .strip()
        .lower()
    )

    optimizer_mappings = {
        'quasinewton': 'bfgs',
        'lbfgs': 'bfgs',
        'bfgs': 'bfgs',
    }

    try:
        ion_dynamics = optimizer_mappings[
            optimizer_key
        ]
    except KeyError:
        raise NotImplementedError(
            "QE geometry optimization currently supports "
            "QuasiNewton and LBFGS only."
        )

    max_force = float(
        max_force
    )

    if max_force <= 0.0:
        raise ValueError(
            "QE geometry force tolerance must be greater than zero."
        )

    max_step = float(
        max_step
    )

    if max_step <= 0.0:
        raise ValueError(
            "QE geometry maximum step must be greater than zero."
        )

    hydrostatic_pressure = float(
        hydrostatic_pressure
    )

    cell_dofree = resolve_qe_cell_dofree(
        relax_cell
    )

    notices = []

    if (
        atoms is not None
        and cell_dofree == 'xyz'
        and not _qe_cell_axes_are_orthogonal(
            atoms
        )
    ):
        cell_dofree = 'all'

        if fix_symmetry:
            notices.append(
                "The normal-strain Relax_cell mask is mapped "
                "to cell_dofree='all' because QE's 'xyz' mode "
                "is unsafe for non-orthogonal cell axes; active "
                "crystal symmetry constrains compatible cell changes."
            )
        else:
            notices.append(
                "The normal-strain Relax_cell mask is mapped "
                "to cell_dofree='all' because QE's 'xyz' mode "
                "is unsafe for non-orthogonal cell axes; this "
                "also enables shear degrees of freedom."
            )

    if (
        cell_dofree is None
        and hydrostatic_pressure != 0.0
    ):
        raise ValueError(
            "QE hydrostatic pressure requires cell relaxation."
        )

    trust_radius_max = (
        max_step
        / Bohr
    )

    settings = {
        'calculation': (
            'vc-relax'
            if cell_dofree is not None
            else 'relax'
        ),
        'notices': notices,
        'control': {
            'forc_conv_thr': ev_to_rydberg(
                max_force
                * Bohr
            ),
        },
        'system': {
            'nosym': not bool(
                fix_symmetry
            ),
        },
        'ions': {
            'ion_dynamics': ion_dynamics,
            'trust_radius_max': trust_radius_max,
            'trust_radius_ini': min(
                0.5,
                trust_radius_max,
            ),
        },
        'cell': None,
    }

    if cell_dofree is not None:
        settings['cell'] = {
            'cell_dynamics': 'bfgs',
            'cell_dofree': cell_dofree,
            'press': (
                hydrostatic_pressure
                * 10.0
            ),
        }

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


def resolve_qe_phonon_qpoint_grid(supercell):
    """Map a diagonal phonon supercell to a QE DFPT q-point grid."""
    if isinstance(supercell, (str, bytes)):
        raise TypeError(
            "Phonon_supercell must be a three-value sequence or "
            "a 3x3 diagonal matrix."
        )

    try:
        outer = tuple(supercell)
    except TypeError as exc:
        raise TypeError(
            "Phonon_supercell must be a three-value sequence or "
            "a 3x3 diagonal matrix."
        ) from exc

    if len(outer) != 3:
        raise ValueError(
            "Phonon_supercell must contain exactly three values "
            "or three matrix rows."
        )

    scalar_input = all(
        not hasattr(value, '__iter__')
        or isinstance(value, (str, bytes))
        for value in outer
    )

    if scalar_input:
        diagonal = outer
    else:
        rows = []

        for row in outer:
            if isinstance(row, (str, bytes)):
                raise TypeError(
                    "Phonon_supercell matrix rows must contain "
                    "three integers."
                )

            try:
                row = tuple(row)
            except TypeError as exc:
                raise TypeError(
                    "Phonon_supercell must not mix scalar values "
                    "and matrix rows."
                ) from exc

            if len(row) != 3:
                raise ValueError(
                    "Phonon_supercell matrix must have shape 3x3."
                )

            rows.append(row)

        for row_index, row in enumerate(rows):
            for column_index, value in enumerate(row):
                if isinstance(value, bool):
                    raise TypeError(
                        "Phonon_supercell matrix values must be integers."
                    )

                try:
                    integer_value = operator.index(value)
                except TypeError as exc:
                    raise TypeError(
                        "Phonon_supercell matrix values must be integers."
                    ) from exc

                if (
                    row_index != column_index
                    and integer_value != 0
                ):
                    raise ValueError(
                        "Native QE phonons currently require a diagonal "
                        "Phonon_supercell matrix; non-diagonal matrices "
                        "cannot be mapped safely to nq1, nq2, nq3."
                    )

        diagonal = tuple(
            rows[index][index]
            for index in range(3)
        )

    if any(isinstance(value, bool) for value in diagonal):
        raise TypeError(
            "Phonon_supercell diagonal values must be integers."
        )

    try:
        grid = tuple(
            operator.index(value)
            for value in diagonal
        )
    except TypeError as exc:
        raise TypeError(
            "Phonon_supercell diagonal values must be integers."
        ) from exc

    if any(value <= 0 for value in grid):
        raise ValueError(
            "Phonon_supercell diagonal values must be positive integers."
        )

    return grid


def build_ph_settings(
    prefix,
    outdir,
    fildyn,
    qpoint_grid,
    tr2_ph=1.0e-12,
):
    """Build the QE ph.x &INPUTPH namelist settings."""
    for name, value in (
        ('prefix', prefix),
        ('outdir', outdir),
        ('fildyn', fildyn),
    ):
        if value is None or not str(value).strip():
            raise ValueError(
                f"QE phonon {name} must not be empty."
            )

    prefix = str(prefix).strip()
    outdir = str(outdir).strip()
    fildyn = str(fildyn).strip()

    try:
        qpoint_grid = tuple(qpoint_grid)
    except TypeError as exc:
        raise TypeError(
            "QE phonon q-point grid must be an iterable "
            "of three positive integers."
        ) from exc

    if len(qpoint_grid) != 3:
        raise ValueError(
            "QE phonon q-point grid must contain exactly 3 values."
        )

    if any(isinstance(value, bool) for value in qpoint_grid):
        raise TypeError(
            "QE phonon q-point grid values must be integers."
        )

    try:
        qpoint_grid = tuple(
            operator.index(value)
            for value in qpoint_grid
        )
    except TypeError as exc:
        raise TypeError(
            "QE phonon q-point grid values must be integers."
        ) from exc

    if any(value <= 0 for value in qpoint_grid):
        raise ValueError(
            "QE phonon q-point grid values must be positive integers."
        )

    try:
        tr2_ph = float(tr2_ph)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            "QE phonon convergence threshold must be a real number."
        ) from exc

    if not math.isfinite(tr2_ph) or tr2_ph <= 0.0:
        raise ValueError(
            "QE phonon convergence threshold must be a positive "
            "finite value."
        )

    nq1, nq2, nq3 = qpoint_grid

    return {
        'prefix': prefix,
        'outdir': outdir,
        'fildyn': fildyn,
        'tr2_ph': tr2_ph,
        'ldisp': True,
        'nq1': nq1,
        'nq2': nq2,
        'nq3': nq3,
    }


def render_ph_input(
    prefix,
    outdir,
    fildyn,
    qpoint_grid,
    tr2_ph=1.0e-12,
    title='Nanoworks native QE phonon calculation',
):
    """Render a complete QE ph.x input for a regular q-point grid."""
    if title is None:
        raise ValueError(
            "QE phonon title must be one non-empty line."
        )

    title = str(title).strip()

    if not title or '\n' in title or '\r' in title:
        raise ValueError(
            "QE phonon title must be one non-empty line."
        )

    settings = build_ph_settings(
        prefix=prefix,
        outdir=outdir,
        fildyn=fildyn,
        qpoint_grid=qpoint_grid,
        tr2_ph=tr2_ph,
    )

    return "\n".join([
        title,
        render_namelist(
            'INPUTPH',
            settings,
        ),
        '',
    ])


def build_q2r_settings(
    fildyn,
    flfrc,
    zasr='no',
):
    """Build the QE q2r.x &INPUT namelist settings."""
    for name, value in (
        ('fildyn', fildyn),
        ('flfrc', flfrc),
    ):
        if value is None or not str(value).strip():
            raise ValueError(
                f"QE q2r {name} must not be empty."
            )

    fildyn = str(fildyn).strip()
    flfrc = str(flfrc).strip()

    if zasr is None:
        raise ValueError(
            "QE q2r zasr must not be empty."
        )

    zasr = str(zasr).strip().lower()

    allowed_zasr = {
        'no',
        'simple',
        'crystal',
        'one-dim',
        'zero-dim',
    }

    if zasr not in allowed_zasr:
        raise ValueError(
            "Unsupported QE q2r zasr setting: "
            f"{zasr!r}. Supported values are: "
            + ", ".join(sorted(allowed_zasr))
            + "."
        )

    return {
        'fildyn': fildyn,
        'flfrc': flfrc,
        'zasr': zasr,
    }


def render_q2r_input(
    fildyn,
    flfrc,
    zasr='no',
):
    """Render a complete QE q2r.x input."""
    settings = build_q2r_settings(
        fildyn=fildyn,
        flfrc=flfrc,
        zasr=zasr,
    )

    return (
        render_namelist(
            'INPUT',
            settings,
        )
        + '\n'
    )


def build_matdyn_band_settings(
    flfrc,
    flfrq,
    acoustic_sum_rule=True,
):
    """Build QE matdyn.x settings for a phonon band path."""
    for name, value in (
        ('flfrc', flfrc),
        ('flfrq', flfrq),
    ):
        if value is None or not str(value).strip():
            raise ValueError(
                f"QE matdyn {name} must not be empty."
            )

    if not isinstance(acoustic_sum_rule, bool):
        raise TypeError(
            "QE matdyn acoustic sum rule setting must be boolean."
        )

    return {
        'flfrc': str(flfrc).strip(),
        'asr': (
            'crystal'
            if acoustic_sum_rule
            else 'no'
        ),
        'dos': False,
        'flfrq': str(flfrq).strip(),
        'q_in_band_form': False,
        'q_in_cryst_coord': True,
    }


def render_matdyn_qpoints(band_path):
    """Render explicit crystalline q-points for matdyn.x."""
    if not isinstance(band_path, dict):
        raise TypeError(
            "QE matdyn band path must be a mapping."
        )

    option = str(
        band_path.get('option', '')
    ).strip().lower()

    if option != 'crystal':
        raise ValueError(
            "QE matdyn band path must use crystal coordinates."
        )

    qpoints = list(
        band_path.get('kpoints', [])
    )

    if not qpoints:
        raise ValueError(
            "QE matdyn band path does not contain any q-points."
        )

    declared_npoints = band_path.get('npoints')

    if (
        declared_npoints is not None
        and int(declared_npoints) != len(qpoints)
    ):
        raise ValueError(
            "QE matdyn band path point count does not match "
            "its metadata."
        )

    lines = [
        str(len(qpoints)),
    ]

    for qpoint in qpoints:
        if len(qpoint) != 3:
            raise ValueError(
                "Each QE matdyn q-point must contain exactly "
                "3 coordinates."
            )

        try:
            x, y, z = (
                float(value)
                for value in qpoint
            )
        except (TypeError, ValueError) as exc:
            raise TypeError(
                "QE matdyn q-point coordinates must be real numbers."
            ) from exc

        if not all(
            math.isfinite(value)
            for value in (x, y, z)
        ):
            raise ValueError(
                "QE matdyn q-point coordinates must be finite."
            )

        lines.append(
            f"{x:.12f} {y:.12f} {z:.12f}"
        )

    return "\n".join(lines)


def render_matdyn_band_input(
    flfrc,
    flfrq,
    band_path,
    acoustic_sum_rule=True,
):
    """Render a QE matdyn.x input for phonon dispersion."""
    settings = build_matdyn_band_settings(
        flfrc=flfrc,
        flfrq=flfrq,
        acoustic_sum_rule=acoustic_sum_rule,
    )

    qpoint_card = render_matdyn_qpoints(
        band_path
    )

    return "\n".join([
        render_namelist(
            'INPUT',
            settings,
        ),
        qpoint_card,
        '',
    ])


def build_matdyn_dos_settings(
    flfrc,
    fldos,
    qpoint_grid,
    acoustic_sum_rule=True,
):
    """Build QE matdyn.x settings for phonon DOS."""
    for name, value in (
        ('flfrc', flfrc),
        ('fldos', fldos),
    ):
        if value is None or not str(value).strip():
            raise ValueError(
                f"QE matdyn {name} must not be empty."
            )

    if not isinstance(acoustic_sum_rule, bool):
        raise TypeError(
            "QE matdyn acoustic sum rule setting must be boolean."
        )

    try:
        qpoint_grid = tuple(qpoint_grid)
    except TypeError as exc:
        raise TypeError(
            "QE matdyn DOS q-point grid must be an iterable "
            "of three positive integers."
        ) from exc

    if len(qpoint_grid) != 3:
        raise ValueError(
            "QE matdyn DOS q-point grid must contain exactly "
            "3 values."
        )

    if any(isinstance(value, bool) for value in qpoint_grid):
        raise TypeError(
            "QE matdyn DOS q-point grid values must be integers."
        )

    try:
        qpoint_grid = tuple(
            operator.index(value)
            for value in qpoint_grid
        )
    except TypeError as exc:
        raise TypeError(
            "QE matdyn DOS q-point grid values must be integers."
        ) from exc

    if any(value <= 0 for value in qpoint_grid):
        raise ValueError(
            "QE matdyn DOS q-point grid values must be positive "
            "integers."
        )

    nk1, nk2, nk3 = qpoint_grid

    return {
        'flfrc': str(flfrc).strip(),
        'asr': (
            'crystal'
            if acoustic_sum_rule
            else 'no'
        ),
        'dos': True,
        'nk1': nk1,
        'nk2': nk2,
        'nk3': nk3,
        'fldos': str(fldos).strip(),
    }


def render_matdyn_dos_input(
    flfrc,
    fldos,
    qpoint_grid,
    acoustic_sum_rule=True,
):
    """Render a QE matdyn.x input for phonon DOS."""
    settings = build_matdyn_dos_settings(
        flfrc=flfrc,
        fldos=fldos,
        qpoint_grid=qpoint_grid,
        acoustic_sum_rule=acoustic_sum_rule,
    )

    return (
        render_namelist(
            'INPUT',
            settings,
        )
        + '\n'
    )

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
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
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
    exx_additional_kpoints=None,
    relaxation_settings=None,
):
    """Render a complete QE pw.x input."""

    calculation = str(
        calculation
    ).strip().lower()

    allowed_calculations = {
        'scf',
        'nscf',
        'bands',
        'relax',
        'vc-relax',
    }

    if calculation not in allowed_calculations:
        raise ValueError(
            f"Unsupported QE pw.x calculation type: {calculation}"
        )

    xc_settings = resolve_qe_xc_settings(
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
    )

    if (
        xc_settings['hybrid']
        and calculation in {
            'nscf',
            'bands',
        }
    ):
        raise NotImplementedError(
            "QE does not support separate NSCF or bands "
            "calculations with hybrid functionals."
        )

    exx_additional_card = None
    exx_qpoint_grid = None

    if exx_additional_kpoints is not None:
        if calculation != 'scf':
            raise ValueError(
                "QE EXX additional k-points can only be used "
                "with calculation='scf'."
            )

        if not xc_settings['hybrid']:
            raise ValueError(
                "QE EXX additional k-points require a hybrid "
                "functional."
            )

        exx_additional_card = (
            render_qe_exx_additional_kpoints(
                exx_additional_kpoints
            )
        )
        exx_qpoint_grid = (
            exx_additional_kpoints.get(
                'qpoint_grid'
            )
        )

    relaxation_calculations = {
        'relax',
        'vc-relax',
    }

    if calculation in relaxation_calculations:
        if relaxation_settings is None:
            raise ValueError(
                f"QE {calculation} calculation requires "
                "relaxation settings."
            )

        resolved_calculation = str(
            relaxation_settings.get(
                'calculation',
                '',
            )
        ).strip().lower()

        if resolved_calculation != calculation:
            raise ValueError(
                "QE relaxation calculation type does not match "
                "the resolved relaxation settings."
            )

    elif relaxation_settings is not None:
        raise ValueError(
            "QE relaxation settings can only be used with "
            "calculation='relax' or calculation='vc-relax'."
        )

    if spinpol:
        if magnetic_moments is None:
            raise ValueError(
                "QE spin-polarized input requires "
                "initial magnetic moments."
            )

        magnetic_model = (
            build_qe_magnetic_species(
                atoms=atoms,
                pseudopotentials=pseudopotentials,
                pseudo_dir=pseudo_dir,
                magnetic_moments=magnetic_moments,
            )
        )

        species = magnetic_model[
            'species'
        ]

        positions = magnetic_model[
            'positions'
        ]

        starting_magnetizations = (
            magnetic_model[
                'starting_magnetizations'
            ]
        )

    else:
        if magnetic_moments is not None:
            raise ValueError(
                "QE initial magnetic moments require "
                "spinpol=True."
            )

        species = build_atomic_species(
            atoms,
            pseudopotentials,
        )

        positions = build_atomic_positions(
            atoms
        )

        starting_magnetizations = {}

    species_by_element = {}

    for element, position in zip(
        atoms.get_chemical_symbols(),
        positions['positions'],
    ):
        species_label = position[0]

        labels = species_by_element.setdefault(
            element,
            [],
        )

        if species_label not in labels:
            labels.append(
                species_label
            )

    hubbard_settings = resolve_qe_hubbard(
        setup_params=setup_params,
        pseudopotentials=pseudopotentials,
        pseudo_dir=pseudo_dir,
    )

    hubbard_card = render_qe_hubbard_card(
        hubbard_settings,
        species_by_element=species_by_element,
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
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
        exx_qpoint_grid=exx_qpoint_grid,
    )
    
    system.update(
        starting_magnetizations
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

    ions = None
    cell_dynamics = None

    if relaxation_settings is not None:
        control.update(
            relaxation_settings['control']
        )

        system.update(
            relaxation_settings['system']
        )

        ions = relaxation_settings['ions']
        cell_dynamics = relaxation_settings['cell']

        if (
            calculation == 'relax'
            and cell_dynamics is not None
        ):
            raise ValueError(
                "QE atomic relaxation cannot contain CELL settings."
            )

        if (
            calculation == 'vc-relax'
            and cell_dynamics is None
        ):
            raise ValueError(
                "QE variable-cell relaxation requires CELL settings."
            )

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
    ]

    if ions is not None:
        lines.append(
            render_namelist(
                'IONS',
                ions,
            )
        )

    if cell_dynamics is not None:
        lines.append(
            render_namelist(
                'CELL',
                cell_dynamics,
            )
        )

    if relaxation_settings is not None:
        for notice in relaxation_settings.get(
            'notices',
            [],
        ):
            lines.append(
                f"! NOTICE: {notice}"
            )

    lines.extend([
        '',
        'ATOMIC_SPECIES',
    ])

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

    if exx_additional_card is not None:
        lines.extend([
            '',
            exx_additional_card,
        ])

    lines.extend([
        '',
        f"CELL_PARAMETERS {cell['option']}",
    ])

    for x, y, z in cell['vectors']:
        lines.append(
            f"{x:.12f} {y:.12f} {z:.12f}"
        )

    if hubbard_card:
        lines.extend([
            '',
            hubbard_card,
        ])

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
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
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
    exx_additional_kpoints=None,
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
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
        exx_additional_kpoints=exx_additional_kpoints,
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
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
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
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
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

def render_relax_input(
    atoms,
    pseudopotentials,
    cutoff_ev,
    kpoint_size,
    optimizer,
    max_force,
    max_step,
    relax_cell,
    hydrostatic_pressure=0.0,
    fix_symmetry=False,
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
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
    """Render a complete QE pw.x relaxation input."""
    relaxation_settings = resolve_qe_relaxation_settings(
        optimizer=optimizer,
        max_force=max_force,
        max_step=max_step,
        relax_cell=relax_cell,
        hydrostatic_pressure=hydrostatic_pressure,
        fix_symmetry=fix_symmetry,
        atoms=atoms,
    )

    return render_pw_input(
        calculation=relaxation_settings['calculation'],
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=kpoint_size,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
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
        relaxation_settings=relaxation_settings,
    )

def render_bands_input(
    atoms,
    pseudopotentials,
    cutoff_ev,
    band_path,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
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
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
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

def render_bands_postprocess_input(
    prefix='nanoworks',
    outdir=None,
    filband='nanoworks.bands',
    lsym=False,
):
    """Render a complete Quantum ESPRESSO bands.x input."""
    if filband is None or not str(filband).strip():
        raise ValueError(
            "QE bands.x output file must not be empty."
        )

    settings = {
        'prefix': str(prefix),
        'filband': str(filband),
        'lsym': bool(lsym),
    }

    if outdir is not None:
        settings['outdir'] = str(outdir)

    return (
        render_namelist(
            'BANDS',
            settings,
        )
        + '\n'
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
    filproj=None,
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

    if filproj is not None:
        settings['filproj'] = str(
            filproj
        )

    return (
        render_namelist(
            'PROJWFC',
            settings,
        )
        + '\n'
    )

def render_pp_input(
    prefix='nanoworks',
    outdir=None,
    filplot='nanoworks.pp',
    fileout='nanoworks.cube',
    plot_num=0,
    spin_component=None,
):
    """Render a QE pp.x input for three-dimensional density data."""
    plot_num = int(
        plot_num
    )

    if plot_num not in {
        0,
        6,
    }:
        raise ValueError(
            "Unsupported QE density plot number: "
            f"{plot_num}"
        )

    input_settings = {
        'prefix': str(prefix),
        'filplot': str(filplot),
        'plot_num': plot_num,
    }

    if outdir is not None:
        input_settings[
            'outdir'
        ] = str(
            outdir
        )

    if spin_component is not None:
        if plot_num != 0:
            raise ValueError(
                "QE spin_component is supported only "
                "for total charge-density output."
            )

        spin_component = int(
            spin_component
        )

        if spin_component not in {
            0,
            1,
            2,
        }:
            raise ValueError(
                "QE charge-density spin_component "
                "must be 0, 1, or 2."
            )

        input_settings[
            'spin_component'
        ] = spin_component

    plot_settings = {
        'nfile': 1,
        'filepp(1)': str(
            filplot
        ),
        'weight(1)': 1.0,
        'iflag': 3,
        'output_format': 6,
        'fileout': str(
            fileout
        ),
    }

    return (
        render_namelist(
            'INPUTPP',
            input_settings,
        )
        + '\n'
        + render_namelist(
            'PLOT',
            plot_settings,
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


def parse_qe_auxiliary_output(
    output,
    expected_program=None,
):
    """Parse common completion metadata from a QE program output."""
    output = Path(
        output
    )

    if not output.is_file():
        raise FileNotFoundError(
            f"Quantum ESPRESSO output was not found: {output}"
        )

    text = output.read_text(
        encoding='utf-8',
        errors='replace',
    )

    program_match = re.search(
        r'Program\s+([A-Za-z0-9_.+-]+)\s+'
        r'v\.(\d+)\.(\d+)(?:\.(\d+))?',
        text,
        flags=re.IGNORECASE,
    )

    program = None
    qe_version = None

    if program_match:
        program = program_match.group(1).upper()
        qe_version = tuple(
            int(value)
            for value in program_match.groups()[1:]
            if value is not None
        )

    if expected_program is not None:
        expected_program = str(
            expected_program
        ).strip().upper()

        if not expected_program:
            raise ValueError(
                "Expected QE program name must not be empty."
            )

        if program is None:
            raise ValueError(
                "Quantum ESPRESSO program and version could not "
                f"be detected in '{output}'."
            )

        if program != expected_program:
            raise ValueError(
                "Unexpected Quantum ESPRESSO program in output: "
                f"expected {expected_program}, found {program}."
            )

    return {
        'program': program,
        'qe_version': qe_version,
        'job_done': 'JOB DONE.' in text,
    }


def parse_matdyn_frequency_file(frequency_file):
    """Parse q-points and phonon frequencies written by matdyn.x."""
    frequency_file = Path(
        frequency_file
    )

    if not frequency_file.is_file():
        raise FileNotFoundError(
            "QE matdyn frequency file was not found: "
            f"{frequency_file}"
        )

    text = frequency_file.read_text(
        encoding='utf-8',
        errors='replace',
    )

    header_match = re.search(
        r'&plot\s+nbnd\s*=\s*(\d+)\s*,\s*'
        r'nks\s*=\s*(\d+)\s*/',
        text,
        flags=re.IGNORECASE,
    )

    if header_match is None:
        raise ValueError(
            "QE matdyn frequency header could not be parsed "
            f"from '{frequency_file}'."
        )

    nmodes = int(
        header_match.group(1)
    )
    nqpoints = int(
        header_match.group(2)
    )

    if nmodes <= 0 or nqpoints <= 0:
        raise ValueError(
            "QE matdyn frequency header contains non-positive "
            "mode or q-point counts."
        )

    number_pattern = re.compile(
        r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)'
        r'(?:[EeDd][+-]?\d+)?'
    )

    data_text = text[
        header_match.end():
    ]

    values = [
        float(
            value.replace('D', 'E').replace('d', 'e')
        )
        for value in number_pattern.findall(data_text)
    ]

    values_per_qpoint = 3 + nmodes
    expected_values = (
        nqpoints * values_per_qpoint
    )

    if len(values) != expected_values:
        raise ValueError(
            "QE matdyn frequency data count does not match "
            f"the header: expected {expected_values} values, "
            f"found {len(values)}."
        )

    qpoints = []
    frequencies_cm1 = []

    for index in range(nqpoints):
        start = index * values_per_qpoint
        qpoints.append(
            tuple(values[start:start + 3])
        )
        frequencies_cm1.append(
            values[
                start + 3:
                start + values_per_qpoint
            ]
        )

    frequencies_thz = [
        [
            value * THZ_PER_CM_MINUS_ONE
            for value in row
        ]
        for row in frequencies_cm1
    ]

    return {
        'qpoints': qpoints,
        'frequencies_cm1': frequencies_cm1,
        'frequencies_thz': frequencies_thz,
        'nqpoints': nqpoints,
        'nmodes': nmodes,
    }


def parse_matdyn_dos_file(dos_file):
    """Parse total and atom-projected phonon DOS from matdyn.x."""
    dos_file = Path(
        dos_file
    )

    if not dos_file.is_file():
        raise FileNotFoundError(
            "QE matdyn DOS file was not found: "
            f"{dos_file}"
        )

    number_pattern = re.compile(
        r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)'
        r'(?:[EeDd][+-]?\d+)?'
    )

    rows = []
    column_count = None

    with dos_file.open(
        'r',
        encoding='utf-8',
        errors='replace',
    ) as fd:
        for line_number, line in enumerate(fd, start=1):
            stripped = line.strip()

            if not stripped or stripped.startswith('#'):
                continue

            matches = number_pattern.findall(stripped)
            remainder = number_pattern.sub('', stripped)

            if not matches or remainder.strip():
                raise ValueError(
                    "QE matdyn DOS row could not be parsed "
                    f"at line {line_number} in '{dos_file}'."
                )

            values = [
                float(
                    value.replace('D', 'E').replace('d', 'e')
                )
                for value in matches
            ]

            if len(values) < 3:
                raise ValueError(
                    "QE matdyn DOS rows must contain frequency, "
                    "total DOS, and at least one projected DOS "
                    f"column in '{dos_file}'."
                )

            if column_count is None:
                column_count = len(values)
            elif len(values) != column_count:
                raise ValueError(
                    "QE matdyn DOS data contains inconsistent "
                    f"column counts in '{dos_file}'."
                )

            rows.append(values)

    if not rows:
        raise ValueError(
            "No QE matdyn DOS data could be parsed from "
            f"'{dos_file}'."
        )

    frequencies_cm1 = [
        row[0]
        for row in rows
    ]
    frequencies_thz = [
        value * THZ_PER_CM_MINUS_ONE
        for value in frequencies_cm1
    ]
    dos = [
        row[1]
        for row in rows
    ]
    natoms = column_count - 2
    atom_projected_dos = [
        [
            row[atom_index + 2]
            for row in rows
        ]
        for atom_index in range(natoms)
    ]

    return {
        'frequencies_cm1': frequencies_cm1,
        'frequencies_thz': frequencies_thz,
        'dos': dos,
        'atom_projected_dos': atom_projected_dos,
        'npoints': len(rows),
        'natoms': natoms,
    }


def write_matdyn_band_data(output_file, band_path, frequencies):
    """Write parsed QE phonon bands as a plot-ready THz table."""
    output_file = Path(
        output_file
    )

    if not isinstance(band_path, dict):
        raise TypeError(
            "QE phonon band path must be a mapping."
        )

    distances = list(
        band_path.get('distances', [])
    )
    qpoints = list(
        frequencies.get('qpoints', [])
    )
    frequency_rows = list(
        frequencies.get('frequencies_thz', [])
    )
    nqpoints = frequencies.get('nqpoints')
    nmodes = frequencies.get('nmodes')

    if not distances or not qpoints or not frequency_rows:
        raise ValueError(
            "QE phonon band data is incomplete."
        )

    row_count = len(distances)

    if not (
        len(qpoints) == row_count
        and len(frequency_rows) == row_count
        and nqpoints == row_count
    ):
        raise ValueError(
            "QE phonon band distances, q-points, and frequency "
            "counts do not match."
        )

    if not isinstance(nmodes, int) or nmodes <= 0:
        raise ValueError(
            "QE phonon band mode count must be a positive integer."
        )

    if any(len(qpoint) != 3 for qpoint in qpoints):
        raise ValueError(
            "Each QE phonon band q-point must contain 3 values."
        )

    if any(len(row) != nmodes for row in frequency_rows):
        raise ValueError(
            "QE phonon band rows do not match the mode count."
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_file.open(
        'w',
        encoding='utf-8',
    ) as fd:
        mode_columns = ' '.join(
            f"Frequency_{index + 1}(THz)"
            for index in range(nmodes)
        )
        print(
            "# Distance(1/Angstrom) qx qy qz "
            + mode_columns,
            file=fd,
        )

        for distance, qpoint, row in zip(
            distances,
            qpoints,
            frequency_rows,
        ):
            values = (
                float(distance),
                *(float(value) for value in qpoint),
                *(float(value) for value in row),
            )
            print(
                ' '.join(
                    f"{value:.10f}"
                    for value in values
                ),
                file=fd,
            )

    return output_file


def write_matdyn_dos_data(output_file, dos_data):
    """Write parsed QE phonon DOS with THz frequencies and units."""
    output_file = Path(
        output_file
    )

    frequencies = list(
        dos_data.get('frequencies_thz', [])
    )
    total_dos_cm1 = list(
        dos_data.get('dos', [])
    )
    projected_dos_cm1 = list(
        dos_data.get('atom_projected_dos', [])
    )
    npoints = dos_data.get('npoints')
    natoms = dos_data.get('natoms')

    if not frequencies or not total_dos_cm1:
        raise ValueError(
            "QE phonon DOS data is incomplete."
        )

    if not (
        len(frequencies) == len(total_dos_cm1)
        and npoints == len(frequencies)
    ):
        raise ValueError(
            "QE phonon DOS frequency and value counts do not match."
        )

    if not isinstance(natoms, int) or natoms <= 0:
        raise ValueError(
            "QE phonon DOS atom count must be a positive integer."
        )

    if (
        len(projected_dos_cm1) != natoms
        or any(
            len(values) != npoints
            for values in projected_dos_cm1
        )
    ):
        raise ValueError(
            "QE atom-projected phonon DOS dimensions do not match."
        )

    total_dos_thz = [
        value / THZ_PER_CM_MINUS_ONE
        for value in total_dos_cm1
    ]
    projected_dos_thz = [
        [
            value / THZ_PER_CM_MINUS_ONE
            for value in atom_values
        ]
        for atom_values in projected_dos_cm1
    ]

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_file.open(
        'w',
        encoding='utf-8',
    ) as fd:
        projected_columns = ' '.join(
            f"Atom_{index + 1}_PDOS(1/THz)"
            for index in range(natoms)
        )
        print(
            "# Frequency(THz) DOS(1/THz) "
            + projected_columns,
            file=fd,
        )

        for point_index, frequency in enumerate(frequencies):
            values = (
                float(frequency),
                float(total_dos_thz[point_index]),
                *(
                    float(atom_values[point_index])
                    for atom_values in projected_dos_thz
                ),
            )
            print(
                ' '.join(
                    f"{value:.10f}"
                    for value in values
                ),
                file=fd,
            )

    return output_file


def calculate_phonon_thermal_properties(
    dos_data,
    t_min=0.0,
    t_max=1000.0,
    t_step=10.0,
    cutoff_frequency_thz=0.0,
):
    """Integrate harmonic thermal properties over a QE phonon DOS."""
    try:
        t_min = float(t_min)
        t_max = float(t_max)
        t_step = float(t_step)
        cutoff_frequency_thz = float(
            cutoff_frequency_thz
        )
    except (TypeError, ValueError) as exc:
        raise TypeError(
            "QE phonon thermal settings must be real numbers."
        ) from exc

    if not all(
        math.isfinite(value)
        for value in (
            t_min,
            t_max,
            t_step,
            cutoff_frequency_thz,
        )
    ):
        raise ValueError(
            "QE phonon thermal settings must be finite."
        )

    if t_min < 0.0 or t_max < t_min or t_step <= 0.0:
        raise ValueError(
            "QE phonon temperatures require 0 <= t_min <= "
            "t_max and t_step > 0."
        )

    if cutoff_frequency_thz < 0.0:
        raise ValueError(
            "QE phonon thermal cutoff frequency must not be negative."
        )

    frequencies = np.asarray(
        dos_data.get('frequencies_thz', []),
        dtype=float,
    )
    dos_cm1 = np.asarray(
        dos_data.get('dos', []),
        dtype=float,
    )

    if (
        frequencies.ndim != 1
        or dos_cm1.ndim != 1
        or frequencies.size < 2
        or frequencies.size != dos_cm1.size
    ):
        raise ValueError(
            "QE phonon DOS data dimensions are inconsistent."
        )

    if not (
        np.all(np.isfinite(frequencies))
        and np.all(np.isfinite(dos_cm1))
    ):
        raise ValueError(
            "QE phonon DOS data must contain only finite values."
        )

    if np.any(np.diff(frequencies) <= 0.0):
        raise ValueError(
            "QE phonon DOS frequencies must be strictly increasing."
        )

    negative_tolerance = max(
        1.0,
        float(np.max(np.abs(dos_cm1))),
    ) * 1.0e-12

    if np.any(dos_cm1 < -negative_tolerance):
        raise ValueError(
            "QE phonon DOS values must not be negative."
        )

    dos_thz = np.maximum(
        dos_cm1,
        0.0,
    ) / THZ_PER_CM_MINUS_ONE
    positive_mask = (
        frequencies > cutoff_frequency_thz
    )

    if np.count_nonzero(positive_mask) < 2:
        raise ValueError(
            "QE phonon DOS does not contain enough positive "
            "frequency points for thermal integration."
        )

    positive_frequencies = frequencies[
        positive_mask
    ]
    positive_dos = dos_thz[
        positive_mask
    ]
    mode_energies = (
        positive_frequencies
        * EV_PER_THZ
    )

    total_mode_weight = float(
        np.trapz(
            dos_thz,
            frequencies,
        )
    )
    integrated_mode_weight = float(
        np.trapz(
            positive_dos,
            positive_frequencies,
        )
    )
    excluded_mode_weight = max(
        0.0,
        total_mode_weight - integrated_mode_weight,
    )
    zero_point_energy_ev = float(
        np.trapz(
            0.5 * mode_energies * positive_dos,
            positive_frequencies,
        )
    )

    temperatures = np.arange(
        t_min,
        t_max + t_step / 2.0,
        t_step,
        dtype=float,
    )
    free_energy_ev = np.empty_like(
        temperatures
    )
    internal_energy_ev = np.empty_like(
        temperatures
    )
    entropy_ev_per_k = np.empty_like(
        temperatures
    )
    heat_capacity_ev_per_k = np.empty_like(
        temperatures
    )

    for index, temperature in enumerate(temperatures):
        if temperature <= 0.0:
            free_energy_ev[index] = zero_point_energy_ev
            internal_energy_ev[index] = zero_point_energy_ev
            entropy_ev_per_k[index] = 0.0
            heat_capacity_ev_per_k[index] = 0.0
            continue

        x = mode_energies / (
            BOLTZMANN_EV_PER_K * temperature
        )
        exp_negative_x = np.exp(-x)
        one_minus_exp_negative_x = -np.expm1(-x)
        occupation = (
            exp_negative_x
            / one_minus_exp_negative_x
        )
        log_bose_factor = np.log(
            one_minus_exp_negative_x
        )

        free_energy_integrand = positive_dos * (
            0.5 * mode_energies
            + BOLTZMANN_EV_PER_K
            * temperature
            * log_bose_factor
        )
        internal_energy_integrand = positive_dos * (
            0.5 * mode_energies
            + mode_energies * occupation
        )
        entropy_integrand = (
            positive_dos
            * BOLTZMANN_EV_PER_K
            * (
                x * occupation
                - log_bose_factor
            )
        )
        heat_capacity_integrand = (
            positive_dos
            * BOLTZMANN_EV_PER_K
            * x ** 2
            * exp_negative_x
            / one_minus_exp_negative_x ** 2
        )

        free_energy_ev[index] = np.trapz(
            free_energy_integrand,
            positive_frequencies,
        )
        internal_energy_ev[index] = np.trapz(
            internal_energy_integrand,
            positive_frequencies,
        )
        entropy_ev_per_k[index] = np.trapz(
            entropy_integrand,
            positive_frequencies,
        )
        heat_capacity_ev_per_k[index] = np.trapz(
            heat_capacity_integrand,
            positive_frequencies,
        )

    energy_conversion = KJ_PER_MOL_PER_EV
    entropy_conversion = (
        KJ_PER_MOL_PER_EV * 1000.0
    )

    return {
        'temperatures_k': temperatures.tolist(),
        'free_energy_kj_mol': (
            free_energy_ev * energy_conversion
        ).tolist(),
        'internal_energy_kj_mol': (
            internal_energy_ev * energy_conversion
        ).tolist(),
        'entropy_j_k_mol': (
            entropy_ev_per_k * entropy_conversion
        ).tolist(),
        'heat_capacity_j_k_mol': (
            heat_capacity_ev_per_k * entropy_conversion
        ).tolist(),
        'zero_point_energy_kj_mol': (
            zero_point_energy_ev * energy_conversion
        ),
        'integrated_mode_weight': integrated_mode_weight,
        'excluded_mode_weight': excluded_mode_weight,
        'cutoff_frequency_thz': cutoff_frequency_thz,
    }


def write_phonon_thermal_properties(output_file, thermal_data):
    """Write QE harmonic phonon thermal properties to CSV."""
    output_file = Path(
        output_file
    )
    columns = (
        thermal_data.get('temperatures_k', []),
        thermal_data.get('free_energy_kj_mol', []),
        thermal_data.get('internal_energy_kj_mol', []),
        thermal_data.get('entropy_j_k_mol', []),
        thermal_data.get('heat_capacity_j_k_mol', []),
    )
    lengths = {
        len(column)
        for column in columns
    }

    if lengths == {0} or len(lengths) != 1:
        raise ValueError(
            "QE phonon thermal property columns do not match."
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with output_file.open(
        'w',
        encoding='utf-8',
    ) as fd:
        print(
            "T(K),Free_Energy(kJ/mol),Internal_Energy(kJ/mol),"
            "Entropy(J/K/mol),Cv(J/K/mol)",
            file=fd,
        )

        for row in zip(*columns):
            print(
                ','.join(
                    f"{float(value):.10f}"
                    for value in row
                ),
                file=fd,
            )

    return output_file


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

def parse_pw_relaxed_structure(
    output,
    reference_atoms,
):
    """Read the final relaxed structure from QE pw.x output."""
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

    coordinate_blocks = re.findall(
        r'Begin\s+final\s+coordinates'
        r'(.*?)'
        r'End\s+final\s+coordinates',
        text,
        flags=(
            re.IGNORECASE
            | re.DOTALL
        ),
    )

    if not coordinate_blocks:
        raise ValueError(
            "No final relaxed coordinates were found "
            f"in the QE output file: {output}"
        )

    block = coordinate_blocks[-1]

    alat_matches = re.findall(
        r'lattice\s+parameter\s+\(alat\)\s*=\s*'
        r'([-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)'
        r'\s*a\.u\.',
        text,
        flags=re.IGNORECASE,
    )

    alat_angstrom = None

    if alat_matches:
        alat_angstrom = (
            float(
                alat_matches[-1]
                .replace('D', 'E')
                .replace('d', 'e')
            )
            * Bohr
        )

    def find_card(name):
        match = re.search(
            rf'{name}\s*'
            r'(?:\(\s*([^)]+?)\s*\)|([A-Za-z_]+))?'
            r'[^\n]*\n',
            block,
            flags=re.IGNORECASE,
        )

        if match is None:
            return None

        option = (
            match.group(1)
            or match.group(2)
            or 'alat'
        )

        return (
            option.strip().lower(),
            block[match.end():],
        )

    def read_rows(card_text, count):
        rows = []

        for line in card_text.splitlines():
            fields = line.split()

            if not fields:
                continue

            if len(fields) < 4:
                break

            try:
                values = [
                    float(
                        value
                        .replace('D', 'E')
                        .replace('d', 'e')
                    )
                    for value in fields[1:4]
                ]
            except ValueError:
                break

            rows.append(
                (
                    fields[0],
                    values,
                )
            )

            if len(rows) == count:
                break

        if len(rows) != count:
            raise ValueError(
                "QE final coordinate block does not contain "
                f"{count} atomic positions."
            )

        return rows

    atoms = reference_atoms.copy()
    natoms = len(atoms)

    cell_card = find_card(
        'CELL_PARAMETERS'
    )

    if cell_card is not None:
        cell_option, cell_text = cell_card
        cell_rows = []

        for line in cell_text.splitlines():
            fields = line.split()

            if not fields:
                continue

            if len(fields) < 3:
                break

            try:
                row = [
                    float(
                        value
                        .replace('D', 'E')
                        .replace('d', 'e')
                    )
                    for value in fields[:3]
                ]
            except ValueError:
                break

            cell_rows.append(
                row
            )

            if len(cell_rows) == 3:
                break

        if len(cell_rows) != 3:
            raise ValueError(
                "QE final coordinate block does not contain "
                "three cell vectors."
            )

        if cell_option == 'angstrom':
            cell_scale = 1.0

        elif cell_option == 'bohr':
            cell_scale = Bohr

        elif cell_option == 'alat':
            if alat_angstrom is None:
                raise ValueError(
                    "QE alat could not be determined for "
                    "the final cell parameters."
                )

            cell_scale = alat_angstrom

        else:
            raise ValueError(
                "Unsupported QE final cell unit: "
                f"{cell_option}"
            )

        atoms.set_cell(
            [
                [
                    value * cell_scale
                    for value in row
                ]
                for row in cell_rows
            ],
            scale_atoms=False,
        )

    position_card = find_card(
        'ATOMIC_POSITIONS'
    )

    if position_card is None:
        raise ValueError(
            "No ATOMIC_POSITIONS card was found in "
            "the QE final coordinate block."
        )

    position_option, position_text = position_card

    position_rows = read_rows(
        position_text,
        natoms,
    )

    symbols = [
        symbol
        for symbol, values in position_rows
    ]

    if symbols != atoms.get_chemical_symbols():
        raise ValueError(
            "QE final atom ordering does not match "
            "the input structure."
        )

    coordinates = [
        values
        for symbol, values in position_rows
    ]

    if position_option == 'crystal':
        atoms.set_scaled_positions(
            coordinates
        )

    elif position_option == 'angstrom':
        atoms.set_positions(
            coordinates
        )

    elif position_option == 'bohr':
        atoms.set_positions([
            [
                value * Bohr
                for value in position
            ]
            for position in coordinates
        ])

    elif position_option == 'alat':
        if alat_angstrom is None:
            raise ValueError(
                "QE alat could not be determined for "
                "the final atomic positions."
            )

        atoms.set_positions([
            [
                value * alat_angstrom
                for value in position
            ]
            for position in coordinates
        ])

    else:
        raise ValueError(
            "Unsupported QE final position unit: "
            f"{position_option}"
        )

    return atoms

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
    """Parse non-spin or collinear-spin eigenvalues from QE pw.x output."""
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

    number_regex = re.compile(
        number_pattern
    )

    def parse_eigenvalue_line(line):
        return [
            float(
                value
                .replace('D', 'E')
                .replace('d', 'e')
            )
            for value in number_regex.findall(
                line
            )
        ]

    def parse_channel(channel_text):
        kpoints = []
        eigenvalues = []

        current_kpoint = None
        current_eigenvalues = None

        for line in channel_text.splitlines():
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
                        parse_eigenvalue_line(
                            remainder
                        )
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
                    parse_eigenvalue_line(
                        stripped
                    )
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
            'kpoints': kpoints,
            'eigenvalues_ev': eigenvalues,
            'nkpoints': len(kpoints),
            'nbands': len(eigenvalues[0]),
        }

    spin_pattern = re.compile(
        r'^\s*[-=]*\s*SPIN\s+(UP|DOWN)\s*[-=]*\s*$',
        flags=(
            re.IGNORECASE
            | re.MULTILINE
        ),
    )

    spin_markers = list(
        spin_pattern.finditer(
            text
        )
    )

    if not spin_markers:
        channel = parse_channel(
            text
        )

        return {
            'spin_polarized': False,
            'nspins': 1,
            'nkpoints': channel['nkpoints'],
            'nbands': channel['nbands'],
            'kpoints': channel['kpoints'],
            'eigenvalues_ev': [
                channel['eigenvalues_ev']
            ],
        }

    channels = {}

    for index, marker in enumerate(
        spin_markers
    ):
        spin_name = (
            marker.group(1)
            .strip()
            .lower()
        )

        if spin_name in channels:
            raise ValueError(
                "QE band output contains duplicate "
                f"SPIN {spin_name.upper()} sections."
            )

        section_end = (
            spin_markers[index + 1].start()
            if index + 1 < len(spin_markers)
            else len(text)
        )

        channels[spin_name] = parse_channel(
            text[
                marker.end():
                section_end
            ]
        )

    if set(channels) != {
        'up',
        'down',
    }:
        raise ValueError(
            "Spin-polarized QE band output must contain "
            "both SPIN UP and SPIN DOWN sections."
        )

    up = channels['up']
    down = channels['down']

    if up['nkpoints'] != down['nkpoints']:
        raise ValueError(
            "QE spin channels contain different "
            "numbers of k-points."
        )

    if up['nbands'] != down['nbands']:
        raise ValueError(
            "QE spin channels contain different "
            "numbers of bands."
        )

    for up_kpoint, down_kpoint in zip(
        up['kpoints'],
        down['kpoints'],
    ):
        if any(
            abs(up_value - down_value) > 1.0e-10
            for up_value, down_value in zip(
                up_kpoint,
                down_kpoint,
            )
        ):
            raise ValueError(
                "QE spin channels use different "
                "k-point grids."
            )

    return {
        'spin_polarized': True,
        'nspins': 2,
        'nkpoints': up['nkpoints'],
        'nbands': up['nbands'],
        'kpoints': up['kpoints'],
        'eigenvalues_ev': [
            up['eigenvalues_ev'],
            down['eigenvalues_ev'],
        ],
    }

def parse_bands_x_output(
    band_file,
    kpoint_indices=None,
):
    """Parse a non-spin bands.x filband file."""
    band_file = Path(
        band_file
    )

    if not band_file.is_file():
        raise FileNotFoundError(
            f"QE bands.x data file was not found: {band_file}"
        )

    text = band_file.read_text(
        encoding='utf-8',
        errors='replace',
    )

    number_pattern = (
        r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)'
        r'(?:[EeDd][-+]?\d+)?'
    )
    header_pattern = re.compile(
        rf'nbnd\s*=\s*({number_pattern})\s*,?\s*'
        rf'nks\s*=\s*({number_pattern})',
        flags=re.IGNORECASE,
    )

    header_match = header_pattern.search(
        text
    )

    if header_match is None:
        raise ValueError(
            "QE bands.x data file does not contain an "
            "nbnd/nks header."
        )

    try:
        nbands = int(
            float(
                header_match.group(1)
            )
        )
        nkpoints = int(
            float(
                header_match.group(2)
            )
        )
    except (TypeError, ValueError) as error:
        raise ValueError(
            "QE bands.x data file has an invalid nbnd/nks header."
        ) from error

    if nbands <= 0 or nkpoints <= 0:
        raise ValueError(
            "QE bands.x data file must contain positive nbnd and nks."
        )

    number_regex = re.compile(
        number_pattern
    )
    body = text[
        header_match.end():
    ]
    values = [
        float(
            value.replace('D', 'E').replace('d', 'e')
        )
        for value in number_regex.findall(body)
    ]
    values_per_kpoint = 3 + nbands
    expected_values = nkpoints * values_per_kpoint

    if len(values) != expected_values:
        raise ValueError(
            "QE bands.x data file contains "
            f"{len(values)} numeric values, expected "
            f"{expected_values}."
        )

    all_kpoints = []
    all_eigenvalues = []
    cursor = 0

    for _ in range(nkpoints):
        all_kpoints.append(
            tuple(
                values[cursor:cursor + 3]
            )
        )
        cursor += 3
        all_eigenvalues.append(
            values[cursor:cursor + nbands]
        )
        cursor += nbands

    if kpoint_indices is None:
        selected_indices = list(
            range(nkpoints)
        )
    else:
        selected_indices = [
            int(index)
            for index in kpoint_indices
        ]

        if not selected_indices:
            raise ValueError(
                "QE bands.x k-point selection must not be empty."
            )

        if any(
            index < 0 or index >= nkpoints
            for index in selected_indices
        ):
            raise ValueError(
                "QE bands.x k-point selection is out of range."
            )

    return {
        'spin_polarized': False,
        'nspins': 1,
        'nkpoints': len(selected_indices),
        'nbands': nbands,
        'kpoints': [
            all_kpoints[index]
            for index in selected_indices
        ],
        'eigenvalues_ev': [[
            all_eigenvalues[index]
            for index in selected_indices
        ]],
    }

def parse_projwfc_band_file(
    projection_file,
):
    """Parse one QE projwfc.x band-projection file."""
    projection_file = Path(
        projection_file
    )

    if not projection_file.is_file():
        raise FileNotFoundError(
            "QE band projection file was not found: "
            f"{projection_file}"
        )

    lines = projection_file.read_text(
        encoding='utf-8',
        errors='replace',
    ).splitlines()

    line_index = 0

    def next_fields(description):
        nonlocal line_index

        while (
            line_index < len(lines)
            and not lines[line_index].strip()
        ):
            line_index += 1

        if line_index >= len(lines):
            raise ValueError(
                "Unexpected end of QE band projection file "
                f"while reading {description}: "
                f"{projection_file}"
            )

        fields = lines[
            line_index
        ].split()

        line_index += 1

        return fields

    def parse_integer(value):
        return int(
            float(
                value
                .replace('D', 'E')
                .replace('d', 'e')
            )
        )

    def parse_float(value):
        return float(
            value
            .replace('D', 'E')
            .replace('d', 'e')
        )

    def parse_logical(value):
        return (
            value
            .strip()
            .strip('.')
            .upper()
            .startswith('T')
        )

    try:
        grid_fields = next_fields(
            'FFT-grid dimensions'
        )

        if len(grid_fields) < 8:
            raise ValueError(
                "The FFT-grid header must contain "
                "eight integer values."
            )

        natoms = parse_integer(
            grid_fields[6]
        )

        ntypes = parse_integer(
            grid_fields[7]
        )

        lattice_fields = next_fields(
            'lattice information'
        )

        if len(lattice_fields) < 7:
            raise ValueError(
                "The lattice header must contain "
                "ibrav and six celldm values."
            )

        ibrav = parse_integer(
            lattice_fields[0]
        )

        if ibrav == 0:
            for vector_index in range(3):
                vector_fields = next_fields(
                    'lattice vector'
                )

                if len(vector_fields) < 3:
                    raise ValueError(
                        "A QE lattice vector must contain "
                        "three values."
                    )

                for value in vector_fields[:3]:
                    parse_float(
                        value
                    )

        cutoff_fields = next_fields(
            'cutoff information'
        )

        if len(cutoff_fields) < 4:
            raise ValueError(
                "The QE cutoff header is incomplete."
            )

        species = {}

        for species_index in range(
            ntypes
        ):
            species_fields = next_fields(
                'atomic species'
            )

            if len(species_fields) < 2:
                raise ValueError(
                    "A QE atomic-species row is incomplete."
                )

            qe_species_index = parse_integer(
                species_fields[0]
            )

            species[
                qe_species_index
            ] = species_fields[1]

        atom_species = {}

        for atom_index in range(
            natoms
        ):
            atom_fields = next_fields(
                'atomic position'
            )

            if len(atom_fields) < 5:
                raise ValueError(
                    "A QE atomic-position row is incomplete."
                )

            qe_atom_index = parse_integer(
                atom_fields[0]
            )

            qe_species_index = parse_integer(
                atom_fields[4]
            )

            if qe_species_index not in species:
                raise ValueError(
                    "QE band projection file references "
                    "an unknown atomic species."
                )

            atom_species[
                qe_atom_index
            ] = species[
                qe_species_index
            ]

        dimension_fields = next_fields(
            'projection dimensions'
        )

        if len(dimension_fields) < 3:
            raise ValueError(
                "The QE projection-dimension header "
                "is incomplete."
            )

        natomic_states = parse_integer(
            dimension_fields[0]
        )

        nkpoints = parse_integer(
            dimension_fields[1]
        )

        nbands = parse_integer(
            dimension_fields[2]
        )

        if (
            natomic_states <= 0
            or nkpoints <= 0
            or nbands <= 0
        ):
            raise ValueError(
                "QE band projection dimensions "
                "must be positive."
            )

        spin_fields = next_fields(
            'spin flags'
        )

        if len(spin_fields) < 2:
            raise ValueError(
                "The QE projection spin header "
                "is incomplete."
            )

        noncollinear = parse_logical(
            spin_fields[0]
        )

        spin_orbit = parse_logical(
            spin_fields[1]
        )

    except (TypeError, ValueError) as exc:
        if isinstance(exc, ValueError) and str(
            exc
        ).startswith(
            (
                'Unexpected end',
                'The ',
                'A QE ',
                'QE band ',
            )
        ):
            raise

        raise ValueError(
            "Invalid QE band projection header in "
            f"{projection_file}: {exc}"
        ) from exc

    if noncollinear or spin_orbit:
        raise NotImplementedError(
            "Noncollinear and spin-orbit QE band "
            "projections are not supported."
        )

    orbital_names = {
        0: 's',
        1: 'p',
        2: 'd',
        3: 'f',
    }

    states = []

    for expected_state in range(
        1,
        natomic_states + 1,
    ):
        state_fields = next_fields(
            'atomic-state header'
        )

        if len(state_fields) < 7:
            raise ValueError(
                "A QE atomic-state header must contain "
                "seven fields."
            )

        try:
            state_index = parse_integer(
                state_fields[0]
            )

            qe_atom_index = parse_integer(
                state_fields[1]
            )

            wfc_index = parse_integer(
                state_fields[4]
            )

            angular_momentum = parse_integer(
                state_fields[5]
            )

            magnetic_component = parse_integer(
                state_fields[6]
            )

        except ValueError as exc:
            raise ValueError(
                "Invalid QE atomic-state header in "
                f"{projection_file}: {exc}"
            ) from exc

        if state_index != expected_state:
            raise ValueError(
                "Unexpected QE atomic-state index "
                f"{state_index}; expected {expected_state}."
            )

        if qe_atom_index not in atom_species:
            raise ValueError(
                "QE band projection state references "
                f"unknown atom index {qe_atom_index}."
            )

        weights = [
            [
                None
                for band_index in range(
                    nbands
                )
            ]
            for kpoint_index in range(
                nkpoints
            )
        ]

        for value_index in range(
            nkpoints * nbands
        ):
            value_fields = next_fields(
                'projection weight'
            )

            if len(value_fields) < 3:
                raise ValueError(
                    "A QE projection-weight row must "
                    "contain three values."
                )

            try:
                raw_kpoint_index = parse_integer(
                    value_fields[0]
                )

                band_index = (
                    parse_integer(
                        value_fields[1]
                    )
                    - 1
                )

                weight = parse_float(
                    value_fields[2]
                )

            except ValueError as exc:
                raise ValueError(
                    "Invalid QE projection weight in "
                    f"{projection_file}: {exc}"
                ) from exc

            if (
                1
                <= raw_kpoint_index
                <= nkpoints
            ):
                kpoint_index = (
                    raw_kpoint_index
                    - 1
                )

            elif (
                nkpoints + 1
                <= raw_kpoint_index
                <= 2 * nkpoints
            ):
                kpoint_index = (
                    raw_kpoint_index
                    - nkpoints
                    - 1
                )

            else:
                raise ValueError(
                    "QE projection k-point index is "
                    f"out of range: {raw_kpoint_index}"
                )

            if not (
                0
                <= band_index
                < nbands
            ):
                raise ValueError(
                    "QE projection band index is "
                    f"out of range: {band_index + 1}"
                )

            if weights[
                kpoint_index
            ][band_index] is not None:
                raise ValueError(
                    "Duplicate QE projection weight for "
                    f"k-point {raw_kpoint_index}, "
                    f"band {band_index + 1}."
                )

            weights[
                kpoint_index
            ][band_index] = weight

        if any(
            weight is None
            for kpoint_weights in weights
            for weight in kpoint_weights
        ):
            raise ValueError(
                "QE band projection file contains "
                "incomplete projection weights."
            )

        states.append({
            'state_index': state_index,
            'atom_index': (
                qe_atom_index
                - 1
            ),
            'qe_atom_index': qe_atom_index,
            'symbol': state_fields[2],
            'species': atom_species[
                qe_atom_index
            ],
            'orbital_label': state_fields[3],
            'wfc_index': wfc_index,
            'l': angular_momentum,
            'm': magnetic_component,
            'orbital': orbital_names.get(
                angular_momentum,
                f'l={angular_momentum}',
            ),
            'weights': weights,
        })

    return {
        'file': projection_file,
        'natoms': natoms,
        'ntypes': ntypes,
        'natomic_states': natomic_states,
        'nkpoints': nkpoints,
        'nbands': nbands,
        'noncollinear': noncollinear,
        'spin_orbit': spin_orbit,
        'states': states,
    }

def prepare_qe_band_projection_data(
    projection_result,
    projections=None,
):
    """Aggregate QE atomic-state weights for requested projections."""
    natoms = projection_result[
        'natoms'
    ]

    nkpoints = projection_result[
        'nkpoints'
    ]

    nbands = projection_result[
        'nbands'
    ]

    states = projection_result[
        'states'
    ]

    if not projections:
        projections = [{
            'atoms': list(
                range(natoms)
            ),
            'orbital': None,
            'color': 'blue',
            'label': 'Total Contribution',
        }]

    allowed_orbitals = {
        's',
        'p',
        'd',
        'f',
    }

    prepared = []

    for projection_index, projection in enumerate(
        projections
    ):
        if not isinstance(
            projection,
            dict,
        ):
            raise TypeError(
                "Each QE band projection must be "
                "a dictionary."
            )

        atom_indices = projection.get(
            'atoms',
            [],
        )

        if not isinstance(
            atom_indices,
            (
                list,
                tuple,
                set,
                range,
            ),
        ):
            raise TypeError(
                "QE band projection 'atoms' must "
                "be a sequence of zero-based indices."
            )

        atom_indices = list(
            atom_indices
        )

        for atom_index in atom_indices:
            if (
                isinstance(
                    atom_index,
                    bool,
                )
                or not isinstance(
                    atom_index,
                    int,
                )
            ):
                raise TypeError(
                    "QE band projection atom indices "
                    "must be integers."
                )

            if not (
                0
                <= atom_index
                < natoms
            ):
                raise ValueError(
                    "QE band projection atom index "
                    f"{atom_index} is out of range for "
                    f"a structure containing {natoms} atoms."
                )

        orbital = projection.get(
            'orbital'
        )

        if orbital is not None:
            if not isinstance(
                orbital,
                str,
            ):
                raise TypeError(
                    "QE band projection 'orbital' must "
                    "be s, p, d, f, or None."
                )

            orbital = orbital.lower()

            if orbital not in allowed_orbitals:
                raise ValueError(
                    "Unsupported QE band projection "
                    f"orbital: {orbital}"
                )

        color = projection.get(
            'color',
            'blue',
        )

        label = projection.get(
            'label',
            f"Atoms: {atom_indices}",
        )

        weights = [
            [
                0.0
                for band_index in range(
                    nbands
                )
            ]
            for kpoint_index in range(
                nkpoints
            )
        ]

        selected_state_count = 0
        atom_index_set = set(
            atom_indices
        )

        for state in states:
            state_weights = state[
                'weights'
            ]

            if len(
                state_weights
            ) != nkpoints:
                raise ValueError(
                    "QE band projection state has "
                    "an inconsistent number of k-points."
                )

            if state[
                'atom_index'
            ] not in atom_index_set:
                continue

            if (
                orbital is not None
                and state['orbital'] != orbital
            ):
                continue

            selected_state_count += 1

            for kpoint_index in range(
                nkpoints
            ):
                if len(
                    state_weights[
                        kpoint_index
                    ]
                ) != nbands:
                    raise ValueError(
                        "QE band projection state has "
                        "an inconsistent number of bands."
                    )

                for band_index in range(
                    nbands
                ):
                    weights[
                        kpoint_index
                    ][
                        band_index
                    ] += state_weights[
                        kpoint_index
                    ][
                        band_index
                    ]

        prepared.append({
            'index': projection_index,
            'atoms': atom_indices,
            'orbital': orbital,
            'color': color,
            'label': label,
            'selected_state_count': (
                selected_state_count
            ),
            'weights': weights,
        })

    return {
        'natoms': natoms,
        'nkpoints': nkpoints,
        'nbands': nbands,
        'projections': prepared,
    }

def prepare_qe_band_data(
    bands,
    band_path,
    reference_energy,
):
    """Prepare non-spin or collinear-spin QE band data."""
    spin_polarized = bool(
        bands.get(
            'spin_polarized',
            False,
        )
    )

    nspins = int(
        bands.get(
            'nspins',
            0,
        )
    )

    expected_nspins = (
        2
        if spin_polarized
        else 1
    )

    if nspins != expected_nspins:
        raise ValueError(
            "QE band spin metadata is inconsistent."
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

    if len(eigenvalues_by_spin) != nspins:
        raise ValueError(
            "QE band data does not contain the reported "
            "number of spin channels."
        )

    for eigenvalues in eigenvalues_by_spin:
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

    shifted_by_spin = [
        [
            [
                float(value)
                - reference_energy
                for value in values
            ]
            for values in eigenvalues
        ]
        for eigenvalues in eigenvalues_by_spin
    ]

    return {
        'spin_polarized': spin_polarized,
        'nspins': nspins,
        'nkpoints': nkpoints,
        'nbands': nbands,
        'distances': distances,
        'special_distances': special_distances,
        'labels': labels,
        'reference_energy_ev': reference_energy,
        'eigenvalues_by_spin_ev': shifted_by_spin,
        'eigenvalues_ev': (
            shifted_by_spin[0]
            if not spin_polarized
            else None
        ),
        'eigenvalues_up_ev': (
            shifted_by_spin[0]
            if spin_polarized
            else None
        ),
        'eigenvalues_down_ev': (
            shifted_by_spin[1]
            if spin_polarized
            else None
        ),
    }

def parse_dos_output(dos_file):
    """Parse non-spin or collinear-spin QE dos.x data."""
    dos_file = Path(
        dos_file
    )

    if not dos_file.is_file():
        raise FileNotFoundError(
            f"QE DOS data file was not found: {dos_file}"
        )

    energies = []
    dos_values = []
    dos_up = []
    dos_down = []
    integrated_dos = []

    spin_polarized = None

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
                values = [
                    float(
                        value
                        .replace('D', 'E')
                        .replace('d', 'e')
                    )
                    for value in parts
                ]
            except ValueError:
                continue

            row_is_spin_polarized = (
                len(values) >= 4
            )

            if spin_polarized is None:
                spin_polarized = (
                    row_is_spin_polarized
                )

            elif (
                row_is_spin_polarized
                != spin_polarized
            ):
                raise ValueError(
                    "QE DOS data contains inconsistent "
                    "spin column counts."
                )

            energies.append(
                values[0]
            )

            if spin_polarized:
                up_value = values[1]
                down_value = values[2]

                dos_up.append(
                    up_value
                )

                dos_down.append(
                    down_value
                )

                dos_values.append(
                    up_value
                    + down_value
                )

                integrated_dos.append(
                    values[3]
                )

            else:
                dos_values.append(
                    values[1]
                )

                integrated_dos.append(
                    values[2]
                )

    if not energies:
        raise ValueError(
            f"No DOS data could be parsed from '{dos_file}'."
        )

    return {
        'energies_ev': energies,
        'dos': dos_values,
        'dos_up': (
            dos_up
            if spin_polarized
            else None
        ),
        'dos_down': (
            dos_down
            if spin_polarized
            else None
        ),
        'integrated_dos': integrated_dos,
        'spin_polarized': spin_polarized,
        'npoints': len(energies),
    }

def parse_projwfc_pdos_file(pdos_file):
    """Parse one collinear-spin or non-spin QE PDOS file."""
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

    names = component_names[
        orbital
    ]

    nonspin_columns = (
        2
        + len(names)
    )

    spin_columns = (
        3
        + 2 * len(names)
    )

    energies = []
    ldos = []
    ldos_up = []
    ldos_down = []

    components = {
        name: []
        for name in names
    }

    components_up = {
        name: []
        for name in names
    }

    components_down = {
        name: []
        for name in names
    }

    spin_polarized = None

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

            try:
                values = [
                    float(
                        value
                        .replace('D', 'E')
                        .replace('d', 'e')
                    )
                    for value in parts
                ]
            except ValueError:
                continue

            if len(values) == nonspin_columns:
                row_is_spin_polarized = False

            elif len(values) == spin_columns:
                row_is_spin_polarized = True

            else:
                continue

            if spin_polarized is None:
                spin_polarized = (
                    row_is_spin_polarized
                )

            elif (
                row_is_spin_polarized
                != spin_polarized
            ):
                raise ValueError(
                    "QE PDOS data contains inconsistent "
                    "spin column counts."
                )

            energies.append(
                values[0]
            )

            if spin_polarized:
                up_value = values[1]
                down_value = values[2]

                ldos_up.append(
                    up_value
                )

                ldos_down.append(
                    down_value
                )

                ldos.append(
                    up_value
                    + down_value
                )

                for index, component_name in enumerate(
                    names
                ):
                    component_up = values[
                        3 + 2 * index
                    ]

                    component_down = values[
                        4 + 2 * index
                    ]

                    components_up[
                        component_name
                    ].append(
                        component_up
                    )

                    components_down[
                        component_name
                    ].append(
                        component_down
                    )

                    components[
                        component_name
                    ].append(
                        component_up
                        + component_down
                    )

            else:
                ldos.append(
                    values[1]
                )

                for index, component_name in enumerate(
                    names
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
        'ldos_up': (
            ldos_up
            if spin_polarized
            else None
        ),
        'ldos_down': (
            ldos_down
            if spin_polarized
            else None
        ),
        'components': components,
        'components_up': (
            components_up
            if spin_polarized
            else None
        ),
        'components_down': (
            components_down
            if spin_polarized
            else None
        ),
        'spin_polarized': spin_polarized,
        'npoints': len(energies),
    }

def aggregate_projwfc_pdos(pdos_prefix):
    """Aggregate collinear-spin or non-spin QE atomic PDOS files."""
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
    
    spin_modes = {
        parsed['spin_polarized']
        for parsed in parsed_files
    }

    if len(spin_modes) != 1:
        raise ValueError(
            "QE PDOS files contain inconsistent "
            "spin configurations."
        )

    spin_polarized = (
        spin_modes.pop()
    )

    spin_up = None
    spin_down = None

    if spin_polarized:
        spin_totals_up = {
            orbital: [0.0] * len(energies)
            for orbital in totals
        }

        spin_totals_down = {
            orbital: [0.0] * len(energies)
            for orbital in totals
        }

        spin_components_up = {
            component: [0.0] * len(energies)
            for component in components
        }

        spin_components_down = {
            component: [0.0] * len(energies)
            for component in components
        }

        for parsed in parsed_files:
            orbital = parsed[
                'orbital'
            ]

            for index, value in enumerate(
                parsed['ldos_up']
            ):
                spin_totals_up[
                    orbital
                ][index] += value

            for index, value in enumerate(
                parsed['ldos_down']
            ):
                spin_totals_down[
                    orbital
                ][index] += value

            channel_components = (
                (
                    parsed['components_up'],
                    spin_components_up,
                ),
                (
                    parsed['components_down'],
                    spin_components_down,
                ),
            )

            for source, target in channel_components:
                for component_name, values in (
                    source.items()
                ):
                    if (
                        orbital == 's'
                        and component_name == 's'
                    ):
                        continue

                    if component_name not in target:
                        raise ValueError(
                            "Unexpected QE PDOS component: "
                            f"{component_name}"
                        )

                    for index, value in enumerate(
                        values
                    ):
                        target[
                            component_name
                        ][index] += value

        def build_spin_projection(
            spin_totals,
            spin_components,
        ):
            projected = []

            for index in range(
                len(energies)
            ):
                projected.append(
                    spin_totals['s'][index]
                    + spin_totals['p'][index]
                    + spin_totals['d'][index]
                    + spin_totals['f'][index]
                )

            return {
                's_total': spin_totals['s'],
                'p_total': spin_totals['p'],
                'pz': spin_components['pz'],
                'px': spin_components['px'],
                'py': spin_components['py'],
                'd_total': spin_totals['d'],
                'd3z2_r2': spin_components[
                    'd3z2_r2'
                ],
                'dxz': spin_components['dxz'],
                'dyz': spin_components['dyz'],
                'dx2_y2': spin_components[
                    'dx2_y2'
                ],
                'dxy': spin_components['dxy'],
                'f_total': spin_totals['f'],
                'total': projected,
            }

        spin_up = build_spin_projection(
            spin_totals_up,
            spin_components_up,
        )

        spin_down = build_spin_projection(
            spin_totals_down,
            spin_components_down,
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
        'spin_polarized': spin_polarized,
        'spin_up': spin_up,
        'spin_down': spin_down,
    }

def run_ph(
    input_file,
    output_file,
    state_dir,
    fildyn,
    qpoint_grid,
    tr2_ph=1.0e-12,
    parallel_cores=1,
    executable='ph.x',
    prefix='nanoworks',
):
    """Render, execute, and validate one QE ph.x grid calculation."""
    input_file = Path(input_file)
    output_file = Path(output_file)
    state_dir = Path(state_dir)
    fildyn = Path(fildyn)

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE ground-state directory is required "
            f"for ph.x: {state_dir}"
        )

    input_text = render_ph_input(
        prefix=prefix,
        outdir=state_dir,
        fildyn=fildyn,
        qpoint_grid=qpoint_grid,
        tr2_ph=tr2_ph,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    fildyn.parent.mkdir(
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

    result = parse_qe_auxiliary_output(
        output_file,
        expected_program='PHONON',
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
            "Quantum ESPRESSO ph.x finished without a "
            f"'JOB DONE.' marker. See '{output_file}'."
        )

    grid_file = Path(
        str(fildyn) + '0'
    )

    if not grid_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO ph.x did not produce the q-grid "
            f"metadata file: {grid_file}"
        )

    dynamical_matrix_files = sorted(
        path
        for path in fildyn.parent.glob(
            fildyn.name + '*'
        )
        if path.is_file() and path != grid_file
    )

    if not dynamical_matrix_files:
        raise RuntimeError(
            "Quantum ESPRESSO ph.x did not produce any "
            "dynamical-matrix files for the requested q-grid."
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'fildyn': fildyn,
        'grid_file': grid_file,
        'dynamical_matrix_files': dynamical_matrix_files,
        'execution': execution,
        'result': result,
    }


def run_q2r(
    input_file,
    output_file,
    fildyn,
    flfrc,
    zasr='no',
    parallel_cores=1,
    executable='q2r.x',
):
    """Render, execute, and validate one QE q2r.x calculation."""
    input_file = Path(input_file)
    output_file = Path(output_file)
    fildyn = Path(fildyn)
    flfrc = Path(flfrc)

    grid_file = Path(
        str(fildyn) + '0'
    )

    if not grid_file.is_file():
        raise FileNotFoundError(
            "QE q2r.x requires the ph.x q-grid metadata file: "
            f"{grid_file}"
        )

    input_text = render_q2r_input(
        fildyn=fildyn,
        flfrc=flfrc,
        zasr=zasr,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    flfrc.parent.mkdir(
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

    result = parse_qe_auxiliary_output(
        output_file,
        expected_program='Q2R',
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
            "Quantum ESPRESSO q2r.x finished without a "
            f"'JOB DONE.' marker. See '{output_file}'."
        )

    if not flfrc.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO q2r.x did not produce the "
            f"real-space force-constant file: {flfrc}"
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'fildyn': fildyn,
        'grid_file': grid_file,
        'flfrc': flfrc,
        'execution': execution,
        'result': result,
    }


def run_matdyn_band(
    input_file,
    output_file,
    flfrc,
    flfrq,
    band_path,
    acoustic_sum_rule=True,
    parallel_cores=1,
    executable='matdyn.x',
):
    """Render, execute, and validate a QE phonon band calculation."""
    input_file = Path(input_file)
    output_file = Path(output_file)
    flfrc = Path(flfrc)
    flfrq = Path(flfrq)

    if not flfrc.is_file():
        raise FileNotFoundError(
            "QE matdyn.x requires a real-space force-constant "
            f"file: {flfrc}"
        )

    input_text = render_matdyn_band_input(
        flfrc=flfrc,
        flfrq=flfrq,
        band_path=band_path,
        acoustic_sum_rule=acoustic_sum_rule,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    flfrq.parent.mkdir(
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

    result = parse_qe_auxiliary_output(
        output_file,
        expected_program='MATDYN',
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
            "Quantum ESPRESSO matdyn.x phonon band calculation "
            "finished without a 'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    if not flfrq.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO matdyn.x did not produce the "
            f"phonon frequency file: {flfrq}"
        )

    frequencies = parse_matdyn_frequency_file(
        flfrq
    )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'flfrc': flfrc,
        'flfrq': flfrq,
        'band_path': band_path,
        'frequencies': frequencies,
        'execution': execution,
        'result': result,
    }


def run_matdyn_dos(
    input_file,
    output_file,
    flfrc,
    fldos,
    qpoint_grid,
    acoustic_sum_rule=True,
    parallel_cores=1,
    executable='matdyn.x',
):
    """Render, execute, and validate a QE phonon DOS calculation."""
    input_file = Path(input_file)
    output_file = Path(output_file)
    flfrc = Path(flfrc)
    fldos = Path(fldos)
    qpoint_grid = tuple(qpoint_grid)

    if not flfrc.is_file():
        raise FileNotFoundError(
            "QE matdyn.x requires a real-space force-constant "
            f"file: {flfrc}"
        )

    input_text = render_matdyn_dos_input(
        flfrc=flfrc,
        fldos=fldos,
        qpoint_grid=qpoint_grid,
        acoustic_sum_rule=acoustic_sum_rule,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    fldos.parent.mkdir(
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

    result = parse_qe_auxiliary_output(
        output_file,
        expected_program='MATDYN',
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
            "Quantum ESPRESSO matdyn.x phonon DOS calculation "
            "finished without a 'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    if not fldos.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO matdyn.x did not produce the "
            f"phonon DOS file: {fldos}"
        )

    dos = parse_matdyn_dos_file(
        fldos
    )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'flfrc': flfrc,
        'fldos': fldos,
        'qpoint_grid': qpoint_grid,
        'dos': dos,
        'execution': execution,
        'result': result,
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
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
    occupation=None,
    parallel_cores=1,
    executable='pw.x',
    prefix='nanoworks',
    exx_additional_kpoints=None,
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
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
        exx_additional_kpoints=exx_additional_kpoints,
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

def run_relax(
    atoms,
    input_file,
    output_file,
    state_dir,
    pseudopotentials,
    pseudo_dir,
    cutoff_ev,
    optimizer,
    max_force,
    max_step,
    relax_cell,
    hydrostatic_pressure=0.0,
    fix_symmetry=False,
    kpoint_density=None,
    kpoint_size=(5, 5, 5),
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
    occupation=None,
    parallel_cores=1,
    executable='pw.x',
    prefix='nanoworks',
):
    """Render, execute, and parse one QE geometry optimization."""
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

    occupation_settings = resolve_qe_occupation(
        occupation
    )

    relaxation_settings = resolve_qe_relaxation_settings(
        optimizer=optimizer,
        max_force=max_force,
        max_step=max_step,
        relax_cell=relax_cell,
        hydrostatic_pressure=hydrostatic_pressure,
        fix_symmetry=fix_symmetry,
        atoms=atoms,
    )

    input_text = render_pw_input(
        calculation=relaxation_settings['calculation'],
        atoms=atoms,
        pseudopotentials=pseudopotentials,
        cutoff_ev=cutoff_ev,
        kpoint_size=mesh,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
        occupations=occupation_settings['occupations'],
        smearing=occupation_settings['smearing'],
        width_ev=occupation_settings['width_ev'],
        prefix=prefix,
        pseudo_dir=pseudo_dir,
        outdir=state_dir,
        relaxation_settings=relaxation_settings,
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
            "Quantum ESPRESSO geometry optimization finished "
            "without a 'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    output_text = output_file.read_text(
        encoding='utf-8',
        errors='replace',
    )

    geometry_converged = bool(
        re.search(
            r'bfgs\s+converged\s+in',
            output_text,
            flags=re.IGNORECASE,
        )
    )

    if not geometry_converged:
        raise RuntimeError(
            "Quantum ESPRESSO geometry optimization did not "
            "report BFGS convergence. "
            f"See '{output_file}'."
        )

    relaxed_atoms = parse_pw_relaxed_structure(
        output_file,
        atoms,
    )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'kpoint_size': mesh,
        'calculation': relaxation_settings['calculation'],
        'relaxation_settings': relaxation_settings,
        'execution': execution,
        'result': result,
        'geometry_converged': geometry_converged,
        'atoms': relaxed_atoms,
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
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
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
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
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
    magnetic_moments=None,
    setup_params=None,
    xc_calc='PBE',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
    occupation=None,
    parallel_cores=1,
    executable='pw.x',
    prefix='nanoworks',
    projected_band=False,
    projections=None,
    projection_input_file=None,
    projection_output_file=None,
    projection_prefix=None,
    projection_executable='projwfc.x',
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
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
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

    band_projections = None

    if projected_band:
        projection_paths = {
            'projection_input_file': (
                projection_input_file
            ),
            'projection_output_file': (
                projection_output_file
            ),
            'projection_prefix': (
                projection_prefix
            ),
        }

        missing_paths = [
            name
            for name, value in (
                projection_paths.items()
            )
            if value is None
        ]

        if missing_paths:
            raise ValueError(
                "QE projected bands require: "
                + ', '.join(
                    missing_paths
                )
            )

        projection_workflow = (
            run_band_projections(
                input_file=(
                    projection_input_file
                ),
                output_file=(
                    projection_output_file
                ),
                state_dir=state_dir,
                projection_prefix=(
                    projection_prefix
                ),
                spinpol=spinpol,
                parallel_cores=parallel_cores,
                executable=(
                    projection_executable
                ),
                prefix=prefix,
            )
        )

        raw_up = parse_projwfc_band_file(
            projection_workflow[
                'projection_up_file'
            ]
        )

        if (
            raw_up['nkpoints']
            != bands['nkpoints']
            or raw_up['nbands']
            != bands['nbands']
        ):
            raise RuntimeError(
                "QE spin-up band projections do not "
                "match the calculated band dimensions."
            )

        prepared_up = (
            prepare_qe_band_projection_data(
                raw_up,
                projections=projections,
            )
        )

        raw_down = None
        prepared_down = None

        if spinpol:
            raw_down = (
                parse_projwfc_band_file(
                    projection_workflow[
                        'projection_down_file'
                    ]
                )
            )

            if (
                raw_down['nkpoints']
                != bands['nkpoints']
                or raw_down['nbands']
                != bands['nbands']
            ):
                raise RuntimeError(
                    "QE spin-down band projections do not "
                    "match the calculated band dimensions."
                )

            prepared_down = (
                prepare_qe_band_projection_data(
                    raw_down,
                    projections=projections,
                )
            )

        band_projections = {
            'workflow': projection_workflow,
            'spin_polarized': bool(
                spinpol
            ),
            'raw_up': raw_up,
            'raw_down': raw_down,
            'up': prepared_up,
            'down': prepared_down,
        }

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'band_path': band_path,
        'execution': execution,
        'result': result,
        'bands': bands,
        'band_projections': band_projections,
    }

def run_bands_postprocess(
    input_file,
    output_file,
    state_dir,
    band_file,
    parallel_cores=1,
    executable='bands.x',
    prefix='nanoworks',
    lsym=False,
):
    """Render and execute one Quantum ESPRESSO bands.x calculation."""
    input_file = Path(
        input_file
    )
    output_file = Path(
        output_file
    )
    state_dir = Path(
        state_dir
    )
    band_file = Path(
        band_file
    )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE electronic state is required "
            f"for bands.x: {state_dir}"
        )

    input_text = render_bands_postprocess_input(
        prefix=prefix,
        outdir=state_dir,
        filband=band_file,
        lsym=lsym,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    band_file.parent.mkdir(
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

    result = parse_qe_auxiliary_output(
        output_file,
        expected_program='BANDS',
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
            "Quantum ESPRESSO bands.x finished without a "
            f"'JOB DONE.' marker. See '{output_file}'."
        )

    if not band_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO bands.x finished but the band data "
            f"file was not created: {band_file}"
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'band_file': band_file,
        'execution': execution,
        'result': result,
    }

def run_hybrid_bands(
    atoms,
    scf_input_file,
    scf_output_file,
    bands_input_file,
    bands_output_file,
    state_dir,
    band_file,
    pseudopotentials,
    pseudo_dir,
    cutoff_ev,
    band_path,
    qpoint_grid,
    kpoint_density=None,
    kpoint_size=(5, 5, 5),
    gamma=False,
    total_charge=0.0,
    nbands=None,
    spinpol=False,
    magnetic_moments=None,
    setup_params=None,
    xc_calc='HSE06',
    pseudo_xc='pbe',
    exx_fraction=None,
    omega=None,
    occupation=None,
    parallel_cores=1,
    scf_executable='pw.x',
    bands_executable='bands.x',
    prefix='nanoworks',
    lsym=False,
):
    """Run a QE hybrid SCF followed by bands.x post-processing."""
    xc_settings = resolve_qe_xc_settings(
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
    )

    if not xc_settings['hybrid']:
        raise ValueError(
            "QE hybrid bands require a hybrid functional."
        )

    additional_kpoints = (
        build_qe_exx_additional_kpoints(
            band_path=band_path,
            qpoint_grid=qpoint_grid,
        )
    )

    scf_workflow = run_scf(
        atoms=atoms,
        input_file=scf_input_file,
        output_file=scf_output_file,
        state_dir=state_dir,
        pseudopotentials=pseudopotentials,
        pseudo_dir=pseudo_dir,
        cutoff_ev=cutoff_ev,
        kpoint_density=kpoint_density,
        kpoint_size=kpoint_size,
        gamma=gamma,
        total_charge=total_charge,
        nbands=nbands,
        spinpol=spinpol,
        magnetic_moments=magnetic_moments,
        setup_params=setup_params,
        xc_calc=xc_calc,
        pseudo_xc=pseudo_xc,
        exx_fraction=exx_fraction,
        omega=omega,
        occupation=occupation,
        parallel_cores=parallel_cores,
        executable=scf_executable,
        prefix=prefix,
        exx_additional_kpoints=additional_kpoints,
    )

    bands_workflow = run_bands_postprocess(
        input_file=bands_input_file,
        output_file=bands_output_file,
        state_dir=state_dir,
        band_file=band_file,
        parallel_cores=parallel_cores,
        executable=bands_executable,
        prefix=prefix,
        lsym=lsym,
    )

    return {
        'scf': scf_workflow,
        'bands': bands_workflow,
        'band_path': band_path,
        'additional_kpoints': additional_kpoints,
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

def run_pp_density(
    input_file,
    output_file,
    state_dir,
    filplot,
    cube_file,
    plot_num=0,
    spin_component=None,
    parallel_cores=1,
    executable='pp.x',
    prefix='nanoworks',
):
    """Render and execute one QE pp.x density calculation."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    filplot = Path(
        filplot
    )

    cube_file = Path(
        cube_file
    )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE ground-state result is required "
            f"for the density calculation: {state_dir}"
        )

    input_text = render_pp_input(
        prefix=prefix,
        outdir=state_dir,
        filplot=filplot,
        fileout=cube_file,
        plot_num=plot_num,
        spin_component=spin_component,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    filplot.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    cube_file.parent.mkdir(
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

    output_text = output_file.read_text(
        encoding='utf-8',
        errors='replace',
    )

    if 'JOB DONE.' not in output_text:
        raise RuntimeError(
            "Quantum ESPRESSO density calculation "
            "finished without a 'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    if not cube_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO pp.x finished but the "
            "density Cube file was not created: "
            f"{cube_file}"
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'filplot': filplot,
        'cube_file': cube_file,
        'plot_num': int(
            plot_num
        ),
        'spin_component': (
            int(spin_component)
            if spin_component is not None
            else None
        ),
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
    
def run_band_projections(
    input_file,
    output_file,
    state_dir,
    projection_prefix,
    pdos_prefix=None,
    spinpol=False,
    parallel_cores=1,
    executable='projwfc.x',
    prefix='nanoworks',
):
    """Run projwfc.x for QE band-projection data."""
    input_file = Path(
        input_file
    )

    output_file = Path(
        output_file
    )

    state_dir = Path(
        state_dir
    )

    projection_prefix = Path(
        projection_prefix
    )

    if pdos_prefix is None:
        pdos_prefix = Path(
            str(projection_prefix)
            + '-pdos'
        )
    else:
        pdos_prefix = Path(
            pdos_prefix
        )

    if not has_qe_state(
        state_dir,
        prefix=prefix,
    ):
        raise FileNotFoundError(
            "A valid QE bands state is required "
            f"for band projections: {state_dir}"
        )

    input_text = render_projwfc_input(
        prefix=prefix,
        outdir=state_dir,
        filpdos=pdos_prefix,
        filproj=projection_prefix,
        lsym=False,
        diag_basis=False,
    )

    input_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    projection_prefix.parent.mkdir(
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
            "Quantum ESPRESSO band projection calculation "
            "finished without a 'JOB DONE.' marker. "
            f"See '{output_file}'."
        )

    projection_up_file = Path(
        str(projection_prefix)
        + '.projwfc_up'
    )

    projection_down_file = Path(
        str(projection_prefix)
        + '.projwfc_down'
    )

    if not projection_up_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO projwfc.x finished but the "
            "band projection file was not created: "
            f"{projection_up_file}"
        )

    if spinpol and not projection_down_file.is_file():
        raise RuntimeError(
            "Quantum ESPRESSO spin-polarized projwfc.x "
            "finished but the spin-down projection file "
            f"was not created: {projection_down_file}"
        )

    projection_files = [
        projection_up_file,
    ]

    if spinpol:
        projection_files.append(
            projection_down_file
        )

    return {
        'input_file': input_file,
        'output_file': output_file,
        'state_dir': state_dir,
        'projection_prefix': projection_prefix,
        'projection_up_file': projection_up_file,
        'projection_down_file': (
            projection_down_file
            if spinpol
            else None
        ),
        'projection_files': projection_files,
        'pdos_prefix': pdos_prefix,
        'spin_polarized': bool(
            spinpol
        ),
        'execution': execution,
    }
