#!/usr/bin/env python3
"""Convert a Quantum ESPRESSO pw.x input into Nanoworks QE files.

The script reads a pw.x style input file, extracts common calculation
parameters and lattice/atomic structure, then produces:
  * a CIF geometry file for use with dftsolve's ``-g`` option.
  * a Python configuration module defining the dftsolve.py variables.

Example
-------
qeconverter --input si.scf.in --output-dir example_folder --system-name SiliconQE
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional
import logging
import nanoworks
from typing import Union
from ase.io import read, write
from ase.units import Bohr

RY_TO_EV = 13.605693009

logger = logging.getLogger(__name__)


@dataclass
class QEInputSettings:
    calculation: Optional[str] = None
    ecutwfc: Optional[float] = None
    occupations: Optional[str] = None
    smearing: Optional[str] = None
    degauss: Optional[float] = None
    nspin: Optional[int] = None
    starting_magnetization: Dict[str, float] = field(default_factory=dict)
    conv_thr: Optional[float] = None
    k_mesh: Optional[List[int]] = None
    k_shift: Optional[List[int]] = None
    total_charge: Optional[float] = None
    nbands: Optional[int] = None
    forc_conv_thr: Optional[float] = None
    ion_dynamics: Optional[str] = None
    trust_radius_max: Optional[float] = None
    cell_dynamics: Optional[str] = None
    cell_dofree: Optional[str] = None
    pressure_kbar: Optional[float] = None
    nosym: Optional[bool] = None
    nat: Optional[int] = None
    ntyp: Optional[int] = None
    pseudo_dir: Optional[str] = None
    species_labels: List[str] = field(
        default_factory=list
    )
    species_pseudopotentials: Dict[str, str] = field(
        default_factory=dict
    )
    atomic_position_labels: List[str] = field(
        default_factory=list
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert Quantum ESPRESSO pw.x input into dftsolve.py inputs.",
    )
    parser.add_argument("--input", type=Path, help="Path to pw.x input file")
    parser.add_argument("--output-dir", type=Path, default=Path.cwd(), help="Directory for generated files")
    parser.add_argument("--system-name", help="System name used for file stems and Outdirname")
    parser.add_argument("--outdirname", help="Override Outdirname value inside the generated input")
    parser.add_argument("--input-filename", help="Override dftsolve.py input filename")
    parser.add_argument("--xc", help="Optional XC functional override (default PBE)")
    parser.add_argument("-v", "--version", action='store_true', help="Show version information")
    return parser.parse_args()


def clean_line(line: str) -> str:
    return line.split("!")[0].split("#")[0].strip()

def _split_qe_assignments(
    line: str,
) -> List[str]:
    """Split comma-separated QE namelist assignments safely."""
    assignments = []
    current = []
    quote = None
    parenthesis_depth = 0

    for character in line:
        if quote is not None:
            current.append(
                character
            )

            if character == quote:
                quote = None

            continue

        if character in {
            "'",
            '"',
        }:
            quote = character
            current.append(
                character
            )
            continue

        if character == '(':
            parenthesis_depth += 1
            current.append(
                character
            )
            continue

        if character == ')':
            parenthesis_depth = max(
                0,
                parenthesis_depth - 1,
            )
            current.append(
                character
            )
            continue

        if (
            character == ','
            and parenthesis_depth == 0
        ):
            assignment = ''.join(
                current
            ).strip()

            if assignment:
                assignments.append(
                    assignment
                )

            current = []
            continue

        current.append(
            character
        )

    assignment = ''.join(
        current
    ).strip()

    if assignment:
        assignments.append(
            assignment
        )

    return assignments


def _parse_qe_float(
    value: str,
) -> float:
    """Parse a QE floating-point value, including Fortran D exponents."""
    return float(
        value
        .strip()
        .replace('D', 'E')
        .replace('d', 'e')
    )

def _parse_qe_bool(
    value: str,
) -> bool:
    """Parse a Quantum ESPRESSO logical value."""
    normalized = (
        value
        .strip()
        .lower()
    )

    true_values = {
        '.true.',
        'true',
        't',
        '1',
    }

    false_values = {
        '.false.',
        'false',
        'f',
        '0',
    }

    if normalized in true_values:
        return True

    if normalized in false_values:
        return False

    raise ValueError(
        "Unable to parse QE logical value: "
        f"{value}"
    )

def _normalize_qe_smearing(
    smearing: Optional[str],
) -> str:
    """Translate QE smearing aliases to Nanoworks names."""
    name = (
        'gaussian'
        if smearing is None
        else smearing.strip().lower()
    )

    aliases = {
        'gauss': 'gaussian',
        'gaussian': 'gaussian',
        'mp': 'methfessel-paxton',
        'methfessel-paxton': 'methfessel-paxton',
        'methfessel_paxton': 'methfessel-paxton',
        'mv': 'marzari-vanderbilt',
        'm-v': 'marzari-vanderbilt',
        'cold': 'marzari-vanderbilt',
        'marzari-vanderbilt': 'marzari-vanderbilt',
        'marzari_vanderbilt': 'marzari-vanderbilt',
        'fd': 'fermi-dirac',
        'fermi-dirac': 'fermi-dirac',
        'fermi_dirac': 'fermi-dirac',
    }

    try:
        return aliases[
            name
        ]
    except KeyError as exc:
        raise ValueError(
            "Unsupported QE smearing scheme: "
            f"{smearing}"
        ) from exc


def _build_occupation_line(
    settings: QEInputSettings,
) -> str:
    """Build a Nanoworks occupation setting from QE input."""
    occupation = (
        'fixed'
        if settings.occupations is None
        else settings.occupations.strip().lower()
    )

    if occupation == 'fixed':
        return "Occupation = 'fixed'"

    tetrahedra = {
        'tetrahedra',
        'tetrahedra_lin',
        'tetrahedra_opt',
    }

    if occupation in tetrahedra:
        return (
            f"Occupation = '{occupation}'"
        )

    if occupation != 'smearing':
        raise ValueError(
            "Unsupported QE occupation scheme: "
            f"{settings.occupations}"
        )

    if settings.degauss is None:
        raise ValueError(
            "QE occupations='smearing' requires "
            "an explicit degauss value for conversion."
        )

    smearing = _normalize_qe_smearing(
        settings.smearing
    )

    width_ev = (
        settings.degauss
        * RY_TO_EV
    )

    return (
        "Occupation = {"
        f"'name': '{smearing}', "
        f"'width': {width_ev:.12g}"
        "}"
    )

def _build_kpoint_lines(
    settings: QEInputSettings,
) -> List[str]:
    """Convert a QE automatic k-point mesh to Nanoworks settings."""
    if settings.k_mesh is None:
        return []

    mesh = list(
        settings.k_mesh
    )

    mesh.extend([
        1,
        1,
        1,
    ])

    mesh = [
        int(
            value
        )
        for value in mesh[:3]
    ]

    input_shift = list(
        settings.k_shift
        or [
            0,
            0,
            0,
        ]
    )

    input_shift.extend([
        0,
        0,
        0,
    ])

    input_shift = [
        int(
            value
        )
        for value in input_shift[:3]
    ]

    notices = []

    if any(
        value not in {
            0,
            1,
        }
        for value in input_shift
    ):
        original_shift = list(
            input_shift
        )

        input_shift = [
            0
            if value == 0
            else 1
            for value in input_shift
        ]

        notices.append(
            "QE k-point shift "
            f"{original_shift} contains values other "
            "than 0 or 1; it is normalized to "
            f"{input_shift}."
        )

    gamma_shift = [
        0,
        0,
        0,
    ]

    shifted_mesh_shift = [
        1
        if value % 2 == 0
        else 0
        for value in mesh
    ]

    if input_shift == gamma_shift:
        ground_gamma = True
        output_shift = gamma_shift

    elif input_shift == shifted_mesh_shift:
        ground_gamma = False
        output_shift = shifted_mesh_shift

    else:
        gamma_distance = sum(
            abs(
                input_value
                - output_value
            )
            for input_value, output_value in zip(
                input_shift,
                gamma_shift,
            )
        )

        shifted_distance = sum(
            abs(
                input_value
                - output_value
            )
            for input_value, output_value in zip(
                input_shift,
                shifted_mesh_shift,
            )
        )

        if gamma_distance <= shifted_distance:
            ground_gamma = True
            output_shift = gamma_shift

        else:
            ground_gamma = False
            output_shift = shifted_mesh_shift

        notices.append(
            "QE k-point shift "
            f"{input_shift} cannot be represented "
            "exactly by Ground_gamma; the nearest "
            f"available shift {output_shift} is used."
        )

    lines = [
        f"Ground_kpts_x = {mesh[0]}",
        f"Ground_kpts_y = {mesh[1]}",
        f"Ground_kpts_z = {mesh[2]}",
    ]

    lines.extend(
        "# NOTICE: "
        + notice
        for notice in notices
    )

    lines.append(
        "Ground_gamma = "
        + repr(
            ground_gamma
        )
    )

    return lines

def _build_workflow_lines(
    settings: QEInputSettings,
) -> List[str]:
    """Translate a QE calculation type to Nanoworks workflow flags."""
    calculation = (
        settings.calculation
        or 'scf'
    ).strip().lower()

    ground_calc = True
    geo_optim = False
    dos_calc = False
    band_calc = False
    notices = []

    if calculation == 'scf':
        pass

    elif calculation in {
        'relax',
        'vc-relax',
    }:
        geo_optim = True

    elif calculation == 'nscf':
        dos_calc = True

        notices.append(
            "QE calculation = 'nscf' is interpreted as "
            "a DOS workflow. Review DOS-specific settings "
            "if the NSCF calculation had another purpose."
        )

    elif calculation == 'bands':
        band_calc = True

        notices.append(
            "QE calculation = 'bands' enables the Nanoworks "
            "band workflow. The QE explicit k-point path is "
            "not converted yet, so review Band_path."
        )

    else:
        notices.append(
            "QE calculation = "
            f"'{calculation}' has no direct Nanoworks "
            "workflow mapping; a ground-state calculation "
            "is enabled as the fallback."
        )

    lines = [
        f"Ground_calc = {ground_calc}",
        f"Geo_optim = {geo_optim}",
        "Elastic_calc = False",
        f"DOS_calc = {dos_calc}",
        f"Band_calc = {band_calc}",
        "Density_calc = False",
        "Optical_calc = False",
    ]

    lines.extend(
        "# NOTICE: "
        + notice
        for notice in notices
    )

    return lines

def _build_relaxation_lines(
    settings: QEInputSettings,
) -> List[str]:
    """Convert QE relaxation settings to Nanoworks settings."""
    calculation = (
        settings.calculation
        or 'scf'
    ).strip().lower()

    if calculation not in {
        'relax',
        'vc-relax',
    }:
        return []

    notices = []
    lines = []

    if settings.ion_dynamics not in {
        None,
        'bfgs',
    }:
        notices.append(
            "QE ion_dynamics = "
            f"'{settings.ion_dynamics}' is approximated "
            "with Nanoworks LBFGS."
        )

    if (
        calculation == 'vc-relax'
        and settings.cell_dynamics not in {
            None,
            'bfgs',
        }
    ):
        notices.append(
            "QE cell_dynamics = "
            f"'{settings.cell_dynamics}' is approximated "
            "with QE BFGS through Nanoworks."
        )

    lines.append(
        "Optimizer = 'LBFGS'"
    )

    if settings.forc_conv_thr is not None:
        max_force_ev_angstrom = (
            settings.forc_conv_thr
            * RY_TO_EV
            / Bohr
        )

        lines.append(
            "Max_F_tolerance = "
            f"{max_force_ev_angstrom:.12g}"
        )

    if settings.trust_radius_max is not None:
        max_step_angstrom = (
            settings.trust_radius_max
            * Bohr
        )

        lines.append(
            "Max_step = "
            f"{max_step_angstrom:.12g}"
        )

    fix_symmetry = (
        not settings.nosym
        if settings.nosym is not None
        else True
    )

    if calculation == 'relax':
        relax_cell = [
            False,
            False,
            False,
            False,
            False,
            False,
        ]

        if settings.cell_dofree is not None:
            notices.append(
                "QE cell_dofree is ignored because "
                "calculation = 'relax' keeps the cell fixed."
            )

        if (
            settings.pressure_kbar is not None
            and settings.pressure_kbar != 0.0
        ):
            notices.append(
                "QE press is ignored because calculation = "
                "'relax' does not relax the cell."
            )

    else:
        cell_dofree = (
            settings.cell_dofree
            or 'all'
        ).strip().lower()

        exact_mappings = {
            'x': [
                True,
                False,
                False,
                False,
                False,
                False,
            ],
            'y': [
                False,
                True,
                False,
                False,
                False,
                False,
            ],
            'z': [
                False,
                False,
                True,
                False,
                False,
                False,
            ],
            'xy': [
                True,
                True,
                False,
                False,
                False,
                False,
            ],
            'xz': [
                True,
                False,
                True,
                False,
                False,
                False,
            ],
            'yz': [
                False,
                True,
                True,
                False,
                False,
                False,
            ],
            'xyz': [
                True,
                True,
                True,
                False,
                False,
                False,
            ],
            '2dxy': [
                True,
                True,
                False,
                False,
                False,
                True,
            ],
            'all': [
                True,
                True,
                True,
                True,
                True,
                True,
            ],
        }

        approximate_mappings = {
            'volume': [
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            'shape': [
                True,
                True,
                True,
                True,
                True,
                True,
            ],
            '2dshape': [
                True,
                True,
                False,
                False,
                False,
                True,
            ],
            'ibrav': [
                True,
                True,
                True,
                True,
                True,
                True,
            ],
        }

        if cell_dofree in exact_mappings:
            relax_cell = exact_mappings[
                cell_dofree
            ]

        elif cell_dofree in approximate_mappings:
            relax_cell = approximate_mappings[
                cell_dofree
            ]

            notices.append(
                "QE cell_dofree = "
                f"'{settings.cell_dofree}' cannot be "
                "represented exactly by Relax_cell; "
                "the nearest available mask is used."
            )

            if cell_dofree in {
                'volume',
                'ibrav',
            }:
                fix_symmetry = True

                notices.append(
                    "Fix_symmetry = True is used to keep "
                    "the approximate cell relaxation close "
                    "to the original QE constraint."
                )

        else:
            relax_cell = [
                True,
                True,
                True,
                True,
                True,
                True,
            ]

            notices.append(
                "Unknown QE cell_dofree = "
                f"'{settings.cell_dofree}'; all cell "
                "components are enabled as a fallback."
            )

    for notice in notices:
        lines.append(
            "# NOTICE: "
            + notice
        )

    lines.append(
        "Relax_cell = "
        + repr(
            relax_cell
        )
    )

    lines.append(
        "Fix_symmetry = "
        + repr(
            fix_symmetry
        )
    )

    if calculation == 'vc-relax':
        hydrostatic_pressure = (
            0.0
            if settings.pressure_kbar is None
            else settings.pressure_kbar / 10.0
        )

        lines.append(
            "Hydrostatic_pressure = "
            f"{hydrostatic_pressure:.12g}"
        )

    return lines

def parse_qe_input(
    path: Path,
) -> QEInputSettings:
    settings = QEInputSettings()

    lines = [
        clean_line(
            line
        )
        for line in path.read_text(
            encoding='utf-8',
        ).splitlines()
    ]

    lines = [
        line
        for line in lines
        if line
    ]

    expect_automatic_kpoints = False
    remaining_species_rows = 0
    remaining_position_rows = 0

    for line in lines:
        upper = line.upper()

        if remaining_species_rows > 0:
            fields = line.split()

            if len(fields) >= 3:
                species_label = fields[
                    0
                ]

                pseudopotential = fields[
                    2
                ]

                settings.species_labels.append(
                    species_label
                )

                settings.species_pseudopotentials[
                    species_label
                ] = pseudopotential

                remaining_species_rows -= 1

                continue

            logger.warning(
                "Unable to parse ATOMIC_SPECIES "
                "row %r; leaving it unresolved",
                line,
            )

            remaining_species_rows = 0

        if remaining_position_rows > 0:
            fields = line.split()

            if len(fields) >= 4:
                settings.atomic_position_labels.append(
                    fields[
                        0
                    ]
                )

                remaining_position_rows -= 1

                continue

            logger.warning(
                "Unable to parse ATOMIC_POSITIONS "
                "row %r; leaving it unresolved",
                line,
            )

            remaining_position_rows = 0

        if upper.startswith('&'):
            continue

        if upper == '/':
            continue

        if upper.startswith(
            'ATOMIC_SPECIES'
        ):
            if settings.ntyp is None:
                logger.warning(
                    "ATOMIC_SPECIES cannot be parsed "
                    "because ntyp was not found."
                )
            else:
                remaining_species_rows = (
                    settings.ntyp
                )

            continue

        if upper.startswith(
            'ATOMIC_POSITIONS'
        ):
            if settings.nat is None:
                logger.warning(
                    "ATOMIC_POSITIONS labels cannot be "
                    "parsed because nat was not found."
                )
            else:
                remaining_position_rows = (
                    settings.nat
                )

            continue
            
        if upper.startswith('K_POINTS'):
            parts = line.split()

            k_mode = (
                parts[1]
                if len(parts) > 1
                else 'automatic'
            )

            k_mode = (
                k_mode
                .strip()
                .lower()
                .strip('{}()')
            )

            if k_mode == 'gamma':
                settings.k_mesh = [
                    1,
                    1,
                    1,
                ]
                settings.k_shift = [
                    0,
                    0,
                    0,
                ]
                expect_automatic_kpoints = False

            elif k_mode == 'automatic':
                expect_automatic_kpoints = True

            else:
                expect_automatic_kpoints = False

            continue

        if expect_automatic_kpoints:
            expect_automatic_kpoints = False

            tokens = line.split()

            if len(tokens) < 3:
                raise ValueError(
                    "QE automatic K_POINTS card requires "
                    "at least three mesh values."
                )

            try:
                settings.k_mesh = [
                    int(
                        float(
                            token
                        )
                    )
                    for token in tokens[:3]
                ]

                settings.k_shift = (
                    [
                        int(
                            float(
                                token
                            )
                        )
                        for token in tokens[3:6]
                    ]
                    if len(tokens) >= 6
                    else [
                        0,
                        0,
                        0,
                    ]
                )

            except ValueError as exc:
                raise ValueError(
                    "Unable to parse QE automatic "
                    f"K_POINTS row: {line}"
                ) from exc

            continue

        if '=' not in line:
            continue

        for assignment in _split_qe_assignments(
            line
        ):
            if '=' not in assignment:
                continue

            key, value = [
                part.strip()
                for part in assignment.split(
                    '=',
                    1,
                )
            ]

            key_lower = key.lower()

            if key_lower.startswith(
                'starting_magnetization'
            ):
                species_match = re.search(
                    r'starting_magnetization'
                    r'\(([^)]+)\)',
                    key_lower,
                )

                species = (
                    species_match.group(1)
                    if species_match
                    else 'default'
                )

                try:
                    settings.starting_magnetization[
                        species
                    ] = _parse_qe_float(
                        value
                    )

                except ValueError:
                    logger.warning(
                        "Unable to parse "
                        "starting_magnetization for "
                        "species %r from value %r; "
                        "leaving default",
                        species,
                        value,
                    )

                continue

            value_clean = (
                value
                .strip()
                .strip("'\"")
                .strip()
            )

            if key_lower == 'calculation':
                settings.calculation = (
                    value_clean.lower()
                )

            elif key_lower == 'ecutwfc':
                try:
                    settings.ecutwfc = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse ecutwfc "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'occupations':
                settings.occupations = (
                    value_clean.lower()
                )

            elif key_lower == 'smearing':
                settings.smearing = (
                    value_clean.lower()
                )

            elif key_lower == 'nat':
                try:
                    settings.nat = int(
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse nat "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'ntyp':
                try:
                    settings.ntyp = int(
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse ntyp "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'pseudo_dir':
                settings.pseudo_dir = (
                    value_clean
                )

            elif key_lower == 'degauss':
                try:
                    settings.degauss = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse degauss "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'nspin':
                try:
                    settings.nspin = int(
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse nspin "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'conv_thr':
                try:
                    settings.conv_thr = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse conv_thr "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'tot_charge':
                try:
                    settings.total_charge = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse tot_charge "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'nbnd':
                try:
                    settings.nbands = int(
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse nbnd "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'forc_conv_thr':
                try:
                    settings.forc_conv_thr = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse forc_conv_thr "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'ion_dynamics':
                settings.ion_dynamics = (
                    value_clean.lower()
                )

            elif key_lower == 'trust_radius_max':
                try:
                    settings.trust_radius_max = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse trust_radius_max "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'cell_dynamics':
                settings.cell_dynamics = (
                    value_clean.lower()
                )

            elif key_lower == 'cell_dofree':
                settings.cell_dofree = (
                    value_clean.lower()
                )

            elif key_lower == 'press':
                try:
                    settings.pressure_kbar = (
                        _parse_qe_float(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse press "
                        "from value %r; leaving default",
                        value_clean,
                    )

            elif key_lower == 'nosym':
                try:
                    settings.nosym = (
                        _parse_qe_bool(
                            value_clean
                        )
                    )
                except ValueError:
                    logger.warning(
                        "Unable to parse nosym "
                        "from value %r; leaving default",
                        value_clean,
                    )

    return settings


def determine_system_name(input_path: Path, override: Optional[str]) -> str:
    if override:
        return _sanitize_name(override)
    return _sanitize_name(input_path.stem)


def _sanitize_name(name: str) -> str:
    safe = "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in name.strip())
    return safe or "system"


def build_config_lines(
    name: str,
    geom_filename: str,
    settings: QEInputSettings,
    args: argparse.Namespace,
) -> List[str]:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    outdirname = args.outdirname or f"{name}-results"
    xc = args.xc or "PBE"
    spin_calc = settings.nspin == 2

    workflow_lines = (
        _build_workflow_lines(
            settings
        )
    )

    lines: List[str] = [
        f"# Auto-generated on {timestamp} by qeconverter.py",
        f"Outdirname = '{outdirname}'",
        "",
        "Engine = 'QE'",
        "Mode = 'PW'",
        *workflow_lines,
        "",
    ]

    relaxation_lines = (
        _build_relaxation_lines(
            settings
        )
    )

    if relaxation_lines:
        lines.extend(
            relaxation_lines
        )
        lines.append(
            ''
        )

    if settings.ecutwfc is not None:
        lines.append(f"Cut_off_energy = {settings.ecutwfc * RY_TO_EV:.1f}")
    else:
        lines.append("Cut_off_energy = 340.0")

    lines.extend(
        _build_kpoint_lines(
            settings
        )
    )

    lines.append(f"XC_calc = '{xc}'")

    lines.append(
        _build_occupation_line(
            settings
        )
    )

    if settings.total_charge is not None:
        lines.append(
            "Total_charge = "
            f"{settings.total_charge:.12g}"
        )

    if settings.nbands is not None:
        calculation = (
            settings.calculation
            or 'scf'
        ).strip().lower()

        if calculation == 'bands':
            band_keyword = (
                'Band_num_of_bands'
            )

        elif calculation == 'nscf':
            band_keyword = (
                'DOS_num_of_bands'
            )

        else:
            band_keyword = (
                'Ground_num_of_bands'
            )

        lines.append(
            f"{band_keyword} = "
            f"{settings.nbands}"
        )

    lines.append(f"Spin_calc = {spin_calc}")

    if settings.conv_thr is not None:
        energy_conv = settings.conv_thr * RY_TO_EV
        lines.append(f"Ground_convergence = {{'energy': {energy_conv}}}")

    lines.extend([
        "MPI_cores = 4",
        "Localisation = 'en_UK'",
        "",
        f"# Geometry file to use with dftsolve: {geom_filename}",
    ])

    return [line.rstrip() for line in lines]


def main() -> None:
    args = parse_args()
    if args.version:
        print(f"nanoworks: qeconverter: version: {nanoworks.__version__}")
        return

    if not args.input:
        print("Error: the following arguments are required: --input")
        return

    input_path = args.input.resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Quantum ESPRESSO input not found: {input_path}")

    atoms = read(input_path, format='espresso-in')
    name = determine_system_name(input_path, args.system_name)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    geom_filename = f"{name}.cif"
    geometry_path = output_dir / geom_filename
    write(geometry_path, atoms, format='cif')

    settings = parse_qe_input(input_path)

    input_filename = args.input_filename or f"{name}.py"
    input_path_out = output_dir / input_filename
    config_lines = build_config_lines(
        name=name,
        geom_filename=geom_filename,
        settings=settings,
        args=args,
    )
    input_path_out.write_text("\n".join(config_lines) + "\n")

    print(f"Wrote geometry to {geometry_path}")
    print(
        f"Wrote Nanoworks input to {input_path_out}"
    )
    print(
        "Run: dftsolve "
        f"-i {input_path_out} "
        f"-g {geometry_path}"
    )


if __name__ == "__main__":  # pragma: no cover
    main()
