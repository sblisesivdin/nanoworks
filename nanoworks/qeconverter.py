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

from ase.io import read, write

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

    for line in lines:
        upper = line.upper()

        if upper.startswith('&'):
            continue

        if upper == '/':
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

    lines: List[str] = [
        f"# Auto-generated on {timestamp} by qeconverter.py",
        f"Outdirname = '{outdirname}'",
        "",
        "Engine = 'QE'",
        "Mode = 'PW'",
        "Ground_calc = True",
        f"Geo_optim = {settings.calculation in {'relax', 'vc-relax'}}",
        "Elastic_calc = False",
        "DOS_calc = False",
        "Band_calc = False",
        "Density_calc = False",
        "Optical_calc = False",
        "",
    ]

    if settings.ecutwfc is not None:
        lines.append(f"Cut_off_energy = {settings.ecutwfc * RY_TO_EV:.1f}")
    else:
        lines.append("Cut_off_energy = 340.0")

    if settings.k_mesh:
        mesh = settings.k_mesh + [1, 1, 1]
        mesh = mesh[:3]
        lines.extend([
            f"Ground_kpts_x = {mesh[0]}",
            f"Ground_kpts_y = {mesh[1]}",
            f"Ground_kpts_z = {mesh[2]}",
        ])
        if settings.k_shift:
            gamma = all(shift == 0 for shift in settings.k_shift[:3])
            lines.append(f"Gamma = {gamma}")

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
        lines.append(
            "Ground_num_of_bands = "
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
