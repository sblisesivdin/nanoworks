#!/usr/bin/env python
# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.


'''
mdsolve.py: Quick Geometric Optimization
              using interatomic potentials.
More information: $ mdsolve.py -h
'''

Description = f''' 
 Usage: 
 $ mdsolve.py <args>
 
 -------------------------------------------------------------
   Some potentials
 -------------------------------------------------------------
 | Name                                                  | Information                           | 
 | ----------------------------------------------------- | ------------------------------------- |
 | LJ_ElliottAkerson_2015_Universal__MO_959249795837_003 | General potential for all elements    |
 
 '''

import getopt
import sys
import os
import time
import textwrap
import shutil
import subprocess
import requests
import nanoworks
from argparse import ArgumentParser, HelpFormatter
from pathlib import Path
from numbers import Number
from itertools import product
from ase import *
from ase.data import atomic_numbers, atomic_masses
from ase.io import read, write
from ase.io.cif import write_cif
from ase.spacegroup.symmetrize import check_symmetry
from asap3 import Atoms, units
from asap3.md.langevin import Langevin
from ase.calculators.kim import KIM


# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------

# Simulation parameters
Engine = 'ASAP'
Ensemble = 'NVT'
OpenKIM_potential = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'
Temperature = 1.0 # K
Time_step = 5.0 # fs
Temperature_damp = 200.0 # fs
Pressure = 0.0 # GPa
Pressure_damp = 1000.0 # fs
Random_seed = 12345

# Optional LAMMPS pre-MD minimization
Minimize = False
Minimize_energy_tolerance = 1.0e-10
Minimize_force_tolerance = 1.0e-6 # eV/Angstrom
Minimize_max_iterations = 10000
Minimize_max_evaluations = 100000

# Optional LAMMPS trajectory analysis
MSD_calc = False
MSD_interval = 1
MSD_remove_com = True

RDF_calc = False
RDF_bins = 100
RDF_interval = 10

VACF_calc = False
VACF_interval = 1

# Molecular dynamics loop configuration
MD_cycles = 25
MD_steps_per_cycle = 10

# Short aliases for verbose OpenKIM potential identifiers
POTENTIAL_ALIASES = {
    'lj': 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003',
}


def _ensure_iterable(value):
    """Normalize schedule-like inputs into a list of floats."""
    if isinstance(value, tuple):
        if len(value) == 3 and all(isinstance(i, Number) for i in value):
            start, stop, step = value
            if step == 0:
                raise ValueError('Step cannot be zero when defining a range.')
            seq = []
            current = start
            comparison = (lambda a, b: a <= b + 1e-12) if step > 0 else (lambda a, b: a >= b - 1e-12)
            while comparison(current, stop):
                seq.append(current)
                current += step
            return seq
        return list(value)
    if isinstance(value, list):
        return list(value)
    if isinstance(value, dict):
        start = value.get('start')
        stop = value.get('stop', start)
        count = value.get('count')
        step = value.get('step')
        if start is None:
            raise ValueError('Range dictionaries must define a "start" key.')
        if count:
            if count <= 1:
                return [float(start)]
            delta = (stop - start) / (count - 1)
            return [start + delta * idx for idx in range(count)]
        if step is None:
            return [float(start), float(stop)] if stop is not None else [float(start)]
        if step == 0:
            raise ValueError('Step cannot be zero when defining a range.')
        seq = []
        current = start
        comparison = (lambda a, b: a <= b + 1e-12) if step > 0 else (lambda a, b: a >= b - 1e-12)
        while comparison(current, stop):
            seq.append(current)
            current += step
        return seq
    if isinstance(value, Number):
        return [value]
    return []


def _build_profile(name, base_value, cycle_count, namespace):
    """Return per-cycle parameter profile with graceful fallbacks."""
    keys = [f"{name}_profile", f"{name}_range"]
    sequence = []
    for key in keys:
        if key in namespace:
            sequence = _ensure_iterable(namespace[key])
            break
    if not sequence:
        sequence = _ensure_iterable(namespace.get(name, base_value))
    if not sequence:
        sequence = [base_value]
    if len(sequence) >= cycle_count:
        return sequence[:cycle_count]
    return sequence + [sequence[-1]] * (cycle_count - len(sequence))

def _resolve_engine(namespace):
    """Resolve and validate the molecular dynamics engine."""

    engine = str(namespace.get('Engine', Engine)).strip().upper()

    supported_engines = ('ASAP', 'LAMMPS')

    if engine not in supported_engines:
        raise ValueError(
            f"Unsupported MD engine: {engine}. "
            f"Supported engines: {', '.join(supported_engines)}"
        )

    return engine

def _resolve_ensemble(namespace, engine):
    """Resolve and validate the molecular dynamics ensemble."""

    ensemble = str(
        namespace.get('Ensemble', Ensemble)
    ).strip().upper()

    supported = {
        'ASAP': ('NVT',),
        'LAMMPS': ('NVT', 'NVE', 'NPT'),
    }

    if ensemble not in supported[engine]:
        raise ValueError(
            f'Unsupported MD ensemble {ensemble} for {engine}. '
            f'Supported ensembles: '
            f'{", ".join(supported[engine])}'
        )

    return ensemble

def _validate_ensemble_settings(namespace, ensemble):
    """Validate ensemble-specific molecular dynamics settings."""

    unsupported = ()

    if ensemble == 'NVE':
        unsupported = (
            'Temperature_profile',
            'Temperature_range',
            'Temperature_damp_profile',
            'Temperature_damp_range',
            'Temperature_damp_values',
            'Pressure_profile',
            'Pressure_range',
            'Pressure_values',
            'Pressure_damp_profile',
            'Pressure_damp_range',
            'Pressure_damp_values',
        )

    elif ensemble == 'NVT':
        unsupported = (
            'Pressure_profile',
            'Pressure_range',
            'Pressure_values',
            'Pressure_damp_profile',
            'Pressure_damp_range',
            'Pressure_damp_values',
        )

    present = [
        name for name in unsupported
        if name in namespace
    ]

    if present:
        raise ValueError(
            f'{ensemble} does not use these settings: '
            + ', '.join(present)
        )

def _check_lammps_available():
    """Check whether the system-wide LAMMPS executable is available."""

    executable = shutil.which('lmp')

    if executable is None:
        raise FileNotFoundError(
            "LAMMPS executable 'lmp' was not found in PATH."
        )

    return executable

def _resolve_potential(namespace, alias):
    """Resolve OpenKIM potential either from alias or explicit id."""
    if alias and alias in POTENTIAL_ALIASES:
        return POTENTIAL_ALIASES[alias]
    if alias and alias not in POTENTIAL_ALIASES:
        raise KeyError(f'Unknown potential alias: {alias}')
    potential = namespace.get('OpenKIM_potential', OpenKIM_potential)
    alias_key = namespace.get('OpenKIM_potential_alias')
    if alias_key:
        if alias_key in POTENTIAL_ALIASES:
            potential = POTENTIAL_ALIASES[alias_key]
        else:
            raise KeyError(f'Unknown potential alias: {alias_key}')
    return potential


def _write_energy_csv(path, records):
    if not records:
        return
    header = (
    "Step,Cycle,PotentialEnergyPerAtom(eV),"
    "KineticEnergyPerAtom(eV),TotalEnergyPerAtom(eV),"
    "Temperature(K),TimeStep(fs),TemperatureDamp(fs),"
    "Pressure(GPa),Volume(A^3)"
)
    with open(path, 'w') as fd:
        fd.write(header + "\n")
        for rec in records:
            temperature_damp = rec.get(
                'temperature_damp'
            )

            if temperature_damp is None:
                temperature_damp_text = ''
            else:
                temperature_damp_text = (
                    f'{temperature_damp:.6f}'
                )

            pressure = rec.get('pressure')
            volume = rec.get('volume')

            pressure_text = (
                '' if pressure is None
                else f'{pressure:.6f}'
            )
            volume_text = (
                '' if volume is None
                else f'{volume:.6f}'
            )

            fd.write(
                f"{rec['step']},{rec['cycle']},"
                f"{rec['epot']:.6f},{rec['ekin']:.6f},"
                f"{rec['total']:.6f},"
                f"{rec['temperature']:.6f},"
                f"{rec['timestep']:.6f},"
                f"{temperature_damp_text},"
                f"{pressure_text},{volume_text}\n"
            )


def _get_run_values(name, default, namespace):
    values_key = f"{name}_values"
    if values_key in namespace:
        seq = _ensure_iterable(namespace[values_key])
        if seq:
            return seq
    value = namespace.get(name, default)
    if isinstance(value, Number):
        return [value]
    seq = _ensure_iterable(value)
    return seq if seq else [default]


def _format_suffix(prefix, value):
    if isinstance(value, Number):
        formatted = ('{:.6g}'.format(value)).replace('-', 'm').replace('.', 'p')
    else:
        formatted = str(value)
    return f"{prefix}{formatted}"

def _export_cif(struct_prefix, atoms):
    write_cif(struct_prefix+'-FinalStructure.cif', atoms)

def _write_lammps_data(atoms, struct_prefix):
    """Write the ASE structure as a LAMMPS data file."""

    data_file = struct_prefix + '-LAMMPS.data'

    species = list(dict.fromkeys(atoms.get_chemical_symbols()))

    write(
        data_file,
        atoms,
        format='lammps-data',
        atom_style='atomic',
        specorder=species,
    )

    return data_file, species

def _write_lammps_input(
    struct_prefix,
    data_file,
    species,
    pbc,
    ensemble,
    openkim_potential,
    temperature_profile,
    timestep_profile,
    temperature_damp_profile,
    pressure_profile,
    pressure_damp_profile,
    minimize,
    minimize_energy_tolerance,
    minimize_force_tolerance,
    minimize_max_iterations,
    minimize_max_evaluations,
    msd_calc,
    msd_interval,
    msd_remove_com,
    rdf_calc,
    rdf_bins,
    rdf_interval,
    vacf_calc,
    vacf_interval,
    random_seed,
    md_cycles,
    md_steps_per_cycle,
):
    """Write a LAMMPS input file for molecular dynamics."""

    input_file = struct_prefix + '-LAMMPS.in'
    dump_file = struct_prefix + '-LAMMPS.dump'
    msd_file = struct_prefix + '-LAMMPS-MSD.dat'
    rdf_file = struct_prefix + '-LAMMPS-RDF.dat'
    vacf_file = struct_prefix + '-LAMMPS-VACF.dat'

    species_string = ' '.join(species)
    
    boundary_string = ' '.join(
        'p' if periodic else 'm'
        for periodic in pbc
    )

    mass_lines = []

    for type_id, symbol in enumerate(species, start=1):
        atomic_number = atomic_numbers[symbol]
        atomic_mass = atomic_masses[atomic_number]

        mass_lines.append(
            f'mass {type_id} {atomic_mass:.8f}'
        )

    initial_temperature = float(
        temperature_profile[0]
    )

    if ensemble == 'NPT' and not all(pbc):
        raise ValueError(
            'LAMMPS NPT currently requires periodic boundary '
            'conditions in all three directions.'
        )

    lines = [
        '# LAMMPS input generated by Nanoworks',
        '',
        f'kim init {openkim_potential} metal',
        '',
        'atom_style atomic',
        f'boundary {boundary_string}',
        f'read_data "{data_file}"',
        '',
        *mass_lines,
        '',
        f'kim interactions {species_string}',
        '',
    ]

    if minimize:
        if minimize_energy_tolerance < 0.0:
            raise ValueError(
                'Minimize_energy_tolerance cannot be negative.'
            )
        if minimize_force_tolerance < 0.0:
            raise ValueError(
                'Minimize_force_tolerance cannot be negative.'
            )
        if minimize_max_iterations <= 0:
            raise ValueError(
                'Minimize_max_iterations must be positive.'
            )
        if minimize_max_evaluations <= 0:
            raise ValueError(
                'Minimize_max_evaluations must be positive.'
            )

        lines.extend(
            [
                '# Pre-MD energy minimization',
                'min_style cg',
                (
                    f'minimize '
                    f'{minimize_energy_tolerance:.12g} '
                    f'{minimize_force_tolerance:.12g} '
                    f'{int(minimize_max_iterations)} '
                    f'{int(minimize_max_evaluations)}'
                ),
                '',
            ]
        )

    lines.extend(
        [
        (
            f'velocity all create '
            f'{initial_temperature:.8f} {int(random_seed)} '
            'mom yes rot yes dist gaussian'
        ),
        '',
        ]
    )

    if msd_calc:
        if int(msd_interval) <= 0:
            raise ValueError(
                'MSD_interval must be a positive integer.'
            )

        com_option = 'yes' if msd_remove_com else 'no'

        lines.extend(
            [
                (
                    f'compute nw_msd all msd '
                    f'com {com_option}'
                ),
                (
                    f'fix nw_msd_output all ave/time '
                    f'{int(msd_interval)} 1 '
                    f'{int(msd_interval)} '
                    f'c_nw_msd[1] c_nw_msd[2] '
                    f'c_nw_msd[3] c_nw_msd[4] '
                    f'file "{msd_file}" mode scalar'
                ),
                '',
            ]
        )

    if rdf_calc:
        if int(rdf_bins) <= 0:
            raise ValueError(
                'RDF_bins must be a positive integer.'
            )
        if int(rdf_interval) <= 0:
            raise ValueError(
                'RDF_interval must be a positive integer.'
            )

        lines.extend(
            [
                f'compute nw_rdf all rdf {int(rdf_bins)}',
                (
                    f'fix nw_rdf_output all ave/time '
                    f'{int(rdf_interval)} 1 '
                    f'{int(rdf_interval)} '
                    f'c_nw_rdf[*] '
                    f'file "{rdf_file}" mode vector'
                ),
                '',
            ]
        )

    if vacf_calc:
        if int(vacf_interval) <= 0:
            raise ValueError(
                'VACF_interval must be a positive integer.'
            )

        lines.extend(
            [
                'compute nw_vacf all vacf',
                (
                    f'fix nw_vacf_output all ave/time '
                    f'{int(vacf_interval)} 1 '
                    f'{int(vacf_interval)} '
                    f'c_nw_vacf[1] c_nw_vacf[2] '
                    f'c_nw_vacf[3] c_nw_vacf[4] '
                    f'file "{vacf_file}" mode scalar'
                ),
                '',
            ]
        )

    lines.extend(
        [
        (
            f'dump nw_dump all custom 1 "{dump_file}" '
            'id type x y z vx vy vz'
        ),
        'dump_modify nw_dump sort id',
        '',
        'thermo 1',
        'thermo_style custom step atoms temp press pe ke etotal vol',
        '',
        ]
    )

    if ensemble in ('NVT', 'NVE'):
        lines.extend(
            [
                'fix nw_integrator all nve',
                '',
            ]
        )

    for cycle in range(md_cycles):
        timestep_fs = float(
            timestep_profile[cycle]
        )

        if timestep_fs <= 0.0:
            raise ValueError(
                'Time_step must be greater than zero.'
            )

        # LAMMPS metal units use picoseconds.
        timestep_ps = timestep_fs / 1000.0

        lines.extend(
            [
                f'# MD cycle {cycle + 1}',
                f'timestep {timestep_ps:.10f}',
            ]
        )

        if ensemble == 'NVT':
            temperature = float(
                temperature_profile[cycle]
            )
            temperature_damp_fs = float(
                temperature_damp_profile[cycle]
            )

            if temperature_damp_fs <= 0.0:
                raise ValueError(
                    'Temperature_damp must be greater than zero.'
                )

            temperature_damp_ps = (
                temperature_damp_fs / 1000.0
            )
            cycle_seed = int(random_seed) + cycle

            lines.extend(
                [
                    (
                        f'fix nw_thermostat all langevin '
                        f'{temperature:.8f} {temperature:.8f} '
                        f'{temperature_damp_ps:.10f} '
                        f'{cycle_seed} zero yes'
                    ),
                    f'run {int(md_steps_per_cycle)}',
                    'unfix nw_thermostat',
                    '',
                ]
            )

        elif ensemble == 'NPT':
            temperature = float(
                temperature_profile[cycle]
            )
            temperature_damp_fs = float(
                temperature_damp_profile[cycle]
            )
            pressure_gpa = float(
                pressure_profile[cycle]
            )
            pressure_damp_fs = float(
                pressure_damp_profile[cycle]
            )

            if temperature_damp_fs <= 0.0:
                raise ValueError(
                    'Temperature_damp must be greater than zero.'
                )

            if pressure_damp_fs <= 0.0:
                raise ValueError(
                    'Pressure_damp must be greater than zero.'
                )

            temperature_damp_ps = (
                temperature_damp_fs / 1000.0
            )
            pressure_damp_ps = (
                pressure_damp_fs / 1000.0
            )

            # LAMMPS metal pressure units are bar.
            pressure_bar = pressure_gpa * 10000.0

            lines.extend(
                [
                    (
                        f'fix nw_barostat all npt '
                        f'temp {temperature:.8f} '
                        f'{temperature:.8f} '
                        f'{temperature_damp_ps:.10f} '
                        f'iso {pressure_bar:.8f} '
                        f'{pressure_bar:.8f} '
                        f'{pressure_damp_ps:.10f}'
                    ),
                    f'run {int(md_steps_per_cycle)}',
                    'unfix nw_barostat',
                    '',
                ]
            )

        else:
            lines.extend(
                [
                    f'run {int(md_steps_per_cycle)}',
                    '',
                ]
            )

    if ensemble in ('NVT', 'NVE'):
        lines.append('unfix nw_integrator')

    if msd_calc:
        lines.extend(
            [
                'unfix nw_msd_output',
                'uncompute nw_msd',
            ]
        )

    if rdf_calc:
        lines.extend(
            [
                'unfix nw_rdf_output',
                'uncompute nw_rdf',
            ]
        )

    if vacf_calc:
        lines.extend(
            [
                'unfix nw_vacf_output',
                'uncompute nw_vacf',
            ]
        )

    lines.extend(
        [
            'undump nw_dump',
            '',
        ]
    )

    with open(input_file, 'w') as fd:
        fd.write('\n'.join(lines))

    return input_file

def _execute_lammps(input_file, struct_prefix):
    """Execute the system-wide LAMMPS binary."""

    log_file = struct_prefix + '-LAMMPS.log'

    command = [
        'lmp',
        '-in',
        input_file,
        '-log',
        log_file,
    ]

    print('')
    print('Starting LAMMPS calculation...')

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
    )

    if result.stdout:
        print(result.stdout, end='')

    if result.returncode != 0:
        if result.stderr:
            print(result.stderr, file=sys.stderr, end='')

        raise RuntimeError(
            f'LAMMPS calculation failed with return code '
            f'{result.returncode}. See {log_file}'
        )

    print('LAMMPS calculation finished.')

    return log_file

def _update_atoms_from_lammps_dump(
    atoms,
    struct_prefix,
    species,
):
    """Update ASE atoms from the final LAMMPS dump frame."""

    dump_file = struct_prefix + '-LAMMPS.dump'

    final_atoms = read(
        dump_file,
        index=-1,
        format='lammps-dump-text',
        specorder=species,
    )

    if len(final_atoms) != len(atoms):
        raise ValueError(
            'LAMMPS final structure contains a different '
            'number of atoms.'
        )

    atoms.set_cell(
        final_atoms.get_cell(),
        scale_atoms=False,
    )

    atoms.set_positions(
        final_atoms.get_positions()
    )

    atoms.set_pbc(
        final_atoms.get_pbc()
    )

    velocities = final_atoms.get_velocities()

    if velocities is not None:
        atoms.set_velocities(velocities)

def _write_lammps_trajectory(
    struct_prefix,
    species,
):
    """Convert the LAMMPS dump trajectory to ASE trajectory format."""

    dump_file = struct_prefix + '-LAMMPS.dump'
    trajectory_file = struct_prefix + '-Results.traj'

    frames = read(
        dump_file,
        index=':',
        format='lammps-dump-text',
        specorder=species,
    )

    write(
        trajectory_file,
        frames,
        format='traj',
    )

    return trajectory_file

def _write_lammps_msd_csv(struct_prefix):
    """Convert LAMMPS MSD output to a Nanoworks CSV file."""

    source_file = struct_prefix + '-LAMMPS-MSD.dat'
    csv_file = struct_prefix + '-MSD.csv'

    if not os.path.isfile(source_file):
        raise FileNotFoundError(
            f'LAMMPS MSD output was not found: {source_file}'
        )

    records = []

    with open(source_file, 'r') as fd:
        for line in fd:
            stripped = line.strip()

            if not stripped or stripped.startswith('#'):
                continue

            fields = stripped.split()

            if len(fields) < 5:
                continue

            try:
                step = int(float(fields[0]))
                msd_x = float(fields[1])
                msd_y = float(fields[2])
                msd_z = float(fields[3])
                msd_total = float(fields[4])
            except ValueError:
                continue

            records.append(
                (
                    step,
                    msd_x,
                    msd_y,
                    msd_z,
                    msd_total,
                )
            )

    with open(csv_file, 'w') as fd:
        fd.write(
            'Step,MSD_X(A^2),MSD_Y(A^2),'
            'MSD_Z(A^2),MSD_Total(A^2)\n'
        )

        for record in records:
            fd.write(
                f'{record[0]},'
                f'{record[1]:.10g},'
                f'{record[2]:.10g},'
                f'{record[3]:.10g},'
                f'{record[4]:.10g}\n'
            )

    return csv_file

def _write_lammps_rdf_csv(struct_prefix):
    """Convert LAMMPS RDF output to a Nanoworks CSV file."""

    source_file = struct_prefix + '-LAMMPS-RDF.dat'
    csv_file = struct_prefix + '-RDF.csv'

    if not os.path.isfile(source_file):
        raise FileNotFoundError(
            f'LAMMPS RDF output was not found: {source_file}'
        )

    records = []

    with open(source_file, 'r') as fd:
        lines = [
            line.strip()
            for line in fd
            if line.strip()
            and not line.lstrip().startswith('#')
        ]

    index = 0

    while index < len(lines):
        header = lines[index].split()

        if len(header) != 2:
            index += 1
            continue

        try:
            step = int(float(header[0]))
            row_count = int(float(header[1]))
        except ValueError:
            index += 1
            continue

        index += 1

        for _ in range(row_count):
            if index >= len(lines):
                break

            fields = lines[index].split()
            index += 1

            if len(fields) < 4:
                continue

            try:
                bin_index = int(float(fields[0]))
                radius = float(fields[1])
                rdf = float(fields[2])
                coordination = float(fields[3])
            except ValueError:
                continue

            records.append(
                (
                    step,
                    bin_index,
                    radius,
                    rdf,
                    coordination,
                )
            )

    with open(csv_file, 'w') as fd:
        fd.write(
            'Step,Bin,R(A),g(r),CoordinationNumber\n'
        )

        for record in records:
            fd.write(
                f'{record[0]},'
                f'{record[1]},'
                f'{record[2]:.10g},'
                f'{record[3]:.10g},'
                f'{record[4]:.10g}\n'
            )

    return csv_file

def _write_lammps_vacf_csv(struct_prefix):
    """Convert LAMMPS VACF output to a Nanoworks CSV file."""

    source_file = struct_prefix + '-LAMMPS-VACF.dat'
    csv_file = struct_prefix + '-VACF.csv'

    if not os.path.isfile(source_file):
        raise FileNotFoundError(
            f'LAMMPS VACF output was not found: {source_file}'
        )

    records = []

    with open(source_file, 'r') as fd:
        for line in fd:
            stripped = line.strip()

            if not stripped or stripped.startswith('#'):
                continue

            fields = stripped.split()

            if len(fields) < 5:
                continue

            try:
                step = int(float(fields[0]))
                vacf_x = float(fields[1])
                vacf_y = float(fields[2])
                vacf_z = float(fields[3])
                vacf_total = float(fields[4])
            except ValueError:
                continue

            records.append(
                (
                    step,
                    vacf_x,
                    vacf_y,
                    vacf_z,
                    vacf_total,
                )
            )

    with open(csv_file, 'w') as fd:
        fd.write(
            'Step,VACF_X(A^2/ps^2),'
            'VACF_Y(A^2/ps^2),'
            'VACF_Z(A^2/ps^2),'
            'VACF_Total(A^2/ps^2)\n'
        )

        for record in records:
            fd.write(
                f'{record[0]},'
                f'{record[1]:.10g},'
                f'{record[2]:.10g},'
                f'{record[3]:.10g},'
                f'{record[4]:.10g}\n'
            )

    return csv_file

def _parse_lammps_thermo(
    log_file,
    ensemble,
    temperature_profile,
    timestep_profile,
    temperature_damp_profile,
    pressure_profile,
    pressure_damp_profile,
    md_steps_per_cycle,
):
    """Parse LAMMPS thermo output into Nanoworks energy records."""

    records = []
    seen_steps = set()
    in_thermo = False

    with open(log_file, 'r') as fd:
        for line in fd:
            stripped = line.strip()

            if stripped.startswith(
                'Step'
            ) and 'PotEng' in stripped and 'KinEng' in stripped:
                in_thermo = True
                continue

            if not in_thermo:
                continue

            if (
                stripped.startswith('Loop time')
                or stripped.startswith('ERROR')
                or not stripped
            ):
                in_thermo = False
                continue

            fields = stripped.split()

            if len(fields) < 8:
                continue

            try:
                step = int(fields[0])
                atoms = int(fields[1])
                instantaneous_temperature = float(fields[2])
                instantaneous_pressure_bar = float(fields[3])
                potential_energy = float(fields[4])
                kinetic_energy = float(fields[5])
                total_energy = float(fields[6])
                volume = float(fields[7])
            except ValueError:
                continue

            if step == 0:
                continue

            if step % md_steps_per_cycle != 0:
                continue

            if step in seen_steps:
                continue

            seen_steps.add(step)

            cycle = step // md_steps_per_cycle
            profile_index = cycle - 1

            if profile_index >= len(temperature_profile):
                continue

            if ensemble == 'NVE':
                recorded_temperature = (
                    instantaneous_temperature
                )
                recorded_temperature_damp = None
            else:
                recorded_temperature = float(
                    temperature_profile[profile_index]
                )
                recorded_temperature_damp = float(
                    temperature_damp_profile[
                        profile_index
                    ]
                )

            recorded_pressure = (
                instantaneous_pressure_bar / 10000.0
            )

            records.append(
                {
                    'cycle': cycle,
                    'step': step,
                    'epot': potential_energy / atoms,
                    'ekin': kinetic_energy / atoms,
                    'total': total_energy / atoms,
                    'temperature': recorded_temperature,
                    'timestep': float(
                        timestep_profile[profile_index]
                    ),
                    'temperature_damp': recorded_temperature_damp,
                    'pressure': recorded_pressure,
                    'volume': volume,
                }
            )

    return records

def _run_asap_langevin(
    atoms,
    struct_prefix,
    openkim_potential,
    temperature_profile,
    timestep_profile,
    temperature_damp_profile,
    md_cycles,
    md_steps_per_cycle,
):
    """Run the ASAP3 Langevin MD workflow."""

    atoms.calc = KIM(
        openkim_potential,
        options={"ase_neigh": False},
    )

    initial_temperature = float(temperature_profile[0])
    initial_timestep = float(timestep_profile[0])
    initial_temperature_damp = float(temperature_damp_profile[0])

    if initial_temperature_damp <= 0.0:
        raise ValueError(
            'Temperature_damp must be greater than zero.'
        )

    initial_friction = 1.0 / (
        initial_temperature_damp * units.fs
    )

    dyn = Langevin(
        atoms,
        timestep=initial_timestep * units.fs,
        trajectory=struct_prefix + '-Results.traj',
        logfile=struct_prefix + '-Log.txt',
        temperature_K=initial_temperature,
        friction=initial_friction,
    )

    print("")
    print("Energy per atom (cycle averages):")
    print(
        "  %6s %15s %15s %15s %12s %12s %15s"
        % (
            "Cycle",
            "Pot. energy",
            "Kin. energy",
            "Total energy",
            "Temp (K)",
            "dt (fs)",
            "T-damp (fs)",
        )
    )

    energy_records = []
    total_steps = 0

    for cycle in range(md_cycles):
        current_temperature = float(
            temperature_profile[cycle]
        )
        current_timestep = float(
            timestep_profile[cycle]
        )
        current_temperature_damp = float(
            temperature_damp_profile[cycle]
        )

        if current_temperature_damp <= 0.0:
            raise ValueError(
                'Temperature_damp must be greater than zero.'
            )

        current_friction = 1.0 / (
            current_temperature_damp * units.fs
        )

        if cycle > 0:
            if abs(
                current_temperature - initial_temperature
            ) > 1e-9:
                dyn.set_temperature(
                    temperature_K=current_temperature
                )

            if abs(
                current_timestep - initial_timestep
            ) > 1e-12:
                dyn.set_timestep(
                    current_timestep * units.fs
                )

            if abs(
                current_temperature_damp
                - initial_temperature_damp
            ) > 1e-12:
                dyn.set_friction(
                    current_friction
                )

        dyn.run(md_steps_per_cycle)

        total_steps += md_steps_per_cycle

        epot = atoms.get_potential_energy() / len(atoms)
        ekin = atoms.get_kinetic_energy() / len(atoms)
        total_energy = epot + ekin

        print(
            "%7d %15.5f %15.5f %15.5f %12.2f %12.3f %15.3f"
            % (
                cycle + 1,
                epot,
                ekin,
                total_energy,
                current_temperature,
                current_timestep,
                current_temperature_damp,
            )
        )

        energy_records.append(
            {
                'cycle': cycle + 1,
                'step': total_steps,
                'epot': epot,
                'ekin': ekin,
                'total': total_energy,
                'temperature': current_temperature,
                'timestep': current_timestep,
                'temperature_damp': current_temperature_damp,
            }
        )

        initial_temperature = current_temperature
        initial_timestep = current_timestep
        initial_temperature_damp = current_temperature_damp

    return energy_records

def _run_md_engine(
    engine,
    ensemble,
    atoms,
    struct_prefix,
    openkim_potential,
    temperature_profile,
    timestep_profile,
    temperature_damp_profile,
    pressure_profile,
    pressure_damp_profile,
    minimize,
    minimize_energy_tolerance,
    minimize_force_tolerance,
    minimize_max_iterations,
    minimize_max_evaluations,
    msd_calc,
    msd_interval,
    msd_remove_com,
    rdf_calc,
    rdf_bins,
    rdf_interval,
    vacf_calc,
    vacf_interval,
    random_seed,
    md_cycles,
    md_steps_per_cycle,
):
    """Dispatch the molecular dynamics run to the selected engine."""

    if engine == 'ASAP':
        return _run_asap_langevin(
            atoms=atoms,
            struct_prefix=struct_prefix,
            openkim_potential=openkim_potential,
            temperature_profile=temperature_profile,
            timestep_profile=timestep_profile,
            temperature_damp_profile=temperature_damp_profile,
            md_cycles=md_cycles,
            md_steps_per_cycle=md_steps_per_cycle,
        )

    if engine == 'LAMMPS':
        data_file, species = _write_lammps_data(
            atoms=atoms,
            struct_prefix=struct_prefix,
        )

        input_file = _write_lammps_input(
            struct_prefix=struct_prefix,
            data_file=data_file,
            species=species,
            pbc=atoms.get_pbc(),
            ensemble=ensemble,
            openkim_potential=openkim_potential,
            temperature_profile=temperature_profile,
            timestep_profile=timestep_profile,
            temperature_damp_profile=temperature_damp_profile,
            pressure_profile=pressure_profile,
            pressure_damp_profile=pressure_damp_profile,
            minimize=minimize,
            minimize_energy_tolerance=minimize_energy_tolerance,
            minimize_force_tolerance=minimize_force_tolerance,
            minimize_max_iterations=minimize_max_iterations,
            minimize_max_evaluations=minimize_max_evaluations,
            msd_calc=msd_calc,
            msd_interval=msd_interval,
            msd_remove_com=msd_remove_com,
            rdf_calc=rdf_calc,
            rdf_bins=rdf_bins,
            rdf_interval=rdf_interval,
            vacf_calc=vacf_calc,
            vacf_interval=vacf_interval,
            random_seed=random_seed,
            md_cycles=md_cycles,
            md_steps_per_cycle=md_steps_per_cycle,
        )

        print(f'LAMMPS data file written: {data_file}')
        print(f'LAMMPS input file written: {input_file}')

        log_file = _execute_lammps(
            input_file=input_file,
            struct_prefix=struct_prefix,
        )

        print(f'LAMMPS log file written: {log_file}')

        _update_atoms_from_lammps_dump(
            atoms=atoms,
            struct_prefix=struct_prefix,
            species=species,
        )

        trajectory_file = _write_lammps_trajectory(
            struct_prefix=struct_prefix,
            species=species,
        )

        print(
            f'LAMMPS trajectory file written: '
            f'{trajectory_file}'
        )

        if msd_calc:
            msd_csv = _write_lammps_msd_csv(
                struct_prefix=struct_prefix,
            )
            print(
                f'LAMMPS MSD file written: {msd_csv}'
            )

        if rdf_calc:
            rdf_csv = _write_lammps_rdf_csv(
                struct_prefix=struct_prefix,
            )
            print(
                f'LAMMPS RDF file written: {rdf_csv}'
            )

        if vacf_calc:
            vacf_csv = _write_lammps_vacf_csv(
                struct_prefix=struct_prefix,
            )
            print(
                f'LAMMPS VACF file written: {vacf_csv}'
            )

        energy_records = _parse_lammps_thermo(
            log_file=log_file,
            ensemble=ensemble,
            temperature_profile=temperature_profile,
            timestep_profile=timestep_profile,
            temperature_damp_profile=temperature_damp_profile,
            pressure_profile=pressure_profile,
            pressure_damp_profile=pressure_damp_profile,
            md_steps_per_cycle=md_steps_per_cycle,
        )

        return energy_records

    raise ValueError(f'Unsupported MD engine: {engine}')
    
Scaled = False # Scaled or Cartesian coordinates
Manual_PBC = False # If you need manual constraint axis

# If Manual_PBC is true then change following:
PBC_constraints = [True, True, False]

# If you do not want to use a CIF file for geometry, please provide
# ASE Atoms object information below. You can use ciftoase.py to
# make your own ASE Atoms object from a CIF file.
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

bulk_configuration = Atoms(
    [
        Atom('Ge', (1.222474e-31, 4.094533468076675e-32, 5.02)),
        Atom('Ge', (-1.9999999993913775e-06, 2.3094022314590417, 4.98)),
    ],
    cell=[
        (4.0, 0.0, 0.0),
        (-1.9999999999999991, 3.464101615137755, 0.0),
        (0.0, 0.0, 20.0)
    ],
    pbc=True,
)
# -------------------------------------------------------------
# ///////   YOU DO NOT NEED TO CHANGE ANYTHING BELOW    \\\\\\\
# -------------------------------------------------------------

def main():
    # Start time
    time0 = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
    # To print Description variable with argparse
    class RawFormatter(HelpFormatter):
        def _fill_text(self, text, width, indent):
            return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

    # Arguments parsing
    parser = ArgumentParser(prog ='mdsolve.py', description=Description, formatter_class=RawFormatter)


    parser.add_argument("-i", "--input", dest = "inputfile", help="Use input file for calculation variables (also you can insert geometry)")
    parser.add_argument("-g", "--geometry",dest ="geometryfile", help="Use CIF file for geometry")
    parser.add_argument("-v", "--version", dest="version", action='store_true')

    args = parser.parse_args()

    if args.version:
        print(f"nanoworks: mdsolve: version: {nanoworks.__version__}")
        sys.exit(1)

    if not args.inputfile or not args.geometryfile:
        print('ERROR: Please provide both input (-i) and geometry (-g) files.')
        sys.exit(1)

    Outdirname = ''
    configpath = None
    config_dir = Path.cwd()
    input_stem = None
    inFile = None

    try:
        configpath = Path(args.inputfile).expanduser()
        if not configpath.is_absolute():
            configpath = (Path.cwd() / configpath).resolve()
        config_dir = configpath.parent
        input_stem = configpath.stem
        if str(config_dir) not in sys.path:
            sys.path.append(str(config_dir))
        conf = __import__(configpath.stem, globals(), locals(), ['*'])
        for k in dir(conf):
            globals()[k] = getattr(conf, k)

        geometry_path = Path(args.geometryfile).expanduser()
        if not geometry_path.is_absolute():
            geometry_path = (config_dir / geometry_path).resolve()
        inFile = str(geometry_path)

    except getopt.error as err:
        print(str(err))

    namespace = globals()

    try:
        resolved_engine = _resolve_engine(namespace)
    except ValueError as exc:
        print(str(exc))
        sys.exit(1)

    Engine = resolved_engine
    namespace['Engine'] = Engine

    try:
        resolved_ensemble = _resolve_ensemble(
            namespace,
            Engine,
        )
    except ValueError as exc:
        print(str(exc))
        sys.exit(1)

    Ensemble = resolved_ensemble
    namespace['Ensemble'] = Ensemble

    try:
        _validate_ensemble_settings(
            namespace,
            Ensemble,
        )
    except ValueError as exc:
        print(str(exc))
        sys.exit(1)

    lammps_only_features = []

    if bool(namespace.get('Minimize', Minimize)):
        lammps_only_features.append('Minimize')

    if bool(namespace.get('MSD_calc', MSD_calc)):
        lammps_only_features.append('MSD_calc')

    if bool(namespace.get('RDF_calc', RDF_calc)):
        lammps_only_features.append('RDF_calc')

    if bool(namespace.get('VACF_calc', VACF_calc)):
        lammps_only_features.append('VACF_calc')

    if Engine != 'LAMMPS' and lammps_only_features:
        print(
            'These settings are currently supported only '
            'by LAMMPS: '
            + ', '.join(lammps_only_features)
        )
        sys.exit(1)

    if Engine == 'LAMMPS':
        try:
            _check_lammps_available()
        except FileNotFoundError as exc:
            print(str(exc))
            sys.exit(1)

    try:
        resolved_potential = _resolve_potential(namespace, None)
    except KeyError as exc:
        print(str(exc))
        sys.exit(1)

    OpenKIM_potential = resolved_potential
    namespace['OpenKIM_potential'] = OpenKIM_potential

    md_cycles = int(namespace.get('MD_cycles', MD_cycles))
    if md_cycles <= 0:
        print('MD_cycles must be a positive integer.')
        sys.exit(1)
    md_steps_per_cycle = int(namespace.get('MD_steps_per_cycle', MD_steps_per_cycle))
    if md_steps_per_cycle <= 0:
        print('MD_steps_per_cycle must be a positive integer.')
        sys.exit(1)
    namespace['MD_cycles'] = md_cycles
    namespace['MD_steps_per_cycle'] = md_steps_per_cycle

    # Geometry input is mandatory
    if inFile is None:
        print('ERROR: Geometry file could not be resolved.')
        sys.exit(1)

    struct_name = Path(inFile).stem
    bulk_configuration = read(inFile, index='-1')
    print("Number of atoms imported from CIF file:"+str(bulk_configuration.get_global_number_of_atoms()))
    try:
        symmetry = check_symmetry(
            bulk_configuration,
            symprec=1e-2,
            verbose=False,
        )

        print(
            "Spacegroup of CIF file (ASE):",
            f"{symmetry.international} "
            f"(No. {symmetry.number})",
        )

    except Exception as exc:
        print(
            f"Could not determine spacegroup: {exc}"
        )

    base_directory = config_dir if configpath else Path.cwd()

    if Outdirname != '':
        structpath = Path(Outdirname)
        if not structpath.is_absolute():
            structpath = (base_directory / structpath).resolve()
    else:
        # Use the geometry file's directory and stem for the output folder name
        if inFile is not None:
            structdir = Path(inFile).resolve().parent
            structpath = (structdir / struct_name).resolve() # Now uses struct_name directly
        else:
            # Fallback if inFile is not available (though it should be mandatory)
            structpath = (base_directory / "results").resolve()

    if not os.path.isdir(structpath):
        os.makedirs(structpath, exist_ok=True)

    # Always use the geometry file's stem for the base name of the files
    # This ensures consistency: output files inside the folder named after the geometry file,
    # also start with the geometry file's stem.
    struct_base = os.path.join(str(structpath), struct_name)

    initial_structure = bulk_configuration.copy()
    temperature_options = [
        float(v)
        for v in _get_run_values(
            'Temperature',
            Temperature,
            namespace,
        )
    ]

    timestep_options = [
        float(v)
        for v in _get_run_values(
            'Time_step',
            Time_step,
            namespace,
        )
    ]

    if Ensemble in ('NVT', 'NPT'):
        temperature_damp_options = [
            float(v)
            for v in _get_run_values(
                'Temperature_damp',
                Temperature_damp,
                namespace,
            )
        ]
    else:
        temperature_damp_options = [None]

    if Ensemble == 'NPT':
        pressure_options = [
            float(v)
            for v in _get_run_values(
                'Pressure',
                Pressure,
                namespace,
            )
        ]
        pressure_damp_options = [
            float(v)
            for v in _get_run_values(
                'Pressure_damp',
                Pressure_damp,
                namespace,
            )
        ]
    else:
        pressure_options = [None]
        pressure_damp_options = [None]

    combinations = list(
        product(
            temperature_options,
            timestep_options,
            temperature_damp_options,
            pressure_options,
            pressure_damp_options,
        )
    )

    varying_lengths = {
        'T': len(set(temperature_options)),
        'dt': len(set(timestep_options)),
        'Tdamp': len(set(temperature_damp_options)),
        'P': len(set(pressure_options)),
        'Pdamp': len(set(pressure_damp_options)),
    }

    original_temperature = namespace.get(
        'Temperature',
        Temperature,
    )

    original_timestep = namespace.get(
        'Time_step',
        Time_step,
    )

    original_temperature_damp = namespace.get(
        'Temperature_damp',
        Temperature_damp,
    )

    original_pressure = namespace.get(
        'Pressure',
        Pressure,
    )

    original_pressure_damp = namespace.get(
        'Pressure_damp',
        Pressure_damp,
    )

    for combo_index, (
        temperature_value,
        timestep_value,
        temperature_damp_value,
        pressure_value,
        pressure_damp_value,
    ) in enumerate(combinations, 1):
        asestruct = initial_structure.copy()
        if Manual_PBC:
            asestruct.set_pbc(PBC_constraints)

        suffix_parts = []
        if len(combinations) > 1:
            if varying_lengths['T'] > 1:
                suffix_parts.append(_format_suffix('T', temperature_value))
            if varying_lengths['dt'] > 1:
                suffix_parts.append(_format_suffix('dt', timestep_value))
            if varying_lengths['Tdamp'] > 1:
                suffix_parts.append(
                    _format_suffix(
                        'Tdamp',
                        temperature_damp_value,
                    )
                )
            if varying_lengths['P'] > 1:
                suffix_parts.append(
                    _format_suffix(
                        'P',
                        pressure_value,
                    )
                )
            if varying_lengths['Pdamp'] > 1:
                suffix_parts.append(
                    _format_suffix(
                        'Pdamp',
                        pressure_damp_value,
                    )
                )
        run_struct = struct_base if not suffix_parts else struct_base + '_' + '_'.join(suffix_parts)
        struct_prefix = run_struct

        namespace['Temperature'] = temperature_value
        namespace['Time_step'] = timestep_value
        if Ensemble in ('NVT', 'NPT'):
            namespace['Temperature_damp'] = (
                temperature_damp_value
            )
        if Ensemble == 'NPT':
            namespace['Pressure'] = pressure_value
            namespace['Pressure_damp'] = (
                pressure_damp_value
            )

        temperature_profile = _build_profile(
            'Temperature',
            temperature_value,
            MD_cycles,
            namespace,
        )

        timestep_profile = _build_profile(
            'Time_step',
            timestep_value,
            MD_cycles,
            namespace,
        )

        if Ensemble in ('NVT', 'NPT'):
            temperature_damp_profile = _build_profile(
                'Temperature_damp',
                temperature_damp_value,
                MD_cycles,
                namespace,
            )
        else:
            temperature_damp_profile = None

        if Ensemble == 'NPT':
            pressure_profile = _build_profile(
                'Pressure',
                pressure_value,
                MD_cycles,
                namespace,
            )
            pressure_damp_profile = _build_profile(
                'Pressure_damp',
                pressure_damp_value,
                MD_cycles,
                namespace,
            )
        else:
            pressure_profile = None
            pressure_damp_profile = None

        if len(combinations) > 1:
            print("")
            message = (
                f"Run {combo_index}/{len(combinations)}: "
                f"T={temperature_value} K, "
                f"dt={timestep_value} fs"
            )
            if Ensemble in ('NVT', 'NPT'):
                message += (
                    f", T-damp={temperature_damp_value} fs"
                )
            if Ensemble == 'NPT':
                message += (
                    f", P={pressure_value} GPa, "
                    f"P-damp={pressure_damp_value} fs"
                )
            print(message)

        energy_records = _run_md_engine(
            engine=Engine,
            ensemble=Ensemble,
            atoms=asestruct,
            struct_prefix=struct_prefix,
            openkim_potential=OpenKIM_potential,
            temperature_profile=temperature_profile,
            timestep_profile=timestep_profile,
            temperature_damp_profile=temperature_damp_profile,
            pressure_profile=pressure_profile,
            pressure_damp_profile=pressure_damp_profile,
            minimize=bool(
                namespace.get('Minimize', Minimize)
            ),
            minimize_energy_tolerance=float(
                namespace.get(
                    'Minimize_energy_tolerance',
                    Minimize_energy_tolerance,
                )
            ),
            minimize_force_tolerance=float(
                namespace.get(
                    'Minimize_force_tolerance',
                    Minimize_force_tolerance,
                )
            ),
            minimize_max_iterations=int(
                namespace.get(
                    'Minimize_max_iterations',
                    Minimize_max_iterations,
                )
            ),
            minimize_max_evaluations=int(
                namespace.get(
                    'Minimize_max_evaluations',
                    Minimize_max_evaluations,
                )
            ),
            msd_calc=bool(
                namespace.get('MSD_calc', MSD_calc)
            ),
            msd_interval=int(
                namespace.get(
                    'MSD_interval',
                    MSD_interval,
                )
            ),
            msd_remove_com=bool(
                namespace.get(
                    'MSD_remove_com',
                    MSD_remove_com,
                )
            ),
            rdf_calc=bool(
                namespace.get('RDF_calc', RDF_calc)
            ),
            rdf_bins=int(
                namespace.get(
                    'RDF_bins',
                    RDF_bins,
                )
            ),
            rdf_interval=int(
                namespace.get(
                    'RDF_interval',
                    RDF_interval,
                )
            ),
            vacf_calc=bool(
                namespace.get('VACF_calc', VACF_calc)
            ),
            vacf_interval=int(
                namespace.get(
                    'VACF_interval',
                    VACF_interval,
                )
            ),
            md_cycles=MD_cycles,
            random_seed=Random_seed,
            md_steps_per_cycle=MD_steps_per_cycle,
        )

        # PRINT TO FILE PART -----------------------------------
        _write_energy_csv(struct_prefix+'-Energy.csv', energy_records)

        with open(struct_prefix+'Results.py', 'w') as f:
            f.write("bulk_configuration = Atoms(\n")
            f.write("    [\n")
            if Scaled == True:
                positions = asestruct.get_scaled_positions()
            else:
                positions = asestruct.get_positions()
            
            for symbol, position in zip(
                asestruct.get_chemical_symbols(),
                positions,
            ):
                f.write(
                    "        Atom("
                    f"'{symbol}', "
                    f"({position[0]}, "
                    f"{position[1]}, "
                    f"{position[2]})),\n"
                )

            f.write("    ],\n")

            cell = asestruct.get_cell()
            f.write("    cell=[\n")
            for vector in cell:
                f.write(
                    "        ("
                    f"{vector[0]}, "
                    f"{vector[1]}, "
                    f"{vector[2]}),\n"
                )
            f.write("    ],\n")

            pbc = asestruct.get_pbc()
            f.write(
                "    pbc=["
                f"{bool(pbc[0])},"
                f"{bool(pbc[1])},"
                f"{bool(pbc[2])}],\n"
            )
            f.write(")\n")

        _export_cif(struct_prefix, asestruct)

    namespace['Temperature'] = original_temperature
    namespace['Time_step'] = original_timestep
    namespace['Temperature_damp'] = original_temperature_damp
    namespace['Pressure'] = original_pressure
    namespace['Pressure_damp'] = original_pressure_damp

if __name__ == "__main__":
    main()
