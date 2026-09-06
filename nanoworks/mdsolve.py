#!/usr/bin/env python

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
from ase.spacegroup import get_spacegroup
from asap3 import Atoms, units
from asap3.md.langevin import Langevin
from ase.calculators.kim import KIM


# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------

# Simulation parameters
Engine = 'ASAP'
OpenKIM_potential = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'
Temperature = 1.0 # K
Time_step = 5.0 # fs
Temperature_damp = 200.0 # fs
Random_seed = 12345

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
    "Temperature(K),TimeStep(fs),TemperatureDamp(fs)"
)
    with open(path, 'w') as fd:
        fd.write(header + "\n")
        for rec in records:
            fd.write(
                f"{rec['step']},{rec['cycle']},"
                f"{rec['epot']:.6f},{rec['ekin']:.6f},{rec['total']:.6f},"
                f"{rec['temperature']:.6f},{rec['timestep']:.6f},"
                f"{rec['temperature_damp']:.6f}\n"
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

def _print_attention_message():
    print("    )")
    print("ATTENTION: If you have double number of atoms, it may be caused by ")
    print("           repeating ASE bug https://gitlab.com/ase/ase/-/issues/169 ")
    print("           Please assign Solve_double_element_problem variable as True in this script if necessary.")

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
    openkim_potential,
):
    """Write a minimal LAMMPS input file for OpenKIM initialization."""

    input_file = struct_prefix + '-LAMMPS.in'

    species_string = ' '.join(species)

    mass_lines = []

    for type_id, symbol in enumerate(species, start=1):
        atomic_number = atomic_numbers[symbol]
        atomic_mass = atomic_masses[atomic_number]

        mass_lines.append(
            f'mass {type_id} {atomic_mass:.8f}'
        )

    lines = [
        '# LAMMPS input generated by Nanoworks',
        '',
        f'kim init {openkim_potential} metal',
        '',
        'atom_style atomic',
        f'read_data "{data_file}"',
        '',
        *mass_lines,
        '',
        f'kim interactions {species_string}',
        '',
        'thermo 1',
        'thermo_style custom step atoms temp pe ke etotal',
        '',
        'run 0',
        '',
    ]

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

    atoms.set_calculator(
        KIM(openkim_potential, options={"ase_neigh": False})
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
    atoms,
    struct_prefix,
    openkim_potential,
    temperature_profile,
    timestep_profile,
    temperature_damp_profile,
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
            openkim_potential=openkim_potential,
        )

        print(f'LAMMPS data file written: {data_file}')
        print(f'LAMMPS input file written: {input_file}')

        log_file = _execute_lammps(
            input_file=input_file,
            struct_prefix=struct_prefix,
        )

        print(f'LAMMPS log file written: {log_file}')

        return []

    raise ValueError(f'Unsupported MD engine: {engine}')
    
Scaled = False # Scaled or Cartesian coordinates
Manual_PBC = False # If you need manual constraint axis

# If Manual_PBC is true then change following:
PBC_constraints = [True, True, False]

# If you have double number of elements in your final file
Solve_double_element_problem = True

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
        spacegroup = get_spacegroup(bulk_configuration, symprec=1e-2)
        print("Spacegroup of CIF file (ASE):", f"{spacegroup.symbol} (No. {spacegroup.no})")
    except Exception as exc:
        print(f"Could not determine spacegroup: {exc}")

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

    temperature_damp_options = [
        float(v)
        for v in _get_run_values(
            'Temperature_damp',
            Temperature_damp,
            namespace,
        )
    ]

    combinations = list(
        product(
            temperature_options,
            timestep_options,
            temperature_damp_options,
        )
    )

    varying_lengths = {
        'T': len(set(temperature_options)),
        'dt': len(set(timestep_options)),
        'Tdamp': len(set(temperature_damp_options)),
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

    for combo_index, (
        temperature_value,
        timestep_value,
        temperature_damp_value,
    ) in enumerate(combinations, 1):
        asestruct = initial_structure.copy()

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
        run_struct = struct_base if not suffix_parts else struct_base + '_' + '_'.join(suffix_parts)
        struct_prefix = run_struct

        namespace['Temperature'] = temperature_value
        namespace['Time_step'] = timestep_value
        namespace['Temperature_damp'] = temperature_damp_value

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

        temperature_damp_profile = _build_profile(
            'Temperature_damp',
            temperature_damp_value,
            MD_cycles,
            namespace,
        )

        if len(combinations) > 1:
            print("")
            print(
                f"Run {combo_index}/{len(combinations)}: "
                f"T={temperature_value} K, "
                f"dt={timestep_value} fs, "
                f"T-damp={temperature_damp_value} fs"
            )

        energy_records = _run_md_engine(
            engine=Engine,
            atoms=asestruct,
            struct_prefix=struct_prefix,
            openkim_potential=OpenKIM_potential,
            temperature_profile=temperature_profile,
            timestep_profile=timestep_profile,
            temperature_damp_profile=temperature_damp_profile,
            md_cycles=MD_cycles,
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
            nn=0
            mm=0

            if Solve_double_element_problem == True:
                for n in asestruct.get_chemical_symbols():
                    nn=nn+1
                    for m in positions:
                        mm=mm+1
                        if mm == nn:
                            f.write("    Atom('"+n+"', ( "+str(m[0])+", "+str(m[1])+", "+str(m[2])+" )),\n")
                    mm=0
            else:
                for n in asestruct.get_chemical_symbols():
                    for m in positions:
                        f.write("    Atom('"+n+"', ( "+str(m[0])+", "+str(m[1])+", "+str(m[2])+" )),\n")
            f.write("    ],\n")
            f.write("    cell=[("+str(asestruct.cell[0,0])+", "+str(asestruct.cell[0,1])+", "+str(asestruct.cell[0,2])+"), ("+str(asestruct.cell[1,0])+", "+str(asestruct.cell[1,1])+", "+str(asestruct.cell[1,2])+"), ("+str(asestruct.cell[2,0])+", "+str(asestruct.cell[2,1])+", "+str(asestruct.cell[2,2])+")],\n")
            if Manual_PBC == False:
                f.write("    pbc=True,\n")
            else:
                f.write("    pbc=["+str(PBC_constraints[0])+","+str(PBC_constraints[1])+","+str(PBC_constraints[2])+"],\n")
            f.write("    )\n")

        _print_attention_message()
        _export_cif(struct_prefix, asestruct)

    namespace['Temperature'] = original_temperature
    namespace['Time_step'] = original_timestep
    namespace['Temperature_damp'] = original_temperature_damp

if __name__ == "__main__":
    main()
