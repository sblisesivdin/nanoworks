#!/usr/bin/env python

'''
dftsolve: High-level Wrapper Script for GPAW
More information: $ dftsolve -h
'''

Description = f'''
 Usage:
 $ dftsolve -p <corenumbers> <args>
'''

import sys
import os, glob
import gc
import importlib.util
import json
import shlex
import shutil
import subprocess

def log_energy_consumption(meter, struct_name, engine):
    """
    Reads pyRAPL results, write energy consumption to a file MPI-safe 
    """
    if meter is None:
        return
        
    from ase.parallel import paropen
    
    meter.end()
    energyresult = meter.result
   
    # Safely handle potential None values from pyRAPL
    duration = energyresult.duration if energyresult.duration is not None else 0.0
    pkg_energy = sum(energyresult.pkg) if energyresult.pkg is not None else 0.0
    dram_energy = sum(energyresult.dram) if energyresult.dram is not None else 0.0
   
    with paropen(
        struct_name
        + f'-ENERGY-{engine}-Log-Energy_consumption.txt',
        'a',
    ) as f1:
        print("Energy measurement:-----------------------------------------", end="\n", file=f1)
        print(1e-6 * duration, " Computation time in seconds", end="\n", file=f1)
        print(1e-6 * pkg_energy, " CPU energy consumption in Joules", end="\n", file=f1)
       
        if energyresult.dram is None:
            print("0.0 DRAM energy consumption in Joules (Hardware domain unsupported)", end="\n", file=f1)
        else:
            print(1e-6 * dram_energy, " DRAM energy consumption in Joules", end="\n", file=f1)
           
        print(2.77777778e-7 * (1e-6 * pkg_energy + 1e-6 * dram_energy), " Total energy consumption in kWh", end="\n", file=f1)

# Parallel execution request.
#
# GPAW runs the whole dftsolve process under MPI.
# Quantum ESPRESSO keeps dftsolve serial and launches the QE
# executable itself under MPI.
def extract_parallel_request():
    parallel = None
    filtered_args = []

    i = 1

    while i < len(sys.argv):
        arg = sys.argv[i]

        if arg in ['-p', '--parallel']:
            if i + 1 >= len(sys.argv):
                print(
                    "Error: -p/--parallel requires an argument"
                )
                sys.exit(1)

            parallel = sys.argv[i + 1]
            i += 2

        elif arg.startswith('--parallel='):
            parallel = arg.split(
                '=',
                1,
            )[1]
            i += 1

        elif arg.startswith('-p') and arg != '-p':
            parallel = arg[2:]
            i += 1

        else:
            filtered_args.append(arg)
            i += 1

    if parallel is not None:
        try:
            parallel = int(parallel)
        except ValueError:
            print(
                "Error: -p/--parallel must be a positive integer"
            )
            sys.exit(1)

        if parallel <= 0:
            print(
                "Error: -p/--parallel must be a positive integer"
            )
            sys.exit(1)

    return parallel, filtered_args


GPAW_STAGE_GROUP_ENV = 'NANOWORKS_GPAW_STAGE_GROUP'


def build_gpaw_process_command(
    parallel,
    filtered_args,
):
    """Build a fresh GPAW process command for the requested core count."""
    script_path = os.path.abspath(
        __file__
    )

    if parallel is None:
        return [
            sys.executable,
            script_path,
            *filtered_args,
        ]

    mpi_exe = (
        shutil.which('mpiexec')
        or shutil.which('mpirun')
        or shutil.which('srun')
    )

    if mpi_exe is None:
        raise RuntimeError(
            "mpiexec, mpirun, or srun not found for "
            "parallel GPAW execution."
        )

    flag = (
        '-n'
        if 'srun' in os.path.basename(mpi_exe)
        else '-np'
    )

    gpaw_exe = shutil.which('gpaw')

    if gpaw_exe is None:
        raise RuntimeError(
            "GPAW command was not found in PATH."
        )
    
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
    os.environ['NUMEXPR_NUM_THREADS'] = '1'

    return [
        mpi_exe,
        flag,
        str(parallel),
        gpaw_exe,
        'python',
        '--',
        script_path,
    ] + filtered_args


def restart_gpaw_with_mpi(
    parallel,
    filtered_args,
):
    """Restart dftsolve under MPI for the GPAW backend."""
    try:
        cmd = build_gpaw_process_command(
            parallel,
            filtered_args,
        )
    except RuntimeError as exc:
        print(
            'Error: ' + str(exc)
        )
        sys.exit(1)

    print(
        f"Restarting GPAW calculation with "
        f"{parallel} cores: {' '.join(cmd)}"
    )

    sys.stdout.flush()
    
    child_env = os.environ.copy()

    # "gpaw python" configures the GPAW MPI backend itself.
    # Do not pass the backend selected for the original serial
    # Nanoworks process to the MPI child processes.
    child_env.pop(
        'GPAW_MPI_BACKEND',
        None,
    )

    os.execvpe(
        cmd[0],
        cmd,
        child_env,
    )


REQUESTED_PARALLEL, FILTERED_ARGS = (
    extract_parallel_request()
)

# argparse must not see Nanoworks' process-count argument.
sys.argv = [
    sys.argv[0],
    *FILTERED_ARGS,
]

# Since gpaw-python is removed, the standard python interpreter needs this 
# environment variable to enable MPI parallelization when imported as a library.
if "GPAW_MPI_BACKEND" not in os.environ:
    os.environ["GPAW_MPI_BACKEND"] = "cgpaw" # or "mpi4py"

import getopt, time
import textwrap
import pickle
import nanoworks
from nanoworks.engine import (
    normalize_engine_name,
    resolve_calculation_stages,
    resolve_initial_magnetic_moments,
    resolve_stage_kpoint_settings,
    resolve_stage_occupation,
    load_engine_module,
)
from nanoworks.pseudos import (
    get_qe_pseudo_dir,
    resolve_qe_pseudopotentials,
)
from argparse import ArgumentParser, HelpFormatter
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any
from ase import *
from ase.spacegroup import get_spacegroup
from ase.dft.kpoints import get_special_points
from ase.parallel import paropen, world, parprint, broadcast
from ase.optimize import QuasiNewton
from ase.io import read, write
from ase.calculators.singlepoint import SinglePointCalculator
from ase.eos import calculate_eos
from ase.units import Bohr, GPa, kJ
import matplotlib.pyplot as plt
from ase.constraints import FixSymmetry
from ase.filters import FrechetCellFilter
from ase.io.cif import write_cif
from pathlib import Path
import numpy as np
from numpy import genfromtxt
import warnings
warnings.filterwarnings('ignore')

DFT_ENGINE_DEFAULTS = {
    'GPAW': {
        'XC_calc': 'LDA',
        'Opt_calc_type': 'BSE',
        'DOS_occupation': None,
        'Fix_symmetry': False,
        'Phonon_PW_cutoff': 400,
        'Phonon_kpts_x': 3,
        'Phonon_kpts_y': 3,
        'Phonon_kpts_z': 3,
    },
    'QE': {
        'XC_calc': 'PBE',
        'Opt_calc_type': 'RPA',
        'DOS_occupation': 'tetrahedra',
        'Fix_symmetry': True,
        'Phonon_PW_cutoff': None,
        'Phonon_kpts_x': None,
        'Phonon_kpts_y': None,
        'Phonon_kpts_z': None,
    },
}

@dataclass
class DFTConfig:
    """
    Configuration dataclass to hold all DFT calculation parameters.
    """
    # Engine, mode and calculation flags
    Engine: str = 'GPAW'
    Mode: str = 'PW'
    Ground_calc: bool = False
    Geo_optim: bool = False
    Elastic_calc: bool = False
    DOS_calc: bool = False
    Band_calc: bool = False
    Density_calc: bool = False
    Phonon_calc: bool = False
    Optical_calc: bool = False
    SOC_calc: bool = False
    vdW_calc: str = 'None'
    
    # Geometry optimization parameters
    Optimizer: str = 'QuasiNewton'
    Max_F_tolerance: float = 0.05
    Max_step: float = 0.1
    Alpha: float = 60.0
    Damping: float = 1.0
    Fix_symmetry: bool = None
    Relax_cell: List[bool] = field(default_factory=lambda: [False, False, False, False, False, False])
    Hydrostatic_pressure: float = 0.0
    
    # Elastic parameters
    Elastic_kpts_density: Optional[float] = None
    Elastic_kpts_x: Optional[int] = None
    Elastic_kpts_y: Optional[int] = None
    Elastic_kpts_z: Optional[int] = None
    Elastic_gamma: Optional[bool] = None
    
    # Ground state parameters
    Cut_off_energy: float = 340
    Ground_num_of_bands: Optional[int] = None
    Ground_gamma: Optional[bool] = None
    Ground_kpts_density: Optional[float] = None
    Ground_kpts_x: int = 5
    Ground_kpts_y: int = 5
    Ground_kpts_z: int = 5
    Ground_gpts_density: Optional[float] = None
    Ground_gpts_x: int = 8
    Ground_gpts_y: int = 8
    Ground_gpts_z: int = 8
    Setup_params: Dict = field(default_factory=dict)
    XC_calc: Any = None
    # Optional hybrid (HSE06/HSE03/PBE0/B3LYP/EXX) tuning. When left as None,
    # GPAW's documented defaults for each functional are used (e.g. HSE06 uses
    # omega=0.11 1/Bohr and 25% exact exchange).
    XC_exx_fraction: Optional[float] = None
    XC_omega: Optional[float] = None
    XC_backend: str = 'pw'
    Ground_convergence: Dict = field(default_factory=dict)
    Occupation: Dict = field(default_factory=lambda: {'name': 'fermi-dirac', 'width': 0.05})
    Mixer_type: Any = None
    Spin_calc: bool = False
    Magmom_per_atom: Any = 1.0
    Magmom_single_atom: Optional[List] = None
    Total_charge: float = 0.0
    
    # DOS parameters
    DOS_npoints: int = 501
    DOS_width: float = 0.1
    DOS_convergence: Dict = field(default_factory=dict)
    DOS_num_of_bands: Optional[int] = None
    DOS_kpts_density: Optional[float] = None
    DOS_kpts_x: Optional[int] = None
    DOS_kpts_y: Optional[int] = None
    DOS_kpts_z: Optional[int] = None
    DOS_gamma: Optional[bool] = None
    DOS_occupation: Any = None
    
    # Band structure parameters
    Gamma: bool = True
    Band_path: str = 'LGL'
    Band_npoints: int = 61
    Band_num_of_bands: Optional[int] = None
    Energy_max: float = 5
    Energy_min: float = -5
    Band_convergence: Dict = field(default_factory=lambda: {'bands': 8})
    Projected_band_plot: bool = False
    Projections: List[Dict[str, Any]] = field(default_factory=list)
    
    # Electron density parameters
    Refine_grid: int = 4
    
    # Phonon parameters
    Phonon_PW_cutoff: Optional[float] = None
    Phonon_kpts_x: Optional[int] = None
    Phonon_kpts_y: Optional[int] = None
    Phonon_kpts_z: Optional[int] = None
    Phonon_supercell: Any = None
    Phonon_displacement: float = 1e-3
    Phonon_path: str = 'LGL'
    Phonon_npoints: int = 61
    Phonon_acoustic_sum_rule: bool = True
    Phonon_qpts_x: int = 20
    Phonon_qpts_y: int = 20
    Phonon_qpts_z: int = 20
    Phonon_thermal_calc: bool = False
    Phonon_T_min: float = 0.0
    Phonon_T_max: float = 1000.0
    Phonon_T_step: float = 10.0
    
    # Optical parameters
    Opt_calc_type: Optional[str] = None
    Opt_shift_en: float = 0.0
    Opt_BSE_valence: Any = None
    Opt_BSE_conduction: Any = None
    Opt_BSE_min_en: float = 0.0
    Opt_BSE_max_en: float = 20.0
    Opt_BSE_num_of_data: int = 1001
    Opt_min_en: Optional[float] = None
    Opt_max_en: Optional[float] = None
    Opt_num_of_data: Optional[int] = None
    Opt_num_of_bands: int = 8
    Opt_kpts_density: Optional[float] = None
    Opt_kpts_x: Optional[int] = None
    Opt_kpts_y: Optional[int] = None
    Opt_kpts_z: Optional[int] = None
    Opt_gamma: Optional[bool] = None
    Opt_FD_smearing: float = 0.05
    Opt_eta: float = 0.05
    Opt_domega0: float = 0.05
    Opt_omega2: float = 5.0
    Opt_cut_of_energy: float = 100
    Opt_nblocks: Any = None
    
    # General parameters
    Localization: str = "en_UK"
    Outdirname: str = ''
    
    # Bulk configuration
    bulk_configuration: Any = None
    
    def __post_init__(self):
        """Initialize default values that depend on other objects."""
        self.Engine = normalize_engine_name(self.Engine)

        try:
            engine_defaults = (
                DFT_ENGINE_DEFAULTS[
                    self.Engine
                ]
            )
        except KeyError:
            raise ValueError(
                "Unsupported DFT engine: "
                f"{self.Engine}"
            )

        for name, value in engine_defaults.items():
            if getattr(
                self,
                name,
            ) is None:
                setattr(
                    self,
                    name,
                    value,
                )

        if self.Mixer_type is None and self.Engine == 'GPAW':
            engine = load_engine_module(self.Engine)
            self.Mixer_type = engine.create_default_mixer()
        if self.Phonon_supercell is None:
            self.Phonon_supercell = np.diag([2, 2, 2])
        if self.Opt_BSE_valence is None:
            self.Opt_BSE_valence = range(0, 3)
        if self.Opt_BSE_conduction is None:
            self.Opt_BSE_conduction = range(4, 7)
        if self.Opt_min_en is None:
            self.Opt_min_en = self.Opt_BSE_min_en
        if self.Opt_max_en is None:
            self.Opt_max_en = self.Opt_BSE_max_en
        if self.Opt_num_of_data is None:
            self.Opt_num_of_data = self.Opt_BSE_num_of_data
        if self.Opt_nblocks is None:
            self.Opt_nblocks = world.size
        
        sanitized_projections = []
        
        for idx, proj in enumerate(self.Projections):
            if isinstance(proj, dict):
                # Rebuild the dictionary with safe defaults if keys are missing
                safe_proj = {
                    'atoms': proj.get('atoms', []),
                    'orbital': proj.get('orbital', None),
                    'color': proj.get('color', 'blue'),  # Default to blue if missing
                    'label': proj.get('label', f"Proj-{idx+1}") # Default label if missing
                }
                sanitized_projections.append(safe_proj)
        
        # Replace the user's raw list with the safely formatted list
        self.Projections = sanitized_projections

class RawFormatter(HelpFormatter):
    """To print Description variable with argparse"""
    def _fill_text(self, text, width, indent):
        return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

def struct_from_file(
    inputfile,
    geometryfile,
    create_output=True,
    report_structure=True,
):
    """Load variables from parse function and return DFTConfig instance."""
    # Works like from FILE import *
    sys.path.append(str(Path(inputfile).parent))
    inputf = __import__(Path(inputfile).stem, globals(), locals(), ['*'])
    
    # Create a config object with loaded parameters
    config_dict = {}
    for k in dir(inputf):
        if not k.startswith('_'):
            config_dict[k] = getattr(inputf, k)

    
    # Create DFTConfig instance
    config = DFTConfig(**{k: v for k, v in config_dict.items() if k in DFTConfig.__dataclass_fields__})
    
    # If there is a CIF input, use it. Otherwise use the bulk configuration provided above.
    if geometryfile is None:
        if config.Outdirname != '':
            struct = config.Outdirname
        else:
            struct = 'results' # All files will get their names from this file
    else:
        struct = Path(geometryfile).stem
        config.bulk_configuration = read(geometryfile, index='-1')
        if report_structure:
            parprint("Number of atoms imported from CIF file:"+str(config.bulk_configuration.get_global_number_of_atoms()))
            parprint("Spacegroup of CIF file:",get_spacegroup(config.bulk_configuration, symprec=1e-2))
            parprint("Special Points usable for this spacegroup:",get_special_points(config.bulk_configuration.get_cell()))

    # Output directory
    input_dir = Path(inputfile).parent
    if config.Outdirname != '':
        structpath = input_dir / config.Outdirname
    else:
        structpath = input_dir / struct

    if create_output and not os.path.isdir(structpath):
        os.makedirs(structpath, exist_ok=True)
    struct = os.path.join(str(structpath), struct)
    return struct, config

def struct_from_auto(
    geometryfile,
    write_output=True,
    report_structure=True,
):
    """Generate configuration automatically from geometry file."""
    struct_path = Path(geometryfile)
    struct_name = struct_path.stem
    
    # Generate Auto Config
    atoms = read(geometryfile)
    config = DFTConfig()
    config.Mode = 'PW'
    config.Ground_calc = True
    config.XC_calc = 'PBE'
    config.Cut_off_energy = 450
    config.Gamma = True
    config.Optimizer = 'LBFGS' 
    
    # ---------------------------------------------------------
    # 1. Magnetism Detection
    # ---------------------------------------------------------
    mag_elements = ['Cr', 'Mn', 'Fe', 'Co', 'Ni']
    symbols = set(atoms.get_chemical_symbols())
    
    if any(s in mag_elements for s in symbols):
        config.Spin_calc = True
        config.Magmom_per_atom = 2.0 

    # ---------------------------------------------------------
    # 2. Dimensionality & K-points (Vacuum detection)
    # ---------------------------------------------------------
    # We want to avoid many k-points in vacuum directions.
    # Heuristic: If cell length is large (>12 A) and atoms occupy < 50% of it, it's vacuum.
    
    atoms_an = atoms.copy()
    atoms_an.center() # Center atoms to get correct span
    positions = atoms_an.get_positions()
    if len(atoms) > 0:
        spans = np.max(positions, axis=0) - np.min(positions, axis=0)
    else:
        spans = np.zeros(3)
    cell_lengths = atoms.cell.lengths()
    
    # Target density (pts/A^-1)
    target_density = 3.0
    kpts = [1, 1, 1]
    
    is_vacuum = [False, False, False]
    
    for i in range(3):
        # Check for vacuum
        if cell_lengths[i] > 12.0 and spans[i] < (cell_lengths[i] * 0.5):
            is_vacuum[i] = True
            kpts[i] = 1
        else:
            # Calculate k-points based on density
            # density = N * 2pi / L  => N = density * L / 2pi
            # We add 0.5 to round nearest, and ensure at least 1
            n_k = int(target_density * cell_lengths[i] / (2 * np.pi) + 0.5)
            kpts[i] = max(1, n_k)
            
    config.Ground_kpts_x = kpts[0]
    config.Ground_kpts_y = kpts[1]
    config.Ground_kpts_z = kpts[2]
    config.Ground_kpts_density = None # Use explicit kpts
    
    # Enable other calculations
    config.DOS_calc = True
    config.Band_calc = True
    
    config.bulk_configuration = atoms
    
    # Determine output path
    input_dir = struct_path.parent
    structpath = input_dir / struct_name
    
    if write_output and not os.path.isdir(structpath):
        os.makedirs(structpath, exist_ok=True)
    struct = os.path.join(str(structpath), struct_name)
    
    # Write input file for future use
    input_filename = (
        struct
        + '-CONFIG-NANOWORKS-Input-Auto.py'
    )
    if write_output:
        with open(input_filename, 'w') as f:
            f.write("from ase.io import read\n")
            f.write("import numpy as np\n\n")
            f.write(f"Mode = '{config.Mode}'\n")
            f.write(f"Ground_calc = {config.Ground_calc}\n")
            f.write(f"XC_calc = '{config.XC_calc}'\n")
            f.write(f"XC_exx_fraction = {getattr(config, 'XC_exx_fraction', None)}\n")
            f.write(f"XC_omega = {getattr(config, 'XC_omega', None)}\n")
            f.write(f"XC_backend = '{getattr(config, 'XC_backend', 'pw')}'\n")
            f.write(f"Cut_off_energy = {config.Cut_off_energy}\n")
            f.write(f"Gamma = {config.Gamma}\n")
            f.write(f"Optimizer = '{config.Optimizer}'\n")
            f.write(f"Spin_calc = {config.Spin_calc}\n")
            if config.Spin_calc:
                f.write(f"Magmom_per_atom = {config.Magmom_per_atom}\n")
            f.write(f"Ground_kpts_x = {config.Ground_kpts_x}\n")
            f.write(f"Ground_kpts_y = {config.Ground_kpts_y}\n")
            f.write(f"Ground_kpts_z = {config.Ground_kpts_z}\n")
            f.write(f"DOS_calc = {config.DOS_calc}\n")
            f.write(f"Band_calc = {config.Band_calc}\n")
            f.write(f"\n# Geometry is handled via command line -g or loaded here if needed\n")

    if report_structure:
        parprint(f"Auto-configured for {struct_name}: PBE, 450eV, Spin={config.Spin_calc}")
        parprint(f"Geometry analysis: Cell {cell_lengths}, Spans {spans}")
        parprint(f"Vacuum detected: {is_vacuum} -> K-points set to {kpts}")
        if write_output:
            parprint(f"Generated input file: {input_filename}")
    
    return struct, config

def autoscale_y(ax, margin=0.1):
    """
    Automatically scales the Y-axis of a matplotlib ax based on the visible X-limits.
    """
    import numpy as np

    def get_bottom_top(line):
        # Force the data into numpy arrays to allow boolean masking
        xd = np.asarray(line.get_xdata())
        yd = np.asarray(line.get_ydata())

        lo, hi = ax.get_xlim()

        # Create a boolean mask for the visible region
        mask = (xd > lo) & (xd < hi)
        y_displayed = yd[mask]

        # Safety guard: if there is no data in the visible range, skip
        if len(y_displayed) == 0:
            return np.nan, np.nan

        h = np.max(y_displayed) - np.min(y_displayed)
        bot = np.min(y_displayed) - margin * h
        top = np.max(y_displayed) + margin * h
        return bot, top

    # Gather limits for all lines plotted in the axis
    lines = ax.lines
    bot_top = [get_bottom_top(line) for line in lines]

    # Filter out empty/invalid bounds
    bot_top = [bt for bt in bot_top if not np.isnan(bt[0])]

    # Apply the new global min and max to the Y-axis
    if len(bot_top) > 0:
        bots, tops = zip(*bot_top)
        ax.set_ylim(min(bots), max(tops))

class dftsolve:
    """
    The dftsolve class is a high-level interaction script for GPAW calculations.
    It handles various types of calculations such as ground state, structure optimization,
    elastic properties, density of states, band structure, density, and optical properties.
    The class takes input parameters from a DFTConfig instance and performs the calculations
    accordingly.
    """
    def __init__(
        self,
        struct: str,
        config: DFTConfig,
        parallel_cores: int = 1,
    ):
        """Initialize dftsolve with struct path and configuration.
        
        Args:
            struct: Path to structure files
            config: DFTConfig instance containing all calculation parameters
        """
        self.struct = struct
        self.config = config
        
        # For backward compatibility, expose config attributes as instance attributes
        self.Engine = config.Engine
        self.engine = load_engine_module(self.Engine)
        self.parallel_cores = int(
            parallel_cores
        )
        self.Mode = config.Mode
        self.Ground_calc = config.Ground_calc
        self.Geo_optim = config.Geo_optim
        self.Elastic_calc = config.Elastic_calc
        self.DOS_calc = config.DOS_calc
        self.Band_calc = config.Band_calc
        self.Density_calc = config.Density_calc
        self.Optical_calc = config.Optical_calc
        self.SOC_calc = config.SOC_calc
        self.vdW_calc = config.vdW_calc
        self.Optimizer = config.Optimizer
        self.Max_F_tolerance = config.Max_F_tolerance
        self.Max_step = config.Max_step
        self.Alpha = config.Alpha
        self.Damping = config.Damping
        self.Fix_symmetry = config.Fix_symmetry
        self.Relax_cell = config.Relax_cell
        self.Hydrostatic_pressure = config.Hydrostatic_pressure
        self.Elastic_kpts_density = config.Elastic_kpts_density
        self.Elastic_kpts_x = config.Elastic_kpts_x
        self.Elastic_kpts_y = config.Elastic_kpts_y
        self.Elastic_kpts_z = config.Elastic_kpts_z
        self.Elastic_gamma = config.Elastic_gamma
        self.Cut_off_energy = config.Cut_off_energy
        self.Ground_gamma = config.Ground_gamma
        self.Ground_kpts_density = config.Ground_kpts_density
        self.Ground_kpts_x = config.Ground_kpts_x
        self.Ground_kpts_y = config.Ground_kpts_y
        self.Ground_kpts_z = config.Ground_kpts_z
        self.Ground_gpts_density = config.Ground_gpts_density
        self.Ground_gpts_x = config.Ground_gpts_x
        self.Ground_gpts_y = config.Ground_gpts_y
        self.Ground_gpts_z = config.Ground_gpts_z
        self.Setup_params = config.Setup_params
        self.XC_calc = config.XC_calc
        self.XC_exx_fraction = getattr(config, 'XC_exx_fraction', None)
        self.XC_omega = getattr(config, 'XC_omega', None)
        self.XC_backend = getattr(config, 'XC_backend', 'pw')
        # Fermi level (eV) of the converged ground state. Stored here so that
        # the DOS and band methods can reference hybrid eigenvalues correctly
        # instead of hard-coding 0.0 eV.
        self.Ground_fermi_level = None
        self.Ground_convergence = config.Ground_convergence
        self.Ground_num_of_bands = config.Ground_num_of_bands
        self.Occupation = config.Occupation
        self.Mixer_type = config.Mixer_type
        self.Spin_calc = config.Spin_calc
        self.Magmom_per_atom = config.Magmom_per_atom
        self.Magmom_single_atom = config.Magmom_single_atom
        self.Total_charge = config.Total_charge
        self.DOS_npoints = config.DOS_npoints
        self.DOS_width = config.DOS_width
        self.DOS_convergence = config.DOS_convergence
        self.DOS_kpts_density = config.DOS_kpts_density
        self.DOS_kpts_x = config.DOS_kpts_x
        self.DOS_kpts_y = config.DOS_kpts_y
        self.DOS_kpts_z = config.DOS_kpts_z
        self.DOS_gamma = config.DOS_gamma
        self.DOS_occupation = config.DOS_occupation
        self.DOS_num_of_bands = config.DOS_num_of_bands
        self.Gamma = config.Gamma
        self.Band_path = config.Band_path
        self.Projected_band_plot = config.Projected_band_plot
        self.Projections = config.Projections
        self.Band_npoints = config.Band_npoints
        self.Band_num_of_bands = config.Band_num_of_bands
        self.Energy_max = config.Energy_max
        self.Energy_min = config.Energy_min
        self.Band_convergence = config.Band_convergence
        self.Refine_grid = config.Refine_grid
        self.Phonon_PW_cutoff = config.Phonon_PW_cutoff
        self.Phonon_kpts_x = config.Phonon_kpts_x
        self.Phonon_kpts_y = config.Phonon_kpts_y
        self.Phonon_kpts_z = config.Phonon_kpts_z
        self.Phonon_supercell = config.Phonon_supercell
        self.Phonon_displacement = config.Phonon_displacement
        self.Phonon_path = config.Phonon_path
        self.Phonon_npoints = config.Phonon_npoints
        self.Phonon_acoustic_sum_rule = config.Phonon_acoustic_sum_rule
        self.Phonon_qpts_x = config.Phonon_qpts_x
        self.Phonon_qpts_y = config.Phonon_qpts_y
        self.Phonon_qpts_z = config.Phonon_qpts_z
        self.Phonon_thermal_calc = config.Phonon_thermal_calc
        self.Phonon_T_min = config.Phonon_T_min
        self.Phonon_T_max = config.Phonon_T_max
        self.Phonon_T_step = config.Phonon_T_step
        self.Opt_calc_type = config.Opt_calc_type
        self.Opt_shift_en = config.Opt_shift_en
        self.Opt_BSE_valence = config.Opt_BSE_valence
        self.Opt_BSE_conduction = config.Opt_BSE_conduction
        self.Opt_BSE_min_en = config.Opt_BSE_min_en
        self.Opt_BSE_max_en = config.Opt_BSE_max_en
        self.Opt_BSE_num_of_data = config.Opt_BSE_num_of_data
        self.Opt_min_en = config.Opt_min_en
        self.Opt_max_en = config.Opt_max_en
        self.Opt_num_of_data = config.Opt_num_of_data
        self.Opt_num_of_bands = config.Opt_num_of_bands
        self.Opt_kpts_density = config.Opt_kpts_density
        self.Opt_kpts_x = config.Opt_kpts_x
        self.Opt_kpts_y = config.Opt_kpts_y
        self.Opt_kpts_z = config.Opt_kpts_z
        self.Opt_gamma = config.Opt_gamma
        self.Opt_FD_smearing = config.Opt_FD_smearing
        self.Opt_eta = config.Opt_eta
        self.Opt_domega0 = config.Opt_domega0
        self.Opt_omega2 = config.Opt_omega2
        self.Opt_cut_of_energy = config.Opt_cut_of_energy
        self.Opt_nblocks = config.Opt_nblocks
        self.bulk_configuration = config.bulk_configuration
        
        from nanoworks.localization import Translator
        
        # Control localisation/localization
        if hasattr(self.config, 'Localization'):
            current_lang = self.config.Localization
        elif hasattr(self.config, 'Localisation'):
            current_lang = self.config.Localisation
        else:
            current_lang = "en"  # İkisi de yoksa İngilizce
            
        # translator function
        self._t = Translator(lang_code=current_lang).get

    def _load_existing_final_structure(self):
        """Load a saved final structure for subsequent calculations."""
        final_structure_file = Path(
            self.struct
            + f'-GROUND-{self.Engine}-Result-Final.cif'
        )

        if not final_structure_file.is_file():
            return False

        try:
            final_structure = read(
                final_structure_file,
                index='-1',
            )
        except Exception as exc:
            raise RuntimeError(
                "The existing final structure could not be loaded "
                f"from '{final_structure_file}': {exc}"
            ) from exc

        self.bulk_configuration = final_structure
        self.config.bulk_configuration = final_structure

        parprint(
            "\033[93mWARNING:\033[0m "
            "Ground_calc = False and an existing final structure "
            "was found. Subsequent calculations will use the "
            "geometry from: "
            + str(final_structure_file)
        )

        return True

    def structurecalc(self):
        """
        This method calculates and writes the spacegroup and special points of the given structure.
        It reads the bulk configuration from the CIF file and prints the number of atoms, spacegroup,
        and special points usable for the spacegroup to a text file.
        """

        # -------------------------------------------------------------
        # STRUCTURE
        # -------------------------------------------------------------

        with paropen(self.struct+'-STRUCTURE-ASE-Result-Spacegroup-and-SpecialPoints.txt', "w") as fd:
            print("Number of atoms imported from CIF file:"+str(self.bulk_configuration.get_global_number_of_atoms()), file=fd)
            print("Spacegroup of CIF file:",get_spacegroup(self.bulk_configuration, symprec=1e-2), file=fd)
            print("Special Points usable for this spacegroup:",get_special_points(self.bulk_configuration.get_cell()), file=fd)

    def groundcalc(self):
        """Run the ground-state workflow using the selected DFT engine."""
        time11 = time.time()

        if self.Engine == 'GPAW':
            self._groundcalc_gpaw()

        elif self.Engine == 'QE':
            self._groundcalc_qe()

        else:
            raise ValueError(
                f"Unsupported DFT engine: {self.Engine}"
            )

        time12 = time.time()

        with paropen(
            self.struct
            + f'-TIMINGS-{self.Engine}-Log-Timings.txt',
            'a',
        ) as f1:
            print(
                'Ground state: ',
                round(
                    time12 - time11,
                    2,
                ),
                end="\n",
                file=f1,
            )

    def _groundcalc_qe(self):
        """Run the Quantum ESPRESSO PW ground-state workflow."""

        # -------------------------------------------------------------
        # GROUND STATE - QE
        # -------------------------------------------------------------

        if self.Mode != 'PW':
            parprint(
                "\033[91mERROR:\033[0m "
                "Quantum ESPRESSO backend currently supports "
                "PW mode only."
            )
            sys.exit(1)

        if self.config.vdW_calc.upper() != 'NONE':
            parprint(
                "\033[91mERROR:\033[0m "
                "vdW corrections are not implemented "
                "for the QE backend yet."
            )
            sys.exit(1)

        try:
            self.engine.validate_qe_xc(
                self.XC_calc,
                pseudo_xc='pbe',
                allow_hybrid=True,
            )
        except ValueError as exc:
            parprint(
                f"\033[91mERROR:\033[0m {exc}"
            )
            sys.exit(1)

        state_dir = Path(
            self.struct
            + '-GROUND-QE-Result-State'
        )

        final_structure_file = Path(
            self.struct
            + f'-GROUND-{self.Engine}-Result-Final.cif'
        )

        if not self.Ground_calc:
            parprint(
                "Passing QE PW ground state calculation..."
            )

            if not self.engine.has_qe_state(
                state_dir,
                prefix='nanoworks',
            ):
                parprint(
                    "\033[91mERROR:\033[0m "
                    + str(state_dir)
                    + " does not contain a valid QE ground-state result. "
                    "It is needed by subsequent QE calculations. "
                    "First complete the ground-state calculation with "
                    "Ground_calc = True. Exiting."
                )
                sys.exit(1)

            final_structure_loaded = (
                self._load_existing_final_structure()
            )

            if (
                self.Geo_optim
                and not final_structure_loaded
            ):
                parprint(
                    "\033[91mERROR:\033[0m "
                    "The optimized QE structure could not be found: "
                    + str(final_structure_file)
                )
                parprint(
                    "Run the ground-state calculation once with "
                    "Ground_calc = True and Geo_optim = True."
                )
                sys.exit(1)

            return

        pseudo_dir = get_qe_pseudo_dir(
            relativistic='scalar',
        )

        try:
            pseudopotentials = (
                resolve_qe_pseudopotentials(
                    self.bulk_configuration,
                    relativistic='scalar',
                )
            )
        except (FileNotFoundError, RuntimeError) as exc:
            parprint(
                f"\033[91mERROR:\033[0m {exc}"
            )
            sys.exit(1)

        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )

        magnetic_moments = None

        if self.Spin_calc:
            magnetic_moments = (
                resolve_initial_magnetic_moments(
                    atoms=self.bulk_configuration,
                    magmom_per_atom=self.Magmom_per_atom,
                    magmom_single_atom=(
                        self.Magmom_single_atom
                    ),
                )
            )

        if self.Geo_optim:
            variable_cell = (
                True in self.Relax_cell
            )

            relaxation_label = (
                'VC-RELAX'
                if variable_cell
                else 'RELAX'
            )

            if variable_cell:
                parprint(
                    "Starting QE PW variable-cell "
                    "geometry optimization..."
                )
            else:
                parprint(
                    "Starting QE PW atomic "
                    "geometry optimization..."
                )

            input_file = Path(
                self.struct
                + f'-GROUND-QE-Input-{relaxation_label}.in'
            )

            output_file = Path(
                self.struct
                + f'-GROUND-{self.Engine}-Log-'
                + f'{relaxation_label}.txt'
            )

            try:
                workflow = self.engine.run_relax(
                    atoms=self.bulk_configuration,
                    input_file=input_file,
                    output_file=output_file,
                    state_dir=state_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    cutoff_ev=self.Cut_off_energy,
                    optimizer=self.Optimizer,
                    max_force=self.Max_F_tolerance,
                    max_step=self.Max_step,
                    relax_cell=self.Relax_cell,
                    hydrostatic_pressure=(
                        self.Hydrostatic_pressure
                    ),
                    fix_symmetry=self.Fix_symmetry,
                    kpoint_density=self.Ground_kpts_density,
                    kpoint_size=(
                        self.Ground_kpts_x,
                        self.Ground_kpts_y,
                        self.Ground_kpts_z,
                    ),
                    gamma=ground_gamma,
                    total_charge=self.Total_charge,
                    nbands=self.Ground_num_of_bands,
                    spinpol=self.Spin_calc,
                    magnetic_moments=magnetic_moments,
                    setup_params=self.Setup_params,
                    xc_calc=self.XC_calc,
                    exx_fraction=self.XC_exx_fraction,
                    omega=self.XC_omega,
                    occupation=self.Occupation,
                    parallel_cores=self.parallel_cores,
                    executable='pw.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    f"QE geometry optimization failed: {exc}"
                )
                raise

            self.bulk_configuration = workflow[
                'atoms'
            ]

            self.config.bulk_configuration = (
                self.bulk_configuration
            )

            result = workflow[
                'result'
            ]

            if variable_cell:
                parprint(
                    "QE variable-cell geometry "
                    "optimization finished."
                )
            else:
                parprint(
                    "QE atomic geometry optimization finished."
                )

        else:
            parprint(
                "Starting QE PW ground state calculation..."
            )

            input_file = Path(
                self.struct
                + '-GROUND-QE-Input-SCF.in'
            )

            output_file = Path(
                self.struct
                + f'-GROUND-{self.Engine}-Log-SCF.txt'
            )

            try:
                workflow = self.engine.run_scf(
                    atoms=self.bulk_configuration,
                    input_file=input_file,
                    output_file=output_file,
                    state_dir=state_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    cutoff_ev=self.Cut_off_energy,
                    kpoint_density=self.Ground_kpts_density,
                    kpoint_size=(
                        self.Ground_kpts_x,
                        self.Ground_kpts_y,
                        self.Ground_kpts_z,
                    ),
                    gamma=ground_gamma,
                    total_charge=self.Total_charge,
                    nbands=self.Ground_num_of_bands,
                    spinpol=self.Spin_calc,
                    magnetic_moments=magnetic_moments,
                    setup_params=self.Setup_params,
                    xc_calc=self.XC_calc,
                    exx_fraction=self.XC_exx_fraction,
                    omega=self.XC_omega,
                    occupation=self.Occupation,
                    parallel_cores=self.parallel_cores,
                    executable='pw.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    f"QE ground-state calculation failed: {exc}"
                )
                raise

            result = workflow[
                'result'
            ]

            parprint(
                "QE ground state calculation finished."
            )

        parprint(
            "Total energy: "
            f"{result['total_energy_ev']:.8f} eV "
            f"({result['total_energy_ry']:.8f} Ry)"
        )

        write_cif(
            final_structure_file,
            self.bulk_configuration,
        )

    def _groundcalc_gpaw(self):
        """
        This method performs ground state calculations for the given structure using various settings
        and parameters specified in the configuration file. It handles different XC functionals,
        spin calculations, and geometry optimizations. The results are saved in appropriate files,
        including the final configuration as a CIF file and the ground state results in a GPW file.
        """

        # -------------------------------------------------------------
        # GROUND STATE - GPAW
        # -------------------------------------------------------------
        
        # Resolve XC functional string and PAW setups
        actual_xc, resolved_setups, is_libxc = self.engine.resolve_xc_and_setups(self.XC_calc, self.Setup_params)
        
        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )

        gpaw_state_file = Path(
            self.struct
            + '-GROUND-GPAW-Result-State.gpw'
        )

        if (
            not self.Ground_calc
            and gpaw_state_file.is_file()
        ):
            self._load_existing_final_structure()
        
        if self.Spin_calc:
            initial_moments = (
                resolve_initial_magnetic_moments(
                    atoms=self.bulk_configuration,
                    magmom_per_atom=self.Magmom_per_atom,
                    magmom_single_atom=(
                        self.Magmom_single_atom
                    ),
                )
            )

            self.bulk_configuration.set_initial_magnetic_moments(
                initial_moments
            )

        if self.Mode == 'PW':
            if self.config.Ground_calc == True:
                # PW Ground State Calculations
                parprint("Starting PW ground state calculation...")
                if self.Geo_optim and True in self.Relax_cell:
                    # Cell relaxation needs the stress tensor, which is not
                    # available for GLLBSC(M) nor for plane-wave hybrids.
                    if self.XC_calc in ['GLLBSC', 'GLLBSCM'] or self.engine.is_hybrid(self.XC_calc):
                        parprint(
                            "\033[91mERROR:\033[0m Cell relaxation can not be used "
                            "with "+self.XC_calc+" xc."
                        )
                        parprint(
                            "Disable Relax_cell, or optimize the structure with PBE "
                            "and use the resulting CIF for this calculation."
                        )
                        parprint("Exiting...")
                        sys.exit(1)
                if self.engine.is_hybrid(actual_xc):
                    parprint('Starting Hybrid XC calculations...')
                    calc = self.engine.create_hybrid_pw_ground_calc(
                        cutoff=self.Cut_off_energy,
                        xc_calc=self.XC_calc,
                        exx_fraction=self.XC_exx_fraction,
                        omega=self.XC_omega,
                        backend=self.XC_backend,
                        mixer=self.Mixer_type,
                        charge=self.Total_charge,
                        spinpol=self.Spin_calc,
                        txt=self.struct+f'-GROUND-{self.Engine}-Log-SCF.txt',
                        convergence=self.Ground_convergence,
                        occupations=self.Occupation,
                        kpoint_density=self.Ground_kpts_density,
                        kpoint_size=(
                            self.Ground_kpts_x,
                            self.Ground_kpts_y,
                            self.Ground_kpts_z,
                        ),
                        gamma=ground_gamma,
                        nbands=self.Ground_num_of_bands,
                    )
                else:
                    parprint(f'Starting calculations with {actual_xc}...')
                    # Fix the spacegroup in the geometric optimization if wanted
                    if self.Geo_optim and self.Fix_symmetry:
                        self.bulk_configuration.set_constraint(
                            FixSymmetry(self.bulk_configuration)
                        )

                    calc = self.engine.create_regular_pw_ground_calc(
                        cutoff=self.Cut_off_energy,
                        xc=actual_xc,
                        setups=resolved_setups,
                        parallel={'domain': world.size},
                        mixer=self.Mixer_type,
                        charge=self.Total_charge,
                        spinpol=self.Spin_calc,
                        txt=self.struct+f'-GROUND-{self.Engine}-Log-SCF.txt',
                        convergence=self.Ground_convergence,
                        occupations=self.Occupation,
                        kpoint_density=self.Ground_kpts_density,
                        kpoint_size=(
                            self.Ground_kpts_x,
                            self.Ground_kpts_y,
                            self.Ground_kpts_z,
                        ),
                        gamma=ground_gamma,
                        nbands=self.Ground_num_of_bands,
                    )
                # Wrapping for vdW
                if hasattr(self.config, 'vdW_calc') and self.config.vdW_calc.upper() == 'D3':
                    from ase.calculators.dftd3 import DFTD3
                    gpaw_calc = calc
                    calc = DFTD3(dft=gpaw_calc, xc=self.config.XC_calc)
                    calc.write = gpaw_calc.write

                    if hasattr(gpaw_calc, 'get_fermi_level'): calc.get_fermi_level = gpaw_calc.get_fermi_level
                    parprint("Applying Grimme DFT-D3 via native ASE Wrapper on legacy GPAW...")
                self.bulk_configuration.calc = calc
                if self.Geo_optim == True:
                    if True in self.Relax_cell:
                        if self.Hydrostatic_pressure > 0.0:
                            uf = FrechetCellFilter(self.bulk_configuration, mask=self.Relax_cell, hydrostatic_strain=True, scalar_pressure=self.Hydrostatic_pressure)
                        else:
                            uf = FrechetCellFilter(self.bulk_configuration, mask=self.Relax_cell)
                        # Optimizer Selection
                        if self.Optimizer == 'FIRE':
                            from ase.optimize.fire import FIRE
                            relax = FIRE(uf, maxstep=self.Max_step, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        elif  self.Optimizer == 'LBFGS':
                            from ase.optimize.lbfgs import LBFGS
                            relax = LBFGS(uf, maxstep=self.Max_step, alpha=self.Alpha, damping=self.Damping, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        elif  self.Optimizer == 'GPMin':
                            from ase.optimize import GPMin
                            relax = GPMin(uf, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        else:
                            relax = QuasiNewton(uf, maxstep=self.Max_step, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                    else:
                        # Optimizer Selection
                        if self.Optimizer == 'FIRE':
                            from ase.optimize.fire import FIRE
                            relax = FIRE(self.bulk_configuration, maxstep=self.Max_step, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        elif  self.Optimizer == 'LBFGS':
                            from ase.optimize.lbfgs import LBFGS
                            relax = LBFGS(self.bulk_configuration, maxstep=self.Max_step, alpha=self.Alpha, damping=self.Damping, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        elif  self.Optimizer == 'GPMin':
                            from ase.optimize import GPMin
                            relax = GPMin(self.bulk_configuration, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        else:
                            relax = QuasiNewton(self.bulk_configuration, maxstep=self.Max_step, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                    relax.run(fmax=self.Max_F_tolerance)  # Consider tighter fmax!

                else:
                    self.bulk_configuration.set_calculator(calc)
                    self.bulk_configuration.get_potential_energy()

                # Store the ground-state Fermi level so that DOS/band methods
                # have a valid energy reference, including for hybrids where
                # get_fermi_level() can not be read back from the NSCF runs.
                try:
                    self.Ground_fermi_level = calc.get_fermi_level()
                except Exception:
                    self.Ground_fermi_level = None

                calc.write(self.struct+'-GROUND-GPAW-Result-State.gpw', mode="all")

                # Writes final configuration as CIF file
                write_cif(self.struct+f'-GROUND-{self.Engine}-Result-Final.cif', self.bulk_configuration)
            else:
                parprint("Passing PW ground state calculation...")
                # Control the ground state GPW file
                if not os.path.exists(self.struct+'-GROUND-GPAW-Result-State.gpw'):
                    parprint('\033[91mERROR:\033[0m'+self.struct+'-GROUND-GPAW-Result-State.gpw file can not be found. It is needed in other calculations. Firstly, finish the ground state calculation. You must have \033[95mGround_calc = True\033[0m line in your input file. Exiting.')
                    sys.exit(1)

            # A little clean-up
            if hasattr(self.config, 'vdW_calc') and self.config.vdW_calc.upper() == 'D3':
                clean_path = os.path.join(self.struct, "../..")
                clean_path = os.path.normpath(clean_path)
                # Erase the dftd3 temp files from input_dir
                temp_files = (glob.glob(os.path.join(clean_path, "dftd3_*")) +
                            glob.glob(os.path.join(clean_path, "ase_dftd3*")) +
                            [os.path.join(clean_path, ".dftd3_info")])

                for f in temp_files:
                    try:
                        if os.path.exists(f):
                            os.remove(f)
                    except Exception:
                        pass

        elif self.Mode == 'LCAO':
            if self.Ground_calc == True:
                parprint("Starting LCAO ground state calculation...")
                # Fix the spacegroup in the geometric optimization if wanted
                if self.Geo_optim and self.Fix_symmetry:
                    self.bulk_configuration.set_constraint(
                        FixSymmetry(self.bulk_configuration)
                    )

                calc = self.engine.create_lcao_ground_calc(
                    setups=self.Setup_params,
                    parallel={'domain': world.size},
                    mixer=self.Mixer_type,
                    charge=self.Total_charge,
                    spinpol=self.Spin_calc,
                    txt=self.struct+f'-GROUND-{self.Engine}-Log-SCF.txt',
                    convergence=self.Ground_convergence,
                    occupations=self.Occupation,
                    kpoint_density=self.Ground_kpts_density,
                    kpoint_size=(
                        self.Ground_kpts_x,
                        self.Ground_kpts_y,
                        self.Ground_kpts_z,
                    ),
                    gamma=ground_gamma,
                    nbands=self.Ground_num_of_bands,
                    grid_spacing=self.Ground_gpts_density,
                    grid_size=(
                        self.Ground_gpts_x,
                        self.Ground_gpts_y,
                        self.Ground_gpts_z,
                    ),
                )

                # Wrapping for vdW
                if hasattr(self.config, 'vdW_calc') and self.config.vdW_calc.upper() == 'D3':
                    from ase.calculators.dftd3 import DFTD3
                    gpaw_calc = calc
                    
                    # wrap it
                    calc = DFTD3(dft=gpaw_calc)
                    calc.write = gpaw_calc.write 
                    if hasattr(gpaw_calc, 'get_fermi_level'): calc.get_fermi_level = gpaw_calc.get_fermi_level
                    parprint("Applying Grimme DFT-D3 via native ASE Wrapper on legacy GPAW...")
                self.bulk_configuration.calc = calc
                if self.Geo_optim == True:
                    if True in self.Relax_cell:
                        #uf = FrechetCellFilter(self.bulk_configuration, mask=self.Relax_cell)
                        #relax = LBFGS(uf, maxstep=self.Max_step, alpha=self.Alpha, damping=self.Damping, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        parprint('\033[91mERROR:\033[0mModifying supercell and atom positions with a filter (Relax_cell keyword) is not implemented in LCAO mode.')
                        sys.exit(1)
                    else:
                        # Optimizer Selection
                        if self.Optimizer == 'FIRE':
                            from ase.optimize.fire import FIRE
                            relax = FIRE(self.bulk_configuration, maxstep=self.Max_step, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        elif self.Optimizer == 'LBFGS':
                            from ase.optimize.lbfgs import LBFGS
                            relax = LBFGS(self.bulk_configuration, maxstep=self.Max_step, alpha=self.Alpha, damping=self.Damping, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        elif self.Optimizer == 'GPMin':
                            from ase.optimize import GPMin
                            relax = GPMin(self.bulk_configuration, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                        else:
                            relax = QuasiNewton(self.bulk_configuration, maxstep=self.Max_step, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                    relax.run(fmax=self.Max_F_tolerance)  # Consider tighter fmax!
                else:
                    self.bulk_configuration.set_calculator(calc)
                    self.bulk_configuration.get_potential_energy()
                #relax = LBFGS(self.bulk_configuration, maxstep=self.Max_step, alpha=self.Alpha, damping=self.Damping, trajectory=self.struct+'-GROUND-GPAW-Result-Trajectory.traj')
                #relax.run(fmax=self.Max_F_tolerance)  # Consider much tighter fmax!
                #self.bulk_configuration.get_potential_energy()
                
                calc.write(self.struct+'-GROUND-GPAW-Result-State.gpw', mode="all")

                # Writes final configuration as CIF file
                write_cif(self.struct+f'-GROUND-{self.Engine}-Result-Final.cif', self.bulk_configuration)
                # Print final spacegroup information
                parprint("Final Spacegroup:",get_spacegroup(self.bulk_configuration, symprec=1e-2))
            else:
                parprint("Passing LCAO ground state calculation...")
                # Control the ground state GPW file
                if not os.path.exists(self.struct+'-GROUND-GPAW-Result-State.gpw'):
                    parprint('\033[91mERROR:\033[0m'+self.struct+'-GROUND-GPAW-Result-State.gpw file can not be found. It is needed in other calculations. Firstly, finish the ground state calculation. You must have \033[95mGround_calc = True\033[0m line in your input file. Exiting.')
                    sys.exit(1)

        elif self.Mode == 'FD':
            parprint("\033[91mERROR:\033[0mFD mode is not implemented in Nanoworks yet...")
            sys.exit(1)
        else:
            parprint("\033[91mERROR:\033[0mPlease enter correct mode information.")
            sys.exit(1)

    def elasticcalc(self, drawfigs=False, strain_n=5, strain_mag=0.01, thickness=None):
        """
        Calculate the full elastic constant tensor and derived moduli.
        - strain_n: Number of strain points (including zero) for each independent strain mode.
        - strain_mag: Maximum strain magnitude (fractional, e.g., 0.01 for 1% strain).
        - thickness: Effective thickness for 2D materials (Angstrom).
        """
        
        from elastic import get_elastic_tensor, get_elementary_deformations
        
        # -------------------------------------------------------------
        # ELASTIC CALCULATION
        # -------------------------------------------------------------
        
        # Start Elastic calc
        time151 = time.time()

        elastic_xc, resolved_setups, elastic_parallel, hybrid = (
            self.engine.resolve_elastic_settings(
                xc_calc=self.XC_calc,
                setups=self.Setup_params,
                world_size=world.size,
                exx_fraction=self.XC_exx_fraction,
                omega=self.XC_omega,
                backend=self.XC_backend,
            )
        )
        
        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )

        elastic_kpoint_density, elastic_kpoint_size, elastic_gamma = (
            resolve_stage_kpoint_settings(
                stage_density=self.Elastic_kpts_density,
                stage_size=(
                    self.Elastic_kpts_x,
                    self.Elastic_kpts_y,
                    self.Elastic_kpts_z,
                ),
                stage_gamma=self.Elastic_gamma,
                ground_density=self.Ground_kpts_density,
                ground_size=(
                    self.Ground_kpts_x,
                    self.Ground_kpts_y,
                    self.Ground_kpts_z,
                ),
                ground_gamma=ground_gamma,
            )
        )

        # Elastic constants rely on the stress tensor. Plane-wave hybrid
        # stress is not reliably available in GPAW, so retain the existing
        # warning when a hybrid functional is requested.
        if hybrid:
            parprint(
                "\033[93mWARNING:\033[0m Elastic constants with hybrid ("
                + self.XC_calc
                + ") XC rely on stress/forces that are not reliable "
                  "in plane-wave GPAW."
            )
            parprint(
                "It is recommended to compute elastic properties with PBE "
                "and use hybrids only for the electronic structure."
            )

        def make_elastic_calc():
            return self.engine.create_elastic_calc(
                cutoff=self.config.Cut_off_energy,
                xc=elastic_xc,
                setups=resolved_setups,
                parallel=elastic_parallel,
                spinpol=self.config.Spin_calc,
                kpoint_density=elastic_kpoint_density,
                kpoint_size=elastic_kpoint_size,
                gamma=elastic_gamma,
                mixer=self.config.Mixer_type,
                txt=self.struct
                    + '-ELASTIC-GPAW-Log-Elastic-deformations.txt',
                charge=self.config.Total_charge,
                convergence=self.config.Ground_convergence,
                occupations=self.config.Occupation,
                hybrid=hybrid,
            )
        
        # Load the optimized (reference) structure
        bulk_atoms = self.bulk_configuration
        ref_calc = self.engine.load_gpaw_calc(
            self.struct + '-GROUND-GPAW-Result-State.gpw',
            hybrid=hybrid,
        )
        bulk_atoms.set_calculator(ref_calc)
        parprint('Optimized (reference) structure is loaded.')
        try:
            ref_stress = bulk_atoms.get_stress(voigt=False)
        except Exception as e:
            raise RuntimeError(f"ERROR: Could not compute reference stress: {e}")
        parprint(f"Reference Stress Tensor (GPa):\n{ref_stress / GPa}")

        # --- Define the six independent strain matrices (Voigt components) ---
        strain_matrices = [
            np.array([[1, 0, 0],
                      [0, 0, 0],
                      [0, 0, 0]]),  # ε_xx
            np.array([[0, 0, 0],
                      [0, 1, 0],
                      [0, 0, 0]]),  # ε_yy
            np.array([[0, 0, 0],
                      [0, 0, 0],
                      [0, 0, 1]]),  # ε_zz
            np.array([[0, 1, 0],
                      [1, 0, 0],
                      [0, 0, 0]]),  # ε_xy
            np.array([[0, 0, 1],
                      [0, 0, 0],
                      [1, 0, 0]]),  # ε_xz
            np.array([[0, 0, 0],
                      [0, 0, 1],
                      [0, 1, 0]])   # ε_yz
        ]

        # Define names for each strain mode (following Voigt notation)
        strain_names = ['ε_xx', 'ε_yy', 'ε_zz', 'ε_xy', 'ε_xz', 'ε_yz']
        # --- Cache file for deformed systems ---
        cache_file = self.struct + '-ELASTIC-GPAW-Result-Elastic-deformations.traj'
        if os.path.exists(cache_file):
            parprint("Loading deformed systems from cache.")
            systems = read(cache_file, index=':')
        else:
            systems = []
            # For each strain mode, sample strain_n values linearly from -strain_mag to +strain_mag.
            for mode, strain in enumerate(strain_matrices):
                mode_label = strain_names[mode]
                strain_values = np.linspace(-strain_mag, strain_mag, strain_n)
                for eps in strain_values:
                    # If eps is effectively zero, use the reference ground state.
                    if abs(eps) < 1e-6:
                        parprint(f"Mode {mode_label}: Using ground-state stress for strain {eps:.4f}")
                        ref_atoms = bulk_atoms.copy()
                        # Use the previously computed reference values:
                        ref_atoms.set_calculator(SinglePointCalculator(ref_atoms, 
                                                energy=bulk_atoms.get_potential_energy(), 
                                                stress=ref_stress, 
                                                forces=bulk_atoms.get_forces()))
                        systems.append(ref_atoms)
                        continue
                    success = False
                    # Try a series of fallback factors to help convergence.
                    for factor in [1.0, 0.5, 0.1]:
                        current_eps = eps * factor
                        deformed_atoms = bulk_atoms.copy()
                        strain_tensor = np.eye(3) + current_eps * strain
                        deformed_atoms.set_cell(deformed_atoms.cell @ strain_tensor, scale_atoms=True)
            
                        # Attach a new GPAW calculator for the deformed structure using PBE.
                        deformed_atoms.set_calculator(
                            make_elastic_calc()
                        )
                        # Wrapping for vdW
                        if hasattr(self.config, 'vdW_calc') and self.config.vdW_calc.upper() == 'D3':
                            from ase.calculators.dftd3 import DFTD3
                            gpaw_calc = deformed_atoms
                            
                            # wrap it
                            calc = DFTD3(dft=gpaw_calc)
                            calc.write = gpaw_calc.write 
                            if hasattr(gpaw_calc, 'get_fermi_level'): calc.get_fermi_level = gpaw_calc.get_fermi_level
                            parprint("Applying Grimme DFT-D3 via native ASE Wrapper on legacy GPAW...")
                        
                        try:
                            deformed_atoms.get_potential_energy()
                            # Force stress calculation.
                            deformed_atoms.calc.calculate(deformed_atoms, properties=['stress'])
                            stress_tensor = deformed_atoms.get_stress(voigt=False)
                            parprint(f"Mode {mode_label}: Computed stress for strain {eps:.4f} (factor {factor}):\n{stress_tensor}")
                            systems.append(deformed_atoms)
                            success = True
                            break
                        except Exception as e:
                            parprint(f"ERROR: Mode {mode_label}: Failed to compute stress for strain {eps:.4f} (factor {factor}): {e}")
                if not success:
                    parprint(f"WARNING: Mode {mode_label}: Skipping strain {eps:.4f} for current mode.")
            # Write deformed systems to cache file so that next time they can be loaded
            write(cache_file, systems)
            parprint(f"Deformed systems saved to cache file: {cache_file}")

        expected = len(strain_matrices) * strain_n
        if len(systems) < 0.5 * expected:
            raise RuntimeError(f"ERROR: Not enough valid deformations (only {len(systems)} out of {expected}) to compute the elastic tensor.")

        # --- Compute the elastic tensor using the deformed systems ---
        Cij, fit_info = get_elastic_tensor(bulk_atoms, systems)
        Cij = np.array(Cij)
        parprint(f"Cij raw shape: {Cij.shape}")
        if Cij.size < 36:
            Cij = reconstruct_full_tensor(Cij, bulk_atoms)
        elif Cij.size == 36:
            Cij = Cij.reshape((6,6))
        else:
            raise ValueError(f"ERROR: Unexpected Cij size: {Cij.size}. Cannot reshape.")
        Cij_GPa = Cij / GPa

        # --- 2D vs. 3D handling (unchanged) ---
        cell_lengths = bulk_atoms.cell.lengths()  # [a, b, c]
        is2D = False
        if thickness or (cell_lengths[2] > 3 * max(cell_lengths[0], cell_lengths[1]) and abs(Cij_GPa[2,2]) < 1e-2):
            is2D = True
            parprint("2D system is detected.")
            t_eff = thickness if thickness else cell_lengths[2]
            t_eff_m = t_eff * 1e-10  # convert Angstrom to meter
            C2D = Cij_GPa * (t_eff_m * 1e9)  # Convert GPa to N/m
            C11 = C2D[0, 0]; C22 = C2D[1, 1]; C12 = C2D[0, 1]
            E2D_x = (C11**2 - C12**2) / C11 if C11 > 0 else 0.0
            E2D = (E2D_x + ((C22**2 - C12**2) / C22 if C22 > 0 else 0.0)) / 2
            nu2D = C12 / C11 if C11 != 0 else 0.0
        else:
            parprint("2D system is not detected.")
            Bv = (Cij_GPa[0,0] + Cij_GPa[1,1] + Cij_GPa[2,2] +
                  2*(Cij_GPa[0,1] + Cij_GPa[0,2] + Cij_GPa[1,2])) / 9.0
            Gv = ((Cij_GPa[0,0] + Cij_GPa[1,1] + Cij_GPa[2,2]) -
                  (Cij_GPa[0,1] + Cij_GPa[0,2] + Cij_GPa[1,2]) +
                  3*(Cij_GPa[3,3] + Cij_GPa[4,4] + Cij_GPa[5,5])) / 15.0
            try:
                Sij = np.linalg.inv(Cij_GPa)
            except np.linalg.LinAlgError:
                Sij = np.linalg.pinv(Cij_GPa)
            Br = 1.0 / (Sij[0,0] + Sij[1,1] + Sij[2,2] +
                        2*(Sij[0,1] + Sij[0,2] + Sij[1,2]))
            Gr = 15.0 / (4*(Sij[0,0] + Sij[1,1] + Sij[2,2]) -
                         4*(Sij[0,1] + Sij[0,2] + Sij[1,2]) +
                         3*(Sij[3,3] + Sij[4,4] + Sij[5,5]))
            B_hill = 0.5 * (Bv + Br)
            G_hill = 0.5 * (Gv + Gr)
            E_hill = (9 * B_hill * G_hill) / (3 * B_hill + G_hill)
            nu_hill = (3 * B_hill - 2 * G_hill) / (2 * (3 * B_hill + G_hill))
        
        with paropen(self.struct + '-ELASTIC-GPAW-Result-Elastic-AllResults.txt', 'w') as fd:
            print("Elastic tensor Cij (GPa):", file=fd)
            print(np.array2string(Cij_GPa, precision=2, floatmode='fixed'), file=fd)
            if is2D:
                print(f"Detected 2D material. Effective thickness = {t_eff:.2f} Å.", file=fd)
                print("In-plane elastic stiffness (N/m):", file=fd)
                print(np.array2string(C2D, precision=2, floatmode='fixed'), file=fd)
                print(f"2D (in-plane) Young's modulus: {E2D:.2f} N/m", file=fd)
                print(f"2D Poisson's ratio (in-plane): {nu2D:.3f}", file=fd)
            else:
                print(f"Bulk modulus (Hill avg): {B_hill:.2f} GPa", file=fd)
                print(f"Shear modulus (Hill avg): {G_hill:.2f} GPa", file=fd)
                print(f"Young's modulus (Hill avg): {E_hill:.2f} GPa", file=fd)
                print(f"Poisson's ratio (Hill avg): {nu_hill:.3f}", file=fd)
        
        # Finish elastic timing
        time152 = time.time()

        # Write timings of calculation
        with paropen(self.struct+f'-TIMINGS-{self.Engine}-Log-Timings.txt', 'a') as f1:
            print('Elastic Calculation: ', round((time152-time151),2), end="\n", file=f1)


    def hybrid_fermi_level(self, calc=None):
        """
        Return a sensible Fermi-level reference (eV) for hybrid calculations.

        Tries, in order: the value stored from the ground-state run, then the
        Fermi level of a freshly read ground-state .gpw, then the Fermi level
        of the calculator passed in. Falls back to 0.0 eV only if nothing is
        available, instead of unconditionally hard-coding 0.0 eV as before.
        """
        if self.Ground_fermi_level is not None and not np.isnan(self.Ground_fermi_level):
            return self.Ground_fermi_level
        # Try reading the converged ground-state .gpw written by groundcalc().
        try:
            ref_calc = self.engine.load_gpaw_calc(
                self.struct+'-GROUND-GPAW-Result-State.gpw',
                hybrid=True,
            )
            ef = ref_calc.get_fermi_level()
            if ef is not None and not np.isnan(ef):
                self.Ground_fermi_level = ef
                return ef
        except Exception as e:
            parprint(f"\033[93mWARNING:\033[0m Could not read Fermi level from ground-state file: {e}")
        if calc is not None:
            try:
                ef = calc.get_fermi_level()
                if ef is not None and not np.isnan(ef):
                    return ef
            except Exception as e:
                parprint("\033[93mWARNING:\033[0m Failed to read Fermi level from provided calculator:", e)
        parprint("\033[93mWARNING:\033[0m Could not determine the Fermi level for the hybrid calculation; using 0.0 eV as reference.")
        return 0.0

    def doscalc(self):
        """Run the DOS workflow using the selected DFT engine."""
        if self.Engine == 'GPAW':
            return self._doscalc_gpaw()

        if self.Engine == 'QE':
            return self._doscalc_qe()

        raise ValueError(
            f"Unsupported DFT engine: {self.Engine}"
        )
        
    def _doscalc_gpaw(self):
        """
        This method performs density of states (DOS) calculations for the given structure using
        the ground state results. It computes the DOS for various energy levels and saves the
        results in appropriate files for further electronic analysis and visualization.
        """

        # -------------------------------------------------------------
        # DOS CALCULATION
        # -------------------------------------------------------------

        time21 = time.time()
        parprint("Starting DOS calculation...")

        hybrid = self.engine.is_hybrid(self.XC_calc)
        
        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )

        if hybrid:
            # Hybrids can NOT use fixed_density(): the exchange operator
            # depends on the occupied orbitals, not only on the density. We
            # therefore read the eigenvalues stored in the converged ground
            # state and reference them to the ground-state Fermi level.
            parprint(
                'Passing DOS NSCF calculations '
                '(using ground-state eigenvalues for hybrid)...'
            )

        dos_kpoint_density, dos_kpoint_size, dos_gamma = (
            resolve_stage_kpoint_settings(
                stage_density=self.DOS_kpts_density,
                stage_size=(
                    self.DOS_kpts_x,
                    self.DOS_kpts_y,
                    self.DOS_kpts_z,
                ),
                stage_gamma=self.DOS_gamma,
                ground_density=self.Ground_kpts_density,
                ground_size=(
                    self.Ground_kpts_x,
                    self.Ground_kpts_y,
                    self.Ground_kpts_z,
                ),
                ground_gamma=ground_gamma,
            )
        )

        dos_occupation = resolve_stage_occupation(
            self.DOS_occupation,
            self.Occupation,
        )

        calc = self.engine.prepare_dos_calc(
            filename=self.struct+'-GROUND-GPAW-Result-State.gpw',
            hybrid=hybrid,
            txt=self.struct+f'-DOS-{self.Engine}-Log-DOS.txt',
            convergence=self.DOS_convergence,
            occupations=dos_occupation,
            kpoint_density=dos_kpoint_density,
            kpoint_size=dos_kpoint_size,
            gamma=dos_gamma,
            nbands=self.DOS_num_of_bands,
        )

        if hybrid:
            ef = self.hybrid_fermi_level(calc)
        else:
            ef = calc.get_fermi_level()

        # Keep GPAW CSV energies on the same E - Ef axis as the DOS graph and QE exports.
        chem_sym = self.bulk_configuration.get_chemical_symbols()

        if self.SOC_calc:
            parprint("Calculating Total DOS with Spin-Orbit Coupling...")
            from gpaw.spinorbit import soc_eigenstates

            soc = soc_eigenstates(calc)

            # Take eigenvalues and use only real parts
            soc_evals_raw = soc.eigenvalues()
            soc_evals = np.real(np.array(soc_evals_raw).flatten())

            # NaN error protection for Fermi level
            ef_safe = ef if not np.isnan(ef) else 0.0
            soc_evals = soc_evals - ef_safe

            energies = np.linspace(self.Energy_min, self.Energy_max, self.DOS_npoints)
            dos = np.zeros_like(energies)

            # Divide by zero protection
            smearing = getattr(self, 'DOS_width', 0.1)
            smearing = max(smearing, 0.01)

            # Gaussian smearing
            for e in soc_evals:
                if not np.isnan(e): # Only valid numbers
                    dos += np.exp(-((energies - e) / smearing)**2) / (smearing * np.sqrt(np.pi))

            # Normalize with dividing to k-points number
            try:
                nkpts = max(1, len(calc.get_ibz_k_points()))
            except:
                nkpts = max(1, len(soc_evals_raw))
            dos /= nkpts

            # If data is empty...
            if np.isnan(dos).all() or np.max(dos) == 0.0:
                parprint("WARNING: Computed SOC-DOS array is empty or NaN. Matplotlib crash prevented.")
                dos[:] = 0.0

            with paropen(self.struct+'-DOS-GPAW-Result-DOS-SOC.csv', "w") as fd:
                for en, d in zip(energies, dos):
                    print(f"{en}, {d}", file=fd)

            # Use explicit figure creation to avoid memory overlap
            if world.rank == 0:
                fig_soc, ax_soc = plt.subplots(figsize=(8, 6))
                ax_soc.plot(energies, dos, 'purple', label='DOS (SOC)')
                ax_soc.set_xlabel(self._t("fig_dos_xlabel"))
                ax_soc.set_ylabel(self._t("fig_dos_ylabel"))
                ax_soc.set_xlim(self.Energy_min, self.Energy_max)
                ax_soc.legend()

                # Draw scaling of DOS
                if np.max(dos) > 0.0 and not np.isnan(dos).all():
                    autoscale_y(ax_soc)
                else:
                    ax_soc.set_ylim(0, 1)

                plt.tight_layout()
                plt.savefig(self.struct+'-DOS-GPAW-Graph-DOS-SOC.png', dpi=300)
                plt.close(fig_soc)

            parprint("Spin-Orbit Total DOS is computed. PDOS is not computed for SOC.")
            return

        if self.Spin_calc == True:
            # ==========================================
            # SPIN DOWN CALCULATIONS
            # ==========================================
            parprint("Calculating and saving Raw PDOS for spin down...")
            
            from gpaw.dos import DOSCalculator
            
            rawdos = DOSCalculator.from_calculator(filename=self.struct+'-GROUND-GPAW-Result-State.gpw', soc=False, theta=0.0, phi=0.0, shift_fermi_level=False)
            energies = rawdos.get_energies(npoints=self.DOS_npoints)
            shifted_energies = energies - ef

            # Weights initialization
            pdossweightsdown = [0.0] * self.DOS_npoints
            pdospweightsdown = [0.0] * self.DOS_npoints
            pdospxweightsdown = [0.0] * self.DOS_npoints
            pdospyweightsdown = [0.0] * self.DOS_npoints
            pdospzweightsdown = [0.0] * self.DOS_npoints
            pdosdweightsdown = [0.0] * self.DOS_npoints
            pdosdxyweightsdown = [0.0] * self.DOS_npoints
            pdosdyzweightsdown = [0.0] * self.DOS_npoints
            pdosd3z2_r2weightsdown = [0.0] * self.DOS_npoints
            pdosdzxweightsdown = [0.0] * self.DOS_npoints
            pdosdx2_y2weightsdown = [0.0] * self.DOS_npoints
            pdosfweightsdown = [0.0] * self.DOS_npoints
            totaldosweightsdown = [0.0] * self.DOS_npoints

            # Writing RawPDOS
            with paropen(self.struct+'-DOS-GPAW-Result-Raw-PDOS-Each-Atom-Down.csv', "w") as fd:
                print("Energy, s-total, p-total, pz, px, py, d-total, d3z2_r2, dzx, dyz, dx2_y2, dxy, f-total, TOTAL", file=fd)
                for j in range(0, self.bulk_configuration.get_global_number_of_atoms()):
                    print(f"Atom no: {j+1}, Atom Symbol: {chem_sym[j]} --------------------", file=fd)
                    pdoss = rawdos.raw_pdos(energies, a=j, l=0, m=None, spin=0, width=self.DOS_width)
                    pdosp = rawdos.raw_pdos(energies, a=j, l=1, m=None, spin=0, width=self.DOS_width)
                    pdospx = rawdos.raw_pdos(energies, a=j, l=1, m=2, spin=0, width=self.DOS_width)
                    pdospy = rawdos.raw_pdos(energies, a=j, l=1, m=0, spin=0, width=self.DOS_width)
                    pdospz = rawdos.raw_pdos(energies, a=j, l=1, m=1, spin=0, width=self.DOS_width)
                    pdosd = rawdos.raw_pdos(energies, a=j, l=2, m=None, spin=0, width=self.DOS_width)
                    pdosdxy = rawdos.raw_pdos(energies, a=j, l=2, m=0, spin=0, width=self.DOS_width)
                    pdosdyz = rawdos.raw_pdos(energies, a=j, l=2, m=1, spin=0, width=self.DOS_width)
                    pdosd3z2_r2 = rawdos.raw_pdos(energies, a=j, l=2, m=2, spin=0, width=self.DOS_width)
                    pdosdzx = rawdos.raw_pdos(energies, a=j, l=2, m=3, spin=0, width=self.DOS_width)
                    pdosdx2_y2 = rawdos.raw_pdos(energies, a=j, l=2, m=4, spin=0, width=self.DOS_width)
                    pdosf = rawdos.raw_pdos(energies, a=j, l=3, m=None, spin=0, width=self.DOS_width)
                    dosspdf = pdoss + pdosp + pdosd + pdosf
                    pdossweightsdown = pdossweightsdown + pdoss
                    pdospweightsdown = pdospweightsdown + pdosp
                    pdospxweightsdown = pdospxweightsdown + pdospx
                    pdospyweightsdown = pdospyweightsdown + pdospy
                    pdospzweightsdown = pdospzweightsdown + pdospz
                    pdosdweightsdown = pdosdweightsdown + pdosd
                    pdosdxyweightsdown = pdosdxyweightsdown + pdosdxy
                    pdosdyzweightsdown = pdosdyzweightsdown + pdosdyz
                    pdosd3z2_r2weightsdown = pdosd3z2_r2weightsdown + pdosd3z2_r2
                    pdosdzxweightsdown = pdosdzxweightsdown + pdosdzx
                    pdosdx2_y2weightsdown = pdosdx2_y2weightsdown + pdosdx2_y2
                    pdosfweightsdown = pdosfweightsdown + pdosf
                    totaldosweightsdown = totaldosweightsdown + dosspdf

                    for x in zip(shifted_energies, pdoss, pdosp, pdospz, pdospx, pdospy, pdosd, pdosd3z2_r2, pdosdzx, pdosdyz, pdosdx2_y2, pdosdxy, pdosf, dosspdf):
                        print(*x, sep=", ", file=fd)

            # Writing DOS
            parprint("Saving DOS for spin down...")
            with paropen(self.struct+f'-DOS-{self.Engine}-Result-DOS-Down.csv', "w") as fd:
                for x in zip(shifted_energies, totaldosweightsdown):
                    print(*x, sep=", ", file=fd)

            # Writing PDOS
            parprint("Saving PDOS for spin down...")
            with paropen(self.struct+f'-DOS-{self.Engine}-Result-PDOS-Down.csv', "w") as fd:
                print("Energy, s-total, p-total, pz, px, py, d-total, d3z2_r2, dzx, dyz, dx2_y2, dxy, f-total, TOTAL", file=fd)
                for x in zip(shifted_energies, pdossweightsdown, pdospweightsdown, pdospzweightsdown, pdospxweightsdown, pdospyweightsdown, pdosdweightsdown, pdosd3z2_r2weightsdown, pdosdzxweightsdown, pdosdyzweightsdown, pdosdx2_y2weightsdown, pdosdxyweightsdown, pdosfweightsdown, totaldosweightsdown):
                    print(*x, sep=", ", file=fd)

            # ==========================================
            # SPIN UP CALCULATIONS
            # ==========================================
            parprint("Calculating and saving Raw PDOS for spin up...")
            rawdos = DOSCalculator.from_calculator(self.struct+'-GROUND-GPAW-Result-State.gpw', soc=False, theta=0.0, phi=0.0, shift_fermi_level=False)
            energies = rawdos.get_energies(npoints=self.DOS_npoints)
            shifted_energies = energies - ef

            # Weights initialization
            pdossweightsup = [0.0] * self.DOS_npoints
            pdospweightsup = [0.0] * self.DOS_npoints
            pdospxweightsup = [0.0] * self.DOS_npoints
            pdospyweightsup = [0.0] * self.DOS_npoints
            pdospzweightsup = [0.0] * self.DOS_npoints
            pdosdweightsup = [0.0] * self.DOS_npoints
            pdosdxyweightsup = [0.0] * self.DOS_npoints
            pdosdyzweightsup = [0.0] * self.DOS_npoints
            pdosd3z2_r2weightsup = [0.0] * self.DOS_npoints
            pdosdzxweightsup = [0.0] * self.DOS_npoints
            pdosdx2_y2weightsup = [0.0] * self.DOS_npoints
            pdosfweightsup = [0.0] * self.DOS_npoints
            totaldosweightsup = [0.0] * self.DOS_npoints

            # Writing RawPDOS
            with paropen(self.struct+'-DOS-GPAW-Result-Raw-PDOS-Each-Atom-Up.csv', "w") as fd:
                print("Energy, s-total, p-total, pz, px, py, d-total, d3z2_r2, dzx, dyz, dx2_y2, dxy, f-total, TOTAL", file=fd)
                for j in range(0, self.bulk_configuration.get_global_number_of_atoms()):
                    print(f"Atom no: {j+1}, Atom Symbol: {chem_sym[j]} --------------------", file=fd)
                    pdoss = rawdos.raw_pdos(energies, a=j, l=0, m=None, spin=1, width=self.DOS_width)
                    pdosp = rawdos.raw_pdos(energies, a=j, l=1, m=None, spin=1, width=self.DOS_width)
                    pdospx = rawdos.raw_pdos(energies, a=j, l=1, m=2, spin=1, width=self.DOS_width)
                    pdospy = rawdos.raw_pdos(energies, a=j, l=1, m=0, spin=1, width=self.DOS_width)
                    pdospz = rawdos.raw_pdos(energies, a=j, l=1, m=1, spin=1, width=self.DOS_width)
                    pdosd = rawdos.raw_pdos(energies, a=j, l=2, m=None, spin=1, width=self.DOS_width)
                    pdosdxy = rawdos.raw_pdos(energies, a=j, l=2, m=0, spin=1, width=self.DOS_width)
                    pdosdyz = rawdos.raw_pdos(energies, a=j, l=2, m=1, spin=1, width=self.DOS_width)
                    pdosd3z2_r2 = rawdos.raw_pdos(energies, a=j, l=2, m=2, spin=1, width=self.DOS_width)
                    pdosdzx = rawdos.raw_pdos(energies, a=j, l=2, m=3, spin=1, width=self.DOS_width)
                    pdosdx2_y2 = rawdos.raw_pdos(energies, a=j, l=2, m=4, spin=1, width=self.DOS_width)
                    pdosf = rawdos.raw_pdos(energies, a=j, l=3, m=None, spin=1, width=self.DOS_width)
                    dosspdf = pdoss + pdosp + pdosd + pdosf
                    pdossweightsup = pdossweightsup + pdoss
                    pdospweightsup = pdospweightsup + pdosp
                    pdospxweightsup = pdospxweightsup + pdospx
                    pdospyweightsup = pdospyweightsup + pdospy
                    pdospzweightsup = pdospzweightsup + pdospz
                    pdosdweightsup = pdosdweightsup + pdosd
                    pdosdxyweightsup = pdosdxyweightsup + pdosdxy
                    pdosdyzweightsup = pdosdyzweightsup + pdosdyz
                    pdosd3z2_r2weightsup = pdosd3z2_r2weightsup + pdosd3z2_r2
                    pdosdzxweightsup = pdosdzxweightsup + pdosdzx
                    pdosdx2_y2weightsup = pdosdx2_y2weightsup + pdosdx2_y2
                    pdosfweightsup = pdosfweightsup + pdosf
                    totaldosweightsup = totaldosweightsup + dosspdf

                    for x in zip(shifted_energies, pdoss, pdosp, pdospz, pdospx, pdospy, pdosd, pdosd3z2_r2, pdosdzx, pdosdyz, pdosdx2_y2, pdosdxy, pdosf, dosspdf):
                        print(*x, sep=", ", file=fd)

            # Writing DOS
            parprint("Saving DOS for spin up...")
            with paropen(self.struct+f'-DOS-{self.Engine}-Result-DOS-Up.csv', "w") as fd:
                for x in zip(shifted_energies, totaldosweightsup):
                    print(*x, sep=", ", file=fd)

            # Writing PDOS
            parprint("Saving PDOS for spin up...")
            with paropen(self.struct+f'-DOS-{self.Engine}-Result-PDOS-Up.csv', "w") as fd:
                print("Energy, s-total, p-total, pz, px, py, d-total, d3z2_r2, dzx, dyz, dx2_y2, dxy, f-total, TOTAL", file=fd)
                for x in zip(shifted_energies, pdossweightsup, pdospweightsup, pdospzweightsup, pdospxweightsup, pdospyweightsup, pdosdweightsup, pdosd3z2_r2weightsup, pdosdzxweightsup, pdosdyzweightsup, pdosdx2_y2weightsup, pdosdxyweightsup, pdosfweightsup, totaldosweightsup):
                    print(*x, sep=", ", file=fd)

        else:
            # ==========================================
            # NON-SPIN CALCULATIONS
            # ==========================================
            parprint("Calculating and saving Raw PDOS...")
            
            from gpaw.dos import DOSCalculator
            
            rawdos = DOSCalculator.from_calculator(self.struct+'-GROUND-GPAW-Result-State.gpw', soc=False, theta=0.0, phi=0.0, shift_fermi_level=False)
            energies = rawdos.get_energies(npoints=self.DOS_npoints)
            shifted_energies = energies - ef

            totaldosweights = [0.0] * self.DOS_npoints
            pdossweights = [0.0] * self.DOS_npoints
            pdospweights = [0.0] * self.DOS_npoints
            pdospxweights = [0.0] * self.DOS_npoints
            pdospyweights = [0.0] * self.DOS_npoints
            pdospzweights = [0.0] * self.DOS_npoints
            pdosdweights = [0.0] * self.DOS_npoints
            pdosdxyweights = [0.0] * self.DOS_npoints
            pdosdyzweights = [0.0] * self.DOS_npoints
            pdosd3z2_r2weights = [0.0] * self.DOS_npoints
            pdosdzxweights = [0.0] * self.DOS_npoints
            pdosdx2_y2weights = [0.0] * self.DOS_npoints
            pdosfweights = [0.0] * self.DOS_npoints

            # Writing RawPDOS
            with paropen(self.struct+'-DOS-GPAW-Result-Raw-PDOS-Each-Atom.csv', "w") as fd:
                print("Energy, s-total, p-total, pz, px, py, d-total, d3z2_r2, dzx, dyz, dx2_y2, dxy, f-total, TOTAL", file=fd)
                for j in range(0, self.bulk_configuration.get_global_number_of_atoms()):
                    print(f"Atom no: {j+1}, Atom Symbol: {chem_sym[j]} ----------------------------------------", file=fd)
                    pdoss = rawdos.raw_pdos(energies, a=j, l=0, m=None, spin=None, width=self.DOS_width)
                    pdosp = rawdos.raw_pdos(energies, a=j, l=1, m=None, spin=None, width=self.DOS_width)
                    pdospx = rawdos.raw_pdos(energies, a=j, l=1, m=2, spin=None, width=self.DOS_width)
                    pdospy = rawdos.raw_pdos(energies, a=j, l=1, m=0, spin=None, width=self.DOS_width)
                    pdospz = rawdos.raw_pdos(energies, a=j, l=1, m=1, spin=None, width=self.DOS_width)
                    pdosd = rawdos.raw_pdos(energies, a=j, l=2, m=None, spin=None, width=self.DOS_width)
                    pdosdxy = rawdos.raw_pdos(energies, a=j, l=2, m=0, spin=None, width=self.DOS_width)
                    pdosdyz = rawdos.raw_pdos(energies, a=j, l=2, m=1, spin=None, width=self.DOS_width)
                    pdosd3z2_r2 = rawdos.raw_pdos(energies, a=j, l=2, m=2, spin=None, width=self.DOS_width)
                    pdosdzx = rawdos.raw_pdos(energies, a=j, l=2, m=3, spin=None, width=self.DOS_width)
                    pdosdx2_y2 = rawdos.raw_pdos(energies, a=j, l=2, m=4, spin=None, width=self.DOS_width)
                    pdosf = rawdos.raw_pdos(energies, a=j, l=3, m=None, spin=None, width=self.DOS_width)

                    dosspdf = pdoss + pdosp + pdosd + pdosf

                    # Weights accumulation (CRITICAL BUG FIXED HERE)
                    pdossweights = pdossweights + pdoss
                    pdospweights = pdospweights + pdosp
                    pdospxweights = pdospxweights + pdospx
                    pdospyweights = pdospyweights + pdospy
                    pdospzweights = pdospzweights + pdospz
                    pdosdweights = pdosdweights + pdosd
                    pdosdxyweights = pdosdxyweights + pdosdxy
                    pdosdyzweights = pdosdyzweights + pdosdyz
                    pdosd3z2_r2weights = pdosd3z2_r2weights + pdosd3z2_r2
                    pdosdzxweights = pdosdzxweights + pdosdzx
                    pdosdx2_y2weights = pdosdx2_y2weights + pdosdx2_y2
                    pdosfweights = pdosfweights + pdosf
                    totaldosweights = totaldosweights + dosspdf

                    for x in zip(shifted_energies, pdoss, pdosp, pdospz, pdospx, pdospy, pdosd, pdosd3z2_r2, pdosdzx, pdosdyz, pdosdx2_y2, pdosdxy, pdosf, dosspdf):
                        print(*x, sep=", ", file=fd)

            # Writing DOS
            parprint("Saving DOS...")
            with paropen(self.struct+f'-DOS-{self.Engine}-Result-DOS.csv', "w") as fd:
                for x in zip(shifted_energies, totaldosweights):
                    print(*x, sep=", ", file=fd)

            # Writing PDOS
            parprint("Saving PDOS...")
            with paropen(self.struct+f'-DOS-{self.Engine}-Result-PDOS.csv', "w") as fd:
                print("Energy, s-total, p-total, pz, px, py, d-total, d3z2_r2, dzx, dyz, dx2_y2, dxy, f-total, TOTAL", file=fd)
                for x in zip(shifted_energies, pdossweights, pdospweights, pdospzweights, pdospxweights, pdospyweights, pdosdweights, pdosd3z2_r2weights, pdosdzxweights, pdosdyzweights, pdosdx2_y2weights, pdosdxyweights, pdosfweights, totaldosweights):
                    print(*x, sep=", ", file=fd)

        # Finish DOS calc
        time22 = time.time()

        # Write timings of calculation
        with paropen(self.struct+f'-TIMINGS-{self.Engine}-Log-Timings.txt', 'a') as f1:
            print(f'DOS calculation: {round((time22-time21),2)}', file=f1)

        if world.rank == 0:
            import pandas as pd
            fig, ax = plt.subplots(figsize=(8, 6))

            if self.Spin_calc == True:
                downf = pd.read_csv(self.struct+f'-DOS-{self.Engine}-Result-DOS-Down.csv', header=None)
                upf = pd.read_csv(self.struct+f'-DOS-{self.Engine}-Result-DOS-Up.csv', header=None)

                ax.plot(downf[0], -1.0*downf[1], 'r', linewidth=1.5, label='Spin Down')
                ax.plot(upf[0], upf[1], 'b', linewidth=1.5, label='Spin Up')

                ax.fill_between(downf[0], 0, -1.0*downf[1], facecolor='red', alpha=0.2)
                ax.fill_between(upf[0], 0, upf[1], facecolor='blue', alpha=0.2)

            else:
                dosf = pd.read_csv(self.struct+f'-DOS-{self.Engine}-Result-DOS.csv', header=None)

                ax.plot(dosf[0], dosf[1], 'b', linewidth=1.5)
                ax.fill_between(dosf[0], 0, dosf[1], facecolor='blue', alpha=0.2)

            # Axis labels (Using multi-language infrastructure)
            ax.set_xlabel(self._t("fig_dos_xlabel"))
            ax.set_ylabel(self._t("fig_dos_ylabel"))

            # Draw a vertical dashed black line exactly at 0 (Fermi Level)
            ax.axvline(x=0, color='k', linestyle='--', linewidth=1)
            ax.axhline(y=0, color='k', linewidth=0.8)

            # X limits should now be pure min and max since we shifted the data
            ax.set_xlim(self.Energy_min, self.Energy_max)

            autoscale_y(ax)

            # Layout and Saving
            plt.tight_layout()
            plt.savefig(self.struct+f'-DOS-{self.Engine}-Graph-DOS.png', dpi=300)

            # Clear memory to prevent interference with upcoming Band/Optical calculations
            plt.close(fig)

    def _doscalc_qe(self):
        """Run the DOS preparation workflow using Quantum ESPRESSO."""
        time21 = time.time()

        parprint(
            "Starting QE DOS calculation..."
        )

        if self.Mode != 'PW':
            parprint(
                "\033[91mERROR:\033[0m "
                "Quantum ESPRESSO DOS calculations "
                "support PW mode only."
            )
            sys.exit(1)

        if self.SOC_calc:
            parprint(
                "\033[91mERROR:\033[0m "
                "Quantum ESPRESSO SOC DOS calculations are not "
                "supported yet."
            )
            sys.exit(1)

        try:
            validated_xc = self.engine.validate_qe_xc(
                self.XC_calc,
                pseudo_xc='pbe',
                allow_hybrid=True,
            )
        except ValueError as exc:
            parprint(
                f"\033[91mERROR:\033[0m {exc}"
            )
            sys.exit(1)

        hybrid = str(
            validated_xc
        ).strip().lower() in {
            'hse06',
            'hse03',
            'pbe0',
        }

        ground_state_dir = Path(
            self.struct
            + '-GROUND-QE-Result-State'
        )

        if not self.engine.has_qe_state(
            ground_state_dir,
            prefix='nanoworks',
        ):
            parprint(
                "\033[91mERROR:\033[0m "
                + str(ground_state_dir)
                + " does not contain a valid QE ground-state result. "
                "Complete the ground-state calculation first."
            )
            sys.exit(1)

        state_dir = (
            Path(
                self.struct
                + '-DOS-QE-Result-State'
            )
            if hybrid
            else ground_state_dir
        )

        pseudo_dir = get_qe_pseudo_dir(
            relativistic='scalar',
        )

        try:
            pseudopotentials = (
                resolve_qe_pseudopotentials(
                    self.bulk_configuration,
                    relativistic='scalar',
                )
            )
        except (FileNotFoundError, RuntimeError) as exc:
            parprint(
                f"\033[91mERROR:\033[0m {exc}"
            )
            sys.exit(1)

        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )

        (
            dos_kpoint_density,
            dos_kpoint_size,
            dos_gamma,
        ) = resolve_stage_kpoint_settings(
            stage_density=self.DOS_kpts_density,
            stage_size=(
                self.DOS_kpts_x,
                self.DOS_kpts_y,
                self.DOS_kpts_z,
            ),
            stage_gamma=self.DOS_gamma,
            ground_density=self.Ground_kpts_density,
            ground_size=(
                self.Ground_kpts_x,
                self.Ground_kpts_y,
                self.Ground_kpts_z,
            ),
            ground_gamma=ground_gamma,
        )

        dos_occupation = resolve_stage_occupation(
            self.DOS_occupation,
            self.Occupation,
        )

        magnetic_moments = None

        if self.Spin_calc:
            magnetic_moments = (
                resolve_initial_magnetic_moments(
                    atoms=self.bulk_configuration,
                    magmom_per_atom=self.Magmom_per_atom,
                    magmom_single_atom=(
                        self.Magmom_single_atom
                    ),
                )
            )

        if hybrid:
            scf_input_file = Path(
                self.struct
                + '-DOS-QE-Input-Hybrid-SCF.in'
            )

            scf_output_file = Path(
                self.struct
                + '-DOS-QE-Log-Hybrid-SCF.txt'
            )
        else:
            input_file = Path(
                self.struct
                + '-DOS-QE-Input-NSCF.in'
            )

            output_file = Path(
                self.struct
                + '-DOS-QE-Log-NSCF.txt'
            )

        workflow = None
        result = None

        if not hybrid:
            try:
                workflow = self.engine.run_nscf(
                    atoms=self.bulk_configuration,
                    input_file=input_file,
                    output_file=output_file,
                    state_dir=state_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    cutoff_ev=self.Cut_off_energy,
                    kpoint_density=dos_kpoint_density,
                    kpoint_size=dos_kpoint_size,
                    gamma=dos_gamma,
                    total_charge=self.Total_charge,
                    nbands=self.DOS_num_of_bands,
                    spinpol=self.Spin_calc,
                    magnetic_moments=magnetic_moments,
                    setup_params=self.Setup_params,
                    xc_calc=self.XC_calc,
                    exx_fraction=self.XC_exx_fraction,
                    omega=self.XC_omega,
                    occupation=dos_occupation,
                    parallel_cores=self.parallel_cores,
                    executable='pw.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    f"QE DOS NSCF calculation failed: {exc}"
                )
                raise

            result = workflow['result']

            parprint(
                "QE DOS NSCF calculation finished."
            )

            if result['fermi_energy_ev'] is not None:
                parprint(
                    "NSCF Fermi energy: "
                    f"{result['fermi_energy_ev']:.8f} eV"
                )

        if self.DOS_npoints is None or int(self.DOS_npoints) < 2:
            raise ValueError(
                "DOS_npoints must be at least 2 for QE DOS calculations."
            )

        fermi_energy = None

        dos_emin_absolute = None
        dos_emax_absolute = None

        if not hybrid:
            fermi_energy = result[
                'fermi_energy_ev'
            ]

            if fermi_energy is None:
                raise RuntimeError(
                    "QE DOS requires a Fermi energy from the NSCF calculation."
                )

            dos_emin_absolute = (
                float(fermi_energy)
                + float(self.Energy_min)
            )

            dos_emax_absolute = (
                float(fermi_energy)
                + float(self.Energy_max)
            )

        delta_e = (
            float(self.Energy_max)
            - float(self.Energy_min)
        ) / (
            int(self.DOS_npoints) - 1
        )

        dos_input_file = Path(
            self.struct
            + '-DOS-QE-Input-DOS.in'
        )

        dos_output_file = Path(
            self.struct
            + f'-DOS-{self.Engine}-Log-DOS.txt'
        )

        dos_data_file = Path(
            self.struct
            + '-DOS-QE-Result-Raw-DOS.dat'
        )

        pdos_input_file = Path(
            self.struct
            + '-DOS-QE-Input-PDOS.in'
        )

        pdos_output_file = Path(
            self.struct
            + '-DOS-QE-Log-PDOS.txt'
        )

        pdos_prefix = Path(
            self.struct
            + '-DOS-QE-Result-Raw-PDOS'
        )

        qe_dos_occupation = (
            self.engine.resolve_qe_occupation(
                dos_occupation
            )
        )

        bz_sum = qe_dos_occupation[
            'occupations'
        ]

        if bz_sum not in {
            'tetrahedra',
            'tetrahedra_lin',
            'tetrahedra_opt',
        }:
            raise NotImplementedError(
                "Quantum ESPRESSO DOS currently supports "
                "tetrahedra occupations only in Nanoworks."
            )

        if hybrid:
            parprint(
                "Starting QE hybrid SCF, DOS and PDOS workflow..."
            )

            try:
                hybrid_workflow = self.engine.run_hybrid_dos(
                    atoms=self.bulk_configuration,
                    scf_input_file=scf_input_file,
                    scf_output_file=scf_output_file,
                    dos_input_file=dos_input_file,
                    dos_output_file=dos_output_file,
                    dos_file=dos_data_file,
                    pdos_input_file=pdos_input_file,
                    pdos_output_file=pdos_output_file,
                    pdos_prefix=pdos_prefix,
                    state_dir=state_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    cutoff_ev=self.Cut_off_energy,
                    kpoint_density=dos_kpoint_density,
                    kpoint_size=dos_kpoint_size,
                    gamma=dos_gamma,
                    total_charge=self.Total_charge,
                    nbands=self.DOS_num_of_bands,
                    spinpol=self.Spin_calc,
                    magnetic_moments=magnetic_moments,
                    setup_params=self.Setup_params,
                    xc_calc=self.XC_calc,
                    exx_fraction=self.XC_exx_fraction,
                    omega=self.XC_omega,
                    occupation=dos_occupation,
                    emin=self.Energy_min,
                    emax=self.Energy_max,
                    delta_e=delta_e,
                    bz_sum=bz_sum,
                    parallel_cores=self.parallel_cores,
                    relative_to_fermi=True,
                    scf_executable='pw.x',
                    dos_executable='dos.x',
                    projwfc_executable='projwfc.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    f"QE hybrid DOS workflow failed: {exc}"
                )
                raise

            workflow = hybrid_workflow['scf']
            result = workflow['result']
            fermi_energy = hybrid_workflow.get(
                'fermi_energy_ev'
            )

            if fermi_energy is None:
                raise RuntimeError(
                    "QE hybrid DOS requires an energy reference from "
                    "the hybrid SCF calculation."
                )

            dos_workflow = hybrid_workflow['dos']
            pdos_workflow = hybrid_workflow['pdos']

            parprint(
                "QE hybrid SCF, DOS and PDOS workflow finished."
            )
        else:
            parprint(
                "Starting QE total DOS calculation..."
            )

            try:
                dos_workflow = self.engine.run_dos(
                    input_file=dos_input_file,
                    output_file=dos_output_file,
                    state_dir=state_dir,
                    dos_file=dos_data_file,
                    emin=dos_emin_absolute,
                    emax=dos_emax_absolute,
                    delta_e=delta_e,
                    bz_sum=bz_sum,
                    parallel_cores=self.parallel_cores,
                    executable='dos.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    f"QE total DOS calculation failed: {exc}"
                )
                raise

            parprint(
                "QE total DOS calculation finished."
            )

        parprint(
            "QE DOS data saved to: "
            f"{dos_workflow['dos_file']}"
        )
        
        dos_result = self.engine.parse_dos_output(
            dos_workflow['dos_file']
        )

        if (
            dos_result['spin_polarized']
            != bool(self.Spin_calc)
        ):
            raise RuntimeError(
                "QE DOS spin channels do not match "
                "the requested Spin_calc setting."
            )

        shifted_energies = [
            energy - fermi_energy
            for energy in dos_result['energies_ev']
        ]

        if self.Spin_calc:
            dos_up_file = Path(
                self.struct
                + f'-DOS-{self.Engine}-Result-DOS-Up.csv'
            )

            dos_down_file = Path(
                self.struct
                + f'-DOS-{self.Engine}-Result-DOS-Down.csv'
            )

            with dos_up_file.open(
                'w',
                encoding='utf-8',
            ) as fd:
                for energy, dos_value in zip(
                    shifted_energies,
                    dos_result['dos_up'],
                ):
                    print(
                        energy,
                        dos_value,
                        sep=', ',
                        file=fd,
                    )

            with dos_down_file.open(
                'w',
                encoding='utf-8',
            ) as fd:
                for energy, dos_value in zip(
                    shifted_energies,
                    dos_result['dos_down'],
                ):
                    print(
                        energy,
                        dos_value,
                        sep=', ',
                        file=fd,
                    )

            parprint(
                "Saving spin-resolved DOS..."
            )

        else:
            csv_file = Path(
                self.struct
                + f'-DOS-{self.Engine}-Result-DOS.csv'
            )

            with csv_file.open(
                'w',
                encoding='utf-8',
            ) as fd:
                for energy, dos_value in zip(
                    shifted_energies,
                    dos_result['dos'],
                ):
                    print(
                        energy,
                        dos_value,
                        sep=', ',
                        file=fd,
                    )

            parprint(
                "Saving DOS..."
            )

        if world.rank == 0:
            fig, ax = plt.subplots(
                figsize=(8, 6)
            )

            energies = np.array(
                shifted_energies
            )

            if self.Spin_calc:
                dos_up = np.array(
                    dos_result['dos_up']
                )

                dos_down = np.array(
                    dos_result['dos_down']
                )

                ax.plot(
                    energies,
                    dos_up,
                    'b',
                    linewidth=1.5,
                    label='Spin Up',
                )

                ax.plot(
                    energies,
                    -dos_down,
                    'r',
                    linewidth=1.5,
                    label='Spin Down',
                )

                ax.fill_between(
                    energies,
                    0,
                    dos_up,
                    facecolor='blue',
                    alpha=0.2,
                )

                ax.fill_between(
                    energies,
                    0,
                    -dos_down,
                    facecolor='red',
                    alpha=0.2,
                )

                ax.legend()

            else:
                dos_values = np.array(
                    dos_result['dos']
                )

                ax.plot(
                    energies,
                    dos_values,
                    'b',
                    linewidth=1.5,
                )

                ax.fill_between(
                    energies,
                    0,
                    dos_values,
                    facecolor='blue',
                    alpha=0.2,
                )

            ax.set_xlabel(
                self._t("fig_dos_xlabel")
            )

            ax.set_ylabel(
                self._t("fig_dos_ylabel")
            )

            ax.axvline(
                x=0,
                color='k',
                linestyle='--',
                linewidth=1,
            )

            ax.axhline(
                y=0,
                color='k',
                linewidth=0.8,
            )

            ax.set_xlim(
                self.Energy_min,
                self.Energy_max,
            )

            autoscale_y(
                ax
            )

            plt.tight_layout()

            plt.savefig(
                self.struct
                + f'-DOS-{self.Engine}-Graph-DOS.png',
                dpi=300,
            )

            plt.close(
                fig
            )

        if not hybrid:
            parprint(
                "Starting QE projected DOS calculation..."
            )

            try:
                pdos_workflow = self.engine.run_projwfc(
                    input_file=pdos_input_file,
                    output_file=pdos_output_file,
                    state_dir=state_dir,
                    pdos_prefix=pdos_prefix,
                    emin=dos_emin_absolute,
                    emax=dos_emax_absolute,
                    delta_e=delta_e,
                    parallel_cores=self.parallel_cores,
                    executable='projwfc.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    f"QE projected DOS calculation failed: {exc}"
                )
                raise

        parprint(
            "QE projected DOS calculation finished."
        )

        parprint(
            "QE PDOS summary saved to: "
            f"{pdos_workflow['pdos_tot_file']}"
        )
        
        pdos_result = self.engine.aggregate_projwfc_pdos(
            pdos_workflow['pdos_prefix']
        )

        if (
            pdos_result['spin_polarized']
            != bool(self.Spin_calc)
        ):
            raise RuntimeError(
                "QE PDOS spin channels do not match "
                "the requested Spin_calc setting."
            )

        pdos_shifted_energies = [
            energy - fermi_energy
            for energy in pdos_result['energies_ev']
        ]

        def write_pdos_csv(
            output_path,
            projection,
        ):
            with output_path.open(
                'w',
                encoding='utf-8',
            ) as fd:
                print(
                    "Energy, "
                    "s-total, "
                    "p-total, "
                    "pz, "
                    "px, "
                    "py, "
                    "d-total, "
                    "d3z2_r2, "
                    "dxz, "
                    "dyz, "
                    "dx2_y2, "
                    "dxy, "
                    "f-total, "
                    "TOTAL",
                    file=fd,
                )

                for values in zip(
                    pdos_shifted_energies,
                    projection['s_total'],
                    projection['p_total'],
                    projection['pz'],
                    projection['px'],
                    projection['py'],
                    projection['d_total'],
                    projection['d3z2_r2'],
                    projection['dxz'],
                    projection['dyz'],
                    projection['dx2_y2'],
                    projection['dxy'],
                    projection['f_total'],
                    projection['total'],
                ):
                    print(
                        *values,
                        sep=', ',
                        file=fd,
                    )

        if self.Spin_calc:
            pdos_up_file = Path(
                self.struct
                + f'-DOS-{self.Engine}-Result-PDOS-Up.csv'
            )

            pdos_down_file = Path(
                self.struct
                + f'-DOS-{self.Engine}-Result-PDOS-Down.csv'
            )

            write_pdos_csv(
                pdos_up_file,
                pdos_result['spin_up'],
            )

            write_pdos_csv(
                pdos_down_file,
                pdos_result['spin_down'],
            )

            parprint(
                "Saving spin-resolved PDOS..."
            )

            parprint(
                "Nanoworks spin-up PDOS data saved to: "
                f"{pdos_up_file}"
            )

            parprint(
                "Nanoworks spin-down PDOS data saved to: "
                f"{pdos_down_file}"
            )

        else:
            pdos_csv_file = Path(
                self.struct
                + f'-DOS-{self.Engine}-Result-PDOS.csv'
            )

            write_pdos_csv(
                pdos_csv_file,
                pdos_result,
            )

            parprint(
                "Saving PDOS..."
            )

            parprint(
                "Nanoworks PDOS data saved to: "
                f"{pdos_csv_file}"
            )
        
        time22 = time.time()

        with paropen(
            self.struct + f'-TIMINGS-{self.Engine}-Log-Timings.txt',
            'a',
        ) as f1:
            print(
                f'DOS calculation: '
                f'{round((time22-time21),2)}',
                file=f1,
            )
    
    def bandcalc(self):
        """Run the band-structure workflow using the selected DFT engine."""
        if self.Engine == 'GPAW':
            return self._bandcalc_gpaw()

        if self.Engine == 'QE':
            return self._bandcalc_qe()

        raise ValueError(
            f"Unsupported DFT engine: {self.Engine}"
        )

    def _bandcalc_qe(self):
        """Run the Quantum ESPRESSO band-structure workflow."""
        time31 = time.time()

        parprint(
            "Starting QE band structure calculation..."
        )

        if self.Mode != 'PW':
            parprint(
                "\033[91mERROR:\033[0m "
                "Quantum ESPRESSO band calculations "
                "support PW mode only."
            )
            sys.exit(1)

        if self.SOC_calc:
            parprint(
                "\033[91mERROR:\033[0m "
                "Quantum ESPRESSO SOC band calculations "
                "are not supported yet."
            )
            sys.exit(1)

        try:
            validated_xc = self.engine.validate_qe_xc(
                self.XC_calc,
                pseudo_xc='pbe',
                allow_hybrid=True,
            )
        except ValueError as exc:
            parprint(
                f"\033[91mERROR:\033[0m {exc}"
            )
            sys.exit(1)

        hybrid = str(
            validated_xc
        ).strip().lower() in {
            'hse06',
            'hse03',
            'pbe0',
        }

        ground_state_dir = Path(
            self.struct
            + '-GROUND-QE-Result-State'
        )

        if not self.engine.has_qe_state(
            ground_state_dir,
            prefix='nanoworks',
        ):
            parprint(
                "\033[91mERROR:\033[0m "
                + str(ground_state_dir)
                + " does not contain a valid QE ground-state result. "
                "Complete the ground-state calculation first."
            )
            sys.exit(1)

        state_dir = (
            Path(
                self.struct
                + '-BAND-QE-Result-State'
            )
            if hybrid
            else ground_state_dir
        )

        pseudo_dir = get_qe_pseudo_dir(
            relativistic='scalar',
        )

        try:
            pseudopotentials = (
                resolve_qe_pseudopotentials(
                    self.bulk_configuration,
                    relativistic='scalar',
                )
            )
        except (FileNotFoundError, RuntimeError) as exc:
            parprint(
                f"\033[91mERROR:\033[0m {exc}"
            )
            sys.exit(1)

        magnetic_moments = None

        if self.Spin_calc:
            magnetic_moments = (
                resolve_initial_magnetic_moments(
                    atoms=self.bulk_configuration,
                    magmom_per_atom=self.Magmom_per_atom,
                    magmom_single_atom=(
                        self.Magmom_single_atom
                    ),
                )
            )

        try:
            band_path = self.engine.build_band_path(
                self.bulk_configuration,
                path=self.Band_path,
                npoints=self.Band_npoints,
            )
        except (TypeError, ValueError) as exc:
            parprint(
                "\033[91mERROR:\033[0m "
                f"QE band path could not be prepared: {exc}"
            )
            raise

        if hybrid:
            input_file = Path(
                self.struct
                + '-BAND-QE-Input-Hybrid-SCF.in'
            )
            output_file = Path(
                self.struct
                + '-BAND-QE-Log-Hybrid-SCF.txt'
            )
            bands_input_file = Path(
                self.struct
                + '-BAND-QE-Input-Bands.x.in'
            )
            bands_output_file = Path(
                self.struct
                + '-BAND-QE-Log-Bands.x.txt'
            )
            band_data_file = Path(
                self.struct
                + '-BAND-QE-Result-Bands.x.dat'
            )
        else:
            input_file = Path(
                self.struct
                + '-BAND-QE-Input-Bands.in'
            )
            output_file = Path(
                self.struct
                + f'-BAND-{self.Engine}-Log-Bands.txt'
            )

        projection_input_file = Path(
            self.struct
            + '-BAND-QE-Input-Projections.in'
        )

        projection_output_file = Path(
            self.struct
            + '-BAND-QE-Log-Projections.txt'
        )

        projection_prefix = Path(
            self.struct
            + '-BAND-QE-Result-Projections'
        )

        try:
            if hybrid:
                hybrid_gamma = (
                    self.Gamma
                    if self.Ground_gamma is None
                    else self.Ground_gamma
                )

                band_mesh = (
                    self.engine.resolve_qe_kpoint_size(
                        self.bulk_configuration,
                        density=self.Ground_kpts_density,
                        size=(
                            self.Ground_kpts_x,
                            self.Ground_kpts_y,
                            self.Ground_kpts_z,
                        ),
                    )
                )

                hybrid_workflow = (
                    self.engine.run_hybrid_bands(
                        atoms=self.bulk_configuration,
                        scf_input_file=input_file,
                        scf_output_file=output_file,
                        bands_input_file=bands_input_file,
                        bands_output_file=bands_output_file,
                        state_dir=state_dir,
                        band_file=band_data_file,
                        pseudopotentials=pseudopotentials,
                        pseudo_dir=pseudo_dir,
                        cutoff_ev=self.Cut_off_energy,
                        band_path=band_path,
                        qpoint_grid=band_mesh,
                        kpoint_size=band_mesh,
                        gamma=hybrid_gamma,
                        total_charge=self.Total_charge,
                        nbands=self.Band_num_of_bands,
                        spinpol=self.Spin_calc,
                        magnetic_moments=magnetic_moments,
                        setup_params=self.Setup_params,
                        xc_calc=self.XC_calc,
                        exx_fraction=self.XC_exx_fraction,
                        omega=self.XC_omega,
                        occupation=self.Occupation,
                        parallel_cores=self.parallel_cores,
                        scf_executable='pw.x',
                        bands_executable='bands.x',
                        prefix='nanoworks',
                        projected_band=(
                            self.Projected_band_plot
                        ),
                        projections=self.Projections,
                        projection_input_file=(
                            projection_input_file
                        ),
                        projection_output_file=(
                            projection_output_file
                        ),
                        projection_prefix=(
                            projection_prefix
                        ),
                        projection_executable='projwfc.x',
                    )
                )
                workflow = {
                    'result': hybrid_workflow[
                        'scf'
                    ][
                        'result'
                    ],
                    'bands': hybrid_workflow[
                        'band_data'
                    ],
                    'band_projections': hybrid_workflow[
                        'band_projections'
                    ],
                    'hybrid_workflow': hybrid_workflow,
                }
            else:
                workflow = self.engine.run_bands(
                    atoms=self.bulk_configuration,
                    input_file=input_file,
                    output_file=output_file,
                    state_dir=state_dir,
                    pseudopotentials=pseudopotentials,
                    pseudo_dir=pseudo_dir,
                    cutoff_ev=self.Cut_off_energy,
                    band_path=band_path,
                    total_charge=self.Total_charge,
                    nbands=self.Band_num_of_bands,
                    spinpol=self.Spin_calc,
                    magnetic_moments=magnetic_moments,
                    setup_params=self.Setup_params,
                    xc_calc=self.XC_calc,
                    exx_fraction=self.XC_exx_fraction,
                    omega=self.XC_omega,
                    occupation=self.Occupation,
                    parallel_cores=self.parallel_cores,
                    executable='pw.x',
                    prefix='nanoworks',
                    projected_band=(
                        self.Projected_band_plot
                    ),
                    projections=self.Projections,
                    projection_input_file=(
                        projection_input_file
                    ),
                    projection_output_file=(
                        projection_output_file
                    ),
                    projection_prefix=(
                        projection_prefix
                    ),
                    projection_executable='projwfc.x',
                )
        except Exception as exc:
            parprint(
                "\033[91mERROR:\033[0m "
                f"QE band calculation failed: {exc}"
            )
            raise

        reference_candidates = []

        ground_output_file = Path(
            self.struct
            + f'-GROUND-{self.Engine}-Log-SCF.txt'
        )

        if hybrid:
            reference_candidates.append(
                workflow['result']
            )

        if ground_output_file.is_file():
            reference_candidates.append(
                self.engine.parse_pw_output(
                    ground_output_file
                )
            )

        if not hybrid:
            reference_candidates.append(
                workflow['result']
            )

        reference = None

        for candidate in reference_candidates:
            try:
                reference = (
                    self.engine.resolve_qe_band_reference(
                        candidate
                    )
                )
                break
            except ValueError:
                continue

        if reference is None:
            raise RuntimeError(
                "QE band reference energy could not be resolved "
                "from either the ground-state or band calculation."
            )

        bands = workflow[
            'bands'
        ]

        band_projections = workflow[
            'band_projections'
        ]

        if self.Projected_band_plot:
            if band_projections is None:
                raise RuntimeError(
                    "QE projected-band calculation "
                    "did not return projection data."
                )

            projection_count = len(
                band_projections[
                    'up'
                ][
                    'projections'
                ]
            )

            parprint(
                "QE band projections prepared: "
                f"{projection_count} selection(s)."
            )

            parprint(
                "QE band projection data prefix: "
                f"{projection_prefix}"
            )

        parprint(
            "QE band calculation finished."
        )

        parprint(
            "QE band k-points: "
            f"{bands['nkpoints']}"
        )

        parprint(
            "QE bands per k-point: "
            f"{bands['nbands']}"
        )

        parprint(
            "QE band reference energy: "
            f"{reference['energy_ev']:.8f} eV "
            f"({reference['source']})"
        )

        band_data = (
            self.engine.prepare_qe_band_data(
                bands=bands,
                band_path=band_path,
                reference_energy=reference[
                    'energy_ev'
                ],
            )
        )

        distances = band_data[
            'distances'
        ]

        if band_data['spin_polarized']:
            band_channels = [
                (
                    'Up',
                    band_data[
                        'eigenvalues_up_ev'
                    ],
                ),
                (
                    'Down',
                    band_data[
                        'eigenvalues_down_ev'
                    ],
                ),
            ]
        else:
            band_channels = [
                (
                    None,
                    band_data[
                        'eigenvalues_ev'
                    ],
                ),
            ]

        for spin_label, eigenvalues in band_channels:
            suffix = (
                f'-{spin_label}'
                if spin_label is not None
                else ''
            )

            band_file = Path(
                self.struct
                + f'-BAND-{self.Engine}-Result-Band'
                + suffix
                + '.dat'
            )

            with paropen(
                band_file,
                'w',
            ) as fd:
                for band_index in range(
                    band_data['nbands']
                ):
                    for kpoint_index in range(
                        band_data['nkpoints']
                    ):
                        print(
                            kpoint_index,
                            eigenvalues[
                                kpoint_index
                            ][
                                band_index
                            ],
                            file=fd,
                        )

                    print(
                        file=fd,
                    )

            xyyy_file = Path(
                self.struct
                + f'-BAND-{self.Engine}-Result-Band'
                + suffix
                + '-XYYY.dat'
            )

            with paropen(
                xyyy_file,
                'w',
            ) as fd:
                for kpoint_index, values in enumerate(
                    eigenvalues
                ):
                    print(
                        kpoint_index,
                        *values,
                        file=fd,
                    )

            parprint(
                "QE band data saved to: "
                f"{band_file}"
            )

            parprint(
                "QE XYYY band data saved to: "
                f"{xyyy_file}"
            )

        time32 = time.time()

        with paropen(
            self.struct
            + f'-TIMINGS-{self.Engine}-Log-Timings.txt',
            'a',
        ) as f1:
            print(
                f'Band calculation: '
                f'{round((time32-time31), 2)}',
                file=f1,
            )

        if world.rank == 0:
            fig, ax = plt.subplots(
                figsize=(8, 6)
            )

            if band_data['spin_polarized']:
                plot_channels = [
                    (
                        band_data[
                            'eigenvalues_up_ev'
                        ],
                        'blue',
                        '-',
                        'Spin up',
                    ),
                    (
                        band_data[
                            'eigenvalues_down_ev'
                        ],
                        'red',
                        '--',
                        'Spin down',
                    ),
                ]
            else:
                plot_channels = [
                    (
                        band_data[
                            'eigenvalues_ev'
                        ],
                        'blue',
                        '-',
                        None,
                    ),
                ]

            for (
                channel_values,
                color,
                linestyle,
                channel_label,
            ) in plot_channels:
                energy_array = np.array(
                    channel_values
                )

                for band_index in range(
                    band_data['nbands']
                ):
                    ax.plot(
                        distances,
                        energy_array[
                            :,
                            band_index,
                        ],
                        color=color,
                        linestyle=linestyle,
                        linewidth=1.0,
                        label=(
                            channel_label
                            if band_index == 0
                            else None
                        ),
                    )

            if band_data['spin_polarized']:
                ax.legend()

            for special_distance in band_data[
                'special_distances'
            ]:
                ax.axvline(
                    x=special_distance,
                    color='black',
                    linewidth=0.6,
                    alpha=0.5,
                )

            ax.axhline(
                y=0.0,
                color='black',
                linestyle='--',
                linewidth=0.8,
            )

            labels = [
                (
                    r'$\Gamma$'
                    if label in {
                        'G',
                        'Gamma',
                        'Γ',
                    }
                    else label
                )
                for label in band_data[
                    'labels'
                ]
            ]

            ax.set_xticks(
                band_data[
                    'special_distances'
                ]
            )

            ax.set_xticklabels(
                labels
            )

            ax.set_ylabel(
                self._t(
                    "fig_band_ylabel"
                )
            )

            ax.set_xlim(
                distances[0],
                distances[-1],
            )

            ax.set_ylim(
                self.Energy_min,
                self.Energy_max,
            )

            plt.tight_layout()

            graph_file = Path(
                self.struct
                + f'-BAND-{self.Engine}-Graph-Band.png'
            )

            plt.savefig(
                graph_file,
                dpi=300,
            )

            plt.close(
                fig
            )

            parprint(
                "QE band graph saved to: "
                f"{graph_file}"
            )

        if self.Projected_band_plot:
            parprint(
                "Drawing QE projected bands..."
            )

            self._draw_projectedband_qe(
                band_data=band_data,
                band_projections=(
                    band_projections
                ),
            )

    def _draw_projectedband_qe(
        self,
        band_data,
        band_projections,
    ):
        """Draw QE projected-band graphs."""
        if world.rank != 0:
            return

        if band_projections is None:
            raise RuntimeError(
                "QE projected-band data is missing."
            )

        from matplotlib.lines import Line2D

        distances = np.asarray(
            band_data['distances'],
            dtype=float,
        )

        labels = [
            (
                r'$\Gamma$'
                if label in {
                    'G',
                    'Gamma',
                    'Γ',
                }
                else label
            )
            for label in band_data[
                'labels'
            ]
        ]

        if band_data['spin_polarized']:
            channels = [
                (
                    'Up',
                    band_data[
                        'eigenvalues_up_ev'
                    ],
                    band_projections[
                        'up'
                    ],
                ),
                (
                    'Down',
                    band_data[
                        'eigenvalues_down_ev'
                    ],
                    band_projections[
                        'down'
                    ],
                ),
            ]
        else:
            channels = [
                (
                    None,
                    band_data[
                        'eigenvalues_ev'
                    ],
                    band_projections[
                        'up'
                    ],
                ),
            ]

        for (
            spin_label,
            eigenvalues,
            projection_data,
        ) in channels:
            if projection_data is None:
                raise RuntimeError(
                    "QE projected-band data is missing "
                    f"for spin channel: {spin_label}"
                )

            energy_array = np.asarray(
                eigenvalues,
                dtype=float,
            )

            expected_shape = (
                band_data['nkpoints'],
                band_data['nbands'],
            )

            if energy_array.shape != expected_shape:
                raise RuntimeError(
                    "QE projected-band energies have "
                    "an inconsistent shape."
                )

            fig, ax = plt.subplots(
                figsize=(8, 6)
            )

            for band_index in range(
                band_data['nbands']
            ):
                ax.plot(
                    distances,
                    energy_array[
                        :,
                        band_index,
                    ],
                    color='gray',
                    linewidth=0.5,
                    zorder=1,
                )

            legend_handles = []

            for projection in projection_data[
                'projections'
            ]:
                weight_array = np.asarray(
                    projection[
                        'weights'
                    ],
                    dtype=float,
                )

                if weight_array.shape != expected_shape:
                    plt.close(
                        fig
                    )

                    raise RuntimeError(
                        "QE projected-band weights have "
                        "an inconsistent shape."
                    )

                marker_sizes = (
                    np.maximum(
                        weight_array,
                        0.0,
                    )
                    * 80.0
                )

                for band_index in range(
                    band_data['nbands']
                ):
                    ax.scatter(
                        distances,
                        energy_array[
                            :,
                            band_index,
                        ],
                        s=marker_sizes[
                            :,
                            band_index,
                        ],
                        color=projection[
                            'color'
                        ],
                        zorder=2,
                        alpha=0.6,
                        edgecolors='none',
                    )

                legend_handles.append(
                    Line2D(
                        [0],
                        [0],
                        marker='o',
                        linestyle='None',
                        markerfacecolor=(
                            projection[
                                'color'
                            ]
                        ),
                        markeredgecolor='none',
                        markersize=8,
                        label=projection[
                            'label'
                        ],
                    )
                )

            ax.set_xticks(
                band_data[
                    'special_distances'
                ]
            )

            ax.set_xticklabels(
                labels
            )

            ax.set_ylabel(
                self._t(
                    "fig_band_ylabel"
                )
            )

            ax.set_ylim(
                self.Energy_min,
                self.Energy_max,
            )

            ax.set_xlim(
                distances[0],
                distances[-1],
            )

            ax.axhline(
                y=0.0,
                color='black',
                linewidth=1.0,
                linestyle='-',
            )

            for special_distance in band_data[
                'special_distances'
            ]:
                ax.axvline(
                    x=special_distance,
                    color='black',
                    linewidth=0.5,
                    linestyle='--',
                )

            if legend_handles:
                ax.legend(
                    handles=legend_handles,
                    loc='lower left',
                    frameon=True,
                    fontsize=10,
                )

            plt.tight_layout()

            if spin_label is None:
                graph_file = Path(
                    self.struct
                    + '-BAND-QE-Graph-Projected-Band.png'
                )
            else:
                graph_file = Path(
                    self.struct
                    + '-BAND-QE-Graph-Projected-Band'
                    + f'-Spin-{spin_label}.png'
                )

            plt.savefig(
                graph_file,
                dpi=300,
            )

            plt.close(
                fig
            )

            parprint(
                "QE projected-band graph saved to: "
                f"{graph_file}"
            )

    def _bandcalc_gpaw(self):
        """
        This method performs band structure calculations for the given structure using the
        ground state results. It computes the electronic band structure along specified
        k-point paths and saves the results in appropriate files for further analysis
        and visualization.
        """

        # -------------------------------------------------------------
        # BAND STRUCTURE CALCULATION
        # -------------------------------------------------------------

        # Start Band calc
        time31 = time.time()
        parprint("Starting band structure calculation...")

        hybrid = self.engine.is_hybrid(self.XC_calc)

        calc = self.engine.prepare_band_calc(
            filename=self.struct+'-GROUND-GPAW-Result-State.gpw',
            hybrid=hybrid,
            path=self.Band_path,
            npoints=self.Band_npoints,
            txt=self.struct+f'-BAND-{self.Engine}-Log-Bands.txt',
            occupations=self.Occupation,
            convergence=self.Band_convergence,
            nbands=self.Band_num_of_bands,
        )

        if hybrid:
            ef = self.hybrid_fermi_level()
        else:
            ef = calc.get_fermi_level()

        calc.get_potential_energy()

        # For hybrids, refine the reference now that the calculator is
        # populated, in case the ground-state value was unavailable above.
        if hybrid and (ef is None or ef == 0.0):
            ef = self.hybrid_fermi_level(calc)

        bs = calc.band_structure()

        calc.get_potential_energy()
        # For hybrids, refine the reference now that the calculator is
        # populated, in case the ground-state value was unavailable above.
        if hybrid and (ef is None or ef == 0.0):
            ef = self.hybrid_fermi_level(calc)
        bs = calc.band_structure()
            
        Band_num_of_bands = calc.get_number_of_bands()
        parprint('Num of bands:'+str(Band_num_of_bands))

        # No need to write an additional gpaw file. Use json file to use with ase band-structure command
        #calc.write(self.struct+'-BAND-GPAW-Result-Band.gpw')
        bs.write(self.struct+'-BAND-GPAW-Result-State.json')

        if self.Spin_calc == True:
            eps_skn = np.array([[calc.get_eigenvalues(k,s)
                                for k in range(self.Band_npoints)]
                                for s in range(2)]) - ef
            parprint(eps_skn.shape)
            with paropen(self.struct+f'-BAND-{self.Engine}-Result-Band-Down.dat', 'w') as f1:
                for n1 in range(Band_num_of_bands):
                    for k1 in range(self.Band_npoints):
                        print(k1, eps_skn[0, k1, n1], end="\n", file=f1)
                    print (end="\n", file=f1)

            with paropen(self.struct+f'-BAND-{self.Engine}-Result-Band-Up.dat', 'w') as f2:
                for n2 in range(Band_num_of_bands):
                    for k2 in range(self.Band_npoints):
                        print(k2, eps_skn[1, k2, n2], end="\n", file=f2)
                    print (end="\n", file=f2)

            # Thanks to Andrej Kesely (https://stackoverflow.com/users/10035985/andrej-kesely) for helping the problem of general XYYY writer
            currentd, all_groupsd = [], []
            with open(self.struct+f'-BAND-{self.Engine}-Result-Band-Down.dat', 'r') as f_in1:
                for line in map(str.strip, f_in1):
                    if line == "" and currentd:
                        all_groupsd.append(currentd)
                        currentd = []
                    else:
                        currentd.append(line.split(maxsplit=1))

            if currentd:
                all_groupsd.append(currentd)

            try:
                with paropen(self.struct+f'-BAND-{self.Engine}-Result-Band-Down-XYYY.dat', 'w') as f1:
                    for g in zip(*all_groupsd):
                        print('{} {} {}'.format(g[0][0], g[0][1], ' '.join(v for _, v in g[1:])), file=f1)
            except Exception as e:
                print("\033[93mWARNING:\033[0m A problem occurred during writing XYYY formatted spin down Band file. Mostly, the file is created without any problem.")
                print(e)

            currentu, all_groupsu = [], []
            with open(self.struct+f'-BAND-{self.Engine}-Result-Band-Up.dat', 'r') as f_in2:
                for line in map(str.strip, f_in2):
                    if line == "" and currentu:
                        all_groupsu.append(currentu)
                        currentu = []
                    else:
                        currentu.append(line.split(maxsplit=1))

            if currentu:
                all_groupsu.append(currentu)
            try:
                with paropen(self.struct+f'-BAND-{self.Engine}-Result-Band-Up-XYYY.dat', 'w') as f2:
                    for g in zip(*all_groupsu):
                        print('{} {} {}'.format(g[0][0], g[0][1], ' '.join(v for _, v in g[1:])), file=f2)
            except Exception as e:
                print("\033[93mWARNING:\033[0m A problem occurred during writing XYYY formatted spin up Band file. Mostly, the file is created without any problem.")
                print(e)

        else:
            if self.SOC_calc:
                from gpaw.spinorbit import soc_eigenstates
                parprint("Applying Spin-Orbit Coupling to Band Structure...")
                # Spin-Orbit perturbation
                soc = soc_eigenstates(calc)
                soc_evals = np.real(soc.eigenvalues()) # Taking only real parts
                # GPAW returns SOC eigenvalues ​​as a single spin channel (nkpts, 2*nbands)
                eps_skn = np.array([soc_evals]) - ef
                Band_num_of_bands = eps_skn.shape[2]
            else:
                eps_skn = np.array([[calc.get_eigenvalues(k,s)
                                for k in range(self.Band_npoints)]
                                for s in range(1)]) - ef
            with paropen(self.struct+f'-BAND-{self.Engine}-Result-Band.dat', 'w') as f:
                for n in range(Band_num_of_bands):
                    for k in range(self.Band_npoints):
                        print(k, eps_skn[0, k, n], end="\n", file=f)
                    print (end="\n", file=f)

            # Thanks to Andrej Kesely (https://stackoverflow.com/users/10035985/andrej-kesely) for helping the problem of general XYYY writer
            current, all_groups = [], []
            with open(self.struct+f'-BAND-{self.Engine}-Result-Band.dat', 'r') as f_in:
                for line in map(str.strip, f_in):
                    if line == "" and current:
                        all_groups.append(current)
                        current = []
                    else:
                        current.append(line.split(maxsplit=1))

            if current:
                all_groups.append(current)
            try:
                with paropen(self.struct+f'-BAND-{self.Engine}-Result-Band-XYYY.dat', 'w') as f1:
                    for g in zip(*all_groups):
                        print('{} {} {}'.format(g[0][0], g[0][1], ' '.join(v for _, v in g[1:])), file=f1)
            except Exception as e:
                print("\033[93mWARNING:\033[0m A problem occurred during writing XYYY formatted Band file. Mostly, the file is created without any problem.")
                print(e)
                
        # Finish Band calc
        time32 = time.time()
        # Write timings of calculation
        with paropen(self.struct+f'-TIMINGS-{self.Engine}-Log-Timings.txt', 'a') as f1:
            print('Band calculation: ', round((time32-time31),2), end="\n", file=f1)

        # Write or draw figures
        # Draw graphs only on master node
        if world.rank == 0:
            # Band Structure
            if self.SOC_calc and not self.Spin_calc:
                from ase.spectrum.band_structure import BandStructure
                bs = BandStructure(path=bs.path, energies=np.array([soc_evals]), reference=ef)
            bs.plot(filename=self.struct+f'-BAND-{self.Engine}-Graph-Band.png', show=False, emax=self.Energy_max + bs.reference, emin=self.Energy_min + bs.reference, ylabel=self._t("fig_band_ylabel"))
            
        # Projected band
        if self.Projected_band_plot == True:
            parprint(f"Drawing Projected Band...")
            self._draw_projectedband(calc)


    def _draw_projectedband(self, calc):
        """
        Internal method for dftsolve class to plot projected band structures.
        """

        from gpaw.utilities.dos import get_angular_projectors
        
        # Get number of spins (1 for non-magnetic, 2 for spin-polarized)
        nspins = calc.get_number_of_spins()
        
        filename = f"{self.struct}-BAND-GPAW-Graph-Projected-Band.png"
        
        Projections = self.Projections

        if not Projections:
            if world.rank == 0:
                print(
                    "Warning: 'Projections' list is missing "
                    "in config. Auto-generating total projection..."
                )

            total_atoms = len(
                calc.get_atoms()
            )

            Projections = [{
                'atoms': list(
                    range(total_atoms)
                ),
                'orbital': None,
                'color': 'blue',
                'label': 'Total Contribution',
            }]
        
        bs = calc.band_structure()

        x, x_special, labels = bs.path.get_linear_kpoint_axis()
        
        # Loop over each spin channel
        for spin_index in range(nspins):
            # Dynamic filename based on spin index
            if nspins == 1:
                filename = f"{self.struct}-BAND-GPAW-Graph-Projected-Band.png"
            else:
                spin_label = "Up" if spin_index == 0 else "Down"
                filename = f"{self.struct}-BAND-GPAW-Graph-Projected-Band-Spin-{spin_label}.png"
            
            energies = bs.energies[spin_index]  
            nbands = energies.shape[1]
            projection_data = []
            
            
            for proj in Projections:
                atom_indices = proj.get('atoms', [])
                orbital = proj.get('orbital', None)
                color = proj.get('color', 'blue')
                label = proj.get('label', f"Atoms: {atom_indices}")
                
                # local weights matrice
                weights_local = np.zeros(energies.shape)
                
                for kpt in calc.wfs.kpt_u:
                    # Skip if the k-point doesn't belong to the current spin channel
                    if kpt.s != spin_index:
                        continue
                    k = kpt.k  # Global k-point index
                    P_ani = kpt.P_ani
                    for a in atom_indices:
                        if a in P_ani:
                            if orbital is not None:
                                setup = calc.wfs.setups[a]
                                proj_indices = get_angular_projectors(setup, orbital, type='bound')
                                w_n = np.sum(np.abs(P_ani[a][:, proj_indices])**2, axis=1)
                            else:
                                w_n = np.sum(np.abs(P_ani[a])**2, axis=1)
                            # write weights to local k-point index
                            weights_local[k, :] += w_n

                world.sum(weights_local)
                
                projection_data.append((weights_local, color, label))

            legend_handles = []
            # draw
            if world.rank == 0:
                import matplotlib.pyplot as plt
                from matplotlib.lines import Line2D
                fig, ax = plt.subplots(figsize=(8, 6))
                
                # background standard bands
                for n in range(nbands):
                    ax.plot(x, energies[:, n], color='gray', lw=0.5, zorder=1)
                    
                # Projections are scatters
                for weights, color, label in projection_data:
                    # Draw all scatter points without legend labels
                    for n in range(nbands):
                        ax.scatter(
                            x,
                            energies[:, n],
                            s=weights[:, n] * 80,
                            color=color,
                            zorder=2,
                            alpha=0.6,
                            edgecolors='none'
                        )

                    # One fixed-size legend entry
                    legend_handles.append(
                        Line2D(
                            [0], [0],
                            marker='o',
                            linestyle='None',
                            markerfacecolor=color,
                            markeredgecolor='none',
                            markersize=8,
                            label=label
                        )
                    )

                ax.set_xticks(x_special)
                ax.set_xticklabels(labels)
                ax.set_ylabel(self._t("fig_band_ylabel"))
                ax.set_ylim(self.Energy_min + bs.reference, self.Energy_max + bs.reference)
                ax.set_xlim(x[0], x[-1])
                ax.axhline(0, color='black', lw=1.0, ls='-')
                for x_loc in x_special:
                    ax.axvline(x_loc, color='black', lw=0.5, ls='--')
                
                if legend_handles:
                    ax.legend(handles=legend_handles,
                              loc='lower left',
                              frameon=True,
                              fontsize=10)
                    
                plt.tight_layout()
                plt.savefig(filename, dpi=300)
                plt.close()

    def densitycalc(self):
        """Run the electron-density workflow for the selected engine."""

        if self.Engine == 'GPAW':
            return self._densitycalc_gpaw()

        if self.Engine == 'QE':
            return self._densitycalc_qe()

        raise ValueError(
            "Unsupported electron-density engine: "
            f"{self.Engine}"
        )

    def _densitycalc_gpaw(self):
        """
        This method performs density calculations for the given structure using the
        ground state results. It computes the electron density distribution and saves
        the results in appropriate files for further analysis and visualization.
        """

        # -------------------------------------------------------------
        # ALL-ELECTRON DENSITY
        # -------------------------------------------------------------

        #Start Density calc
        time41 = time.time()
        parprint("Starting All-electron density calculation...")
        calc = self.engine.load_gpaw_calc(self.struct+'-GROUND-GPAW-Result-State.gpw', txt=self.struct+'-EDENSITY-GPAW-Log-Calculation.txt')
        self.bulk_configuration.calc = calc
        if self.Spin_calc == True:
            np = calc.get_pseudo_density()
            n = calc.get_all_electron_density(gridrefinement=self.Refine_grid)
            # For spins
            npdown = calc.get_pseudo_density(spin=0)
            ndown = calc.get_all_electron_density(spin=0, gridrefinement=self.Refine_grid)
            npup = calc.get_pseudo_density(spin=1)
            nup = calc.get_all_electron_density(spin=1, gridrefinement=self.Refine_grid)
            # Zeta
            nzeta = (nup - ndown) / (nup + ndown)
            npzeta = (npup - npdown) / (npup + npdown)
            # Writing spin down pseudo and all electron densities to cube file with Bohr unit
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_nall-Down.cube', self.bulk_configuration, data=ndown * Bohr**3)
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_npseudo-Down.cube', self.bulk_configuration, data=npdown * Bohr**3)
            # Writing spin up pseudo and all electron densities to cube file with Bohr unit
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_nall-Up.cube', self.bulk_configuration, data=nup * Bohr**3)
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_npseudo-Up.cube', self.bulk_configuration, data=npup * Bohr**3)
            # Writing total pseudo and all electron densities to cube file with Bohr unit
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_nall-Total.cube', self.bulk_configuration, data=n * Bohr**3)
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_npseudo-Total.cube', self.bulk_configuration, data=np * Bohr**3)
            # Writing zeta pseudo and all electron densities to cube file with Bohr unit
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_nall-Zeta.cube', self.bulk_configuration, data=nzeta * Bohr**3)
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_npseudo-Zeta.cube', self.bulk_configuration, data=npzeta * Bohr**3)
        else:
            np = calc.get_pseudo_density()
            n = calc.get_all_electron_density(gridrefinement=self.Refine_grid)
            # Writing pseudo and all electron densities to cube file with Bohr unit
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_nall-Total.cube', self.bulk_configuration, data=n * Bohr**3)
            write(self.struct+'-EDENSITY-GPAW-Result-All-electron_npseudo-Total.cube', self.bulk_configuration, data=np * Bohr**3)
            
        # Finish Density calc
        time42 = time.time()
        # Write timings of calculation
        with paropen(self.struct+f'-TIMINGS-{self.Engine}-Log-Timings.txt', 'a') as f1:
            print('Density calculation: ', round((time42-time41),2), end="\n", file=f1)

    def _densitycalc_qe(self):
        """Generate QE pseudo-valence electron-density Cube files."""
        self.engine.validate_qe_xc(
            self.XC_calc,
            pseudo_xc='pbe',
            allow_hybrid=True,
        )

        state_dir = Path(
            self.struct
            + '-GROUND-QE-Result-State'
        )

        jobs = [
            {
                'label': 'Pseudo-Total',
                'description': 'total pseudo-valence density',
                'plot_num': 0,
                'spin_component': (
                    0
                    if self.Spin_calc
                    else None
                ),
            },
        ]

        if self.Spin_calc:
            jobs.extend([
                {
                    'label': 'Pseudo-Up',
                    'description': 'spin-up pseudo-valence density',
                    'plot_num': 0,
                    'spin_component': 1,
                },
                {
                    'label': 'Pseudo-Down',
                    'description': 'spin-down pseudo-valence density',
                    'plot_num': 0,
                    'spin_component': 2,
                },
                {
                    'label': 'Spin-Density',
                    'description': 'spin density',
                    'plot_num': 6,
                    'spin_component': None,
                },
            ])

        outputs = {}

        for job in jobs:
            label = job[
                'label'
            ]

            parprint(
                "Starting QE "
                + job['description']
                + " calculation..."
            )

            input_file = Path(
                self.struct
                + f'-EDENSITY-QE-Input-{label}.in'
            )

            output_file = Path(
                self.struct
                + f'-EDENSITY-QE-Log-{label}.txt'
            )

            filplot = Path(
                self.struct
                + f'-EDENSITY-QE-Result-{label}.dat'
            )

            cube_file = Path(
                self.struct
                + f'-EDENSITY-QE-Result-{label}.cube'
            )

            try:
                workflow = self.engine.run_pp_density(
                    input_file=input_file,
                    output_file=output_file,
                    state_dir=state_dir,
                    filplot=filplot,
                    cube_file=cube_file,
                    plot_num=job['plot_num'],
                    spin_component=(
                        job['spin_component']
                    ),
                    parallel_cores=self.parallel_cores,
                    executable='pp.x',
                    prefix='nanoworks',
                )
            except Exception as exc:
                parprint(
                    "\033[91mERROR:\033[0m "
                    "QE electron-density calculation "
                    f"failed for {label}: {exc}"
                )
                raise

            outputs[
                label
            ] = workflow

            parprint(
                "QE "
                + job['description']
                + " written to: "
                + str(cube_file)
            )

        return outputs

    def _plot_qe_phonon_results(
        self,
        output_file,
        band_path,
        frequencies,
        dos_data,
    ):
        """Plot native QE phonon bands and total DOS in THz."""
        output_file = Path(
            output_file
        )
        distances = np.asarray(
            band_path.get('distances', []),
            dtype=float,
        )
        band_frequencies = np.asarray(
            frequencies.get('frequencies_thz', []),
            dtype=float,
        )
        dos_frequencies = np.asarray(
            dos_data.get('frequencies_thz', []),
            dtype=float,
        )
        total_dos = np.asarray(
            dos_data.get('dos', []),
            dtype=float,
        ) / self.engine.THZ_PER_CM_MINUS_ONE

        if (
            distances.ndim != 1
            or distances.size == 0
            or band_frequencies.ndim != 2
            or band_frequencies.shape[0] != distances.size
        ):
            raise ValueError(
                "QE phonon band data dimensions are inconsistent."
            )

        if (
            dos_frequencies.ndim != 1
            or total_dos.ndim != 1
            or dos_frequencies.size == 0
            or dos_frequencies.size != total_dos.size
        ):
            raise ValueError(
                "QE phonon DOS data dimensions are inconsistent."
            )

        special_distances = list(
            band_path.get('special_distances', [])
        )
        labels = list(
            band_path.get('labels', [])
        )

        if len(special_distances) != len(labels):
            raise ValueError(
                "QE phonon high-symmetry labels and positions "
                "do not match."
            )

        display_labels = [
            r'$\Gamma$'
            if label == 'G'
            else label
            for label in labels
        ]

        output_file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        fig, (band_ax, dos_ax) = plt.subplots(
            1,
            2,
            figsize=(9, 6),
            sharey=True,
            gridspec_kw={
                'width_ratios': [3, 1],
                'wspace': 0.08,
            },
        )

        try:
            for mode_values in band_frequencies.T:
                band_ax.plot(
                    distances,
                    mode_values,
                    color='tab:blue',
                    linewidth=1.0,
                )

            for position in special_distances:
                band_ax.axvline(
                    position,
                    color='0.65',
                    linewidth=0.7,
                    linestyle='--',
                )

            band_ax.axhline(
                0.0,
                color='black',
                linewidth=0.8,
            )
            band_ax.set_xlim(
                distances[0],
                distances[-1],
            )
            band_ax.set_xticks(
                special_distances
            )
            band_ax.set_xticklabels(
                display_labels
            )
            band_ax.set_xlabel(
                'Wave vector'
            )
            band_ax.set_ylabel(
                'Frequency (THz)'
            )

            dos_ax.plot(
                total_dos,
                dos_frequencies,
                color='tab:red',
                linewidth=1.2,
            )
            dos_ax.fill_betweenx(
                dos_frequencies,
                0.0,
                total_dos,
                color='tab:red',
                alpha=0.2,
            )
            dos_ax.axhline(
                0.0,
                color='black',
                linewidth=0.8,
            )
            dos_ax.set_xlim(
                left=0.0
            )
            dos_ax.set_xlabel(
                'Phonon DOS (1/THz)'
            )
            dos_ax.tick_params(
                axis='y',
                labelleft=False,
            )

            fig.savefig(
                output_file,
                dpi=300,
                bbox_inches='tight',
            )
        finally:
            plt.close(fig)

        return output_file

    def phononcalc(self):
        """Run the phonon workflow using the selected DFT engine."""
        if self.Engine == 'GPAW':
            return self._phononcalc_gpaw()

        if self.Engine == 'QE':
            return self._phononcalc_qe()

        raise ValueError(
            f"Unsupported phonon engine: {self.Engine}"
        )

    def _phononcalc_qe(self):
        """Run the native Quantum ESPRESSO DFPT phonon workflow."""
        time51 = time.time()

        parprint(
            "Starting native QE phonon calculations..."
        )

        if self.Mode != 'PW':
            raise ValueError(
                "Quantum ESPRESSO phonon calculations support "
                "PW mode only."
            )

        self.engine.validate_qe_xc(
            self.XC_calc,
            pseudo_xc='pbe',
        )

        state_dir = Path(
            self.struct
            + '-GROUND-QE-Result-State'
        )

        if not self.engine.has_qe_state(
            state_dir,
            prefix='nanoworks',
        ):
            raise FileNotFoundError(
                f"{state_dir} does not contain a valid QE "
                "ground-state result. Complete the ground-state "
                "calculation before running phonons."
            )

        electronic_overrides = {
            'Phonon_PW_cutoff': self.Phonon_PW_cutoff,
            'Phonon_kpts_x': self.Phonon_kpts_x,
            'Phonon_kpts_y': self.Phonon_kpts_y,
            'Phonon_kpts_z': self.Phonon_kpts_z,
        }
        supplied_overrides = {
            name: value
            for name, value in electronic_overrides.items()
            if value is not None
        }

        if supplied_overrides:
            formatted_overrides = ', '.join(
                f"{name}={value}"
                for name, value in supplied_overrides.items()
            )
            parprint(
                "NOTICE: Native QE DFPT reuses the converged "
                "ground-state cutoff and electronic k-point mesh; "
                "the following GPAW finite-displacement settings "
                "are ignored: "
                + formatted_overrides
            )

        if self.Phonon_displacement != 1.0e-3:
            parprint(
                "NOTICE: Phonon_displacement is ignored by native "
                "QE DFPT because no displaced supercells are used."
            )

        qpoint_grid = (
            self.engine.resolve_qe_phonon_qpoint_grid(
                self.Phonon_supercell
            )
        )
        dos_qpoint_grid = (
            self.Phonon_qpts_x,
            self.Phonon_qpts_y,
            self.Phonon_qpts_z,
        )
        band_path = self.engine.build_band_path(
            self.bulk_configuration,
            path=self.Phonon_path,
            npoints=self.Phonon_npoints,
        )

        fildyn = Path(
            self.struct
            + '-PHONON-QE-Result-Dynamical-Matrix'
        )
        flfrc = Path(
            self.struct
            + '-PHONON-QE-Result-Force-Constants.fc'
        )
        flfrq = Path(
            self.struct
            + '-PHONON-QE-Result-Band.freq'
        )
        fldos = Path(
            self.struct
            + '-PHONON-QE-Result-DOS.dat'
        )

        ph_workflow = self.engine.run_ph(
            input_file=Path(
                self.struct
                + '-PHONON-QE-Input-PH.in'
            ),
            output_file=Path(
                self.struct
                + '-PHONON-QE-Log-PH.txt'
            ),
            state_dir=state_dir,
            fildyn=fildyn,
            qpoint_grid=qpoint_grid,
            parallel_cores=self.parallel_cores,
            executable='ph.x',
            prefix='nanoworks',
        )

        q2r_workflow = self.engine.run_q2r(
            input_file=Path(
                self.struct
                + '-PHONON-QE-Input-Q2R.in'
            ),
            output_file=Path(
                self.struct
                + '-PHONON-QE-Log-Q2R.txt'
            ),
            fildyn=fildyn,
            flfrc=flfrc,
            zasr='no',
            parallel_cores=self.parallel_cores,
            executable='q2r.x',
        )

        band_workflow = self.engine.run_matdyn_band(
            input_file=Path(
                self.struct
                + '-PHONON-QE-Input-Matdyn-Band.in'
            ),
            output_file=Path(
                self.struct
                + '-PHONON-QE-Log-Matdyn-Band.txt'
            ),
            flfrc=flfrc,
            flfrq=flfrq,
            band_path=band_path,
            acoustic_sum_rule=self.Phonon_acoustic_sum_rule,
            parallel_cores=self.parallel_cores,
            executable='matdyn.x',
        )

        dos_workflow = self.engine.run_matdyn_dos(
            input_file=Path(
                self.struct
                + '-PHONON-QE-Input-Matdyn-DOS.in'
            ),
            output_file=Path(
                self.struct
                + '-PHONON-QE-Log-Matdyn-DOS.txt'
            ),
            flfrc=flfrc,
            fldos=fldos,
            qpoint_grid=dos_qpoint_grid,
            acoustic_sum_rule=self.Phonon_acoustic_sum_rule,
            parallel_cores=self.parallel_cores,
            executable='matdyn.x',
        )

        band_data_file = self.engine.write_matdyn_band_data(
            output_file=Path(
                self.struct
                + '-PHONON-QE-Result-Band-THz.dat'
            ),
            band_path=band_path,
            frequencies=band_workflow[
                'frequencies'
            ],
        )
        dos_data_file = self.engine.write_matdyn_dos_data(
            output_file=Path(
                self.struct
                + '-PHONON-QE-Result-DOS-THz.dat'
            ),
            dos_data=dos_workflow[
                'dos'
            ],
        )
        graph_file = self._plot_qe_phonon_results(
            output_file=Path(
                self.struct
                + '-PHONON-QE-Graph-Phonon.png'
            ),
            band_path=band_path,
            frequencies=band_workflow[
                'frequencies'
            ],
            dos_data=dos_workflow[
                'dos'
            ],
        )

        thermal_data = None
        thermal_data_file = None

        if self.Phonon_thermal_calc:
            thermal_data = (
                self.engine.calculate_phonon_thermal_properties(
                    dos_workflow[
                        'dos'
                    ],
                    t_min=self.Phonon_T_min,
                    t_max=self.Phonon_T_max,
                    t_step=self.Phonon_T_step,
                )
            )
            thermal_data_file = (
                self.engine.write_phonon_thermal_properties(
                    output_file=Path(
                        self.struct
                        + '-PHONON-QE-Result-Thermal-Properties.csv'
                    ),
                    thermal_data=thermal_data,
                )
            )
            parprint(
                "QE phonon thermal properties written to: "
                + str(thermal_data_file)
            )
            parprint(
                "QE positive-frequency mode weight used for "
                "thermal integration: "
                f"{thermal_data['integrated_mode_weight']:.8f}"
            )

            if thermal_data['excluded_mode_weight'] > 1.0e-8:
                parprint(
                    "NOTICE: QE phonon thermal integration "
                    "excluded non-positive mode weight: "
                    f"{thermal_data['excluded_mode_weight']:.8f}"
                )

        time52 = time.time()

        with paropen(
            self.struct
            + '-TIMINGS-QE-Log-Timings.txt',
            'a',
        ) as fd:
            print(
                'Phonon calculation: ',
                round(time52 - time51, 2),
                file=fd,
            )

        parprint(
            "Native QE phonon calculations finished."
        )

        return {
            'qpoint_grid': qpoint_grid,
            'dos_qpoint_grid': dos_qpoint_grid,
            'band_path': band_path,
            'ph': ph_workflow,
            'q2r': q2r_workflow,
            'band': band_workflow,
            'dos': dos_workflow,
            'band_data_file': band_data_file,
            'dos_data_file': dos_data_file,
            'graph_file': graph_file,
            'thermal_data': thermal_data,
            'thermal_data_file': thermal_data_file,
        }

    def _phononcalc_gpaw(self):
        """
        This method performs a phonon calculation for the given structure using the ground state results. 
        It generates atomic displacements, computes force constants, and calculates phonon dispersion and phonon DOS.
        The results are saved as PNG file for now.
        """
        
        from phonopy import Phonopy
        from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
        import phonopy
        
        # -------------------------------------------------------------
        # PHONON CALCULATION
        # -------------------------------------------------------------

        time51 = time.time()
        parprint("Starting phonon calculations.")

        if self.engine.is_hybrid(self.XC_calc):
            parprint("\033[93mWARNING:\033[0m Phonon calculations use finite-difference forces; hybrid ("+self.XC_calc+") forces are expensive and unreliable in plane-wave GPAW.")
            parprint("It is recommended to compute phonons with PBE and use hybrids only for the electronic structure.")
            sys.exit(1)

        calc = self.engine.load_gpaw_calc(self.struct+'-GROUND-GPAW-Result-State.gpw')
        self.bulk_configuration.calc = calc

        # Pre-process
        bulk_configuration_ph = convert_atoms_to_phonopy(self.bulk_configuration)
        phonon = Phonopy(bulk_configuration_ph, self.Phonon_supercell, log_level=1)
        phonon.generate_displacements(distance=self.Phonon_displacement)
        with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
            print("[Phonopy] Atomic displacements:", end="\n", file=f2)
            disps = phonon.displacements
            for d in disps:
                print("[Phonopy] %d %s" % (d[0], d[1:]), end="\n", file=f2)

        calc = self.engine.create_phonon_calc(
            cutoff=self.Phonon_PW_cutoff,
            kpoint_size=(
                self.Phonon_kpts_x,
                self.Phonon_kpts_y,
                self.Phonon_kpts_z,
            ),
            txt=self.struct+'-PHONON-GPAW-Log-Calculation.txt',
        )

        self.bulk_configuration.calc = calc

        path = get_band_path(self.bulk_configuration, self.Phonon_path, self.Phonon_npoints)

        phonon_path = self.struct+'-PHONON-GPAW-Result-Force-Constants.npy'
        sum_rule=self.Phonon_acoustic_sum_rule

        if os.path.exists(phonon_path):
            with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
                print('Reading FCs from {!r}'.format(phonon_path), end="\n", file=f2)
            phonon.force_constants = np.load(phonon_path)

        else:
            with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
                print('Computing FCs',end="\n", file=f2)
                #os.makedirs('force-sets', exist_ok=True)
            supercells = list(phonon.supercells_with_displacements)
            fnames = [self.struct+'-PHONON-GPAW-Result-Supercell-{:04}.npy'.format(i) for i in range(len(supercells))]
            set_of_forces = [
                self.load_or_compute_force(fname, calc, supercell)
                for (fname, supercell) in zip(fnames, supercells)
            ]
            with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
                print('Building FC matrix', end="\n", file=f2)
            phonon.produce_force_constants(forces=set_of_forces, calculate_full_force_constants=False)
            if sum_rule:
                phonon.symmetrize_force_constants()
            with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
                print('Writing FCs to {!r}'.format(phonon_path), end="\n", file=f2)
            np.save(phonon_path, phonon.force_constants)
            #shutil.rmtree('force-sets')

        with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
            print('', end="\n", file=f2)
            print("[Phonopy] Phonon frequencies at Gamma:", end="\n", file=f2)
            for i, freq in enumerate(phonon.get_frequencies((0, 0, 0))):
                print("[Phonopy] %3d: %10.5f THz" %  (i + 1, freq), end="\n", file=f2) # THz

            # DOS
            print("[Phonopy] Initializing mesh...", end="\n", file=f2)
            phonon.init_mesh([self.Phonon_qpts_x, self.Phonon_qpts_y, self.Phonon_qpts_z])
            print("[Phonopy] Running total DOS calculation...", end="\n", file=f2)
            phonon.run_total_dos()
            print("[Phonopy] DOS calculation completed. Type of total_dos:", type(phonon.total_dos), end="\n", file=f2)
            if phonon.total_dos is not None:
                dos_array = np.array(phonon.total_dos)
                print("[Phonopy] DOS array shape:", dos_array.shape, "ndim:", dos_array.ndim, end="\n", file=f2)
            print('', end="\n", file=f2)
            print("[Phonopy] Phonon DOS:", end="\n", file=f2)
            
            # Check if total_dos is properly calculated and has the expected structure
            try:
                if phonon.total_dos is not None:
                    # Check if total_dos is a TotalDos object (new Phonopy format)
                    if hasattr(phonon.total_dos, 'frequency_points') and hasattr(phonon.total_dos, 'dos'):
                        frequencies = phonon.total_dos.frequency_points
                        dos_values = phonon.total_dos.dos
                        print("[Phonopy] Found TotalDos object with %d frequency points" % len(frequencies), end="\n", file=f2)
                        for omega, dos in zip(frequencies, dos_values):
                            print("%15.7f%15.7f" % (omega, dos), end="\n", file=f2)
                    # Check if total_dos is a tuple (freq, dos) as in older Phonopy versions
                    elif isinstance(phonon.total_dos, (tuple, list)) and len(phonon.total_dos) == 2:
                        frequencies, dos_values = phonon.total_dos
                        freq_array = np.array(frequencies)
                        dos_array = np.array(dos_values)
                        if freq_array.ndim == 1 and dos_array.ndim == 1 and len(freq_array) == len(dos_array):
                            for omega, dos in zip(freq_array, dos_array):
                                print("%15.7f%15.7f" % (omega, dos), end="\n", file=f2)
                        else:
                            print("[Phonopy] Warning: DOS frequency/values arrays have incompatible shapes", end="\n", file=f2)
                            print("[Phonopy] freq shape: %s, dos shape: %s" % (freq_array.shape, dos_array.shape), end="\n", file=f2)
                    else:
                        # Try the old method in case the format is different
                        dos_array = np.array(phonon.total_dos)
                        print("[Phonopy] DOS array shape:", dos_array.shape, "ndim:", dos_array.ndim, end="\n", file=f2)
                        if dos_array.ndim >= 2:  # Check if it's at least 2D
                            for omega, dos in dos_array.T:
                                print("%15.7f%15.7f" % (omega, dos), end="\n", file=f2)
                        else:
                            print("[Phonopy] Warning: DOS data has unexpected structure (ndim=%d)" % dos_array.ndim, end="\n", file=f2)
                            print("[Phonopy] DOS calculation may have failed or incomplete", end="\n", file=f2)
                            print("[Phonopy] Available attributes: %s" % [attr for attr in dir(phonon.total_dos) if not attr.startswith('_')], end="\n", file=f2)
                else:
                    print("[Phonopy] Warning: DOS data is None, calculation may have failed", end="\n", file=f2)
            except Exception as e:
                print("[Phonopy] Error processing DOS data: %s" % str(e), end="\n", file=f2)
                print("[Phonopy] Skipping DOS output to log file", end="\n", file=f2)

        qpoints, labels, connections = path
        phonon.run_band_structure(qpoints, path_connections=connections, labels=labels)

        # without DOS
        # fig = phonon.plot_band_structure()

        # with DOS
        #phonon.run_mesh([self.Phonon_qpts_x, self.Phonon_qpts_y, self.Phonon_qpts_z])
        phonon.run_total_dos()
        fig = phonon.plot_band_structure_and_dos()

        # with PDOS
        # phonon.run_mesh([20, 20, 20], with_eigenvectors=True, is_mesh_symmetry=False)
        # fig = phonon.plot_band_structure_and_dos(pdoc_indices=[[0], [1]])

        fig.savefig(self.struct+'-PHONON-GPAW-Graph-Phonon.png', dpi=300)
        # Standard Phonopy YAML format export
        try:
            phonon.write_yaml_band_structure(filename=self.struct+'-PHONON-GPAW-Result-Band.yaml')
        except Exception as e:
            parprint("YAML band data can not be saved:", e)

        # Saving in .dat format
        try:
            band_dict = phonon.get_band_structure_dict()
            distances = band_dict['distances']
            frequencies = band_dict['frequencies']
            
            with paropen(self.struct+'-PHONON-GPAW-Result-Band.dat', 'w') as f_band:
                f_band.write("Distance(1/A)    Frequencies(THz)...\n")
                # Return for every k-way segment
                for dist_path, freq_path in zip(distances, frequencies):
                    for d, f in zip(dist_path, freq_path):
                        # Write all frequency branches at the same q point
                        freq_str = "    ".join([f"{x:.6f}" for x in f])
                        f_band.write(f"{d:.6f}    {freq_str}\n")
                    f_band.write("\n") # Give space between segments
        except Exception as e:
            parprint("DAT band data can not be saved:", e)
        #Thermodynamic calculations
        if self.Phonon_thermal_calc:
            parprint(f"Thermal properties calculation (T: {self.Phonon_T_min}K - {self.Phonon_T_max}K)...")
                      
            # Termal properties calc
            phonon.run_thermal_properties(t_min=self.Phonon_T_min, 
                                          t_max=self.Phonon_T_max, 
                                          t_step=self.Phonon_T_step)
            
            # Phonopy standard YAML output
            phonon.write_yaml_thermal_properties(filename=self.struct+'-PHONON-GPAW-Result-Thermal-Properties.yaml')
            
            # convert to CSV format
            tp_dict = phonon.get_thermal_properties_dict()
            temperatures = tp_dict['temperatures']
            free_energy = tp_dict['free_energy']
            entropy = tp_dict['entropy']
            heat_capacity = tp_dict['heat_capacity']
            
            with paropen(self.struct+"-PHONON-GPAW-Result-Thermal-Properties.csv", "w") as f:
                f.write("T(K),Free_Energy(kJ/mol),Entropy(J/K/mol),Cv(J/K/mol)\n")
                for i in range(len(temperatures)):
                    f.write(f"{temperatures[i]:.2f},{free_energy[i]:.6f},{entropy[i]:.6f},{heat_capacity[i]:.6f}\n")
                    
            parprint(f"Thermal calculations finished!")
        
        time52 = time.time()
        # Write timings of calculation
        with paropen(self.struct+f'-TIMINGS-{self.Engine}-Log-Timings.txt', 'a') as f1:
            print('Phonon calculation: ', round((time52-time51),2), end="\n", file=f1)

    def opticalcalc(self):
        """Run the optical workflow using the selected DFT engine."""
        if self.Engine == 'GPAW':
            return self._opticalcalc_gpaw()

        if self.Engine == 'QE':
            return self._opticalcalc_qe()

        raise ValueError(
            f"Unsupported optical engine: {self.Engine}"
        )

    def _opticalcalc_qe(self):
        """Run the native QE NSCF and epsilon.x optical workflow."""
        time61 = time.time()

        parprint(
            "Starting native QE RPA optical calculation..."
        )

        if self.Mode != 'PW':
            raise ValueError(
                "Quantum ESPRESSO optical calculations support "
                "PW mode only."
            )

        if self.SOC_calc:
            raise NotImplementedError(
                "Quantum ESPRESSO SOC optical calculations are "
                "not supported yet."
            )

        calculation_type = str(
            self.Opt_calc_type
        ).strip().upper()

        if calculation_type != 'RPA':
            raise NotImplementedError(
                "Native Quantum ESPRESSO optical calculations "
                "currently support Opt_calc_type = 'RPA' only; "
                "epsilon.x does not provide the BSE workflow."
            )

        self.engine.validate_qe_xc(
            self.XC_calc,
            pseudo_xc='pbe',
        )

        state_dir = Path(
            self.struct
            + '-GROUND-QE-Result-State'
        )

        if not self.engine.has_qe_state(
            state_dir,
            prefix='nanoworks',
        ):
            raise FileNotFoundError(
                f"{state_dir} does not contain a valid QE "
                "ground-state result. Complete the ground-state "
                "calculation before running optics."
            )

        pseudo_dir = get_qe_pseudo_dir(
            relativistic='scalar',
        )
        pseudopotentials = resolve_qe_pseudopotentials(
            self.bulk_configuration,
            relativistic='scalar',
        )
        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )
        (
            optical_kpoint_density,
            optical_kpoint_size,
            optical_gamma,
        ) = resolve_stage_kpoint_settings(
            stage_density=self.Opt_kpts_density,
            stage_size=(
                self.Opt_kpts_x,
                self.Opt_kpts_y,
                self.Opt_kpts_z,
            ),
            stage_gamma=self.Opt_gamma,
            ground_density=self.Ground_kpts_density,
            ground_size=(
                self.Ground_kpts_x,
                self.Ground_kpts_y,
                self.Ground_kpts_z,
            ),
            ground_gamma=ground_gamma,
        )
        magnetic_moments = None

        if self.Spin_calc:
            magnetic_moments = resolve_initial_magnetic_moments(
                atoms=self.bulk_configuration,
                magmom_per_atom=self.Magmom_per_atom,
                magmom_single_atom=self.Magmom_single_atom,
            )

        nscf_workflow = self.engine.run_nscf(
            atoms=self.bulk_configuration,
            input_file=Path(
                self.struct
                + '-OPTICAL-QE-Input-NSCF.in'
            ),
            output_file=Path(
                self.struct
                + '-OPTICAL-QE-Log-NSCF.txt'
            ),
            state_dir=state_dir,
            pseudopotentials=pseudopotentials,
            pseudo_dir=pseudo_dir,
            cutoff_ev=self.Cut_off_energy,
            kpoint_density=optical_kpoint_density,
            kpoint_size=optical_kpoint_size,
            gamma=optical_gamma,
            total_charge=self.Total_charge,
            nbands=self.Opt_num_of_bands,
            spinpol=self.Spin_calc,
            magnetic_moments=magnetic_moments,
            setup_params=self.Setup_params,
            xc_calc=self.XC_calc,
            exx_fraction=self.XC_exx_fraction,
            omega=self.XC_omega,
            occupation={
                'name': 'fermi-dirac',
                'width': self.Opt_FD_smearing,
            },
            nosym=True,
            parallel_cores=self.parallel_cores,
            executable='pw.x',
            prefix='nanoworks',
        )

        parprint(
            "QE optical NSCF calculation finished."
        )

        epsilon_workflow = self.engine.run_epsilon(
            input_file=Path(
                self.struct
                + '-OPTICAL-QE-Input-Epsilon.in'
            ),
            output_file=Path(
                self.struct
                + '-OPTICAL-QE-Log-Epsilon.txt'
            ),
            state_dir=state_dir,
            result_dir=Path(
                self.struct
                + '-OPTICAL-QE-Result-Raw'
            ),
            calculation='eps',
            smeartype='gauss',
            intersmear=self.Opt_eta,
            intrasmear=0.0,
            wmin=self.Opt_min_en,
            wmax=self.Opt_max_en,
            nw=self.Opt_num_of_data,
            shift=self.Opt_shift_en,
            parallel_cores=self.parallel_cores,
            executable='epsilon.x',
            prefix='nanoworks',
        )

        parprint(
            "QE epsilon.x calculation finished."
        )

        optical_data = epsilon_workflow[
            'optical_data'
        ]
        table_files = self.engine.write_epsilon_optical_data(
            optical_data,
            Path(
                self.struct
                + '-OPTICAL-QE-Result-Calculation-RPA'
            ),
        )

        for direction, table_file in table_files.items():
            parprint(
                "QE optical "
                + direction
                + "-direction data saved to: "
                + str(table_file)
            )

        figure_files = {}

        if world.rank == 0:
            for direction, table_file in table_files.items():
                plot_data = np.loadtxt(
                    table_file,
                    skiprows=1,
                )
                figure_prefix = (
                    self.struct
                    + '-OPTICAL-QE-Graph-RPA-'
                    + direction
                    + 'direction'
                )
                self._generate_optical_figures(
                    plot_data,
                    figure_prefix,
                    f"QE RPA ({direction})",
                )
                figure_files[direction] = figure_prefix
                parprint(
                    "QE optical figures saved with prefix: "
                    + figure_prefix
                )

        time62 = time.time()

        with paropen(
            self.struct
            + f'-TIMINGS-{self.Engine}-Log-Timings.txt',
            'a',
        ) as fd:
            print(
                'Optical calculation: ',
                round(time62 - time61, 2),
                file=fd,
            )

        return {
            'nscf': nscf_workflow,
            'epsilon': epsilon_workflow,
            'table_files': table_files,
            'figure_prefixes': figure_files,
        }

    def _opticalcalc_gpaw(self):
        """
        This method performs optical property calculations for the given structure using the
        ground state results. It computes the dielectric function, absorption spectrum, and
        other related optical properties. The results are saved in appropriate files for
        further analysis and visualization.
        """

        # -------------------------------------------------------------
        # OPTICAL CALCULATION
        # -------------------------------------------------------------

        #Start Optical calc
        time61 = time.time()
        
        ground_gamma = (
            self.Gamma
            if self.Ground_gamma is None
            else self.Ground_gamma
        )
        
        if self.Mode == 'PW':
            parprint("Starting optical calculation...")
            try:
                hybrid = self.engine.is_hybrid(self.XC_calc)
                
                opt_kpoint_density, opt_kpoint_size, opt_gamma = (
                    resolve_stage_kpoint_settings(
                        stage_density=self.Opt_kpts_density,
                        stage_size=(
                            self.Opt_kpts_x,
                            self.Opt_kpts_y,
                            self.Opt_kpts_z,
                        ),
                        stage_gamma=self.Opt_gamma,
                        ground_density=self.Ground_kpts_density,
                        ground_size=(
                            self.Ground_kpts_x,
                            self.Ground_kpts_y,
                            self.Ground_kpts_z,
                        ),
                        ground_gamma=ground_gamma,
                    )
                )

                if hybrid:
                    parprint(
                        "Hybrid XC detected: reading ground state directly "
                        "for optical calculation..."
                    )

                calc = self.engine.prepare_optical_calc(
                    filename=self.struct+'-GROUND-GPAW-Result-State.gpw',
                    hybrid=hybrid,
                    txt=self.struct+'-OPTICAL-GPAW-Log-Calculation.txt',
                    nbands=self.Opt_num_of_bands,
                    smearing=self.Opt_FD_smearing,
                    kpoint_density=opt_kpoint_density,
                    kpoint_size=opt_kpoint_size,
                    gamma=opt_gamma,
                )
            except FileNotFoundError as err:
                # output error, and return with an error code
                parprint('\033[91mERROR:\033[0mOptical computations must be done separately. Please do ground calculations first.')
                sys.exit(1)

            calc.get_potential_energy()

            calc.diagonalize_full_hamiltonian(nbands=self.Opt_num_of_bands)  # diagonalize Hamiltonian
            calc.write(self.struct+'-OPTICAL-GPAW-Result-State.gpw', mode= 'all')  # write wavefunctions

            #from mpi4py import MPI
            if self.Opt_calc_type == 'BSE':
                from gpaw.response.bse import BSE
                
                if self.Spin_calc == True:
                   parprint('\033[91mERROR:\033[0mBSE calculations can not run with spin dependent data.')
                   sys.exit(1)
                parprint('Starting BSE calculations')
                bse = BSE(calc= self.struct+'-OPTICAL-GPAW-Result-State.gpw', ecut=self.Opt_cut_of_energy,
                             valence_bands=self.Opt_BSE_valence,
                             conduction_bands=self.Opt_BSE_conduction,
                             nbands=self.Opt_num_of_bands,
                             mode='BSE',
                             integrate_gamma='sphere', txt=self.struct+'-OPTICAL-GPAW-Log-Calculation-BSE.txt')

                # Getting dielectric function spectrum
                parprint("Starting dielectric function calculation...")
                # Writing to files
                bse.get_dielectric_function(filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec.csv',
                                            eta=self.Opt_eta, w_w=np.linspace(self.Opt_min_en, self.Opt_max_en, self.Opt_num_of_data),
                                            write_eig=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_eig.dat')
                # Loading dielectric function spectrum to numpy
                dielec = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec.csv', delimiter=',')
                # dielec.shape[0] will give us the length of data.
                c_opt = 29979245800
                h_opt = 6.58E-16
                # Initialize arrays
                opt_n_bse = np.array ([1e-6,]*dielec.shape[0])
                opt_k_bse = np.array ([1e-6,]*dielec.shape[0])
                opt_abs_bse = np.array([1e-6,]*dielec.shape[0])
                opt_ref_bse = np.array([1e-6,]*dielec.shape[0])
                # Calculation of other optical data
                for n in range(dielec.shape[0]):
                    opt_n_bse[n] = np.sqrt((np.sqrt(np.square(dielec[n][1])+np.square(dielec[n][2]))+dielec[n][1])/2.0)
                    opt_k_bse[n] = np.sqrt((np.sqrt(np.square(dielec[n][1])+np.square(dielec[n][2]))-dielec[n][1])/2.0)
                    opt_abs_bse[n] = 2*dielec[n][0]*opt_k_bse[n]/(h_opt*c_opt)
                    opt_ref_bse[n] = (np.square(1-opt_n_bse[n])+np.square(opt_k_bse[n]))/(np.square(1+opt_n_bse[n])+np.square(opt_k_bse[n]))
                
                # Saving other data
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE-AllData.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec.shape[0]):
                        print(dielec[n][0], dielec[n][1], dielec[n][2], opt_n_bse[n], opt_k_bse[n], opt_abs_bse[n], opt_ref_bse[n], end="\n", file=f1)
                    print (end="\n", file=f1)
                    
                '''
                # DIRECTION IS NOT WORKING FOR A WHILE, IN FUTURE THESE LINES CAN BE USED
                bse.get_dielectric_function(filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec_xdirection.csv',
                                            q_c = [0.0, 0.0, 0.0], direction=0, eta=self.Opt_eta,
                                            w_w=np.linspace(self.Opt_min_en, self.Opt_max_en, self.Opt_num_of_data),
                                            write_eig=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_eig_xdirection.dat')
                bse.get_dielectric_function(q_c = [0.0, 0.0, 0.0], direction=1, eta=self.Opt_eta,
                                            w_w=np.linspace(self.Opt_min_en, self.Opt_max_en, self.Opt_num_of_data),
                                            filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec_ydirection.csv',
                                            write_eig=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_eig_ydirection.dat')
                bse.get_dielectric_function(q_c = [0.0, 0.0, 0.0], direction=2, eta=self.Opt_eta,
                                            w_w=np.linspace(self.Opt_min_en, self.Opt_max_en, self.Opt_num_of_data),
                                            filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec_zdirection.csv',
                                            write_eig=self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_eig_zdirection.dat')

                # Loading dielectric function spectrum to numpy
                dielec_x = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec_xdirection.csv', delimiter=',')
                dielec_y = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec_ydirection.csv', delimiter=',')
                dielec_z = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE_dielec_zdirection.csv', delimiter=',')
                # dielec_x.shape[0] will give us the length of data.
                # c and h
                c_opt = 29979245800
                h_opt = 6.58E-16
                #c_opt = 1
                #h_opt = 1

                # Initialize arrays
                opt_n_bse_x = np.array ([1e-6,]*dielec_x.shape[0])
                opt_k_bse_x = np.array ([1e-6,]*dielec_x.shape[0])
                opt_abs_bse_x = np.array([1e-6,]*dielec_x.shape[0])
                opt_ref_bse_x = np.array([1e-6,]*dielec_x.shape[0])
                opt_n_bse_y = np.array ([1e-6,]*dielec_y.shape[0])
                opt_k_bse_y = np.array ([1e-6,]*dielec_y.shape[0])
                opt_abs_bse_y = np.array([1e-6,]*dielec_y.shape[0])
                opt_ref_bse_y = np.array([1e-6,]*dielec_y.shape[0])
                opt_n_bse_z = np.array ([1e-6,]*dielec_z.shape[0])
                opt_k_bse_z = np.array ([1e-6,]*dielec_z.shape[0])
                opt_abs_bse_z = np.array([1e-6,]*dielec_z.shape[0])
                opt_ref_bse_z = np.array([1e-6,]*dielec_z.shape[0])

                # Calculation of other optical data
                for n in range(dielec_x.shape[0]):
                    # x-direction
                    opt_n_bse_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))+dielec_x[n][1])/2.0)
                    opt_k_bse_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))-dielec_x[n][1])/2.0)
                    opt_abs_bse_x[n] = 2*dielec_x[n][0]*opt_k_bse_x[n]/(h_opt*c_opt)
                    opt_ref_bse_x[n] = (np.square(1-opt_n_bse_x[n])+np.square(opt_k_bse_x[n]))/(np.square(1+opt_n_bse_x[n])+np.square(opt_k_bse_x[n]))
                    # y-direction
                    opt_n_bse_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))+dielec_y[n][1])/2.0)
                    opt_k_bse_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))-dielec_y[n][1])/2.0)
                    opt_abs_bse_y[n] = 2*dielec_y[n][0]*opt_k_bse_y[n]/(h_opt*c_opt)
                    opt_ref_bse_y[n] = (np.square(1-opt_n_bse_y[n])+np.square(opt_k_bse_y[n]))/(np.square(1+opt_n_bse_y[n])+np.square(opt_k_bse_y[n]))
                    # z-direction
                    opt_n_bse_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))+dielec_z[n][1])/2.0)
                    opt_k_bse_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))-dielec_z[n][1])/2.0)
                    opt_abs_bse_z[n] = 2*dielec_z[n][0]*opt_k_bse_z[n]/(h_opt*c_opt)
                    opt_ref_bse_z[n] = (np.square(1-opt_n_bse_z[n])+np.square(opt_k_bse_z[n]))/(np.square(1+opt_n_bse_z[n])+np.square(opt_k_bse_z[n]))

                # Saving other data for x-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE-AllData_xdirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_x.shape[0]):
                        print(dielec_x[n][0], dielec_x[n][1], dielec_x[n][2], opt_n_bse_x[n], opt_k_bse_x[n], opt_abs_bse_x[n], opt_ref_bse_x[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # Saving other data for y-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE-AllData_ydirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_y.shape[0]):
                        print(dielec_y[n][0], dielec_y[n][1], dielec_y[n][2], opt_n_bse_y[n], opt_k_bse_y[n], opt_abs_bse_y[n], opt_ref_bse_y[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # Saving other data for z-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-BSE-AllData_zdirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_z.shape[0]):
                        print(dielec_z[n][0], dielec_z[n][1], dielec_z[n][2], opt_n_bse_z[n], opt_k_bse_z[n], opt_abs_bse_z[n], opt_ref_bse_z[n], end="\n", file=f1)
                    print (end="\n", file=f1)
               '''
            
            elif self.Opt_calc_type == 'RPA':
                parprint('Starting RPA calculations')
                from gpaw.response.df import DielectricFunction
                df = DielectricFunction(calc=self.struct+'-OPTICAL-GPAW-Result-State.gpw',
                                        frequencies={'type': 'nonlinear', 'domega0': self.Opt_domega0, 'omega2': self.Opt_omega2},
                                        eta=self.Opt_eta, intraband=False, hilbert=False,
                                        ecut=self.Opt_cut_of_energy, txt=self.struct+'-OPTICAL-GPAW-Log-Calculation-RPA.txt')
                # Writing to files as: omega, nlfc.real, nlfc.imag, lfc.real, lfc.imag
                # Here lfc is local field correction
                # Getting dielectric function spectrum
                parprint("Starting dielectric function calculation...")
                df.get_dielectric_function(direction='x', filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA_dielec_xdirection.csv')
                df.get_dielectric_function(direction='y', filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA_dielec_ydirection.csv')
                df.get_dielectric_function(direction='z', filename=self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA_dielec_zdirection.csv')

                # Loading dielectric function spectrum to numpy
                dielec_x = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA_dielec_xdirection.csv', delimiter=',')
                dielec_y = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA_dielec_ydirection.csv', delimiter=',')
                dielec_z = genfromtxt(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA_dielec_zdirection.csv', delimiter=',')
                # dielec_x.shape[0] will give us the length of data.
                # c and h
                c_opt = 29979245800
                h_opt = 6.58E-16
                #c_opt = 1
                #h_opt = 1
                # ---- NLFC ----
                # Initialize arrays for NLFC
                opt_n_nlfc_x = np.array ([1e-6,]*dielec_x.shape[0])
                opt_k_nlfc_x = np.array ([1e-6,]*dielec_x.shape[0])
                opt_abs_nlfc_x = np.array([1e-6,]*dielec_x.shape[0])
                opt_ref_nlfc_x = np.array([1e-6,]*dielec_x.shape[0])
                opt_n_nlfc_y = np.array ([1e-6,]*dielec_y.shape[0])
                opt_k_nlfc_y = np.array ([1e-6,]*dielec_y.shape[0])
                opt_abs_nlfc_y = np.array([1e-6,]*dielec_y.shape[0])
                opt_ref_nlfc_y = np.array([1e-6,]*dielec_y.shape[0])
                opt_n_nlfc_z = np.array ([1e-6,]*dielec_z.shape[0])
                opt_k_nlfc_z = np.array ([1e-6,]*dielec_z.shape[0])
                opt_abs_nlfc_z = np.array([1e-6,]*dielec_z.shape[0])
                opt_ref_nlfc_z = np.array([1e-6,]*dielec_z.shape[0])

                # Calculation of other optical spectrum for NLFC
                for n in range(dielec_x.shape[0]):
                    # NLFC-x
                    opt_n_nlfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))+dielec_x[n][1])/2.0)
                    opt_k_nlfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))-dielec_x[n][1])/2.0)
                    opt_abs_nlfc_x[n] = 2*dielec_x[n][0]*opt_k_nlfc_x[n]/(h_opt*c_opt)
                    opt_ref_nlfc_x[n] = (np.square(1-opt_n_nlfc_x[n])+np.square(opt_k_nlfc_x[n]))/(np.square(1+opt_n_nlfc_x[n])+np.square(opt_k_nlfc_x[n]))
                    # NLFC-y
                    opt_n_nlfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))+dielec_y[n][1])/2.0)
                    opt_k_nlfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))-dielec_y[n][1])/2.0)
                    opt_abs_nlfc_y[n] = 2*dielec_y[n][0]*opt_k_nlfc_y[n]/(h_opt*c_opt)
                    opt_ref_nlfc_y[n] = (np.square(1-opt_n_nlfc_y[n])+np.square(opt_k_nlfc_y[n]))/(np.square(1+opt_n_nlfc_y[n])+np.square(opt_k_nlfc_y[n]))
                    # NLFC-z
                    opt_n_nlfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))+dielec_z[n][1])/2.0)
                    opt_k_nlfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))-dielec_z[n][1])/2.0)
                    opt_abs_nlfc_z[n] = 2*dielec_z[n][0]*opt_k_nlfc_z[n]/(h_opt*c_opt)
                    opt_ref_nlfc_z[n] = (np.square(1-opt_n_nlfc_z[n])+np.square(opt_k_nlfc_z[n]))/(np.square(1+opt_n_nlfc_z[n])+np.square(opt_k_nlfc_z[n]))

                # Saving NLFC other optical spectrum for x-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA-NLFC-AllData_xdirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_x.shape[0]):
                        print(dielec_x[n][0], dielec_x[n][1], dielec_x[n][2], opt_n_nlfc_x[n], opt_k_nlfc_x[n], opt_abs_nlfc_x[n], opt_ref_nlfc_x[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # Saving NLFC other optical spectrum for y-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA-NLFC-AllData_ydirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_y.shape[0]):
                        print(dielec_y[n][0], dielec_y[n][1], dielec_y[n][2], opt_n_nlfc_y[n], opt_k_nlfc_y[n], opt_abs_nlfc_y[n], opt_ref_nlfc_y[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # Saving NLFC other optical spectrum for z-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA-NLFC-AllData_zdirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_z.shape[0]):
                        print(dielec_z[n][0], dielec_z[n][1], dielec_z[n][2], opt_n_nlfc_z[n], opt_k_nlfc_z[n], opt_abs_nlfc_z[n], opt_ref_nlfc_z[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # ---- LFC ----
                # Initialize arrays for LFC
                opt_n_lfc_x = np.array ([1e-6,]*dielec_x.shape[0])
                opt_k_lfc_x = np.array ([1e-6,]*dielec_x.shape[0])
                opt_abs_lfc_x = np.array([1e-6,]*dielec_x.shape[0])
                opt_ref_lfc_x = np.array([1e-6,]*dielec_x.shape[0])
                opt_n_lfc_y = np.array ([1e-6,]*dielec_y.shape[0])
                opt_k_lfc_y = np.array ([1e-6,]*dielec_y.shape[0])
                opt_abs_lfc_y = np.array([1e-6,]*dielec_y.shape[0])
                opt_ref_lfc_y = np.array([1e-6,]*dielec_y.shape[0])
                opt_n_lfc_z = np.array ([1e-6,]*dielec_z.shape[0])
                opt_k_lfc_z = np.array ([1e-6,]*dielec_z.shape[0])
                opt_abs_lfc_z = np.array([1e-6,]*dielec_z.shape[0])
                opt_ref_lfc_z = np.array([1e-6,]*dielec_z.shape[0])

                # Calculation of other optical spectrum for LFC
                for n in range(dielec_x.shape[0]):
                    # LFC-x
                    opt_n_lfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][3])+np.square(dielec_x[n][4]))+dielec_x[n][3])/2.0)
                    opt_k_lfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][3])+np.square(dielec_x[n][4]))-dielec_x[n][3])/2.0)
                    opt_abs_lfc_x[n] = 2*dielec_x[n][0]*opt_k_nlfc_x[n]/(h_opt*c_opt)
                    opt_ref_lfc_x[n] = (np.square(1-opt_n_lfc_x[n])+np.square(opt_k_lfc_x[n]))/(np.square(1+opt_n_lfc_x[n])+np.square(opt_k_lfc_x[n]))
                    # LFC-y
                    opt_n_lfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][3])+np.square(dielec_y[n][4]))+dielec_y[n][3])/2.0)
                    opt_k_lfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][3])+np.square(dielec_y[n][4]))-dielec_y[n][3])/2.0)
                    opt_abs_lfc_y[n] = 2*dielec_y[n][0]*opt_k_lfc_y[n]/(h_opt*c_opt)
                    opt_ref_lfc_y[n] = (np.square(1-opt_n_lfc_y[n])+np.square(opt_k_lfc_y[n]))/(np.square(1+opt_n_lfc_y[n])+np.square(opt_k_lfc_y[n]))
                    # LFC-z
                    opt_n_lfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][3])+np.square(dielec_z[n][4]))+dielec_z[n][3])/2.0)
                    opt_k_lfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][3])+np.square(dielec_z[n][4]))-dielec_z[n][3])/2.0)
                    opt_abs_lfc_z[n] = 2*dielec_z[n][0]*opt_k_lfc_z[n]/(h_opt*c_opt)
                    opt_ref_lfc_z[n] = (np.square(1-opt_n_lfc_z[n])+np.square(opt_k_lfc_z[n]))/(np.square(1+opt_n_lfc_z[n])+np.square(opt_k_lfc_z[n]))

                # Saving LFC other optical spectrum for x-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA-LFC-AllData_xdirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_x.shape[0]):
                        print(dielec_x[n][0], dielec_x[n][3], dielec_x[n][4], opt_n_lfc_x[n], opt_k_lfc_x[n], opt_abs_lfc_x[n], opt_ref_lfc_x[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # Saving LFC other optical spectrum for y-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA-LFC-AllData_ydirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_y.shape[0]):
                        print(dielec_y[n][0], dielec_y[n][3], dielec_y[n][4], opt_n_lfc_y[n], opt_k_lfc_y[n], opt_abs_lfc_y[n], opt_ref_lfc_y[n], end="\n", file=f1)
                    print (end="\n", file=f1)

                # Saving LFC other optical spectrum for z-direction
                with paropen(self.struct+'-OPTICAL-GPAW-Result-Calculation-RPA-LFC-AllData_zdirection.dat', 'w') as f1:
                    print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                    for n in range(dielec_z.shape[0]):
                        print(dielec_z[n][0], dielec_z[n][3], dielec_z[n][4], opt_n_lfc_z[n], opt_k_lfc_z[n], opt_abs_lfc_z[n], opt_ref_lfc_z[n], end="\n", file=f1)
                    print (end="\n", file=f1)

            else:
                parprint('\033[91mERROR:\033[0mUnknown optical calculation type.')
                sys.exit(1)

        elif self.Mode == 'LCAO':
            parprint('\033[91mERROR:\033[0mNot implemented in LCAO mode yet.')
        else:
            parprint('\033[91mERROR:\033[0mNot implemented in FD mode yet.')
        # Finish Optical calc
        time62 = time.time()
        # Write timings of calculation
        with paropen(self.struct+f'-TIMINGS-{self.Engine}-Log-Timings.txt', 'a') as f1:
            print('Optical calculation: ', round((time62-time61),2), end="\n", file=f1)
        
        # Plot generated optical data
        self._plot_optical_results()
    
    def _generate_optical_figures(self, data, file_prefix, title_suffix):
        """
        Helper method to generate 4 standard optical figures from parsed array data.
        Forces the 'Agg' backend for matplotlib to prevent X11 server errors on clusters.
        
        Args:
            data (numpy.ndarray): 2D array containing the optical data.
                                  Col 0: Energy, Col 1: Eps_real, Col 2: Eps_img,
                                  Col 3: Refractive_Index, Col 4: Extinction_Index,
                                  Col 5: Absorption, Col 6: Reflectivity
            file_prefix (str): Prefix string for the saved .png files.
            title_suffix (str): String appended to the plot titles (e.g., 'BSE' or 'RPA LFC (x)').
        """
        import matplotlib
        matplotlib.use('Agg') # Essential for headless cluster/HPC environments
        import matplotlib.pyplot as plt

        # Extract columns based on the .dat file structure generated in opticalcalc
        energy = data[:, 0]
        eps_real = data[:, 1]
        eps_img = data[:, 2]
        n_idx = data[:, 3]
        k_idx = data[:, 4]
        absorption = data[:, 5]
        reflectivity = data[:, 6]

        # ---------------------------------------------------------
        # 1. Plot Complex Dielectric Function (Epsilon Real & Imag)
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(energy, eps_real, label=r'$\epsilon_1$ (Re)', color='blue', linewidth=2)
        ax.plot(energy, eps_img, label=r'$\epsilon_2$ (Im)', color='red', linestyle='--', linewidth=2)
        ax.set_xlabel(self._t("fig_optical_xlabel"))
        ax.set_ylabel('Dielectric Function')
        ax.set_xlim(left=0) 
        ax.axhline(y=0, color='black', linestyle='-', linewidth=0.8) 
        ax.legend(loc='upper right', frameon=True)
        ax.grid(True, linestyle=':', alpha=0.7)
        fig.tight_layout()
        fig.savefig(f"{file_prefix}-Dielectric.png", dpi=300)
        plt.close(fig)

        # ---------------------------------------------------------
        # 2. Plot Complex Refractive Index (n & k)
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(energy, n_idx, label='n', color='darkgreen', linewidth=2)
        ax.plot(energy, k_idx, label='k', color='orange', linestyle='--', linewidth=2)
        ax.set_xlabel(self._t("fig_optical_xlabel"))
        ax.set_ylabel('n,k')
        ax.set_xlim(left=0)
        ax.legend(loc='upper right', frameon=True)
        ax.grid(True, linestyle=':', alpha=0.7)
        fig.tight_layout()
        fig.savefig(f"{file_prefix}-n-and-k.png", dpi=300)
        plt.close(fig)

        # ---------------------------------------------------------
        # 3. Plot Absorption Coefficient
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(energy, absorption, label=r'Abs ($\alpha$)', color='purple', linewidth=2)
        ax.set_xlabel(self._t("fig_optical_xlabel"))
        ax.set_ylabel(self._t("fig_optical_abs"))
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0) 
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
        ax.legend(loc='upper right', frameon=True)
        ax.grid(True, linestyle=':', alpha=0.7)
        fig.tight_layout()
        fig.savefig(f"{file_prefix}-Absorption.png", dpi=300)
        plt.close(fig)

        # ---------------------------------------------------------
        # 4. Plot Reflectivity
        # ---------------------------------------------------------
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.plot(energy, reflectivity, label='Ref (R)', color='teal', linewidth=2)
        ax.set_xlabel(self._t("fig_optical_xlabel"))
        ax.set_ylabel(self._t("fig_optical_ref"))
        ax.set_xlim(left=0)
        ax.set_ylim(bottom=0)
        ax.legend(loc='upper right', frameon=True)
        ax.grid(True, linestyle=':', alpha=0.7)
        fig.tight_layout()
        fig.savefig(f"{file_prefix}-Reflectivity.png", dpi=300)
        plt.close(fig)

    def _plot_optical_results(self):
        """
        Scans for generated optical .dat files and plots them.
        Execution is strictly isolated to the master node (rank 0).
        """
        from ase.parallel import world, parprint
        import os
        import numpy as np

        # Safety check: Ensure only the master core executes plotting
        if world.rank == 0:
            parprint("Generating optical property figures...")
            
            # 1. Process BSE Data
            bse_filename = f"{self.struct}-OPTICAL-GPAW-Result-Calculation-BSE-AllData.dat"
            if os.path.exists(bse_filename):
                try:
                    bse_data = np.loadtxt(bse_filename, skiprows=1)
                    # Create prefix aligning with Nanoworks standard
                    file_prefix = f"{self.struct}-OPTICAL-GPAW-Graph-BSE"
                    self._generate_optical_figures(bse_data, file_prefix, "BSE")
                except Exception as e:
                    print(f"Error plotting BSE data: {e}")

            # 2. Process RPA Data (LFC and NLFC across x, y, z directions)
            rpa_types = ["LFC", "NLFC"]
            directions = ["xdirection", "ydirection", "zdirection"]

            for rtype in rpa_types:
                for direction in directions:
                    rpa_filename = f"{self.struct}-OPTICAL-GPAW-Result-Calculation-RPA-{rtype}-AllData_{direction}.dat"
                    if os.path.exists(rpa_filename):
                        try:
                            rpa_data = np.loadtxt(rpa_filename, skiprows=1)
                            # Create prefix aligning with Nanoworks standard
                            file_prefix = f"{self.struct}-OPTICAL-GPAW-Graph-RPA-{rtype}-{direction}"
                            title_suffix = f"RPA {rtype} ({direction})"
                            self._generate_optical_figures(rpa_data, file_prefix, title_suffix)
                        except Exception as e:
                            print(f"Error plotting RPA data ({rtype}, {direction}): {e}")

    def run_gpaw(self, calc, cell):
        cell = convert_atoms_to_ase(cell)
        cell.set_calculator(calc)
        forces = cell.get_forces()
        drift_force = forces.sum(axis=0)
        with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
            print(("[Phonopy] Drift force:" + "%11.5f" * 3) % tuple(drift_force), end="\n", file=f2)
        # Simple translational invariance
        for force in forces:
            force -= drift_force / forces.shape[0]
        return forces

    def load_or_compute_force(self, path, calc, atoms):
        if os.path.exists(path):
            with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
                print('Reading {!r}'.format(path), end="\n", file=f2)
            return np.load(path)

        else:
            with paropen(self.struct+'-PHONON-GPAW-Log-Phonopy.txt', 'a') as f2:
                print('Computing {!r}'.format(path), end="\n", file=f2)
            force_set = self.run_gpaw(calc, atoms)
            np.save(path, force_set)
            return force_set

# Elastic related functions
from ase.spacegroup import get_spacegroup
import numpy as np

def reconstruct_full_tensor(independent, atoms):
    """
    Reconstruct the full 6x6 elastic tensor from the independent elastic constants
    returned by get_elastic_tensor(), based on the spacegroup determined by ASE.
    
    Parameters:
      independent : NumPy array of independent elastic constants.
      atoms       : ASE Atoms object.
      
    Returns:
      A 6x6 NumPy array representing the full elastic tensor.
      
    Implemented cases:
      - 3 independent constants: assumed cubic or isotropic, 
        mapped as:
            C11  C12  C12   0    0    0
            C12  C11  C12   0    0    0
            C12  C12  C11   0    0    0
             0    0    0   C44   0    0
             0    0    0    0   C44   0
             0    0    0    0    0   (C11-C12)/2
      - 5 independent constants: assumed hexagonal, mapped as:
            C11  C12  C13   0    0    0
            C12  C11  C13   0    0    0
            C13  C13  C33   0    0    0
             0    0    0   C44   0    0
             0    0    0    0   C44   0
             0    0    0    0    0   (C11-C12)/2
      - 21 independent constants: assumed triclinic; these are
        taken as the lower triangular elements (in order) and then mirrored.
    """
    sg = get_spacegroup(atoms, symprec=1e-2)
    symbol = sg.symbol
    parprint(f"[DEBUG] Reconstructing tensor from {independent.size} independent constants; spacegroup: {symbol}")
    if independent.size == 3:
        # Assume cubic or isotropic system.
        C11, C12, C44 = independent
        full = np.array([
            [C11, C12, C12, 0,   0,   0],
            [C12, C11, C12, 0,   0,   0],
            [C12, C12, C11, 0,   0,   0],
            [0,   0,   0,   C44, 0,   0],
            [0,   0,   0,   0,   C44, 0],
            [0,   0,   0,   0,   0,   (C11-C12)/2]
        ])
        return full
    elif independent.size == 5:
        # Assume hexagonal symmetry.
        C11, C12, C13, C33, C44 = independent
        C66 = (C11 - C12) / 2
        full = np.array([
            [C11, C12, C13, 0,   0,   0],
            [C12, C11, C13, 0,   0,   0],
            [C13, C13, C33, 0,   0,   0],
            [0,   0,   0,   C44, 0,   0],
            [0,   0,   0,   0,   C44, 0],
            [0,   0,   0,   0,   0,   C66]
        ])
        return full
    elif independent.size == 21:
        full = np.zeros((6,6))
        idx = 0
        for i in range(6):
            for j in range(i+1):
                full[i,j] = independent[idx]
                idx += 1
        # Mirror the lower triangle to the upper triangle.
        full = full + full.T - np.diag(np.diag(full))
        return full
    else:
        raise ValueError("Reconstruction for symmetry with {} independent constants is not implemented.".format(independent.size))


# Phonon related functions
# The remaining functions related to phonon calculations in this file are MIT-licensed by (C) 2020 Michael Lamparski

def get_band_path(atoms, path_str, npoints, path_frac=None, labels=None):
    from ase.dft.kpoints import bandpath

    atoms = convert_atoms_to_ase(atoms)
    if path_str is None:
        path_str = atoms.get_cell().bandpath().path

    # Commas are part of ase's supported syntax, but we'll take care of them
    # ourselves to make it easier to get things the way phonopy wants them
    if path_frac is None:
        path_frac = []
        for substr in path_str.split(','):
            path = bandpath(substr, atoms.get_cell()[...], npoints=1)
            path_frac.append(path.kpts)

    if labels is None:
        labels = []
        for substr in path_str.split(','):
            path = bandpath(substr, atoms.get_cell()[...], npoints=1)

            _, _, substr_labels = path.get_linear_kpoint_axis()
            labels.extend(['$\\Gamma$' if s == 'G' else s for s in substr_labels])

    qpoints, connections = get_band_qpoints_and_path_connections(path_frac, npoints=npoints)
    return qpoints, labels, connections

def convert_atoms_to_ase(atoms):
    if hasattr(atoms, 'get_chemical_symbols'):
        return Atoms(
            symbols=atoms.get_chemical_symbols(),
            scaled_positions=atoms.get_scaled_positions(),
            cell=atoms.get_cell(),
            pbc=True
        )
    else:
        return Atoms(
            symbols=atoms.symbols,
            scaled_positions=atoms.scaled_positions,
            cell=atoms.cell,
            pbc=True
        )

def convert_atoms_to_phonopy(atoms):
    from phonopy.structure.atoms import PhonopyAtoms

    return PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        scaled_positions=atoms.get_scaled_positions(),
        cell=atoms.get_cell()
    )


# End of phonon related functions------------------------------

# Projected Band Structure related functions-------------------

def projected_weights(calc):
    ns = calc.get_number_of_spins()
    atom_num = calc.atoms.get_atomic_numbers()
    atoms = calc.atoms
    
    # Defining the atoms and angular momentum to project onto.
    lan = range(58, 72)
    act = range(90, 104)
    atom_num = np.asarray(atom_num)
    ang_mom_a = {}
    atoms = Atoms(numbers=atom_num)
    magnetic_elements = {'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In',
                'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl'}

    for a, (z, magn) in enumerate(zip(atom_num, magnetic_elements)):
        if z in lan or z in act:
            ang_mom_a[a] = 'spdf'
        else:
            ang_mom_a[a] = 'spd' if magn else 'sp'

    # For each unique atom
    a_x = [a for a in ang_mom_a for _ in ang_mom_a[a]]
    ang_mom_x = [ang_mom for a in ang_mom_a for ang_mom in ang_mom_a[a]]

    # Get i index for each unique symbol
    sym_ang_mom_i = []
    i_x = []
    for a, ang_mom in zip(a_x, ang_mom_x):
        symbol = atoms.symbols[a]
        sym_ang_mom = '.'.join([str(symbol), str(ang_mom)])
        if sym_ang_mom in sym_ang_mom_i:
            i = sym_ang_mom_i.index(sym_ang_mom)
        else:
            i = len(sym_ang_mom_i)
            sym_ang_mom_i.append(sym_ang_mom)
        i_x.append(i)

    nk, nb = len(calc.get_ibz_k_points()), calc.get_number_of_bands()
    projector_weight_skni = np.zeros((ns, nk, nb, len(sym_ang_mom_i)))
    ali_x = [(a, ang_mom, i) for (a, ang_mom, i) in zip(a_x, ang_mom_x, i_x)]
    
    from gpaw.utilities.dos import raw_orbital_LDOS
    
    for _, (a, ang_mom, i) in enumerate(ali_x):
        # Extract weights
        for s in range(ns):
            __, weights = raw_orbital_LDOS(calc, a, s, ang_mom)
            projector_weight_kn = weights.reshape((nk, nb))
            projector_weight_kn /= calc.wfs.kd.weight_k[:, np.newaxis]
            projector_weight_skni[s, :, :, i] += projector_weight_kn

    return projector_weight_skni, sym_ang_mom_i

# End of Projected Band Structure related functions----------------


def release_stage_resources(solver):
    """Release calculator references and synchronize calculation ranks."""
    atoms = getattr(
        solver,
        'bulk_configuration',
        None,
    )

    if atoms is not None:
        atoms.calc = None

    gc.collect()
    world.barrier()


def should_split_gpaw_optical(config):
    """Return whether a mixed GPAW workflow needs an optical process boundary."""
    if (
        config.Engine != 'GPAW'
        or not getattr(config, 'Optical_calc', False)
    ):
        return False

    return any((
        getattr(config, 'Ground_calc', False),
        getattr(config, 'Elastic_calc', False),
        getattr(config, 'DOS_calc', False),
        getattr(config, 'Band_calc', False),
        getattr(config, 'Density_calc', False),
        getattr(config, 'Phonon_calc', False),
    ))


def run_gpaw_stage_processes(
    parallel,
    filtered_args,
):
    """Run GPAW electronic and optical stages in fresh processes."""
    try:
        command = build_gpaw_process_command(
            parallel,
            filtered_args,
        )
    except RuntimeError as exc:
        parprint(
            'ERROR: ' + str(exc)
        )
        return 1

    for stage_group in ('electronic', 'optical'):
        child_env = os.environ.copy()
        child_env[GPAW_STAGE_GROUP_ENV] = stage_group
        child_env['OMP_NUM_THREADS'] = '1'
        child_env['OPENBLAS_NUM_THREADS'] = '1'
        child_env['MKL_NUM_THREADS'] = '1'
        child_env['VECLIB_MAXIMUM_THREADS'] = '1'
        child_env['NUMEXPR_NUM_THREADS'] = '1'
        child_env['OMP_DYNAMIC'] = 'FALSE'

        if parallel is not None:
            child_env.pop(
                'GPAW_MPI_BACKEND',
                None,
            )

        parprint(
            "Starting GPAW "
            + stage_group
            + " stage process..."
        )
        completed = subprocess.run(
            command,
            env=child_env,
            check=False,
        )

        if completed.returncode != 0:
            parprint(
                "ERROR: GPAW "
                + stage_group
                + " stage process failed with exit code "
                + str(completed.returncode)
                + "."
            )
            return completed.returncode

    return 0


def run_calculation_stages(
    solver,
    config,
    stages=None,
):
    """Run all requested DFT stages in dependency-safe order."""
    if stages is None:
        stages = resolve_calculation_stages(config)

    for stage in stages:
        if stage == 'optical':
            release_stage_resources(
                solver
            )

        getattr(
            solver,
            f'{stage}calc',
        )()

        if stage == 'optical':
            release_stage_resources(
                solver
            )


def required_dft_executables(config):
    """Return external executables required by the selected workflow."""
    if config.Engine != 'QE':
        return ()

    executables = set()

    if config.Ground_calc:
        executables.add('pw.x')

    if config.DOS_calc:
        executables.update({
            'pw.x',
            'dos.x',
            'projwfc.x',
        })

    if config.Band_calc:
        executables.update({
            'pw.x',
            'bands.x',
        })

        if config.Projected_band_plot:
            executables.add('projwfc.x')

    if config.Density_calc:
        executables.add('pp.x')

    if config.Phonon_calc:
        executables.update({
            'ph.x',
            'q2r.x',
            'matdyn.x',
        })

    if config.Optical_calc:
        executables.update({
            'pw.x',
            'epsilon.x',
        })

    return tuple(
        sorted(executables)
    )


def check_dft_configuration(
    config,
    struct,
    parallel_cores=1,
    check_executables=True,
    check_saved_state=True,
):
    """Check workflow capabilities and runtime dependencies without running DFT."""
    checks = []

    def add(status, name, detail):
        checks.append({
            'status': status,
            'name': name,
            'detail': detail,
        })

    stages = resolve_calculation_stages(
        config
    )
    add(
        'ok',
        'workflow',
        ' -> '.join(stages),
    )

    try:
        parallel_cores = int(
            parallel_cores
        )
    except (TypeError, ValueError):
        parallel_cores = 0

    if parallel_cores <= 0:
        add(
            'error',
            'parallel',
            'Parallel core count must be a positive integer.',
        )
    else:
        add(
            'ok',
            'parallel',
            f'{parallel_cores} process(es)',
        )

    if config.bulk_configuration is None:
        add(
            'error',
            'structure',
            'No atomic structure was supplied.',
        )
    else:
        add(
            'ok',
            'structure',
            str(
                config.bulk_configuration
                .get_global_number_of_atoms()
            )
            + ' atom(s)',
        )

    if config.Cut_off_energy <= 0:
        add(
            'error',
            'cutoff',
            'Cut_off_energy must be positive.',
        )
    else:
        add(
            'ok',
            'cutoff',
            f'{config.Cut_off_energy:g} eV',
        )

    if config.Energy_max <= config.Energy_min:
        add(
            'error',
            'energy-window',
            'Energy_max must be greater than Energy_min.',
        )

    engine = load_engine_module(
        config.Engine
    )

    if config.Engine == 'QE':
        if config.Mode != 'PW':
            add(
                'error',
                'mode',
                'The QE backend supports PW mode only.',
            )
        else:
            add(
                'ok',
                'mode',
                'PW',
            )

        if str(config.vdW_calc).strip().upper() != 'NONE':
            add(
                'error',
                'vdW',
                'QE vdW corrections are not supported yet.',
            )

        if config.SOC_calc:
            add(
                'error',
                'soc',
                'QE SOC workflows are not supported yet.',
            )

        if config.Elastic_calc:
            add(
                'error',
                'elastic',
                'Native QE elastic calculations are not supported yet.',
            )

        hybrid = str(
            config.XC_calc
        ).strip().lower() in {
            'hse06',
            'hse03',
            'pbe0',
        }

        if hybrid and config.Geo_optim:
            add(
                'error',
                'hybrid-geometry',
                'QE hybrid geometry optimization is not supported.',
            )

        if hybrid and config.Phonon_calc:
            add(
                'error',
                'hybrid-phonon',
                'QE hybrid phonons are not supported.',
            )

        if hybrid and config.Optical_calc:
            add(
                'error',
                'hybrid-optical',
                'QE hybrid optical calculations are not supported.',
            )

        if (
            config.Optical_calc
            and str(config.Opt_calc_type).strip().upper() != 'RPA'
        ):
            add(
                'error',
                'optical-method',
                "Native QE optics requires Opt_calc_type = 'RPA'.",
            )

        if config.Optical_calc:
            if config.Opt_max_en <= config.Opt_min_en:
                add(
                    'error',
                    'optical-grid',
                    'Opt_max_en must be greater than Opt_min_en.',
                )

            if int(config.Opt_num_of_data) < 2:
                add(
                    'error',
                    'optical-grid',
                    'Opt_num_of_data must be at least 2.',
                )

        if config.DOS_calc:
            try:
                dos_occupation = engine.resolve_qe_occupation(
                    config.DOS_occupation
                )
            except Exception as exc:
                add(
                    'error',
                    'dos-occupation',
                    str(exc),
                )
            else:
                if dos_occupation['occupations'] not in {
                    'tetrahedra',
                    'tetrahedra_lin',
                    'tetrahedra_opt',
                }:
                    add(
                        'error',
                        'dos-occupation',
                        'QE DOS requires a tetrahedron occupation.',
                    )

        if check_executables:
            for executable in required_dft_executables(config):
                resolved = shutil.which(
                    executable
                )

                if resolved is None:
                    add(
                        'error',
                        f'executable:{executable}',
                        f'{executable} was not found in PATH.',
                    )
                else:
                    add(
                        'ok',
                        f'executable:{executable}',
                        resolved,
                    )

        if check_executables and parallel_cores > 1:
            mpi_executable = (
                shutil.which('mpiexec')
                or shutil.which('mpirun')
                or shutil.which('srun')
            )

            if mpi_executable is None:
                add(
                    'error',
                    'mpi-launcher',
                    'mpiexec, mpirun, or srun was not found.',
                )
            else:
                add(
                    'ok',
                    'mpi-launcher',
                    mpi_executable,
                )

        pseudo_required = any((
            config.Ground_calc,
            config.DOS_calc,
            config.Band_calc,
            config.Optical_calc,
        ))

        if pseudo_required and config.bulk_configuration is not None:
            try:
                pseudo_dir = get_qe_pseudo_dir(
                    relativistic='scalar',
                )
                pseudopotentials = resolve_qe_pseudopotentials(
                    config.bulk_configuration,
                    relativistic='scalar',
                )
            except Exception as exc:
                add(
                    'error',
                    'pseudopotentials',
                    str(exc),
                )
            else:
                add(
                    'ok',
                    'pseudopotentials',
                    f'{len(pseudopotentials)} species in {pseudo_dir}',
                )

        if check_saved_state and not config.Ground_calc:
            state_dir = Path(
                struct
                + '-GROUND-QE-Result-State'
            )

            if engine.has_qe_state(
                state_dir,
                prefix='nanoworks',
            ):
                add(
                    'ok',
                    'ground-state',
                    str(state_dir),
                )
            else:
                add(
                    'error',
                    'ground-state',
                    'Ground_calc is False and no valid QE state was found at '
                    + str(state_dir),
                )

    elif config.Engine == 'GPAW':
        if config.Mode not in {
            'PW',
            'LCAO',
        }:
            add(
                'error',
                'mode',
                'The GPAW backend supports PW and LCAO modes.',
            )
        else:
            add(
                'ok',
                'mode',
                config.Mode,
            )

        if (
            check_executables
            and importlib.util.find_spec('gpaw') is None
        ):
            add(
                'error',
                'python:gpaw',
                'The gpaw Python package is not installed.',
            )
        elif check_executables:
            add(
                'ok',
                'python:gpaw',
                'installed',
            )

        for stage, module_name in (
            ('elastic', 'elastic'),
            ('phonon', 'phonopy'),
        ):
            if (
                getattr(config, f'{stage.capitalize()}_calc')
                and importlib.util.find_spec(module_name) is None
            ):
                add(
                    'error',
                    f'python:{module_name}',
                    f'{module_name} is required for the {stage} stage.',
                )

        if config.Optical_calc and config.Mode != 'PW':
            add(
                'error',
                'optical-mode',
                'GPAW optical calculations require PW mode.',
            )

        if check_executables and parallel_cores > 1:
            mpi_executable = (
                shutil.which('mpiexec')
                or shutil.which('mpirun')
                or shutil.which('srun')
            )

            if mpi_executable is None:
                add(
                    'error',
                    'mpi-launcher',
                    'mpiexec, mpirun, or srun was not found.',
                )

            if shutil.which('gpaw') is None:
                add(
                    'error',
                    'executable:gpaw',
                    'The gpaw command was not found in PATH.',
                )

        if check_saved_state and not config.Ground_calc:
            state_file = Path(
                struct
                + '-GROUND-GPAW-Result-State.gpw'
            )

            if state_file.is_file():
                add(
                    'ok',
                    'ground-state',
                    str(state_file),
                )
            else:
                add(
                    'error',
                    'ground-state',
                    'Ground_calc is False and no GPAW state was found at '
                    + str(state_file),
                )

    errors = [
        check
        for check in checks
        if check['status'] == 'error'
    ]

    return {
        'ok': not errors,
        'engine': config.Engine,
        'stages': stages,
        'checks': checks,
        'errors': errors,
    }


def format_dft_preflight_report(report):
    """Format a deterministic human-readable preflight report."""
    lines = [
        'dftsolve preflight check',
        f"Engine: {report['engine']}",
        'Stages: ' + ' -> '.join(report['stages']),
        '',
    ]

    for check in report['checks']:
        label = check['status'].upper()
        lines.append(
            f"[{label}] {check['name']}: {check['detail']}"
        )

    lines.extend([
        '',
        (
            'Result: READY'
            if report['ok']
            else f"Result: BLOCKED ({len(report['errors'])} error(s))"
        ),
    ])

    return '\n'.join(lines)


def format_dft_preflight_json(report):
    """Format a versioned machine-readable preflight report."""
    payload = {
        'schema_version': 1,
        'ok': bool(report['ok']),
        'engine': report['engine'],
        'stages': list(report['stages']),
        'checks': report['checks'],
        'errors': report['errors'],
        'error_count': len(report['errors']),
    }

    return json.dumps(
        payload,
        indent=2,
        sort_keys=True,
    )


def prepare_qe_dry_run(
    config,
    struct,
    parallel_cores=1,
):
    """Write QE input files and a command plan without executing QE."""
    if config.Engine != 'QE':
        raise NotImplementedError(
            "dftsolve --dry-run currently supports Engine = 'QE' only."
        )

    engine = load_engine_module('QE')
    validated_xc = engine.validate_qe_xc(
        config.XC_calc,
        pseudo_xc='pbe',
        allow_hybrid=True,
    )
    hybrid = str(validated_xc).strip().lower() in {
        'hse06',
        'hse03',
        'pbe0',
    }

    parallel_cores = int(parallel_cores)

    if parallel_cores <= 0:
        raise ValueError(
            "Dry-run parallel core count must be a positive integer."
        )

    struct = Path(struct).expanduser().resolve()
    struct.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    atoms = config.bulk_configuration
    jobs = []
    setup_directories = set()
    notes = []

    if config.Geo_optim and any((
        config.DOS_calc,
        config.Band_calc,
        config.Density_calc,
        config.Phonon_calc,
        config.Optical_calc,
    )):
        notes.append(
            "Downstream inputs contain the supplied geometry. Regenerate the "
            "dry-run deck after relaxation if the optimized geometry must be "
            "embedded in those inputs."
        )

    def command_for(executable, input_file):
        command = []

        if parallel_cores > 1:
            launcher = (
                shutil.which('mpiexec')
                or shutil.which('mpirun')
                or shutil.which('srun')
                or 'mpiexec'
            )
            flag = (
                '-n'
                if Path(launcher).name == 'srun'
                else '-np'
            )
            command.extend([
                str(launcher),
                flag,
                str(parallel_cores),
            ])

        command.extend([
            executable,
            '-i',
            str(input_file),
        ])
        return command

    def add_job(
        job_id,
        stage,
        executable,
        input_file,
        output_file,
        input_text,
        depends_on=None,
        working_directory=None,
        metadata=None,
    ):
        input_file = Path(input_file).expanduser().resolve()
        output_file = Path(output_file).expanduser().resolve()
        input_file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        input_file.write_text(
            input_text,
            encoding='utf-8',
        )

        command = command_for(
            executable,
            input_file,
        )
        working_directory = (
            Path(working_directory).expanduser().resolve()
            if working_directory is not None
            else None
        )

        if working_directory is not None:
            setup_directories.add(
                str(working_directory)
            )

        setup_directories.add(
            str(output_file.parent)
        )
        jobs.append({
            'id': job_id,
            'stage': stage,
            'executable': executable,
            'input_file': str(input_file),
            'output_file': str(output_file),
            'working_directory': (
                str(working_directory)
                if working_directory is not None
                else None
            ),
            'command': command,
            'depends_on': list(depends_on or []),
            'metadata': dict(metadata or {}),
        })

    ground_state_dir = Path(
        str(struct) + '-GROUND-QE-Result-State'
    ).resolve()
    setup_directories.add(
        str(ground_state_dir)
    )

    pseudo_dir = None
    pseudopotentials = None
    magnetic_moments = None
    needs_pw_input = any((
        config.Ground_calc,
        config.DOS_calc,
        config.Band_calc,
        config.Optical_calc,
    ))

    if needs_pw_input:
        pseudo_dir = get_qe_pseudo_dir(
            relativistic='scalar',
        )
        pseudopotentials = resolve_qe_pseudopotentials(
            atoms,
            relativistic='scalar',
        )

        if config.Spin_calc:
            magnetic_moments = resolve_initial_magnetic_moments(
                atoms=atoms,
                magmom_per_atom=config.Magmom_per_atom,
                magmom_single_atom=config.Magmom_single_atom,
            )

    ground_gamma = (
        config.Gamma
        if config.Ground_gamma is None
        else config.Ground_gamma
    )
    common_pw = {
        'atoms': atoms,
        'pseudopotentials': pseudopotentials,
        'cutoff_ev': config.Cut_off_energy,
        'total_charge': config.Total_charge,
        'spinpol': config.Spin_calc,
        'magnetic_moments': magnetic_moments,
        'setup_params': config.Setup_params,
        'xc_calc': config.XC_calc,
        'exx_fraction': config.XC_exx_fraction,
        'omega': config.XC_omega,
        'prefix': 'nanoworks',
        'pseudo_dir': pseudo_dir,
        'outdir': ground_state_dir,
    }

    if config.Ground_calc:
        ground_mesh = engine.resolve_qe_kpoint_size(
            atoms,
            density=config.Ground_kpts_density,
            size=(
                config.Ground_kpts_x,
                config.Ground_kpts_y,
                config.Ground_kpts_z,
            ),
        )
        ground_occupation = engine.resolve_qe_occupation(
            config.Occupation
        )
        ground_kwargs = {
            **common_pw,
            'kpoint_size': ground_mesh,
            'gamma': ground_gamma,
            'nbands': config.Ground_num_of_bands,
            'occupations': ground_occupation['occupations'],
            'smearing': ground_occupation['smearing'],
            'width_ev': ground_occupation['width_ev'],
        }

        if config.Geo_optim:
            variable_cell = True in config.Relax_cell
            label = 'VC-RELAX' if variable_cell else 'RELAX'
            input_text = engine.render_relax_input(
                optimizer=config.Optimizer,
                max_force=config.Max_F_tolerance,
                max_step=config.Max_step,
                relax_cell=config.Relax_cell,
                hydrostatic_pressure=config.Hydrostatic_pressure,
                fix_symmetry=config.Fix_symmetry,
                **ground_kwargs,
            )
        else:
            label = 'SCF'
            input_text = engine.render_scf_input(
                **ground_kwargs
            )

        add_job(
            'ground',
            'ground',
            'pw.x',
            Path(str(struct) + f'-GROUND-QE-Input-{label}.in'),
            Path(str(struct) + f'-GROUND-QE-Log-{label}.txt'),
            input_text,
        )
    else:
        notes.append(
            "Ground_calc is False; generated post-processing commands expect "
            f"an existing QE state in {ground_state_dir}."
        )

    ground_dependency = ['ground'] if config.Ground_calc else []

    if config.DOS_calc:
        (
            dos_density,
            dos_size,
            dos_gamma,
        ) = resolve_stage_kpoint_settings(
            stage_density=config.DOS_kpts_density,
            stage_size=(
                config.DOS_kpts_x,
                config.DOS_kpts_y,
                config.DOS_kpts_z,
            ),
            stage_gamma=config.DOS_gamma,
            ground_density=config.Ground_kpts_density,
            ground_size=(
                config.Ground_kpts_x,
                config.Ground_kpts_y,
                config.Ground_kpts_z,
            ),
            ground_gamma=ground_gamma,
        )
        dos_mesh = engine.resolve_qe_kpoint_size(
            atoms,
            density=dos_density,
            size=dos_size,
        )
        dos_occupation = resolve_stage_occupation(
            config.DOS_occupation,
            config.Occupation,
        )
        qe_dos_occupation = engine.resolve_qe_occupation(
            dos_occupation
        )
        delta_e = (
            float(config.Energy_max)
            - float(config.Energy_min)
        ) / (int(config.DOS_npoints) - 1)

        if hybrid:
            dos_state_dir = Path(
                str(struct) + '-DOS-QE-Result-State'
            ).resolve()
            setup_directories.add(
                str(dos_state_dir)
            )
            dos_pw = {
                **common_pw,
                'outdir': dos_state_dir,
            }
            dos_electronic_job = 'dos-hybrid-scf'
            dos_input_file = Path(
                str(struct) + '-DOS-QE-Input-Hybrid-SCF.in'
            )
            dos_output_file = Path(
                str(struct) + '-DOS-QE-Log-Hybrid-SCF.txt'
            )
            dos_input_text = engine.render_scf_input(
                **dos_pw,
                kpoint_size=dos_mesh,
                gamma=dos_gamma,
                nbands=config.DOS_num_of_bands,
                occupations=qe_dos_occupation['occupations'],
                smearing=qe_dos_occupation['smearing'],
                width_ev=qe_dos_occupation['width_ev'],
            )
        else:
            dos_state_dir = ground_state_dir
            dos_electronic_job = 'dos-nscf'
            dos_input_file = Path(
                str(struct) + '-DOS-QE-Input-NSCF.in'
            )
            dos_output_file = Path(
                str(struct) + '-DOS-QE-Log-NSCF.txt'
            )
            dos_input_text = engine.render_nscf_input(
                **common_pw,
                kpoint_size=dos_mesh,
                gamma=dos_gamma,
                nbands=config.DOS_num_of_bands,
                occupations=qe_dos_occupation['occupations'],
                smearing=qe_dos_occupation['smearing'],
                width_ev=qe_dos_occupation['width_ev'],
            )

        add_job(
            dos_electronic_job,
            'dos',
            'pw.x',
            dos_input_file,
            dos_output_file,
            dos_input_text,
            depends_on=ground_dependency,
            metadata={
                'calculation': (
                    'scf'
                    if hybrid
                    else 'nscf'
                ),
                'hybrid': hybrid,
                'kpoint_size': list(dos_mesh),
            },
        )
        add_job(
            'dos-total',
            'dos',
            'dos.x',
            Path(str(struct) + '-DOS-QE-Input-DOS.in'),
            Path(str(struct) + '-DOS-QE-Log-DOS.txt'),
            engine.render_dos_input(
                prefix='nanoworks',
                outdir=dos_state_dir,
                fildos=Path(str(struct) + '-DOS-QE-Result-Raw-DOS.dat'),
                bz_sum=qe_dos_occupation['occupations'],
                delta_e=delta_e,
            ),
            depends_on=[dos_electronic_job],
        )
        add_job(
            'dos-projected',
            'dos',
            'projwfc.x',
            Path(str(struct) + '-DOS-QE-Input-PDOS.in'),
            Path(str(struct) + '-DOS-QE-Log-PDOS.txt'),
            engine.render_projwfc_input(
                prefix='nanoworks',
                outdir=dos_state_dir,
                filpdos=Path(str(struct) + '-DOS-QE-Result-Raw-PDOS'),
                delta_e=delta_e,
            ),
            depends_on=[dos_electronic_job],
        )
        notes.append(
            "DOS and PDOS dry-run inputs omit Emin/Emax because the absolute "
            "energy window depends on the electronic-stage Fermi energy."
        )

    if config.Band_calc:
        band_path = engine.build_band_path(
            atoms,
            path=config.Band_path,
            npoints=config.Band_npoints,
        )
        band_occupation = engine.resolve_qe_occupation(
            config.Occupation
        )

        if hybrid:
            band_state_dir = Path(
                str(struct) + '-BAND-QE-Result-State'
            ).resolve()
            setup_directories.add(
                str(band_state_dir)
            )
            band_mesh = engine.resolve_qe_kpoint_size(
                atoms,
                density=config.Ground_kpts_density,
                size=(
                    config.Ground_kpts_x,
                    config.Ground_kpts_y,
                    config.Ground_kpts_z,
                ),
            )
            additional_kpoints = (
                engine.build_qe_exx_additional_kpoints(
                    band_path=band_path,
                    qpoint_grid=band_mesh,
                )
            )
            helper_indices = list(
                additional_kpoints['helper_indices']
            )
            hybrid_index_metadata = {
                'band_indices': list(
                    additional_kpoints['band_indices']
                ),
                'helper_count': len(helper_indices),
                'helper_index_range': (
                    [
                        helper_indices[0],
                        helper_indices[-1] + 1,
                    ]
                    if helper_indices
                    else []
                ),
            }
            band_pw = {
                **common_pw,
                'outdir': band_state_dir,
            }
            add_job(
                'band-hybrid-scf',
                'band',
                'pw.x',
                Path(
                    str(struct)
                    + '-BAND-QE-Input-Hybrid-SCF.in'
                ),
                Path(
                    str(struct)
                    + '-BAND-QE-Log-Hybrid-SCF.txt'
                ),
                engine.render_scf_input(
                    **band_pw,
                    kpoint_size=band_mesh,
                    gamma=ground_gamma,
                    nbands=config.Band_num_of_bands,
                    occupations=band_occupation['occupations'],
                    smearing=band_occupation['smearing'],
                    width_ev=band_occupation['width_ev'],
                    exx_additional_kpoints=additional_kpoints,
                ),
                depends_on=ground_dependency,
                metadata={
                    'calculation': 'scf',
                    'hybrid': True,
                    'qpoint_grid': list(
                        additional_kpoints['qpoint_grid']
                    ),
                    **hybrid_index_metadata,
                },
            )
            add_job(
                'band-postprocess',
                'band',
                'bands.x',
                Path(
                    str(struct)
                    + '-BAND-QE-Input-Bands.x.in'
                ),
                Path(
                    str(struct)
                    + '-BAND-QE-Log-Bands.x.txt'
                ),
                engine.render_bands_postprocess_input(
                    prefix='nanoworks',
                    outdir=band_state_dir,
                    filband=Path(
                        str(struct)
                        + '-BAND-QE-Result-Bands.x.dat'
                    ),
                    lsym=False,
                ),
                depends_on=['band-hybrid-scf'],
                metadata={
                    **hybrid_index_metadata,
                },
            )
            band_projection_dependency = 'band-postprocess'
        else:
            band_state_dir = ground_state_dir
            add_job(
                'band',
                'band',
                'pw.x',
                Path(str(struct) + '-BAND-QE-Input-Bands.in'),
                Path(str(struct) + '-BAND-QE-Log-Bands.txt'),
                engine.render_bands_input(
                    **common_pw,
                    band_path=band_path,
                    nbands=config.Band_num_of_bands,
                    occupations=band_occupation['occupations'],
                    smearing=band_occupation['smearing'],
                    width_ev=band_occupation['width_ev'],
                ),
                depends_on=ground_dependency,
                metadata={
                    'calculation': 'bands',
                    'hybrid': False,
                },
            )
            band_projection_dependency = 'band'

        if config.Projected_band_plot:
            add_job(
                'band-projections',
                'band',
                'projwfc.x',
                Path(str(struct) + '-BAND-QE-Input-Projections.in'),
                Path(str(struct) + '-BAND-QE-Log-Projections.txt'),
                engine.render_projwfc_input(
                    prefix='nanoworks',
                    outdir=band_state_dir,
                    filpdos=Path(
                        str(struct)
                        + '-BAND-QE-Result-Projections-pdos'
                    ),
                    filproj=Path(
                        str(struct)
                        + '-BAND-QE-Result-Projections'
                    ),
                    lsym=False,
                    diag_basis=False,
                ),
                depends_on=[band_projection_dependency],
                metadata={
                    'hybrid': hybrid,
                    **(
                        hybrid_index_metadata
                        if hybrid
                        else {}
                    ),
                },
            )

    if config.Density_calc:
        density_jobs = [
            ('Pseudo-Total', 0, 0 if config.Spin_calc else None),
        ]

        if config.Spin_calc:
            density_jobs.extend([
                ('Pseudo-Up', 0, 1),
                ('Pseudo-Down', 0, 2),
                ('Spin-Density', 6, None),
            ])

        for label, plot_num, spin_component in density_jobs:
            add_job(
                'density-' + label.lower(),
                'density',
                'pp.x',
                Path(str(struct) + f'-EDENSITY-QE-Input-{label}.in'),
                Path(str(struct) + f'-EDENSITY-QE-Log-{label}.txt'),
                engine.render_pp_input(
                    prefix='nanoworks',
                    outdir=ground_state_dir,
                    filplot=Path(
                        str(struct)
                        + f'-EDENSITY-QE-Result-{label}.dat'
                    ),
                    fileout=Path(
                        str(struct)
                        + f'-EDENSITY-QE-Result-{label}.cube'
                    ),
                    plot_num=plot_num,
                    spin_component=spin_component,
                ),
                depends_on=ground_dependency,
            )

    if config.Phonon_calc:
        qpoint_grid = engine.resolve_qe_phonon_qpoint_grid(
            config.Phonon_supercell
        )
        phonon_path = engine.build_band_path(
            atoms,
            path=config.Phonon_path,
            npoints=config.Phonon_npoints,
        )
        fildyn = Path(
            str(struct) + '-PHONON-QE-Result-Dynamical-Matrix'
        )
        flfrc = Path(
            str(struct) + '-PHONON-QE-Result-Force-Constants.fc'
        )
        add_job(
            'phonon-grid',
            'phonon',
            'ph.x',
            Path(str(struct) + '-PHONON-QE-Input-PH.in'),
            Path(str(struct) + '-PHONON-QE-Log-PH.txt'),
            engine.render_ph_input(
                prefix='nanoworks',
                outdir=ground_state_dir,
                fildyn=fildyn,
                qpoint_grid=qpoint_grid,
            ),
            depends_on=ground_dependency,
        )
        add_job(
            'phonon-force-constants',
            'phonon',
            'q2r.x',
            Path(str(struct) + '-PHONON-QE-Input-Q2R.in'),
            Path(str(struct) + '-PHONON-QE-Log-Q2R.txt'),
            engine.render_q2r_input(
                fildyn=fildyn,
                flfrc=flfrc,
                zasr='no',
            ),
            depends_on=['phonon-grid'],
        )
        add_job(
            'phonon-band',
            'phonon',
            'matdyn.x',
            Path(str(struct) + '-PHONON-QE-Input-Matdyn-Band.in'),
            Path(str(struct) + '-PHONON-QE-Log-Matdyn-Band.txt'),
            engine.render_matdyn_band_input(
                flfrc=flfrc,
                flfrq=Path(
                    str(struct) + '-PHONON-QE-Result-Band.freq'
                ),
                band_path=phonon_path,
                acoustic_sum_rule=config.Phonon_acoustic_sum_rule,
            ),
            depends_on=['phonon-force-constants'],
        )
        add_job(
            'phonon-dos',
            'phonon',
            'matdyn.x',
            Path(str(struct) + '-PHONON-QE-Input-Matdyn-DOS.in'),
            Path(str(struct) + '-PHONON-QE-Log-Matdyn-DOS.txt'),
            engine.render_matdyn_dos_input(
                flfrc=flfrc,
                fldos=Path(
                    str(struct) + '-PHONON-QE-Result-DOS.dat'
                ),
                qpoint_grid=(
                    config.Phonon_qpts_x,
                    config.Phonon_qpts_y,
                    config.Phonon_qpts_z,
                ),
                acoustic_sum_rule=config.Phonon_acoustic_sum_rule,
            ),
            depends_on=['phonon-force-constants'],
        )

    if config.Optical_calc:
        (
            optical_density,
            optical_size,
            optical_gamma,
        ) = resolve_stage_kpoint_settings(
            stage_density=config.Opt_kpts_density,
            stage_size=(
                config.Opt_kpts_x,
                config.Opt_kpts_y,
                config.Opt_kpts_z,
            ),
            stage_gamma=config.Opt_gamma,
            ground_density=config.Ground_kpts_density,
            ground_size=(
                config.Ground_kpts_x,
                config.Ground_kpts_y,
                config.Ground_kpts_z,
            ),
            ground_gamma=ground_gamma,
        )
        optical_mesh = engine.resolve_qe_kpoint_size(
            atoms,
            density=optical_density,
            size=optical_size,
        )
        optical_occupation = engine.resolve_qe_occupation({
            'name': 'fermi-dirac',
            'width': config.Opt_FD_smearing,
        })
        add_job(
            'optical-nscf',
            'optical',
            'pw.x',
            Path(str(struct) + '-OPTICAL-QE-Input-NSCF.in'),
            Path(str(struct) + '-OPTICAL-QE-Log-NSCF.txt'),
            engine.render_nscf_input(
                **common_pw,
                kpoint_size=optical_mesh,
                gamma=optical_gamma,
                nbands=config.Opt_num_of_bands,
                occupations=optical_occupation['occupations'],
                smearing=optical_occupation['smearing'],
                width_ev=optical_occupation['width_ev'],
                nosym=True,
            ),
            depends_on=ground_dependency,
        )
        optical_result_dir = Path(
            str(struct) + '-OPTICAL-QE-Result-Raw'
        ).resolve()
        add_job(
            'optical-epsilon',
            'optical',
            'epsilon.x',
            Path(str(struct) + '-OPTICAL-QE-Input-Epsilon.in'),
            Path(str(struct) + '-OPTICAL-QE-Log-Epsilon.txt'),
            engine.render_epsilon_input(
                prefix='nanoworks',
                outdir=ground_state_dir,
                calculation='eps',
                smeartype='gauss',
                intersmear=config.Opt_eta,
                intrasmear=0.0,
                wmin=config.Opt_min_en,
                wmax=config.Opt_max_en,
                nw=config.Opt_num_of_data,
                shift=config.Opt_shift_en,
            ),
            depends_on=['optical-nscf'],
            working_directory=optical_result_dir,
        )

    plan_file = Path(
        str(struct) + '-DRYRUN-QE-Plan.json'
    )
    script_file = Path(
        str(struct) + '-DRYRUN-QE-Run.sh'
    )
    plan = {
        'schema_version': 1,
        'engine': 'QE',
        'dry_run': True,
        'scheduler': 'local',
        'parallel_cores': parallel_cores,
        'spin_polarized': bool(config.Spin_calc),
        'stages': list(resolve_calculation_stages(config)),
        'jobs': jobs,
        'notes': notes,
        'setup_directories': sorted(setup_directories),
        'plan_file': str(plan_file),
        'script_file': str(script_file),
    }
    plan_file.write_text(
        json.dumps(
            plan,
            indent=2,
            sort_keys=True,
        ) + '\n',
        encoding='utf-8',
    )

    script_lines = [
        '#!/usr/bin/env bash',
        'set -euo pipefail',
        '',
        '# Generated by dftsolve --dry-run; inspect before execution.',
    ]

    for directory in sorted(setup_directories):
        script_lines.append(
            'mkdir -p ' + shlex.quote(directory)
        )

    script_lines.append('')

    for job in jobs:
        command = shlex.join(job['command'])
        output_file = shlex.quote(job['output_file'])

        if job['working_directory'] is None:
            script_lines.append(
                f"{command} > {output_file}"
            )
        else:
            script_lines.append(
                '(cd '
                + shlex.quote(job['working_directory'])
                + ' && '
                + command
                + ' > '
                + output_file
                + ')'
            )

    script_file.write_text(
        '\n'.join(script_lines) + '\n',
        encoding='utf-8',
    )
    script_file.chmod(
        script_file.stat().st_mode | 0o111
    )

    return plan


SLURM_PROFILE_KEYS = {
    'time',
    'memory',
    'partition',
    'account',
    'qos',
    'modules',
    'job_name',
}


def resolve_slurm_profile_path(profile):
    """Resolve a Slurm profile path or a named user profile."""
    profile = str(profile).strip()

    if not profile:
        raise ValueError(
            "Slurm profile name must not be empty."
        )

    requested = Path(profile).expanduser()

    if requested.is_file():
        return requested.resolve()

    if requested.parent != Path('.'):
        raise FileNotFoundError(
            f"Slurm profile was not found: {requested}"
        )

    profile_name = requested.name

    if not profile_name.endswith('.json'):
        profile_name += '.json'

    named_profile = (
        Path.home()
        / '.config'
        / 'nanoworks'
        / 'clusters'
        / profile_name
    )

    if not named_profile.is_file():
        raise FileNotFoundError(
            "Slurm profile was not found as a file or named profile: "
            + profile
        )

    return named_profile.resolve()


def load_slurm_profile(profile):
    """Load and validate a reusable JSON Slurm profile."""
    profile_file = resolve_slurm_profile_path(
        profile
    )

    try:
        document = json.loads(
            profile_file.read_text(
                encoding='utf-8'
            )
        )
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid JSON in Slurm profile '{profile_file}': {exc}"
        ) from exc

    if not isinstance(document, dict):
        raise ValueError(
            "Slurm profile must contain a JSON object."
        )

    settings = document.get(
        'slurm',
        document,
    )

    if not isinstance(settings, dict):
        raise ValueError(
            "The Slurm profile 'slurm' value must be an object."
        )

    unknown = sorted(
        set(settings) - SLURM_PROFILE_KEYS
    )

    if unknown:
        raise ValueError(
            "Unsupported Slurm profile setting(s): "
            + ', '.join(unknown)
        )

    resolved = dict(settings)
    modules = resolved.get('modules')

    if modules is not None:
        if (
            not isinstance(modules, list)
            or not all(
                isinstance(module, str)
                for module in modules
            )
        ):
            raise ValueError(
                "Slurm profile modules must be a list of strings."
            )

        resolved['modules'] = list(modules)

    for key in SLURM_PROFILE_KEYS - {'modules'}:
        value = resolved.get(key)

        if value is not None and not isinstance(value, str):
            raise ValueError(
                f"Slurm profile {key} must be a string."
            )

    resolved['profile_file'] = str(
        profile_file
    )
    return resolved


def write_qe_slurm_script(
    plan,
    wall_time='24:00:00',
    memory=None,
    partition=None,
    account=None,
    qos=None,
    modules=None,
    job_name=None,
    profile_file=None,
):
    """Write a Slurm batch script from a QE dry-run plan."""
    if plan.get('engine') != 'QE' or not plan.get('dry_run'):
        raise ValueError(
            "A native QE dry-run plan is required for Slurm generation."
        )

    jobs = list(plan.get('jobs', []))

    if not jobs:
        raise ValueError(
            "The QE dry-run plan does not contain any executable jobs."
        )

    parallel_cores = int(plan.get('parallel_cores', 0))

    if parallel_cores <= 0:
        raise ValueError(
            "The QE dry-run plan has an invalid parallel core count."
        )

    def validate_token(name, value, extra_characters=''):
        if value is None:
            return None

        value = str(value).strip()

        if not value:
            raise ValueError(
                f"Slurm {name} must not be empty."
            )

        allowed = set(
            'abcdefghijklmnopqrstuvwxyz'
            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            '0123456789._-'
            + extra_characters
        )

        if any(character not in allowed for character in value):
            raise ValueError(
                f"Slurm {name} contains unsupported characters."
            )

        return value

    wall_time = str(wall_time).strip()
    day_count = None
    clock = wall_time

    if '-' in wall_time:
        day_text, clock = wall_time.split('-', 1)

        if not day_text.isdigit():
            raise ValueError(
                "Slurm time must use HH:MM:SS or D-HH:MM:SS."
            )

        day_count = int(day_text)

    clock_parts = clock.split(':')

    if (
        len(clock_parts) != 3
        or not all(part.isdigit() for part in clock_parts)
    ):
        raise ValueError(
            "Slurm time must use HH:MM:SS or D-HH:MM:SS."
        )

    hours, minutes, seconds = (
        int(part)
        for part in clock_parts
    )

    if (
        minutes >= 60
        or seconds >= 60
        or (day_count is not None and hours >= 24)
        or (
            (day_count is None or day_count == 0)
            and hours == 0
            and minutes == 0
            and seconds == 0
        )
    ):
        raise ValueError(
            "Slurm time contains an invalid clock value."
        )

    memory = validate_token('memory', memory)
    partition = validate_token('partition', partition)
    account = validate_token('account', account)
    qos = validate_token('qos', qos)
    modules = [
        validate_token(
            'module name',
            module,
            extra_characters='+/@:',
        )
        for module in (modules or [])
    ]

    if job_name is None:
        plan_name = Path(
            plan['plan_file']
        ).name
        job_name = plan_name.split(
            '-DRYRUN-',
            1,
        )[0]

    raw_job_name = str(job_name).strip()
    job_name = ''.join(
        character
        if character.isalnum() or character in '-_'
        else '-'
        for character in raw_job_name
    ).strip('-_')[:64]

    if not job_name:
        job_name = 'nanoworks-qe'

    local_script = Path(
        plan['script_file']
    )
    slurm_script = local_script.with_suffix(
        '.slurm'
    )
    output_prefix = str(
        local_script.with_suffix('')
    )
    slurm_output = output_prefix + '-SLURM-%j.out'
    slurm_error = output_prefix + '-SLURM-%j.err'
    lines = [
        '#!/usr/bin/env bash',
        f'#SBATCH --job-name={job_name}',
        '#SBATCH --nodes=1',
        f'#SBATCH --ntasks={parallel_cores}',
        '#SBATCH --cpus-per-task=1',
        f'#SBATCH --time={wall_time}',
        f'#SBATCH --output={slurm_output}',
        f'#SBATCH --error={slurm_error}',
    ]

    if partition is not None:
        lines.append(
            f'#SBATCH --partition={partition}'
        )

    if account is not None:
        lines.append(
            f'#SBATCH --account={account}'
        )

    if qos is not None:
        lines.append(
            f'#SBATCH --qos={qos}'
        )

    if memory is not None:
        lines.append(
            f'#SBATCH --mem={memory}'
        )

    lines.extend([
        '',
        'set -euo pipefail',
        '',
    ])

    if modules:
        for module in modules:
            lines.append(
                'module load ' + shlex.quote(module)
            )
        lines.append('')
    else:
        lines.extend([
            '# Load the site-specific Quantum ESPRESSO module here if needed.',
        ])

    lines.extend([
        'export OMP_NUM_THREADS=1',
        'export OPENBLAS_NUM_THREADS=1',
        'export MKL_NUM_THREADS=1',
        'export VECLIB_MAXIMUM_THREADS=1',
        'export NUMEXPR_NUM_THREADS=1',
        'export OMP_DYNAMIC=FALSE',
        '',
    ])

    for directory in plan.get('setup_directories', []):
        lines.append(
            'mkdir -p ' + shlex.quote(str(directory))
        )

    if plan.get('setup_directories'):
        lines.append('')

    for job in jobs:
        command = shlex.join([
            'srun',
            '-n',
            str(parallel_cores),
            str(job['executable']),
            '-i',
            str(job['input_file']),
        ])
        output_file = shlex.quote(
            str(job['output_file'])
        )
        working_directory = job.get(
            'working_directory'
        )

        if working_directory is None:
            lines.append(
                f"{command} > {output_file}"
            )
        else:
            lines.append(
                '(cd '
                + shlex.quote(str(working_directory))
                + ' && '
                + command
                + ' > '
                + output_file
                + ')'
            )

    slurm_script.write_text(
        '\n'.join(lines) + '\n',
        encoding='utf-8',
    )
    slurm_script.chmod(
        slurm_script.stat().st_mode | 0o111
    )
    plan['scheduler'] = 'slurm'
    plan['slurm'] = {
        'script_file': str(slurm_script),
        'job_name': job_name,
        'nodes': 1,
        'ntasks': parallel_cores,
        'cpus_per_task': 1,
        'time': wall_time,
        'memory': memory,
        'partition': partition,
        'account': account,
        'qos': qos,
        'modules': modules,
        'profile_file': profile_file,
        'output': slurm_output,
        'error': slurm_error,
    }
    plan['slurm_script_file'] = str(
        slurm_script
    )
    Path(plan['plan_file']).write_text(
        json.dumps(
            plan,
            indent=2,
            sort_keys=True,
        ) + '\n',
        encoding='utf-8',
    )

    return slurm_script


def format_dft_dry_run_report(plan):
    """Format a concise summary of generated QE dry-run artifacts."""
    lines = [
        'dftsolve QE dry run',
        'Stages: ' + ' -> '.join(plan['stages']),
        f"Jobs: {len(plan['jobs'])}",
        f"Plan: {plan['plan_file']}",
        f"Script: {plan['script_file']}",
    ]

    if plan.get('slurm_script_file'):
        lines.append(
            f"Slurm: {plan['slurm_script_file']}"
        )

    if plan.get('slurm', {}).get('profile_file'):
        lines.append(
            "Slurm profile: "
            + plan['slurm']['profile_file']
        )

    lines.append(
        'Result: PREPARED (no calculations executed)'
    )
    return '\n'.join(lines)


def main():
    meter = None
    parser = ArgumentParser(prog ='dftsolve.py', description=Description, formatter_class=RawFormatter)
    parser.add_argument("-i", "--input", dest = "inputfile", help="Use input file for calculation variables (also you can insert geometry)")
    parser.add_argument("-g", "--geometry",dest ="geometryfile", help="Use CIF file for geometry")
    parser.add_argument("-a", "--auto", dest="auto", action='store_true', help="Automatically generate input parameters based on geometry")
    parser.add_argument("-v", "--version", dest="version", action='store_true')
    parser.add_argument("-e", "--energy", dest="energymeas", action='store_true')
    parser.add_argument("-p", "--parallel", dest="parallel", type=int, help="Number of cores to run in parallel")
    parser.add_argument("--check", dest="check", action='store_true', help="Validate the selected workflow without running calculations")
    parser.add_argument("--json", dest="json", action='store_true', help="Print --check results as machine-readable JSON")
    parser.add_argument("--dry-run", dest="dry_run", action='store_true', help="Write QE inputs and an execution plan without running calculations")
    parser.add_argument("--scheduler", choices=('local', 'slurm'), help="Execution script type generated by --dry-run")
    parser.add_argument("--cluster-profile", help="Slurm profile JSON file or named profile from ~/.config/nanoworks/clusters")
    parser.add_argument("--slurm-time", help="Slurm wall time (HH:MM:SS or D-HH:MM:SS)")
    parser.add_argument("--slurm-memory", help="Optional Slurm memory request, for example 32G")
    parser.add_argument("--slurm-partition", help="Optional Slurm partition name")
    parser.add_argument("--slurm-account", help="Optional Slurm account/project name")
    parser.add_argument("--slurm-qos", help="Optional Slurm quality-of-service name")
    parser.add_argument("--slurm-module", action='append', help="Module to load in the Slurm script; repeat for multiple modules")
    parser.add_argument("--slurm-job-name", help="Optional Slurm job name")

    args = None

    # Parse arguments
    try:
        if world.rank == 0:
            args = parser.parse_args()
    finally:
        args = broadcast(args, root=0, comm=world)

    if args is None:
        parprint("No arguments used.")
        sys.exit(1)

    if args.json and not args.check:
        parprint("ERROR: --json requires --check.")
        return 2

    if args.check and args.dry_run:
        parprint("ERROR: --check and --dry-run cannot be used together.")
        return 2

    slurm_profile = {}

    if args.cluster_profile is not None:
        try:
            slurm_profile = load_slurm_profile(
                args.cluster_profile
            )
        except (OSError, ValueError) as exc:
            parprint(
                "ERROR: Could not load Slurm profile: "
                + str(exc)
            )
            return 2

    scheduler = args.scheduler

    if scheduler is None:
        scheduler = (
            'slurm'
            if slurm_profile
            else 'local'
        )

    if slurm_profile and scheduler != 'slurm':
        parprint(
            "ERROR: --cluster-profile requires the Slurm scheduler."
        )
        return 2

    if scheduler != 'local' and not args.dry_run:
        parprint("ERROR: --scheduler requires --dry-run.")
        return 2

    energymeas = False
    inFile = None
    configpath = None


    try:
        if args.version == True:
            parprint(f"nanoworks: dftsolve: version: {nanoworks.__version__}")
            sys.exit(1)

        if args.auto:
            if args.geometryfile is None:
                parprint("\033[91mERROR:\033[0m Auto mode (-a) requires a geometry file (-g).")
                sys.exit(1)
            inFile = os.path.join(os.getcwd(), args.geometryfile)
        elif args.inputfile is not None:
            configpath = os.path.join(os.getcwd(),args.inputfile)
            sys.path.append(os.getcwd())
        else:
            parprint("\033[91mERROR:\033[0m Please use an input file with -i argument or use auto mode -a.")
            sys.exit(1)


        if args.geometryfile :
            inFile = os.path.join(os.getcwd(),args.geometryfile)

        if (
            args.energymeas == True
            and not args.check
            and not args.dry_run
        ):
            try:
                import pyRAPL
                energymeas = True
                # Start energy consumption calculation.
                pyRAPL.setup()
                meter = pyRAPL.Measurement('dftsolve')
                meter.begin()
            except:
                parprint("\033[91mERROR:\033[0m Unexpected error while using -e argument.")
                parprint("-e works only with Intel CPUs after Sandy Bridge generation. Do not use with AMD CPUs")
                parprint("You also need to install pymongo and pandas libraries.")
                parprint("If you got permission error, try: sudo chmod -R a+r /sys/class/powercap/intel-rapl")
                parprint("More information about the error:")
                parprint(sys.exc_info()[0])
                sys.exit(1)

    except getopt.error as err:
        # output error, and return with an error code
        parprint (str(err))

    # Start time
    time0 = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))

    # Load struct and config
    if args.auto:
        struct, config = struct_from_auto(
            inFile,
            write_output=(
                not args.check
                and not args.dry_run
            ),
            report_structure=(
                not args.check
                and not args.dry_run
            ),
        )
    else:
        struct, config = struct_from_file(
            inputfile=configpath,
            geometryfile=inFile,
            create_output=(
                not args.check
                and not args.dry_run
            ),
            report_structure=(
                not args.check
                and not args.dry_run
            ),
        )

    if REQUESTED_PARALLEL is not None:
        parallel_cores = REQUESTED_PARALLEL
    else:
        parallel_cores = world.size

    if args.check:
        report = check_dft_configuration(
            config,
            struct=struct,
            parallel_cores=parallel_cores,
        )
        if args.json:
            output = format_dft_preflight_json(
                report
            )
        else:
            output = format_dft_preflight_report(
                report
            )

        parprint(output)
        return 0 if report['ok'] else 2

    if args.dry_run:
        report = check_dft_configuration(
            config,
            struct=struct,
            parallel_cores=parallel_cores,
            check_executables=False,
            check_saved_state=False,
        )

        if not report['ok']:
            parprint(
                format_dft_preflight_report(
                    report
                )
            )
            return 2

        try:
            plan = prepare_qe_dry_run(
                config,
                struct=struct,
                parallel_cores=parallel_cores,
            )

            if scheduler == 'slurm':
                write_qe_slurm_script(
                    plan,
                    wall_time=(
                        args.slurm_time
                        if args.slurm_time is not None
                        else slurm_profile.get(
                            'time',
                            '24:00:00',
                        )
                    ),
                    memory=(
                        args.slurm_memory
                        if args.slurm_memory is not None
                        else slurm_profile.get('memory')
                    ),
                    partition=(
                        args.slurm_partition
                        if args.slurm_partition is not None
                        else slurm_profile.get('partition')
                    ),
                    account=(
                        args.slurm_account
                        if args.slurm_account is not None
                        else slurm_profile.get('account')
                    ),
                    qos=(
                        args.slurm_qos
                        if args.slurm_qos is not None
                        else slurm_profile.get('qos')
                    ),
                    modules=(
                        args.slurm_module
                        if args.slurm_module is not None
                        else slurm_profile.get('modules')
                    ),
                    job_name=(
                        args.slurm_job_name
                        if args.slurm_job_name is not None
                        else slurm_profile.get('job_name')
                    ),
                    profile_file=slurm_profile.get(
                        'profile_file'
                    ),
                )
        except Exception as exc:
            parprint(
                "ERROR: QE dry-run preparation failed: "
                + str(exc)
            )
            return 2

        parprint(
            format_dft_dry_run_report(
                plan
            )
        )
        return 0

    gpaw_stage_group = os.environ.get(
        GPAW_STAGE_GROUP_ENV
    )

    if gpaw_stage_group not in (None, 'electronic', 'optical'):
        parprint(
            "ERROR: Invalid internal GPAW stage group: "
            + gpaw_stage_group
        )
        return 2

    if (
        gpaw_stage_group is None
        and should_split_gpaw_optical(config)
    ):
        if world.size != 1:
            parprint(
                "ERROR: Automatic GPAW optical process separation must "
                "start from the serial dftsolve command. Use "
                "'dftsolve -p N' instead of launching dftsolve under "
                "MPI directly."
            )
            return 2

        parprint(
            "GPAW optical calculation will run in a fresh process "
            "after the electronic stages to release memory."
        )
        return run_gpaw_stage_processes(
            REQUESTED_PARALLEL,
            FILTERED_ARGS,
        )

    # Parallel execution is backend-specific.
    #
    # GPAW requires the complete Python calculation to run under MPI.
    # QE keeps Nanoworks serial and launches pw.x with the requested
    # number of MPI processes.
    if (
        config.Engine == 'GPAW'
        and REQUESTED_PARALLEL is not None
    ):
        restart_gpaw_with_mpi(
            REQUESTED_PARALLEL,
            FILTERED_ARGS,
        )
    
    # Write timings of calculation
    with paropen(
        struct
        + f'-TIMINGS-{config.Engine}-Log-Timings.txt',
        'a',
    ) as f1:
        print("dftsolve.py execution timings (seconds):", end="\n", file=f1)
        print("Execution started:", time0, end="\n", file=f1)

    # Load dftsolve() class with config and core numbers
    dftsolver = dftsolve(
        struct,
        config,
        parallel_cores=parallel_cores,
    )

    # Run structure calculation
    if gpaw_stage_group != 'optical':
        dftsolver.structurecalc()

    # Run ground state and every requested downstream calculation.
    # Optical stays last because its response step can require much more
    # memory than the other post-processing stages.
    selected_stages = None

    if gpaw_stage_group == 'electronic':
        selected_stages = tuple(
            stage
            for stage in resolve_calculation_stages(config)
            if stage != 'optical'
        )
    elif gpaw_stage_group == 'optical':
        selected_stages = ('optical',)

    run_calculation_stages(
        dftsolver,
        config,
        stages=selected_stages,
    )

    # Ending of timings
    with paropen(
        struct
        + f'-TIMINGS-{config.Engine}-Log-Timings.txt',
        'a',
    ) as f1:
        print("---------------------------------------", end="\n", file=f1)

    if args.energymeas == True and meter is not None:
        log_energy_consumption(
            meter,
            struct,
            config.Engine,
        )

    return 0

if __name__ == "__main__":
    sys.exit(main())
