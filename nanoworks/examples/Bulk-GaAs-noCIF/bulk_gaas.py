# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase import Atoms,Atom

Outdirname = 'bulk-gaas-results'

bulk_configuration = Atoms(
    [
    Atom('Ga', ( 0.0, 0.0, 0.0 )),
    Atom('As', ( 2.033, 1.174, 0.830 ))
    ],
    cell=[(4.066, 0.0, 0.0),
          (2.033, 3.521, 0.0),
          (2.033, 1.174, 3.320)],
    pbc=True,
    )

# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Ground_calc = True     # Ground state calculations
Geo_optim = True       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = True         # DOS calculation
Band_calc = True        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# GEOMETRY
Geometry_optimizer = 'LBFGS'      # QuasiNewton, GPMin, LBFGS or FIRE
Geometry_force_tolerance = 0.05 	 # Maximum force tolerance in LBFGS geometry optimization. Unit is eV/Ang.
Geometry_max_step = 0.2    # How far is a single atom allowed to move. Default is 0.2 Ang.
Fix_symmetry = True    # True for preserving the spacegroup symmetry during optimisation
# Which components of strain will be relaxed: EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
Relax_cell=[True, True, True, False, False, False]
Hydrostatic_pressure=0.0 #GPa

# ELECTRONIC
Wavefunction_cutoff = 300 	# eV
Ground_kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
Ground_kpts_x = 5 	        # kpoints in x direction
Ground_kpts_y = 5	 	# kpoints in y direction
Ground_kpts_z = 5		# kpoints in z direction
Gamma = True
Band_path = 'LGXG'	    # Brillouin zone high symmetry points
Band_npoints = 40		# Number of points between high symmetry points
Hubbard_U = {}            # Can be used like {'N-2p': 6.0}, for none use {}

XC_calc = 'LDA'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0, EXX, B3LYP

SCF_accuracy = 'normal'  # Portable across GPAW and QE
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05

DOS_npoints = 501        # Number of points
DOS_integration = 'smearing'
DOS_width = 0.1          # Gaussian DOS broadening in eV

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
Refine_grid = 4             # refine grid for all electron density (1, 2 [=default] and 4)
Total_charge = 0.0       # Total charge. Normally 0.0 for a neutral system.

#GENERAL
Energy_min = -5 		# eV. It is the minimum energy value for band structure and DOS figures.
Energy_max = 10  		# eV. It is the maximum energy value for band structure and DOS figures.
# Localization = "tr_TR" # Localization setting for figures. en_UK is default.
