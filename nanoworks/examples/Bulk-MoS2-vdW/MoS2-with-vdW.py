# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Outdirname = 'MoS2-with-vdW'

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
SOC_calc = False              # Calculate Spin Orbit Coupling Effects
vdW_calc = 'D3'        # Calculate with vdW correction. 'D3' for Grimme-D3 or 'None' for none.
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
Wavefunction_cutoff = 400 	# eV
Ground_kpts_x = 6 	        # kpoints in x direction
Ground_kpts_y = 6	 	# kpoints in y direction
Ground_kpts_z = 4		# kpoints in z direction
Gamma = True
Band_path = 'GMKG'	    # Brillouin zone high symmetry points

XC_calc = 'PBE'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0, EXX, B3LYP

SCF_accuracy = 'normal'  # Portable across GPAW and QE
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05
