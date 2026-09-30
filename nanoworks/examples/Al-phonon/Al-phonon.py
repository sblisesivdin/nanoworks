import numpy as np

# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Ground_calc = True     # Ground state calculations
Geo_optim = True       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = False         # DOS calculation
Band_calc = False        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Phonon_calc = True
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# GEOMETRY
Geometry_optimizer = 'LBFGS'     # QuasiNewton, GPMin, LBFGS or FIRE
Geometry_force_tolerance = 0.05 	# Maximum force tolerance in LBFGS geometry optimization. Unit is eV/Ang.
Geometry_max_step = 0.2          # How far is a single atom allowed to move. Default is 0.2 Ang.
Fix_symmetry = True    # True for preserving the spacegroup symmetry during optimisation
# Which components of strain will be relaxed: EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
Relax_cell=[True, True, True, False, False, False]
Hydrostatic_pressure=0.0 #GPa

# ELECTRONIC
Wavefunction_cutoff = 700 	# eV
#Ground_kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
Ground_kpts_x = 5			    # kpoints in x direction
Ground_kpts_y = 5				# kpoints in y direction
Ground_kpts_z = 5				# kpoints in z direction
Gamma = True
Band_path = 'LGXWK'	    # Brillouin zone high symmetry points
Band_npoints = 40		# Number of points between high symmetry points
Hubbard_U = {}            # Can be used like {'N-2p': 6.0}, for none use {}
Total_charge = 0.0       # Total charge. Normally 0.0 for a neutral system.

XC_calc = 'PBE'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0, EXX, B3LYP

SCF_accuracy = 'normal'  # Portable across GPAW and QE
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05

#PHONON
Phonon_PW_cutoff = 700
Phonon_kpts_x = 5
Phonon_kpts_y = 5
Phonon_kpts_z = 5
Phonon_qpts_x = 21
Phonon_qpts_y = 21
Phonon_qpts_z = 21
Phonon_supercell = np.diag([2, 2, 2])
Phonon_displacement = 1e-3
Phonon_path = 'GXKGL'
Phonon_npoints = 61
Phonon_acoustic_sum_rule = True
Phonon_thermal_calc = True    # Thermodynamic calculations
Phonon_T_min = 0.0            
Phonon_T_max = 1000.0         
Phonon_T_step = 10.0


#GENERAL
Energy_min = -5 		# eV. It is the minimum energy value for band structure and DOS figures.
Energy_max = 5  		# eV. It is the maximum energy value for band structure and DOS figures.
Localization = "en_UK"  # Localization setting for figures. en_UK is default.
