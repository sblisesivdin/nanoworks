# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Ground_calc = True     # Ground state calculations
Geo_optim = False       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = True         # DOS calculation
Band_calc = True        # Band structure calculation
Density_calc = True    # Calculate the all-electron density?
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
Relax_cell=[False, False, False, False, False, False]
Hydrostatic_pressure=0.0 #GPa

# ELECTRONIC
Cut_off_energy = 340 	# eV
Ground_kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
Ground_kpts_x = 3			    # kpoints in x direction
Ground_kpts_y = 3				# kpoints in y direction
Ground_kpts_z = 3				# kpoints in z direction
Gamma = True
Band_path = 'GXWKG'	    # Brillouin zone high symmetry points
Band_npoints = 40		# Number of points between high symmetry points
Setup_params = {}            # Can be used like {'N': ':p,6.0'}, for none use {}

XC_calc = 'PBE'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0, EXX, B3LYP

SCF_accuracy = 'normal'  # Portable across GPAW and QE
Occupation_scheme = 'fermi-dirac'
Smearing_width = 0.05

DOS_npoints = 301        # Number of points
DOS_width = 0.0          # Width of Gaussian smearing.  Use 0.0 for linear tetrahedron interpolation

Spin_calc = True        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
Refine_grid = 4             # refine grid for all electron density (1, 2 [=default] and 4)
Total_charge = 0.0       # Total charge. Normally 0.0 for a neutral system.
Projected_band_plot = True  # Projected Band Configuration for Cr2O
# Assumption: Atom 0 and Atom 1 = Cr, Atom 2 = O (same with .cif file)

Projections = [
    # ==========================================
    # CHROMIUM (Cr) COMPONENTS
    # ==========================================
    # As a transition metal, the d-band is the most critical for magnetism and conductivity.
    {'atoms': [0, 1], 'orbital': 'd', 'color': 'red',   'label': 'Cr-d'},
    
    # Chromium's s-orbitals generally form wider, more dispersed bands.
    {'atoms': [0, 1], 'orbital': 's', 'color': 'orange','label': 'Cr-s'},

    # ==========================================
    # OXYGEN (O) COMPONENTS
    # ==========================================
    # Oxygen's valence electrons are in p-orbitals, which strongly hybridize with Cr-d.
    {'atoms': [2],    'orbital': 'p', 'color': 'blue',  'label': 'O-p'},
    
    # Oxygen's s-orbitals are usually located in deep valence bands, far below the Fermi level.
    {'atoms': [2],    'orbital': 's', 'color': 'cyan',  'label': 'O-s'}
]

#GENERAL
Energy_min = -5 		# eV. It is the minimum energy value for band structure and DOS figures.
Energy_max = 5  		# eV. It is the maximum energy value for band structure and DOS figures.
Localization = "en_UK"  # Localization setting for figures. en_UK is default.
