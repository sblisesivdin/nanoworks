# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Density_calc = True  # Export the electron density.
Wavefunction_cutoff = 340  # Set the plane-wave basis cutoff (eV).
Ground_kpts_density = 2.5  # Choose reciprocal-space sampling density for ground state.
Band_path = 'GXWKG'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 301  # Set the number of DOS energy samples.
DOS_integration = 'tetrahedron'  # Choose how to integrate electronic states.
Spin_calc = True  # Allow different up- and down-spin states.
Magmom_per_atom = 1.0  # Set initial atomic moments (Bohr magnetons).
Projected_band_plot = True  # Show orbital weights on the bands.
Projections = [{'atoms': [0, 1], 'orbital': 'd', 'color': 'red', 'label': 'Cr-d'},  # Select atoms and orbitals for band weights (zero-based atoms).
 {'atoms': [0, 1], 'orbital': 's', 'color': 'orange', 'label': 'Cr-s'},
 {'atoms': [2], 'orbital': 'p', 'color': 'blue', 'label': 'O-p'},
 {'atoms': [2], 'orbital': 's', 'color': 'cyan', 'label': 'O-s'}]
