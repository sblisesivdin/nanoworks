# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase.build import bulk  # Build the structure or supercell.

Engine = 'QE'  # Choose the calculation engine.
Outdirname = 'Fe-QE-spin-smoke-results'  # Keep the output in this directory.
bulk_configuration = bulk('Fe', 'bcc', a=2.87)  # Build the atomic structure without a CIF.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Projected_band_plot = True  # Show orbital weights on the bands.
Density_calc = True  # Export the electron density.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Wavefunction_cutoff = 500  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 4  # Sample ground state along reciprocal x.
Ground_kpts_y = 4  # Sample ground state along reciprocal y.
Ground_kpts_z = 4  # Sample ground state along reciprocal z.
Ground_num_of_bands = 16  # Include this many electronic bands.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.1  # Set occupation broadening (eV).
Spin_calc = True  # Allow different up- and down-spin states.
Magmom_per_atom = 2.2  # Set initial atomic moments (Bohr magnetons).
DOS_integration = 'tetrahedron'  # Choose how to integrate electronic states.
DOS_npoints = 101  # Set the number of DOS energy samples.
Band_path = 'GHNGPH'  # Choose the band-structure symmetry path.
Band_npoints = 16  # Sample this many points along the path.
Projections = [{'atoms': [0], 'orbital': 'd', 'color': 'red', 'label': 'Fe-d'},  # Select atoms and orbitals for band weights (zero-based atoms).
 {'atoms': [0], 'orbital': 's', 'color': 'orange', 'label': 'Fe-s'}]
Energy_min = -8  # Set the lower plot energy (eV relative to Ef).
Energy_max = 8  # Set the upper plot energy (eV relative to Ef).
