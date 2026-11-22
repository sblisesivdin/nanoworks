# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

from ase.build import bulk  # Build the structure or supercell.

Engine = 'QE'  # Choose the calculation engine.
Outdirname = 'ZnO-QE-withHubbard-results'  # Keep the output in this directory.
bulk_configuration = bulk('ZnO', 'wurtzite', a=3.25, c=5.2)  # Build the atomic structure without a CIF.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Wavefunction_cutoff = 800  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 5  # Sample ground state along reciprocal x.
Ground_kpts_y = 5  # Sample ground state along reciprocal y.
Ground_kpts_z = 5  # Sample ground state along reciprocal z.
Ground_num_of_bands = 40  # Include this many electronic bands.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
Hubbard_U = {'O-2p': 7.0, 'Zn-3d': 10.0}  # Apply orbital-specific U corrections (eV).
DOS_kpts_x = 7  # Sample DOS along reciprocal x.
DOS_kpts_y = 7  # Sample DOS along reciprocal y.
DOS_kpts_z = 7  # Sample DOS along reciprocal z.
DOS_gamma = True  # Choose Gamma-centered DOS sampling.
DOS_integration = 'tetrahedron'  # Choose how to integrate electronic states.
DOS_npoints = 501  # Set the number of DOS energy samples.
Band_path = 'ALMGAHKG'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
Energy_min = -10  # Set the lower plot energy (eV relative to Ef).
Energy_max = 10  # Set the upper plot energy (eV relative to Ef).
