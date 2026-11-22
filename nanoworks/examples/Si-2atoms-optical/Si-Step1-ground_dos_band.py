# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'Si-results'  # Keep the output in this directory.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Density_calc = True  # Export the electron density.
Wavefunction_cutoff = 340  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 4  # Sample ground state along reciprocal x.
Ground_kpts_y = 4  # Sample ground state along reciprocal y.
Ground_kpts_z = 4  # Sample ground state along reciprocal z.
Band_path = 'GXWKL'  # Choose the band-structure symmetry path.
Band_npoints = 40  # Sample this many points along the path.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 501  # Set the number of DOS energy samples.
DOS_integration = 'smearing'  # Choose how to integrate electronic states.
DOS_width = 0.3  # Set Gaussian DOS broadening (eV).
