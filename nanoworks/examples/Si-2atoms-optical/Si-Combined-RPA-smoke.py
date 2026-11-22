# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'Si-combined-RPA-smoke-results'  # Keep the output in this directory.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
Band_calc = True  # Calculate bands along a symmetry path.
Density_calc = True  # Export the electron density.
Optical_calc = True  # Calculate the optical response.
Wavefunction_cutoff = 340  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 2  # Sample ground state along reciprocal x.
Ground_kpts_y = 2  # Sample ground state along reciprocal y.
Ground_kpts_z = 2  # Sample ground state along reciprocal z.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
DOS_npoints = 101  # Set the number of DOS energy samples.
DOS_integration = 'smearing'  # Choose how to integrate electronic states.
DOS_width = 0.3  # Set Gaussian DOS broadening (eV).
Band_path = 'GXWKL'  # Choose the band-structure symmetry path.
Band_npoints = 8  # Sample this many points along the path.
Refine_grid = 2  # Refine the exported electron-density grid.
Opt_calc_type = 'RPA'  # Choose RPA or the excitonic BSE method.
Opt_num_of_bands = 8  # Include occupied and empty optical states.
Opt_FD_smearing = 0.05  # Set optical occupation broadening (eV).
Opt_eta = 0.1  # Set optical spectral broadening (eV).
Opt_domega0 = 0.1  # Set low-energy frequency spacing (eV).
Opt_omega2 = 5.0  # Control the nonlinear frequency grid (eV).
Opt_cut_of_energy = 30  # Set the response-function basis cutoff (eV).
