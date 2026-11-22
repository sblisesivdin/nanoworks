# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'Si-results'  # Keep the output in this directory.
Optical_calc = True  # Calculate the optical response.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Opt_calc_type = 'RPA'  # Choose RPA or the excitonic BSE method.
Opt_num_of_bands = 8  # Include occupied and empty optical states.
Opt_FD_smearing = 0.05  # Set optical occupation broadening (eV).
Opt_eta = 0.05  # Set optical spectral broadening (eV).
Opt_domega0 = 0.05  # Set low-energy frequency spacing (eV).
Opt_omega2 = 5.0  # Control the nonlinear frequency grid (eV).
Opt_cut_of_energy = 100  # Set the response-function basis cutoff (eV).
Opt_nblocks = None  # Distribute response matrices automatically when None.
Opt_kpts_x = 4  # Sample optics along reciprocal x.
Opt_kpts_y = 4  # Sample optics along reciprocal y.
Opt_kpts_z = 4  # Sample optics along reciprocal z.
Opt_gamma = True  # Choose Gamma-centered optical sampling.
