# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'Si-results'  # Keep the output in this directory.
Optical_calc = True  # Calculate the optical response.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Opt_calc_type = 'BSE'  # Choose RPA or the excitonic BSE method.
Opt_shift_en = 0.0  # Set a rigid optical energy shift (eV).
Opt_BSE_valence = range(0, 4)  # Select occupied bands for BSE (zero-based).
Opt_BSE_conduction = range(4, 8)  # Select empty bands for BSE (zero-based).
Opt_num_of_bands = 8  # Include occupied and empty optical states.
Opt_min_en = 0.0  # Start the photon-energy range (eV).
Opt_max_en = 20.0  # End the photon-energy range (eV).
Opt_num_of_data = 1001  # Set the number of photon-energy samples.
Opt_FD_smearing = 0.05  # Set optical occupation broadening (eV).
Opt_eta = 0.05  # Set optical spectral broadening (eV).
Opt_cut_of_energy = 100  # Set the response-function basis cutoff (eV).
Opt_kpts_x = 4  # Sample optics along reciprocal x.
Opt_kpts_y = 4  # Sample optics along reciprocal y.
Opt_kpts_z = 4  # Sample optics along reciprocal z.
Opt_gamma = True  # Choose Gamma-centered optical sampling.
