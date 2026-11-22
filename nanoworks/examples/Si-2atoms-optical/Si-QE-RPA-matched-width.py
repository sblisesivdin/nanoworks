# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'QE'  # Choose the calculation engine.
Outdirname = 'Si-QE-RPA-matched-width-results'  # Keep the output in this directory.
Ground_calc = True  # Calculate the electronic ground state.
Optical_calc = True  # Calculate the optical response.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Wavefunction_cutoff = 340  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 2  # Sample ground state along reciprocal x.
Ground_kpts_y = 2  # Sample ground state along reciprocal y.
Ground_kpts_z = 2  # Sample ground state along reciprocal z.
Ground_num_of_bands = 8  # Include this many electronic bands.
Gamma = True  # Center the k-point mesh at Gamma.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
Total_charge = 0.0  # Remove this many electrons when positive.
Opt_calc_type = 'RPA'  # Choose RPA or the excitonic BSE method.
Opt_num_of_bands = 8  # Include occupied and empty optical states.
Opt_kpts_x = 2  # Sample optics along reciprocal x.
Opt_kpts_y = 2  # Sample optics along reciprocal y.
Opt_kpts_z = 2  # Sample optics along reciprocal z.
Opt_gamma = True  # Choose Gamma-centered optical sampling.
Opt_FD_smearing = 0.05  # Set optical occupation broadening (eV).
Opt_eta = 0.1  # Set optical spectral broadening (eV).
Opt_min_en = 0.0  # Start the photon-energy range (eV).
Opt_max_en = 10.0  # End the photon-energy range (eV).
Opt_num_of_data = 1001  # Set the number of photon-energy samples.
Opt_shift_en = 0.0  # Set a rigid optical energy shift (eV).
