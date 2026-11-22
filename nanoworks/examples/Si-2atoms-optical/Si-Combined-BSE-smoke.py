# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Outdirname = 'Si-BSE-smoke-results'  # Keep the output in this directory.
Ground_calc = True  # Calculate the electronic ground state.
Optical_calc = True  # Calculate the optical response.
XC_calc = 'PBE'  # Choose the exchange-correlation functional.
Wavefunction_cutoff = 250  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 2  # Sample ground state along reciprocal x.
Ground_kpts_y = 2  # Sample ground state along reciprocal y.
Ground_kpts_z = 2  # Sample ground state along reciprocal z.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.05  # Set occupation broadening (eV).
Opt_calc_type = 'BSE'  # Choose RPA or the excitonic BSE method.
Opt_num_of_bands = 8  # Include occupied and empty optical states.
Opt_BSE_valence = range(0, 4)  # Select occupied bands for BSE (zero-based).
Opt_BSE_conduction = range(4, 8)  # Select empty bands for BSE (zero-based).
Opt_cut_of_energy = 50  # Set the response-function basis cutoff (eV).
Opt_min_en = 0.0  # Start the photon-energy range (eV).
Opt_max_en = 10.0  # End the photon-energy range (eV).
Opt_num_of_data = 101  # Set the number of photon-energy samples.
Opt_eta = 0.1  # Set optical spectral broadening (eV).
Opt_FD_smearing = 0.05  # Set optical occupation broadening (eV).
Opt_shift_en = 0.0  # Set a rigid optical energy shift (eV).
