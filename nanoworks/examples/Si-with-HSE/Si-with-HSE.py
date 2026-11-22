# SPDX-FileCopyrightText: Sefer Bora Lisesivdin and Beyza Lisesivdin
# SPDX-License-Identifier: MIT
# See LICENSE.md in the project root for license terms.

Engine = 'GPAW'  # Choose the calculation engine.
Ground_calc = True  # Calculate the electronic ground state.
DOS_calc = True  # Calculate the density of electronic states.
XC_calc = 'HSE06'  # Choose the exchange-correlation functional.
Wavefunction_cutoff = 400  # Set the plane-wave basis cutoff (eV).
Ground_kpts_x = 3  # Sample ground state along reciprocal x.
Ground_kpts_y = 3  # Sample ground state along reciprocal y.
Ground_kpts_z = 3  # Sample ground state along reciprocal z.
SCF_accuracy = 'loose'  # Set electronic convergence accuracy.
Occupation_scheme = 'fermi-dirac'  # Choose how electronic states are occupied.
Smearing_width = 0.2  # Set occupation broadening (eV).
DOS_npoints = 301  # Set the number of DOS energy samples.
DOS_integration = 'smearing'  # Choose how to integrate electronic states.
DOS_width = 0.2  # Set Gaussian DOS broadening (eV).
