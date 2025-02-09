import numpy as np

# Define physical constants

c = 299792458                      # Speed of light in m/s
mu_0 = 1.25663706212 * 1e-6        # Permeability of free space in m kg s^-2 A^-2
eps_0 = 1 / (np.square(c) * mu_0)  # Permittivity of free space in m^-3 kg^-1 s^4 A^2
k_B = 1.380649 * 1e-23             # Boltzmann constant in m^2 kg s^-2 K^-1
h = 6.62607015 * 1e-34             # Planck constant in m^2 kg s^-1
hbar = h / (2 * np.pi)             # Reduced Planck constant in m^2 kg s^-1

# Define conversion factors

Angst2Bohr = 1.889725989      # Angstrom to Bohr
Bohr2Angst = 1 / Angst2Bohr   # Bohr to Angstrom
factor2cm = 3739.4256800756   # Convert from phononpy to cm-1
nm2cm = 1e-7                  # Convert from nm to cm-1

header = r"""
 *******************************************************************************
 **                       ______     _____     ______                         **
 **                      |  ____|   / ____|   |  ____|                        **
 **                      | |__     | (___     | |__                           **
 **                      |  __|     \___ \    |  __|                          **
 **                      | |____    ____) |   | |____                         **
 **                      |______|  |_____/    |______|                        **
 **                                                                           **
 **                     Electronic Structure Experiments                      **
 **           A automated python tool to run DFT Simulations with CP2K        **
 **                                                                           **
 **      {:^65}    **
 **                                                                           **
 **          "Simulating molecules and materials, one electron at a time!"    **
 **                                                                           **
 *******************************************************************************
"""
