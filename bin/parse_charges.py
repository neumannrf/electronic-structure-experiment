#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import argparse

from modules.calculate_properties import (get_AtomicPositions,
                                          get_CellParameters,
                                          get_CM5AtomicCharges,
                                          get_DDECAtomicCharges)
from modules.io_files import saveCIF

# Required parameters
parser = argparse.ArgumentParser(description='Create the Chargemol simulation input.')
parser.add_argument('output_folder',
                    type=str,
                    action='store',
                    metavar='OUTPUT_FOLDER',
                    help='Directory for storing output files.')
parser.add_argument('--FrameworkName',
                    type=str,
                    required=True,
                    action='store',
                    metavar='FRAMEWORK_NAME',
                    help='Name of the CIF file describing the nanoporous material structure.')
parser.add_argument('--CM5',
                    action='store_true',
                    required=False,
                    help='Get the CM5 charges calculated by Chargemol.')

arg = parser.parse_args()

# Check if the optimization was successful.
with open(os.path.join(arg.output_folder, 'simulation_SCF.out'), 'r') as f:
    lines = f.readlines()

converged = False
for line in lines:
    if 'SCF run converged in' in line:
        converged = True
        break

if not converged:
    print('SCF failed to converge! Exiting...')
    exit(1)

# Get the DDEC charges
DDEC_Charges = get_DDECAtomicCharges(os.path.join(arg.output_folder, 'DDEC6_even_tempered_net_atomic_charges.xyz'))

# Get the CM5 charges
CM5_Charges = get_CM5AtomicCharges(os.path.join(arg.output_folder, 'valence_cube_DDEC_analysis.output'))

# Get the lattice parameters and atomic positions
cif_filename = os.path.join(arg.output_folder, arg.FrameworkName + '.cif')

CellParameters = get_CellParameters(cif_filename)

AtomicTypes, PosX, PosY, PosZ = get_AtomicPositions(cif_filename)

# Write the DDEC charges to file
saveCIF(arg.FrameworkName + '_DDEC',
        CellParameters,
        AtomicTypes,
        PosX,
        PosY,
        PosZ,
        DDEC_Charges,
        arg.output_folder)

if arg.CM5:
    # Write the CM5 charges to file
    saveCIF(arg.FrameworkName + '_CM5',
            CellParameters,
            AtomicTypes,
            PosX,
            PosY,
            PosZ,
            CM5_Charges,
            arg.output_folder)
