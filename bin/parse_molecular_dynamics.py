#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import argparse

from ase import Atoms

from modules.parse_cp2k import (getCellParameters, getStructures, getEnergies, getForces, getStress)

from modules.io_files import save_axsf

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
# Optional parameters
parser.add_argument('--SaveHistory',
                    action='store_true',
                    required=False,
                    help='Save each step of the optimization as a cif file.')

arg = parser.parse_args()

# Check if the optimization was successful.
with open(os.path.join(arg.output_folder, 'simulation_MolecularDynamics.out'), 'r') as f:
    lines = f.readlines()

CellMatrixList, CellParametersList = getCellParameters(arg.output_folder, arg.FrameworkName)
atomLabelList, atomPosList = getStructures(arg.output_folder, arg.FrameworkName)
ForcesList = getForces(arg.output_folder, arg.FrameworkName)
EnergyList = getEnergies(arg.output_folder, arg.FrameworkName)
StressList = getStress(arg.output_folder, arg.FrameworkName)

# Save the optimization history
if arg.SaveHistory:

    # Save the axsf file with the optimization history
    save_axsf(arg.output_folder, arg.FrameworkName, CellMatrixList, atomLabelList, atomPosList, ForcesList)

    # Create a directory to store the optimization history as cif files
    save_path = os.path.join(arg.output_folder, 'MolecularDynamicsHistory')
    os.makedirs(save_path, exist_ok=True)
    for i in range(len(CellParametersList)):
        tempStructure = Atoms(atomLabelList[i],
                              cell=CellParametersList[i],
                              pbc=(1, 1, 1),
                              positions=atomPosList[i])

        print('Saving structure ' + str(i + 1) + ' of ' + str(len(CellParametersList)))

        tempStructure.write(os.path.join(save_path, arg.FrameworkName + '_MD_' + str(i + 1) + '.cif'))


# Write the final structure to file
tempStructure = Atoms(atomLabelList[-1],
                      cell=CellParametersList[-1],
                      pbc=(1, 1, 1),
                      positions=atomPosList[-1])

print('Saving optimized structure.')

tempStructure.write(os.path.join(arg.output_folder, arg.FrameworkName + '_last' + '.cif'))
