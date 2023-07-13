#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import argparse

from ase import Atoms

from modules.calculate_properties import getCellParametersFromOptimization, getStructuresFromOptimization

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
with open(os.path.join(arg.output_folder, 'simulation_Optimization.out'), 'r') as f:
    lines = f.readlines()

optimized = False
for line in lines:
    if 'GEOMETRY OPTIMIZATION COMPLETED' in line:
        optimized = True
        break

if not optimized:
    print('Optimization failed! Exiting...')
    exit(1)

CellParametersList = getCellParametersFromOptimization(arg.output_folder, arg.FrameworkName)
StructureList = getStructuresFromOptimization(arg.output_folder, arg.FrameworkName)

if arg.SaveHistory:
    # Create a directory to store the optimization history
    save_path = os.path.join(arg.output_folder, 'OptimizationHistory')
    os.makedirs(save_path, exist_ok=True)
    for i in range(len(CellParametersList)):
        tempStructure = Atoms(StructureList[i][0],
                              cell=CellParametersList[i],
                              pbc=(1, 1, 1),
                              positions=StructureList[i][1].T)

        print('Saving structure ' + str(i) + ' of ' + str(len(CellParametersList)))

        tempStructure.write(os.path.join(save_path, arg.FrameworkName + '_Optimization_' + str(i) + '.cif'))


# Write the final structure to file
tempStructure = Atoms(StructureList[-1][0],
                      cell=CellParametersList[-1],
                      pbc=(1, 1, 1),
                      positions=StructureList[i][-1].T)

print('Saving optimized structure.')

tempStructure.write(os.path.join(arg.output_folder, arg.FrameworkName + '_optimized' + '.cif'))
