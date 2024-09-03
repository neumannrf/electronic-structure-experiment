#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import argparse

from ase import Atoms

from modules.parse_cp2k import (getCellParameters, getStructures, getEnergies, getForces, getStress)

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

CellParametersList = getCellParameters(arg.output_folder, arg.FrameworkName)
StructureList = getStructures(arg.output_folder, arg.FrameworkName)
EnergyList = getEnergies(arg.output_folder, arg.FrameworkName)
ForcesList = getForces(arg.output_folder, arg.FrameworkName)
StressList = getStress(arg.output_folder, arg.FrameworkName)

# To-Do: Save the history of the optimization


# Write the final structure to file
tempStructure = Atoms(StructureList[-1][0],
                      cell=CellParametersList[-1],
                      pbc=(1, 1, 1),
                      positions=StructureList[-1][1].T)

print('Saving optimized structure.')

tempStructure.write(os.path.join(arg.output_folder, arg.FrameworkName + '_optimized' + '.cif'))
