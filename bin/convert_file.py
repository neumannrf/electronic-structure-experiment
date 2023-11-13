#!/usr/bin/env -S python -B

# SPDX-License-Identifier: Apache2.0

import argparse

from modules.io_files import (readCIF, saveCIF,
                              readChemicalJSON, saveChemicalJSON,
                              readGJF, saveGJF,
                              readXSF, saveXSF)

# Required parameters
parser = argparse.ArgumentParser(description='Create the CP2K simulation input.')
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
parser.add_argument('--InputFormat',
                    type=str,
                    required=True,
                    choices=['cif', 'cjson', 'xsf', 'gjf'],
                    action='store',
                    metavar='INPUT_FORMAT',
                    help='Format of the input file.')
parser.add_argument('--OutputFormat',
                    type=str,
                    required=True,
                    action='store',
                    choices=['cif', 'cjson', 'xsf', 'gjf'],
                    metavar='OUTPUT_FORMAT',
                    help='Format of the output file.')
parser.add_argument('--ChargeType',
                    type=str,
                    default='none',
                    action='store',
                    metavar='CHARGE_TYPE',
                    help='Type of partial charges.')

# Parse the arguments
arg = parser.parse_args()

# Dictionary containing the function to read files
read_dict = {
    'cif': readCIF,
    'gjf': readGJF,
    'cjson': readChemicalJSON,
    'xsf': readXSF,
}

# Dictionary containing the function to save files
save_dict = {
    'cif': saveCIF,
    'gjf': saveGJF,
    'cjson': saveChemicalJSON,
    'xsf': saveXSF,
}

CellParameters, labels, frac_x, frac_y, frac_z, charges, charge_type = read_dict[arg.InputFormat](arg.FrameworkName,
                                                                                                  arg.output_folder)

save_dict[arg.OutputFormat](FrameworkName=arg.FrameworkName,
                            CellParameters=CellParameters,
                            labels=labels,
                            frac_x=frac_x,
                            frac_y=frac_y,
                            frac_z=frac_z,
                            charges=charges,
                            charge_type=arg.ChargeType,
                            OutputFolder=arg.output_folder)
