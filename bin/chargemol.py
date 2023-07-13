#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import argparse
import os
from textwrap import dedent

from modules.calculate_properties import get_CellParameters, get_AtomicPositions
from modules.atom_data import ATOMIC_NUMBER, CORE_NUM

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
parser.add_argument('--NetCharge',
                    type=float,
                    default=0.0,
                    action='store',
                    required=False,
                    metavar='NET_CHARGE',
                    help='Net charge of the unit cell.')
parser.add_argument('--DataFolder',
                    type=str,
                    default=os.environ.get("CHARGEMOL_DATA_FOLDER"),
                    action='store',
                    required=False,
                    metavar='DATA_FOLDER',
                    help='Directory containing the atomic density files for Chargemol.')

arg = parser.parse_args()

# Calculate self-consistent properties
cif_filename = os.path.join(arg.output_folder, arg.FrameworkName + '.cif')

# Get the lattice parameters
arg.Cell = get_CellParameters(cif_filename)

# Get the atomic positions
AtomLabels, _, _, _ = get_AtomicPositions(cif_filename)

# Get the unique elements in the unit cell
Kinds = list(set(AtomLabels))

# Create string with the number of core electrons for each element
CoreTxt = [f'{ATOMIC_NUMBER[atom]} {CORE_NUM[atom]}' for atom in Kinds]
arg.CoreTxt = '\n'.join(CoreTxt)

# Create input file for CP2K Energy calculation
inputfile = dedent("""\
<net charge>
{NetCharge} <-- specifies the net charge of the unit cell (defaults to 0.0 if nothing specified)
</net charge>

<periodicity along A, B, and C vectors>
.true.
.true.
.true.
</periodicity along A, B, and C vectors>

<atomic densities directory complete path>
{DataFolder}
</atomic densities directory complete path>

<charge type>
DDEC6 <-- specifies the charge type (DDEC3 or DDEC6)
</charge type

<compute BOs>
.false. <-- specifies whether to compute bond orders or not
</compute BOs>

<number of core electrons>
{CoreTxt}
</number of core electrons>
""").format(**arg.__dict__)

# Write string to file
with open(os.path.join(arg.output_folder, 'job_control.txt'), 'w') as f:
    f.write(inputfile)
