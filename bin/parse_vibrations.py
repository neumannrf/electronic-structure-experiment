#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import argparse
import numpy as np

from modules.calculate_properties import get_vibrational_data, lorentzian, saveVibrationalVectors

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
parser.add_argument('--HalfWidth',
                    type=float,
                    default=5.0,
                    action='store',
                    metavar='HALF_WIDTH',
                    help='Half width of the Lorentzian function used to broaden the vibrational frequencies.')
parser.add_argument('--SaveVibrations',
                    action='store_true',
                    required=False,
                    help='Save the vibrational modes as AXSF files.')

arg = parser.parse_args()

frequency, IR_intensity, RAMAN_intensity = get_vibrational_data('simulation_Vibrations.out')

if len(IR_intensity) == 0:
    IR_intensity = np.zeros_like(frequency)

if len(RAMAN_intensity) == 0:
    RAMAN_intensity = np.zeros_like(frequency)

# Save the vibrational frequencies and intensities as a numpy csv file
np.savetxt(os.path.join(arg.output_folder, f'{arg.FrameworkName}_VibrationsTable.csv'),
           np.transpose([frequency, IR_intensity, RAMAN_intensity]),
           header='Frequency (cm-1), IR Intensity (km/mol), Raman Intensity (A2/amu)',
           delimiter=',')

X = np.linspace(min(frequency), max(frequency)*1.2, 10000)
IR_curve = np.zeros_like(X)
RAMAN_curve = np.zeros_like(X)

for i, freq in enumerate(frequency):
    IR_curve += lorentzian(X, freq, arg.HalfWidth) * IR_intensity[i]
    RAMAN_curve += lorentzian(X, freq, arg.HalfWidth) * RAMAN_intensity[i]

# Save as a numpy csv file
np.savetxt(os.path.join(arg.output_folder, f'{arg.FrameworkName}_RAMAN_IR_Curve.csv'),
           np.transpose([X, IR_curve, RAMAN_curve]),
           header='Frequency (cm-1), IR Intensity (km/mol), Raman Intensity (A2/amu)',
           delimiter=',')

if arg.SaveVibrations:
    saveVibrationalVectors(arg.FrameworkName, arg.output_folder)
