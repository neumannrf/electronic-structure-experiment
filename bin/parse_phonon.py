#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2024
# SPDX-License-Identifier: Apache2.0

import argparse
import os

import numpy as np
from ase.cell import Cell
from modules.calculate_properties import (calculate_UnitCells,
                                          get_AtomicPositions,
                                          get_CellParameters,
                                          get_forces,
                                          get_pol_tensor,
                                          get_spg_class,
                                          diff_cross_section,
                                          lorentzian)
from modules.constants import factor2cm
from modules.io_files import save_axsf, saveVibrationalChemicalJSON
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.units import CP2KToTHz

# Required parameters
parser = argparse.ArgumentParser(description='Create symmetric shifts for the small displacement method.')
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
parser.add_argument('--UnitCells',
                    type=str,
                    default=None,
                    metavar='X,Y,Z',
                    required=False,
                    action='store',
                    help='Number of unit cell replications in the supercell (comma-separated).')
parser.add_argument('--PrimitiveMatrix',
                    type=str,
                    default=None,
                    action='store',
                    required=False,
                    metavar='PRIMITIVE_MATRIX',
                    help='Primitive matrix for unit cell creation. (comma-separated row-wise flattened 3x3 matrix)')
parser.add_argument('--dR',
                    type=float,
                    default=0.001,
                    action='store',
                    required=False,
                    metavar='DR',
                    help='Finite difference step size for numerical differentiation in angstrom.')
parser.add_argument('--symmPrec',
                    type=float,
                    default=1e-3,
                    action='store',
                    required=False,
                    metavar='SYMMPREC',
                    help='Symmetry precision for the force constants calculation.')
parser.add_argument('--LaserWaveLength',
                    type=float,
                    default=532,
                    action='store',
                    required=False,
                    metavar='LASER_WAVELENGTH',
                    help='Laser wavelength in nm. Usually 532 or 785')
parser.add_argument('--ExternalTemperature',
                    type=float,
                    default=300,
                    action='store',
                    required=False,
                    metavar='TEMPERATURE',
                    help='Temperature in Kelvin for the Raman cross-section calculation.')
parser.add_argument('--CalculateRaman',
                    action='store_true',
                    required=False,
                    help='Calculate the Raman spectrum.')
parser.add_argument('--CalculateIR',
                    action='store_true',
                    required=False,
                    help='Calculate the IR spectrum.')
parser.add_argument('--HalfWidth',
                    type=float,
                    default=5.0,
                    action='store',
                    metavar='HALF_WIDTH',
                    help='Half width of the Lorentzian function used to broaden the vibrational frequencies.')
parser.add_argument('--Resolution',
                    type=float,
                    default=0.1,
                    action='store',
                    metavar='RESOLUTION',
                    help='Resolution of the Raman/IR spectrum in cm-1.')
parser.add_argument('--CurveLimits',
                    type=str,
                    default='0,4000',
                    action='store',
                    metavar='CURVE_LIMITS',
                    help='Limits for the Raman/IR spectrum plot separated by comma. Ex. 0,4000')
parser.add_argument('--SaveVecs',
                    action='store_true',
                    required=False,
                    help='Save the normal mode displacement vectors.')

# Parse the arguments
arg = parser.parse_args()

# Read the cif file and get the lattice parameters and atomic positions
cif_filename = os.path.join(arg.output_folder, arg.FrameworkName + '.cif')

CellParameters = get_CellParameters(cif_filename)
AtomicTypes, PosX, PosY, PosZ = get_AtomicPositions(cif_filename)

fracPos = np.array([PosX, PosY, PosZ]).T
aseCell = Cell.fromcellpar(CellParameters)

if arg.UnitCells is None:
    arg.UnitCells = calculate_UnitCells(cif_filename, 6).replace(' ', ',')
    print('Supercell size not specified. Using a default value of 6A for the supercell size.')
    print('Calculated supercell size:', arg.UnitCells)

arg.UnitCells = np.array([int(i) for i in arg.UnitCells.split(',')])

if arg.PrimitiveMatrix is None:
    arg.PrimitiveMatrix = 'auto'
else:
    arg.PrimitiveMatrix = np.array([int(i) for i in arg.PrimitiveMatrix.split(',')]).shape(3, 3)

# Get the cell matrix
CellMatrix = aseCell.tolist()

# Create the phonon object
unitcell = PhonopyAtoms(symbols=AtomicTypes,
                        cell=CellMatrix,
                        scaled_positions=fracPos)

PhononCalc = Phonopy(unitcell,
                     supercell_matrix=np.eye(3),
                     primitive_matrix=arg.PrimitiveMatrix,
                     factor=CP2KToTHz,
                     symprec=1e-5)

atomTypes = {
    'primitive': PhononCalc.primitive.get_chemical_symbols(),
    'supercell': PhononCalc.supercell.get_chemical_symbols()
    }

atomNumbers = {
    'primitive': PhononCalc.primitive.get_atomic_numbers(),
    'supercell': PhononCalc.supercell.get_atomic_numbers()
    }

atomMasses = {
    'primitive': PhononCalc.primitive.get_masses(),
    'supercell': PhononCalc.supercell.get_masses()
    }

nAtoms = {
    'primitive': PhononCalc.primitive.get_number_of_atoms(),
    'supercell': PhononCalc.supercell.get_number_of_atoms()
    }

cellMatrix = {
    'primitive': PhononCalc.primitive.get_cell(),
    'supercell': PhononCalc.supercell.get_cell()
    }

cellParameters = {
    'primitive': Cell(cellMatrix['primitive']).cellpar(),
    'supercell': Cell(cellMatrix['supercell']).cellpar()
    }

fracPos = {
    'primitive': PhononCalc.primitive.get_scaled_positions(),
    'supercell': PhononCalc.supercell.get_scaled_positions()
    }

cartPos = {
    'primitive': PhononCalc.primitive.get_positions(),
    'supercell': PhononCalc.supercell.get_positions()
    }

cellVolume = {
    'primitive': np.dot(cellMatrix['primitive'][0], np.cross(cellMatrix['primitive'][1], cellMatrix['primitive'][2])),
    'supercell': np.dot(cellMatrix['supercell'][0], np.cross(cellMatrix['supercell'][1], cellMatrix['supercell'][2]))
    }

IndAtoms = {
    'primitive': [int(i) for i in
                  (PhononCalc.symmetry.get_independent_atoms() / (nAtoms['supercell'] / nAtoms['primitive']))],
    'supercell': PhononCalc.symmetry.get_independent_atoms()
}

invAtomicMass = {
    'primitive': np.sqrt(np.reciprocal([atomMasses['primitive'][i] for i in range(nAtoms['primitive'])])),
    'supercell': np.sqrt(np.reciprocal([atomMasses['supercell'][i] for i in range(nAtoms['supercell'])]))
}

# Get the space group
spaceGroupString, spaceGrounNumber = PhononCalc.symmetry.get_international_table().split()
spaceGrounNumber = int(spaceGrounNumber[1:-1])
spaceGroupClass = get_spg_class(spaceGrounNumber)

# Print the results
cell_txt = "{:7.4f}A  {:7.4f}A {:7.4f}A {:5.2f}° {:5.2f}° {:5.2f}°"
print(f"Primitive cell with {nAtoms['primitive']} atoms:", cell_txt.format(*cellParameters['primitive']))
print(f"Supercell with {nAtoms['supercell']} atoms:", cell_txt.format(*cellParameters['supercell']))
print(f'Found space group: {spaceGroupClass} {spaceGroupString} with number {spaceGrounNumber}')
print(f"{len(IndAtoms['supercell'])} independent atoms on supercell:")

for atom in IndAtoms['supercell']:
    print("    Atom {:3} with type {:2} at position {:7.4f}  {:7.4f}  {:7.4f}".format(atom,
                                                                                      atomTypes['supercell'][atom],
                                                                                      *fracPos['supercell'][atom]))

print('Reading the symmetric shifts for the small displacement method...')

displacement_dict = {
    'natom': nAtoms['supercell'],
    'first_atoms': []
    }

for i in IndAtoms['primitive']:
    for label, dirVec in [['x', [1, 0, 0]], ['y', [0, 1, 0]], ['z', [0, 0, 1]]]:
        for sl, s in [['+', +1], ['-', -1]]:
            dirVec = np.array(dirVec) * s

            outName = f'{arg.FrameworkName}_{i}_{sl}{label}-forces-1_0.xyz'
            forces = get_forces(outName, f'{i}_{sl}{label}')

            displacement_dict['first_atoms'].append({
                'number': i,
                'displacement': arg.dR * dirVec,
                'forces': forces
                })

# Displacement dict
PhononCalc.dataset = displacement_dict

# Set the force constants
PhononCalc.produce_force_constants(calculate_full_force_constants=True)

# Symmetrize the force constants
PhononCalc.symmetrize_force_constants()

# Perform irreducible representation analysis at Gamma
ir_labels = [""] * nAtoms['primitive'] * 3
PhononCalc.supercell.set_magnetic_moments(None)
PhononCalc.set_irreps([0.0, 0.0, 0.0])

for (deg_set, ir) in zip(PhononCalc.get_irreps()._degenerate_sets, PhononCalc.get_irreps()._ir_labels):
    for j in deg_set:
        if (ir is None):
            ir_labels[j] = "Non"
        else:
            ir_labels[j] = ir

# Get the dynamical matrix at Gamma
dymMatrix_Gamma = PhononCalc.get_dynamical_matrix_at_q([0, 0, 0])

# Calculate the eigenvalues (frequencies) and eigenvectors (normal displacements)
eigenValues, eigenVectors = np.linalg.eigh(dymMatrix_Gamma)

# Calculate sqrt(omega) keeping the sign of the eigenvalues and converting to cm-1
frequencies = np.sqrt(np.abs(eigenValues.real)) * np.sign(eigenValues) * factor2cm

# Check if there is any negative frequency
if np.any(frequencies < -5e-3):
    print('WARNING: Negative frequencies found!')
    for i, freq in enumerate(frequencies):
        if freq < -5e-3:
            print(f'Mode {i:3} {ir_labels[i]:4}: {freq:8.2f} cm-1')

if arg.SaveVecs:
    shiftVecs = np.zeros((len(frequencies), nAtoms['primitive'], 3))
    for i, mode in enumerate(eigenVectors.real.T):
        for j, atom in enumerate(mode.reshape(-1, 3)):
            # Normalize the displacement vectors by the square root of the atomic mass
            shiftVecs[i, j] = atom / np.sqrt(atomMasses['primitive'][j])

    os.makedirs(os.path.join(arg.output_folder, 'VIBRATION_FILES'), exist_ok=True)

    # Save independend files for each mode
    for i, freq in enumerate(frequencies):
        save_axsf(os.path.join(arg.output_folder, 'VIBRATION_FILES'),
                  f'{arg.FrameworkName}_{i}_{ir_labels[i]}_{freq}',
                  [cellMatrix['primitive']],
                  [atomTypes['primitive']],
                  [cartPos['primitive']],
                  [shiftVecs[i]])

    # Save all modes in a single file
    save_axsf(arg.output_folder,
              f'{arg.FrameworkName}_all',
              [cellMatrix['primitive'] for _ in range(len(shiftVecs))],
              [atomTypes['primitive'] for _ in range(len(shiftVecs))],
              [cartPos['primitive'] for _ in range(len(shiftVecs))],
              shiftVecs)

# Calculating dP/dR: N atoms, 3 directions (x, y, z), 2 polarizations (+, -), 3x3 tensor
d_polarizability = np.zeros((nAtoms['primitive'], 3, 3, 3))

# Iterate over atoms on primitive cell
for i in range(nAtoms['primitive']):
    for d, label, _ in [[0, 'x', [1, 0, 0]], [1, 'y', [0, 1, 0]], [2, 'z', [0, 0, 1]]]:
        # Get the polarizability tensor in Angstrom^2
        polTensor_plus = get_pol_tensor(f'{arg.FrameworkName}_{i}_+{label}-raman-1_0.data', f'{i}_+{label}')[1]
        polTensor_minus = get_pol_tensor(f'{arg.FrameworkName}_{i}_-{label}-raman-1_0.data', f'{i}_-{label}')[1]

        # Use the two-point finite difference formula to calculate the polarizability tensor derivatives
        # f'(x) = (f(x + h) - f(x - h)) / (2 * h)
        d_polarizability[i][d] = (polTensor_plus - polTensor_minus) / (2 * arg.dR)

# Get the phonon eigendisplacements vectors as (3 * N_atoms) x N_modes array
phonon_eigendisplacements = np.zeros((nAtoms['primitive'], len(frequencies), 3))

# Fill the phonon eigendisplacements array
for i, mode in enumerate(eigenVectors.real.T):
    for j, atom in enumerate(mode.reshape(-1, 3)):
        phonon_eigendisplacements[j, i] = atom

# Calculating Raman tensor from polarizability tensor derivatives (dP/dR), phonon eigendisplacements, and atomic masses
alpha = np.einsum('ad...,akd,a->k...',
                  d_polarizability,
                  phonon_eigendisplacements,
                  invAtomicMass['primitive']) * np.sqrt(cellVolume['primitive'])

# Calculating Raman Tensor Placzek Invariants following:
# The Raman Effect: A Unified Treatment of the Theory of Raman Scattering by Molecules
# by Derek A. Long, 2002
# This does not require the tensor to be symmetric.

# Calculate the mean polarizability squared
a_sq = np.square(np.trace(alpha, 0, 2) / 3).reshape((-1, 1))

# Create an empty vector for the anisotropy
gamma_sq = np.zeros((len(frequencies), 1))

# Create an empty vector for asymmetric anisotropy
delta_sq = np.zeros_like(gamma_sq)

for k in range(len(frequencies)):
    delta_sq[k] = 3/4 * (np.square(alpha[k][0][1] - alpha[k][1][0])
                         + np.square(alpha[k][1][2] - alpha[k][2][1])
                         + np.square(alpha[k][2][0] - alpha[k][0][2]))

    gamma_sq[k] = 1/2 * (np.square(alpha[k][0][0] - alpha[k][1][1])
                         + np.square(alpha[k][1][1] - alpha[k][2][2])
                         + np.square(alpha[k][2][2] - alpha[k][0][0])) \
        + 3/4 * (np.square(alpha[k][0][1] + alpha[k][1][0])
                 + np.square(alpha[k][0][2] + alpha[k][2][0])
                 + np.square(alpha[k][1][2] + alpha[k][2][1]))

# Create the Raman Intensities vector
I_raman = np.zeros((len(frequencies), 3))

# Calculate absolute Raman intensity: Total, Perpendicular, and Parallel considering
# incident linear polarized radiation
for k in range(len(frequencies)):
    I_total = 45 * a_sq[k] + 7 * gamma_sq[k] + 5 * delta_sq[k]
    I_parallel = 45 * a_sq[k] + 4 * gamma_sq[k]
    I_perpendicular = 3 * gamma_sq[k] + 5 * delta_sq[k]

    I_raman[k] = np.array([I_total, I_perpendicular, I_parallel]).flatten() / 45

cs = diff_cross_section(frequencies, arg.LaserWaveLength, arg.ExternalTemperature)

# Create the Raman cross section vector
raman_cross_section = np.zeros((len(frequencies), 3))

for k in range(len(frequencies)):
    raman_cross_section[k] = cs[k] * I_raman[k]

# Prepare the Raman data to save as a csv file
raman_data = [[i, ir_labels[i], freq, *I_raman[i], *raman_cross_section[i]] for i, freq in enumerate(frequencies)]

header_list = [
    'Mode',
    'Symmetry',
    'Frequency (cm-1)',
    'Total Raman Int (a.u)',
    'Perpendicular Raman Int (a.u)',
    'Parallel Raman Int (a.u)',
    'Total Cross Section (Å^4.amu^-1)',
    'Perpendicular Cross Section (Å^4.amu^-1)',
    'Parallel Cross Section (Å^4.amu^-1)']

# Save raman_data as a csv file
np.savetxt(os.path.join(arg.output_folder, f'{arg.FrameworkName}_RamanTable.csv'),
           np.array(raman_data, dtype=object),
           header=','.join(header_list),
           delimiter=',',
           fmt='%5d,%4s,%10.2f,%15.5e,%15.5e,%15.5e,%15.5e,%15.5e,%15.5e')


curve_limits = [int(i) for i in arg.CurveLimits.split(',')]

# Calculate the Raman spectrum
X = np.arange(round(min(frequencies)) - 100, max(frequencies) + 100, arg.Resolution)
I_tot = np.zeros_like(X)
I_perp = np.zeros_like(X)
I_par = np.zeros_like(X)
Cs_tot = np.zeros_like(X)
Cs_perp = np.zeros_like(X)
Cs_par = np.zeros_like(X)


for i, freq in enumerate(frequencies):
    # Skip the frequencies outside the curve limits
    if freq < curve_limits[0] or freq > curve_limits[1]:
        continue

    I_tot += lorentzian(X, freq, arg.HalfWidth) * I_raman[i][0]
    I_perp += lorentzian(X, freq, arg.HalfWidth) * I_raman[i][1]
    I_par += lorentzian(X, freq, arg.HalfWidth) * I_raman[i][2]
    Cs_tot += lorentzian(X, freq, arg.HalfWidth) * raman_cross_section[i][0]
    Cs_perp += lorentzian(X, freq, arg.HalfWidth) * raman_cross_section[i][1]
    Cs_par += lorentzian(X, freq, arg.HalfWidth) * raman_cross_section[i][2]

# Normalize the Raman intensities
norm_factor = np.max(I_tot)

I_tot /= norm_factor
I_perp /= norm_factor
I_par /= norm_factor

norm_factor = np.max(Cs_tot)

Cs_tot /= norm_factor
Cs_perp /= norm_factor
Cs_par /= norm_factor

header_list = ['Frequency (cm-1)',
               'Total Raman Int (a.u)',
               'Perpendicular Raman Int (a.u)',
               'Parallel Raman Int (a.u)',
               'Total Cross Section (a.u.)',
               'Perpendicular Cross Section (a.u.)',
               'Parallel Cross Section (a.u.)']

# Save as a numpy csv file
np.savetxt(os.path.join(arg.output_folder, f'{arg.FrameworkName}_RAMAN_Curve.csv'),
           np.transpose([X, I_tot, I_perp, I_par, Cs_tot, Cs_perp, Cs_par]),
           header=','.join(header_list),
           delimiter=',',
           fmt='%15.7f')

# Save the vibrations as cjson file
saveVibrationalChemicalJSON(OutputFolder=arg.output_folder,
                            Frameworkname=arg.FrameworkName,
                            CellParameters=cellParameters['primitive'],
                            atomTypes=atomTypes['primitive'],
                            cartPos=cartPos['primitive'],
                            modes=ir_labels,
                            eigenVectors=np.array([i.flatten() for i in shiftVecs]).tolist(),
                            freqList=frequencies,
                            IR_intensity=np.zeros(len(frequencies)),
                            RAMAN_intensity=I_raman.T[0])
