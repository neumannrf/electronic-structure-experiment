#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2024
# SPDX-License-Identifier: Apache2.0

import argparse
import os
from copy import deepcopy

import numpy as np
from ase.cell import Cell
from modules.constants import header
from modules.calculate_properties import (calculate_UnitCells,
                                          create_input_file,
                                          get_AtomicPositions,
                                          get_CellParameters, get_spg_class)
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
                    default=1e-5,
                    action='store',
                    required=False,
                    metavar='SYMM_PREC',
                    help='Symmetry precision used in the phonopy calculations.')
parser.add_argument('--Charge',
                    type=int,
                    default=0,
                    action='store',
                    required=False,
                    metavar='CHARGE',
                    help='Total charge of the unit cell.')
parser.add_argument('--Multiplicity',
                    type=int,
                    default=1,
                    action='store',
                    required=False,
                    metavar='MULTIPLICITY',
                    help='Total multiplicity of the unit cell.')
parser.add_argument('--MaxSCFcycles',
                    type=int,
                    default=50,
                    action='store',
                    required=False,
                    metavar='MAX_SCF_CYCLES',
                    help='Maximum number of SCF cycles.')
parser.add_argument('--UseOT',
                    action='store_true',
                    required=False,
                    help='Use the Orbital Transformation (OT) method.')
parser.add_argument('--MaxOuterSCFycles',
                    type=int,
                    default=5,
                    action='store',
                    required=False,
                    metavar='MAX_OUTER_SCF_CYCLES',
                    help='Maximum number of Outer SCF cycles for OT simulations.')
parser.add_argument('--SCFGuess',
                    type=str,
                    default='atomic',
                    action='store',
                    required=False,
                    choices=['atomic',
                             'restart',
                             'core',
                             'random',
                             'sparse',
                             'mopac'],
                    metavar='SCF_GUESS',
                    help='Initial guess for the SCF cycle.')
parser.add_argument('--MixingMethod',
                    type=str,
                    default='direct_p_mixing',
                    action='store',
                    required=False,
                    choices=['direct_p_mixing',
                             'broyden_mixing',
                             'broyden_mixing_new',
                             'kerker_mixing'],
                    metavar='MIXING_METHOD',
                    help='Method for mixing the density matrix.')
parser.add_argument('--MixingAlpha',
                    type=float,
                    default=0.2,
                    action='store',
                    required=False,
                    metavar='MIXING_ALPHA',
                    help='Mixing parameter for the density matrix.')
parser.add_argument('--EPSDefault',
                    type=float,
                    default=1e-10,
                    action='store',
                    required=False,
                    metavar='EPS_DEFAULT',
                    help='Default value for the electronic density convergence threshold.')
parser.add_argument('--PWCutoff',
                    type=float,
                    default=800,
                    action='store',
                    required=False,
                    metavar='PW_CUTOFF',
                    help='Plane wave cutoff energy in Ry.')
parser.add_argument('--NGrid',
                    type=int,
                    default=5,
                    action='store',
                    required=False,
                    metavar='N_GRID',
                    help='Number of grids for the multigrid method.')
parser.add_argument('--RelativeCutOff',
                    type=float,
                    default=60,
                    action='store',
                    required=False,
                    metavar='RELATIVE_CUTOFF',
                    help='Relative cutoff for the multigrid method.')
parser.add_argument('--Functional',
                    type=str,
                    default='PBE',
                    action='store',
                    required=False,
                    choices=['PBE', 'XTB'],
                    metavar='FUNCTIONAL',
                    help='Functional used to calculate the total energy.')
parser.add_argument('--Parametrization',
                    type=str,
                    default='ORIG',
                    action='store',
                    required=False,
                    choices=['ORIG', 'PBESOL', 'REVPBE'],
                    metavar='PARAMETRIZATION',
                    help='PBE functional parametrization used to calculate the total energy.')
parser.add_argument('--DispersionCorrection',
                    type=str,
                    default='DFTD3',
                    action='store',
                    required=False,
                    choices=['DFTD3', 'DFTD3(BJ)'],
                    metavar='DISPERSION_CORRECTION',
                    help='Dispersion correction used to calculate the total energy')
parser.add_argument('--BasisSet',
                    type=str,
                    default='DZVP',
                    action='store',
                    required=False,
                    choices=['SZV', 'DZVP', 'TZVP', 'TZV2P'],
                    metavar='BASIS_SET',
                    help='Gaussian basis set type.')
parser.add_argument('--SCFConvergence',
                    type=float,
                    default=1e-8,
                    action='store',
                    required=False,
                    metavar='SCF_CONVERGENCE',
                    help='SCF convergence threshold.')
parser.add_argument('--CP2KDataDir',
                    type=str,
                    default=os.environ.get("CP2K_DATA_DIR"),
                    action='store',
                    required=False,
                    metavar='CP2K_DATA_DIR',
                    help='Directory containing the Basis set and pseudopotential files.')
parser.add_argument('--CalculateRaman',
                    action='store_true',
                    required=False,
                    help='Calculate the Raman spectrum.')
parser.add_argument('--CalculateIR',
                    action='store_true',
                    required=False,
                    help='Calculate the IR spectrum.')
parser.add_argument('--UseScalapack',
                    action='store_true',
                    required=False,
                    help='Use Scalapack as preferred diagonalization library')


# Parse the arguments
arg = parser.parse_args()

print(header.format('Small Displacement input creation'))

# Read the cif file and get the lattice parameters and atomic positions
cif_filename = arg.FrameworkName + '.cif'

print(cif_filename)
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
                     supercell_matrix=np.eye(3) * arg.UnitCells,
                     primitive_matrix=arg.PrimitiveMatrix,
                     factor=CP2KToTHz,
                     symprec=arg.symmPrec)

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

print('Generating the symmetric shifts for the small displacement method...')
print('{} shited structures will be generated.'.format(nAtoms['supercell'] * 6))

for i in range(nAtoms['supercell']):
    for label, dirVec in [['x', [1, 0, 0]], ['y', [0, 1, 0]], ['z', [0, 0, 1]]]:
        for sl, s in [['+', +1], ['-', -1]]:
            dirVec = np.array(dirVec) * s

            print('Shift of {} on atom {}{} along {}{} direction: {:2} {:2} {:2}'.format(arg.dR,
                                                                                         atomTypes['supercell'][i],
                                                                                         i,
                                                                                         sl,
                                                                                         label,
                                                                                         *dirVec))

            # Shift in the plus direction
            shifted_cartPos = deepcopy(cartPos['supercell'])
            shifted_cartPos[i] = shifted_cartPos[i] + arg.dR * dirVec

            file_name = f'{arg.FrameworkName}_{i}_{sl}{label}'

            os.makedirs(f'{i}_{sl}{label}', exist_ok=True)

            create_input_file(FrameworkName=file_name,
                              output_folder=os.path.join(os.getcwd(), f'{i}_{sl}{label}'),
                              CalcType='energy_force',
                              CellMatrix=cellMatrix['supercell'],
                              AtomicTypes=atomTypes['supercell'],
                              CartX=shifted_cartPos[:, 0],
                              CartY=shifted_cartPos[:, 1],
                              CartZ=shifted_cartPos[:, 2],
                              Charge=arg.Charge,
                              Multiplicity=arg.Multiplicity,
                              UseOT=arg.UseOT,
                              MaxSCFcycles=arg.MaxSCFcycles,
                              MaxOuterSCFycles=arg.MaxOuterSCFycles,
                              SCFGuess=arg.SCFGuess,
                              MixingMethod=arg.MixingMethod,
                              MixingAlpha=arg.MixingAlpha,
                              EPSDefault=arg.EPSDefault,
                              PWCutoff=arg.PWCutoff,
                              NGrid=arg.NGrid,
                              RelativeCutOff=arg.RelativeCutOff,
                              Functional=arg.Functional,
                              Parametrization=arg.Parametrization,
                              DispersionCorrection=arg.DispersionCorrection,
                              BasisSet=arg.BasisSet,
                              CalculateRaman=arg.CalculateRaman,
                              CalculateIR=arg.CalculateIR,
                              UseScalapack=arg.UseScalapack
                              )

print('Done!')
