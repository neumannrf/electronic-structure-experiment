#!/usr/bin/env -S python -B

# SPDX-License-Identifier: Apache2.0

import os
import numpy as np
from ase.cell import Cell


def getCellParameters(outputfolder: str, FrameworkName: str) -> list[float]:
    """
    Get the cell parameters from CP2K optimization.

    Parameters
    ----------
    outputfolder : str
        Path to the output folder.
    FrameworkName : str
        Name of the framework.

    Returns
    -------
    cellParameters : list
        List of the cell parameters.
    """

    if not os.path.isfile(os.path.join(outputfolder, FrameworkName + '-1.cell')):
        return [], []

    # Open the FrameworkName-1.cell file
    with open(os.path.join(outputfolder, FrameworkName + '-1.cell'), 'r') as f:
        lines = f.read().splitlines()[1:]

    cellMatrix = np.array([np.array(line.split()[2:-1]).astype(float).reshape(3, 3) for line in lines])
    # Use ase library to convert the cell matrix to cell parameters
    cellParameters = np.array([Cell(i).cellpar() for i in cellMatrix])

    return cellMatrix.tolist(), cellParameters.tolist()


def getStructures(outputfolder: str, FrameworkName: str) -> list[float]:
    """
    Get the structures from CP2K optimization.

    Parameters
    ----------
    outputfolder : str
        Path to the output folder.
    FrameworkName : str
        Name of the framework.

    Returns
    -------
    structure_list : list
        List of the structures.
    """

    if not os.path.isfile(os.path.join(outputfolder, FrameworkName + '-pos-1.xyz')):
        return []

    # Open the FrameworkName-pos-1 file
    with open(os.path.join(outputfolder, FrameworkName + '-pos-1.xyz'), 'r') as f:
        lines = f.read().splitlines()

    n_atoms = int(lines[0])

    # Reshape lines to have the shape (n_atoms + 2, -1)
    lines = [i[2:] for i in np.array(lines).reshape((-1, n_atoms + 2))]

    atomLabelList = []
    atomPosList = []

    for structure in lines:
        atomLabelList.append([i.split()[0] for i in structure])
        atomPosList.append(np.array([np.array(i.split()[1:]).astype(float) for i in structure]))

    return atomLabelList, atomPosList


def getEnergies(outputfolder: str, FrameworkName: str) -> list[float]:
    """
    Get the forces from CP2K optimization.

    Parameters
    ----------
    outputfolder : str
        Path to the output folder.
    FrameworkName : str
        Name of the framework.

    Returns
    -------
    energyList : list
        List of the energies.
    """

    if not os.path.isfile(os.path.join(outputfolder, FrameworkName + '-pos-1.xyz')):
        return []

    # Open the FrameworkName-pos-1 file
    with open(os.path.join(outputfolder, FrameworkName + '-pos-1.xyz'), 'r') as f:
        lines = f.read().splitlines()

    energyList = np.array([line.split()[-1] for line in lines if 'E = ' in line]).astype(float)

    return energyList.tolist()


def getForces(outputfolder: str, FrameworkName: str) -> list:
    """
    Get the forces from CP2K optimization.

    Parameters
    ----------
    outputfolder : str
        Path to the output folder.
    FrameworkName : str
        Name of the framework.

    Returns
    -------
    force_list : list
        List of the forces.
    """

    if not os.path.isfile(os.path.join(outputfolder, FrameworkName + '-frc-1.xyz')):
        return []

    # Open the FrameworkName-pos-1 file
    with open(os.path.join(outputfolder, FrameworkName + '-frc-1.xyz'), 'r') as f:
        lines = f.read().splitlines()

    n_atoms = int(lines[0])

    # Reshape lines to have the shape (n_atoms + 2, -1)
    lines = [i[2:] for i in np.array(lines).reshape((-1, n_atoms + 2))]

    force_list = []

    for structure in lines:
        force_list.append(np.array([np.array(i.split()[1:]).astype(float) for i in structure]))

    return force_list


def getStress(outputfolder: str, FrameworkName: str) -> list:
    """
    Get the stress from CP2K optimization.

    Parameters
    ----------
    outputfolder : str
        Path to the output folder.
    FrameworkName : str
        Name of the framework.

    Returns
    -------
    stressList : list
        List of the stress tensors.
    """

    # Check if the FrameworkName-1.stress file exists
    if not os.path.isfile(os.path.join(outputfolder, FrameworkName + '-1.stress')):
        return []
    # Open the FrameworkName-1.stress file
    with open(os.path.join(outputfolder, FrameworkName + '-1.stress'), 'r') as f:
        lines = f.read().splitlines()[1:]

    stressList = np.array([np.array(line.split()[2:]).astype(float).reshape(3, 3) for line in lines])

    return stressList.tolist()


def get_MoldenData(OutputFolder, Frameworkname):
    """
    Get the vibrational information from the molden file.

    Parameters
    ----------
    Frameworkname : str
        Name of the framework.

    Returns
    -------
    eigenVectors : list
        List of the vibrational vectors.
    modes : list
        List of the vibrational modes.
    atom_labels : list
        List of the atomic labels.
    atom_pos : list
        List of the atomic positions.
    """
    # Read the molden file
    with open(os.path.join(OutputFolder, f'{Frameworkname}-VIBRATIONS-1.mol')) as f:
        molden_file = f.read().splitlines()

    position = {' [FREQ]': None,
                ' [FR-COORD]': None,
                ' [FR-NORM-COORD]': None,
                ' [INT]': None
                }

    for i, line in enumerate(molden_file):
        if line in position.keys():
            position[line] = i

    atom_list = molden_file[position[' [FR-COORD]'] + 1: position[' [FR-NORM-COORD]']]
    atom_labels = [atom.split()[0] for atom in atom_list]
    atom_pos = np.array([atom.split()[1:] for atom in atom_list]).astype(float)

    freq_list = molden_file[position[' [FREQ]'] + 1: position[' [FR-COORD]']]

    # convert atom_pos from bohr to angstrom
    atom_pos *= 0.529177

    # Create the string combinin the atom labels and positions
    atom_list = []

    for i, atom in enumerate(atom_pos):
        atom_list.append(f"{atom_labels[i]:3}     {atom[0]:15.9f}   {atom[1]:15.9f} {atom[2]:15.9f} ")

    vibrations = molden_file[position[' [FR-NORM-COORD]'] + 1: position[' [INT]']]

    # Reshape vibrations to the shape of (len(atom_list) + 1, -1)
    vibrations = [vibrations[i + 1:i + len(atom_list) + 1] for i in range(0, len(vibrations), len(atom_list) + 1)]

    eigenVectors = [[i.split() for i in mode] for mode in vibrations]
    eigenVectors = [[float(i) for sublist in mode for i in sublist] for mode in eigenVectors]

    intensity = [float(i) for i in molden_file[position[' [INT]'] + 1:]]

    modes = [i + 1 for i in range(len(intensity))]

    return atom_labels, atom_pos, vibrations, modes, eigenVectors, intensity, freq_list


def get_vibrational_data(CP2K_output_name) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get the vibrational data from the CP2K output file.

    Parameters
    ----------
    CP2K_output_name : str
        Path to the CP2K output file.

    Returns
    -------
    frequency : np.ndarray
        Array of the vibrational frequencies in cm^-1
    IR_intensity : np.ndarray
        Array of the IR intensities in KM/Mole
    RAMAN_intensity : np.ndarray
        Array of the RAMAN intensities in A^4/AMU
    """
    output_file = open(CP2K_output_name, 'r').read().splitlines()

    # Find the line with the text: "NORMAL MODES - CARTESIAN DISPLACEMENTS"
    normal_modes = None

    for i, line in enumerate(output_file):
        if 'NORMAL MODES - CARTESIAN DISPLACEMENTS' in line:
            normal_modes = output_file[i:]

    frequency = np.array([])
    IR_intensity = np.array([])
    RAMAN_intensity = np.array([])

    for line in normal_modes:
        if ' VIB|Frequency (cm^-1)' in line:
            freq = np.array([n.replace('*', '0') for n in line.split()[2:]]).astype(float)
            frequency = np.append(frequency, freq)
        if 'VIB|IR int (KM/Mole)' in line:
            ir_int = np.array([n.replace('*', '0') for n in line.split()[3:]]).astype(float)
            IR_intensity = np.append(IR_intensity, ir_int)
        if ' VIB|Raman (A^4/amu)' in line:
            raman_int = np.array([n.replace('*', '0') for n in line.split()[2:]]).astype(float)
            RAMAN_intensity = np.append(RAMAN_intensity, raman_int)

    return frequency, IR_intensity, RAMAN_intensity