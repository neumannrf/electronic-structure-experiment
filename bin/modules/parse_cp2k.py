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

    structure_list = []

    for structure in lines:
        atom_labels = [i.split()[0] for i in structure]
        atom_pos = np.array([np.array(i.split()[1:]).astype(float) for i in structure]).T
        structure_list.append([atom_labels, atom_pos.tolist()])

    return structure_list


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
        atom_labels = [i.split()[0] for i in structure]
        atom_pos = np.array([np.array(i.split()[1:]).astype(float) for i in structure]).T
        force_list.append([atom_labels, atom_pos.tolist()])

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
