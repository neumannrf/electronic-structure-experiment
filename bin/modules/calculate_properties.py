#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import gemmi

import numpy as np

from ase.cell import Cell


def calculate_Perpendicular_Widths(cif_filename: str) -> tuple[float, float, float]:
    """
    Calculate the perpendicular widths of the unit cell.
    RASPA considers the perpendicular directions as the directions perpendicular to the `ab`,
    `bc`, and `ca` planes. Thus, the directions depend on the crystallographic vectors `a`, `b`,
    and `c`.
    The length in the perpendicular directions are the projections of the crystallographic vectors
    on the vectors `a x b`, `b x c`, and `c x a`. (here `x` means cross product)
    """
    # Read data from CIF file
    cif = gemmi.cif.read_file(cif_filename).sole_block()
    a = float(cif.find_value('_cell_length_a').split('(')[0])
    b = float(cif.find_value('_cell_length_b').split('(')[0])
    c = float(cif.find_value('_cell_length_c').split('(')[0])
    beta = float(cif.find_value('_cell_angle_beta').split('(')[0]) * np.pi / 180.0
    gamma = float(cif.find_value('_cell_angle_gamma').split('(')[0]) * np.pi / 180.0
    alpha = float(cif.find_value('_cell_angle_alpha').split('(')[0]) * np.pi / 180.0

    # Calculate the nu value
    nu = (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)

    # Build the transformation matrix as a numpy array
    CellBox = np.array([[a, 0.0, 0.0],
                        [b * np.cos(gamma), b * np.sin(gamma), 0.0],
                        [c * np.cos(beta), c * nu, c * np.sqrt(1.0 - np.cos(beta)**2 - nu**2)]])

    # Calculate the cross products
    axb = np.cross(CellBox[0], CellBox[1])
    bxc = np.cross(CellBox[1], CellBox[2])
    cxa = np.cross(CellBox[2], CellBox[0])

    # Calculates the volume of the unit cell
    V = np.dot(np.cross(CellBox[0], CellBox[1]), CellBox[2])

    # Calculate perpendicular widths
    p_width_1 = V / np.linalg.norm(bxc)
    p_width_2 = V / np.linalg.norm(cxa)
    p_width_3 = V / np.linalg.norm(axb)

    return p_width_1, p_width_2, p_width_3


def calculate_UnitCells(cif_filename: str, cutoff: float) -> str:
    """
    Calculate the number of unit cell repetitions so that all supercell lengths are larger than
    twice the interaction potential cut-off radius.
    """

    # Calculate the perpendicular widths
    p_width_1, p_width_2, p_width_3 = calculate_Perpendicular_Widths(cif_filename)

    # Calculate UnitCells string
    uc_array = np.ceil(2.0 * cutoff / np.array([p_width_1, p_width_2, p_width_3])).astype(int)
    unit_cells = ' '.join(map(str, uc_array))

    return unit_cells


def get_CellParameters(cif_filename: str) -> tuple[float, float, float, float, float, float]:
    """
    Calculate the perpendicular widths of the unit cell.
    RASPA considers the perpendicular directions as the directions perpendicular to the `ab`,
    `bc`, and `ca` planes. Thus, the directions depend on the crystallographic vectors `a`, `b`,
    and `c`.
    The length in the perpendicular directions are the projections of the crystallographic vectors
    on the vectors `a x b`, `b x c`, and `c x a`. (here `x` means cross product)

    Parameters
    ----------
    cif_filename : str
        Path to the CIF file.

    Returns
    -------
    a : float
        Length of the `a` vector.
    b : float
        Length of the `b` vector.
    c : float
        Length of the `c` vector.
    alpha : float
        Angle between the `b` and `c` vectors.
    beta : float
        Angle between the `a` and `c` vectors.
    gamma : float
        Angle between the `a` and `b` vectors.
    """
    # Read data from CIF file
    cif = gemmi.cif.read_file(cif_filename).sole_block()
    a = float(cif.find_value('_cell_length_a').split('(')[0])
    b = float(cif.find_value('_cell_length_b').split('(')[0])
    c = float(cif.find_value('_cell_length_c').split('(')[0])
    beta = float(cif.find_value('_cell_angle_beta').split('(')[0])
    gamma = float(cif.find_value('_cell_angle_gamma').split('(')[0])
    alpha = float(cif.find_value('_cell_angle_alpha').split('(')[0])

    return a, b, c, alpha, beta, gamma


def get_AtomicPositions(cif_filename: str) -> tuple[list[str], list[float], list[float], list[float]]:
    """
    Get the atomic positions of the unit cell.

    Parameters
    ----------
    cif_filename : str
        Path to the CIF file.

    Returns
    -------
    atom_site_type_symbol : list
        List of the atomic symbols.
    atom_site_fract_x : list
        List of the fractional coordinates of the atoms along the `a` vector.
    atom_site_fract_y : list
        List of the fractional coordinates of the atoms along the `b` vector.
    atom_site_fract_z : list
        List of the fractional coordinates of the atoms along the `c` vector.
    """
    # Read data from CIF file
    cif = gemmi.cif.read_file(cif_filename).sole_block()
    atom_site_type_symbol = cif.find_values('_atom_site_type_symbol')
    atom_site_fract_x = np.array(cif.find_values('_atom_site_fract_x')).astype(float)
    atom_site_fract_y = np.array(cif.find_values('_atom_site_fract_y')).astype(float)
    atom_site_fract_z = np.array(cif.find_values('_atom_site_fract_z')).astype(float)

    return atom_site_type_symbol, atom_site_fract_x, atom_site_fract_y, atom_site_fract_z


def get_DDECAtomicCharges(xyz_filename: str) -> list[float]:
    """
    Get the atomic charges from the xyz file.

    Parameters
    ----------
    xyz_filename : str
        Path to the xyz file.

    Returns
    -------
    charges : list
        List of the atomic charges.
    """
    # Read the xyz file
    with open(xyz_filename, 'r') as f:
        # Read the lines
        lines = f.readlines()

        # Get the number of atoms
        n_atoms = int(lines[0].split()[0])

        # Get the atomic charges
        charges = []
        for line in lines[2:2 + n_atoms]:
            charges.append(float(line.split()[-1]))

    return charges


def get_CM5AtomicCharges(output_filename: str) -> list[float]:
    """
    Get the atomic charges from the xyz file.

    Parameters
    ----------
    xyz_filename : str
        Path to the xyz file.

    Returns
    -------
    charges : list
        List of the atomic charges.
    """
    # Read the xyz file
    with open(output_filename, 'r') as f:
        # Read the lines
        lines = f.readlines()

        start_pos = [i for i in range(len(lines)) if 'The computed CM5 net atomic charges are:' in lines[i]][0]
        end_pos = [i for i in range(len(lines)) if 'Hirshfeld and CM5 analysis finished' in lines[i]][0]

        charge_lines = ' '.join(lines[start_pos + 1: end_pos]).split()

        charges = [float(i) for i in charge_lines]

    return charges


def get_vibrational_data(CP2K_output_name) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
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
    '''
    output_file = open(CP2K_output_name, 'r').read().splitlines()

    # Find the line with the text: "NORMAL MODES - CARTESIAN DISPLACEMENTS"
    normal_modes = None

    for i, line in enumerate(output_file):
        if 'NORMAL MODES - CARTESIAN DISPLACEMENTS' in line:
            normal_modes = output_file[i:]

    frequency = np.array([])
    IR_intensity = np.array([])
    RAMAN_intensity = np.array([])

    for i, line in enumerate(normal_modes):
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


def get_MoldenData(OutputFolder, Frameworkname):
    '''
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

    '''
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


def lorentzian(x, x0, gamma) -> np.ndarray[float]:
    return 1/np.pi * gamma / ((x-x0)**2 + gamma**2)


def getCellParametersFromOptimization(outputfolder, FrameworkName) -> list[np.ndarray[float]]:

    # Open the FrameworkName-1.cell file
    with open(os.path.join(outputfolder, FrameworkName + '-1.cell'), 'r') as f:
        lines = f.read().splitlines()[1:]

    cellList = [np.array(line.split()[2:-1]).astype(float).reshape(3, 3) for line in lines]
    # Use ase library to convert the cell matrix to cell parameters
    cellParameters = [Cell(i).cellpar() for i in cellList]

    return cellParameters


def getStructuresFromOptimization(outputfolder, FrameworkName) -> list:

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
        structure_list.append([atom_labels, atom_pos])

    return structure_list
