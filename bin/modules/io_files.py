#!/usr/bin/env -S python -B

# SPDX-License-Identifier: Apache2.0

import os
import json
import gemmi

import numpy as np
from textwrap import dedent

from ase.cell import Cell

from modules.atom_data import ATOMIC_NUMBER
from modules.calculate_properties import (get_MoldenData,
                                          get_vibrational_data,
                                          get_CellParameters)


def readChemicalJSON(FrameworkName: str, OutputFolder: str = '.', **kwargs):
    """
    Read the chemical JSON file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    OutputFolder : str
        Path to the output folder. Default: `.`

    A chargeType can be passed as kwargs to read the charges from the cjson file.
    """

    # Read the cif file and get the lattice parameters and atomic positions
    cjson_filename = os.path.join(OutputFolder, FrameworkName + '.cjson')

    with open(cjson_filename, 'r') as f:
        ChemJSON = json.load(f)

    # Get the cell parameters from cif file
    if 'a' in ChemJSON['unitCell']:
        CellParameters = [ChemJSON['unitCell']['a'],
                          ChemJSON['unitCell']['b'],
                          ChemJSON['unitCell']['c'],
                          ChemJSON['unitCell']['alpha'],
                          ChemJSON['unitCell']['beta'],
                          ChemJSON['unitCell']['gamma']]

    elif 'cellVectors' in ChemJSON['unitCell']:
        CellMatrix = ChemJSON['unitCell']['cellVectors']
        aseCell = Cell(CellMatrix)
        CellParameters = aseCell.cellpar()

    else:
        print('Could not find the cell parameters in the cjson file.')
        CellParameters = None

    # Get the atomic labels
    if 'type' in ChemJSON['atoms']['elements']:
        labels = ChemJSON['atoms']['elements']['type']
    elif 'number' in ChemJSON['atoms']['elements']:
        labels = [gemmi.Element(i).name for i in ChemJSON['atoms']['elements']['number']]

    # Get the fractional coordinates
    if '3dFractional' in ChemJSON['atoms']['coords']:
        frac_x = ChemJSON['atoms']['coords']['3dFractional'][0::3]
        frac_y = ChemJSON['atoms']['coords']['3dFractional'][1::3]
        frac_z = ChemJSON['atoms']['coords']['3dFractional'][2::3]

    elif '3d' in ChemJSON['atoms']['coords']:
        cart_x = ChemJSON['atoms']['coords']['3dFractional'][0::3]
        cart_y = ChemJSON['atoms']['coords']['3dFractional'][1::3]
        cart_z = ChemJSON['atoms']['coords']['3dFractional'][2::3]

        # Convert fractional coordinates to cartesian
        cartPositions = np.array([cart_x, cart_y, cart_z]).T
        frac_x, frac_y, frac_z = aseCell.scaled_positions(cartPositions).T

    # Get the charges
    if 'partialCharges' in ChemJSON:
        # Check if some charge type was passed as kwargs
        if 'charge_type' in kwargs:
            chargeType = kwargs['charge_type']
        else:
            chargeType = list(ChemJSON['partialCharges'].keys())[0]
        charges = ChemJSON['partialCharges'][chargeType]
    else:
        charges = None
        chargeType = None

    return CellParameters, labels, frac_x, frac_y, frac_z, charges, chargeType


def saveChemicalJSON(FrameworkName: str,
                     CellParameters: float,
                     labels: list[str],
                     frac_x: list[float],
                     frac_y: list[float],
                     frac_z: list[float],
                     charges: list[float] = None,
                     charge_type: str = 'none',
                     OutputFolder: str = '.'):
    """
    Save the chemical JSON file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    CellParameters : list
        List of the cell parameters.
    labels : list
        List of the atomic labels.
    frac_x : list
        List of the atomic positions along the `a` vector.
    frac_y : list
        List of the atomic positions along the `b` vector.
    frac_z : list
        List of the atomic positions along the `c` vector.
    charges : list
        List of the atomic charges.
    charge_type : str
        Type of the charges.
    OutputFolder : str
        Path to the output folder. Default: `.`
    """

    # Convert atom_labels to atom_number
    atom_number = [gemmi.Element(atom).atomic_number for atom in labels]

    # Flatten the atom_pos list
    atom_pos = np.array([frac_x, frac_y, frac_z]).T

    aseCell = Cell.fromcellpar(CellParameters)

    # Get the cell matrix
    CellMatrix = aseCell.tolist()

    FracCoords = atom_pos.flatten().tolist()
    CartCoords = aseCell.cartesian_positions(atom_pos).flatten().tolist()

    formula = ' '.join([f'{atom}{labels.count(atom)}' for atom in set(labels)])

    ChemJSON = {
        "chemicalJson": 1,
        "name": FrameworkName,
        "formula": formula,
        "unitCell": {
            "a": CellParameters[0],
            "b": CellParameters[1],
            "c": CellParameters[2],
            "alpha": CellParameters[3],
            "beta":  CellParameters[4],
            "gamma": CellParameters[5],
            "cellVectors": CellMatrix
        },
        "atoms": {
            "elements": {
                "type": labels,
                "number": atom_number
                },
            "coords": {
                "3d": CartCoords,
                "3dFractional": FracCoords
                }
        }
    }

    if charges is not None:
        ChemJSON["partialCharges"] = {
            charge_type: charges.tolist()
        }

    # Save ChemJSON as a json file
    with open(os.path.join(OutputFolder, f'{FrameworkName}.cjson'), 'w') as f:
        json.dump(ChemJSON, f, indent=4)


def readXSF(FrameworkName: str, OutputFolder: str = '.', **kwargs):
    """
    Read the XSF file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    OutputFolder : str
        Path to the output folder. Default: `.`
    """

    with open(os.path.join(OutputFolder, FrameworkName + '.xsf'), 'r') as f:
        lines = f.read().splitlines()

    # Get the cell parameters
    cellMatrix = np.array([lines[2].split(), lines[3].split(), lines[4].split()]).astype(float)

    aseCell = Cell(cellMatrix)

    cellParameters = aseCell.cellpar()

    # Get the atomic labels and positions
    n_atoms = int(lines[6].split()[0])
    atoms = lines[7:7 + n_atoms]

    atomLabels = [i.split()[0] for i in atoms]

    cart_x = [float(i.split()[1]) for i in atoms]
    cart_y = [float(i.split()[2]) for i in atoms]
    cart_z = [float(i.split()[3]) for i in atoms]

    cartPositions = np.array([cart_x, cart_y, cart_z]).T

    # Convert fractional coordinates to cartesian
    frac_x, frac_y, frac_z = aseCell.scaled_positions(cartPositions).T

    return cellParameters, atomLabels, frac_x, frac_y, frac_z, None, None


def saveXSF(FrameworkName: str,
            CellParameters: float,
            labels: list[str],
            frac_x: list[float],
            frac_y: list[float],
            frac_z: list[float],
            OutputFolder: str = '.',
            **kwargs):
    '''
    Save the XSF file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    CellParameters : list
        List of the cell parameters.
    labels : list
        List of the atomic labels.
    frac_x : list
        List of the atomic positions along the `a` vector.
    frac_y : list
        List of the atomic positions along the `b` vector.
    frac_z : list
        List of the atomic positions along the `c` vector.
    OutputFolder : str
        Path to the output folder. Default: `.`
    '''

    aseCell = Cell.fromcellpar(CellParameters)
    # Get the cell parameters from cif file
    CellMatrix = aseCell.tolist()

    # Convert fractional coordinates to cartesian
    cart_coords = aseCell.cartesian_positions(np.array([frac_x, frac_y, frac_z]).T)

    xsf_file = dedent(f'''CRYSTAL
PRIMVEC
  {CellMatrix[0][0]:15.10f}    {CellMatrix[0][1]:15.10f}    {CellMatrix[0][2]:15.10f}
  {CellMatrix[1][0]:15.10f}    {CellMatrix[1][1]:15.10f}    {CellMatrix[1][2]:15.10f}
  {CellMatrix[2][0]:15.10f}    {CellMatrix[2][1]:15.10f}    {CellMatrix[2][2]:15.10f}
PRIMCOORD    1
      {len(labels)}   1
''')
    for j, atom in enumerate(labels):
        xsf_file += ' {} {:15.10f}  {:15.10f}  {:15.10f}\n'.format(atom,
                                                                   cart_coords[j][0],
                                                                   cart_coords[j][1],
                                                                   cart_coords[j][2])

    with open(os.path.join(OutputFolder, f'{FrameworkName}.xsf'), 'w') as f:
        f.write(xsf_file)


def readCIF(FrameworkName: str,
            OutputFolder: str = '.',
            **kwargs):
    """
    Read the CIF file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    OutputFolder : str
        Path to the output folder. Default: `.`
    """

    # Read the cif file and get the lattice parameters and atomic positions
    cif_filename = os.path.join(OutputFolder, FrameworkName + '.cif')

    cif = gemmi.cif.read_file(cif_filename).sole_block()

    a = float(cif.find_value('_cell_length_a').split('(')[0])
    b = float(cif.find_value('_cell_length_b').split('(')[0])
    c = float(cif.find_value('_cell_length_c').split('(')[0])
    beta = float(cif.find_value('_cell_angle_beta').split('(')[0])
    gamma = float(cif.find_value('_cell_angle_gamma').split('(')[0])
    alpha = float(cif.find_value('_cell_angle_alpha').split('(')[0])

    CellParameters = [a, b, c, alpha, beta, gamma]

    AtomicTypes = list(cif.find_values('_atom_site_type_symbol'))
    PosX = np.array(cif.find_values('_atom_site_fract_x')).astype(float)
    PosY = np.array(cif.find_values('_atom_site_fract_y')).astype(float)
    PosZ = np.array(cif.find_values('_atom_site_fract_z')).astype(float)
    try:
        charges = np.array(cif.find_values('_atom_site_charge')).astype(float)
        charge_type = 'DDEC'
    except Exception:
        charges = None
        charge_type = None

    return CellParameters, AtomicTypes, PosX, PosY, PosZ, charges, charge_type


def saveCIF(FrameworkName: str,
            CellParameters: float,
            labels: list[str],
            frac_x: list[float],
            frac_y: list[float],
            frac_z: list[float],
            charges: list[float] = None,
            OutputFolder: str = '.',
            **kwargs) -> None:
    """
    Save the CIF file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    cell : list
        List of the cell parameters.
    labels : list
        List of the atomic labels.
    frac_x : list
        List of the atomic positions along the `a` vector.
    frac_y : list
        List of the atomic positions along the `b` vector.
    frac_z : list
        List of the atomic positions along the `c` vector.
    charges : list
        List of the atomic charges.
    OutputFolder : str
        Path to the output folder. Default: `.`
    """
    cif_file = dedent(f"""\
data_{FrameworkName}
_chemical_name_common                  '{FrameworkName}'
_cell_length_a                          {CellParameters[0]:10.5f}
_cell_length_b                          {CellParameters[1]:10.5f}
_cell_length_c                          {CellParameters[2]:10.5f}
_cell_angle_alpha                       {CellParameters[3]:10.5f}
_cell_angle_beta                        {CellParameters[4]:10.5f}
_cell_angle_gamma                       {CellParameters[5]:10.5f}

_symmetry_cell_setting          triclinic
_symmetry_space_group_name_Hall 'P 1'
_symmetry_space_group_name_H-M  'P 1'
_symmetry_Int_Tables_number     1

_symmetry_equiv_pos_as_xyz 'x,y,z'

loop_
   _atom_site_type_symbol
   _atom_site_label
   _atom_site_fract_x
   _atom_site_fract_y
   _atom_site_fract_z
""")

    if charges is not None:
        cif_file += '   _atom_site_charge\n'

    # Unique atoms in labes
    symbols_number = {i: 0 for i in set(labels)}
    symbols = []
    for label in labels:
        symbols_number[label] += 1
        symbols.append(label + str(symbols_number[label]))

    if charges is None:
        for label, symbol, x, y, z in zip(labels, symbols, frac_x, frac_y, frac_z):
            cif_file += "   {:3s}   {:6s}   {:15.10f}   {:15.10f}   {:15.10f}\n".format(
                label, symbol, x, y, z
                )
    else:
        for label, symbol, x, y, z, charge in zip(labels, symbols, frac_x, frac_y, frac_z, charges):
            cif_file += "   {:3s}   {:6s}   {:15.10f}   {:15.10f}   {:15.10f}   {:10.7f}\n".format(
                label, symbol, x, y, z, charge
                )

    with open(os.path.join(OutputFolder, f'{FrameworkName}.cif'), 'w') as f:
        f.write(cif_file)


def readGJF(FrameworkName: str, OutputFolder: str, **kwargs):
    """
    Read a GJF file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    OutputFolder : str
        Path to the output folder. Default: `.`
    """

    with open(os.path.join(OutputFolder, FrameworkName + '.gjf'), 'r') as f:
        lines = f.read().splitlines()

    # Remove empty lines
    lines = [line for line in lines if len(line.split()) > 0]

    # Get the cell parameters
    cellMatrix = []

    for line in lines:
        if 'Tv' in line:
            cellMatrix.append(line.split()[1:])

    cellMatrix = np.array(cellMatrix).astype(float)

    atomLabels = []
    cartPositions = []
    # Get the atomic labels
    for line in lines:
        if line.split()[0] in ATOMIC_NUMBER.keys():
            atomLabels.append(line.split()[0])
            cartPositions.append(line.split()[1:])

    cartPositions = np.array(cartPositions).astype(float)

    # Convert cartesian coordinates to fractional
    aseCell = Cell(cellMatrix)
    frac_x, frac_y, frac_z = aseCell.scaled_positions(cartPositions).T

    # Get the cell matrix
    cellParameters = aseCell.cellpar()

    return cellParameters, atomLabels, frac_x, frac_y, frac_z, None, None


def saveGJF(FrameworkName: str,
            CellParameters: float,
            labels: list[str],
            frac_x: list[float],
            frac_y: list[float],
            frac_z: list[float],
            OutputFolder: str = '.',
            **kwargs):
    """
    Save the GJF file.

    Parameters
    ----------
    FrameworkName : str
        Name of the framework.
    cell : list
        List of the cell parameters.
    labels : list
        List of the atomic labels.
    frac_x : list
        List of the atomic positions along the `a` vector.
    frac_y : list
        List of the atomic positions along the `b` vector.
    frac_z : list
        List of the atomic positions along the `c` vector.
    OutputFolder : str
        Path to the output folder. Default: `.`
    """

    aseCell = Cell.fromcellpar(CellParameters)
    cellMatrix = aseCell.tolist()

    # Get the cartesian coordinates
    carCoords = aseCell.cartesian_positions(np.array([frac_x, frac_y, frac_z]).T)

    gjf_file = dedent(f"""%chk={FrameworkName}.chk
# pbepbe/3-21g/auto

Title Card Required

0 1
""")
    for i, atom in enumerate(labels):
        gjf_file += f' {atom:3} {carCoords[i][0]:15.10f} {carCoords[i][1]:15.10f} {carCoords[i][2]:15.10f}\n'

    for i in range(3):
        gjf_file += f' Tv {cellMatrix[i][0]:15.10f} {cellMatrix[i][1]:15.10f} {cellMatrix[i][2]:15.10f}\n'

    with open(os.path.join(OutputFolder, f'{FrameworkName}.gjf'), 'w') as f:
        f.write(gjf_file)


def saveVibrationalVectors(OutputFolder, Frameworkname):
    '''
    Get the vibrational vectors from the CP2K output file.
    '''

    atom_labels, _, vibrations, _, _, _, freq_list = get_MoldenData(OutputFolder, Frameworkname)

    # Get the cell parameters from cif file
    CellParameters = get_CellParameters(Frameworkname + '.cif')
    CellMatrix = Cell.fromcellpar(CellParameters).tolist()

    os.makedirs(os.path.join(OutputFolder, 'VIBRATION_FILES'), exist_ok=True)

    for i, vib in enumerate(vibrations):

        axsf_filename = os.path.join(os.path.join(OutputFolder, 'VIBRATION_FILES'),
                                     f'VIBRATIONS-{i}-{float(freq_list[i])}.axsf')
        axsf_file = open(axsf_filename, 'w')

        axsf_file.write('CRYSTAL\n')
        axsf_file.write('PRIMVEC  \n')
        axsf_file.write(f'  {CellMatrix[0][0]:15.10f}    {CellMatrix[0][1]:15.10f}    {CellMatrix[0][2]:15.10f}\n')
        axsf_file.write(f'  {CellMatrix[1][0]:15.10f}    {CellMatrix[1][1]:15.10f}    {CellMatrix[1][2]:15.10f}\n')
        axsf_file.write(f'  {CellMatrix[2][0]:15.10f}    {CellMatrix[2][1]:15.10f}    {CellMatrix[2][2]:15.10f}\n')
        axsf_file.write('PRIMCOORD    1\n')
        axsf_file.write(f'      {len(atom_labels)}   1\n')
        for j, atom in enumerate(atom_labels):
            axsf_file.write(f' {atom} {vib[j]}\n')

        axsf_file.close()

    return None


def saveVibrationalChemicalJSON(OutputFolder, Frameworkname):

    frequency, IR_intensity, RAMAN_intensity = get_vibrational_data(
        os.path.join(OutputFolder, 'simulation_Vibrations.out')
        )

    atom_labels, atom_pos, _, modes, eigenVectors, _, _ = get_MoldenData(OutputFolder,
                                                                         Frameworkname)

    # Convert atom_labels to atom_number
    atom_number = [gemmi.Element(atom).atomic_number for atom in atom_labels]

    # Flatten the atom_pos list
    atom_pos = [i for sublist in atom_pos for i in sublist]

    # Get the cell parameters from cif file
    CellParameters = get_CellParameters(Frameworkname + '.cif')
    CellMatrix = Cell.fromcellpar(CellParameters).flatten().tolist()

    formula = ' '.join([f'{atom}{atom_labels.count(atom)}' for atom in set(atom_labels)])

    ChemJSON = {
        "chemicalJson": 1,
        "name": Frameworkname,
        "formula": formula,
        "unitCell": {
            "a": CellParameters[0],
            "b": CellParameters[1],
            "c": CellParameters[2],
            "alpha": CellParameters[3],
            "beta":  CellParameters[4],
            "gamma": CellParameters[5],
            "cellVectors": CellMatrix
        },
        "atoms": {
            "elements": {
                "type": atom_labels,
                "number": atom_number
                },
            "coords": {
                "3d": atom_pos
                }
        },
        'vibrations': {
            'eigenVectors': eigenVectors,
            'frequencies': list(frequency),
            'intensities': list(IR_intensity),
            'ramanIntensities': list(RAMAN_intensity),
            'modes': modes
        }
        }

    # Save ChemJSON as a json file
    with open(os.path.join(OutputFolder, f'{Frameworkname}.cjson'), 'w') as f:
        json.dump(ChemJSON, f, indent=4)
