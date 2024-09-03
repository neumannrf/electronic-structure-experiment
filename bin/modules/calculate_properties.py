#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
from types import SimpleNamespace

import gemmi
import numpy as np
from ase.cell import Cell
from cp2k_input_tools.generator import CP2KInputGenerator
from modules.atom_data import BASIS_SET, PSEUDO_POTENTIALS
from phonopy.harmonic.force_constants import similarity_transformation
from modules.constants import (c, k_B, h)


def calculate_Perpendicular_Widths(cif_filename: str) -> tuple[float, float, float]:
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
    p_width_1 : float
        Perpendicular width in the direction perpendicular to the `ab` plane.
    p_width_2 : float
        Perpendicular width in the direction perpendicular to the `bc` plane.
    p_width_3 : float
        Perpendicular width in the direction perpendicular to the `ca` plane.
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

    Parameters
    ----------
    cif_filename : str
        Path to the CIF file.
    cutoff : float
        Interaction potential cut-off radius in angstrom.

    Returns
    -------
    unit_cells : str
        String containing the number of unit cell repetitions in the `a`, `b`, and `c` directions.
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


def getForcesFromOptimization(outputfolder, FrameworkName) -> list:
    """
    Get the forces from the optimization output files.

    Parameters
    ----------
    outputfolder : str
        Path to the output folder.
    FrameworkName : str
        Name of the framework.

    Returns
    -------
    forces_list : list
        List of the forces.
    """

    # List the files in the output folder with name {FrameworkName}-forces-1_1.xyz
    files = [i for i in os.listdir(outputfolder) if f'{FrameworkName}-forces-1' in i]

    forces_list = []

    for i in range(len(files)):
        forces_list.append(get_forces(FileName=f'{FrameworkName}-forces-1_{i + 1}',
                                      output_folder=outputfolder))

    return forces_list

def get_spg_class(spgnum) -> str:
    """
    Get the space group class from the space group number.
    Parameters
    ----------
    spgnum : int
        Space group number

    Returns
    -------
    spgclass : str
        Space group class
    """
    # Triclinic or monoclinic or orthorhombic
    spgclass = ""
    if (spgnum < 3):
        spgclass = "triclinic"
    elif (spgnum > 3 and spgnum < 16):
        spgclass = "monoclinic"
    elif (spgnum > 15 and spgnum < 75):
        spgclass = "orthorhombic"
    elif (spgnum > 74 and spgnum < 143):
        spgclass = "tetragonal"
    elif (spgnum > 142 and spgnum < 168):
        spgclass = "trigonal"
    elif (spgnum > 167 and spgnum < 195):
        spgclass = "hexagonal"
    elif (spgnum > 194):
        spgclass = "cubic"

    return spgclass


def get_reciprocal_vectors(CellMatrix) -> tuple[float, float, float]:
    """
    Get the reciprocal vectors of a cell given in cell parameters of cell vectors
    ----------
    CellMatrix : array
        (3,1) array for cell vectors
    Returns
    -------
    b1 : array
        (3,1) array containing b_1 vector in the reciprocal space
    b2 : array
        (3,1) array containing b_2 vector in the reciprocal space
    b3 : array
        (3,1) array containing b_3 vector in the reciprocal space
    """

    v1, v2, v3 = CellMatrix

    vol = np.dot(v1, np.cross(v2, v3))

    b1 = 2 * np.pi * np.cross(v2, v3) / vol
    b2 = 2 * np.pi * np.cross(v3, v1) / vol
    b3 = 2 * np.pi * np.cross(v1, v2) / vol

    return b1, b2, b3


def get_kgrid(cell, dist=0.3) -> tuple[float, float, float]:
    """Get the k-points grid in the reciprocal space with a given distance for a
    cell given in cell parameters of cell vectors.
    ----------
    cell : array
        (3,1) array for cell vectors or (6,1) array for cell parameters
    distance : float
        distance between the points in the reciprocal space
    Returns
    -------
    kx : int
        Number of points in the x direction on reciprocal space
    ky : int
        Number of points in the y direction on reciprocal space
    kz : int
        Number of points in the z direction on reciprocal space
    """

    b1, b2, b3 = get_reciprocal_vectors(cell)

    b = np.array([np.linalg.norm(b1),
                  np.linalg.norm(b2),
                  np.linalg.norm(b3)])

    kx = np.ceil(b[0]/dist).astype(int)
    ky = np.ceil(b[1]/dist).astype(int)
    kz = np.ceil(b[2]/dist).astype(int)

    return kx, ky, kz


def create_input_file(FrameworkName: str,
                      output_folder: str,
                      **kwargs):
    """
    Create the input file for CP2K

    Parameters
    ----------
    FrameworkName : str
        Name of the framework
    output_folder : str
        Path to the output folder
    **kwargs : dict
        Dictionary with the parameters to be used in the CP2K input file creation.
    """

    CalcDict = {
        'FrameworkName': FrameworkName.split('.')[0],
        'Charge': 0,
        'Multiplicity': 1,
        'CalcType': 'energy_force',  # Can be 'energy_force', 'cell_opt', 'geo_opt', 'md', or 'normal_modes'
        'UseOT': False,
        'UseSmearing': False,
        'SmearingMethod': 'fermi_dirac',  # Can be 'fermi_dirac' or 'energy_window'
        'ElectronicTemperature': 300,
        'WindowSize': 0.1,
        'AddedMOs': 0,
        'MixingMethod': 'broyden_mixing',  # Can be 'direct_p_mixing', 'broyden_mixing_new', or 'kerker_mixing'
        'MixingAlpha': 0.2,
        'MaxSCFycles': 30,
        'MaxOuterSCFycles': 10,
        'EPSDefault': 1e-8,
        'PWCutoff': 1200,
        'NGrid': 5,
        'RelativeCutOff': 60,
        'Functional': 'PBE',  # Can be 'PBE', 'XTB', or 'PBE0'
        'Parametrization': 'ORIG',  # Can be 'ORIG', 'PBESOL', or 'REVPBE'
        'DispersionCorrection': 'DFTD3',  # Can be None, 'DFTD2', 'DFTD3', or 'DFTD3(BJ)'
        'CheckAtomicCharges': True,
        'BasisSet': 'DZVP',  # Can be 'DZVP', 'TZVP', or 'TZV2P'
        'SCFGuess': 'atomic',  # Can be 'atomic', 'restart', 'core', 'random', 'sparse', or 'mopac'
        'SCFConvergence': 1e-8,
        'CP2KDataDir': os.environ.get("CP2K_DATA_DIR"),
        'KeepSymmetry': False,
        'KeepSpaceGroup': False,
        'KeepAngles': False,
        'MaxIterations': 100,
        'Restart': False,
        'MaxDR': 1e-3,
        'MaxForce': 1e-3,
        'RMSDR': 1e-3,
        'RMSForce': 1e-3,
        'UseScalapack': False,
        'CellParameters': None,
        'CellMatrix': None,
        'AtomicTypes': None,
        'FracX': None,
        'FracY': None,
        'FracZ': None,
        'CartX': None,
        'CartY': None,
        'CartZ': None,
        'ProcsPerReplica': 4,
        'dX': 0.001,
        'CalculateRaman': False,
        'CalculateIR': False,
        'Ensemble': 'NPT_F',
        'Temperature': 400,
        'TimeStep': 0.5,
        'MDSteps': 100,
        'Pressure': 1,
        'TimeCon': 1000,
        'KPoints': False,
        'RecDist': 0.3
    }

    # TO-DO: Add conversion from frac to cart and vice versa

    # Update the dictionary with the user input
    CalcDict.update(kwargs)

    calcPar = SimpleNamespace(**CalcDict)

    Coord_Dict = {
        'scaled': False,
        '*': ['{:3} {:11.6f} {:11.6f} {:11.6f}'.format(calcPar.AtomicTypes[i],
                                                       calcPar.CartX[i],
                                                       calcPar.CartY[i],
                                                       calcPar.CartZ[i]) for i in range(len(calcPar.AtomicTypes))]
                    }

    Kind_List = []

    for specie in set(calcPar.AtomicTypes):
        Kind_List.append(
            {
                "_": specie,
                'element': specie,
                'potential': PSEUDO_POTENTIALS[specie],
                'basis_set': BASIS_SET[calcPar.BasisSet][specie]
            }
        )

    if calcPar.CellParameters is not None:
        Cell_Dict = {
            'abc': [calcPar.CellParameters[0], calcPar.CellParameters[1], calcPar.CellParameters[2]],
            'alpha_beta_gamma': [calcPar.CellParameters[3], calcPar.CellParameters[4], calcPar.CellParameters[5]],
            'periodic': 'XYZ'
            }
    elif calcPar.CellMatrix is not None:
        Cell_Dict = {
            'a': [calcPar.CellMatrix[0][0], calcPar.CellMatrix[0][1], calcPar.CellMatrix[0][2]],
            'b': [calcPar.CellMatrix[1][0], calcPar.CellMatrix[1][1], calcPar.CellMatrix[1][2]],
            'c': [calcPar.CellMatrix[2][0], calcPar.CellMatrix[2][1], calcPar.CellMatrix[2][2]],
            'periodic': 'XYZ'
            }
    else:
        raise ValueError('Either the cell parameters or the cell matrix must be provided')

    Global_Dict = {
        "project_name": calcPar.FrameworkName,
        "run_type": calcPar.CalcType.lower(),
    }

    if calcPar.UseScalapack:
        Global_Dict["preferred_diag_library"] = "scalapack"

    Vibrational_Analysis_Dict = {
        'print': {'program_run_info': {'_': 'ON'}},
        'nproc_rep': calcPar.ProcsPerReplica,
        'dx': calcPar.dX,
        'fully_periodic': True,
        'intensities': True
        }

    Force_Eval_Dict = {
                "+dft": {
                    "+qs": {
                        'eps_default': calcPar.EPSDefault,
                        },
                    "+print": {
                        "+hirshfeld": {"_": "OFF"},
                        "+lowdin": {"_": "OFF"},
                        "+mulliken": {"_": "OFF"},
                    },
                    "+scf": {
                        "scf_guess": calcPar.SCFGuess,
                        "max_scf": calcPar.MaxSCFycles,
                        "eps_scf": calcPar.SCFConvergence,
                        "+mixing": {"method": calcPar.MixingMethod,
                                    "alpha": calcPar.MixingAlpha},
                        "+outer_scf": {"max_scf": calcPar.MaxOuterSCFycles,
                                       "eps_scf": calcPar.SCFConvergence}
                    },
                    "charge": calcPar.Charge,
                    "multiplicity": calcPar.Multiplicity
                },
                "+subsys": {
                    "+cell": Cell_Dict,
                    "+coord": Coord_Dict,
                    "+print": {'+symmetry': {'symmetry_elements': True}},
                },
                "stress_tensor": "analytical"
            }

    if calcPar.KPoints:
        if calcPar.KPoints is True:
            calcPar.KPoints = get_kgrid(calcPar.CellMatrix, dist=calcPar.RecDist)

        Force_Eval_Dict["+dft"]['+kpoints'] = {
            "scheme": ('MONKHORST-PACK', str(calcPar.KPoints[0]), str(calcPar.KPoints[1]), str(calcPar.KPoints[2])),
            "symmetry": True,
            "full_grid": True,
            "verbose": True,
            "parallel_group_size": -1,
            "eps_geo": 1e-3,
            }

        if calcPar.KPoints == 'auto':
            calcPar.KPoints = get_kgrid(calcPar.CellMatrix, dist=calcPar.RecDist)

    if calcPar.CalcType.lower() == 'energy_force':
        Force_Eval_Dict['+print'] = {
            "+forces": {"filename": "forces", "_": "ON"},
            "+stress_tensor": {"_": "ON"}
            }

    if calcPar.Functional == 'XTB':
        Force_Eval_Dict['+dft']['+qs'] = {
                        'method': 'XTB',
                        '+XTB': {
                            'check_atomic_charges': calcPar.CheckAtomicCharges,
                            'do_ewald': True,
                            '+parameter': {'dispersion_parameter_file': 'dftd3.dat'},
                        },
                    }

    if calcPar.Functional == 'PBE':
        Force_Eval_Dict["+dft"]['+xc'] = {
                        "+xc_functional": {
                            "+pbe": {"parametrization": calcPar.Parametrization}
                            },
                        "+vdw_potential": {
                            "potential_type": "pair_potential",
                            "+pair_potential": {
                                "type": calcPar.DispersionCorrection,
                                "reference_functional": calcPar.Functional,
                                "r_cutoff": 16,
                                "parameter_file_name": "dftd3.dat"
                                }
                            }
                        }
        Force_Eval_Dict["+dft"]['+mgrid'] = {
            'cutoff': calcPar.PWCutoff,
            'ngrids': calcPar.NGrid,
            'rel_cutoff': calcPar.RelativeCutOff
            }

        Force_Eval_Dict["+dft"]["basis_set_file_name"] = [
            "BASIS_MOLOPT",
            "BASIS_MOLOPT_UZH"
            ]

        Force_Eval_Dict["+dft"]["potential_file_name"] = "GTH_POTENTIALS"

        Force_Eval_Dict["+subsys"]["+kind"] = Kind_List

    if calcPar.Functional == 'PBE0':
        Force_Eval_Dict["+dft"]['+xc'] = {
                        "+xc_functional": {
                            "_": calcPar.Functional
                        },
                        "+vdw_potential": {
                            "potential_type": "pair_potential",
                            "+pair_potential": {
                                "type": calcPar.DispersionCorrection,
                                "reference_functional": calcPar.Functional,
                                "r_cutoff": 16,
                                "parameter_file_name": "dftd3.dat"
                                }
                            }
                        }
        Force_Eval_Dict["+dft"]['+mgrid'] = {
            'cutoff': calcPar.PWCutoff,
            'ngrids': calcPar.NGrid,
            'rel_cutoff': calcPar.RelativeCutOff
            }

        Force_Eval_Dict["+dft"]["basis_set_file_name"] = [
            "BASIS_MOLOPT",
            "BASIS_MOLOPT_UZH"
            ]

        Force_Eval_Dict["+dft"]["potential_file_name"] = "GTH_POTENTIALS"

        Force_Eval_Dict["+subsys"]["+kind"] = Kind_List

    if calcPar.UseOT:
        Force_Eval_Dict["+dft"]['+scf']["+ot"] = {"minimizer": "DIIS",
                                                  "n_diis": 7,
                                                  "preconditioner": "FULL_SINGLE_INVERSE"}

    if calcPar.UseSmearing:
        if calcPar.SmearingMethod == 'fermi_dirac':
            Force_Eval_Dict["+dft"]['+scf']['+smear'] = {
                "method": 'FERMI_DIRAC',
                "electronic_temperature": calcPar.ElectronicTemperature
            }
        elif calcPar.SmearingMethod == 'energy_window':
            Force_Eval_Dict["+dft"]['+scf']['+smear'] = {
                "method": 'energy_window',
                "width": calcPar.WindowSize
            }
        if calcPar.AddedMOs == 0:
            calcPar.AddedMOs = 50

        Force_Eval_Dict["+dft"]['+scf']['added_mos'] = calcPar.AddedMOs

    motion_dict = {
        "+print": [
            {
                "+forces": {"+each": {"cell_opt": 1, "geo_opt": 1, "md": 1}},
                "+cell": {"+each": {"cell_opt": 1, "geo_opt": 1, "md": 1}},
                "+trajectory": {"+each": {"cell_opt": 1, "geo_opt": 1, "md": 1}},
                "+velocities": {"+each": {"cell_opt": 1, "geo_opt": 1, "md": 1}},
                "+stress": {"+each": {"cell_opt": 1, "geo_opt": 1, "md": 1}},
                "+restart": {"+each": {"cell_opt": 1, "geo_opt": 1, "md": 1}, "backup_copies": 0},
                "+restart_history": {"_": "OFF"}
            }
        ]
    }

    if calcPar.CalcType.lower() == 'cell_opt':
        motion_dict['+cell_opt'] = {
            "+lbfgs": {"trust_radius": 0.25},
            "optimizer": "lbfgs",
            "max_iter": calcPar.MaxIterations,
            "max_dr": calcPar.MaxDR,
            "max_force": calcPar.MaxForce,
            "rms_dr": calcPar.RMSDR,
            "rms_force": calcPar.RMSForce
        }

        if calcPar.KeepSymmetry:
            motion_dict['+cell_opt']['keep_symmetry'] = True
            motion_dict['+cell_opt']['keep_space_group'] = True
            motion_dict['+cell_opt']['keep_angles'] = True

    if calcPar.CalcType.lower() == 'geo_opt':
        motion_dict['+geo_opt'] = {
            "+bfgs": {"trust_radius": 0.25},
            "max_iter": calcPar.MaxIterations,
            "max_dr": calcPar.MaxDR,
            "max_force": calcPar.MaxForce,
            "rms_dr": calcPar.RMSDR,
            "rms_force": calcPar.RMSForce
        }

    if calcPar.CalcType.lower() == 'md':
        motion_dict['+md'] = {
            "ensemble": calcPar.Ensemble,
            "temperature": calcPar.Temperature,
            "timestep": calcPar.TimeStep,
            "steps": calcPar.MDSteps,
            "+barostat": {
                "pressure": calcPar.Pressure,
                "timecon": calcPar.TimeCon
            },
            "+thermostat": {
                "type": 'CSVR',
                "+csvr": {'timecon': 0.1},
            }
        }

    if calcPar.CalculateRaman:
        Force_Eval_Dict["+properties"] = {
            'linres': {'polar': {'do_raman': True},
                       'max_iter': 200,
                       'preconditioner': 'full_all',
                       'eps': 1e-08
                       },
            }

    if calcPar.CalculateIR:
        Force_Eval_Dict['+dft']['+print']['+moments'] = {"periodic": True}

    input_dict = {
        "+global": Global_Dict,
        "+force_eval": [Force_Eval_Dict]
    }

    if calcPar.CalcType.lower() in ['cell_opt', 'geo_opt', 'md']:
        input_dict['+motion'] = motion_dict

    if calcPar.CalcType.lower() == 'normal_modes':
        input_dict['+vibrational_analysis'] = Vibrational_Analysis_Dict

    if calcPar.Restart:
        input_dict['+ext_restart'] = {
            "restart_file_name": f"{calcPar.FrameworkName}-1.restart"
        }

    generator = CP2KInputGenerator()

    with open(os.path.join(output_folder, FrameworkName), "w") as fhandle:
        for line in generator.line_iter(input_dict):
            fhandle.write(f"{line}\n")


def get_forces(FrameworkName, output_folder):
    """ Get the CP2K forces from the output file in atomic units [Hartree/a.u.]

    Parameters
    ----------
    FrameworkName : str
        Name of the framework
    output_folder : str
        Path to the output folder

    Returns
    -------
    forces : np.ndarray
        Nx3 Array of the forces in atomic units [a.u.]
    """

    with open(os.path.join(output_folder, f"{FrameworkName}-forces-1_0.xyz"), "r") as f:
        lines = f.read().splitlines()

    forces = []

    for line in lines[4:-1]:
        forces.append([float(i) for i in line.split()[3:]])

    if np.any(np.isnan(forces)) or np.any(np.isinf(forces)):
        print(f'Warning: Found NaN or Inf values on forces for {output_folder}')

    return np.array(forces)


def get_pol_tensor(file_name, output_folder, symmetrize=False):
    """
    Get the polarizability tensor from the CP2K output file.

    Parameters
    ----------
    file_name : str
        Name of the output file.
    output_folder : str
        Path to the output folder.
    symmetrize : bool, optional
        Symmetrize the polarizability tensor.

    Returns
    -------
    pol_au_order : np.ndarray
        Polarizability tensor (3,3) in atomic units [a.u.^3].
    po_angs_order : np.ndarray
        Polarizability tensor (3,3) in [angs^3].
    """
    with open(os.path.join(output_folder, file_name), "r") as f:
        lines = f.read().splitlines()

    pol_au = np.zeros((3, 3))
    pol_angs = np.zeros((3, 3))

    for i, line in enumerate(lines):
        if "POLARIZABILITY TENSOR (atomic units):" in line:
            t1 = lines[i+1].split()[1:]
            t2 = lines[i+2].split()[1:]
            t3 = lines[i+3].split()[1:]
            pol_au = np.array([t1, t2, t3], dtype=float)

        if "POLARIZABILITY TENSOR (Angstrom^3):" in line:

            t1 = lines[i+1].split()[1:]
            t2 = lines[i+2].split()[1:]
            t3 = lines[i+3].split()[1:]
            pol_angs = np.array([t1, t2, t3], dtype=float)

    # The output of CP2K is out of order. We need to reorder it
    pol_au_order = np.zeros((3, 3))
    pol_au_order[0, 0] = pol_au[0, 0]
    pol_au_order[1, 1] = pol_au[0, 1]
    pol_au_order[2, 2] = pol_au[0, 2]
    pol_au_order[0, 1] = pol_au[1, 0]
    pol_au_order[0, 2] = pol_au[1, 1]
    pol_au_order[1, 2] = pol_au[1, 2]
    pol_au_order[1, 0] = pol_au[2, 0]
    pol_au_order[2, 0] = pol_au[2, 1]
    pol_au_order[2, 1] = pol_au[2, 2]

    po_angs_order = np.zeros((3, 3))

    po_angs_order[0, 0] = pol_angs[0, 0]
    po_angs_order[1, 1] = pol_angs[0, 1]
    po_angs_order[2, 2] = pol_angs[0, 2]
    po_angs_order[0, 1] = pol_angs[1, 0]
    po_angs_order[0, 2] = pol_angs[1, 1]
    po_angs_order[1, 2] = pol_angs[1, 2]
    po_angs_order[1, 0] = pol_angs[2, 0]
    po_angs_order[2, 0] = pol_angs[2, 1]
    po_angs_order[2, 1] = pol_angs[2, 2]

    if symmetrize:
        pol_au_order = 0.5 * (pol_au_order + pol_au_order.T)
        po_angs_order = 0.5 * (po_angs_order + po_angs_order.T)

    if np.linalg.norm(pol_au_order) == 0:
        print(f'Warning: Found zero values on polarizability tensor of {output_folder}')

    if np.any(np.isnan(pol_au_order)) or np.any(np.isinf(pol_au_order)):
        print(f'Warning: Found NaN or Inf values on polarizability tensor of {output_folder}')
        pol_au_order = np.zeros((3, 3))
        po_angs_order = np.zeros((3, 3))

    return pol_au_order, po_angs_order


def diff_cross_section(w, laser_wl=532, T=298):
    """
    Calculate the differential cross section of a Raman scattering process for
    a given frequency, laser wavelength, and temperature.

    Taken from: The Raman Effect: A Unified Treatment of the Theory of Raman Scattering
    Equation 5.7.16

    Parameters
    ----------
    w : float
        Frequency of the Raman scattering process in cm^-1.
    laser_wl : float, optional
        Wavelength of the laser in nm. Default is 532 nm.
    T : float, optional
        Temperature in Kelvin. Default is 298 K.

    Returns
    -------
    cross_section : float
        Differential cross section of the Raman scattering process.
    """

    # Convert the laser wavelength to wavenumber
    wl = np.reciprocal(laser_wl * 1e-7)

    n_m = np.reciprocal(1 - np.exp(-h * w * 1e2 * c / (k_B * T)))

    return ((wl - w)**4 / w) * (n_m + 1)


def expand_tensor_by_symm(tensor, primitive, prim_symmetry):
    """
    Expand the tensor to all atoms in the primitive cell.

    Parameters
    ----------
    tensor : np.ndarray
        Tensor to be expanded.
    primitive : Atoms
        Primitive cell.
    prim_symmetry : Symmetry
        Symmetry of the primitive cell.

    Returns
    -------
    tensor : np.ndarray
        Expanded tensor.
    """

    # Expand tensor to all atoms in the primitive cell
    rotations = prim_symmetry.get_symmetry_operations()['rotations']
    map_operations = prim_symmetry.get_map_operations()
    map_atoms = prim_symmetry.get_map_atoms()

    for na in range(primitive.get_number_of_atoms()):
        # R_cart = L R L^-1
        rotc = similarity_transformation(primitive.get_cell().transpose(),
                                         rotations[map_operations[na]])
        rotct = rotc.transpose()
        # R_cart^T B R_cart^-1 (inverse rotation is required to transform)
        # 3d rotational matrix is R_{kij}=sum_{lmn}rotc_{kl}*rotc_{im}*rotc_{jn}
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    for p in range(3):
                        for m in range(3):
                            for n in range(3):
                                tensor[na][k][i][j] += (
                                    tensor[map_atoms[na]][p][m][n]
                                    * rotct[k][p]
                                    * rotct[i][m]
                                    * rotct[j][n])

    return tensor
