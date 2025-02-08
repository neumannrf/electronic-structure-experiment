#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import os
import warnings

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

    with open(os.path.join(outputfolder, FrameworkName + '-frc-1.xyz'), 'r') as f:
        lines = f.read().splitlines()

    n_atoms = int(lines[0])

    # Reshape lines to be (n_atoms + 2, -1)
    forces_list = [i[2:] for i in np.array(lines).reshape((-1, n_atoms + 2))]

    # Remove the atom label (first column) and convert the forces to float
    forces_list = [np.array([np.array(i.split()[1:]).astype(float) for i in structure]) for structure in forces_list]

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
                      CalcType: str = 'energy_force',
                      Charge: int = 0,
                      Multiplicity: int = 1,
                      UseOT: bool = False,
                      UseSmearing: bool = False,
                      SmearingMethod: str = 'fermi_dirac',
                      ElectronicTemperature: int = 300,
                      WindowSize: float = 0.1,
                      AddedMOs: int = 0,
                      MixingMethod: str = 'broyden_mixing',
                      MixingAlpha: float = 0.2,
                      MaxSCFcycles: int = 30,
                      MaxOuterSCFycles: int = 10,
                      EPSDefault: float = 1e-8,
                      PWCutoff: int = 1200,
                      NGrid: int = 5,
                      RelativeCutOff: int = 60,
                      Functional: str = 'PBE',
                      Parametrization: str = 'ORIG',
                      DispersionCorrection: str = 'DFTD3',
                      CheckAtomicCharges: bool = True,
                      BasisSet: str = 'DZVP',
                      SCFGuess: str = 'atomic',
                      SCFConvergence: float = 1e-8,
                      CP2KDataDir: str = None,
                      KeepSymmetry: bool = False,
                      KeepSpaceGroup: bool = False,
                      KeepAngles: bool = False,
                      MaxIterations: int = 100,
                      Restart: bool = False,
                      MaxDR: float = 1e-3,
                      MaxForce: float = 1e-3,
                      RMSDR: float = 1e-3,
                      RMSForce: float = 1e-3,
                      UseScalapack: bool = False,
                      CellParameters: str = None,
                      CellMatrix: str = None,
                      AtomicTypes: str = None,
                      FracX: float = None,
                      FracY: float = None,
                      FracZ: float = None,
                      CartX: float = None,
                      CartY: float = None,
                      CartZ: float = None,
                      ProcsPerReplica: int = 4,
                      dX: float = 0.001,
                      CalculateRaman: bool = False,
                      CalculateIR: bool = False,
                      Ensemble: str = 'NPT_F',
                      Temperature: int = 400,
                      TimeStep: float = 0.5,
                      MDSteps: int = 100,
                      Pressure: int = 1,
                      TimeCon: int = 1000,
                      KPoints: bool = False,
                      RecDist: float = 0.3) -> None:
    """
    Create the input file for CP2K

    Parameters
    ----------
    FrameworkName : str
        Name of the framework
    output_folder : str
        Path to the output folder
    CalcType : str, optional
        Type of calculation. Can be 'energy_force', 'cell_opt', 'geo_opt', 'md', or 'normal_modes'
    Charge : int, optional
        Charge of the system. Default is 0
    Multiplicity : int, optional
        Multiplicity of the system. Default is 1
    UseOT : bool, optional
        Use the OT method. Default is False
    UseSmearing : bool, optional
        Use smearing method. Default is False
    SmearingMethod : str, optional
        Smearing method to be used. Can be 'fermi_dirac' or 'energy_window'. Default is 'fermi_dirac'
    ElectronicTemperature : int, optional
        Electronic temperature in Kelvin. Default is 300
    WindowSize : float, optional
        Window size for smearing. Default is 0.1
    AddedMOs : int, optional
        Number of added MOs. Default is 0
    MixingMethod : str, optional
        Mixing method. Can be 'direct_p_mixing', 'broyden_mixing_new', or 'kerker_mixing'. Default is 'broyden_mixing'
    MixingAlpha : float, optional
        Mixing alpha parameter. Default is 0.2
    MaxSCFcycles : int, optional
        Maximum number of SCF cycles. Default is 30
    MaxOuterSCFycles : int, optional
        Maximum number of outer SCF cycles. Default is 10
    EPSDefault : float, optional
        Default convergence criterion for SCF. Default is 1e-8
    PWCutoff : int, optional
        Plane wave cutoff. Default is 1200
    NGrid : int, optional
        Grid density. Default is 5
    RelativeCutOff : int, optional
        Relative cutoff for potentials. Default is 60
    Functional : str, optional
        Exchange-correlation functional. Can be 'PBE', 'XTB', or 'PBE0'. Default is 'PBE'
    Parametrization : str, optional
        Functional parametrization. Can be 'ORIG', 'PBESOL', or 'REVPBE'. Default is 'ORIG'
    DispersionCorrection : str, optional
        Dispersion correction method. Can be None, 'DFTD2', 'DFTD3', or 'DFTD3(BJ)'. Default is 'DFTD3'
    CheckAtomicCharges : bool, optional
        Whether to check atomic charges. Default is True
    BasisSet : str, optional
        Basis set to use. Can be 'SZV', 'DZVP', 'TZVP', or 'TZV2P'. Default is 'DZVP'
    SCFGuess : str, optional
        SCF guess method. Can be 'atomic', 'restart', 'core', 'random', 'sparse', or 'mopac'. Default is 'atomic'
    SCFConvergence : float, optional
        SCF convergence criterion. Default is 1e-8
    CP2KDataDir : str, optional
        Directory for CP2K data. If None, will use environment variable 'CP2K_DATA_DIR'
    KeepSymmetry : bool, optional
        Whether to keep symmetry. Default is False
    KeepSpaceGroup : bool, optional
        Whether to keep space group information. Default is False
    KeepAngles : bool, optional
        Whether to keep angles. Default is False
    MaxIterations : int, optional
        Maximum number of iterations. Default is 100
    Restart : bool, optional
        Whether to restart from a previous calculation. Default is False
    MaxDR : float, optional
        Maximum change in coordinates. Default is 1e-3
    MaxForce : float, optional
        Maximum force on atoms. Default is 1e-3
    RMSDR : float, optional
        RMS deviation of coordinates. Default is 1e-3
    RMSForce : float, optional
        RMS force on atoms. Default is 1e-3
    UseScalapack : bool, optional
        Whether to use ScaLAPACK. Default is False
    CellParameters : str, optional
        Cell parameters. Default is None
    CellMatrix : str, optional
        Cell matrix. Default is None
    AtomicTypes : str, optional
        Atomic types. Default is None
    FracX : float, optional
        Fractional coordinate in X direction. Default is None
    FracY : float, optional
        Fractional coordinate in Y direction. Default is None
    FracZ : float, optional
        Fractional coordinate in Z direction. Default is None
    CartX : float, optional
        Cartesian coordinate in X direction. Default is None
    CartY : float, optional
        Cartesian coordinate in Y direction. Default is None
    CartZ : float, optional
        Cartesian coordinate in Z direction. Default is None
    ProcsPerReplica : int, optional
        Number of processors per replica. Default is 4
    dX : float, optional
        Increment for coordinate adjustments. Default is 0.001
    CalculateRaman : bool, optional
        Whether to calculate Raman spectra. Default is False
    CalculateIR : bool, optional
        Whether to calculate IR spectra. Default is False
    Ensemble : str, optional
        MD ensemble type. Can be 'NVE', 'NVT', NPT_I', and 'NPT_F'. Default is 'NPT_F'
    Temperature : int, optional
        Temperature for MD. Default is 400
    TimeStep : float, optional
        Time step for MD simulations. Default is 0.5
    MDSteps : int, optional
        Number of MD steps. Default is 100
    Pressure : int, optional
        Pressure for MD. Default is 1
    TimeCon : int, optional
        Time constant for MD. Default is 1000
    KPoints : bool, optional
        Whether to use k-points. Default is False
    RecDist : float, optional
        Reciprocal distance for k-points. Default is 0.3
    """

    Coord_Dict = {
        'scaled': False,
        '*': ['{:3} {:11.6f} {:11.6f} {:11.6f}'.format(AtomicTypes[i],
                                                       CartX[i],
                                                       CartY[i],
                                                       CartZ[i]) for i in range(len(AtomicTypes))]
                    }

    Kind_List = []

    for specie in set(AtomicTypes):
        Kind_List.append(
            {
                "_": specie,
                'element': specie,
                'potential': PSEUDO_POTENTIALS[specie],
                'basis_set': BASIS_SET[BasisSet][specie]
            }
        )

    if CellParameters is not None:
        Cell_Dict = {
            'abc': [CellParameters[0], CellParameters[1], CellParameters[2]],
            'alpha_beta_gamma': [CellParameters[3], CellParameters[4], CellParameters[5]],
            'periodic': 'XYZ'
            }
    elif CellMatrix is not None:
        Cell_Dict = {
            'a': [CellMatrix[0][0], CellMatrix[0][1], CellMatrix[0][2]],
            'b': [CellMatrix[1][0], CellMatrix[1][1], CellMatrix[1][2]],
            'c': [CellMatrix[2][0], CellMatrix[2][1], CellMatrix[2][2]],
            'periodic': 'XYZ'
            }
    else:
        raise ValueError('Either the cell parameters or the cell matrix must be provided')

    Global_Dict = {
        "project_name": FrameworkName,
        "run_type": CalcType.lower(),
    }

    if UseScalapack:
        Global_Dict["preferred_diag_library"] = "scalapack"

    Vibrational_Analysis_Dict = {
        'print': {'program_run_info': {'_': 'ON'}},
        'nproc_rep': ProcsPerReplica,
        'dx': dX,
        'fully_periodic': True,
        'intensities': True
        }

    Force_Eval_Dict = {
                "+dft": {
                    "+qs": {
                        'eps_default': EPSDefault,
                        },
                    "+print": {
                        "+hirshfeld": {"_": "OFF"},
                        "+lowdin": {"_": "OFF"},
                        "+mulliken": {"_": "OFF"},
                    },
                    "+scf": {
                        "scf_guess": SCFGuess,
                        "max_scf": MaxSCFcycles,
                        "eps_scf": SCFConvergence,
                        "+mixing": {"method": MixingMethod,
                                    "alpha": MixingAlpha},
                        "+outer_scf": {"max_scf": MaxOuterSCFycles,
                                       "eps_scf": SCFConvergence}
                    },
                    "charge": Charge,
                    "multiplicity": Multiplicity
                },
                "+subsys": {
                    "+cell": Cell_Dict,
                    "+coord": Coord_Dict,
                    "+print": {'+symmetry': {'symmetry_elements': True}},
                },
                "stress_tensor": "analytical"
            }

    if KPoints:
        if KPoints is True:
            KPoints = get_kgrid(CellMatrix, dist=RecDist)

        Force_Eval_Dict["+dft"]['+kpoints'] = {
            "scheme": ('MONKHORST-PACK', str(KPoints[0]), str(KPoints[1]), str(KPoints[2])),
            "symmetry": True,
            "full_grid": True,
            "verbose": True,
            "parallel_group_size": -1,
            "eps_geo": 1e-3,
            }

        if KPoints == 'auto':
            KPoints = get_kgrid(CellMatrix, dist=RecDist)

    if CalcType.lower() == 'energy_force':
        Force_Eval_Dict['+print'] = {
            "+forces": {"filename": "forces", "_": "ON"},
            "+stress_tensor": {"_": "ON"}
            }

    if Functional == 'XTB':
        Force_Eval_Dict['+dft']['+qs'] = {
                        'method': 'XTB',
                        '+XTB': {
                            'check_atomic_charges': CheckAtomicCharges,
                            'do_ewald': True,
                            '+parameter': {'dispersion_parameter_file': 'dftd3.dat'},
                        },
                    }

    if Functional == 'PBE':
        Force_Eval_Dict["+dft"]['+xc'] = {
                        "+xc_functional": {
                            "+pbe": {"parametrization": Parametrization}
                            },
                        "+vdw_potential": {
                            "potential_type": "pair_potential",
                            "+pair_potential": {
                                "type": DispersionCorrection,
                                "reference_functional": Functional,
                                "r_cutoff": 16,
                                "parameter_file_name": "dftd3.dat"
                                }
                            }
                        }
        Force_Eval_Dict["+dft"]['+mgrid'] = {
            'cutoff': PWCutoff,
            'ngrids': NGrid,
            'rel_cutoff': RelativeCutOff
            }

        Force_Eval_Dict["+dft"]["basis_set_file_name"] = [
            "BASIS_MOLOPT",
            "BASIS_MOLOPT_UZH"
            ]

        Force_Eval_Dict["+dft"]["potential_file_name"] = "GTH_POTENTIALS"

        Force_Eval_Dict["+subsys"]["+kind"] = Kind_List

    if Functional == 'PBE0':
        Force_Eval_Dict["+dft"]['+xc'] = {
                        "+xc_functional": {
                            "_": Functional
                        },
                        "+vdw_potential": {
                            "potential_type": "pair_potential",
                            "+pair_potential": {
                                "type": DispersionCorrection,
                                "reference_functional": Functional,
                                "r_cutoff": 16,
                                "parameter_file_name": "dftd3.dat"
                                }
                            }
                        }
        Force_Eval_Dict["+dft"]['+mgrid'] = {
            'cutoff': PWCutoff,
            'ngrids': NGrid,
            'rel_cutoff': RelativeCutOff
            }

        Force_Eval_Dict["+dft"]["basis_set_file_name"] = [
            "BASIS_MOLOPT",
            "BASIS_MOLOPT_UZH"
            ]

        Force_Eval_Dict["+dft"]["potential_file_name"] = "GTH_POTENTIALS"

        Force_Eval_Dict["+subsys"]["+kind"] = Kind_List

    if UseOT:
        Force_Eval_Dict["+dft"]['+scf']["+ot"] = {"minimizer": "DIIS",
                                                  "n_diis": 7,
                                                  "preconditioner": "FULL_SINGLE_INVERSE"}

    if UseSmearing:
        if SmearingMethod == 'fermi_dirac':
            Force_Eval_Dict["+dft"]['+scf']['+smear'] = {
                "method": 'FERMI_DIRAC',
                "electronic_temperature": ElectronicTemperature
            }
        elif SmearingMethod == 'energy_window':
            Force_Eval_Dict["+dft"]['+scf']['+smear'] = {
                "method": 'energy_window',
                "width": WindowSize
            }
        if AddedMOs == 0:
            AddedMOs = 50

        Force_Eval_Dict["+dft"]['+scf']['added_mos'] = AddedMOs

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

    if CalcType.lower() == 'cell_opt':
        motion_dict['+cell_opt'] = {
            "+lbfgs": {"trust_radius": 0.25},
            "optimizer": "lbfgs",
            "max_iter": MaxIterations,
            "max_dr": MaxDR,
            "max_force": MaxForce,
            "rms_dr": RMSDR,
            "rms_force": RMSForce,
            "keep_symmetry": KeepAngles,
            "keep_space_group": KeepSpaceGroup,
            "keep_angles": KeepAngles
        }

    if CalcType.lower() == 'geo_opt':
        motion_dict['+geo_opt'] = {
            "+bfgs": {"trust_radius": 0.25},
            "max_iter": MaxIterations,
            "max_dr": MaxDR,
            "max_force": MaxForce,
            "rms_dr": RMSDR,
            "rms_force": RMSForce
        }

    if CalcType.lower() == 'md':
        motion_dict['+md'] = {
            "ensemble": Ensemble,
            "temperature": Temperature,
            "timestep": TimeStep,
            "steps": MDSteps,
            "+barostat": {
                "pressure": Pressure,
                "timecon": TimeCon
            },
            "+thermostat": {
                "type": 'CSVR',
                "+csvr": {'timecon': 0.1},
            }
        }

    if CalculateRaman:
        Force_Eval_Dict["+properties"] = {
            'linres': {'polar': {'do_raman': True},
                       'max_iter': 200,
                       'preconditioner': 'full_all',
                       'eps': 1e-08
                       },
            }

    if CalculateIR:
        Force_Eval_Dict['+dft']['+print']['+moments'] = {"periodic": True}

    input_dict = {
        "+global": Global_Dict,
        "+force_eval": [Force_Eval_Dict]
    }

    if CalcType.lower() in ['cell_opt', 'geo_opt', 'md']:
        input_dict['+motion'] = motion_dict

    if CalcType.lower() == 'normal_modes':
        input_dict['+vibrational_analysis'] = Vibrational_Analysis_Dict

    if Restart:
        input_dict['+ext_restart'] = {
            "restart_file_name": f"{FrameworkName}-1.restart"
        }

    generator = CP2KInputGenerator()

    with open(os.path.join(output_folder, FrameworkName), "w") as fhandle:
        for line in generator.line_iter(input_dict):
            fhandle.write(f"{line}\n")


def get_forces(FileName, output_folder):
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

    with open(os.path.join(output_folder, FileName), "r") as f:
        lines = f.read().splitlines()

    forces_list = []

    for line in lines[4:-1]:
        try:
            forces = np.array(line.split()[3:]).astype(float)
        except Exception:
            print(f'Error: Could not read forces for {FileName}')
            forces = np.zeros(3)
        if np.any(np.isnan(forces)) or np.any(np.isinf(forces)):
            print(f'Warning: Found NaN or Inf values on forces for {FileName}')

        forces_list.append(forces)

    return np.array(forces_list)


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
        warnings.warn(f'Warning: Found zero values on polarizability tensor of {output_folder}')

    if np.any(np.isnan(pol_au_order)) or np.any(np.isinf(pol_au_order)):
        warnings.warn(f'Warning: Found NaN or Inf values on polarizability tensor of {output_folder}')
        pol_au_order = np.zeros((3, 3))
        po_angs_order = np.zeros((3, 3))

    return pol_au_order, po_angs_order


def calc_placzek_invariants(frequencies, alpha):
    """
    Calculate the Raman Tensor Placzek Invariants of the polarizability tensor.
    This method follows the procedure described in:
    The Raman Effect: A Unified Treatment of the Theory of Raman Scattering by Molecules
    by Derek A. Long, 2002, Section A14.7.5, page 490.

    This method does not require the tensor to be symmetric.

    Parameters
    ----------
    alpha : np.ndarray
        Polarizability tensor in atomic units [a.u.^3].

    Returns
    -------
    a_sq : np.ndarray
        Mean polarizability squared.
    gamma_sq : np.ndarray
        Anisotropy squared.
    delta_sq : np.ndarray
        Asymmetric anisotropy squared
    """

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

    return a_sq, gamma_sq, delta_sq


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
