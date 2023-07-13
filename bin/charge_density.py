#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

import argparse
import os

from cp2k_input_tools.generator import CP2KInputGenerator

from modules.calculate_properties import get_CellParameters, get_AtomicPositions
from modules.atom_data import BASIS_SET, PSEUDO_POTENTIALS

# Required parameters
parser = argparse.ArgumentParser(description='Create the CP2K simulation input.')
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
parser.add_argument('--CP2KDataDir',
                    type=str,
                    default=os.environ.get("CP2K_DATA_DIR"),
                    action='store',
                    required=False,
                    metavar='CP2K_DATA_DIR',
                    help='Directory containing the Basis set and pseudopotential files.')
parser.add_argument('--NetCharge',
                    type=int,
                    default=0,
                    action='store',
                    required=False,
                    metavar='NET_CHARGE',
                    help='Net charge of the unit cell.')
parser.add_argument('--EPSDefault',
                    type=float,
                    default=1e-7,
                    action='store',
                    required=False,
                    metavar='EPS_DEFAULT',
                    help='Default value for the electronic density convergence threshold.')
parser.add_argument('--PWCutoff',
                    type=float,
                    default=500,
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
parser.add_argument('--MaxSCFycles',
                    type=int,
                    default=50,
                    action='store',
                    required=False,
                    metavar='MAX_SCF_CYCLES',
                    help='Maximum number of SCF cycles.')
parser.add_argument('--MaxOuterSCFycles',
                    type=int,
                    default=10,
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
parser.add_argument('--UseOT',
                    action='store_true',
                    required=False,
                    help='Use the OT method for the SCF cycle.')
parser.add_argument('--SCFConvergence',
                    type=float,
                    default=1e-7,
                    action='store',
                    required=False,
                    metavar='SCF_CONVERGENCE',
                    help='SCF convergence threshold.')
parser.add_argument('--UseSmearing',
                    action='store_true',
                    required=False,
                    help='Use smearing for the electronic density.')
parser.add_argument('--SmearingMethod',
                    type=str,
                    default='fermi_dirac',
                    action='store',
                    required=False,
                    choices=['fermi_dirac',
                             'energy_window'],
                    metavar='SMEARING_METHOD',
                    help='Smearing method for the electronic density.')
parser.add_argument('--ElectronicTemperature',
                    type=float,
                    default=300.0,
                    action='store',
                    required=False,
                    metavar='ELECTRONIC_TEMPERATURE',
                    help='Electronic temperature in the case of Fermi-Dirac smearing.')
parser.add_argument('--WindowSize',
                    type=float,
                    default=0.1,
                    action='store',
                    required=False,
                    metavar='WINDOW_SIZE',
                    help='Size of the energy window centred at the Fermi level.')
parser.add_argument('--AddedMOs',
                    type=int,
                    default=0,
                    action='store',
                    required=False,
                    metavar='ADDED_MOS',
                    help='Number of unnocupied molecular orbitals added.')
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

# Parse the arguments
arg = parser.parse_args()

# Read the cif file and get the lattice parameters and atomic positions
cif_filename = os.path.join(arg.output_folder, arg.FrameworkName + '.cif')

CellParameters = get_CellParameters(cif_filename)
AtomicTypes, PosX, PosY, PosZ = get_AtomicPositions(cif_filename)

Coord_Dict = {
    'scaled': True,
    '*': [f'{AtomicTypes[i]:3} {PosX[i]:11.6f} {PosY[i]:11.6f} {PosZ[i]:11.6f}' for i in range(len(AtomicTypes))]
                }

Kind_List = []

for i, specie in enumerate(set(AtomicTypes)):
    Kind_List.append({
        "_": specie,
        'element': specie,
        'potential': PSEUDO_POTENTIALS[specie],
        'basis_set': BASIS_SET[specie]
    })

Cell_Dict = {'abc': [CellParameters[0], CellParameters[1], CellParameters[2]],
             'alpha_beta_gamma': [CellParameters[3], CellParameters[4], CellParameters[5]],
             'periodic': 'XYZ'}

Global_Dict = {
    "project_name": arg.FrameworkName,
    "run_type": "energy",
    "preferred_diag_library": "scalapack",
}

Force_Eval_Dict = {
        "+dft": {
            "+qs": {'eps_default': arg.EPSDefault},
            '+mgrid': {
                'cutoff': arg.PWCutoff,
                'ngrids': arg.NGrid,
                'rel_cutoff': arg.RelativeCutOff
                },
            "+print": {
                "+hirshfeld": {"_": "OFF"},
                "+lowdin": {"_": "OFF"},
                "+mulliken": {"_": "OFF"},
                "+E_DENSITY_CUBE": {"_": "ON", "filename": "valence_density", "stride": 1},
            },
            "+scf": {
                "scf_guess": arg.SCFGuess,
                "max_scf": arg.MaxSCFycles,
                "eps_scf": arg.SCFConvergence,
                "+mixing": {"method": arg.MixingMethod, "alpha": arg.MixingAlpha},
                "+outer_scf": {"max_scf": arg.MaxOuterSCFycles, "eps_scf": arg.SCFConvergence}
            },
            "+xc": {
                "+xc_functional": {
                    "_": "PBE"
                }
            },
            "charge": 0,
            "multiplicity": 1,
            "basis_set_file_name": f"{arg.CP2KDataDir}/BASIS_MOLOPT",
            "potential_file_name": f"{arg.CP2KDataDir}/GTH_POTENTIALS",
        },
        "+print": {
            "+forces": {"filename": "forces", "_": "ON"},
            "+stress_tensor": {"_": "ON"}
        },
        "+subsys": {
            "+cell": Cell_Dict,
            "+coord": Coord_Dict,
            "+kind": Kind_List
        },
        "method": "quickstep",
        "stress_tensor": "analytical"
    }

if arg.UseOT:
    Force_Eval_Dict["+dft"]['+scf']["+ot"] = {"minimizer": "DIIS", "n_diis": 7, "preconditioner": "FULL_SINGLE_INVERSE"}

if arg.UseSmearing:
    if arg.SmearingMethod == 'fermi_dirac':
        Force_Eval_Dict["+dft"]['+scf']['+smear'] = {
            "method": 'FERMI_DIRAC',
            "electronic_temperature": arg.ElectronicTemperature
        }
    elif arg.SmearingMethod == 'energy_window':
        Force_Eval_Dict["+dft"]['+scf']['+smear'] = {
            "method": 'energy_window',
            "width": arg.WindowSize
        }
    Force_Eval_Dict["+dft"]['+scf']['added_mos'] = arg.AddedMOs

input = {"+global": Global_Dict, "+force_eval": [Force_Eval_Dict]}

generator = CP2KInputGenerator()

with open("simulation_SCF.inp", "w") as fhandle:
    for line in generator.line_iter(input):
        fhandle.write(f"{line}\n")
