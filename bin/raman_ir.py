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
parser.add_argument('--MaxSCFycles',
                    type=int,
                    default=50,
                    action='store',
                    required=False,
                    metavar='MAX_SCF_CYCLES',
                    help='Maximum number of SCF cycles.')
parser.add_argument('--NetCharge',
                    type=int,
                    default=0,
                    action='store',
                    required=False,
                    metavar='NET_CHARGE',
                    help='Net charge of the unit cell.')
parser.add_argument('--MaxOuterSCFycles',
                    type=int,
                    default=5,
                    action='store',
                    required=False,
                    metavar='MAX_OUTER_SCF_CYCLES',
                    help='Maximum number of Outer SCF cycles for OT simulations.')
parser.add_argument('--DataFolder',
                    type=str,
                    default=os.getenv('CP2K_DATA'),
                    action='store',
                    required=False,
                    metavar='DATA_FOLDER',
                    help='Directory containing the Basis set and pseudopotential files.')
parser.add_argument('--ProcsPerReplica',
                    type=int,
                    default=1,
                    action='store',
                    required=False,
                    metavar='PROCS_PER_REPLICA',
                    help='Number of processors per replica. There is a total of 6N replicas (N = number of atoms).')
parser.add_argument('--dX',
                    type=float,
                    default=0.001,
                    action='store',
                    required=False,
                    metavar='DX',
                    help='Finite difference step size for numerical differentiation.')
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

arg = parser.parse_args()

# Calculate self-consistent properties
cif_filename = os.path.join(arg.output_folder, arg.FrameworkName + '.cif')

# Get the lattice parameters
CellParameters = get_CellParameters(cif_filename)

# Get the atomic positions
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
    "run_type": "normal_modes"
}

if arg.UseScalapack:
    Global_Dict["preferred_diag_library"] = "scalapack"

Vibrational_Analysis_Dict = {
    'print': {'program_run_info': {'_': 'ON'}},
    'nproc_rep': arg.ProcsPerReplica,
    'dx': arg.dX,
    'fully_periodic': True,
    'intensities': True
    }

Force_Eval_Dict = {
        "+dft": {
            "+qs": {'eps_default': 1e-10},
            '+mgrid': {
                'cutoff': 600.0,
                'ngrids': 5,
                'rel_cutoff': 60.0
                },
            "+print": {
                "+hirshfeld": {"_": "OFF"},
                "+lowdin": {"_": "OFF"},
                "+mulliken": {"_": "OFF"},
            },
            "+scf": {
                "+ot": {"preconditioner": "full_single_inverse", "minimizer": "diis"},
                "+mixing": {"method": "direct_p_mixing"},
                "scf_guess": "atomic",
                "max_scf": arg.MaxSCFycles,
                "eps_scf": 1e-10
            },
            'xc': {
                'xc_functional': {'pbe': {'parametrization': 'orig'}},
                'vdw_potential': {
                        'pair_potential': {
                            'parameter_file_name': 'dftd3.dat',
                            'type': 'dftd3(bj)',
                            'reference_functional': 'PBE',
                            'r_cutoff': 16.0
                            },
                        'potential_type': 'pair_potential'
                        }
                },
            "charge": arg.NetCharge,
            "multiplicity": 1,
            "basis_set_file_name": ['BASIS_MOLOPT', 'BASIS_MOLOPT_UCL', 'BASIS_ADMM'],
            "potential_file_name": "GTH_POTENTIALS",
        },
        "+print": {
            "+forces": {"filename": "forces", "_": "ON"},
            "+stress_tensor": {"_": "ON"}
        },
        "+subsys": {
            "+cell": Cell_Dict,
            "+coord": Coord_Dict,
            "+kind": Kind_List,
            "+print": {'+symmetry': {'symmetry_elements': True}},
        },
        "method": "quickstep",
        "stress_tensor": "analytical"
    }

if arg.CalculateRaman:
    Force_Eval_Dict["+properties"] = {
        'linres': {'polar': {'do_raman': True},
                   'max_iter': 200,
                   'preconditioner': 'full_all',
                   'eps': 1e-08
                   },
        }

if arg.CalculateIR:
    Force_Eval_Dict['+dft']['+print']['+moments'] = {"periodic": True}

input = {"+global": Global_Dict, "+force_eval": [Force_Eval_Dict], '+vibrational_analysis': Vibrational_Analysis_Dict}

generator = CP2KInputGenerator()

with open("simulation_Vibrations.inp", "w") as fhandle:
    for line in generator.line_iter(input):
        fhandle.write(f"{line}\n")
