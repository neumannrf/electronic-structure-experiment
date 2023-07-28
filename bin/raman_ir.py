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
parser.add_argument('--Charge',
                    type=int,
                    default=0,
                    action='store',
                    required=False,
                    metavar='CHARGE',
                    help='Total charge of the unit cell.')
parser.add_argument('--Multipliticy',
                    type=int,
                    default=1,
                    action='store',
                    required=False,
                    metavar='MULTIPLICITY',
                    help='Total multiplicity of the unit cell.')
parser.add_argument('--MaxSCFycles',
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
                    choices=['PBE', 'xTB'],
                    metavar='FUNCTIONAL',
                    help='Functional used to calculate the total energy.')
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
                    choices=['DZVP', 'TZV2P'],
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
        'basis_set': BASIS_SET[arg.BasisSet][specie]
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
                    "_": arg.Functional
                },
                "+vdw_potential": {
                    "potential_type": "pair_potential",
                    "+pair_potential": {
                        "type": arg.DispersionCorrection,
                        "reference_functional": arg.Functional,
                        "r_cutoff": 16,
                        "parameter_file_name": f"{arg.CP2KDataDir}/dftd3.dat"
                        }
                    }
                },
            "charge": arg.Charge,
            "multiplicity": arg.Multipliticy,
            "basis_set_file_name": [
                f"{arg.CP2KDataDir}/BASIS_MOLOPT",
                f"{arg.CP2KDataDir}/BASIS_MOLOPT_UZH"],
            "potential_file_name": f"{arg.CP2KDataDir}/GTH_POTENTIALS",
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

if arg.UseOT:
    Force_Eval_Dict["+dft"]['+scf']["+ot"] = {"minimizer": "DIIS",
                                              "n_diis": 7,
                                              "preconditioner": "FULL_SINGLE_INVERSE"}

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
