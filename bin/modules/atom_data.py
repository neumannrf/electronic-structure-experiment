#!/usr/bin/env -S python -B

# Copyright IBM Corp. 2023
# SPDX-License-Identifier: Apache2.0

# valence number for PBE-GTH potentials and DZVP-MOLOPT-SR-GTH basis
VALENCE_NUM = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 3, 'C': 4, 'N': 5, 'O': 6, 'F': 7, 'Ne': 8, 'Na': 9, 'Mg': 10,
    'Al': 3, 'Si': 4, 'P': 5, 'S': 6, 'Cl': 7, 'Ar': 8, 'K': 9, 'Ca': 10, 'Sc': 11, 'Ti': 12, 'V': 13, 'Cr': 14,
    'Mn': 15, 'Fe': 16, 'Co': 17, 'Ni': 18, 'Cu': 11, 'Zn': 12, 'Ga': 13, 'Ge': 4, 'As': 5, 'Se': 6, 'Br': 7, 'Kr': 8,
    'Rb': 9, 'Sr': 10, 'Y': 11, 'Zr': 12, 'Nb': 13, 'Mo': 14, 'Tc': 15, 'Ru': 16, 'Rh': 17, 'Pd': 18, 'Ag': 11,
    'Cd': 12, 'In': 13, 'Sn': 4, 'Sb': 5, 'Te': 6, 'I': 7, 'Xe': 8, 'Cs': 9, 'Ba': 10, 'Hf': 12, 'Ta': 13, 'W': 14,
    'Re': 15, 'Os': 16, 'Ir': 17, 'Pt': 18, 'Au': 11, 'Hg': 12, 'Tl': 13, 'Pb': 4, 'Bi': 5, 'Po': 6, 'At': 7, 'Rn': 8}

# core number for PBE-GTH potentials and DZVP-MOLOPT-SR-GTH basis
CORE_NUM = {'H': 0, 'He': 0, 'Li': 0, 'Be': 0, 'B': 2, 'C': 2, 'N': 2, 'O': 2, 'F': 2, 'Ne': 2, 'Na': 2, 'Mg': 2,
            'Al': 10, 'Si': 10, 'P': 10, 'S': 10, 'Cl': 10, 'Ar': 10, 'K': 10, 'Ca': 10, 'Sc': 10, 'Ti': 10, 'V': 10,
            'Cr': 10, 'Mn': 10, 'Fe': 10, 'Co': 10, 'Ni': 10, 'Cu': 18, 'Zn': 18, 'Ga': 18, 'Ge': 28, 'As': 28,
            'Se': 28, 'Br': 28, 'Kr': 28, 'Rb': 28, 'Sr': 28, 'Y': 28, 'Zr': 28, 'Nb': 28, 'Mo': 28, 'Tc': 28,
            'Ru': 28, 'Rh': 28, 'Pd': 28, 'Ag': 36, 'Cd': 36, 'In': 36, 'Sn': 46, 'Sb': 46, 'Te': 46, 'I': 46,
            'Xe': 46, 'Cs': 46, 'Ba': 46, 'Hf': 60, 'Ta': 60, 'W': 60, 'Re': 60, 'Os': 60, 'Ir': 60, 'Pt': 60,
            'Au': 68, 'Hg': 68, 'Tl': 68, 'Pb': 78, 'Bi': 78, 'Po': 78, 'At': 78, 'Rn': 78}

ATOMIC_NUMBER = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11,
                 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti':
                 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32,
                 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42,
                 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52,
                 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62,
                 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72,
                 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
                 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92,
                 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102,
                 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111,
                 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116, 'Ts': 117, 'Og': 118}

# Basis set names for DZVP-MOLOPT-SR-GTH basis
BASIS_SET = {'DZVP': {
    'H': 'DZVP-MOLOPT-SR-GTH-q1', 'He': 'DZVP-MOLOPT-SR-GTH-q2', 'Li': 'DZVP-MOLOPT-SR-GTH-q3',
    'Be': 'DZVP-MOLOPT-SR-GTH-q4', 'B': 'DZVP-MOLOPT-SR-GTH-q3', 'C': 'DZVP-MOLOPT-SR-GTH-q4',
    'N': 'DZVP-MOLOPT-SR-GTH-q5', 'O': 'DZVP-MOLOPT-SR-GTH-q6', 'F': 'DZVP-MOLOPT-SR-GTH-q7',
    'Ne': 'DZVP-MOLOPT-SR-GTH-q8', 'Na': 'DZVP-MOLOPT-SR-GTH-q9', 'Mg': 'DZVP-MOLOPT-SR-GTH-q10',
    'Al': 'DZVP-MOLOPT-SR-GTH-q3', 'Si': 'DZVP-MOLOPT-SR-GTH-q4', 'P': 'DZVP-MOLOPT-SR-GTH-q5',
    'S': 'DZVP-MOLOPT-SR-GTH-q6', 'Cl': 'DZVP-MOLOPT-SR-GTH-q7', 'Ar': 'DZVP-MOLOPT-SR-GTH-q8',
    'K': 'DZVP-MOLOPT-SR-GTH-q9', 'Ca': 'DZVP-MOLOPT-SR-GTH-q10', 'Sc': 'DZVP-MOLOPT-SR-GTH-q11',
    'Ti': 'DZVP-MOLOPT-SR-GTH-q12', 'V': 'DZVP-MOLOPT-SR-GTH-q13', 'Cr': 'DZVP-MOLOPT-SR-GTH-q14',
    'Mn': 'DZVP-MOLOPT-SR-GTH-q15', 'Fe': 'DZVP-MOLOPT-SR-GTH-q16', 'Co': 'DZVP-MOLOPT-SR-GTH-q17',
    'Ni': 'DZVP-MOLOPT-SR-GTH-q18', 'Cu': 'DZVP-MOLOPT-SR-GTH-q11', 'Zn': 'DZVP-MOLOPT-SR-GTH-q12',
    'Ga': 'DZVP-MOLOPT-SR-GTH-q13', 'Ge': 'DZVP-MOLOPT-SR-GTH-q4', 'As': 'DZVP-MOLOPT-SR-GTH-q5',
    'Se': 'DZVP-MOLOPT-SR-GTH-q6', 'Br': 'DZVP-MOLOPT-SR-GTH-q7', 'Kr': 'DZVP-MOLOPT-SR-GTH-q8',
    'Rb': 'DZVP-MOLOPT-SR-GTH-q9', 'Sr': 'DZVP-MOLOPT-SR-GTH-q10', 'Y': 'DZVP-MOLOPT-SR-GTH-q11',
    'Zr': 'DZVP-MOLOPT-SR-GTH-q12', 'Nb': 'DZVP-MOLOPT-SR-GTH-q13', 'Mo': 'DZVP-MOLOPT-SR-GTH-q14',
    'Tc': 'DZVP-MOLOPT-SR-GTH-q15', 'Ru': 'DZVP-MOLOPT-SR-GTH-q16', 'Rh': 'DZVP-MOLOPT-SR-GTH-q17',
    'Pd': 'DZVP-MOLOPT-SR-GTH-q18', 'Ag': 'DZVP-MOLOPT-SR-GTH-q11', 'Cd': 'DZVP-MOLOPT-SR-GTH-q12',
    'In': 'DZVP-MOLOPT-SR-GTH-q13', 'Sn': 'DZVP-MOLOPT-SR-GTH-q4', 'Sb': 'DZVP-MOLOPT-SR-GTH-q5',
    'Te': 'DZVP-MOLOPT-SR-GTH-q6', 'I': 'DZVP-MOLOPT-SR-GTH-q7', 'Xe': 'DZVP-MOLOPT-SR-GTH-q8',
    'Cs': 'DZVP-MOLOPT-SR-GTH-q9', 'Ba': 'DZVP-MOLOPT-SR-GTH-q10', 'Hf': 'DZVP-MOLOPT-SR-GTH-q12',
    'Ta': 'DZVP-MOLOPT-SR-GTH-q13', 'W': 'DZVP-MOLOPT-SR-GTH-q14', 'Re': 'DZVP-MOLOPT-SR-GTH-q15',
    'Ir': 'DZVP-MOLOPT-SR-GTH-q17', 'Pt': 'DZVP-MOLOPT-SR-GTH-q18', 'Au': 'DZVP-MOLOPT-SR-GTH-q11',
    'Hg': 'DZVP-MOLOPT-SR-GTH-q12', 'Tl': 'DZVP-MOLOPT-SR-GTH-q13', 'Pb': 'DZVP-MOLOPT-SR-GTH-q4',
    'Bi': 'DZVP-MOLOPT-SR-GTH-q5', 'Po': 'DZVP-MOLOPT-SR-GTH-q6', 'At': 'DZVP-MOLOPT-SR-GTH-q7',
    'Rn': 'DZVP-MOLOPT-SR-GTH-q8', 'U': 'DZVP-MOLOPT-GTH-q14'
    },
        'TZV2P': {
    'H': 'TZV2P-MOLOPT-PBE-GTH-q1', 'He': 'TZV2P-MOLOPT-PBE-GTH-q2', 'Li': 'TZV2P-MOLOPT-PBE-GTH-q3',
    'Be': 'TZV2P-MOLOPT-PBE-GTH-q4', 'B': 'TZV2P-MOLOPT-PBE-GTH-q3', 'C': 'TZV2P-MOLOPT-PBE-GTH-q4',
    'N': 'TZV2P-MOLOPT-PBE-GTH-q5', 'O': 'TZV2P-MOLOPT-PBE-GTH-q6', 'F': 'TZV2P-MOLOPT-PBE-GTH-q7',
    'Ne': 'TZV2P-MOLOPT-PBE-GTH-q8', 'Na': 'TZV2P-MOLOPT-PBE-GTH-q9', 'Mg': 'TZV2P-MOLOPT-PBE-GTH-q10',
    'Al': 'TZV2P-MOLOPT-PBE-GTH-q3', 'Si': 'TZV2P-MOLOPT-PBE-GTH-q4', 'P': 'TZV2P-MOLOPT-PBE-GTH-q5',
    'S': 'TZV2P-MOLOPT-PBE-GTH-q6', 'Cl': 'TZV2P-MOLOPT-PBE-GTH-q7', 'Ar': 'TZV2P-MOLOPT-PBE-GTH-q8',
    'K': 'TZV2P-MOLOPT-PBE-GTH-q9', 'Ca': 'TZV2P-MOLOPT-PBE-GTH-q10', 'Sc': 'TZV2P-MOLOPT-PBE-GTH-q11',
    'Ti': 'TZV2P-MOLOPT-PBE-GTH-q12', 'V': 'TZV2P-MOLOPT-PBE-GTH-q13', 'Cr': 'TZV2P-MOLOPT-PBE-GTH-q14',
    'Mn': 'TZV2P-MOLOPT-PBE-GTH-q15', 'Fe': 'TZV2P-MOLOPT-PBE-GTH-q16', 'Co': 'TZV2P-MOLOPT-PBE-GTH-q17',
    'Ni': 'TZV2P-MOLOPT-PBE-GTH-q18', 'Cu': 'TZV2P-MOLOPT-PBE-GTH-q11', 'Zn': 'TZV2P-MOLOPT-PBE-GTH-q12',
    'Ga': 'TZV2P-MOLOPT-PBE-GTH-q13', 'Ge': 'TZV2P-MOLOPT-PBE-GTH-q4', 'As': 'TZV2P-MOLOPT-PBE-GTH-q5',
    'Se': 'TZV2P-MOLOPT-PBE-GTH-q6', 'Br': 'TZV2P-MOLOPT-PBE-GTH-q7', 'Kr': 'TZV2P-MOLOPT-PBE-GTH-q8',
    'Rb': 'TZV2P-MOLOPT-PBE-GTH-q9', 'Sr': 'TZV2P-MOLOPT-PBE-GTH-q10', 'Y': 'TZV2P-MOLOPT-PBE-GTH-q11',
    'Zr': 'TZV2P-MOLOPT-PBE-GTH-q12', 'Nb': 'TZV2P-MOLOPT-PBE-GTH-q13', 'Mo': 'TZV2P-MOLOPT-PBE-GTH-q14',
    'Tc': 'TZV2P-MOLOPT-PBE-GTH-q15', 'Ru': 'TZV2P-MOLOPT-PBE-GTH-q16', 'Rh': 'TZV2P-MOLOPT-PBE-GTH-q17',
    'Pd': 'TZV2P-MOLOPT-PBE-GTH-q18', 'Ag': 'TZV2P-MOLOPT-PBE-GTH-q11', 'Cd': 'TZV2P-MOLOPT-PBE-GTH-q12',
    'In': 'TZV2P-MOLOPT-PBE-GTH-q13', 'Sn': 'TZV2P-MOLOPT-PBE-GTH-q4', 'Sb': 'TZV2P-MOLOPT-PBE-GTH-q5',
    'Te': 'TZV2P-MOLOPT-PBE-GTH-q6', 'I': 'TZV2P-MOLOPT-PBE-GTH-q7', 'Xe': 'TZV2P-MOLOPT-PBE-GTH-q8',
    'Cs': 'TZV2P-MOLOPT-PBE-GTH-q9', 'Ba': 'TZV2P-MOLOPT-PBE-GTH-q10', 'Hf': 'TZV2P-MOLOPT-PBE-GTH-q12',
    'Ta': 'TZV2P-MOLOPT-PBE-GTH-q13', 'W': 'TZV2P-MOLOPT-PBE-GTH-q14', 'Re': 'TZV2P-MOLOPT-PBE-GTH-q15',
    'Ir': 'TZV2P-MOLOPT-PBE-GTH-q17', 'Pt': 'TZV2P-MOLOPT-PBE-GTH-q18', 'Au': 'TZV2P-MOLOPT-PBE-GTH-q11',
    'Hg': 'TZV2P-MOLOPT-PBE-GTH-q12', 'Tl': 'TZV2P-MOLOPT-PBE-GTH-q13', 'Pb': 'TZV2P-MOLOPT-PBE-GTH-q4',
    'Bi': 'TZV2P-MOLOPT-PBE-GTH-q5', 'Po': 'TZV2P-MOLOPT-PBE-GTH-q6', 'At': 'TZV2P-MOLOPT-PBE-GTH-q7',
    'Rn': 'TZV2P-MOLOPT-PBE-GTH-q8', 'U': 'TZV2P-MOLOPT-GTH-q14'
    }
    }

# Potential names for the PBE-GTH potentials
PSEUDO_POTENTIALS = {
    'H': 'GTH-PBE-q1', 'He': 'GTH-PBE-q2', 'Li': 'GTH-PBE-q3', 'Be': 'GTH-PBE-q4', 'B': 'GTH-PBE-q3',
    'C': 'GTH-PBE-q4', 'N': 'GTH-PBE-q5', 'O': 'GTH-PBE-q6', 'F': 'GTH-PBE-q7', 'Ne': 'GTH-PBE-q8',
    'Na': 'GTH-PBE-q9', 'Mg': 'GTH-PBE-q10', 'Al': 'GTH-PBE-q3', 'Si': 'GTH-PBE-q4', 'P': 'GTH-PBE-q5',
    'S': 'GTH-PBE-q6', 'Cl': 'GTH-PBE-q7', 'Ar': 'GTH-PBE-q8', 'K': 'GTH-PBE-q9', 'Ca': 'GTH-PBE-q10',
    'Sc': 'GTH-PBE-q11', 'Ti': 'GTH-PBE-q12', 'V': 'GTH-PBE-q13', 'Cr': 'GTH-PBE-q14', 'Mn': 'GTH-PBE-q15',
    'Fe': 'GTH-PBE-q16', 'Co': 'GTH-PBE-q17', 'Ni': 'GTH-PBE-q18', 'Cu': 'GTH-PBE-q11', 'Zn': 'GTH-PBE-q12',
    'Ga': 'GTH-PBE-q13', 'Ge': 'GTH-PBE-q4', 'As': 'GTH-PBE-q5', 'Se': 'GTH-PBE-q6', 'Br': 'GTH-PBE-q7',
    'Kr': 'GTH-PBE-q8', 'Rb': 'GTH-PBE-q9', 'Sr': 'GTH-PBE-q10', 'Y': 'GTH-PBE-q11', 'Zr': 'GTH-PBE-q12',
    'Nb': 'GTH-PBE-q13', 'Mo': 'GTH-PBE-q14', 'Tc': 'GTH-PBE-q15', 'Ru': 'GTH-PBE-q16', 'Rh': 'GTH-PBE-q17',
    'Pd': 'GTH-PBE-q18', 'Ag': 'GTH-PBE-q11', 'Cd': 'GTH-PBE-q12', 'In': 'GTH-PBE-q13', 'Sn': 'GTH-PBE-q4',
    'Sb': 'GTH-PBE-q5', 'Te': 'GTH-PBE-q6', 'I': 'GTH-PBE-q7', 'Xe': 'GTH-PBE-q8', 'Cs': 'GTH-PBE-q9',
    'Ba': 'GTH-PBE-q10', 'Hf': 'GTH-PBE-q12', 'Ta': 'GTH-PBE-q13', 'W': 'GTH-PBE-q14', 'Re': 'GTH-PBE-q15',
    'Ir': 'GTH-PBE-q17', 'Pt': 'GTH-PBE-q18', 'Au': 'GTH-PBE-q11', 'Hg': 'GTH-PBE-q12', 'Tl': 'GTH-PBE-q13',
    'Pb': 'GTH-PBE-q4', 'Bi': 'GTH-PBE-q5', 'Po': 'GTH-PBE-q6', 'At': 'GTH-PBE-q7', 'Rn': 'GTH-PBE-q8',
    'U': 'GTH-PBE-q14'}
