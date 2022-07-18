import sys
import argparse

#############################################################################
#   This file contains a function which is the input parser for Vibham
#
#   To read all of the options and information for each input variable, type:
#       "python3 VibHam.py -h"
#
#   To read adjust an input variable, type:
#       "python3 VibHam -$VAR -$OPTION"
#       
#   where $VAR is the input variable you wish to change to $OPTION
############################################################################

def Input():

    parser = argparse.ArgumentParser(allow_abbrev=False,
                                     description='Dyanamical Simulations of a Water Molecule in an HCP ParaHydrogen Matrix'
                                     )

    parser.version='1.0'

    parser.add_argument('-shells', 
                       action='store', 
                       help='Number of pH2 shells to simulate', 
                       type=int,
                       default=1,
                       choices=[1, 2, 3]
                       )

    parser.add_argument('-distance', 
                        action='store', 
                        help='R value to separate the pH2 matrix',
                        type=float,
                        default=3.79
                        )

    parser.add_argument('-J',
                        action='store',
                        help='Maximum J-value for the rotation of the H2O molecule',
                        type=int,
                        default=1
                        )

    parser.add_argument('-state',
                        action='store',
                        help='Vibrational state to use for the H2O molecule (0-000, 1-100, 2-010, 3-001)',
                        type=int,
                        choices=[0, 1, 2, 3],
                        default=0
                        )

    parser.add_argument('-polyDeg', 
                        action='store', 
                        help='Degree of the polynomial used to fit the Lebedev coefficients',
                        type=int,
                        default=12
                        )

    parser.add_argument('-threshold',
                        action='store',
                        help='Cutoff threshold for the Hamiltonian matrix elements',
                        type=float,
                        default=1e-8
                        )

    parser.add_argument('-lebThresh',
                        action='store',
                        help='Cutoff threshold for the Lebedev Interpolation',
                        type=float,
                        default=1e-4
                        )

    parser.add_argument('-rotQuad',
                        action='store',
                        help='Maximum rotational quadrature values for the H2O molecule',
                        type=int,
                        nargs=3,
                        default=[4,4,4]
                        )

    parser.add_argument('-smearH',
                        action='store',
                        help='Maximum Gaussian quadrature value for the ZPM of an H-atom',
                        type=int,
                        default=3
                        )

    parser.add_argument('-smearH2',
                        action='store',
                        help='Maximum Gaussian quadrature value for the ZPM of a pH2 molecule',
                        type=int,
                        default=3
                        )

    parser.add_argument('-sigH',
                        action='store',
                        help='Sigma value for the gaussian quadrature of the ZPM of an H-atom',
                        type=float,
                        default=0.5563
                        )

    parser.add_argument('-sigH2',
                        action='store',
                        help='Sigma value for the gaussian quadrature of the ZPM of a pH2 molecule',
                        type=float,
                        default=0.4875265
                        )

    parser.add_argument('-h_atoms',
                        action='store',
                        help='Positions of H-atoms to replace pH2 molecules',
                        type=int,
                        nargs="+",
                        default=[]
                        )

    parser.add_argument('-holes',
                        action='store',
                        help='Positions of holes to replace pH2 molecules',
                        type=int,
                        nargs="+",
                        default=[]
                        )

    parser.add_argument('-print',
                        action='store',
                        help='Print level (1 or 2)',
                        type=int,
                        nargs="+",
                        default=[1, 10, 10]
                        )

    parser.add_argument('-stark',
                        action='store_true',
                        help='Perform a Stark calculation to test the methodology',
                        default=False
                        )

    parser.add_argument('-dipole',
                        action='store',
                        help='Dipole moment of water (D)',
                        type=float,
                        default=1.8546
                        )

    parser.add_argument('-environment',
                        action='store',
                        help='Perform a simulation using a given lattice environment',
                        choices=['fcc', 'hcp'],
                        default='hcp',
                        type=str.lower
                        )

    parser.add_argument('-time',
                        action='store_true',
                        help='Perform a runtime estimation for a simulation',
                        default=False
                        )

    parser.add_argument('-h_test',
                        action='store_true',
                        help="Don't read in Hamiltonian file",
                        default=False
                        )

    parser.add_argument('-effect',
                        action='store',
                        help='Amount of Potential Energy Matrix',
                        type=float,
                        default=1.
                        )

    parser.add_argument('-load',
                        action='store_true',
                        help="Attempt to load in previous datafiles",
                        default=False
                        )

    parser.add_argument('-method',
                        action='store',
                        help='Choose the code to use (python or fortran77)',
                        choices=['fortran', 'python'],
                        default='python',
                        type=str.lower
                        )

    
    return parser.parse_args()
