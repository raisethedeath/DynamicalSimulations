import time
import numpy as np
from datetime import timedelta

def PrintSimulationInfo(*args):
    arg = args[0]
    print ()
    print ("\tSystem Information:\n")
    print ("\tEnvironment                          - ", arg.environment.upper())
    print ("\tMaximum J-value                      - ", arg.J)
    print ("\tDipole Moment                        - ", arg.dipole, "D")
    print ("\tpH2-pH2 Distance                     - ", arg.distance, "A")
    print ("\tLebedev Cutoff Threshold             - ", arg.lebThresh, "cm^-1")
    print ("\tLebedev Polynomial Fit Degree        - ", arg.polyDeg)
    print ("\tPrint Level                          - ", arg.print)
    print ("\tRotational Quadrature                - ", str(arg.rotQuad)[1:-1])
    print ("\tNumber of HCP Shells                 - ", arg.shells)
    print ("\tSigma Value for H-atom               - ", arg.sigH)
    print ("\tSigma Value fo pH2 molecule          - ", arg.sigH2)
    print ("\tGaussian Quadrature for H-atom       - ", arg.smearH)
    print ("\tGaussian Quadrature for pH2 molecule - ", arg.smearH2)
    print ("\tVibrational State                    - ", arg.state)
    print ("\tHamiltonian Cutoff Threshold         - ", arg.threshold, "cm^-1")
    print ("\n")


def PrintXYZ(*args):

    ParaH   = args[0]
    H_atoms = args[1]

    print ("\tPosition of pH2 molecules:")
    for i, ph2 in enumerate(ParaH):
        print ("\t", i+1, " - ", end=' ')
        for j in ph2:
            print (round(j,4), end='\t')
        print ()
    print()

    if len(H_atoms) != 0:
        print ("\tPosition of H-atoms:")
        for i, h in enumerate(H_atoms):
            print ("\t", i+1, " - ", end=' ')
            for j in h:
                print (round(j,4), end='\t')
            print ()
        print()

def PrintFinalInfo(*args):
    val = args[0][0]

    if val == 1:
        # PRINT N LOWEST EIGENVALUES

        eigval = args[1]
        numval = min(len(eigval), args[0][1])-1

        print ("\tThe following are the eigenvalues for the first", numval+1, "states:")
        print ()
        print ("\t{:>6s}{:>16s}{:>16s}".format("State", "Value", "dE"))

        for j in range(numval+1):
            if j == 0:
                print ("\t{:>6d}{:>16.8f}{:>16.8f}".format(j, eigval[j], 0))
            else:
                print ("\t{:>6d}{:>16.8f}{:>16.8f}".format(j, eigval[j], eigval[j] - eigval[0]))

    elif val == 2:
        # PRINT N LOWEST EIGENVALUES AND EIGENVECTORS

        eigval = args[1]
        eigvec = args[2]
        numval = min(len(eigval), args[0][1])-1
        numvec = min(len(eigval), args[0][2])-1

        vecval = len(eigval)-1

        print ("\tThe following are the eigenvalues and eigenvectors for the first", numval+1, "states:")
        print ()
        print ("\t{:>6s}{:>16s}{:>16s}".format("State", "Value", "dE"))
        print ("\t= = = = = = = = = = = = = = = = = = = =")
        print ("\t{:>6s}{:>16s}{:>16s}".format("Vector", "Contribution", "Weight"))
        print ()

        for j in range(numval+1):
            if j == 0:
                print ("\t{:>6d}{:>16.8f}{:>16.8f}".format(j, eigval[j], 0))
            else:
                print ("\t{:>6d}{:>16.8f}{:>16.8f}".format(j, eigval[j], eigval[j] - eigval[0]))
            print ("\t- - - - - - - - - - - - - - - - - - - -")
            for jj in range(numvec+1):
                print ("\t{:>6d}{:>16.8f}{:>16.8f}".format(jj, eigvec[jj][j], eigvec[jj][j]**2))
            print ()



def PrintEnd(*args):
    start = args[0]
    string = "Total Excecution Time:"
    print ("\n")
    print ("\tEnd of Simulation")
    e_Time = str(timedelta(seconds=start)).split(":")
    print ("\t" + string + " - ", "%s hours, %s minutes, %s seconds" %(e_Time[0], e_Time[1], e_Time[2]), "\n")
    
