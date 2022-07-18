import numpy as np
import pathlib
import sys
import time
import os

from Input import Input
from Environment import HCP, FCC
from Print import PrintSimulationInfo, PrintXYZ, PrintFinalInfo, PrintEnd
from InteractionEnergies import InteractionEnergies
from Hamiltonian import Potential, Kinetic

import Dynamical

class DynamicalSimulations:
    def __init__(self, args):
        self.args = args

        PrintSimulationInfo(args)

        self.start = time.time()


    def Diagonalize(self, M):
        '''Function used to diagonalize square matrices and
            sort the eigenvalues according to magnitude.

            Variables:
                val - Eigenvalues
                vec - Eigenvectors
            '''

        val, vec = np.linalg.eig(M)
        idx = val.argsort()
        val = val[idx]
        vec = vec[:,idx]
        return val, vec

    def Environment(self):
        '''Function used to control the environment around the 
            rotating H2O molecule'''

        if self.args.environment == 'fcc':
            self.paraH = FCC(self.args.shells, self.args.distance)
        elif self.args.environment == 'hcp':
            self.paraH = HCP(self.args.shells, self.args.distance)

        print ("\n\n")
        print ("\tGenerating the pH2", self.args.environment.upper(), "lattice.")
        print ("\tLattice contains", self.args.shells, "shell(s) at a distance of", self.args.distance, "Angstrom")
        print ()

        if bool(set(self.args.h_atoms) & set(self.args.holes)) == True:
            # Check to make sure no H-atoms and holes occupy in the same space
            print ("\tCannot have an H-atom and a hole at the same place")
            print ("\tTerminating pogram")
            exit()

        if len(self.args.h_atoms) > 0:
            try:
                while max(self.args.h_atoms) > self.paraH.shape[0]-1:
                    self.args.h_atoms = sorted(self.args.h_atoms)[:-1]
            except:
                pass
            print ("\tH-atoms   found at positions - ", *self.args.h_atoms)
        if len(self.args.holes) > 0:
            try:
                while max(self.args.holes) > self.paraH.shape[0]-1:
                    self.args.holes = sorted(self.args.holes)[:-1]
            except:
                pass
            print ("\tVacancies found at positions - ", *self.args.holes)
        print ()


        # Modify lattice environment to account for h-atoms and vacancies
        self.hAtoms = self.paraH[sorted(self.args.h_atoms)]
        self.args.holes[:] = [x - len(self.args.h_atoms) for x in self.args.holes]
        self.paraH = np.delete(self.paraH, self.args.h_atoms, 0)
        self.paraH = np.delete(self.paraH, self.args.holes, 0)

        # Print the array of the pH2 matrix and any H-atoms
        # Terminates program if args.Print == 'XYZ'
        PrintXYZ(self.paraH, self.hAtoms)

    def Hamiltonian(self):

        IE = InteractionEnergies(self.args, 
                                 self.paraH, 
                                 self.hAtoms,
                                 script_path)

        self.interactionEnergy = IE.angleArray


        HAMIL = Potential(self.args, 
                          self.interactionEnergy)

        self.potentialHamil = HAMIL.Hamil


        lowvals = abs(self.potentialHamil) < self.args.threshold
        self.potentialHamil[lowvals] = 0

        # Calculate Kinetic energy of water molecule in the absence of an external field
        # First Parameter is the Jmax value
        # Second parameter is the vibrational excited state

        print ("\tCalculating the Kinetic energy of the free water molecule for state", self.args.state, " with Jmax of", self.args.J, "\n")
        KINETIC = Kinetic(self.args)

        self.kineticHamil = KINETIC.kineticEnergyMatrix

        self.val, self.vec = self.Diagonalize(self.potentialHamil + self.kineticHamil)

        PrintFinalInfo(self.args.print, self.val, self.vec)
        PrintEnd(self.start)



def Main():
    '''Main Driver for the Dynamical Simulations'''

    args = Input()

    RunDynam = DynamicalSimulations(args)

    RunDynam.Environment()

    RunDynam.Hamiltonian()


if __name__ == "__main__":

    script_path = str(pathlib.Path(__file__).parent.resolve()) + "/"

    
    Main()
