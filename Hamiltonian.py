import numpy as np
import os
import time
import math

class Potential:
    '''Class designed to construct and populate the rotational Hamiltonian matrix



    '''

    def __init__(self, args, angleEnergyArray):
        self.args = args
        self.angleEnergyArray = angleEnergyArray

        self.dim = int((4.0/3.0) * self.args.J**3 + 4*self.args.J**2 + (11.0/3.0) * self.args.J + 1.0000001)

        if self.args.load:
            self.lowJHamil = self.loadHamil()

            if self.lowJHamil.size == self.dim**2:
                self.Hamil = self.lowJHamil
            else:
                self.GeneratePreHamil()
                self.GenerateHamil()

        else:
            self.lowJHamil = np.zeros((1,1))
            self.GeneratePreHamil()
            self.GenerateHamil()

        lowvals = abs(self.Hamil) < self.args.threshold
        self.Hamil[lowvals] = 0

    def loadHamil(self):
        for j in range(self.args.J, -1, -1):
            if os.path.exists(str(j) + ".hamil"):
                print ("\n\tLower J Hamiltonian found with J = ", j)
                return np.loadtxt(str(j) + ".hamil")
        return np.zeros((1,1))

    def GeneratePreHamil(self):

        print ("\n\tConstructing pre-Hamiltonian arrays\n")

        self.dim = int((4.0/3.0) * self.args.J**3 + 4*self.args.J**2 + (11.0/3.0) * self.args.J + 1.0000001)

        self.PreHamil = np.zeros((self.dim, self.dim, 12))
    
        col = 0

        for J in range(0, self.args.J+1):
            Bra_Norm = np.sqrt((2 * J + 1) / (8 * np.pi**2))
            for K in range(-J, J+1):
                for M in range(-J, J+1):
                    Bra_Nmin = max(0, K-M)
                    Bra_Nmax = min((J-M), (J+K))

                    self.PreHamil[col, :, 0] = J
                    self.PreHamil[col, :, 1] = K
                    self.PreHamil[col, :, 2] = M
                    self.PreHamil[col, :, 3] = Bra_Norm
                    self.PreHamil[col, :, 4] = Bra_Nmin
                    self.PreHamil[col, :, 5] = Bra_Nmax

                    self.PreHamil[:, col, 6 ] = J
                    self.PreHamil[:, col, 7 ] = K
                    self.PreHamil[:, col, 8 ] = M
                    self.PreHamil[:, col, 9 ] = Bra_Norm
                    self.PreHamil[:, col, 10] = Bra_Nmin
                    self.PreHamil[:, col, 11] = Bra_Nmax

                    col += 1

    def GenerateHamil(self):

        start = time.time()

        time_check = False

        self.Hamil = np.zeros((self.dim, self.dim))

        tmpHamil = np.zeros((self.dim, self.dim))

        if self.args.load:
            #TODO

            for col in range(self.lowJHamil.shape[0]):
                for row in range(self.lowJHamil.shape[1]):
                    if self.lowJHamil[col, row] != 0.0:
                        tmpHamil[col, row] = self.lowJHamil[col, row]


        if self.args.method == 'python':

            def WignerD(J, K, M, nmin, nmax, beta):
                def fac(x):
                    return math.factorial(x)
                def cos(x):
                    return math.cos(x)
                def sin(x):
                    return math.sin(x)
                def sqrt(x):
                    return math.sqrt(x)

                lilD = 0

                for N in range(nmin, nmax+1):

                    lilw = sqrt(fac(J+M) * fac(J-M) * fac(J+K) * fac(J-K)) / (fac(J-M-N) * fac(J+K-N) * fac(N+M-K) * fac(N))
                    BigW = lilw * ((cos(beta/2.0))**(2*J+K-M-2*N)) * ((-sin(beta/2.0))**(M-K+2*N))

                    lilD += ((-1)**N) * BigW

                return lilD

            def PrintTime(val, string):
                '''Function used to print time value'''
                minutes, seconds = divmod(val, 60)
                hours,   minutes = divmod(minutes, 60)
                days,    hours   = divmod(hours, 24)
                print ("\t" + string + " - ", "%d days, %d hours, %d minutes, %d seconds" %(days, hours, minutes, seconds), "\n")
                return

            for col in range(self.dim):
                J        = int(self.PreHamil[col, 0, 0])
                K        = int(self.PreHamil[col, 0, 1])
                M        = int(self.PreHamil[col, 0, 2])
                braNorm  =     self.PreHamil[col, 0, 3]
                braNmin = int(self.PreHamil[col, 0, 4])
                braNmax = int(self.PreHamil[col, 0, 5])

                for row in range(col, self.dim):
                    JJ          = int(self.PreHamil[0, row, 6])
                    KK          = int(self.PreHamil[0, row, 7])
                    MM          = int(self.PreHamil[0, row, 8])
                    ketNorm     =     self.PreHamil[0, row, 9]
                    ketNmin    = int(self.PreHamil[0, row, 10])
                    ketNmax    = int(self.PreHamil[0, row, 11])

                    if tmpHamil[col, row] != 0.:
                        self.Hamil[col, row] = tmpHamil[col, row]
                        self.Hamil[row, col] = tmpHamil[col, row]

                    else:
                        start_ = time.time()

                        H = 0.

                        for angle in range(self.angleEnergyArray.shape[0]):
                            alpha  = self.angleEnergyArray[angle, 0]
                            beta   = self.angleEnergyArray[angle, 1]
                            gamma  = self.angleEnergyArray[angle, 2]

                            totalWeight = self.angleEnergyArray[angle, 3] * self.angleEnergyArray[angle, 4] * self.angleEnergyArray[angle, 5]

                            energy = self.angleEnergyArray[angle, 6]

                            braLilD = WignerD(J,  K,  M,  braNmin, braNmax, beta)
                            ketLilD = WignerD(JJ, KK, MM, ketNmin, ketNmax, beta)

                            braBigD = np.exp(-(M*alpha + K*gamma) * complex(0,1)) * braLilD
                            ketBigD = np.exp(-(MM*alpha + KK*gamma) * complex(0,1)) * ketLilD

                            bra = braBigD.conjugate() * braNorm
                            ket = ketBigD * ketNorm

                            H += totalWeight * energy * bra * ket

                        self.Hamil[col, row] = H.real
                        self.Hamil[row, col] = H.real

                        if time_check == False:
                            stop = time.time()
                            val = (stop - start_)

                            corner = int(self.dim**2 / 2. + self.dim / 2.)

                            precount = np.count_nonzero(tmpHamil)

                            val *= self.dim**2 - precount
                            val *= (corner / self.dim**2)

                            PrintTime(val, "Estimated time for Populating Rotational Hamiltonian Matrix")

                            if self.args.time == True:
                                exit()
                            time_check = True

                print ("\t%4d out of %4d entries in the Rotational Hamiltonian Matrix Populated" % (row + col*self.dim + 1, self.dim**2))

            print ("\n\n")
            print ("\tConstruction of Hamiltonian matrix completed.")

            np.savetxt(str(self.args.J) + ".hamil", self.Hamil)

            stop = time.time()
            val = round((stop - start), 3)
            PrintTime(val, "Total time")
            print ()

        elif self.args.method == 'fortran':
            #TODO

            self


class Kinetic:

    def __init__(self, args):
        self.args = args

        self.dim = int((4.0/3.0) * self.args.J**3 + 4*self.args.J**2 + (11.0/3.0) * self.args.J + 1.0000000001)

        self.Energy()

    def Energy(self):

        if self.args.state == 0 :
        # v1, v2, v3 = 0, 0, 0
            A = 27.8572554962
            B = 14.5144875043
            C = 9.27986250769

        elif self.args.state == 1:
        # v1, v2, v3 = 1, 0, 0
            A = 27.1388599924
            B = 14.2991000050
            C = 9.10149945389

        elif self.args.state == 2:
        # v1, v2, v3 = 0, 1, 0
            A = 31.0847300062
            B = 14.6748099834
            C = 9.13606001363

        elif self.args.state == 3:
        # v1, v2, v3 = 0, 0, 1
            A = 26.6303599946
            B = 14.4225400054
            C = 9.14182999037

        IULC = 0

        self.kineticEnergyMatrix = np.zeros((self.dim, self.dim))

        for J in range(0,self.args.J+1):
            for M in range(-J, J+1):
                for K in range(-J, J+1):
                    for KK in range(-J, J+1):

                        val = 0.0

                        if K == KK:
                            val = ((0.50) * (B+C) * (J*(J+1)-K**2) + A*K**2)
                        elif K + 2 == KK:
                            val = ((0.25) * (B-C) * np.sqrt((J*(J+1) - K * (K+1))) * np.sqrt(J*(J+1) - (K+1) * (K+2)))
                        elif K - 2 == KK:
                            val = ((0.25) * (B-C) * np.sqrt((J*(J+1) - K * (K-1))) * np.sqrt(J*(J+1) - (K-1) * (K-2)))

                        try:
                            self.kineticEnergyMatrix[IULC+KK+J, IULC+K+J] = val
                        except IndexError:
                            pass

                IULC = IULC + 2*J+1

