import numpy as np
import os
import scipy.special
import time

class InteractionEnergies:
    '''Class desinged to calculate the interaction energies

        Functions:
            QuadValues          - Calculate Euler angle quadrature values and weights
            FitSurface          - Fit the potential energy surface using Lebedev quadrature
            GuassSmearing       - Function used to apply Gaussian smearing to the zero-point vibrational motion
                                   of the pH2 and/or H-atoms
            GetEnergyValues     - Calculate the interaction energy for all points of the simulation

    '''

    def __init__(self, args, paraH, hAtoms, path):
        self.args   = args
        self.paraH  = paraH
        self.hAtoms = hAtoms
        self.path = path
        self.jMax = 8

        if self.args.load and os.path.exists("Energy.en"):
            self.angleArray = np.loadtxt("Energy.en")

            print ("\tEnergy file loaded successfully")
            print ()
            print ("\tSkipping calculation related to:")
            print ("\t\t-Fitting of the potential energy surface")
            print ("\t\t-Gaussian quadrature (smearing)")
            print ("\t\t-Interaction energies")

        else:
            self.QuadValues(*args.rotQuad)
            self.FitSurface()
            self.GaussianSmearing()
            self.GetEnergyValues()


    def QuadValues(self, *rotQuad):
        '''Function used to calculate the Euler angle 
            quadrature values and weights

            Variables:
            numAValues  - Number of alpha Euler angle values
            numBValues  - Number of beta  Euler angle values
            numGValues  - Number of gamma Euler angle values
            angleArray  - Matrix of all Euler angle values and weights

        '''
        
        numAValues = rotQuad[0]
        numBValues = rotQuad[1]
        numGValues = rotQuad[2]

        betaVals, betaWeights = np.polynomial.legendre.leggauss(numBValues)
        betaVals  = np.arccos(betaVals)

        total_ang = numAValues * numBValues * numGValues

        self.angleArray = np.zeros((total_ang, 7))

        angle_count = 0

        for a in range(numAValues):
            alphaVal    = (2 * np.pi * a) / numAValues
            alphaWeight = ((2. * np.pi) / numAValues)

            for g in range(numGValues):
                gammaVal    = (2 * np.pi * g) / numGValues
                gammaWeight = ((2. * np.pi) / numGValues)

                for b in range(numBValues):
                    betaVal    = betaVals[b]
                    betaWeight = betaWeights[b]

                    self.angleArray[angle_count, 0] = alphaVal
                    self.angleArray[angle_count, 1] = betaVal
                    self.angleArray[angle_count, 2] = gammaVal
                    self.angleArray[angle_count, 3] = alphaWeight
                    self.angleArray[angle_count, 4] = betaWeight
                    self.angleArray[angle_count, 5] = gammaWeight

                    angle_count += 1


    def FitSurface(self):
        '''Function used to fit the potential energy surface using a Lebedev
            quadrature scheme.

            Functions:
                Coef        - Used to calculate the Lebedev expansion coefficients
        
            Vairables:
                self.h2Coef     - Array of Lebedev expansion coefficients for the H2O-H2 interaction
                self.hCoef      - Array of Lebedev expansion coefficients for the H2O-H  interaction
                R               - Array of monomer distances values
                WTP             - Matrix of weight, theta, and phi values for Lebedev Quadrature
                W               - Array of weight values for Lebedev Quadrature
                T               - Array of theta values for Lebedev Quadrature
                P               - Array of phi values for Lebedev Quadrature
                E               - Array of energy values for Lebedev Quadrature
        '''

        def Coef(T, P, W, E):

            count = 0
            C_jm = np.zeros(((self.jMax+1)**2))
            
            for j in range(0, self.jMax+1):
                for m in range(-j, j+1):
                    C = 0.0
                    for pts in range(len(E)):
                        theta, phi = T[pts], P[pts]
                        Sphere_harm = np.conj(scipy.special.sph_harm(m, j, phi, theta))
                        V_abinit    = E[pts]
                        Weight      = W[pts]

                        C += Sphere_harm * V_abinit * Weight

                    C_jm[count] = C.real
                    count += 1

            return C_jm

        if self.args.load:
            if os.path.exists('H2O-pH2.leb'):
                print ("\n\tPrevious H2O-pH2.leb file found for Lebedev expansion coefficients")
                self.h2Coef = np.loadtxt('H2O-pH2.leb')

            if os.path.exists('H2O-H.leb'):
                print ("\tPrevious H2O-H.leb   file found for Lebedev expansion coefficients")
                self.hCoef = np.loadtxt('H2O-H.leb')

        else:
            R = np.loadtxt(self.path + "/PES/R").transpose()

            WTP = np.loadtxt(self.path + "/PES/WTP").transpose()
            W, T, P = WTP[0], WTP[1], WTP[2]
            
            h2Data = np.loadtxt(self.path + "PES/H2O-pH2").transpose()[self.args.state]

            self.h2Coef = np.zeros((len(R), (self.jMax+1)**2))

            print ("\n\tFitting the Ab initio datapoints to a Lebedev Quadrature scheme for H2O-pH2")

            for j in range(len(R)):
                E = h2Data[j*110:(j+1)*110]
                self.h2Coef[j,:] = Coef(T, P, W, E)

            self.h2Coef = self.h2Coef.T

            print ("\tSaving H2O-pH2 Lebedev expansion coefficients to file 'H2O-pH2.leb'")

            np.savetxt("H2O-pH2.leb", self.h2Coef)

            self.h2PolyFit = np.zeros(((self.jMax+1)**2, self.args.polyDeg+1))

            for j in range((self.jMax+1)**2):
                self.h2PolyFit[j,:] = np.polyfit(R, self.h2Coef[j,:], self.args.polyDeg)



            if len(self.hAtoms) > 0:

                print ("\n\tFitting the Ab initio datapoints to a Lebedev Quadrature Scheme for H2O-H")
                hData = np.loadtxt(self.path + "PES/H2O-H").transpose()[self.args.state]

                self.hCoef = np.zeros((len(R), (self.jMax+1)**2))

                for j in range(len(R)):
                    E = hData[j*110:(j+1)*110]
                    self.hCoef[j,:] = Coef(T, P, W, E)

                self.hCoef = self.hCoef.T

                print ("\tSaving H2O-H Lebedev expansion coefficients to file 'H2O-pH2.leb'")

                np.savetxt("H2O-H.leb", self.hCoef)

                self.hPolyFit = np.zeros(((self.jMax+1)**2, self.args.polyDeg+1))

                for j in range((self.jMax+1)**2):
                    self.hPolyFit[j,:] = np.polyfit(R, self.hCoef[j,:], self.args.polyDeg)


    def GaussianSmearing(self):
        '''Function used to apply Gaussian smearing to the zero-point vibrational motion
            of the pH2 and/or H-atoms.

            Variables:
                self.pH2AbscVal     - Abscissae values for pH2
                self.pH2AbscVal     - Abscissae values for H
                self.pH2AbscWeight  - Abscissae weights for pH2
                self.pH2AbscWeight  - Abscissae weights for H

        '''


        print ("\n\tGenerating weights and abscissae for Gaussian quadrature (smearing) for:")
        print ("\tA pH2 molecule with a sigma value of", self.args.sigH2)

        self.pH2AbscVal, self.pH2AbscWeight = np.polynomial.legendre.leggauss(self.args.smearH2)

        if len(self.hAtoms) > 0:
            print ("\tAn H atom      with a sigma value of", self.args.sigH)
            self.hAbscVal, self.hAbscWeight = np.polynomial.legendre.leggauss(self.args.smearH)
        print ()


    def GetEnergyValues(self):
        '''Functions used to calculate the interaction energy for all points of the simulation

            Functions:
                Rotate      - Generate a rotation matrix based on given Euler angles
                IE_python   - Python-based approach to calculate interaction energies

            Variables::
                totalMol            - Total number of molecules in simulation
                totalAng            - Total number of rotations in simulation
                timeCheck           - boolean used to print predicted time for calculation of all interaction energies
                Rot                 - Rotation matrix
                E                   - Energy value
                paraH_rot           - Rotated matrix of pH2 positions
                hAtoms_rot          - Rotated matrix of H positions
                self.angleArray     - Matrix of all Euler angle values and weights and energies
        
        '''

        def Rotate(A, B, G):
            '''Function to calculate rotation matrix based on given Euler angles'''
            cosA = np.cos(A)
            cosB = np.cos(B)
            cosG = np.cos(G)

            sinA = np.sin(A)
            sinB = np.sin(B)
            sinG = np.sin(G)

            Rot = np.array([[ cosA * cosB * cosG - sinA * sinG, -cosA * cosB * sinG - sinA * cosG, cosA * sinB],
                            [ sinA * cosB * cosG + cosA * sinG, -sinA * cosB * sinG + cosA * cosG, sinA * sinB],
                            [-sinB * cosG,                       sinB * sinG,                      cosB       ]])


            return Rot

        def IE_python(pos, FIT, absc, weight, quad, sig, cutoff):
            '''Python-based approach to calculate interaction energies'''

            def Leb(C, theta, phi):
                '''Function that uses Lebedev quadrature to calculate interaction energies

                    Variables:
                        Sphere_harm         - Spherical Harmonic function
                        jx, jy, jz          - Quadrature values for Guassian smearing
                        weight_tot          - Total weight for the given Gaussian quadrature points
                        move                - Amount to shift position of the point
                        new_pos             - New position of point
                        R                   - Distance of point
                        Coef                - Lebedev expansion coefficients
                        T                   - Theta value of point
                        P                   - Phi value of point
                        E                   - Energy value
                '''
                count = 0
                jmax  = 8

                V_e = 0.0
                for j in range(0, jmax+1):
                    for m in range(-j, j+1):
                        if abs(C[count]) > cutoff:
                            Sphere_harm = scipy.special.sph_harm(m, j, phi, theta).real
                            V_e += Sphere_harm * C[count]
                        else:
                            pass
                        count += 1

                return V_e

            E = 0.

            for jx in range(quad):
                for jy in range(quad):
                    for jz in range(quad):

                        weight_tot = (weight[jx] * weight[jy] * weight[jz])

                        move = sig * np.array(([absc[jx], absc[jy], absc[jz]]))

                        new_pos = pos + move

                        R = np.linalg.norm(new_pos)
                        new_pos = new_pos / R

                        Coef = np.zeros((81))

                        for j in range (81):
                            Coef[j] = np.poly1d(FIT[j,:])(R)

                        T = np.arccos(new_pos[2]/R)
                        P = np.arctan2(new_pos[1], new_pos[0])

                        if P < 0.:
                            P += 2 * np.pi

                        E_ = Leb(Coef, T, P)

                        E += E_ * weight_tot

            return E / sum(weight)**3

        def PrintTime(val, string):
            '''Function used to print time value'''
            minutes, seconds = divmod(val, 60)
            hours,   minutes = divmod(minutes, 60)
            days,    hours   = divmod(hours, 24)
            print ("\t" + string + " - ", "%d days, %d hours, %d minutes, %d seconds" %(days, hours, minutes, seconds), "\n")
            return


        print ()
        
        totalMol = self.paraH.shape[0] + self.hAtoms.shape[0]
        totalAng = self.angleArray.shape[0]

        angCounter = np.linspace(0, totalAng, 10)
        aCount = 0

        timeCheck = False

        start = time.time()

        for angle in range(self.angleArray.shape[0]):
            alpha  = self.angleArray[angle, 0]
            beta   = self.angleArray[angle, 1]
            gamma  = self.angleArray[angle, 2]

            E = 0.

            Rot = Rotate(alpha, beta, gamma)

            paraH_rot = np.matmul(self.paraH, Rot)

            count = 0

            for item in paraH_rot:

                if self.args.method == 'python':

                    E += IE_python(item, 
                                   self.h2PolyFit, 
                                   self.pH2AbscVal,  
                                   self.pH2AbscWeight, 
                                   self.args.smearH2, 
                                   self.args.sigH2, 
                                   self.args.lebThresh)
                else:
                    #TODO
                        # Add fortran code
                    E += 0

                if angle == 0 and timeCheck == False:
                    stop = time.time()
                    val = round((stop - start) * totalMol * totalAng, 3)

                    PrintTime(val, "Estimated time for Calculating all Interaction Energies")

                    timeCheck = True

                if self.args.time == True:
                    exit()

                count += 1

            if len(self.hAtoms) != 0:
                hAtoms_rot = np.matmul(self.hAtoms, Rot)

                for item in hAtoms_rot:
                    if self.args.method == 'python':
                        E += IE_python(item,
                                       self.hPolyFit,
                                       self.hAbscVal,
                                       self.hAbscWeight,
                                       self.args.smearH,
                                       self.args.sigH,
                                       self.args.lebThresh)

                    else:
                        #TODO
                            # Add fortran code
                        E += 0

                    count += 1

            self.angleArray[angle, 6] = E

            if angle > angCounter[aCount]:
                print ("\t%4d out of %4d Interaction Energies Calculated (%4d" % (angle*totalMol+count, totalAng*totalMol, (aCount+1)*10) + "%)")
                aCount += 1

        if self.angleArray.shape == (5,1):
            self.angleArray = np.array(([self.angleArray]))

        print ("\n\n")
        print ("\tCalculation of interaction energies completed.")
        
        np.savetxt("Energy.en", self.angleArray)

        stop = time.time()
        val = round((stop - start), 3)
        PrintTime(val, "Total time")
        print ()

