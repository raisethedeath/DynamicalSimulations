import numpy as np
import sys

def HCP(SHELL, Dis):

    NUM = SHELL * 2 + 1
    MAX = NUM ** 3

    if SHELL == 1:
        H2Pos = np.zeros((13,3))
    elif SHELL == 2:
        H2Pos = np.zeros((57,3))
    elif SHELL == 3:
        H2Pos = np.zeros((157,3))
    elif SHELL >= 4 or SHELL == 0:
        print ("This code can only handle shells between 1 and up to 3")
        sys.exit()

    Pos   = np.zeros((3))
    Atoms = np.zeros((MAX, 3))
    Base  = np.zeros((NUM**2, 3))

    X = 1/2
    Y = np.sqrt(3) / 2
    Z = np.sqrt(33) / 6

    count = 0

    for j in range(0, NUM):
        Pos[0] = j * X
        Pos[1] = j * Y
        for k in range(0, NUM):
            Base[count] = Pos
            Pos[0] +=  1.0
            count += 1

    Mid = int((MAX - 1) / 2)

    for i in range(0, NUM):
        if i % 2 == 0:
            for j in range(len(Base)):
                Atoms[i * NUM**2 + j][0] = Base[j][0]
                Atoms[i * NUM**2 + j][1] = Base[j][1]
                Atoms[i * NUM**2 + j][2] = Base[j][2] + i * np.sqrt(2.0/3.0)

        else:
            for j in range(len(Base)):
                Atoms[i * NUM**2 + j][0] = Base[j][0]  
                Atoms[i * NUM**2 + j][1] = Base[j][1] + np.sqrt(3)/3.0
                Atoms[i * NUM**2 + j][2] = Base[j][2] + i * np.sqrt(2.0/3.0)


    Atoms = Atoms * Dis            

    Atoms -= Atoms[Mid]

    count = 0 

    for i in Atoms:
        R = np.linalg.norm(i)
        if R <= 1.01 * Dis * SHELL and R != 0:
            H2Pos[count] = i

            count += 1

    return H2Pos[:-1]


def FCC(SHELL, Dis):

    print ("HERE")

    return np.loadtxt("FCC.xyz")
