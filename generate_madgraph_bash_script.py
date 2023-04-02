#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + Type-II Seesaw Term
INFO: Intended    -> Generate a bash script with MadGraph instructions
                     to calculate cross sections.
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> December 31, 2022.
"""

###########
# MODULES #
###########

import numpy as np

#############
# CONSTANTS #
#############

mh = 125.18  # mass of the 125-Higgs in [GeV]
v = 246.22  # vev of the 125-Higgs in [GeV]

##################
# ABOUT THE RUNS #
##################

sA = 0
vDelta = 1e-4  # [GeV]
M1 = 1e-4  # [GeV]
M2 = 0  # [GeV]

m5_max = 1000  # [GeV]
m5_min = 45  # [GeV]

ebeam1 = 7000  # [GeV]
ebeam2 = 7000  # [GeV]

N = 30  # Number of runs

massenergy = 'm3Em5_14TeV'
process = 'ppTOWpmH3pmTOH5ppmmH3mp'

#############
# FUNCTIONS #
#############

def tanH(vDelta):
    sH = (2 * np.sqrt(2) * vDelta) / v
    cH = np.sqrt(1 - sH ** 2)
    return sH / cH


def vphi(vDelta):
    vD = vDelta
    return np.sqrt(v ** 2 - 8 * vD ** 2)


def lambda2_inv(sA, vDelta, M1, m3, mh, mH):
    vD = vDelta
    vp = vphi(vDelta)
    s2A = np.sin(2 * np.arcsin(sA))
    numerator = np.sqrt(3) * v ** 2 * s2A * (mH ** 2 - mh ** 2) - 3 * M1 * vp * v ** 2 + 24 * m3 ** 2 * vp * vD
    denominator = 24 * vp * vD * v ** 2
    return numerator / denominator


def lambda3_inv(vDelta, M1, M2, m3, m5):
    vp = vphi(vDelta)
    vD = vDelta
    numerator = m5 ** 2 + (M1 * vp ** 2) / (2 * vD) - 12 * M2 * vD - 3 * m3 ** 2 * vp ** 2 / v ** 2
    denominator = 8 * vD ** 2
    return numerator / denominator


def lambda4_inv(sA, vDelta, M1, M2, m3, m5, mh, mH):
    vp = vphi(vDelta)
    vD = vDelta
    cA = np.sqrt(1 - sA ** 2)
    numerator = mh ** 2 * sA ** 2 + mH ** 2 * cA ** 2 - m5 ** 2 + 3 * m3 ** 2 * vp ** 2 / v ** 2 - 3 * M1 * vp ** 2 / (4 * vD) + 18 * M2 * vD
    denominator = 24 * vD ** 2
    return numerator / denominator


def lambda5_inv(vDelta, M1, m3):
    vD = vDelta
    numerator = 8 * vD * m3 ** 2 - 2 * M1 * v ** 2
    denominator = 4 * vD * v ** 2
    return numerator / denominator


def mass5(i):
    return ((m5_max - m5_min) / (N - 1)) * i + m5_min


def mass3(i):
    return mass5(i)


def massH(i):
    return 200

######################
# WRITING THE SCRIPT #
######################

with open(f'MG5_bash_script_{massenergy}_{process}.sh', 'w') as f:

    f.write('#! /bin/bash' + '\n' * 2)

    for i in range(N):
        
        m5 = mass5(i)
        m3 = mass3(i)
        mH = massH(i)

        lam2 = lambda2_inv(sA, vDelta, M1, m3, mh, mH)
        lam3 = lambda3_inv(vDelta, M1, M2, m3, m5)
        lam4 = lambda4_inv(sA, vDelta, M1, M2, m3, m5, mh, mH)
        lam5 = lambda5_inv(vDelta, M1, m3)

        f.write('# --run = ' + str(i+1) + '\n'
                'python3.8 /home/sebastian/Projects/Thesis/MG5_aMC_v3_4_0/bin/mg5_aMC <<EOD' + '\n'
                f'  launch Processes/{massenergy}/{process}' + '\n'
                '   0' + '\n'
                '   set tanth ' + str(tanH(vDelta)) + '\n'
                '   set m1coeff ' + str(M1) + '\n'
                '   set m2coeff ' + str(M2) + '\n'
                '   set lam2 ' + str(lam2) + '\n'
                '   set lam3 ' + str(lam3) + '\n'
                '   set lam4 ' + str(lam4) + '\n'
                '   set lam5 ' + str(lam5) + '\n'
                '   set ebeam1 ' + str(ebeam1) + '\n'
                '   set ebeam2 ' + str(ebeam2) + '\n'
                '   0' + '\n'
                '   quit()' + '\n'
                'EOD' + '\n' * 2)


print('DONE!')