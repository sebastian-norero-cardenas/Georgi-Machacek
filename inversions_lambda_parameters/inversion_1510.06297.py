###############
# INFORMATION #
###############

""" PAPER: https://arxiv.org/abs/1510.06297v2 """

###########
# MODULES #
###########

from scipy.optimize import fsolve
from math import exp
import numpy as np

##############
# PARAMETERS #
##############

# CONSTANTS
v = 246.22  # [GeV]
mh = 125.280865  # [GeV]

# INDEPENDENT PARAMETERS
sA = 0.502369549610856
vDelta = 33.2992327030555  # [GeV]
mu1 = -55.1761999832967  # [GeV]
mu2 = -0.0000000038214549  # [GeV]
m3 = 133.14803  # [GeV]
m5 = 66.996471  # [GeV]
mH = 161.253739  # [GeV]

# (SOME) DEPENDENT PARAMETERS
vphi = np.sqrt(v ** 2 - 8 * vDelta ** 2)
cA = np.sqrt(1 - sA ** 2)
s2A = 2 * sA * cA

####################
# INPUT PARAMETERS #
####################

print('---------------------------------------------------')
print('--- INPUT PARAMETERS ------------------------------')
print('---------------------------------------------------')
print('sA: ' + str(sA), 'vDelta: ' + str(vDelta), 'mu_1: ' + str(mu1), 'mu_2: ' + str(mu2), 'm3: ' + str(m3), 'm5: ' + str(m5), 'mh: ' + str(mh), 'mH: ' + str(mH), sep='\n')

##################################
# CONSTRUCTING LAMBDA PARAMETERS #
##################################

tH = 2 * np.sqrt(2) * vDelta / vphi
sH = tH / np.sqrt(1 + tH ** 2)
cH = 1 / np.sqrt(1 + tH ** 2)
s2H = 2 * sH * cH

M1sq = - (v * mu1) / (np.sqrt(2) * sH)
M2sq = - 3 * np.sqrt(2) * sH * v * mu2

def lambda1():
    numerator = mh ** 2 * cA ** 2 + mH ** 2 * sA ** 2
    denominator = 8 * v ** 2 * cH ** 2
    return numerator / denominator


def lambda2():
    numerator = 2 * mH ** 2 * cA ** 2 + 2 * mh ** 2 * sA ** 2 + 3 * M2sq - 2 * m5 ** 2 + 6 * cH ** 2 * (m3 ** 2 - M1sq)
    denominator = 6 * v ** 2 * sH ** 2
    return numerator / denominator


def lambda3():
    numerator = cH ** 2 * (2 * M1sq - 3 * m3 ** 2) + m5 ** 2 - M2sq
    denominator = v ** 2 * sH ** 3
    return numerator / denominator


def lambda4():
    numerator = np.sqrt(6) * s2A * (mh ** 2 - mH ** 2) / 2 + 3 * sH * cH * (2 * m3 ** 2 - M1sq)
    denominator = 6 * v ** 2 * sH * cH
    return numerator / denominator


def lambda5():
    numerator = 2 * (M1sq - m3 ** 2)
    denominator = v ** 2
    return numerator / denominator

print('---------------------------------------------------')
print('--- arXiv:1510.06297 ------------------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lambda1()), 'lam2: ' + str(lambda2()), 'lam3: ' + str(lambda3()), 'lam4: ' + str(lambda4()), 'lam5: ' + str(lambda5()), sep='\n')
print('---------------------------------------------------')


# let's try to recover the masses and alpha-mixing angle
M11sq = 8 * cH ** 2 * lambda1() * v ** 2
M22sq = sH ** 2 * (3 * lambda2() + lambda3()) * v ** 2 + cH ** 2 * M1sq - M2sq / 2
M12sq = np.sqrt(3 / 2) * sH * cH * ((2 * lambda4() + lambda5()) * v ** 2 - M1sq)

tan2A = (2 * M12sq) / (M22sq - M11sq)
def sinA():
    if 0 < sA < 1 / np.sqrt(2):
        return np.sqrt((np.sqrt(tan2A ** 2 + 1) - 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if 1 / np.sqrt(2) < sA < 1:
        return np.sqrt((np.sqrt(tan2A ** 2 + 1) + 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if -1 / np.sqrt(2) < sA < 0:
        return - np.sqrt((np.sqrt(tan2A ** 2 + 1) - 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if -1 < sA < -1 / np.sqrt(2):
        return - np.sqrt((np.sqrt(tan2A ** 2 + 1) + 1) / (2 * np.sqrt(tan2A ** 2 + 1)))

mass3 = np.sqrt(- lambda5() * v ** 2 / 2 + M1sq)
mass5 = np.sqrt((sH ** 2 * lambda3() - 3 * cH ** 2 * lambda5() / 2) * v ** 2 + cH ** 2 * M1sq + M2sq)
massh = np.sqrt(M11sq * cA ** 2 + M22sq * sA ** 2 - 2 * M12sq * sA * cA)
massH = np.sqrt(M11sq * sA ** 2 + M22sq * cA ** 2 + 2 * M12sq * sA * cA)
def sinA():
    if 0 < sA < 1 / np.sqrt(2):
        return np.sqrt((np.sqrt(tan2A ** 2 + 1) - 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if 1 / np.sqrt(2) < sA < 1:
        return np.sqrt((np.sqrt(tan2A ** 2 + 1) + 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if -1 / np.sqrt(2) < sA < 0:
        return - np.sqrt((np.sqrt(tan2A ** 2 + 1) - 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if -1 < sA < -1 / np.sqrt(2):
        return - np.sqrt((np.sqrt(tan2A ** 2 + 1) + 1) / (2 * np.sqrt(tan2A ** 2 + 1)))

print('[..] and trying to recover the input parameters: ')
print('sA: ' + str(sinA()), 'm3: ' + str(mass3), 'm5: ' + str(mass5), 'm_h: ' + str(massh), 'm_H: ' + str(massH), sep='\n')

