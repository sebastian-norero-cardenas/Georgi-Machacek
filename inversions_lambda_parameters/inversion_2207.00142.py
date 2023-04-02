###############
# INFORMATION #
###############

""" PAPER: https://arxiv.org/abs/2207.00142 """

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
sA = -0.502369549610856
vxi = 33.2992327030555  # [GeV]
mu1 = 55.1761999832967  # [GeV]
mu2 = 0.0000000038214549  # [GeV]
m3 = 133.14803  # [GeV]
m5 = 66.996471  # [GeV]
meta = 161.253739  # [GeV]

# (SOME) DEPENDENT PARAMETERS
vphi = np.sqrt(v ** 2 - 8 * vxi ** 2)
cA = np.sqrt(1 - sA ** 2)
s2A = 2 * sA * cA

# In this paper, v_xi = v_Delta and v_chi = sqrt(2) * vDelta

tB = 2 * np.sqrt(2) * vxi / vphi
sB = tB / np.sqrt(1 + tB ** 2)
cB = 1 / np.sqrt(1 + tB ** 2)

####################
# INPUT PARAMETERS #
####################

print('---------------------------------------------------')
print('--- INPUT PARAMETERS ------------------------------')
print('---------------------------------------------------')
print('sA: ' + str(sA), 'm3: ' + str(m3), 'm5: ' + str(m5), 'm_eta: ' + str(meta), sep='\n')


##################################
# CONSTRUCTING LAMBDA PARAMETERS #
##################################

def lambda1():
    numerator = mh ** 2 * cA ** 2 + meta ** 2 * sA ** 2
    denominator = 8 * v ** 2 * cB ** 2
    return numerator / denominator

def lambda2():
    numerator = 2 * np.sqrt(6) * (mh ** 2 - meta ** 2) * sA * cA + 12 * m3 ** 2 * cB  * cB - 3 * np.sqrt(2) * v * cB * mu1
    denominator = 12 * v ** 2 * sB * cB
    return numerator / denominator

def lambda3():
    numerator = m5 ** 2 - 3 * m3 ** 2 * cB ** 2 + np.sqrt(2) * v * mu1 * cB ** 2 / sB - 3 * np.sqrt(2) * v * sB * mu2
    denominator = v ** 2 * sB ** 2
    return numerator / denominator

def lambda4():
    numerator = 2 * meta ** 2 * cA ** 2 + 2 * mh ** 2 * sA ** 2 - 2 * m5 ** 2 + 6 * cB ** 2 * m3 ** 2 - 3 * np.sqrt(2) * v * mu1 * cB ** 3 / sB ** 2 + 9 * np.sqrt(2) * v * mu2 * sB
    denominator = 6 * v ** 2 * sB ** 2
    return numerator / denominator

def lambda5():
    numerator = 2 * m3 ** 2 * v * sB - np.sqrt(2) * mu1 * v ** 2
    denominator = v ** 3 * sB
    return numerator / denominator

print('---------------------------------------------------')
print('--- arXiv:2207.00142 ------------------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lambda1()), 'lam2: ' + str(lambda2()), 'lam3: ' + str(lambda3()), 'lam4: ' + str(lambda4()), 'lam5: ' + str(lambda5()), sep='\n')

# let's try to recover the masses and alpha-mixing angle
M11sq = 8 * lambda1() * cB ** 2 * v ** 2
M22sq = (mu1 * cB ** 2 * v) / (np.sqrt(2) * sB) - 3 * mu2 * sB * v / np.sqrt(2) + (lambda3() + 3 * lambda4()) * sB ** 2 * v ** 2
M12sq = (np.sqrt(3) * cB * v / 2) * (-mu1 + np.sqrt(2) * (2 * lambda2() - lambda5()) * sB * v)

tan2A = 2 * M12sq / (M22sq - M11sq)
def sinA():
    if 0 < sA < 1 / np.sqrt(2):
        return np.sqrt((np.sqrt(tan2A ** 2 + 1) - 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if 1 / np.sqrt(2) < sA < 1:
        return np.sqrt((np.sqrt(tan2A ** 2 + 1) + 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if -1 / np.sqrt(2) < sA < 0:
        return - np.sqrt((np.sqrt(tan2A ** 2 + 1) - 1) / (2 * np.sqrt(tan2A ** 2 + 1)))
    if -1 < sA < -1 / np.sqrt(2):
        return - np.sqrt((np.sqrt(tan2A ** 2 + 1) + 1) / (2 * np.sqrt(tan2A ** 2 + 1)))

mass3 = np.sqrt((mu1 / (np.sqrt(2) * sB * v) + lambda5() / 2) * v ** 2)
mass5 = np.sqrt((mu1 * cB ** 2 * v) / (np.sqrt(2) * sB) + 6 * mu2 * sB * v / np.sqrt(2) + 3 * lambda5() * cB ** 2 * v ** 2 / 2 + lambda3() * sB ** 2 * v ** 2)
# massh = (1 / np.sqrt(2)) * np.sqrt(M11sq + M22sq - np.sqrt((M11sq - M22sq) ** 2 + 4 * M12sq ** 2))
masseta = (1 / np.sqrt(2)) * np.sqrt(M11sq + M22sq + np.sqrt((M11sq - M22sq) ** 2 + 4 * M12sq ** 2))

print('[..] and trying to recover the input parameters: ')
print('sA: ' + str(sinA()), 'm3: ' + str(mass3), 'm5: ' + str(mass5), 'm_eta: ' + str(masseta), sep='\n')
