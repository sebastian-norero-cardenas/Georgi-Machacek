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
vDelta = 33.2992327030555  # [GeV]
M1 = 55.1761999832967  # [GeV]
M2 = 0.0000000038214549  # [GeV]0.2463699
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
print('sA: ' + str(sA), 'vDelta: ' + str(vDelta), 'M1: ' + str(M1), 'M2: ' + str(M2), 'm3: ' + str(m3), 'm5: ' + str(m5), 'mH: ' + str(mH), sep='\n')


#######################
# NUMERICAL INVERSION #
#######################

def equations(vars):
    lam1, lam2, lam3, lam4, lam5 = vars
    M11sq = 8 * lam1 * vphi ** 2
    M22sq = (M1 * vphi ** 2) / (4 * vDelta) - 6 * M2 * vDelta + 8 * (lam3 + 3 * lam4) * vDelta ** 2
    M12sq = (np.sqrt(3) * vphi / 2) * (-M1 + 4 * (2 * lam2 - lam5) * vDelta)
    eq1 = (M1 / (4 * vDelta) + lam5 / 2) * v ** 2 - m3 ** 2
    eq2 = (M1 * vphi ** 2) / (4 * vDelta) + 12 * vDelta * M2 + 8 * vDelta ** 2 * lam3 + (3/2) * vphi ** 2 * lam5 - m5 ** 2
    eq3 = M11sq + M22sq - np.sqrt((M11sq - M22sq) ** 2 + 4 * (M12sq) ** 2) - 2 * mh ** 2
    eq4 = M11sq + M22sq + np.sqrt((M11sq - M22sq) ** 2 + 4 * (M12sq) ** 2) - 2 * mH ** 2
    eq5 = (1 / (8 * vphi ** 2)) * (mh ** 2 + (M12sq) ** 2 / (M22sq - mh ** 2)) - lam1
    return [eq1, eq2, eq3, eq4, eq5]

lam1_num, lam2_num, lam3_num, lam4_num, lam5_num = fsolve(equations, (0.5, 0.5, 0.5, 0.5, 0.5))

print('---------------------------------------------------')
print('--- NUMERICAL INVERTION ---------------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lam1_num), 'lam2: ' + str(lam2_num), 'lam3: ' + str(lam3_num), 'lam4: ' + str(lam4_num), 'lam5: ' + str(lam5_num), sep='\n')

#################################################
# INVERSION BY https://arxiv.org/abs/2207.00142 #
#################################################

# In this paper, v_xi = v_Delta and v_chi = sqrt(2) * vDelta

tB = 2 * np.sqrt(2) * vDelta / vphi
sB = tB / np.sqrt(1 + tB ** 2)
cB = np.sqrt(1 - sB ** 2)

def lambda1_2207():
    numerator = mh ** 2 * cA ** 2 + mH ** 2 * sA ** 2
    denominator = 8 * v ** 2 * cB ** 2
    return numerator / denominator

def lambda2_2207():
    numerator = 2 * np.sqrt(6) * (mh ** 2 - mH ** 2) * sA * cA + 12 * m3 ** 2 * cB * cB - 3 * np.sqrt(2) * v * cB * M1
    denominator = 12 * v ** 2 * sB * cB
    return numerator / denominator

def lambda3_2207():
    numerator = m5 ** 2 - 3 * m3 ** 2 * cB ** 2 + np.sqrt(2) * v * M1 * cB ** 2 / sB - 3 * np.sqrt(2) * v * sB * M2
    denominator = v ** 2 * sB ** 2
    return numerator / denominator

def lambda4_2207():
    numerator = 2 * mH ** 2 * cB ** 2 + 2 * mh ** 2 * sA ** 2 - 2 * m5 ** 2 + 6 * cB ** 2 * m3 ** 2 - 3 * np.sqrt(2) * v * M1 * cB ** 3 / sB ** 2 + 9 * np.sqrt(2) * v * M2 * sB
    denominator = 6 * v ** 2 * sB ** 2
    return numerator / denominator

def lambda5_2207():
    numerator = 2 * m3 ** 2 * v * sB - np.sqrt(2) * M1 * v ** 2
    denominator = v ** 3 * sB
    return numerator / denominator

print('---------------------------------------------------')
print('--- arXiv:2207.00142 ------------------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lambda1_2207()), 'lam2: ' + str(lambda2_2207()), 'lam3: ' + str(lambda3_2207()), 'lam4: ' + str(lambda4_2207()), 'lam5: ' + str(lambda5_2207()), sep='\n')

# let's try to recover the masses and alpha-mixing angle
M11sq_2207 = 8 * lambda1_2207() * cB ** 2 * v ** 2
M22sq_2207 = M1 * cB ** 2 * v / np.sqrt(2) - 3 * M2 * sB * v / np.sqrt(2) + (lambda3_2207() + 3 * lambda4_2207()) * sB ** 2 * v ** 2
M12sq_2207 = (np.sqrt(3) * cB * v / 2) * (-M1 + np.sqrt(2) * (2 * lambda2_2207() - lambda5_2207()) * sB * v)

m3_2207 = np.sqrt((M1 / (np.sqrt(2) * sB * v) + lambda5_2207() / 2) * v ** 2)
m5_2207 = np.sqrt((M1 * cB ** 2 * v) / (np.sqrt(2) * sB) + 6 * M2 * sB * v / np.sqrt(2) + 3 * lambda5_2207() * cB ** 2 * v ** 2 / 2 + lambda3_2207() * sB ** 2 * v ** 2)
mh_2207 = (1 / np.sqrt(2)) * np.sqrt(np.abs(M11sq_2207 + M22sq_2207 - np.sqrt((M11sq_2207 - M22sq_2207) ** 2 + 4 * M12sq_2207 ** 2)))
mH_2207 = (1 / np.sqrt(2)) * np.sqrt(M11sq_2207 + M22sq_2207 + np.sqrt((M11sq_2207 - M22sq_2207) ** 2 + 4 * M12sq_2207 ** 2))
sA_2207 = np.sin(np.arctan(2 * M12sq_2207 / (M22sq_2207 - M11sq_2207)) / 2)

print('--- which imply: ')
print('m3: ' + str(m3_2207), 'm5: ' + str(m5_2207), 'mh: ' + str(mh_2207), 'mH: ' + str(mH_2207), 'sA: ' + str(sA_2207), sep='\n')

###################################################
# INVERSION BY https://arxiv.org/abs/1510.06297v2 #
###################################################

tH = 2 * np.sqrt(2) * vDelta / vphi
sH = tB / np.sqrt(1 + tB ** 2)
cH = np.sqrt(1 - sB ** 2)
s2H = 2 * sH * cH

M1sq_1510 = - (v * (-M1)) / (np.sqrt(2) * sH)
M2sq_1510 = - 3 * np.sqrt(2) * sH *v * (-M2)

def lambda1_1510():
    numerator = mh ** 2 * cA ** 2 + mH ** 2 * sA ** 2
    denominator = 8 * cH ** 2 * v ** 2
    return numerator / denominator


def lambda2_1510():
    numerator = np.sqrt(6) * s2A * (mh ** 2 - mH ** 2) / 2 + 3 * sH * cH * (2 * m3 ** 2 - M1sq_1510)
    denominator = 6 * v ** 2 * sH * cH
    return numerator / denominator


def lambda3_1510():
    numerator = cH ** 2 * (2 * M1sq_1510 - 3 * m3 ** 2) + m5 ** 2 - M2sq_1510
    denominator = v ** 2 * sH ** 3
    return numerator / denominator


def lambda4_1510():
    numerator = 2 * mH ** 2 * cA ** 2 + 2 * mh ** 2 * sA ** 2 + 3 * M2sq_1510 - 2 * m5 ** 2 + 6 * cH ** 2 * (m3 ** 2 - M1sq_1510)
    denominator = 6 * v ** 2 * sH ** 2
    return numerator / denominator


def lambda5_1510():
    numerator = 2 * (M1sq_1510 - m3 ** 2)
    denominator = v ** 2
    return - numerator / denominator

print('---------------------------------------------------')
print('--- arXiv:1510.06297 ------------------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lambda1_1510()), 'lam2: ' + str(lambda2_1510()), 'lam3: ' + str(lambda3_1510()), 'lam4: ' + str(lambda4_1510()), 'lam5: ' + str(lambda5_1510()), sep='\n')
print('---------------------------------------------------')
