###############
# INFORMATION #
###############

""" PAPER: https://arxiv.org/abs/2205.03896 """

###########
# MODULES #
###########

import numpy as np

##############
# PARAMETERS #
##############

# CONSTANTS
v = 246.22  # [GeV]
mh = 125.280865  # [GeV]

# INDEPENDENT PARAMETERS
sA = -0.502369549610856
v2 = 33.2992327030555  # [GeV]
M1 = 55.1761999832967  # [GeV]
M2 = 0.0000000038214549  # [GeV]
m3 = 133.14803  # [GeV]
m5 = 66.996471  # [GeV]
mH = 161.253739  # [GeV]

# (SOME) DEPENDENT PARAMETERS
v1 = np.sqrt(v ** 2 - 8 * v2 ** 2)
cA = np.sqrt(1 - sA ** 2)
s2A = 2 * sA * cA

tB = 2 * np.sqrt(2) * v2 / v1
sB = tB / np.sqrt(1 + tB ** 2)
cB = 1 / np.sqrt(1 + tB ** 2)

####################
# INPUT PARAMETERS #
####################

print('---------------------------------------------------')
print('--- INPUT PARAMETERS ------------------------------')
print('---------------------------------------------------')
print('sA: ' + str(sA), 'm3: ' + str(m3), 'm5: ' + str(m5), 'm_h: ' + str(mh), 'm_H: ' + str(mH), sep='\n')


##################################
# CONSTRUCTING LAMBDA PARAMETERS #
##################################

def lambda1():
    numerator = mh ** 2 * cA ** 2 + mH ** 2 * sA ** 2
    denominator = 8 * v1 ** 2
    return numerator / denominator

def lambda2():
    numerator = - (mH ** 2 - mh ** 2) * (1 / (np.sqrt(3) * v1)) * s2A - M1 + 8 * m3 ** 2 * v2 / v ** 2
    denominator = 8 * v2
    return numerator / denominator

def lambda3():
    numerator = m5 ** 2 - 3 * m3 ** 2 * v1 ** 2 / v ** 2 + (M1 * v1 ** 2) / (2 * v2) - 12 * M2 * v2
    denominator = 8 * v2 ** 2
    return numerator / denominator

def lambda4():
    numerator = mh ** 2 * sA ** 2 + mH ** 2 * cA ** 2 - m5 ** 2 + 3 * m3 ** 2 * v1 ** 2 / v ** 2 - (3 * M1 * v1 ** 2) / (4 * v2) + 18 * M2 * v2
    denominator = 24 * v2 ** 2
    return numerator / denominator

def lambda5():
    numerator = 2 * (4 * v2 * m3 ** 2 - M1 * v ** 2)
    denominator = 4 * v ** 2 * v2
    return numerator / denominator

print('---------------------------------------------------')
print('--- arXiv:2207.00142 ------------------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lambda1()), 'lam2: ' + str(lambda2()), 'lam3: ' + str(lambda3()), 'lam4: ' + str(lambda4()), 'lam5: ' + str(lambda5()), sep='\n')

# let's try to recover the masses and alpha-mixing angle
M11sq = 8 * lambda1() * v1 ** 2
M22sq = (M1 * v1 ** 2) / (4 * v2)
M12sq = (np.sqrt(3) / 2) * (-M1 + 4 * (2 * lambda2() - lambda5()) * v2) * v1

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

mass3 = np.sqrt((M1 / (4 * v2) + lambda5() / 2) * v ** 2)
mass5 = np.sqrt((M1 * v1 ** 2) / (4 * v2) + 12 * M2 * v2 + 3 * lambda5() * v1 ** 2 / 2 + 8 * lambda3() * v2 ** 2)
massh = (1 / np.sqrt(2)) * np.sqrt(M11sq + M22sq - np.sqrt((M11sq - M22sq) ** 2 + 4 * M12sq ** 2))
massH = (1 / np.sqrt(2)) * np.sqrt(M11sq + M22sq + np.sqrt((M11sq - M22sq) ** 2 + 4 * M12sq ** 2))

print('[..] and trying to recover the input parameters: ')
print('sA: ' + str(sinA()), 'm3: ' + str(mass3), 'm5: ' + str(mass5), 'm_h: ' + str(massh), 'm_H: ' + str(massH), sep='\n')
