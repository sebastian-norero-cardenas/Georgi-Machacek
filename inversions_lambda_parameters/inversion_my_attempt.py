###############
# INFORMATION #
###############

""" MY ATTEMPT """

###########
# MODULES #
###########

import numpy as np

##############
# PARAMETERS #
##############

# CONSTANTS
v = 246.22  # [GeV]

# INDEPENDENT PARAMETERS
sA = -0.0212657990268057
vDelta = 1.40005394559093  # [GeV]
M1 = 16.2673715453875  # [GeV]
M2 = -339.401615596681  # [GeV]
m3 = 419.7434837462  # [GeV]
m5 = 413.0458127676  # [GeV]
mh = 125.2696585732  # [GeV]
mH = 423.0590866122  # [GeV]

# (SOME) DEPENDENT PARAMETERS
vphi = np.sqrt(v ** 2 - 8 * vDelta ** 2)
cA = np.sqrt(1 - sA ** 2)
s2A = 2 * sA * cA

tH = 2 * np.sqrt(2) * vDelta / vphi
sH = tH / np.sqrt(1 + tH ** 2)
cH = 1 / np.sqrt(1 + tH ** 2)

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
    denominator = 8 * vphi ** 2
    return numerator / denominator

def lambda2():
    numerator = np.sqrt(3) * v ** 2 * s2A * (mH ** 2 - mh ** 2) - 3 * M1 * vphi * v ** 2 + 24 * m3 ** 2 * vphi * vDelta
    denominator = 24 * vphi * vDelta * v ** 2
    return numerator / denominator

def lambda3():
    numerator = m5 ** 2 + (M1 * vphi ** 2) / (2 * vDelta) - 12 * M2 * vDelta - 3 * m3 ** 2 * vphi ** 2 / v ** 2
    denominator = 8 * vDelta ** 2
    return numerator / denominator

def lambda4():
    numerator = mh ** 2 * sA ** 2 + mH ** 2 * cA ** 2 - m5 ** 2 + 3 * m3 ** 2 * vphi ** 2 / v ** 2 - 3 * M1 * vphi ** 2 / (4 * vDelta) + 18 * M2 * vDelta
    denominator = 24 * vDelta ** 2
    return numerator / denominator

def lambda5():
    numerator = 8 * vDelta * m3 ** 2 - 2 * M1 * v ** 2
    denominator = 4 * vDelta * v ** 2
    return numerator / denominator

print('---------------------------------------------------')
print('--- ACCORDING TO MY FORMULAS ----------------------')
print('---------------------------------------------------')
print('lam1: ' + str(lambda1()), 'lam2: ' + str(lambda2()), 'lam3: ' + str(lambda3()), 'lam4: ' + str(lambda4()), 'lam5: ' + str(lambda5()), sep='\n')

# let's try to recover the masses and alpha-mixing angle
M11sq = 8 * lambda1() * cH ** 2 * v ** 2
M22sq = (M1 * cH ** 2 * v) / (np.sqrt(2) * sH) - 3 * M2 * sH * v / np.sqrt(2) + (lambda3() + 3 * lambda4()) * sH ** 2 * v ** 2
M12sq = (np.sqrt(3) * cH * v / 2) * (-M1 + np.sqrt(2) * (2 * lambda2() - lambda5()) * sH * v)

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

mass3 = np.sqrt((M1 / (np.sqrt(2) * sH * v) + lambda5() / 2) * v ** 2)
mass5 = np.sqrt((M1 * cH ** 2 * v) / (np.sqrt(2) * sH) + 6 * M2 * sH * v / np.sqrt(2) + 3 * lambda5() * cH ** 2 * v ** 2 / 2 + lambda3() * sH ** 2 * v ** 2)
massh = (1 / np.sqrt(2)) * np.sqrt(M11sq + M22sq - np.sqrt((M11sq - M22sq) ** 2 + 4 * M12sq ** 2))
massH = (1 / np.sqrt(2)) * np.sqrt(M11sq + M22sq + np.sqrt((M11sq - M22sq) ** 2 + 4 * M12sq ** 2))

print('[..] and trying to recover the input parameters: ')
print('sA: ' + str(sinA()), 'm3: ' + str(mass3), 'm5: ' + str(mass5), 'm_h: ' + str(massh), 'm_H: ' + str(massH), sep='\n')
