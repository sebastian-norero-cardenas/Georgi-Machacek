# -*- coding: utf-8 -*-
#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Print a benchmark point
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Version     -> 1.0.0
INFO: Last update -> November 08, 2022.
"""

###########
# MODULES #
###########
import pandas as pd
import numpy as np
import math
import cmath
import random as rd

###############
# IMPORT DATA #
###############

# Directory to the .csv with random points
DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterData20000_GMSEESAW.csv'

df = pd.read_csv(DirRandomScan)

#############
# CONSTANTS #
#############

# Independent physical constants
# Ref.: R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
alphaEW = 1 / 127.9  # EW coupling constant at energy = mZ
GF = 1.166378 * 1E-5  # Fermi constant [GeV]

# Masses
mZ = 91.1876  # Z boson mass [GeV]

me = 5.11 * 1E-4  # Electron mass [GeV]
mmu = 0.10566  # Muon mass [GeV]
mtau = 1.777  # Tau mass [GeV]

mup = 2.55 * 1E-3  # Up quark mass [GeV]
mcharm = 1.42  # Charm quark mass [GeV]
mtop = 172  # Top quark mass [GeV]
mdown = 5.04 * 1E-3  # Down quark mass [GeV]
mstrange = 0.101  # Strange quark mass [GeV]
mbottom = 4.7  # Bottom quark mass [GeV]

# Dependent parameters
v = 1 / np.sqrt(np.sqrt(2) * GF)  # SM Higgs vev ~ 246.22 [GeV]
mW = np.sqrt(mZ ** 2 / 2 + np.sqrt((mZ ** 2 / 2) ** 2 - (np.pi / np.sqrt(2)) * (alphaEW / GF) * mZ ** 2))  # W boson mass

# Internal parameters
r1 = 6  # round to r1 decimals
r2 = 16  # round to r2 decimals (used for neutrino yukawas and other very small parameters)

################################
# INDEPENDENT PARAMETER RANGES #
################################

log_vDelta = 0.937618

M1 = 100.0
M2 = 500.0

lam2 = 0.1
lam3 = 0.1
lam4 = 0.1
lam5 = 3.1

# NEUTRINO PHYSICS
# All values assume a normal neutrino mass ordering

# MASS SPLITTINGS
# Solar mass splitting
Dm21 = 8.14 * 1E-5  # Delta m^2_21 maximum value (3-sigma) [eV^2]

# Atmospheric mass splitting
Dm31 = 2.63 * 1E-3  # Delta m^2_31 maximum value (3-sigma) [eV^2]

# NEUTRINO MIXING ANGLES
PMNSs12 = 5.640000 * 1E-1  # PMNS: sin(theta_12) maximum value (3-sigma)
PMNSs23 = 7.580000 * 1E-1  # PMNS: sin(theta_23) maximum value (3-sigma)
PMNSs13 = 1.480000 * 1E-1  # PMNS: sin(theta_13) maximum value (3-sigma)
PMNSphase = 0  # PMNS: CP-violating phase maximum value (3-sigma)

# NEUTRINO MASSES
mv1 = 8.600000 * 1E-2  # maximum mass of the lightest neutrino [eV]

# HIGGS PHYSICS
mh = 1.2500 * 1E+2  # maximum mass of the 125-Higgs boson [GeV]

# QUARK MIXING PHYSICS
CKMs12 = 2.265000 * 1E-1
CKMs23 = 4.053000 * 1E-2
CKMs13 = 3.610000 * 1E-3
CKMphase = 1.196000


########################
# DEPENDENT PARAMETERS #
########################

def vphi(log_vDelta):
    cH = costh(log_vDelta)
    return v * cH


def vDelta(log_vDelta):
    return 10 ** log_vDelta


def sinth(log_vDelta):
    vD = vDelta(log_vDelta)
    return 2 * np.sqrt(2) * vD / v


def costh(log_vDelta):
    sH = sinth(log_vDelta)
    return cmath.sqrt(1 - sH ** 2)


def tanth(log_vDelta):
    sH = sinth(log_vDelta)
    cH = costh(log_vDelta)
    return sH / cH


def mass3(log_vDelta, M1, lam5):
    vD = vDelta(log_vDelta)
    return cmath.sqrt(M1 / (4 * vD) +  lam5 / 2) * v


def mass5(log_vDelta, M1, M2, lam3, lam5):
    vp = vphi(log_vDelta)
    vD = vDelta(log_vDelta)
    x = M1 * (vp ** 2) / (4 * vD)
    y = 12 * M2 * vD
    z = 8 * lam3 * (vD ** 2)
    w = (3 / 2) * lam5 * (vp ** 2)
    return cmath.sqrt(x + y + z + w)


def M12square(log_vDelta, M1, lam2, lam5):
    vp = vphi(log_vDelta)
    vD = vDelta(log_vDelta)
    p = (np.sqrt(3) / 2) * vp
    return p * (-M1 + 4 * vD * (2 * lam2 - lam5))


def M22square(log_vDelta, M1, M2, lam3, lam4):
    vp = vphi(log_vDelta)
    vD = vDelta(log_vDelta)
    x = M1 * (vp ** 2) / (4 * vD)
    y = -6 * M2 * vD
    z = 8 * (vD ** 2) * (lam3 + 3 * lam4)
    return x + y + z


def lambda1(log_vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(log_vDelta)
    M12sq = M12square(log_vDelta, M1, lam2, lam5)
    M22sq = M22square(log_vDelta, M1, M2, lam3, lam4)
    p = 1 / (8 * vp ** 2)
    x = mh ** 2
    y = M12sq ** 2 / (M22sq - mh ** 2)
    return p * (x + y)


def M11square(log_vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(log_vDelta)
    lam1 = lambda1(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return 8 * lam1 * (vp ** 2)


def musquarePhi(log_vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(log_vDelta)
    vD = vDelta(log_vDelta)
    lam1 = lambda1(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)
    x = -4 * lam1 * (vp ** 2)
    y = -3 * (2 * lam2 - lam5) * (vD ** 2)
    z = (3 / 2) * M1 * vD
    return x + y + z


def musquareDelta(log_vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(log_vDelta)
    vD = vDelta(log_vDelta)
    x = -(2 * lam2 - lam5) * (vp ** 2)
    y = -4 * (lam3 + 3 * lam4) * (vD ** 2)
    z = M1 * (vp ** 2) / (4 * vD)
    w = 6 * M2 * vD
    return x + y + z + w


def massH(log_vDelta, M1, M2, lam2, lam3, lam4, lam5):
    M11sq = M11square(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)
    M22sq = M22square(log_vDelta, M1, M2, lam3, lam4)
    x = M11sq + M22sq
    y = np.square(mh)
    return cmath.sqrt(x - y)


def alpha(log_vDelta, M1, M2, lam2, lam3, lam4, lam5):
    M12sq = M12square(log_vDelta, M1, lam2, lam5)
    mH = massH(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)
    numerator = 2 * M12sq
    denominator = mH ** 2 - mh ** 2
    return 0.5 * np.arcsin(numerator / denominator)


# Mass of the second neutrino
def massv2(mv1, Dm21):
    return np.sqrt(mv1 ** 2 + Dm21)


# Mass of the third neutrino
def massv3(mv1, Dm31):
    return np.sqrt(mv1 ** 2 + Dm31)


# CKM matrix
def CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase):
    CKMc12 = math.sqrt(1 - CKMs12 ** 2)
    CKMc23 = math.sqrt(1 - CKMs23 ** 2)
    CKMc13 = math.sqrt(1 - CKMs13 ** 2)
    return np.array([[CKMc12 * CKMc13, CKMs12 * CKMc13, CKMs13 * np.exp(-1j*CKMphase)],
                    [-CKMs12 * CKMc23 - CKMc12 * CKMs23 * CKMs13 * np.exp(1j*CKMphase), CKMc12 * CKMc23 - CKMs12 * CKMs23 * CKMs13 * np.exp(1j*CKMphase), CKMs23 * CKMc13],
                    [CKMs12 * CKMs23 - CKMc12 * CKMc23 * CKMs13 * np.exp(1j*CKMphase), -CKMc12 * CKMs23 - CKMs12 * CKMc23 * CKMs13 * np.exp(1j*CKMphase), CKMc23 * CKMc13]])


# PMNS matrix
def PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase):
    """We are assuming that PMNSphase = 0"""
    PMNSc12 = math.sqrt(1 - PMNSs12 ** 2)
    PMNSc23 = math.sqrt(1 - PMNSs23 ** 2)
    PMNSc13 = math.sqrt(1 - PMNSs13 ** 2)
    return np.array([[PMNSc12 * PMNSc13, PMNSs12 * PMNSc13, PMNSs13 * np.exp(-1j*PMNSphase)],
                    [-PMNSs12 * PMNSc23 - PMNSc12 * PMNSs23 * PMNSs13 * np.exp(1j*PMNSphase), PMNSc12 * PMNSc23 - PMNSs12 * PMNSs23 * PMNSs13 * np.exp(1j*PMNSphase), PMNSs23 * PMNSc13],
                    [PMNSs12 * PMNSs23 - PMNSc12 * PMNSc23 * PMNSs13 * np.exp(1j*PMNSphase), -PMNSc12 * PMNSs23 - PMNSs12 * PMNSc23 * PMNSs13 * np.exp(1j*PMNSphase), PMNSc23 * PMNSc13]])


# Diagonalized charged lepton yukawa matrix
def Yl_DIAG_MATRIX(log_vDelta):
    vp = vphi(log_vDelta)
    return np.array([[(np.sqrt(2) * me) / vp, 0, 0],
                    [0, (np.sqrt(2) * mmu) / vp, 0],
                    [0, 0, (np.sqrt(2) * mtau) / vp]])


# Diagonalized neutrino yukawa matrix
"""Note the factor 10^-9 with neutrino masses; this is because they are
expected to be given in unit of [eV], while v_Delta is given in units of [GeV]"""
def Yv_DIAG_MATRIX(log_vDelta, mv1, Dm21, Dm31):
    vD = vDelta(log_vDelta)
    mv2 = massv2(mv1, Dm21)
    mv3 = massv3(mv1, Dm31)
    return np.array([[(mv1 * 1E-9) / vD, 0, 0],
                    [0, (mv2 * 1E-9) / vD, 0],
                    [0, 0, (mv3 * 1E-9) / vD]])


# Diagonalized up-quark yukawa matrix
def Yu_DIAG_MATRIX(log_vDelta):
    vp = vphi(log_vDelta)
    return np.array([[(np.sqrt(2) * mup) / vp, 0, 0],
                    [0, (np.sqrt(2) * mcharm) / vp, 0],
                    [0, 0, (np.sqrt(2) * mtop) / vp]])


# Diagonalized down-quark yukawa matrix
def Yd_DIAG_MATRIX(log_vDelta):
    vp = vphi(log_vDelta)
    return np.array([[(np.sqrt(2) * mdown) / vp, 0, 0],
                    [0, (np.sqrt(2) * mstrange) / vp, 0],
                    [0, 0, (np.sqrt(2) * mbottom) / vp]])



# DEPENDENT PARAMETERS

lam1 = lambda1(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)  # lambda_1 parameter

tH = tanth(log_vDelta)  # tan(theta_H)
vp = vphi(log_vDelta)  # v_phi, i.e., vacuum expectation of the Higgs doublet
vD = vDelta(log_vDelta)  # v_Delta, i.e., vacuum expectation of the triplet

m3 = mass3(log_vDelta, M1, lam5)  # mass of the H3 scalars
m5 = mass5(log_vDelta, M1, M2, lam3, lam5)  # mass of the H5 scalars
mH = massH(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)  # mass of the H scalar

M12sq = M12square(log_vDelta, M1, lam2, lam5)  # M^2_12 auxiliary variable
M22sq = M22square(log_vDelta, M1, M2, lam3, lam4)  # M^2_22 auxiliary variable
M11sq = M11square(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)  # M^2_11 auxiliary variable

musqPhi = musquarePhi(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)  # mu^2_Phi parameter
musqDelta = musquareDelta(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)  # mu^2_Delta parameter

A = alpha(log_vDelta, M1, M2, lam2, lam3, lam4, lam5)  # alpha mixing angle
sA = np.sin(A)  # sin(alpha)

mv2 = massv2(mv1, Dm21)  # mass of the 2nd neutrino
mv3 = massv3(mv1, Dm31)  # mass of the 3rd neutrino

# CKM MATRIX
VCKM = PMNS_MATRIX(CKMs12, CKMs23, CKMs23, CKMphase)  # CKM matrix
VCKM_DAGGER = np.conjugate(np.transpose(VCKM))  # conjugate transpose of CKM matrix

# PMNS MATRIX
VPMNS = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs23, PMNSphase)  # PMNS matrix
VPMNS_DAGGER = np.conjugate(np.transpose(VPMNS))  # conjugate transpose of PMNS matrix

# CHARGED LEPTON YUKAWA MATRIX
YlD = Yl_DIAG_MATRIX(log_vDelta)  # diagonal charged lepton yukawa matrix

# NEUTRINO YUKAWA MATRIX
YvD = Yv_DIAG_MATRIX(log_vDelta, mv1, Dm21, Dm31)  # diagonal neutrino yukawa matrix
YvND = np.dot(np.dot(VPMNS_DAGGER, YvD), VPMNS)  # non-diagonal neutrino yukawa matrix

# CHARGED UP QUARK YUKAWA MATRIX
YuD = Yu_DIAG_MATRIX(log_vDelta)  # diagonal charged up quark yukawa matrix

# CHARGED DOWN QUARK YUKAWA MATRIX
YdD = Yd_DIAG_MATRIX(log_vDelta)  # diagonal charged down quark yukawa matrix

#####################
# PRINTING THE DATA #
#####################

print('sin(alpha): ' + str(round(np.real(sA), r1)))

print('tan(theta_H): '+ str(round(np.real(tH), r1)))

print('M1: ' + str(round(np.real(M1), r1)))

print('M2: ' + str(round(np.real(M2), r1)))

if 0 < lam1 < (1 / 3) * np.pi:
    print('lam1: ' + str(round(np.real(lam1), r1)))
else:
    print('lam1: ' + str(round(np.real(lam1), r1)) + '\033[91m [WARNING: OUT OF RANGE!]\033[00m')

if -(2 / 3) * np.pi < lam2 < (2 / 3) * np.pi:
    print('lam2: ' + str(round(np.real(lam2), r1)))
else:
    print('lam2: ' + str(round(np.real(lam2), r1)) + '\033[91m [WARNING: OUT OF RANGE!]\033[00m')

if -(1 / 2) * np.pi < lam3 < (3 / 5) * np.pi:
    print('lam3: ' + str(round(np.real(lam3), r1)))
else:
    print('lam3: ' + str(round(np.real(lam3), r1)) + '\033[91m [WARNING: OUT OF RANGE!]\033[00m')

if -(1 / 5) * np.pi < lam4 < (1 / 2) * np.pi:
    print('lam4: ' + str(round(np.real(lam4), r1)))
else:
    print('lam4: ' + str(round(np.real(lam4), r1)) + '\033[91m [WARNING: OUT OF RANGE!]\033[00m')

if -(8 / 3) * np.pi < lam5 < (8 / 3) * np.pi:
    print('lam5: ' + str(round(np.real(lam5), r1)))
else:
    print('lam5: ' + str(round(np.real(lam5), r1)) + '\033[91m [WARNING: OUT OF RANGE!]\033[00m')

print('v_phi: ' + str(round(np.real(vp), r1)))

print('v_Delta: ' + str(round(np.real(vD), r1)))

print('m_3: ' + str(round(np.real(m3), r1)))

print('m_5: ' + str(round(np.real(m5), r1)))

print('m_h: ' + str(round(np.real(mh), r1)))

print('m_H: ' + str(round(np.real(mH), r1)))

print('mu^2_Phi: ' + str(round(np.real(musqPhi), r1)))

print('mu^2_Delta: ' + str(round(np.real(musqDelta), r1)))

print('Dm21: ' + str(round(np.real(Dm21), r2)))

print('Dm31: ' + str(round(np.real(Dm31), r2)))

print('PMNSs12: ' + str(round(np.real(PMNSs12), r1)))

print('PMNSs23: ' + str(round(np.real(PMNSs23), r1)))

print('PMNSs13: ' + str(round(np.real(PMNSs13), r1)))

print('PMNSphase: ' + str(round(np.real(PMNSphase), r1)))

print('m_v1: ' + str(round(np.real(mv1), r1)))

print('m_v2: ' + str(round(np.real(mv2), r1)))

print('m_v3: ' + str(round(np.real(mv3), r1)))

print('|VCKM_11|: ' + str(round(np.abs(VCKM[0, 0]), r1)))

print('|VCKM_12|: ' + str(round(np.abs(VCKM[0, 1]), r1)))

print('|VCKM_13|: ' + str(round(np.abs(VCKM[0, 2]), r1)))

print('|VCKM_21|: ' + str(round(np.abs(VCKM[1, 0]), r1)))

print('|VCKM_22|: ' + str(round(np.abs(VCKM[1, 1]), r1)))

print('|VCKM_23|: ' + str(round(np.abs(VCKM[1, 2]), r1)))

print('|VCKM_31|: ' + str(round(np.abs(VCKM[2, 0]), r1)))

print('|VCKM_32|: ' + str(round(np.abs(VCKM[2, 1]), r1)))

print('|VCKM_33|: ' + str(round(np.abs(VCKM[2, 2]), r1)))

print('|VPMNS_11|: ' + str(round(np.abs(VPMNS[0, 0]), r1)))

print('|VPMNS_12|: ' + str(round(np.abs(VPMNS[0, 1]), r1)))

print('|VPMNS_13|: ' + str(round(np.abs(VPMNS[0, 2]), r1)))

print('|VPMNS_21|: ' + str(round(np.abs(VPMNS[1, 0]), r1)))

print('|VPMNS_22|: ' + str(round(np.abs(VPMNS[1, 1]), r1)))

print('|VPMNS_23|: ' + str(round(np.abs(VPMNS[1, 2]), r1)))

print('|VPMNS_31|: ' + str(round(np.abs(VPMNS[2, 0]), r1)))

print('|VPMNS_32|: ' + str(round(np.abs(VPMNS[2, 1]), r1)))

print('|VPMNS_33|: ' + str(round(np.abs(VPMNS[2, 2]), r1)))

print('YvND_11: ' + str(round(np.real(YvND[0, 0]), r2)))

print('YvND_12: ' + str(round(np.real(YvND[0, 1]), r2)))

print('YvND_13: ' + str(round(np.real(YvND[0, 2]), r2)))

print('YvND_21: ' + str(round(np.real(YvND[1, 0]), r2)))

print('YvND_22: ' + str(round(np.real(YvND[1, 1]), r2)))

print('YvND_23: ' + str(round(np.real(YvND[1, 2]), r2)))

print('YvND_31: ' + str(round(np.real(YvND[2, 0]), r2)))

print('YvND_32: ' + str(round(np.real(YvND[2, 1]), r2)))

print('YvND_33: ' + str(round(np.real(YvND[2, 2]), r2)))

print('y_v1: ' + str(round(np.real(YvD[0, 0]), r2)))

print('y_v2: ' + str(round(np.real(YvD[1, 1]), r2)))

print('y_v3: ' + str(round(np.real(YvD[2, 2]), r2)))

print('y_electron: ' + str(round(np.real(YlD[0, 0]), r2)))

print('y_muon: ' + str(round(np.real(YlD[1, 1]), r2)))

print('y_tau: ' + str(round(np.real(YlD[2, 2]), r2)))

print('y_up: ' + str(round(np.real(YuD[0, 0]), r2)))

print('y_charm: ' + str(round(np.real(YuD[1, 1]), r2)))

print('y_top: ' + str(round(np.real(YuD[2, 2]), r2)))

print('y_down: ' + str(round(np.real(YdD[0, 0]), r2)))

print('y_strange: ' + str(round(np.real(YdD[1, 1]), r2)))

print('y_bottom: ' + str(round(np.real(YdD[2, 2]), r2)))
