# -*- coding: utf-8 -*-
#!/usr/bin/env python3.8

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + type-II seesaw mechanism.
INFO: Intended    -> Random scan of the parameter space.
                     (!) Logarithmic scan.
INFO: Language    -> Python 3.8
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> November 23, 2022.
"""


###########
# MODULES #
###########

import logging
logging.basicConfig(
    #filename = 'file.log',
    format = '%(levelname)s:%(message)s',
    level = logging.DEBUG)

logging.info('Loading libraries')

import numpy as np
import pandas as pd
import math
import cmath
import random as rd
from time import perf_counter
from tqdm import tqdm

# timer
start = perf_counter()

#############
# CONSTANTS #
#############

logging.info('Setting parameters')

from ConstantsModule import *

# Internal parameters
N = 12  # Number of points to be found
n = 0  # Initial number of points

r1 = 10  # round to r1 decimals
r2 = 20  # round to r2 decimals (used for neutrino yukawas and other very small parameters)


################################
# INDEPENDENT PARAMETER RANGES #
################################

logging.info('Setting parameter rangers')

# Maximum ranges for vDelta and M1,M2 coefficients
vDelta_max = 80  # [GeV]
vDelta_min = 1e-12  # [GeV]

M1coeff_max = 1000  # [GeV]
M1coeff_min = 1e-09  # [GeV]
# M1coeff_lower_bound = 1e-09  # [GeV]

M2coeff_max = 1000  # [GeV]
M2coeff_min = -1000  # [GeV]
M2coeff_lower_bound = 1e-09  # [GeV]


# Maximum ranges for lambda parameters
# Ref.: https://arxiv.org/abs/1404.2640
lam1_max = (1 / 3) * np.pi
lam1_min = 0

lam2_max = (2 / 3) * np.pi
lam2_min = -(2 / 3) * np.pi
lam2_lower_bound = 1e-09  # this means that lam2 will be chosen in the interval [-2pi/3,-1e-09] U [1e-09,2pi/3]

lam3_max = (3 / 5) * np.pi
lam3_min = -(1 / 2) * np.pi
lam3_lower_bound = 1e-09

lam4_max = (1 / 2) * np.pi
lam4_min = -(1 / 5) * np.pi
lam4_lower_bound = 1e-09

lam5_max = (8 / 3) * np.pi
lam5_min = -(8 / 3) * np.pi
lam5_lower_bound = 1e-09


# NEUTRINO PHYSICS
# All values assume a normal neutrino mass ordering
# Ref.: https://link.springer.com/article/10.1007/JHEP02(2021)071

# MASS SPLITTINGS
# Solar mass splitting
Dm21_max = 8.14 * 1E-5  # Delta m^2_21 maximum value (3-sigma) [eV^2]
Dm21_min = 6.94 * 1E-5  # Delta m^2_21 minimum value (3-sigma) [eV^2]

# Atmospheric mass splitting
Dm31_max = 2.63 * 1E-3  # Delta m^2_31 maximum value (3-sigma) [eV^2]
Dm31_min = 2.47 * 1E-3  # Delta m^2_31 minimum value (3-sigma) [eV^2]

# NEUTRINO MIXING ANGLES
PMNSs12_max = 0.607454  # PMNS: sin(theta_12) maximum value (3-sigma)
PMNSs12_min = 0.520577  # PMNS: sin(theta_12) minimum value (3-sigma)

PMNSs23_max = 0.781025  # PMNS: sin(theta_23) maximum value (3-sigma)
PMNSs23_min = 0.658787  # PMNS: sin(theta_23) minimum value (3-sigma)

PMNSs13_max = 0.155081  # PMNS: sin(theta_13) maximum value (3-sigma)
PMNSs13_min = 0.141421  # PMNS: sin(theta_13) minimum value (3-sigma)

PMNSphase_max = 0  # PMNS: CP-violating phase maximum value (3-sigma)
PMNSphase_min = 0  # PMNS: CP-violating phase minimum value (3-sigma)

# NEUTRINO MASSES
massv1_max = 0.1  # maximum mass of the lightest neutrino [eV]
massv1_min = 0  # minimum value of the lightest neutrino [eV]

# QUARK MIXING PHYSICS
# Ref.: R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
CKMs12_max = 0.22698
CKMs12_min = 0.22602

CKMs23_max = 0.04136
CKMs23_min = 0.03992

CKMs13_max = 0.00372
CKMs13_min = 0.00352

CKMphase_max = 1.241
CKMphase_min = 1.153


########################
# DEPENDENT PARAMETERS #
########################

logging.info('Setting dependent parameters')

from FunctionsModule import *


#############################
# OPEN FILE AND SET COLUMNS #
#############################

logging.info('Creating the data list')

points_dict = {
        'sA': [],
        'log_sA': [],
        'sA_inv': [],
        'sA/sA_inv': [],
        'tH': [],
        'log_tH': [],
        'M1 [GeV]': [],
        'log_M1': [],
        'M2 [GeV]': [],
        'log_M2': [],
        'lam1': [],
        'log_lam1': [],
        'lam1_inv': [],
        'lam1/lam1_inv': [],
        'lam2': [],
        'log_lam2': [],
        'lam2_inv': [],
        'lam2/lam2_inv': [],
        'lam3': [],
        'log_lam3': [],
        'lam3_inv': [],
        'lam3/lam3_inv': [],
        'lam4': [],
        'log_lam4': [],
        'lam4_inv': [],
        'lam4/lam4_inv': [],
        'lam5': [],
        'log_lam5': [],
        'lam5_inv': [],
        'lam5/lam5_inv': [],
        'vphi [GeV]': [],
        'vDelta [GeV]': [],
        'log_vDelta': [],
        'm3 [GeV]': [],
        'm3_inv [GeV]': [],
        'm3/m3_inv': [],
        'm5 [GeV]': [],
        'm5_inv [GeV]': [],
        'm5/m5_inv': [],
        'mh [GeV]': [],
        'mH [GeV]': [],
        'mH_inv [GeV]': [],
        'mH/mH_inv': [],
        'm3/m5': [],
        'm3/mH': [],
        'm5/mH': [],
        'mu2Phi [GeV^2]': [],
        'mu2Delta [GeV^2]': [],
        'Dm21 [eV^2]': [],
        'Dm31 [eV^2]': [],
        'PMNSs12': [],
        'PMNSs23': [],
        'PMNSs13': [],
        'PMNSphase': [],
        'mv1 [eV]': [],
        'mv2 [eV]': [],
        'mv3 [eV]': [],
        'CKMs12': [],
        'CKMs23': [],
        'CKMs13': [],
        'CKMphase': [],
        'VCKM11': [],
        'VCKM12': [],
        'VCKM13': [],
        'VCKM21': [],
        'VCKM22': [],
        'VCKM23': [],
        'VCKM31': [],
        'VCKM32': [],
        'VCKM33': [],
        'VPMNS11': [],
        'VPMNS12': [],
        'VPMNS13': [],
        'VPMNS21': [],
        'VPMNS22': [],
        'VPMNS23': [],
        'VPMNS31': [],
        'VPMNS32': [],
        'VPMNS33': [],
        'YvND11': [],
        'log_YvND11': [],
        'YvND12': [],
        'log_YvND12': [],
        'YvND13': [],
        'log_YvND13': [],
        'YvND21': [],
        'log_YvND21': [],
        'YvND22': [],
        'log_YvND22': [],
        'YvND23': [],
        'log_YvND23': [],
        'YvND31': [],
        'log_YvND31': [],
        'YvND32': [],
        'log_YvND32': [],
        'YvND33': [],
        'log_YvND33': [],
        'yv1': [],
        'log_yv1': [],
        'yv2': [],
        'log_yv2': [],
        'yv3': [],
        'log_yv3': [],
        'yelectron': [],
        'log_yelectron': [],
        'ymuon': [],
        'log_ymuon': [],
        'ytau': [],
        'log_ytau': [],
        'yup': [],
        'log_yup': [],
        'ycharm': [],
        'log_ycharm': [],
        'ytop': [],
        'log_ytop': [],
        'ydown': [],
        'log_ydown': [],
        'ystrange': [],
        'log_ystrange': [],
        'ybottom': [],
        'log_ybottom': [],
        'kappaV': [],
        'kappaf': [],
        'RbGM': [],
        'RC(RbGM)%': [],
        'RGMDeltam': [],
        'RC(RGMDeltam)%': [],
        'RGMsmu': [],
        'RC(RGMsmu)%': []
    }


###################################
# MAIN LOOP: LOOKING VALID POINTS #
###################################

logging.info('Initiating the main loop')

pbar = tqdm(total=N, desc='COMPATIBLE POINTS FOUND', unit='points')

while n < N:

    ################################
    # INITIAL SET OF RANDOM VALUES #
    ################################

    # RANDOM INDEPENDENT PARAMETERS
    log_vDelta, vDelta = rnd_log_possitive(vDelta_min, vDelta_max)
    log_M1, M1 = rnd_log_possitive(M1coeff_min, M1coeff_max)
    log_M2, M2 = rnd_log_possitivenegative(M2coeff_min, M2coeff_max, M2coeff_lower_bound)

    log_lam2, lam2 = rnd_log_possitivenegative(lam2_min, lam2_max, lam2_lower_bound)
    log_lam3, lam3 = rnd_log_possitivenegative(lam3_min, lam3_max, lam3_lower_bound)
    log_lam4, lam4 = rnd_log_possitivenegative(lam4_min, lam4_max, lam4_lower_bound)
    log_lam5, lam5 = rnd_log_possitivenegative(lam5_min, lam5_max, lam5_lower_bound)

    Dm21 = GRN(Dm21_min, Dm21_max)
    Dm31 = GRN(Dm31_min, Dm31_max)
    PMNSs12 = GRN(PMNSs12_min, PMNSs12_max)
    PMNSs23 = GRN(PMNSs23_min, PMNSs23_max)
    PMNSs13 = GRN(PMNSs13_min, PMNSs13_max)
    PMNSphase = GRN(PMNSphase_min, PMNSphase_max)

    mv1 = GRN(massv1_min, massv1_max)

    CKMs12 = GRN(CKMs12_min, CKMs12_max)
    CKMs23 = GRN(CKMs23_min, CKMs23_max)
    CKMs13 = GRN(CKMs13_min, CKMs13_max)
    CKMphase = GRN(CKMphase_min, CKMphase_max)

    # DEPENDENT PARAMETERS
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)  # lambda_1 parameter

    tH = tanth(vDelta)  # tan(theta_H)
    vp = vphi(vDelta)  # v_phi, i.e., vacuum expectation of the Higgs doublet

    m3 = mass3(vDelta, M1, lam5)  # mass of the H3 scalars
    m5 = mass5(vDelta, M1, M2, lam3, lam5)  # mass of the H5 scalars
    mH = massH(vDelta, M1, M2, lam2, lam3, lam4, lam5)  # mass of the H scalar

    M12sq = M12square(vDelta, M1, lam2, lam5)  # M^2_12 auxiliary variable
    M22sq = M22square(vDelta, M1, M2, lam3, lam4)  # M^2_22 auxiliary variable
    M11sq = M11square(vDelta, M1, M2, lam2, lam3, lam4, lam5)  # M^2_11 auxiliary variable

    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)  # mu^2_Phi parameter
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)  # mu^2_Delta parameter

    sA = sinA(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    mv2 = massv2(mv1, Dm21)  # mass of the 2nd neutrino
    mv3 = massv3(mv1, Dm31)  # mass of the 3rd neutrino

    # RECOVERED BY INVERSION
    lam1_inv = lambda1_inv(sA, vDelta, mH)
    lam2_inv = lambda2_inv(sA, vDelta, M1, m3, mH)
    lam3_inv = lambda3_inv(vDelta, M1, M2, m3, m5)
    lam4_inv = lambda4_inv(sA, vDelta, M1, M2, m3, m5, mH)
    lam5_inv = lambda5_inv(vDelta, M1, m3)

    m3_inv = mass3(vDelta, M1, lam5_inv)  # mass of the H3 scalars from the lambda-inverted parameters
    m5_inv = mass5(vDelta, M1, M2, lam3_inv, lam5_inv)  # mass of the H5 scalars from the lambda-inverted parameters
    mH_inv = massH(vDelta, M1, M2, lam2_inv, lam3_inv, lam4_inv, lam5_inv)  # mass of the H scalar from the lambda-inverted parameters

    sA_inv = sinA(vDelta, M1, M2, lam2_inv, lam3_inv, lam4_inv, lam5_inv)

    # CKM MATRIX
    VCKM = PMNS_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)  # CKM matrix
    VCKM_DAGGER = np.conjugate(np.transpose(VCKM))  # conjugate transpose of CKM matrix

    # PMNS MATRIX
    VPMNS = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)  # PMNS matrix
    VPMNS_DAGGER = np.conjugate(np.transpose(VPMNS))  # conjugate transpose of PMNS matrix

    # CHARGED LEPTON YUKAWA MATRIX
    YlD = Yl_DIAG_MATRIX(vDelta)  # diagonal charged lepton yukawa matrix

    # NEUTRINO YUKAWA MATRIX
    YvD = Yv_DIAG_MATRIX(vDelta, mv1, Dm21, Dm31)  # diagonal neutrino yukawa matrix
    YvND = np.dot(np.dot(VPMNS_DAGGER, YvD), VPMNS)  # non-diagonal neutrino yukawa matrix

    # CHARGED UP QUARK YUKAWA MATRIX
    YuD = Yu_DIAG_MATRIX(vDelta)  # diagonal charged up quark yukawa matrix
    
    # CHARGED DOWN QUARK YUKAWA MATRIX
    YdD = Yd_DIAG_MATRIX(vDelta)  # diagonal charged down quark yukawa matrix

    
    # ALL THE MINIMA
    Vmin_real = VScalarPotential(0, vp, 0, vDelta, 0, 0, vDelta, vDelta, M1, M2, lam2, lam3, lam4, lam5)

    # IN THE CP-EVEN SUBSPACE
    Vmin_V0p1p = V0p1p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p1m = V0p1m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p2p = V0p2p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p2m = V0p2m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p3p = V0p3p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p3m = V0p3m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p4 = V0p4(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p5p = V0p5p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p5m = V0p5m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p6p = V0p6p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p6m = V0p6m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p7p = V0p7p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p7m = V0p7m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p8p = V0p8p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0p8m = V0p8m(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    # IN THE CP-ODD SUBSPACE
    Vmin_V0m1p = V0m1p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0m1m = V0m1m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0m2p = V0m2p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0m2m = V0m2m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_V0m3 = V0m3(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    # IN THE SINGLY-CHARGED SUBSPACE
    Vmin_Vpm1p = Vpm1p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm1m = Vpm1m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm2p = Vpm2p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm2m = Vpm2m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm3p = Vpm3p(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm3m = Vpm3m(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm4 = Vpm4(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm5 = Vpm5(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    Vmin_Vpm6 = Vpm6(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    # IN THE DOUBLY-CHARGED SUBSPACE
    Vmin_Vppmm = Vppmm(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    ###############
    # CONSTRAINTS #
    ###############

    # EXCLUDING NON-ZERO DIVISION
    if tH > 0.0:
        pass
    else:
        continue

    # MASS CONSTRAINTS

    # CHECKING THAT m3 > 0
    if np.isreal(m3) and np.real(m3) >= 0.0 and np.real(m3) <= 1000:
        pass
    else:
        continue

    # CHECKING THAT m5 > 0
    if np.isreal(m5) and np.real(m5) >= 0.0 and np.real(m5) <= 1000:
        pass
    else:
        continue

    # CHECKING THAT mH > 0
    if np.isreal(mH) and np.real(mH) >= Mh and np.real(mH) <= 1000:
        pass
    else:
        continue


    # UNITARITY CONSTRAINTS
    # Ref.: https://arxiv.org/abs/1404.2640

    # 1rt UNITARITY CONSTRAINT
    if np.real(cmath.sqrt(np.square(6 * lam1 - 7 * lam3 - 11 * lam4)
                          + 36 * np.square(lam2)) + abs(6 * lam1 + 7 * lam3 + 11 * lam4)) < 4 * np.pi:
        pass
    else:
        continue

    # 2nd UNITARITY CONSTRAINT
    if np.real(cmath.sqrt(np.square(2 * lam1 + lam3 - 2 * lam4)
                          + np.square(lam5)) + abs(2 * lam1 - lam3 + 2 * lam4)) < 4 * np.pi:
        pass
    else:
        continue

    # 3rd UNITARITY CONSTRAINT
    if np.real(abs(2 * lam3 + lam4)) < np.pi:
        pass
    else:
        continue

    # 4th UNITARITY CONSTRAINT
    if np.real(abs(lam2 - lam5)) < 2 * np.pi:
        pass
    else:
        continue


    # POTENTIAL BOUNDED-FROM-BELOW CONSTRAINTS
    # Ref.: https://arxiv.org/abs/1404.2640

    # 1st BOUNDED-FROM-BELOW CONSTRAIN
    if lam1 > 0.0:
        pass
    else:
        continue

    # 2nd BOUNDED-FROM-BELOW CONSTRAIN
    if lam3 >= 0:
        if lam4 > (-1/3) * lam3:
            pass
        else:
            continue
    else:
        if lam4 > -lam3:
            pass
        else:
            continue

    # 3rd BOUNDED-FROM-BELOW CONSTRAIN
    def omega_plus(zeta):
        # Inner auxiliary function B
        def B(zeta):
            return np.real(cmath.sqrt((3/2) * (zeta - 1/3)))
        # Omega_plus function
        return np.real((1/6) * (1 - B(zeta)) + (cmath.sqrt(2)/3) * cmath.sqrt((1 - B(zeta)) * (1/2 + B(zeta))))


    def omega_minus(zeta):
        # Inner auxiliary function B
        def B(zeta):
            return np.real(cmath.sqrt((3/2) * (zeta - 1/3)))
        # Omega_plus function
        return np.real((1/6) * (1 - B(zeta)) - (cmath.sqrt(2)/3) * cmath.sqrt((1 - B(zeta)) * (1/2 + B(zeta))))


    if lam5 >= 0 and lam3 >= 0:
        if lam2 > np.real((1/2) * lam5 - 2 * cmath.sqrt(lam1 * ((1/3) * lam3 + lam4))):
            pass
        else:
            continue
    elif lam5 >= 0 and lam3 < 0:
        if any(lam2 < np.real(omega_plus(zeta) * lam5 - 2 * cmath.sqrt(lam1 * (zeta * lam3 + lam4))) for zeta in np.linspace(1/3, 1, 40)):
            continue
        else:
            pass
    else:
        if any(lam2 < np.real(omega_minus(zeta) * lam5 - 2 * cmath.sqrt(lam1 * (zeta * lam3 + lam4))) for zeta in np.linspace(1/3, 1, 40)):
            continue
        else:
            pass


    # AVOIDING WRONG MINIMA
    # Ref.: https://arxiv.org/abs/2207.00142
    
    if Vmin_real < Vmin_V0p1p and Vmin_real < Vmin_V0p1m and Vmin_real < Vmin_V0p2p and Vmin_real < Vmin_V0p2m and Vmin_real < Vmin_V0p3p and Vmin_real < Vmin_V0p3m and Vmin_real < Vmin_V0p4 and Vmin_real < Vmin_V0p5p and Vmin_real < Vmin_V0p5m and Vmin_real < Vmin_V0p6p and Vmin_real < Vmin_V0p6m and Vmin_real < Vmin_V0p7p and Vmin_real < Vmin_V0p7m and Vmin_real < Vmin_V0p8p and Vmin_real < Vmin_V0p8m:
        pass
    else:
        continue


    if Vmin_real < Vmin_V0m1p and Vmin_real < Vmin_V0m1m and Vmin_real < Vmin_V0m2p and Vmin_real < Vmin_V0m2m and Vmin_real < Vmin_V0m3:
        pass
    else:
        continue


    # if Vmin_real < Vmin_Vpm1p and Vmin_real < Vmin_Vpm1m and Vmin_real < Vmin_Vpm2p and Vmin_real < Vmin_Vpm2m and Vmin_real < Vmin_Vpm3p and Vmin_real < Vmin_Vpm3m and Vmin_real < Vmin_Vpm4 and Vmin_real < Vmin_Vpm5 and Vmin_real < Vmin_Vpm6:
    #     pass
    # else:
    #     continue


    if Vmin_real < Vmin_Vppmm:
        pass
    else:
        continue
    

    # NEUTRINO CONSTRAINTS

    # EFFECTIVE ELECTRON ANTI-NEUTRINO MASS
    # Ref.: Nat. Phys. 18, 160-166 (2022)
    meff_ve = np.sqrt((mv1 ** 2) * (np.abs(VPMNS[0, 0]) ** 2) + (mv2 ** 2) * (np.abs(VPMNS[0, 1]) ** 2) + (mv3 ** 2) * (np.abs(VPMNS[0, 2]) ** 2))
    if meff_ve < 0.8:  # [eV] at 90% CL
        pass
    else:
        continue

    # EFFECTIVE MUON ANTI-NEUTRINO MASS
    # Ref.: R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
    meff_vmu = np.sqrt((mv1 ** 2) * (np.abs(VPMNS[1, 0]) ** 2) + (mv2 ** 2) * (np.abs(VPMNS[1, 1]) ** 2) + (mv3 ** 2) * (np.abs(VPMNS[1, 2]) ** 2))
    if meff_vmu < 190 * 1E3:  # [eV] at 90% CL
        pass
    else:
        continue

    # EFFECTIVE TAU NEUTRINO MASS
    # Ref.: R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
    meff_vtau = np.sqrt((mv1 ** 2) * (np.abs(VPMNS[2, 0]) ** 2) + (mv2 ** 2) * (np.abs(VPMNS[2, 1]) ** 2) + (mv3 ** 2) * (np.abs(VPMNS[2, 2]) ** 2))
    if meff_vtau < 18.2 * 1E6:  # [eV] at 95% CL
        pass
    else:
        continue

    # EFFECTIVE MAJORANA MASS OF v_e
    m_ee = np.abs(mv1 * (np.abs(VPMNS[0, 0]) ** 2) + mv2 * (np.abs(VPMNS[0, 1]) ** 2) + mv3 * (np.abs(VPMNS[0, 2]) ** 2))
    if m_ee < 0.165:
        pass
    else:
        continue

    # COSMOLOGICAL OBSERVATIONS
    # Ref.1: Phys. Rev. D 103 no. 8, (2021) 083533
    # Ref.2: Astron. Astrophys. 641 (2020) A6
    Sigma = mv1 + mv2 + mv3
    if Sigma < 0.12:  # at 95% CL
        pass
    else:
        continue

    # ACCESIBLE RANGES
    # Ref.: https://arxiv.org/abs/1404.2640
    if lam1_min < lam1 < lam1_max:
        pass
    else:
        continue

    # COLLIDER CONSTRAINTS

    # DECAYS OF THE Z
    if np.real(m5) >= MZ/2:
        pass
    else:
        continue

    # ELECTROWEAK PARAMETER TESTS
    # Ref.: https://arxiv.org/abs/1410.5538
    
    # if 0.05 - 0.11 < S_parameter_function(vDelta, M1, M2, lam2, lam3, lam4, lam5) < 0.05 + 0.11:
    # print(chi_square_S_parameter(vDelta, M1, M2, lam2, lam3, lam4, lam5))
    # if chi_square_S_parameter(vDelta, M1, M2, lam2, lam3, lam4, lam5) < 4:
    #     pass
    # else:
    #     continue

    # Z-POLE OBSERVABLE Rb
    # Ref.: https://arxiv.org/abs/1410.5538
    if RbGM(vDelta, M1, lam5) > 0.21443:
        pass
    else:
        continue

    # B0s-barB0s MIXING
    # Ref.: https://arxiv.org/abs/1410.5538
    if RGMDeltam(vDelta, M1, lam5) < 1.46:
        pass
    else:
        continue

    
    # BR(B0s to mu+ mu-)
    # Ref.: https://arxiv.org/abs/1410.5538
    if RGMsmu(vDelta, M1, lam5) < 1.43:
        pass
    else:
        continue


    # new point confirmed!

    n += 1

    pbar.update(1)

    ################
    # SAVING POINT #
    ################

    points_dict['sA'].append(round(np.real(sA), r2))
    points_dict['log_sA'].append(round(np.log10(np.abs(np.real(sA))), r2))
    points_dict['sA_inv'].append(round(np.real(sA_inv), r2))
    points_dict['sA/sA_inv'].append(round(np.real(sA) / np.real(sA_inv), r2))
    points_dict['tH'].append(round(np.real(tH), r2))
    points_dict['log_tH'].append(round(np.log10(np.abs(np.real(tH))), r2))
    points_dict['M1 [GeV]'].append(round(np.real(M1), r2))
    points_dict['log_M1'].append(round(np.real(log_M1), r2))
    points_dict['M2 [GeV]'].append(round(np.real(M2), r2))
    points_dict['log_M2'].append(round(np.real(log_M2), r2))
    points_dict['lam1'].append(round(np.real(lam1), r2))
    points_dict['log_lam1'].append(round(np.log10(np.abs(np.real(lam1))), r2))
    points_dict['lam1_inv'].append(round(np.real(lam1_inv), r2))
    points_dict['lam1/lam1_inv'].append(round(np.real(lam1) / np.real(lam1_inv), r2))
    points_dict['lam2'].append(round(np.real(lam2), r2))
    points_dict['log_lam2'].append(round(np.real(log_lam2), r2))
    points_dict['lam2_inv'].append(round(np.real(lam2_inv), r2))
    points_dict['lam2/lam2_inv'].append(round(np.real(lam2) / np.real(lam2_inv), r2))
    points_dict['lam3'].append(round(np.real(lam3), r2))
    points_dict['log_lam3'].append(round(np.real(log_lam3), r2))
    points_dict['lam3_inv'].append(round(np.real(lam3_inv), r2))
    points_dict['lam3/lam3_inv'].append(round(np.real(lam3) / np.real(lam3_inv), r2))
    points_dict['lam4'].append(round(np.real(lam4), r2))
    points_dict['log_lam4'].append(round(np.real(log_lam4), r2))
    points_dict['lam4_inv'].append(round(np.real(lam4_inv), r2))
    points_dict['lam4/lam4_inv'].append(round(np.real(lam4) / np.real(lam4_inv), r2))
    points_dict['lam5'].append(round(np.real(lam5), r2))
    points_dict['log_lam5'].append(round(np.real(log_lam5), r2))
    points_dict['lam5_inv'].append(round(np.real(lam5_inv), r2))
    points_dict['lam5/lam5_inv'].append(round(np.real(lam5) / np.real(lam5_inv), r2))
    points_dict['vphi [GeV]'].append(round(np.real(vp), r2))
    points_dict['vDelta [GeV]'].append(round(np.real(vDelta), r2))
    points_dict['log_vDelta'].append(round(np.real(log_vDelta), r2))
    points_dict['m3 [GeV]'].append(round(np.real(m3), r1))
    points_dict['m3_inv [GeV]'].append(round(np.real(m3_inv), r1))
    points_dict['m3/m3_inv'].append(round(np.real(m3) / np.real(m3_inv), r1))
    points_dict['m5 [GeV]'].append(round(np.real(m5), r1))
    points_dict['m5_inv [GeV]'].append(round(np.real(m5_inv), r1))
    points_dict['m5/m5_inv'].append(round(np.real(m5) / np.real(m5_inv), r1))
    points_dict['mh [GeV]'].append(round(np.real(Mh), r1))
    points_dict['mH [GeV]'].append(round(np.real(mH), r1))
    points_dict['mH_inv [GeV]'].append(round(np.real(mH_inv), r1))
    points_dict['mH/mH_inv'].append(round(np.real(mH) / np.real(mH_inv), r1))
    points_dict['m3/m5'].append(round(np.real(m3)/np.real(m5), r1))
    points_dict['m3/mH'].append(round(np.real(m3)/np.real(mH), r1))
    points_dict['m5/mH'].append(round(np.real(m5)/np.real(mH), r1))
    points_dict['mu2Phi [GeV^2]'].append(round(np.real(mu2Phi), r1))
    points_dict['mu2Delta [GeV^2]'].append(round(np.real(mu2Delta), r1))
    points_dict['Dm21 [eV^2]'].append(round(np.real(Dm21), r2))
    points_dict['Dm31 [eV^2]'].append(round(np.real(Dm31), r2))
    points_dict['PMNSs12'].append(round(np.real(PMNSs12), r1))
    points_dict['PMNSs23'].append(round(np.real(PMNSs23), r1))
    points_dict['PMNSs13'].append(round(np.real(PMNSs13), r1))
    points_dict['PMNSphase'].append(round(np.real(PMNSphase), r1))
    points_dict['mv1 [eV]'].append(round(np.real(mv1), r1))
    points_dict['mv2 [eV]'].append(round(np.real(mv2), r1))
    points_dict['mv3 [eV]'].append(round(np.real(mv3), r1))
    points_dict['CKMs12'].append(round(np.real(CKMs12), r1))
    points_dict['CKMs23'].append(round(np.real(CKMs23), r1))
    points_dict['CKMs13'].append(round(np.real(CKMs13), r1))
    points_dict['CKMphase'].append(round(np.real(CKMphase), r1))
    points_dict['VCKM11'].append(round(np.abs(VCKM[0, 0]), r1))
    points_dict['VCKM12'].append(round(np.abs(VCKM[0, 1]), r1))
    points_dict['VCKM13'].append(round(np.abs(VCKM[0, 2]), r1))
    points_dict['VCKM21'].append(round(np.abs(VCKM[1, 0]), r1))
    points_dict['VCKM22'].append(round(np.abs(VCKM[1, 1]), r1))
    points_dict['VCKM23'].append(round(np.abs(VCKM[1, 2]), r1))
    points_dict['VCKM31'].append(round(np.abs(VCKM[2, 0]), r1))
    points_dict['VCKM32'].append(round(np.abs(VCKM[2, 1]), r1))
    points_dict['VCKM33'].append(round(np.abs(VCKM[2, 2]), r1))
    points_dict['VPMNS11'].append(round(np.abs(VPMNS[0, 0]), r1))
    points_dict['VPMNS12'].append(round(np.abs(VPMNS[0, 1]), r1))
    points_dict['VPMNS13'].append(round(np.abs(VPMNS[0, 2]), r1))
    points_dict['VPMNS21'].append(round(np.abs(VPMNS[1, 0]), r1))
    points_dict['VPMNS22'].append(round(np.abs(VPMNS[1, 1]), r1))
    points_dict['VPMNS23'].append(round(np.abs(VPMNS[1, 2]), r1))
    points_dict['VPMNS31'].append(round(np.abs(VPMNS[2, 0]), r1))
    points_dict['VPMNS32'].append(round(np.abs(VPMNS[2, 1]), r1))
    points_dict['VPMNS33'].append(round(np.abs(VPMNS[2, 2]), r1))
    points_dict['YvND11'].append(round(np.real(YvND[0, 0]), r2))
    points_dict['log_YvND11'].append(round(np.log10(np.abs(np.real(YvND[0, 0]))), r2))
    points_dict['YvND12'].append(round(np.real(YvND[0, 1]), r2))
    points_dict['log_YvND12'].append(round(np.log10(np.abs(np.real(YvND[0, 1]))), r2))
    points_dict['YvND13'].append(round(np.real(YvND[0, 2]), r2))
    points_dict['log_YvND13'].append(round(np.log10(np.abs(np.real(YvND[0, 2]))), r2))
    points_dict['YvND21'].append(round(np.real(YvND[1, 0]), r2))
    points_dict['log_YvND21'].append(round(np.log10(np.abs(np.real(YvND[1, 0]))), r2))
    points_dict['YvND22'].append(round(np.real(YvND[1, 1]), r2))
    points_dict['log_YvND22'].append(round(np.log10(np.abs(np.real(YvND[1, 1]))), r2))
    points_dict['YvND23'].append(round(np.real(YvND[1, 2]), r2))
    points_dict['log_YvND23'].append(round(np.log10(np.abs(np.real(YvND[1, 2]))), r2))
    points_dict['YvND31'].append(round(np.real(YvND[2, 0]), r2))
    points_dict['log_YvND31'].append(round(np.log10(np.abs(np.real(YvND[2, 0]))), r2))
    points_dict['YvND32'].append(round(np.real(YvND[2, 1]), r2))
    points_dict['log_YvND32'].append(round(np.log10(np.abs(np.real(YvND[2, 1]))), r2))
    points_dict['YvND33'].append(round(np.real(YvND[2, 2]), r2))
    points_dict['log_YvND33'].append(round(np.log10(np.abs(np.real(YvND[2, 2]))), r2))
    points_dict['yv1'].append(round(np.real(YvD[0, 0]), r2))
    points_dict['log_yv1'].append(round(np.log10(np.abs(np.real(YvD[0, 0]))), r2))
    points_dict['yv2'].append(round(np.real(YvD[1, 1]), r2))
    points_dict['log_yv2'].append(round(np.log10(np.abs(np.real(YvD[1, 1]))), r2))
    points_dict['yv3'].append(round(np.real(YvD[2, 2]), r2))
    points_dict['log_yv3'].append(round(np.log10(np.abs(np.real(YvD[2, 2]))), r2))
    points_dict['yelectron'].append(round(np.real(YlD[0, 0]), r2))
    points_dict['log_yelectron'].append(round(np.log10(np.abs(np.real(YlD[0, 0]))), r2))
    points_dict['ymuon'].append(round(np.real(YlD[1, 1]), r2))
    points_dict['log_ymuon'].append(round(np.log10(np.abs(np.real(YlD[1, 1]))), r2))
    points_dict['ytau'].append(round(np.real(YlD[2, 2]), r2))
    points_dict['log_ytau'].append(round(np.log10(np.abs(np.real(YlD[2, 2]))), r2))
    points_dict['yup'].append(round(np.real(YuD[0, 0]), r2))
    points_dict['log_yup'].append(round(np.log10(np.abs(np.real(YuD[0, 0]))), r2))
    points_dict['ycharm'].append(round(np.real(YuD[1, 1]), r2))
    points_dict['log_ycharm'].append(round(np.log10(np.abs(np.real(YuD[1, 1]))), r2))
    points_dict['ytop'].append(round(np.real(YuD[2, 2]), r2))
    points_dict['log_ytop'].append(round(np.log10(np.abs(np.real(YuD[2, 2]))), r2))
    points_dict['ydown'].append(round(np.real(YdD[0, 0]), r2))
    points_dict['log_ydown'].append(round(np.log10(np.abs(np.real(YdD[0, 0]))), r2))
    points_dict['ystrange'].append(round(np.real(YdD[1, 1]), r2))
    points_dict['log_ystrange'].append(round(np.log10(np.abs(np.real(YdD[1, 1]))), r2))
    points_dict['ybottom'].append(round(np.real(YdD[2, 2]), r2))
    points_dict['log_ybottom'].append(round(np.log10(np.abs(np.real(YdD[2, 2]))), r2))
    points_dict['kappaV'].append(round(np.real(kappaV(vDelta, M1, M2, lam2, lam3, lam4, lam5)), r2))
    points_dict['kappaf'].append(round(np.real(kappaf(vDelta, M1, M2, lam2, lam3, lam4, lam5)), r2))
    points_dict['RbGM'].append(round(np.real(RbGM(vDelta, M1, lam5)), r2))
    points_dict['RC(RbGM)%'].append(round((np.real(RbGM(vDelta, M1, lam5)) - RbExp) * 100 / RbExp, r2))
    points_dict['RGMDeltam'].append(round(np.real(RGMDeltam(vDelta, M1, lam5)), r2))
    points_dict['RC(RGMDeltam)%'].append(round((np.real(RGMDeltam(vDelta, M1, lam5)) - RExpDeltam) * 100 / RExpDeltam, r2))
    points_dict['RGMsmu'].append(round(np.real(RGMsmu(vDelta, M1, lam5)), r2))
    points_dict['RC(RGMsmu)%'].append(round((np.real(RGMsmu(vDelta, M1, lam5)) - RExpsmu) * 100 / RExpsmu, r2))


pbar.close()

####################
# WRITING THE DATA #
####################

logging.info('Creating DataFrame with all data')
df = pd.DataFrame(data = points_dict)

logging.info('Saving information in .csv file')
df.to_csv(f'MasterDataPython{N}_GMSEESAW.csv', index=False)

logging.info('DONE!')

##########################
# INFORMATION ABOUT SCAN #
##########################

end = perf_counter()

runs_in_the_script = f"You will find {N} runs in the script."
average_time_per_point = "Average time to find a compatible point: " + '\n' + str(round((end - start) / N, 4)) + " seconds or " + '\n' + str(round(((end - start) / 60) / N, 4)) + " minutes."
code_execution_time = "Code execution time: " + '\n' + str(round(end - start, 4)) + " seconds or " + '\n' + str(round((end - start) / 60, 4)) + " minutes or" + '\n' + str(round((end - start) / 3600, 4)) + " hours."

logging.info(code_execution_time)
logging.info(average_time_per_point)
logging.info(runs_in_the_script)

###########################
# SEND EMAIL NOTIFICATION #
###########################

# import smtplib
# from email.mime.text import MIMEText

# sender = 'sebastian.norero.cardenas@gmail.com'
# receivers = ['sebastian.norero.cardenas@gmail.com']
# body_of_email = 'This is an automatic notification: the results of the random scan are ready!' + '\n' + runs_in_the_script + '\n' + average_time_per_point + '\n' + code_execution_time

# msg = MIMEText(body_of_email, 'plain')
# # html_body_of_email = '<h1>The html body of the email</h1>'
# # msg = MIMEText(html_body_of_email, 'html')
# msg['Subject'] = 'AUTOMATIC NOTIFICATION: Results from random scan'
# msg['From'] = sender
# msg['To'] = ','.join(receivers)

# s = smtplib.SMTP_SSL(host='smtp.gmail.com', port=465)
# s.login(user='sebastian.norero.cardenas', password='qsMc08735wkL232vUbdnq')
# s.sendmail(sender, receivers, msg.as_string())
# s.quit()

# print('INFO: Email notification sent.')
# print('INFO: DONE!')
