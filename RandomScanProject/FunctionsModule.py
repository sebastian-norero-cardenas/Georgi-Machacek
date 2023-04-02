# -*- coding: utf-8 -*-
#!/usr/bin/env python3.8

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + type-II seesaw mechanism.
INFO: Intended    -> Module storing functions used in the random scan.
INFO: Language    -> Python 3.8
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> January 07, 2023.
"""


###########
# MODULES #
###########

import numpy as np
import math
import cmath
import random as rd

from ConstantsModule import *


########################
# DEPENDENT PARAMETERS #
########################

def sinth(vDelta):
    vD = vDelta
    return 2 * np.sqrt(2) * vD / v


def costh(vDelta):
    sH = sinth(vDelta)
    return cmath.sqrt(1 - sH ** 2)


def tanth(vDelta):
    sH = sinth(vDelta)
    cH = costh(vDelta)
    return sH / cH


def vphi(vDelta):
    vD = vDelta
    return np.sqrt(v ** 2 - 8 * vD ** 2)


def mass3(vDelta, M1, lam5):
    vD = vDelta
    return cmath.sqrt(M1 / (4 * vD) +  lam5 / 2) * v


def mass5(vDelta, M1, M2, lam3, lam5):
    vD = vDelta
    vp = vphi(vDelta)
    return cmath.sqrt((M1 * vp ** 2) / (4 * vD) + 12 * M2 * vD + 8 * lam3 * vD ** 2 + (3 / 2) * (vp ** 2) * lam5)


def M12square(vDelta, M1, lam2, lam5):
    vD = vDelta
    vp = vphi(vDelta)
    p = (np.sqrt(3) / 2) * vp
    return p * (-M1 + 4 * (2 * lam2 - lam5) * vD)


def M22square(vDelta, M1, M2, lam3, lam4):
    vD = vDelta
    vp = vphi(vDelta)
    return (M1 * vp ** 2) / (4 * vD) - 6 * M2 * vD + 8 * (lam3 + 3 * lam4) * vD ** 2


def lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(vDelta)
    M12sq = M12square(vDelta, M1, lam2, lam5)
    M22sq = M22square(vDelta, M1, M2, lam3, lam4)
    p = 1 / (8 * vp ** 2)
    return p * (Mh ** 2 + M12sq ** 2 / (M22sq - Mh ** 2))


def M11square(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(vDelta)
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return 8 * lam1 * vp ** 2


def musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vD = vDelta
    vp = vphi(vDelta)
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return - 4 * lam1 * vp ** 2 - 3 * (2 * lam2 - lam5) * vD ** 2 + (3 / 2) * M1 * vD


def musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vD = vDelta
    vp = vphi(vDelta)
    return - (2 * lam2 - lam5) * vp ** 2 - 4 * (lam3 + 3 * lam4) * vD ** 2 + (M1 * vp ** 2) / (4 * vD) + 6 * M2 * vD


def massH(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    M11sq = M11square(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    M22sq = M22square(vDelta, M1, M2, lam3, lam4)
    return cmath.sqrt(M11sq + M22sq - Mh ** 2)


def sin2A(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mH = massH(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    M12sq = M12square(vDelta, M1, lam2, lam5)
    return (2 * M12sq) / (mH ** 2 - Mh ** 2)


def sinA(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    s2A = sin2A(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return np.sin(np.arcsin(s2A) / 2)


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
def Yl_DIAG_MATRIX(vDelta):
    vp = vphi(vDelta)
    return np.array([[(np.sqrt(2) * Me) / vp, 0, 0],
                    [0, (np.sqrt(2) * Mmu) / vp, 0],
                    [0, 0, (np.sqrt(2) * Mtau) / vp]])


# Diagonalized neutrino yukawa matrix
"""Note the factor 10^-9 with neutrino masses; this is because they are
expected to be given in unit of [eV], while v_Delta is given in units of [GeV]"""
def Yv_DIAG_MATRIX(vDelta, mv1, Dm21, Dm31):
    vD = vDelta
    mv2 = massv2(mv1, Dm21)
    mv3 = massv3(mv1, Dm31)
    return np.array([[(mv1 * 1E-9) / vD, 0, 0],
                    [0, (mv2 * 1E-9) / vD, 0],
                    [0, 0, (mv3 * 1E-9) / vD]])


# Diagonalized up-quark yukawa matrix
def Yu_DIAG_MATRIX(vDelta):
    vp = vphi(vDelta)
    return np.array([[(np.sqrt(2) * Mup) / vp, 0, 0],
                    [0, (np.sqrt(2) * Mcharm) / vp, 0],
                    [0, 0, (np.sqrt(2) * Mtop) / vp]])


# Diagonalized down-quark yukawa matrix
def Yd_DIAG_MATRIX(vDelta):
    vp = vphi(vDelta)
    return np.array([[(np.sqrt(2) * Mdown) / vp, 0, 0],
                    [0, (np.sqrt(2) * Mstrange) / vp, 0],
                    [0, 0, (np.sqrt(2) * Mbottom) / vp]])


# SCALING FACTOR kappa_V
def kappaV(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vD = vDelta
    vp = vphi(vDelta)
    sA = sinA(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    cA = np.sqrt(1 - sA ** 2)
    return cA * vp / v + 8 * sA * vD / (np.sqrt(3) * v)


# SCALING FACTOR kappa_f
def kappaf(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vp = vphi(vDelta)
    sA = sinA(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    cA = np.sqrt(1 - sA ** 2)
    return cA * v / vp


# INVERSE RELATIONS
def lambda1_inv(sA, vDelta, mH):
    vp = vphi(vDelta)
    cA = np.sqrt(1 - sA ** 2)
    numerator = Mh ** 2 * cA ** 2 + mH ** 2 * sA ** 2
    denominator = 8 * vp ** 2
    return numerator / denominator


def lambda2_inv(sA, vDelta, M1, m3, mH):
    vD = vDelta
    vp = vphi(vDelta)
    s2A = np.sin(2 * np.arcsin(sA))
    numerator = np.sqrt(3) * v ** 2 * s2A * (mH ** 2 - Mh ** 2) - 3 * M1 * vp * v ** 2 + 24 * m3 ** 2 * vp * vD
    denominator = 24 * vp * vD * v ** 2
    return numerator / denominator


def lambda3_inv(vDelta, M1, M2, m3, m5):
    vp = vphi(vDelta)
    vD = vDelta
    numerator = m5 ** 2 + (M1 * vp ** 2) / (2 * vD) - 12 * M2 * vD - 3 * m3 ** 2 * vp ** 2 / v ** 2
    denominator = 8 * vD ** 2
    return numerator / denominator


def lambda4_inv(sA, vDelta, M1, M2, m3, m5, mH):
    vp = vphi(vDelta)
    vD = vDelta
    cA = np.sqrt(1 - sA ** 2)
    numerator = Mh ** 2 * sA ** 2 + mH ** 2 * cA ** 2 - m5 ** 2 + 3 * m3 ** 2 * vp ** 2 / v ** 2 - 3 * M1 * vp ** 2 / (4 * vD) + 18 * M2 * vD
    denominator = 24 * vD ** 2
    return numerator / denominator


def lambda5_inv(vDelta, M1, m3):
    vD = vDelta
    numerator = 8 * vD * m3 ** 2 - 2 * M1 * v ** 2
    denominator = 4 * vD * v ** 2
    return numerator / denominator

# ELECTROWEAK PARAMETER TESTS
# Ref.: https://arxiv.org/abs/1410.5538
def f1(x, y):
    numerator = 5 * (y ** 6 - x ** 6) + 27 * (x ** 4 * y ** 2 - x ** 2 * y ** 4) + 12 * (x ** 6 - 3 * x ** 4 * y ** 2) * np.log(x) + 12 * (3 * x ** 2 * y ** 4 - y ** 6) * np.log(y)
    denominator = 36 * (y ** 2 - x ** 2) ** 3
    return numerator / denominator


def f3(x, y):
    numerator = x ** 4 - y ** 4 + 2 * x ** 2 * y ** 2 * (np.log(y ** 2) - np.log(x ** 2))
    denominator = 2 * (x ** 2 - y ** 2) ** 3
    return numerator / denominator


def S_parameter_function(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    vD = vDelta
    vp = vphi(vDelta)
    sA = sinA(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    cA = np.sqrt(1 - sA ** 2)

    m3 = mass3(vDelta, M1, lam5)
    m5 = mass5(vDelta, M1, M2, lam3, lam5)
    mH = massH(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    g_ZhH30 = - 1j * np.sqrt(2 / 3) * (ee / (sW * cW)) * (sA * vp / v + np.sqrt(3) * cA * vD / v)
    g_ZHH30 = 1j * np.sqrt(2 / 3) * (ee / (sW * cW)) * (cA * vp / v - np.sqrt(3) * sA * vD / v)
    g_ZH50H30 = - 1j * np.sqrt(1 / 3) * ((ee * vp) / (sW * cW * v))
    g_ZH5pH3m = (ee * vp) / (2 * sW * cW * v)
    g_ZZh = ((ee ** 2 * v) / (2 * sW ** 2 * cW ** 2)) * (cA * vp / v - (8 * sA * vD) / (np.sqrt(3) * v))
    g_ZZH = ((ee ** 2 * v) / (2 * sW ** 2 * cW ** 2)) * (sA * vp / v + (8 * cA * vD) / (np.sqrt(3) * v))
    g_ZZH50 = - np.sqrt(8 / 3) * (ee ** 2 / (sW ** 2 * cW ** 2)) * vD
    g_ZWpH5m = - (np.sqrt(2) * ee ** 2 * vD) / (cW * sW ** 2)
    gSM_ZZh = (ee ** 2 * v) / (2 * sW ** 2 * cW ** 2)
    mhSM = 125.18  # [GeV]

    p = (sW ** 2 * cW ** 2) / (ee ** 2 * np.pi)
    x1 = - ((ee ** 2) / (12 * sW ** 2 * cW ** 2)) * (np.log(m3 ** 2) + 5 * np.log(m5 ** 2))
    x2 = 2 * np.abs(g_ZhH30) ** 2 * f1(Mh, m3)
    x3 = 2 * np.abs(g_ZHH30) ** 2 * f1(mH, m3)
    x4 = 2 * (np.abs(g_ZH50H30) ** 2 + 2 * np.abs(g_ZH5pH3m) ** 2) * f1(m5, m3)
    x5 = np.abs(g_ZZh) ** 2 * (f1(MZ, Mh) / (2 * MZ ** 2) - f3(MZ, Mh))
    x6 = -np.abs(gSM_ZZh) ** 2 * (f1(MZ, mhSM) / (2 * MZ ** 2) - f3(MZ, mhSM))
    x7 = np.abs(g_ZZH) ** 2 * (f1(MZ, mH) / (2 * MZ ** 2) - f3(MZ, mH))
    x8 = np.abs(g_ZZH50) ** 2 * (f1(MZ, m5) / (2 * MZ ** 2) - f3(MZ, m5))
    x9 = 2 * np.abs(g_ZWpH5m) ** 2 * (f1(mW, m5) / (2 * mW ** 2) - f3(mW, m5))

    return np.real(p * (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9))


def chi_square_S_parameter(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    S = S_parameter_function(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    T = T_exp + rho_ST * (S - S_exp) * (DeltaT_exp / DeltaS_exp)

    p = 1 / (1 - rho_ST ** 2)
    x1 = (S - S_exp) ** 2 / (DeltaS_exp ** 2)
    x2 = (T - T_exp) ** 2 / (DeltaT_exp ** 2)
    x3 = - (2 * rho_ST * (S - S_exp) * (T - T_exp)) / (DeltaS_exp * DeltaT_exp)

    return p * (x1 + x2 + x3)

# Z-POLE OBSERVABLE Rb
# Ref.: https://arxiv.org/abs/1410.5538
def RbGM(vDelta, M1, lam5):
    tH = tanth(vDelta)
    m3 = mass3(vDelta, M1, lam5)
    xtW = mtop ** 2 / mW ** 2
    xt3 = mtop ** 2 / m3 ** 2
    deltaRbGM = (0.7785 / (64 * np.pi ** 2)) * (ee ** 3 / (sW ** 3 * cW)) * tH ** 2 * xtW * (xt3 / (1 - xt3) + (xt3 * np.log(xt3)) / (1 - xt3) ** 2)
    return RbSM + deltaRbGM


# B0s-barB0s MIXING
# Ref.: https://arxiv.org/abs/1410.5538 (See Eq.(35) from Ref.)
def IWW(x):
    return 1 + 9 / (1 - x) - 6 / (1 - x) ** 2 - (6 * x ** 2 * np.log(x)) / (1 - x) ** 3


def IHH(x):
    return x * ((1 + x) / (1 - x) ** 2 + (2 * x * np.log(x)) / (1 - x) ** 3)


def IWH(x, y, z):
    return y * (((2 * z - 8) * np.log(y)) / ((1 - y) ** 2 * (1 - z)) + (6 * z * np.log(x)) / ((1 - x) ** 2 * (1 - z)) - (8 - 2 * x) / ((1 - y) * (1 - x)))


def RGMDeltam(vDelta, M1, lam5):
    tH = tanth(vDelta)
    m3 = mass3(vDelta, M1, lam5)
    xtW = mtop ** 2 / mW ** 2
    xt3 = mtop ** 2 / m3 ** 2
    x3W = m3 ** 2 / mW ** 2
    return 1 + (tH ** 2 * IWH(xtW, xt3, x3W) + tH ** 4 * IHH(xt3)) / IWW(xtW)


# BR(B0s to mu+ mu-)
# Ref.: https://arxiv.org/abs/1410.5538 (See Eq.(35) from Ref.)
def RGMsmu(vDelta, M1, lam5):
    tH = tanth(vDelta)
    m3 = mass3(vDelta, M1, lam5)
    xtW = mtop ** 2 / mW ** 2
    xt3 = mtop ** 2 / m3 ** 2
    CSM10 = -0.9389 * (Mtop / 173.1) ** 1.53 * (alpha_s / 0.1184) ** (-0.09)
    CGM10 = CSM10 + tH ** 2 * (xtW / 8) * (xt3 / (1 - xt3) + (xt3 * np.log(xt3)) / (1 - xt3) ** 2)
    return np.abs(CGM10 / CSM10) ** 2


# GENERATION OF RANDOM NUMBERS
# Generate randon number in a given interval
def GRN(a, b):
    return rd.uniform(a, b)


# Generate a number between "min < number < max" whose
# logarithm was randomly generated. Here, "min, max > 0".
def rnd_log_possitive(min, max):
    log_min = np.log10(min)
    log_max = np.log10(max)
    log = rd.uniform(log_min, log_max)
    return log, 10 ** log


# Generate a number between "min < number < max" whose
# logarithm was randomly generated. Here, "min < 0" and "max > 0".
def rnd_log_possitivenegative(min, max, lower_bound):
    abs_min = np.abs(min)
    abs_max = np.abs(max)
    log_min = np.log10(abs_min)
    log_max = np.log10(abs_max)
    log_lower_bound = np.log10(lower_bound)
    r = rd.random()
    if r < abs_max / (abs_min + abs_max):
        log = rd.uniform(log_lower_bound, log_max)
        return log, 10 ** log
    else:
        log = rd.uniform(log_lower_bound, log_min)
        return log, - 10 ** log


# THE SCALAR POTENTIAL

# 2x2 SU(2) GENERATORS
tau1 = (1 / 2) * np.array([[0, 1], [1, 0]])
tau2 = (1 / 2) * np.array([[0, -1j], [1j, 0]])
tau3 = (1 / 2) * np.array([[1, 0], [0, -1]])

def SU2Generators2x2(n):
    return (tau1, tau2, tau3)[n]


# 3x3 SU(2) GENERATORS
T1 = (1 / np.sqrt(2)) * np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
T2 = (1 / np.sqrt(2)) * np.array([[0, -1j, 0], [1j, 0, -1j], [0, 1j, 0]])
T3 = (1 / np.sqrt(2)) * np.array([[1, 0, 0], [0, 0, 0], [0, 0, -1]])

def SU2Generators3x3(n):
    return (T1, T2, T3)[n]


# SIMILARITY TRANSFORMATION
U = (1 / np.sqrt(2)) * np.array([[-1, 0, 1], [-1j, 0, -1j], [0, np.sqrt(2), 0]])
UDagger = np.conjugate(np.transpose(U))

def PhiBiDoublet(phip, phi0):
    return np.array([[np.conjugate(phi0), phip],
                    [-np.conjugate(phip), phi0]])


def DeltaBiTriplet(xip, xi0, chipp, chip, chi0):
    return np.array([[np.conjugate(chi0), xip, chipp],
                    [-np.conjugate(chip), xi0, chip],
                    [np.conjugate(chipp), -np.conjugate(xip), chi0]])


def VScalarPotential(phip, phi0, xip, xi0, chipp, chip, chi0, vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)

    Phi = PhiBiDoublet(phip, phi0)
    PhiDagger = np.conjugate(np.transpose(Phi))
    Delta = DeltaBiTriplet(xip, xi0, chipp, chip, chi0)
    DeltaDagger = np.conjugate(np.transpose(Delta))

    V1 = (mu2Phi / 2) * np.trace(np.dot(PhiDagger, Phi))
    V2 = lam1 * (np.trace(np.dot(PhiDagger, Phi))) ** 2
    V3 = (mu2Delta / 2) * np.trace(np.dot(DeltaDagger, Delta))
    V4 = lam4 * (np.trace(np.dot(DeltaDagger, Delta))) ** 2 
    V5 = lam3 * np.trace(np.dot(np.dot(DeltaDagger, Delta), np.dot(DeltaDagger, Delta)))
    V6 = lam2 * np.trace(np.dot(PhiDagger, Phi)) * np.trace(np.dot(DeltaDagger, Delta))
    V7 = -lam5 * sum([np.trace(np.dot(PhiDagger, np.dot(SU2Generators2x2(a), np.dot(Phi, SU2Generators2x2(b))))) * np.trace(np.dot(DeltaDagger, np.dot(SU2Generators3x3(a), np.dot(Delta, SU2Generators3x3(b))))) for a in range(3) for b in range(3)])
    V8 = -M1 * sum([np.trace(np.dot(PhiDagger, np.dot(SU2Generators2x2(a), np.dot(Phi, SU2Generators2x2(b))))) * (np.dot(U, np.dot(Delta, UDagger))[a, b]) for a in range(3) for b in range(3)])
    V9 = -M2 * sum([np.trace(np.dot(DeltaDagger, np.dot(SU2Generators3x3(a), np.dot(Delta, SU2Generators3x3(b))))) * (np.dot(U, np.dot(Delta, UDagger))[a, b]) for a in range(3) for b in range(3)])

    return V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9


# AVOIDING WRONG MINIMA
# Ref.: https://arxiv.org/abs/2207.00142

# VScalarPotential(phip, phi0, xip, xi0, chipp, chip, chi0, vDelta, M1, M2, lam2, lam3, lam4, lam5)

# IN THE CP-EVEN SUBSPACE
def V0p1p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, +cmath.sqrt(-lam1 * mu2Phi) / (2 * lam1), 0, 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p1m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, -cmath.sqrt(-lam1 * mu2Phi) / (2 * lam1), 0, 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p2p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, 0, 0, +cmath.sqrt(-2 * mu2Delta * (2 * lam4 + lam3)) / (2 * (2 * lam4 + lam3)), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p2m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, 0, 0, -cmath.sqrt(-2 * mu2Delta * (2 * lam4 + lam3)) / (2 * (2 * lam4 + lam3)), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p3p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, +cmath.sqrt(-mu2Delta * (lam4 + lam3)) / (2 * (lam4 + lam3)), 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p3m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, -cmath.sqrt(-mu2Delta * (lam4 + lam3)) / (2 * (lam4 + lam3)), 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p4(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, - (3 * M2) / (2 * lam3), 0, 0, (1 / lam3) * cmath.sqrt((-mu2Delta * lam3 ** 2 - 9 * M2 ** 2 * lam3 - 9 * M2 * lam4) / (2 * lam3 + 4 * lam4)), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p5p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, (3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, +cmath.sqrt(3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p5m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, (3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, -cmath.sqrt(3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p6p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, (3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, +cmath.sqrt(-3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p6m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, (3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, -cmath.sqrt(-3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p7p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, -(-3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, +cmath.sqrt(3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p7m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, -(-3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, -cmath.sqrt(3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p8p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, -(-3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, +cmath.sqrt(-3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0p8m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, -(-3 * M2 + cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2)) / (4 * lam3 + 12 * lam4), 0, 0, -cmath.sqrt(-3 * M2 * cmath.sqrt(-4 * mu2Delta * lam3 - 12 * mu2Delta * lam4 + 9 * M2 ** 2) - 2 * mu2Delta * lam3 - 6 * mu2Delta * lam4 + 9 * M2 ** 2) / (2 * lam3 + 6 * lam4), vDelta, M1, M2, lam2, lam3, lam4, lam5)


# IN THE CP-ODD SUBSPACE
def V0m1p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, +1j * cmath.sqrt(-lam1 * mu2Phi) / (2 * lam1), 0, 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0m1m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, -1j * cmath.sqrt(-lam1 * mu2Phi) / (2 * lam1), 0, 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0m2p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, 0, +1j * cmath.sqrt(-2 * mu2Delta * (2 * lam4 + lam3)) / (2 * (2 * lam4 + lam3)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0m2m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, 0, -1j * cmath.sqrt(-2 * mu2Delta * (2 * lam4 + lam3)) / (2 * (2 * lam4 + lam3)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def V0m3(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 1j * cmath.sqrt((-8 * mu2Phi * lam3 - 16 * mu2Phi * lam4 + 8 * mu2Delta * lam2 - 2 * mu2Delta * lam5) / (32 * lam1 * lam3 + 64 * lam1 * lam4 - 16 * lam2 ** 2 + 8 * lam2 * lam5 - lam5 ** 2)), 0, 0, 0, 1j * cmath.sqrt((8 * mu2Phi * lam2 - 2 * mu2Phi * lam5 - 16 * mu2Delta * lam1) / (32 * lam1 * lam3 + 64 * lam1 * lam4 - 16 * lam2 ** 2 + 8 * lam2 * lam5 - lam5 ** 2)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


# IN THE SINGLY-CHARGED SUBSPACE
def Vpm1p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(+cmath.sqrt(-2 * lam1 * mu2Phi) / (4 * lam1), 0, 0, 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm1m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(-cmath.sqrt(-2 * lam1 * mu2Phi) / (4 * lam1), 0, 0, 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm2p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, 0, +cmath.sqrt(-2 * mu2Delta * (lam4 + lam3)) / (4 * (lam4 + lam3)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm2m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, 0, -cmath.sqrt(-2 * mu2Delta * (lam4 + lam3)) / (4 * (lam4 + lam3)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm3p(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, +cmath.sqrt(-2 * mu2Delta * (lam4 + lam3)) / (4 * (lam4 + lam3)), 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm3m(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, +cmath.sqrt(-2 * mu2Delta * (lam4 + lam3)) / (4 * (lam4 + lam3)), 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm4(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(cmath.sqrt((lam2 * mu2Phi - 2 * lam1 * mu2Delta) / (16 * lam1 * lam3 + 16 * lam1 * lam4 - 4 * lam2 ** 2)), 0, 0, 0, 0, cmath.sqrt((-2 * lam3 * mu2Phi - 2 * lam4 * mu2Phi + lam2 * mu2Delta) / (16 * lam1 * lam3 + 16 * lam1 * lam4 - 4 * lam2 ** 2)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm5(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(cmath.sqrt((-mu2Delta) / (16 * lam4 + 8 * lam3)), 0, cmath.sqrt((-mu2Delta) / (16 * lam4 + 8 * lam3)), 0, 0, 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


def Vpm6(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    lam1 = lambda1(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    mu2Delta = musquareDelta(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, cmath.sqrt((lam2 * mu2Phi - 2 * lam1 * mu2Delta) / (16 * lam1 * lam3 + 16 * lam1 * lam4 - 4 * lam2 ** 2)), 0, 0, cmath.sqrt((-2 * lam3 * mu2Phi - 2 * lam4 * mu2Phi + lam2 * mu2Delta) / (16 * lam1 * lam3 + 16 * lam1 * lam4 - 4 * lam2 ** 2)), 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)


# IN THE DOUBLY-CHARGED SUBSPACE
def Vppmm(vDelta, M1, M2, lam2, lam3, lam4, lam5):
    mu2Phi = musquarePhi(vDelta, M1, M2, lam2, lam3, lam4, lam5)
    return VScalarPotential(0, 0, 0, 0, cmath.sqrt(-mu2Phi * (2 * lam4 + lam3)) / (2 * (2 * lam4 + lam3)), 0, 0, vDelta, M1, M2, lam2, lam3, lam4, lam5)





