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
INFO: Last update -> November 11, 2022.
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
ee = (4 * np.pi) / 127.9
v = 1 / np.sqrt(np.sqrt(2) * GF)  # SM Higgs vev ~ 246.22 [GeV]
mW = np.sqrt(mZ ** 2 / 2 + np.sqrt((mZ ** 2 / 2) ** 2 - (np.pi / np.sqrt(2)) * (alphaEW / GF) * mZ ** 2))  # W boson mass
sW = np.sqrt(1 - (mW / mZ) ** 2)

# Internal parameters
r1 = 6  # round to r1 decimals
r2 = 16  # round to r2 decimals (used for neutrino yukawas and other very small parameters)

n = 4042  # param_dard_n

################################
# INDEPENDENT PARAMETER RANGES #
################################

sA = float(df.iloc[[n-11.40911338343]]['sA'])
tH = float(df.iloc[[n-1]]['tH'])
M1 = float(df.iloc[[n-1]]['M1 [GeV]'])
M2 = float(df.iloc[[n-1]]['M2 [GeV]'])
lam2 = float(df.iloc[[n-1]]['lam2'])
lam3 = float(df.iloc[[n-1]]['lam3'])
lam4 = float(df.iloc[[n-1]]['lam4'])
lam5 = float(df.iloc[[n-1]]['lam5'])
vphi = float(df.iloc[[n-1]]['vphi [GeV]'])
vDelta = float(df.iloc[[n-1]]['vDelta [GeV]'])
m3 = float(df.iloc[[n-1]]['m3 [GeV]'])
m5 = float(df.iloc[[n-1]]['m5 [GeV]'])
mh = float(df.iloc[[n-1]]['mh [GeV]'])
mH = float(df.iloc[[n-1]]['mH [GeV]'])
mu2Phi = float(df.iloc[[n-1]]['mu2Phi [GeV^2]'])
mu2Delta = float(df.iloc[[n-1]]['mu2Delta [GeV^2]'])
Dm21 = float(df.iloc[[n-1]]['Dm21 [eV^2]'] * 1E-18)
Dm31 = float(df.iloc[[n-1]]['Dm31 [eV^2]'] * 1E-18)
PMNSs12 = float(df.iloc[[n-1]]['PMNSs12'])
PMNSs23 = float(df.iloc[[n-1]]['PMNSs23'])
PMNSs13 = float(df.iloc[[n-1]]['PMNSs13'])
PMNSphase = float(df.iloc[[n-1]]['PMNSphase'])
mv1 = float(df.iloc[[n-1]]['mv1 [eV]'] * 1E-9)
mv2 = float(df.iloc[[n-1]]['mv2 [eV]'] * 1E-9)
mv3 = float(df.iloc[[n-1]]['mv3 [eV]'] * 1E-9)
CKMs12 = float(df.iloc[[n-1]]['CKMs12'])
CKMs23 = float(df.iloc[[n-1]]['CKMs23'])
CKMs13 = float(df.iloc[[n-1]]['CKMs13'])
CKMphase = float(df.iloc[[n-1]]['CKMphase'])

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


VCKM11 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[0,0]
VCKM12 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[0, 1]
VCKM13 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[0, 2]
VCKM21 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[1, 0]
VCKM22 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[1, 1]
VCKM23 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[1, 2]
VCKM31 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[2, 0]
VCKM32 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[2, 1]
VCKM33 = CKM_MATRIX(CKMs12, CKMs23, CKMs13, CKMphase)[2, 2]
VPMNS11 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[0,0]
VPMNS12 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[0,1]
VPMNS13 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[0,2]
VPMNS21 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[1,0]
VPMNS22 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[1,1]
VPMNS23 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[1,2]
VPMNS31 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[2,0]
VPMNS32 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[2,1]
VPMNS33 = PMNS_MATRIX(PMNSs12, PMNSs23, PMNSs13, PMNSphase)[2,2]

YvND11 = float(df.iloc[[n-1]]['YvND11'])
YvND12 = float(df.iloc[[n-1]]['YvND12'])
YvND13 = float(df.iloc[[n-1]]['YvND13'])
YvND21 = float(df.iloc[[n-1]]['YvND21'])
YvND22 = float(df.iloc[[n-1]]['YvND22'])
YvND23 = float(df.iloc[[n-1]]['YvND23'])
YvND31 = float(df.iloc[[n-1]]['YvND31'])
YvND32 = float(df.iloc[[n-1]]['YvND32'])
YvND33 = float(df.iloc[[n-1]]['YvND33'])
yv1 = float(df.iloc[[n-1]]['yv1'])
yv2 = float(df.iloc[[n-1]]['yv2'])
yv3 = float(df.iloc[[n-1]]['yv3'])
yelectron = float(df.iloc[[n-1]]['yelectron'])
ymuon = float(df.iloc[[n-1]]['ymuon'])
ytau = float(df.iloc[[n-1]]['ytau'])
yup = float(df.iloc[[n-1]]['yup'])
ycharm = float(df.iloc[[n-1]]['ycharm'])
ytop = float(df.iloc[[n-1]]['ytop'])
ydown = float(df.iloc[[n-1]]['ydown'])
ystrange = float(df.iloc[[n-1]]['ystrange'])
ybottom = float(df.iloc[[n-1]]['ybottom'])


#####################
# DECAY WIDTH: H5pp #
#####################

# H5pp -> e+ e+
def GammaH5ppTOee():
    f1 = (1/2 * 2) / (8 * np.pi)
    g = (mv1 * VPMNS11 * np.conjugate(VPMNS11) + mv2 * VPMNS12 * np.conjugate(VPMNS12) + mv3 * VPMNS13 * np.conjugate(VPMNS13)) / (2 * vDelta)
    f2 = cmath.sqrt(1 - (4 * me ** 2) /m5 ** 2)
    h = np.heaviside(m5 - 2 * me, 0)
    return np.real(f1 * m5 * np.abs(g) ** 2 * f2 * h)

# H5pp -> e+ mu+
def GammaH5ppTOemu():
    f1 = (1 * 2) / (8 * np.pi)
    g = (mv1 * VPMNS11 * np.conjugate(VPMNS21) + mv2 * VPMNS12 * np.conjugate(VPMNS22) + mv3 * VPMNS13 * np.conjugate(VPMNS23)) / (2 * vDelta)
    f2 = 1 - ((me - mmu) / m5) ** 2
    f3 = cmath.sqrt(1 - 2 * ((me ** 2 + mmu ** 2) / m5 ** 2) + ((me ** 2 - mmu ** 2) / m5 ** 2) ** 2)
    h = np.heaviside(m5 - (me + mmu), 0)
    return np.real(f1 * m5 * np.abs(g) ** 2 * f2 * f3 * h)

# H5pp -> e+ tau+
def GammaH5ppTOetau():
    f1 = (1 * 2) / (8 * np.pi)
    g = (mv1 * VPMNS11 * np.conjugate(VPMNS31) + mv2 * VPMNS12 * np.conjugate(VPMNS32) + mv3 * VPMNS13 * np.conjugate(VPMNS33)) / (2 * vDelta)
    f2 = 1 - ((me - mtau) / m5) ** 2
    f3 = cmath.sqrt(1 - 2 * ((me ** 2 + mtau ** 2) / m5 ** 2) + ((me ** 2 - mtau ** 2) / m5 ** 2) ** 2)
    h = np.heaviside(m5 - (me + mtau), 0)
    return np.real(f1 * m5 * np.abs(g) ** 2 * f2 * f3 * h)

# H5pp -> mu+ mu+
def GammaH5ppTOmumu():
    f1 = (1/2 * 2) / (8 * np.pi)
    g = (mv1 * VPMNS21 * np.conjugate(VPMNS21) + mv2 * VPMNS22 * np.conjugate(VPMNS22) + mv3 * VPMNS23 * np.conjugate(VPMNS23)) / (2 * vDelta)
    f2 = cmath.sqrt(1 - (4 * mmu ** 2) /m5 ** 2)
    h = np.heaviside(m5 - 2 * mmu, 0)
    return np.real(f1 * m5 * np.abs(g) ** 2 * f2 * h)

# H5pp -> mu+ tau+
def GammaH5ppTOmutau():
    f1 = (1 * 2) / (8 * np.pi)
    g = (mv1 * VPMNS21 * np.conjugate(VPMNS31) + mv2 * VPMNS22 * np.conjugate(VPMNS32) + mv3 * VPMNS23 * np.conjugate(VPMNS33)) / (2 * vDelta)
    f2 = 1 - ((mmu - mtau) / m5) ** 2
    f3 = cmath.sqrt(1 - 2 * ((mmu ** 2 + mtau ** 2) / m5 ** 2) + ((mmu ** 2 - mtau ** 2) / m5 ** 2) ** 2)
    h = np.heaviside(m5 - (mmu + mtau), 0)
    return np.real(f1 * m5 * np.abs(g) ** 2 * f2 * f3 * h)

# H5pp -> tau+ tau+
def GammaH5ppTOtautau():
    f1 = (1/2 * 2) / (8 * np.pi)
    g = (mv1 * VPMNS31 * np.conjugate(VPMNS31) + mv2 * VPMNS32 * np.conjugate(VPMNS32) + mv3 * VPMNS33 * np.conjugate(VPMNS33)) / (2 * vDelta)
    f2 = cmath.sqrt(1 - (4 * mtau ** 2) /m5 ** 2)
    h = np.heaviside(m5 - 2 * mtau, 0)
    return np.real(f1 * m5 * np.abs(g) ** 2 * f2 * h)

# H5pp -> W+ W+
def GammaH5ppTOWW():
    f1 = 1 / (32 * np.pi)
    f2 = 1 / m5
    g = (2 * ee * vDelta) / sW ** 2
    f3 = 3 - m5 ** 2 / mW ** 2 + m5 ** 4 / (4 * mW ** 4)
    f4 = cmath.sqrt(1 - (4 * mW ** 2) /m5 ** 2)
    h = np.heaviside(m5 - 2 * mW, 0)
    return np.real(f1 * f2 * np.abs(g) ** 2 * f3 * f4 * h)

# H5pp -> H3+ H3+
def GammaH5ppTOH3pH3p():
    f1 = 1 / (32 * np.pi)
    f2 = 1 / m5
    g = (-2 / v ** 2) * (4 * M1 * vDelta ** 2 + 3 * M2 * vphi ** 2 + 2 * vDelta * vphi ** 2 * lam3 - 8 * vDelta ** 3 * lam5 - 4 * vDelta * vphi ** 2 * lam5)
    f3 = cmath.sqrt(1 - (4 * m3 ** 2) /m5 ** 2)
    h = np.heaviside(m5 - 2 * m3, 0)
    return np.real(f1 * f2 * np.abs(g) ** 2 * f3 * h)

# H5pp -> H3+ W+
def GammaH5ppTOH3pW():
    f1 = 1 / (16 * np.pi)
    f2 = mW ** 2 / m5
    g = (ee * vphi) / (np.sqrt(2) * sW * v)
    f3 = 1 - 2 * ((m5 ** 2 + m3 ** 2) / mW ** 2) + ((m5 ** 2 - m3 ** 2) / mW ** 2) ** 2
    f4 = cmath.sqrt(1 - 2 * ((m3 ** 2 + mW ** 2) / m5 ** 2) + ((m3 ** 2 - mW ** 2) / m5 ** 2) ** 2)
    h = np.heaviside(m5 - (m3 + mW), 0)
    return np.real(f1 * f2 * np.abs(g) ** 2 * f3 * f4 * h)


def MadWidth():
    ch = vphi / v
    MH3 = m3
    MH5 = m5
    sw = sW
    MW = mW
    vchi = vDelta
    M1coeff = M1
    M2coeff = M1
    sh = np.sqrt(1 - ch ** 2)
    return np.real((((12*ee**4*vchi**2)/sw**4 + (ee**4*MH5**4*vchi**2)/(MW**4*sw**4) - (4*ee**4*MH5**2*vchi**2)/(MW**2*sw**4))*cmath.sqrt(MH5**4 - 4*MH5**2*MW**2))/(32.*cmath.pi*abs(MH5)**3))


#################
# PRINT RESULTS #
#################
gamma_unicode = "\u0393"

print('--- MASSES ---')
print(f'm3 = {m3}')
print(f'm5 = {m5}')
print(f'mH = {mH}')
print('--- DECAYS WIDTHS (THEORETICAL VALUES) ---')
print(f'{gamma_unicode}(H5++ -> e+ e+): {GammaH5ppTOee()} GeV')
print(f'{gamma_unicode}(H5++ -> e+ mu+): {GammaH5ppTOemu()} GeV')
print(f'{gamma_unicode}(H5++ -> e+ tau+): {GammaH5ppTOetau()} GeV')
print(f'{gamma_unicode}(H5++ -> mu+ mu+): {GammaH5ppTOmumu()} GeV')
print(f'{gamma_unicode}(H5++ -> mu+ tau+): {GammaH5ppTOmutau()} GeV')
print(f'{gamma_unicode}(H5++ -> tau+ tau+): {GammaH5ppTOtautau()} GeV')
print('---')
print(f'{gamma_unicode}(H5++ -> W+ W+): {GammaH5ppTOWW()} GeV')
print(f'{gamma_unicode}(H5++ -> W+ W+): {MadWidth()} GeV (MadWidth)')
print('---')
print(f'{gamma_unicode}(H5++ -> H3+ H3+): {GammaH5ppTOH3pH3p()} GeV')
print('---')
print(f'{gamma_unicode}(H5++ -> H3+ W+): {GammaH5ppTOH3pW()} GeV')