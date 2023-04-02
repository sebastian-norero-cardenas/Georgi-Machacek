# -*- coding: utf-8 -*-
#!/usr/bin/env python3.8

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + type-II seesaw mechanism.
INFO: Intended    -> Module storing constants used in the random scan.
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

#############
# CONSTANTS #
#############

# Independent physical constants
# Ref.: R.L. Workman et al. (Particle Data Group), Prog. Theor. Exp. Phys. 2022, 083C01 (2022)
alpha_s = 0.1179  # Strong coupling at Q^2 = (M_Z)^2; Central value
alphaEW = 1 / 127.9  # EW coupling at Q^2 = (M_Z)^2
GF = 1.166378e-5  # Fermi constant [GeV]; Central value

# Z-pole observable Rb
# Ref.: https://arxiv.org/abs/1410.5538
RbSM = 0.21577  # Z-pole observable Rb: SM prediction; central value
RbExp = 0.21629  # Z-pole observable Rb: measured value; central value

# B0s-barB0s MIXING
# Ref.: https://arxiv.org/abs/1410.5538
RExpDeltam = 1.02  # See Eq.(38) from Ref.

# BR(B0s to mu+ mu-)
# Ref.: https://arxiv.org/abs/1410.5538
RExpsmu =  0.79  # See Eq.(44) from Ref.

# POLE MASSES
MZ = 91.1876  # Z boson mass [GeV]; Pole mass; Central value

Me = 5.11e-04  # Electron mass [GeV]; Pole mass; Central value
Mmu = 0.10566  # Muon mass [GeV]; Pole mass; Central value
Mtau = 1.777  # Tau mass [GeV]; Pole mass; Central value

Mup = 2.55e-03  # Up quark mass [GeV]; Pole mass; Central value
Mcharm = 1.42  # Charm quark mass [GeV]; Pole mass; Central value
Mtop = 172.69  # Top quark mass [GeV]; Pole mass; Central value
Mdown = 5.04e-03  # Down quark mass [GeV]; Pole mass; Central value
Mstrange = 0.101  # Strange quark mass [GeV]; Pole mass; Central value
Mbottom = 4.7  # Bottom quark mass [GeV]; Pole mass; Central value

Mh = 125.18  # 125-Higgs mass [GeV]; Pole mass; Central value

# RUNNING MASSES
mtop = 163.30  # Top quark mass [GeV]; SMbar running mass at Q^2 = (M_top)^2; Central value

# Dependent parameters
v = 1 / np.sqrt(np.sqrt(2) * GF)  # SM Higgs vev ~ 246.22 [GeV]
mW = np.sqrt(MZ ** 2 / 2 + np.sqrt((MZ ** 2 / 2) ** 2 - (np.pi / np.sqrt(2)) * (alphaEW / GF) * MZ ** 2))  # W boson mass
sW = np.sqrt(1 - (mW / MZ) ** 2)  # Sine of Weinberg angle
cW = mW / MZ  # Cosine of Weinberg angle
ee = np.sqrt(4 * np.pi * alphaEW)  # Unit of electric charge

# 'S' OBLIQUE PARAMETER
S_exp = 0.05  # experimental central value of the 'S' oblique parameter
T_exp = 0.10  # experimental central value of the 'T' oblique parameter
rho_ST = 0.91  # relative correlation between the 'S' and 'T' oblique parameters
DeltaS_exp = 0.09  # 1-sigma experimental uncertainty of the 'S' oblique parameter
DeltaT_exp = 0.07  # 1-sigma experimental uncertainty of the 'T' oblique parameter