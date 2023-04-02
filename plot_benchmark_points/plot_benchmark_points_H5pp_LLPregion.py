#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Scatter plots with benchmark points
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> December 27, 2022.
"""

###########
# MODULES #
###########

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cls

#######################
# EXTERNAL PARAMETERS #
#######################

mh = 125.18
mW = 79.82436
mZ = 91.1876

##################
# IMPORTING DATA #
##################

# Directory to the .csv with random points
DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterData100000_GMSEESAW.csv'

df = pd.read_csv(DirRandomScan)

#####################
# GLOBAL PLOT STYLE #
#####################

plt.style.use('ggplot')

global_params = {
        'font.family': "monospace",
        'legend.fontsize': 'medium',
        # 'figure.figsize': (15, 5),
        'axes.labelsize': 'medium',
        'axes.titlesize':'medium',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'}

plt.rcParams.update(global_params)

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['m5 [GeV]'],
        y=30 * df['ctau_H5pp [cm]'],  # <βγ> = 30
        # c=df['vDelta [GeV]'],
        # cmap='jet',
        marker='+',
        s=0.1,
        alpha=1.0
        # norm=cls.LogNorm()
)

# BENCHMARK POINT 1
plt.plot(145.420659, 20.8427105790815, color='blue', marker='.')
plt.annotate(
        text='BP 1',
        xy=(145.420659, 20.8427105790815), xytext=(50, 10),
        textcoords='offset points', ha='right', va='bottom',
        arrowprops=dict(arrowstyle='->', color='black'))

# # BENCHMARK POINT 2
# plt.plot(98.100759, 5.18497842207607, color='blue', marker='.')
# plt.annotate(
#         text='BP 2',
#         xy=(98.100759, 5.18497842207607), xytext=(70, 10),
#         textcoords='offset points', ha='right', va='bottom',
#         arrowprops=dict(arrowstyle='->', color='black'))

# # BENCHMARK POINT 3
# plt.plot(145.420659, 20.8427105790815, color='blue', marker='.')
# plt.annotate(
#         text='BP 3',
#         xy=(145.420659, 20.8427105790815), xytext=(50, 10),
#         textcoords='offset points', ha='right', va='bottom',
#         arrowprops=dict(arrowstyle='->', color='black'))

# # BENCHMARK POINT 4
# plt.plot(602.426507, 0.0000508477363830789, color='blue', marker='.')
# plt.annotate(
#         text='BP 4',
#         xy=(602.426507, 0.0000508477363830789), xytext=(50, 20),
#         textcoords='offset points', ha='right', va='bottom',
#         arrowprops=dict(arrowstyle='->', color='black'))

# -------------

# m5 = mW line
plt.axvline(x=mW, color="black", linestyle="--", linewidth=0.5)
ax.text(x=mW*(1 - 0.22), y=2e+03, s=r'$m_5 = m_W$', color="black", fontsize=7, rotation=90)

# m5 = mh line
plt.axvline(x=mh, color="black", linestyle="--", linewidth=0.5)
ax.text(x=mh*(1 - 0.15), y=1.8e+03, s=r'$m_5 = m_h$', color="black", fontsize=7, rotation=90)

# m5 = 2mW line
plt.axvline(x=2*mW, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*mW*(1 + 0.06), y=1.5e+03, s=r'$m_5 = 2m_W$', color="black", fontsize=7, rotation=90)

# -------------

plt.axhline(y=1100, color="black", linestyle="--", linewidth=0.5) # ctau = 11 [m] line

ax.text(x=450, y=620, s=r'MUON SPECTROMETER', color="black", fontsize=7, rotation=0)  # MUON SPECTROMETER

plt.axhline(y=450, color="black", linestyle="--", linewidth=0.5) # ctau = 4.5 [m] line

ax.text(x=450, y=240, s=r'HADRONIC CALORIMETER', color="black", fontsize=7, rotation=0)  # HCal

plt.axhline(y=155, color="black", linestyle="--", linewidth=0.5) # ctau = 1550 [mm] line

ax.text(x=450, y=60, s=r'ELECTROMAGNETIC CALORIMETER', color="black", fontsize=7, rotation=0)  # EMCal

plt.axhline(y=108.2, color="black", linestyle="--", linewidth=0.5)  # ctau = 1082 [mm] line

ax.text(x=450, y=6, s=r'INNER DETECTOR', color="black", fontsize=7, rotation=0)  # Inner detector

plt.axhline(y=3.3, color="black", linestyle="--", linewidth=0.5)  # ctau = 33 [mm] line

plt.axhline(y=0.2, color="black", linestyle="--", linewidth=0.5)  # ctau = 2 [mm] line
ax.text(x=450, y=3e-01, s=r'$|d_0| = 2$  [mm]', color="black", fontsize=7, rotation=0)  # |d0| > 2 [mm]

ax.axhspan(3.3, 108.2, facecolor='cyan', alpha=0.2)
ax.axhspan(108.2, 155, facecolor='green', alpha=0.2)
ax.axhspan(155, 450, facecolor='red', alpha=0.2)
ax.axhspan(450, 1100, facecolor='purple', alpha=0.2)

######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
# cbar = fig.colorbar(sc, ax=ax)
# cbar.set_label(r'$v_\Delta$  [GeV]')

# AXES LABELS
ax.set_xlabel(r'$m_5$  [GeV]')
ax.set_ylabel(r'$\beta\gamma c\tau(H^{++}_5)$  [cm]')

# AXES LIMITS
# ax.set_xlim(left=1e-05, right=712)
ax.set_ylim(bottom=4e-02, top=1e+04)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
ax.set_yscale('log')

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
ax.set_title('GM MODEL + TYPE-II SEESAW\n' + r'RANDOM SCAN; N=100.000; $\langle\beta\gamma\rangle = 30$')

# SAVE PLOT
plt.savefig('BENCHMARK_POINTS_LLPsREGION_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()