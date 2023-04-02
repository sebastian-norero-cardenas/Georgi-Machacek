#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Generate random scatter plots
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> October 21, 2022.
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

# Absolute path to the MasterData .csv file
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
        x=df['mH [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H [cm]'],  # second variable scattered over the y-axis
        c=df['vDelta [GeV]'],  # third variable represented as a color gradient
        cmap='jet',  # colormap for the third variable
        marker='+',  # type of marker for the scattered points
        s=0.1,  # size of the markers
        alpha=1.0,  # transparency of the markers
        norm=cls.LogNorm()  # color bar in logarithmic scale
)

# m5=mh line
plt.axvline(x=1.2518e+02, color="black", linestyle="--", linewidth=0.5)
ax.text(x=1.2518e+02*(1 - 0.15), y=1e-15, s=r'$m_H = m_h$', color="black", fontsize=7, rotation=90)

# m5=2mW line
plt.axvline(x=2*7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*7.982436e+01*(1 + 0.05), y=1e-18, s=r'$m_H = 2m_W$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label(r'$v_\Delta$  [GeV]')

# AXES LABELS
ax.set_xlabel(r'$m_H$  [GeV]')
ax.set_ylabel(r'$c\tau(H)$  [cm]')

# AXES LIMITS
ax.set_xlim(left=50, right=750)
ax.set_ylim(bottom=1e-19, top=1e+1)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
ax.set_yscale('log')

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
ax.set_title('GM MODEL + TYPE-II SEESAW\n' + 'RANDOM SCAN; N=100.000')

# SAVE PLOT
plt.savefig('mH_ctauH_vDelta_GMSEESAW.png',
        bbox_inches='tight',
        dpi=100)

# SHOW PLOT
plt.show()