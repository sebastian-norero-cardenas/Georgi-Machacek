#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Plot random scan
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> Dcember 30, 2022.
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
        x=df['m5 [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H5p [cm]'],  # second variable scattered over the y-axis
        c=df['vDelta [GeV]'],  # third variable represented as a color gradient
        cmap='jet',  # colormap for the third variable
        marker='+',  # type of marker for the scattered points
        s=0.1,  # size of the markers
        alpha=1.0,  # transparency of the markers
        norm=cls.LogNorm()  # color bar in logarithmic scale
)

# m5=mW line
plt.axvline(x=7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=7.982436e+01*(1 - 0.30), y=1e+00, s=r'$m_5 = m_W$', color="black", fontsize=7, rotation=90)

# m5=mh line
plt.axvline(x=1.2518e+02, color="black", linestyle="--", linewidth=0.5)
ax.text(x=1.2518e+02*(1 - 0.12), y=1e+02, s=r'$m_5 = m_h$', color="black", fontsize=7, rotation=90)

# m5=mW+mZ line
plt.axvline(x=7.982436e+01 + 9.11876e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=(7.982436e+01 + 9.11876e+01)*(1 - 0.11), y=1e+03, s=r'$m_5 = m_W + m_Z$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label(r'$v_\Delta$  [GeV]')

# AXES LABELS
ax.set_xlabel(r'$m_5$  [GeV]')
ax.set_ylabel(r'$c\tau(H^{+}_5)$  [cm]')

# AXES LIMITS
#ax.set_xlim(left=80, right=712)
#ax.set_ylim(bottom=-250000, top=20000.0)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
ax.set_yscale('log')

# USE A GRAY BACKGROUND
# ax.set_facecolor('#E6E6E6')
# ax.set_axisbelow(True)

# GRID
# ax.grid(
#         True,
#         color='black',
#         linestyle='--'
# )

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
ax.set_title('GM MODEL + TYPE-II SEESAW\n' + 'RANDOM SCAN; N=100.000')

# SAVE PLOT
plt.savefig('m5_ctauH5p_vDelta_GMSEESAW.png',
        bbox_inches='tight',
        dpi=100)

# SHOW PLOT
plt.show()