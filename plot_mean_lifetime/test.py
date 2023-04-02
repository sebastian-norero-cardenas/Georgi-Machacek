#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Generate random scatter plots
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> Oct. 18, 2022.
"""

###########
# MODULES #
###########

# import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##################
# IMPORTING DATA #
##################

# Directory to the .csv with random points
DirRandomScan = '/home/sebastian/Proyects/Thesis/random_scans_madwidth/gather_decay_width/MasterData_10000P.csv'

df = pd.read_csv(DirRandomScan)

#####################
# GLOBAL PLOT STYLE #
#####################

plt.style.use('ggplot')

global_params = {
        'font.family': "sans",
        'legend.fontsize': 'medium',
        # 'figure.figsize': (15, 5),
        'axes.labelsize': 'medium',
        'axes.titlesize':'medium',
        'xtick.labelsize':'medium',
        'ytick.labelsize':'medium'
        }

plt.rcParams.update(global_params)

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['lam4'],  # first variable scattered over the x-axis
        y=df['ctau_H5pp'],  # second variable scattered over the y-axis
        c=df['vDelta'],  # third variable represented as a color gradient
        cmap='jet',  # colormap for the third variable
        marker='.',  # type of marker for the scattered points
        s=5,  # size of the markers
        alpha=1.0  # transparency of the markers
        # norm=cls.LogNorm()  # color bar in logarithmic scale
)
"""
# m3=m5 line
ax.plot(
        np.arange(80, 712, 0.1),
        np.arange(80, 712, 0.1),
        linestyle='--',
        color='black'
        )
"""

######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label(r'$v_\Delta$  [GeV]')

# AXES LABELS
ax.set_xlabel(r'$M_1$')
ax.set_ylabel(r'$c\tau$  [cm]')

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
#ax.set_title(r'$H \to\ ZZ$')

# SAVE PLOT
plt.savefig('test.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()