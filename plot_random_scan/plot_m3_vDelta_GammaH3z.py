###########
# MODULES #
###########

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cls

#######################
# EXTERNAL PARAMETERS #
#######################

mh = 125.18

##################
# IMPORTING DATA #
##################

# Directory to the .csv with random points
DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW_LOGSCAN/MasterData100000_GMSEESAW_LOGSCAN.csv'

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
        x=df['m3 [GeV]'],  # first variable scattered over the x-axis
        y=df['vDelta [GeV]'],  # second variable scattered over the y-axis
        c=df['ctau_H3z [cm]'],  # third variable represented as a color gradient
        cmap='jet',  # colormap for the third variable
        marker='.',  # type of marker for the scattered points
        s=1,  # size of the markers
        alpha=1.0,  # transparency of the markers
        norm=cls.LogNorm()  # color bar in logarithmic scale
)

# m5=mh line
plt.axvline(x=mh, color="black", linestyle="--", linewidth=0.5)
ax.text(x=mh*(1 - 0.20), y=5e-07, s=r'$m_5 = m_h$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label(r'$c\tau(H^{0}_3)$  [cm]')

# AXES LABELS
ax.set_xlabel(r'$m_3$  [GeV]')
ax.set_ylabel(r'$v_\Delta$  [GeV]')
# ax.set_zlabel('z')

# AXES LIMITS
# ax.set_xlim(left=-0.5, right=0.5)
# ax.set_ylim(bottom=-0.5, top=0.5)

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
ax.set_title('GM MODEL + TYPE-II SEESAW\n' + '(LOGARITHMIC SCAN; N=100.000)')

# SAVE PLOT
plt.savefig('m3_vDelta_gammaH3z_GMSEESAW.png',
            bbox_inches='tight',
            dpi=100)

# SHOW PLOT
plt.show()
