###########
# MODULES #
###########

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cls

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
        x=df['vDelta [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H5pp [cm]'],  # second variable scattered over the y-axis
        c=df['m5 [GeV]'],  # third variable represented as a color gradient
        cmap='jet',  # colormap for the third variable
        marker='.',  # type of marker for the scattered points
        s=1,  # size of the markers
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
cbar.set_label(r'$m_5$  [GeV]')

# AXES LABELS
ax.set_xlabel(r'$v_\Delta$  [GeV]')
ax.set_ylabel(r'$\Gamma(H^{++}_5)$  [cm]')
# ax.set_zlabel('z')

# AXES LIMITS
# ax.set_xlim(left=-0.5, right=0.5)
# ax.set_ylim(bottom=-0.5, top=0.5)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
ax.set_xscale('log')
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
ax.set_title(r'GM MODEL + TYPE-II SEESAW (N=100.000)')

# SAVE PLOT
plt.savefig('vDelta_gammaH5pp_GMSEESAW.png',
            bbox_inches='tight',
            dpi=100)

# SHOW PLOT
plt.show()
