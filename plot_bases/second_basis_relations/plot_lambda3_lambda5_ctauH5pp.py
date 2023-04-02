###########
# MODULES #
###########

# import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cls

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
        x=df['lam3'],  # first variable scattered over the x-axis
        y=df['lam5'],  # second variable scattered over the y-axis
        c=df['ctau_H5pp [cm]'],  # third variable represented as a color gradient
        cmap='tab20b',  # colormap for the third variable
        marker='+',  # type of marker for the scattered points
        s=0.1,  # size of the markers
        alpha=1.0,  # transparency of the markers
        norm=cls.LogNorm()  # color bar in logarithmic scale
)


######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
cbar = fig.colorbar(sc, ax=ax)
cbar.set_label(r'$c\tau(H^{++}_5)$  [cm]')

# AXES LABELS
ax.set_xlabel(r'$\lambda_3$')
ax.set_ylabel(r'$\lambda_5$')
# ax.set_zlabel('z')

# AXES LIMITS
#ax.set_xlim(left=80, right=712)
#ax.set_ylim(bottom=-250000, top=20000.0)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
# ax.set_yscale('log')

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
plt.savefig('lambda3_lambda5_ctauH5pp_GMSEESAW.png',
            bbox_inches='tight',
            dpi=100)

# SHOW PLOT
plt.show()
