#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Generate random scatter plots
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> October 31, 2022.
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
        x=df['m5 [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H5pp [cm]'],  # second variable scattered over the y-axis
        # c=df['vDelta [GeV]'],  # third variable represented as a color gradient
        # cmap='jet',  # colormap for the third variable
        marker='+',  # type of marker for the scattered points
        s=0.1,  # size of the markers
        alpha=1.0,  # transparency of the markers
        # norm=cls.LogNorm()  # color bar in logarithmic scale
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
#         xy=(98.100759, 5.18497842207607), xytext=(70, 0),
#         textcoords='offset points', ha='right', va='bottom',
#         arrowprops=dict(arrowstyle='->', color='black'))

# # BENCHMARK POINT 3
# plt.plot(145.420659, 20.8427105790815, color='blue', marker='.')
# plt.annotate(
#         text='BP 3',
#         xy=(145.420659, 20.8427105790815), xytext=(50, 20),
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
ax.text(x=mW*(1 - 0.22), y=1e+00, s=r'$m_5 = m_W$', color="black", fontsize=7, rotation=90)

# m5 = mh line
plt.axvline(x=mh, color="black", linestyle="--", linewidth=0.5)
ax.text(x=mh*(1 - 0.20), y=1e+01, s=r'$m_5 = m_h$', color="black", fontsize=7, rotation=90)

# m5 = 2mW line
plt.axvline(x=2*mW, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*mW*(1 - 0.12), y=1e+04, s=r'$m_5 = 2m_W$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# COLOR BAR
# cbar = fig.colorbar(sc, ax=ax)
# cbar.set_label(r'$v_\Delta$  [GeV]')

# AXES LABELS
ax.set_xlabel(r'$m_5$  [GeV]')
ax.set_ylabel(r'$c\tau(H^{++}_5)$  [cm]')

# AXES LIMITS
#ax.set_xlim(left=80, right=712)
#ax.set_ylim(bottom=-250000, top=20000.0)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
ax.set_yscale('log')

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
ax.set_title('GM MODEL + TYPE-II SEESAW\n' + '(RANDOM SCAN; N=100.000)')

# SAVE PLOT
plt.savefig('BENCHMARK_POINTS_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()