#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Generate random scatter plot with main decay channel in color
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> October 25, 2022.
"""

###########
# MODULES #
###########

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cls
import matplotlib.patches as mpatches

##################
# IMPORTING DATA #
##################

# Directory to the .csv with random points
DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterData20000_GMSEESAW.csv'

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

# USE THIS TO KNOW EVERY DIFFERENT MAIN DECAY CHANNEL
# print(set(df['H3p_maindecaychannel']))
# print(len(set(df['H3p_maindecaychannel'])))

# COLOR ASSIGNMENT TO DECAY CHANNELS

colors = {
        "['4', '-3']": 'purple',
        "['6', '-5']": 'green',
        "['25', '24']": 'blue',
        "['252', '24']": 'darkturquoise',
        "['255', '-24']": 'red',
        "['-5', '5', '24']": 'pink',
        "['-1', '2', '25']": 'orange',
        "['-3', '4', '25']": 'peru',
        "['-1', '2', '252']": 'brown',
        "['-3', '4', '252']": 'olive',
        "['-2', '1', '255']": 'lime',
        "['-4', '3', '255']": 'magenta',
        "['-24', '256', '256']": 'black'}


#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['m3 [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H3p [cm]'],  # second variable scattered over the y-axis
        c=df['H3p_maindecaychannel'].map(colors),  # third variable represented as a color gradient
        marker='.',  # type of marker for the scattered points
        s=1,  # size of the markers
        alpha=1.0  # transparency of the markers
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

# LEGEND
purple_patch = mpatches.Patch(color='purple', label=r'$c\bar{s}$')
green_patch = mpatches.Patch(color='green', label=r'$t\bar{b}$')
blue_patch = mpatches.Patch(color='blue', label=r'$Zw^{+}$')
darkturquoise_patch = mpatches.Patch(color='darkturquoise', label=r'$HW^{+}$')
red_patch = mpatches.Patch(color='red', label=r'$H^{++}_5W^{-}$')
pink_patch = mpatches.Patch(color='pink', label=r'$\bar{b}bW^{+}$')
orange_patch = mpatches.Patch(color='orange', label=r'$\bar{d}uh$')
peru_patch = mpatches.Patch(color='peru', label=r'$\bar{s}ch$')
brown_patch = mpatches.Patch(color='brown', label=r'$\bar{d}uH$')
olive_patch = mpatches.Patch(color='olive', label=r'$\bar{s}cH$')
lime_patch = mpatches.Patch(color='lime', label=r'$\bar{u}dH^{++}_5$')
magenta_patch = mpatches.Patch(color='magenta', label=r'$\bar{c}sH^{++}_5$')
black_patch = mpatches.Patch(color='black', label=r'$W^{-}H^{+}_5H^{+}_5$')

ax.legend(handles=[
        purple_patch,
        green_patch,
        blue_patch,
        darkturquoise_patch,
        red_patch,
        pink_patch,
        orange_patch,
        peru_patch,
        brown_patch,
        olive_patch,
        lime_patch,
        magenta_patch,
        black_patch], loc='upper right', ncol=3, fontsize='small')

# AXES LABELS
ax.set_xlabel(r'$m_3$  [GeV]')
ax.set_ylabel(r'$c\tau(H^{+}_3)$  [cm]')

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
ax.set_title('MAIN DECAY CHANNELS\n' + 'GM MODEL + TYPE-II SEESAW (N=20.000)')

# SAVE PLOT
plt.savefig('ctau_H3p_wMDC_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()
