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
# print(set(df['H_maindecaychannel']))
# print(len(set(df['H_maindecaychannel'])))

# COLOR ASSIGNMENT TO DECAY CHANNELS

colors = {
        "['5', '-5']": 'lime',
        "['6', '-6']": 'darkturquoise',
        "['-24', '24']": 'blue',
        "['253', '-24']": 'purple',
        "['25', '25']": 'brown',
        "['-253', '253']": 'red',
        "['-255', '255']": 'orange',
        "['-24', '-1', '2']": 'green',
        "['-24', '-3', '4']": 'olive' ,
        "['-2', '1', '24']": 'black',
        "['-4', '3', '24']": 'peru',
        "['-255', '24', '24']": 'pink',
        "['-2', '1', '253']": 'magenta',
        "['-4', '3', '253']": 'navy',
        "['-24', '-24', '255']": 'tan',
        "['23', '23', '257']": 'coral'}

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['mH [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H [cm]'],  # second variable scattered over the y-axis
        c=df['H_maindecaychannel'].map(colors),  # third variable represented as a color gradient
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
lime_patch = mpatches.Patch(color='lime', label=r'$b\bar{b}$')
darkturquoise_patch = mpatches.Patch(color='darkturquoise', label=r'$t\bar{t}$')
blue_patch = mpatches.Patch(color='blue', label=r'$W^{-}W^{+}$')
purple_patch = mpatches.Patch(color='purple', label=r'$H^{+}_3W^{-}$')
brown_patch = mpatches.Patch(color='brown', label=r'$hh$')
red_patch = mpatches.Patch(color='red', label=r'$H^{+}_3H^{-}_3$')
orange_patch = mpatches.Patch(color='orange', label=r'$H^{++}_5H^{--}_5$')
green_patch = mpatches.Patch(color='green', label=r'$\bar{d}uW^{-}$')
olive_patch = mpatches.Patch(color='olive', label=r'$\bar{s}cW^{-}$')
black_patch = mpatches.Patch(color='black', label=r'$\bar{d}uW^{+}$')
peru_patch = mpatches.Patch(color='peru', label=r'$\bar{c}sW^{+}$')
pink_patch= mpatches.Patch(color='pink', label=r'$H^{--}_5W^{+}W^{+}$')
magenta_patch = mpatches.Patch(color='magenta', label=r'$\bar{u}dH^{+}_3$')
navy_patch = mpatches.Patch(color='navy', label=r'$\bar{c}sH^{+}_3$')
tan_patch = mpatches.Patch(color='tan', label=r'$W^{-}W^{-}H^{++}_5$')
coral_patch = mpatches.Patch(color='coral', label=r'$ZZH^{0}_5$')

ax.legend(handles=[
        lime_patch,
        darkturquoise_patch,
        blue_patch,
        purple_patch,
        brown_patch,
        red_patch,
        orange_patch,
        green_patch,
        olive_patch,
        black_patch,
        peru_patch,
        pink_patch,
        magenta_patch,
        navy_patch,
        tan_patch,
        coral_patch], loc='upper right', ncol=3, fontsize='small')

# AXES LABELS
ax.set_xlabel(r'$m_H$  [GeV]')
ax.set_ylabel(r'$c\tau(H)$  [cm]')

# AXES LIMITS
# ax.set_xlim(left=80, right=712)
# ax.set_ylim(bottom=3e-18, top=1e-04)

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
plt.savefig('ctau_H_wMDC_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()
