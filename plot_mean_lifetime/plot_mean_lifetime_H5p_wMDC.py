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
# print(set(df['H5p_maindecaychannel']))
# print(len(set(df['H5p_maindecaychannel'])))

# COLOR ASSIGNMENT TO DECAY CHANNELS
colors = {
        "['16', '-13']": 'blue',
        "['16', '-15']": 'orange',
        "['24', '23']": 'red',
        "['254', '24']": 'darkturquoise',
        "['253', '254']": 'hotpink',
        "['-1', '2', '23']": 'brown',
        "['-3', '4', '23']": 'purple',
        "['12', '12', '24']": 'green',
        "['14', '14', '24']": 'olive',
        "['16', '16', '24']": 'black',
        "['-1', '2', '254']": 'peru',
        "['-3', '4', '254']": 'sandybrown',
        "['-1', '-1', '1', '2']": 'tomato',
        "['-1', '2', '12', '12']": 'khaki',
        "['-3', '4', '12', '12']": 'palegreen',
        "['-1', '2', '14', '14']": 'teal',
        "['-3', '4', '14', '14']": 'steelblue',
        "['-1', '2', '16', '16']": 'mediumorchid',
        "['-3', '4', '16', '16']": 'plum'}

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['m5 [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H5p [cm]'],  # second variable scattered over the y-axis
        c=df['H5p_maindecaychannel'].map(colors),  # third variable represented as a color gradient
        marker='.',  # type of marker for the scattered points
        s=1,  # size of the markers
        alpha=1.0  # transparency of the markers
)

# m5=mW line
plt.axvline(x=7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=7.982436e+01*(1 - 0.22), y=1e+00, s=r'$m_5 = m_W$', color="black", fontsize=7, rotation=90)

# m5=mW+mZ line
plt.axvline(x=7.982436e+01 + 9.11876e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=(7.982436e+01 + 9.11876e+01)*(1 - 0.11), y=1e-02, s=r'$m_5 = m_W + m_Z$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# LEGEND
blue_patch = mpatches.Patch(color='blue', label=r'$\nu_\tau\mu^{+}$')
orange_patch = mpatches.Patch(color='orange', label=r'$\nu_\tau\tau^{+}$')
red_patch = mpatches.Patch(color='red', label=r'$W^{+}Z$')
darkturquoise_patch = mpatches.Patch(color='darkturquoise', label=r'$H^{0}_3W^{+}$')
hotpink_patch = mpatches.Patch(color='hotpink', label=r'$H^{+}_3H^{0}_3$')
brown_patch = mpatches.Patch(color='brown', label=r'$\bar{d}uZ$')
purple_patch = mpatches.Patch(color='purple', label=r'$\bar{s}cZ$')
green_patch = mpatches.Patch(color='green', label=r'$\nu_e\nu_eW^{+}$')
olive_patch = mpatches.Patch(color='olive', label=r'$\nu_\mu\nu_\mu W^{+}$')
black_patch = mpatches.Patch(color='black', label=r'$\nu_\tau\nu_\tau W^{+}$')
peru_patch = mpatches.Patch(color='peru', label=r'$\bar{d}uH^{0}_3$')
sandybrown_patch = mpatches.Patch(color='sandybrown', label=r'$\bar{s}cH^{0}_3$')
tomato_patch = mpatches.Patch(color='tomato', label=r'$\bar{d}\bar{d}du$')
khaki_patch = mpatches.Patch(color='khaki', label=r'$\nu_e\nu_e\bar{d}u$')
palegreen_patch = mpatches.Patch(color='palegreen', label=r'$\nu_e\nu_e\bar{s}c$')
teal_patch = mpatches.Patch(color='teal', label=r'$\nu_\mu\nu_\mu\bar{d}u$')
steelblue_patch = mpatches.Patch(color='steelblue', label=r'$\nu_\mu\nu_\mu\bar{s}c$')
mediumorchid_patch = mpatches.Patch(color='mediumorchid', label=r'$\nu_\tau\nu_\tau\bar{d}u$')
plum_patch = mpatches.Patch(color='plum', label=r'$\bar{s}c\nu_\tau\nu_\tau$')

ax.legend(handles=[
        red_patch,
        darkturquoise_patch,
        blue_patch,
        orange_patch,
        hotpink_patch,
        brown_patch,
        purple_patch,
        green_patch,
        olive_patch,
        black_patch,
        peru_patch,
        sandybrown_patch,
        tomato_patch,
        khaki_patch,
        palegreen_patch,
        teal_patch,
        steelblue_patch,
        mediumorchid_patch,
        plum_patch], loc='upper right', ncol=3, fontsize='small')

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
ax.set_title('MAIN DECAY CHANNELS\n' + 'GM MODEL + TYPE-II SEESAW (N=20.000)')

# SAVE PLOT
plt.savefig('ctau_H5p_wMDC_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()
