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
# print(set(df['H5z_maindecaychannel']))
# print(len(set(df['H5z_maindecaychannel'])))

# COLOR ASSIGNMENT TO DECAY CHANNELS
colors = {
        "['16', '16']": 'orange',
        "['23', '23']": 'red',
        "['254', '23']": 'yellowgreen',
        "['-24', '24']": 'darkturquoise',
        "['253', '-24']": 'pink',
        "['254', '254']": 'brown',
        "['-24', '-1', '2']": 'peru',
        "['-24', '-3', '4']": 'sandybrown',
        "['-2', '1', '24']": 'olive',
        "['-4', '3', '24']": 'blue',
        "['12', '12', '254']": 'black',
        "['14', '14', '254']": 'teal',
        "['16', '16', '254']": 'green',
        "['-3', '-1', '1', '3']": 'khaki',
        "['-4', '-3', '3', '4']": 'purple',
        "['-2', '2', '12', '12']": 'palegreen',
        "['-4', '4', '12', '12']": 'steelblue',
        "['-1', '1', '14', '14']": 'slateblue',
        "['-2', '2', '14', '14']": 'hotpink',
        "['-4', '4', '14', '14']": 'mediumorchid',
        "['-1', '1', '16', '16']": 'bisque',
        "['-2', '2', '16', '16']": 'rosybrown'}

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['m5 [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H5z [cm]'],  # second variable scattered over the y-axis
        c=df['H5z_maindecaychannel'].map(colors),  # third variable represented as a color gradient
        marker='.',  # type of marker for the scattered points
        s=1,  # size of the markers
        alpha=1.0  # transparency of the markers
)

# m5=mW line
plt.axvline(x=7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=7.982436e+01*(1 - 0.22), y=1e+00, s=r'$m_5 = m_W$', color="black", fontsize=7, rotation=90)

# m5=2mW line
plt.axvline(x=2*7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*7.982436e+01*(1 - 0.12), y=1e-02, s=r'$m_5 = 2m_W$', color="black", fontsize=7, rotation=90)

# m5=2mZ line
plt.axvline(x=2*91.1876, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*9.11876e+01*(1 - 0.1), y=1e-05, s=r'$m_5 = 2m_Z$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# LEGEND
orange_patch = mpatches.Patch(color='orange', label=r'$\nu_\tau\nu_\tau$')
red_patch = mpatches.Patch(color='red', label=r'$ZZ$')
yellowgreen_patch = mpatches.Patch(color='yellowgreen', label=r'$H^{0}_3Z$')
darkturquoise_patch = mpatches.Patch(color='darkturquoise', label=r'$W^{-}W^{+}$')
pink_patch = mpatches.Patch(color='pink', label=r'$H^{+}_3W^{-}$')
brown_patch = mpatches.Patch(color='brown', label=r'$H^{0}_3H^{0}_3$')
peru_patch = mpatches.Patch(color='peru', label=r'$\bar{d}uW^{-}$')
sandybrown_patch = mpatches.Patch(color='sandybrown', label=r'$\bar{s}cW^{-}$')
olive_patch = mpatches.Patch(color='olive', label=r'$\bar{u}dW^{+}$')
blue_patch = mpatches.Patch(color='blue', label=r'$\bar{c}sW^{+}$')
black_patch = mpatches.Patch(color='black', label=r'$\nu_e\nu_eH^{0}_3$')
teal_patch = mpatches.Patch(color='teal', label=r'$\nu_\mu\nu_\mu H^{0}_3$')
green_patch = mpatches.Patch(color='green', label=r'$\nu_\tau\nu_\tau H^{0}_3$')
khaki_patch = mpatches.Patch(color='khaki', label=r'$\bar{s}s\bar{d}d$')
purple_patch = mpatches.Patch(color='purple', label=r'$\bar{c}\bar{s}sc$')
palegreen_patch = mpatches.Patch(color='palegreen', label=r'$\bar{u}u\nu_e\nu_e$')
steelblue_patch = mpatches.Patch(color='steelblue', label=r'$\bar{c}c\nu_e\nu_e$')
slateblue_patch = mpatches.Patch(color='slateblue', label=r'$\bar{d}d\nu_\mu\nu_\mu$')
hotpink_patch = mpatches.Patch(color='hotpink', label=r'$\bar{u}u\nu_\mu\nu_\mu$')
mediumorchid_patch = mpatches.Patch(color='mediumorchid', label=r'$\bar{c}c\nu_\mu\nu_\mu$')
bisque_patch = mpatches.Patch(color='bisque', label=r'$\bar{d}d\nu_\tau\nu_\tau$')
rosybrown_patch = mpatches.Patch(color='rosybrown', label=r'$\bar{u}u\nu_\tau\nu_\tau$')

ax.legend(handles=[
        orange_patch,
        red_patch,
        yellowgreen_patch,
        darkturquoise_patch,
        pink_patch,
        brown_patch,
        peru_patch,
        sandybrown_patch,
        olive_patch,
        blue_patch,
        black_patch,
        teal_patch,
        green_patch,
        khaki_patch,
        purple_patch,
        palegreen_patch,
        steelblue_patch,
        slateblue_patch,
        hotpink_patch,
        mediumorchid_patch,
        bisque_patch,
        rosybrown_patch], loc='upper right', ncol=3, fontsize='small')

# AXES LABELS
ax.set_xlabel(r'$m_5$  [GeV]')
ax.set_ylabel(r'$c\tau(H^{0}_5)$  [cm]')

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
plt.savefig('ctau_H5z_wMDC_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()
