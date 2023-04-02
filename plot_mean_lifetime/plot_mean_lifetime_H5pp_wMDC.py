#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Generate random scatter plot with main decay channel in color
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

# PLOT STYLE
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
# print(set(df['H5pp_maindecaychannel']))
# print(len(set(df['H5pp_maindecaychannel'])))

# COLOR ASSIGNMENT TO DECAY CHANNELS
colors = {
        "['-13', '-13']": 'coral',
        "['-13', '-15']": 'deepskyblue',
        "['-15', '-15']": 'black',
        "['24', '24']": 'blue',
        "['253', '24']": 'deeppink',
        "['253', '253']": 'limegreen',
        "['-1', '2', '24']": 'olive',
        "['-3', '4', '24']": 'red',
        "['-1', '2', '253']": 'pink',
        "['-3', '4', '253']": 'purple',
        "['-3', '-1', '2', '4']": 'orange',
        "['-5', '-3', '4', '4']": 'darkturquoise',
        "['-15', '-1', '2', '16']": 'peru',
        "['-15', '-3', '4', '16']": 'snow'}

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PLOT
sc = ax.scatter(
        x=df['m5 [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H5pp [cm]'],  # second variable scattered over the y-axis
        c=df['H5pp_maindecaychannel'].map(colors),  # third variable represented as a color gradient
        marker='.',  # type of marker for the scattered points
        s=1,  # size of the markers
        alpha=1.0  # transparency of the markers
)

# m5=mW line
plt.axvline(x=7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=7.982436e+01*(1 - 0.30), y=1e+00, s=r'$m_5 = m_W$', color="black", fontsize=7, rotation=90)

# m5=2mW line
plt.axvline(x=2*7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*7.982436e+01*(1 - 0.15), y=1e-02, s=r'$m_5 = 2m_W$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

# LEGEND
coral_patch = mpatches.Patch(color='coral', label=r'$\mu^{+}\mu^{+}$')
deepskyblue_patch = mpatches.Patch(color='deepskyblue', label=r'$\mu^{+}\tau^{+}$')
black_patch = mpatches.Patch(color='black', label=r'$\tau^{+}\tau^{+}$')
blue_patch = mpatches.Patch(color='blue', label=r'$W^{+}W^{+}$')
deeppink_patch = mpatches.Patch(color='deeppink', label=r'$H^{+}_3W^{+}$')
limegreen_patch = mpatches.Patch(color='limegreen', label=r'$H^{+}_3H^{+}_3$')
olive_patch = mpatches.Patch(color='olive', label=r'$\bar{d}uW^{+}$')
red_patch = mpatches.Patch(color='red', label=r'$\bar{s}cW^{+}$')
pink_patch = mpatches.Patch(color='pink', label=r'$\bar{d}uH^{+}_3$')
purple_patch = mpatches.Patch(color='purple', label=r'$\bar{s}cH^{+}_3$')
orange_patch = mpatches.Patch(color='orange', label=r'$\bar{s}c\bar{d}u$')
darkturquoise_patch = mpatches.Patch(color='darkturquoise', label=r'$\bar{b}c\bar{s}c$')
peru_patch = mpatches.Patch(color='peru', label=r'$\tau^{+}\nu_\tau\bar{d}u$')
snow_patch = mpatches.Patch(color='snow', label=r'$\tau^{+}\nu_\tau\bar{s}c$')

ax.legend(handles=[
        coral_patch,
        deepskyblue_patch,
        black_patch,
        blue_patch,
        deeppink_patch,
        limegreen_patch,
        olive_patch,
        red_patch,
        pink_patch,
        purple_patch,
        orange_patch,
        darkturquoise_patch,
        peru_patch,
        snow_patch], loc='upper right', ncol=3, fontsize='small')

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
plt.savefig('ctau_H5pp_wMDC_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()
