#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw.
INFO: Intended    -> Plot minimum value of ctau with the promptest particle in color.
INFO: Language    -> Python 3.8
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> November 26, 2022.
"""

###########
# MODULES #
###########

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from random import randint

##################
# IMPORTING DATA #
##################

# Absolute path to the MasterData .csv file
DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterData100000_GMSEESAW.csv'

# Absolute path to the PDG_DICTIONARY .csv file
DirPDGDictionary = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW/PDG_DICTIONARY.csv'

# Read files
df = pd.read_csv(DirRandomScan)
df_PDG = pd.read_csv(DirPDGDictionary)

# PDG particles to dictionary
PDG_DICTIONARY = df_PDG.to_dict('records')[0]

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

#################
# COLOR MAPPING #
#################

# USE THIS TO KNOW EVERY DIFFERENT MAIN DECAY CHANNEL
# print(set(df['ctau_minimum_key']))
# print(len(set(df['ctau_minimum_key'])))

# COLOR ASSIGNMENT TO DECAY CHANNELS
colors = {
        'ctau_H5pp [cm]': 'blue',
        'ctau_H5p [cm]': 'orange',
        'ctau_H5z [cm]': 'darkturquoise',
        'ctau_H3p [cm]': 'hotpink',
        'ctau_H3z [cm]': 'brown',
        'ctau_H [cm]': 'green'}

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PROT
sc = ax.scatter(
        x=df['ctau_minimum [cm]'],  # first variable scattered over the x-axis
        y=df['ctau_maximum [cm]'],  # second variable scattered over the y-axis
        c=df['ctau_minimum_key'].map(colors),  # third variable represented as a color gradient
        marker='+',  # type of marker for the scattered points
        s=0.1,  # size of the markers
        alpha=1.0  # transparency of the markers
)

######################
# PLOT CUSTOMIZATION #
######################

# LEGEND
blue_patch = mpatches.Patch(color='blue', label=r'$c\tau(H^{++}_5)$')
orange_patch = mpatches.Patch(color='orange', label=r'$c\tau(H^{+}_5)$')
darkturquoise_patch = mpatches.Patch(color='darkturquoise', label=r'$c\tau(H^{0}_5)$')
hotpink_patch = mpatches.Patch(color='hotpink', label=r'$c\tau(H^{+}_3)$')
brown_patch = mpatches.Patch(color='brown', label=r'$c\tau(H^{0}_3)$')
green_patch = mpatches.Patch(color='green', label=r'$c\tau(H)$')

ax.legend(handles=[
        blue_patch,
        orange_patch,
        darkturquoise_patch,
        hotpink_patch,
        brown_patch,
        green_patch], loc='upper right', ncol=2, fontsize='small')

# AXES LABELS
ax.set_xlabel(r'$\min\{c\tau(X)\}$  [GeV]')
ax.set_ylabel(r'$\max\{c\tau(X)\}$  [cm]')

# AXES LIMITS
#ax.set_xlim(left=80, right=712)
#ax.set_ylim(bottom=-250000, top=20000.0)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
ax.set_xscale('log')
ax.set_yscale('log')

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
ax.set_title('GM MODEL + TYPE-II SEESAW\n' + 'RANDOM SCAN; N=100.000\n' + 'COLOR: PROMPTEST PARTICLE')

# SAVE PLOT
plt.savefig('ctaumin_ctaumax_mincolor_GMSEESAW.png',
        bbox_inches='tight',
        dpi=100)

# SHOW PLOT
plt.show()
