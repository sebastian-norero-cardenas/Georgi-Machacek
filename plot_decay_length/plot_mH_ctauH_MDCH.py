#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw.
INFO: Intended    -> Generate random scatter plot with main decay channel in color.
INFO: Language    -> Python 3.8
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> November 25, 2022.
"""

###########
# MODULES #
###########

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from random import randint

#######################
# EXTERNAL PARAMETERS #
#######################

mh = 125.18
mW = 79.82436
mZ = 91.1876

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

unique_maindecaychannel = list(set(df['H_maindecaychannel'])) 
step_size = 256 ** 3 // len(unique_maindecaychannel)

maindecaychannel_colors = {}

for i, maindecaychannel in enumerate(unique_maindecaychannel):
        maindecaychannel_colors[maindecaychannel] = '#'+ f'{"%06x" % randint(0, 0xFFFFFF)}'

#########
# PLOTS #
#########

fig, ax = plt.subplots(nrows=1, ncols=1)

# SCATTERING PROT
sc = ax.scatter(
        x=df['mH [GeV]'],  # first variable scattered over the x-axis
        y=df['ctau_H [cm]'],  # second variable scattered over the y-axis
        c=df['H_maindecaychannel'].map(maindecaychannel_colors),  # third variable represented as a color gradient
        marker='+',  # type of marker for the scattered points
        s=0.1,  # size of the markers
        alpha=1.0  # transparency of the markers
)

# m5=mh line
plt.axvline(x=1.2518e+02, color="black", linestyle="--", linewidth=0.5)
ax.text(x=1.2518e+02*(1 - 0.15), y=1e-15, s=r'$m_H = m_h$', color="black", fontsize=7, rotation=90)

# m5=2mW line
plt.axvline(x=2*7.982436e+01, color="black", linestyle="--", linewidth=0.5)
ax.text(x=2*7.982436e+01*(1 + 0.05), y=1e-18, s=r'$m_H = 2m_W$', color="black", fontsize=7, rotation=90)

######################
# PLOT CUSTOMIZATION #
######################

handles_list = []

for maindecaychannel, color in maindecaychannel_colors.items():
        decay_label = ''.join([PDG_DICTIONARY[k] for k in eval(maindecaychannel)])
        handles_list.append(mpatches.Patch(color=f'{color}', label=fr'${decay_label}$'))


ax.legend(handles=handles_list, loc='upper right', ncol=3, fontsize='small')

# AXES LABELS
ax.set_xlabel(r'$m_H$  [GeV]')
ax.set_ylabel(r'$c\tau(H)$  [cm]')

# AXES LIMITS
ax.set_xlim(left=50, right=750)
ax.set_ylim(bottom=1e-19, top=1e+1)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
ax.set_yscale('log')

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
ax.set_title('MAIN DECAY CHANNELS\n' + 'GM MODEL + TYPE-II SEESAW\n' + 'RANDOM SCAN; N=100.000')

# SAVE PLOT
plt.savefig('mH_ctauH_MDCH_GMSEESAW.png',
        bbox_inches='tight',
        dpi=100)

# SHOW PLOT
plt.show()
