#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Plot cross section
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> January 02, 2023.
"""

###########
# MODULES #
###########

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cls

####################
# PLOT INFORMATION #
####################

energy = '13TeV'
mass_relation = 'm3E2m5'

##################
# IMPORTING DATA #
##################

# Directory to the .csv with cross section data
DirRandomScan = f'/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterDataCrossSection{energy}_H5pp_GMSEESAW.csv'

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

ax.plot(df['m5 [GeV]'], df[f'sigma({mass_relation}_{energy}_ppTOH5ppmmH5mmpp) [pb]'], color='#d62728', linestyle='dashed', label=r'$pp \to H^{++}_5H^{--}_5$')
ax.plot(df['m5 [GeV]'], df[f'sigma({mass_relation}_{energy}_ppTOH5ppmmH5mp) [pb]'], color='green', linestyle='dashed', label=r'$pp \to H^{\pm\pm}_5H^{\mp}_5$')
ax.plot(df['m5 [GeV]'], df[f'sigma({mass_relation}_{energy}_ppTOWpmH3pmTOH5ppmmH3mp) [pb]'], color='blue', linestyle='dashed', label=r'$pp \to H^{\pm\pm}_5H^{\mp}_3$')

######################
# PLOT CUSTOMIZATION #
######################

# AXES LABELS
ax.set_xlabel(r'$m_5$  [GeV]')
ax.set_ylabel(r'$\sigma$  [pb]')

# AXES LIMITS
#ax.set_xlim(left=80, right=712)
#ax.set_ylim(bottom=-250000, top=20000.0)

# LOG AXES
"""You can use 'symlog' instead of 'log' to set the scale to a
symmetric log scale; this allows the use of negative numbers as well"""
# ax.set_xscale('log')
ax.set_yscale('log')

# LEGEND
plt.legend(loc='upper right')

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
plt.title('GM MODEL + TYPE-II SEESAW\n' + r'$m_3 = 2m_5$; $m_H = 200$ [GeV]; $M_1 = 10^{-4}$ [GeV]; $M_2 = 0$ [GeV];' + '\n' r'$v_\Delta = 10^{-4}$ [GeV]; $\sin(\alpha) = 0$; $\sqrt{s} = 13$ [TeV]')

# SAVE PLOT
plt.savefig(f'CrossSection{energy}_{mass_relation}_H5pp_GMSEESAW.png',
        bbox_inches='tight',
        dpi=300)

# SHOW PLOT
plt.show()