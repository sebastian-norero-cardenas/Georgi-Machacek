###########
# MODULES #
###########

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

fig, ax1 = plt.subplots(nrows=1, ncols=1)

# HISTOGRAM
ax1.hist(
        df['vDelta [GeV]'],
        bins=30,  # Use 'auto' for an automatic selection of bins
        #density=True,
        rwidth=0.85,
        color='#0504aa',
        alpha=0.7
        #label='Data'
)


ax2 = ax1.twinx()

ax2.hist(
        df['vDelta [GeV]'],
        bins=30,  # Use 'auto' for an automatic selection of bins
        # cumulative=True,
        density=True,
        rwidth=0.85,
        color='red',
        alpha=0.1
        #label='Data'
)

######################
# PLOT CUSTOMIZATION #
######################

# AXES LABELS
ax1.set_xlabel(r'$v_\Delta$  [GeV]')

ax1.set_ylabel(r'Frequency')
ax2.set_ylabel(r'Frequency')

# LEGEND
#plt.legend(loc='upper left')

# AXES LIMITS
# ax.set_xlim(left=0, right=712)
# ax.set_ylim(bottom=0, top=700)

# LOG AXES
# ax1.set_xscale('log')
ax1.set_yscale('log')

# ax2.set_xscale('log')
# ax2.set_yscale('log')

# USE A GRAY BACKGROUND
# ax.set_facecolor('#E6E6E6')
# ax.set_axisbelow(True)

########################
# SAVE AND PRINT PLOTS #
########################

# TITLE OF PLOT
plt.title('HISTOGRAM\n' + 'GM MODEL + TYPE-II SEESAW\n' + 'RANDOM SCAN; N=100.000')

# SAVE PLOT
plt.savefig('TEST.png',
            bbox_inches='tight',
            dpi=100)

# SHOW PLOT
plt.show()
