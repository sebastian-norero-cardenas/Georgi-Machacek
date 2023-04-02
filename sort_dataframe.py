#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Generate random scatter plots
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> October 21, 2022.
"""

###########
# MODULES #
###########

import numpy as np
import pandas as pd
import random as  rd
from tqdm import tqdm

##################
# IMPORTING DATA #
##################

# Absolute path to the MasterData .csv file
df = pd.read_csv('/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterData100000_GMSEESAW.csv')

#############
# SORT DATA #
#############

df.sort_values(by='ctau_H5pp [cm]', ascending=False, inplace=True)

# print(df.sort_values(by='m5 [GeV]', ascending=True)['m5 [GeV]'])

#############
# SAVE DATA #
#############

df.to_csv('ctau_H5pp_descending_order.csv', index=False)
