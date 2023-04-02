#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Gather decay width information to the .csv random scan
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> December 11, 2022.
"""

###########
# MODULES #
###########

import logging
logging.basicConfig(
    #filename = 'file.log',
    format = '%(levelname)s:%(message)s',
    level = logging.DEBUG)

logging.info('Loading libraries')

import os
from natsort import os_sorted
import pandas as pd
from tqdm import tqdm

##################
# IMPORTING DATA #
##################

# Directory to the .csv with random points
DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW_LOGSCAN/MasterDataPython100000_GMSEESAW_LOGSCAN.csv'

# Path to the param_card.dat's
DirParamCards = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW_LOGSCAN/param_cards_GMSEESAW_LOGSCAN'

df = pd.read_csv(DirRandomScan)

# The following is just because we have analyze 80000 points so far
# df.drop(df.index[80000:100000], inplace=True)

###########################
# COLLECTING DECAY WIDTHS #
###########################

logging.info("Collecting param_cards's paths")

# Dictionary of decay widths
decay_width_dict = {'TDW_H [GeV]': [], 'TDW_H3p [GeV]': [], 'TDW_H3z [GeV]': [], 'TDW_H5pp [GeV]': [], 'TDW_H5p [GeV]': [], 'TDW_H5z [GeV]': []}

# Let's collect the path of every param_card.dat file
path_list_txt = []

for root, dirs, files in os.walk(DirParamCards):
    for file in files:
        if file.endswith('.dat'):
            path_list_txt.append(os.path.join(root, file))

path_list_txt = os_sorted(path_list_txt)

# Let's extract de total decay width (TDW) of every particle
for path in tqdm(path_list_txt, desc="COLLECTING TOTAL DECAY WIDTHS"):
    file = open(path, 'r')
    lines = [x.strip() for x in file]
    if any([line.startswith('DECAY  252') for line in lines]):
        pass
    else:
        print(path)
    for line in lines:
        if line.startswith('DECAY  252'):
            decay_width_dict['TDW_H [GeV]'].append(float(line.split()[2]))
        if line.startswith('DECAY  253'):
            decay_width_dict['TDW_H3p [GeV]'].append(float(line.split()[2]))
        if line.startswith('DECAY  254'):
            decay_width_dict['TDW_H3z [GeV]'].append(float(line.split()[2]))
        if line.startswith('DECAY  255'):
            decay_width_dict['TDW_H5pp [GeV]'].append(float(line.split()[2]))
        if line.startswith('DECAY  256'):
            decay_width_dict['TDW_H5p [GeV]'].append(float(line.split()[2]))
        if line.startswith('DECAY  257'):
            decay_width_dict['TDW_H5z [GeV]'].append(float(line.split()[2]))


##################################
# COLLECTING MAIN DECAY CHANNELS #
##################################

maindecaychannel_dict = {
    'H_maindecaychannel': [],
    'H3p_maindecaychannel': [],
    'H3z_maindecaychannel': [],
    'H5pp_maindecaychannel': [],
    'H5p_maindecaychannel': [],
    'H5z_maindecaychannel': []}

# Let's extract de main decay channel of every particle
for path in tqdm(path_list_txt, desc="COLLECTING MAIN DECAY CHANNELS"):
    with open(path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('DECAY  252'):  # H
                data = lines[lines.index(line)+2].split()  # line with main decay channel
                # BR = data[0]
                final_states = []
                for i in range(2, 6):
                    if data[i] != '#':
                        final_states.append(data[i])
                    else:
                        final_states = final_states
                        break
                maindecaychannel_dict['H_maindecaychannel'].append(final_states)
            if line.startswith('DECAY  253'):  # H3p
                data = lines[lines.index(line)+2].split()  # line with main decay channel
                # BR = data[0]
                final_states = []
                for i in range(2, 6):
                    if data[i] != '#':
                        final_states.append(data[i])
                    else:
                        final_states = final_states
                        break
                maindecaychannel_dict['H3p_maindecaychannel'].append(final_states)
            if line.startswith('DECAY  254'):  # H3z
                data = lines[lines.index(line)+2].split()  # line with main decay channel
                # BR = data[0]
                final_states = []
                for i in range(2, 6):
                    if data[i] != '#':
                        final_states.append(data[i])
                    else:
                        final_states = final_states
                        break
                maindecaychannel_dict['H3z_maindecaychannel'].append(final_states)
            if line.startswith('DECAY  255'):  # H5pp
                data = lines[lines.index(line)+2].split()  # line with main decay channel
                # BR = data[0]
                final_states = []
                for i in range(2, 6):
                    if data[i] != '#':
                        final_states.append(data[i])
                    else:
                        final_states = final_states
                        break
                maindecaychannel_dict['H5pp_maindecaychannel'].append(final_states)
            if line.startswith('DECAY  256'):  # H5p
                data = lines[lines.index(line)+2].split()  # line with main decay channel
                # BR = data[0]
                final_states = []
                for i in range(2, 6):
                    if data[i] != '#':
                        final_states.append(data[i])
                    else:
                        final_states = final_states
                        break
                maindecaychannel_dict['H5p_maindecaychannel'].append(final_states)
            if line.startswith('DECAY  257'):  # H5z
                data = lines[lines.index(line)+2].split()  # line with main decay channel
                # BR = data[0]
                final_states = []
                for i in range(2, 6):
                    if data[i] != '#':
                        final_states.append(data[i])
                    else:
                        final_states = final_states
                        break
                maindecaychannel_dict['H5z_maindecaychannel'].append(final_states)

###########################
# COLLECTING DECAY WIDTHS #
###########################

logging.info('Constructing proper lifetimes dictionaries')
# Dictionary of proper lifetimes [seconds]
proper_lifetime_dict = {
    'tau_H [s]': [6.58e-25 * (1/TDW) for TDW in decay_width_dict['TDW_H [GeV]']],
    'tau_H3p [s]': [6.58e-25 * (1/TDW) for TDW in decay_width_dict['TDW_H3p [GeV]']],
    'tau_H3z [s]': [6.58e-25 * (1/TDW) for TDW in decay_width_dict['TDW_H3z [GeV]']],
    'tau_H5pp [s]': [6.58e-25 * (1/TDW) for TDW in decay_width_dict['TDW_H5pp [GeV]']],
    'tau_H5p [s]': [6.58e-25 * (1/TDW) for TDW in decay_width_dict['TDW_H5p [GeV]']],
    'tau_H5z [s]': [6.58e-25 * (1/TDW) for TDW in decay_width_dict['TDW_H5z [GeV]']]}


logging.info('Constructing proper decay length dictionaries')
# Dictionary of proper decay lengths [centimeters]
proper_decay_length_dict = {
    'ctau_H [cm]': [(299792458 * tau) * 1e2 for tau in proper_lifetime_dict['tau_H [s]']],
    'ctau_H3p [cm]': [(299792458 * tau) * 1e2 for tau in proper_lifetime_dict['tau_H3p [s]']],
    'ctau_H3z [cm]': [(299792458 * tau) * 1e2 for tau in proper_lifetime_dict['tau_H3z [s]']],
    'ctau_H5pp [cm]': [(299792458 * tau) * 1e2 for tau in proper_lifetime_dict['tau_H5pp [s]']],
    'ctau_H5p [cm]': [(299792458 * tau) * 1e2 for tau in proper_lifetime_dict['tau_H5p [s]']],
    'ctau_H5z [cm]': [(299792458 * tau) * 1e2 for tau in proper_lifetime_dict['tau_H5z [s]']]}

##########################################
# APPEND DECAY WIDTHS AND MEAN LIFETIMES #
##########################################

logging.info('Adding decay widths to .csv')
for key in decay_width_dict:
    df[key] = decay_width_dict[key]

logging.info('Adding proper mean lifetimes to .csv')
for key in proper_lifetime_dict:
    df[key] = proper_lifetime_dict[key]

logging.info('Adding proper decay lengths to .csv')
for key in proper_decay_length_dict:
    df[key] = proper_decay_length_dict[key]

logging.info('Adding main decay channels to .csv')
for key in maindecaychannel_dict:
    df[key] = maindecaychannel_dict[key]

###########################
# APPEND MIN AND MAX CTAU #
###########################

logging.info('Adding minimum ctau')
df['ctau_minimum [cm]'] = df.loc[:, ['ctau_H [cm]', 'ctau_H3p [cm]', 'ctau_H3z [cm]', 'ctau_H5pp [cm]', 'ctau_H5p [cm]', 'ctau_H5z [cm]']].min(axis=1)

logging.info('Adding minimum ctau (row key)')
df['ctau_minimum_key'] = df.loc[:, ['ctau_H [cm]', 'ctau_H3p [cm]', 'ctau_H3z [cm]', 'ctau_H5pp [cm]', 'ctau_H5p [cm]', 'ctau_H5z [cm]']].idxmin(axis=1)

logging.info('Adding maximum ctau')
df['ctau_maximum [cm]'] = df.loc[:, ['ctau_H [cm]', 'ctau_H3p [cm]', 'ctau_H3z [cm]', 'ctau_H5pp [cm]', 'ctau_H5p [cm]', 'ctau_H5z [cm]']].max(axis=1)

logging.info('Adding maximum ctau (row key)')
df['ctau_maximum_key'] = df.loc[:, ['ctau_H [cm]', 'ctau_H3p [cm]', 'ctau_H3z [cm]', 'ctau_H5pp [cm]', 'ctau_H5p [cm]', 'ctau_H5z [cm]']].idxmax(axis=1)

###############
# SAVING .CSV #
###############

logging.info('Saving .csv')
df.to_csv(f'MasterData{len(df.index)}_GMSEESAW_LOGSCAN.csv', index=False)

logging.info('DONE!')
