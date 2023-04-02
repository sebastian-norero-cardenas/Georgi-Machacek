#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek + typeII seesaw
INFO: Intended    -> Gather cross sections and create a .csv
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> January 02, 2023.
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
# ABOUT THE RUNS #
##################

m5_max = 1000  # [GeV]
m5_min = 45  # [GeV]

N = 30  # Number of runs

sqrt_s = 14  # [TeV]

def mass5(n):
    return ((m5_max - m5_min) / (N - 1)) * n + m5_min

##################
# IMPORTING DATA #
##################

massenergy_list = [f'm3E05m5_{sqrt_s}TeV', f'm3Em5_{sqrt_s}TeV', f'm3E2m5_{sqrt_s}TeV']

processes_list = ['ppTOH5ppmmH5mmpp', 'ppTOH5ppmmH5mp', 'ppTOWpmH3pmTOH5ppmmH3mp']

# Path to the run_tag_1_banner.txt's
def DirRuns(massenergy, process): return f'/home/sebastian/Projects/Thesis/MG5_aMC_v3_4_0/Processes/{massenergy}/{process}/Events'

###########################
# COLLECTING DECAY WIDTHS #
###########################

logging.info("Collecting run's paths")

# Dictionary of cross sections
cross_section_dict = {'m5 [GeV]': []}

for n in range(N):
    cross_section_dict['m5 [GeV]'].append(mass5(n))


for massenergy in massenergy_list:
    for process in processes_list:
        cross_section_dict[f'sigma({massenergy}_{process}) [pb]'] = []


for massenergy in massenergy_list:
    for process in processes_list:
        # Let's collect the path of every run_tag_1_banner.txt file
        path_list_txt = []

        for root, dirs, files in os.walk(DirRuns(massenergy, process)):
            for file in files:
                if file.endswith('.txt'):
                    path_list_txt.append(os.path.join(root, file))

        path_list_txt = os_sorted(path_list_txt)

        # Let's extract de total decay width (TDW) of every particle
        for path in tqdm(path_list_txt, desc="COLLECTING CROSS SECTIONS"):
            file = open(path, 'r')
            lines = [x.strip() for x in file]
            for line in lines:
                if line.startswith('#  Integrated weight (pb)  :'):
                    cross_section_dict[f'sigma({massenergy}_{process}) [pb]'].append(float(line.split()[5]))


####################
# WRITING THE DATA #
####################

logging.info('Creating DataFrame with all data')
df = pd.DataFrame(data = cross_section_dict)

logging.info('Saving information in .csv file')
df.to_csv(f'MasterDataCrossSection{sqrt_s}TeV_H5pp_GMSEESAW.csv', index=False)

logging.info('DONE!')
