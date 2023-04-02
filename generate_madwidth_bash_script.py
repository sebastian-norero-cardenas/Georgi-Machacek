#!/usr/bin/env python3

###############
# INFORMATION #
###############

"""
INFO: Model       -> Georgi-Machacek
INFO: Intended    -> Generate a bash script with MadWidth instructions
INFO: Language    -> Python 3
INFO: Author      -> Sebastian Norero C.
INFO: Last update -> October 21, 2022.
"""

###########
# MODULES #
###########

import pandas as pd
import time

start = time.time()

##################
# ABOUT THE RUNS #
##################

max_num_decays = '4'  # Maximum number of particles in a final state
tolerated_error = '.1'  # Tolerated error (.1 = 10%, .01 = 1%, ...) Do NOT include integer part here
new_particles = 'H H3z H3p H5z H5p H5pp'  # Compute decay width of these particles

N_initial = 1  # Initial run
N_final = 3  # Final run
num_runs = N_final - N_initial + 1  # Number of runs

###############
# DIRECTORIES #
###############

DirRandomScan = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW/MasterDataPython3.csv'  # Directory to the .csv with random points
DirMadEvent = f'/home/sebastian/Projects/Thesis/MG5_aMC_v3_4_0/Processes/test/bin/madevent'  # Directory to MadWidth
DirParamCard = f'/home/sebastian/Projects/Thesis/MG5_aMC_v3_4_0/Processes/test/Cards/param_card.dat'  # Location of the param_card.dat that MadWidth uses
DirResults = f'/home/sebastian/Projects/Thesis/python_files_GMSEESAW/param_cards_GMSEESAW{N_initial}TO{N_final}'  # Here we save the resulting param_card.dat's

def CopyParamCard(i): return f'cp {DirParamCard} {DirResults}/param_card_{i}.dat'

##################
# IMPORTING DATA #
##################

df = pd.read_csv(DirRandomScan)

######################
# WRITING THE SCRIPT #
######################

f = open(f'madwidth_bash_script_{N_initial}TO{N_final}_points.sh', 'w')

f.write('#! /bin/bash' + '\n' * 2)

for i in range(N_initial-1, N_final):
    f.write('# --benchmark point = ' + str(i+1) + '\n'
            'python ' + DirMadEvent + ' compute_widths ' + new_particles + ' --body_decay=' + max_num_decays + tolerated_error + ' <<EOD' + '\n'
            '   set tanth ' + str(df['tH'][i]) + '\n'
            '   set m1coeff ' + str(df['M1 [GeV]'][i]) + '\n'
            '   set m2coeff ' + str(df['M2 [GeV]'][i]) + '\n'
            '   set lam2 ' + str(df['lam2'][i]) + '\n'
            '   set lam3 ' + str(df['lam3'][i]) + '\n'
            '   set lam4 ' + str(df['lam4'][i]) + '\n'
            '   set lam5 ' + str(df['lam5'][i]) + '\n'
            '   set mh__2 ' + str(df['mh [GeV]'][i]) + '\n'
            '   set dm21 ' + str(df['Dm21 [eV^2]'][i] * 1E-18) + '\n'
            '   set dm31 ' + str(df['Dm31 [eV^2]'][i] * 1E-18) + '\n'
            '   set pmnss12 ' + str(df['PMNSs12'][i]) + '\n'
            '   set pmnss23 ' + str(df['PMNSs23'][i]) + '\n'
            '   set pmnss13 ' + str(df['PMNSs13'][i]) + '\n'
            '   set pmnsphase ' + str(df['PMNSphase'][i]) + '\n'
            '   set mve ' + str(df['mv1 [eV]'][i] * 1E-9) + '\n'
            '   set ckms12 ' + str(df['CKMs12'][i]) + '\n'
            '   set ckms23 ' + str(df['CKMs23'][i]) + '\n'
            '   set ckms13 ' + str(df['CKMs13'][i]) + '\n'
            '   set ckmphase ' + str(df['CKMphase'][i]) + '\n'
            '   0' + '\n'
            '   quit()' + '\n'
            'EOD' + '\n'
            + CopyParamCard(i+1) + '\n' * 2)

f.close()

##########################
# INFORMATION ABOUT SCAN #
##########################

end = time.time()

runs_in_the_script = "INFO: You will find " + str(num_runs) + " runs in the script."
average_time_per_point = "INFO: Average time to find a compatible point: " + '\n' + str(round((end - start) / num_runs, 4)) + " seconds or " + '\n' + str(round(((end - start) / 60) / num_runs, 4)) + " minutes."
code_execution_time = "INFO: Code execution time: " + '\n' + str(round(end - start, 4)) + " seconds or " + '\n' + str(round((end - start) / 60, 4)) + " minutes or" + '\n' + str(round((end - start) / 3600, 4)) + " hours."

print(code_execution_time)
print(average_time_per_point)
print(runs_in_the_script)