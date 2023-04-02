import pandas as pd

DirPythonRandomScan1 = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW_LOG/MasterDataPython20000_GMSEESAW_LOGSCAN.csv'
DirPythonRandomScan2 = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW_LOG/MasterDataPython20000_GMSEESAW_LOGSCAN_2to4.csv'
DirPythonRandomScan3 = '/home/sebastian/Projects/Thesis/python_files_GMSEESAW_LOG/MasterDataPython60000_GMSEESAW_LOGSCAN_4to10.csv'

df1 = pd.read_csv(DirPythonRandomScan1)
df2 = pd.read_csv(DirPythonRandomScan2)
df3 = pd.read_csv(DirPythonRandomScan3)

df = pd.concat([df1, df2, df3])

nrows = df.shape[0]

df.to_csv(f'MasterDataPython{nrows}.csv', index=False)
