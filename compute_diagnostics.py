# Simulation parameters and library paths
from input_data import *

import glob
import sys, os
import pandas as pd
os.environ['MPLCONFIGDIR'] = base
from scanf import scanf

#from IPython import embed
import numpy as np
import tools


# compute all terms in all equations, average over processes. 
# Store in an excel spreadsheet


#====== GENERATE DATA ===========

# changed name not to delete saved files that start with sim_
# Read the bin files and compute averages of calcium traces in process and soma as a function of sigma 
binFiles = glob.glob("diagnostic*.bin")

#Store in arrays panda list
# Record: run, sig, der, dc, mean(P1->P5), mean (P1_ER->P5_ER), mean(S), mean(S_ER)
# 15 columns

columns = ["run", "sig", "der", "dc"]
columns.extend(["mean_serca", "mean_cicr", "mean_linear", "mean_ca_removal",  "mean_diffc", "mean_differ"]) 
columns.extend(["mean_serca_soma", "mean_cicr_soma", "mean_linear_soma", "mean_ca_removal_soma", 
    "mean_diffc_soma", "mean_differ_soma"])

df_data = []

for filen in binFiles:
    t, serca_current, cicr_current, linear_current, ca_removal_c, diffc_c, differ_c = tools.loadDiagnosticsFast([filen], numProc)

    serca = {}
    cicr = {}
    linear = {}
    ca_removal = {}
    differ = {}
    diffc = {}

    for i in range(numProc+1):
       cicr[i]       = np.asarray(cicr_current["%d" % i])
       serca[i]      = np.asarray(serca_current["%d" % i])
       linear[i]     = np.asarray(linear_current["%d" % i])
       ca_removal[i] = np.asarray(ca_removal_c["%d" % i])
       differ[i]     = np.asarray(differ_c["%d" % i])
       diffc[i]      = np.asarray(diffc_c["%d" % i])

       cicr[i]       = np.mean(cicr[i]);
       serca[i]      = np.mean(serca[i]);
       linear[i]     = np.mean(linear[i]);
       ca_removal[i] = np.mean(ca_removal[i]);
       differ[i]     = np.mean(differ[i]);
       diffc[i]      = np.mean(diffc[i]);

    mean_cicr        = np.mean([cicr[0],   cicr[1],   cicr[2],   cicr[3],   cicr[4]])
    mean_serca       = np.mean([serca[0],  serca[1],  serca[2],  serca[3],  serca[4]])
    mean_linear      = np.mean([linear[0], linear[1], linear[2], linear[3], linear[4]])
    mean_ca_removal  = np.mean([ca_removal[0], ca_removal[1], ca_removal[2], ca_removal[3], ca_removal[4]])
    mean_differ      = np.mean([differ[0], differ[1], differ[2], differ[3], differ[4]])
    mean_diffc       = np.mean([diffc[0],  diffc[1],  diffc[2],  diffc[3],  diffc[4]])

    mean_cicr_soma       = cicr[numProc]
    mean_serca_soma      = serca[numProc]
    mean_linear_soma     = linear[numProc]
    mean_ca_removal_soma = ca_removal[numProc]
    mean_differ_soma     = differ[numProc]
    mean_diffc_soma      = diffc[numProc]

    [run, sig, dc, der] = scanf("diagnostic%d_sig%f_dc%f_der%f.bin",  filen)

    df_data.append([run, sig, der, dc, mean_serca, mean_cicr, mean_linear, mean_ca_removal,  mean_diffc, mean_differ,
         mean_serca_soma, mean_cicr_soma, mean_linear_soma, mean_ca_removal_soma, mean_diffc_soma, mean_differ_soma])

df = pd.DataFrame(df_data, columns=columns)
df.to_csv("balance_data.csv", columns=columns)


