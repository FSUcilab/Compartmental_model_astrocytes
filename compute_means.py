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


#def newRecord():


#====== GENERATE DATA ===========

# changed name not to delete saved files that start with sim_
# Read the bin files and compute averages of calcium traces in process and soma as a function of sigma 
binFiles = glob.glob("sim_*.bin")

#Store in arrays panda list
# Record: run, sig, der, dc, mean(P1->P5), mean (P1_ER->P5_ER), mean(S), mean(S_ER)
# 15 columns

columns = ["run", "sig", "der", "dc"]
columns.extend(["mean_P%d" % i for i in range(1,6)])
columns.extend(["mean_P%dER" %i for i in range(1,6)])
columns.extend(["mean_S", "mean_SER"])
#df = pd.DataFrame(index=np.arange(0,5), columns=columns)
#df = pd.DataFrame(columns=columns)
print "columns: ", columns
#print df

df_data = []

for filen in binFiles:
    t, stimuli, soma, randomInput,stim_er = tools.loadDataFast([filen], numProc)

    """
    print type(t)
    print stimuli.keys()
    print soma.keys()
    print randomInput.keys()
    print stim_er.keys()
    """

    """
    <type 'list'>   # t
    ['caP4', 'caP5', 'caP1', 'caP2', 'caP3']  # stimuli
    ['caS', 'caS_ER']  # soma
    ['dr1', 'dr3', 'dr2', 'dr5', 'dr4']  # randomInput
    ['caP_ER1', 'caP_ER3', 'caP_ER2', 'caP_ER5', 'caP_ER4']  # stim_er
    """

    caS    = np.asarray(soma["caS"])
    caS_ER = np.asarray(soma["caS_ER"])
    
    caP = {}
    caP_ER = {}
    caP_mean = {}
    caP_ER_mean = {}

    for i in range(numProc):
       caP[i]    = np.asarray(stimuli["caP%d" % (i+1)])
       caP_ER[i] = np.asarray(stim_er["caP_ER%d" % (i+1)])
       caP_mean[i] = np.mean(caP[i])
       caP_ER_mean[i] = np.mean(caP_ER[i])
    caS_mean = np.mean(caS)
    caS_ER_mean = np.mean(caS_ER)

    print "%s: <caS>= %f, <caS_ER>= %f" % (filen, caS_mean, caS_ER_mean)
    print "%s: <caP1>= %f, <caP1_ER>= %f" % (filen, caP_mean[0], caP_ER_mean[0])

    # Extract parameters from file name
    # output_file = "sim_run%03d_sig%06.4f_dc%06.3f_der%07.4f.bin" 
    # sim_run000_sig0.0000_dc00.050_der00.0000.bin
    #[run, sig, dc, der] = scanf("sim_run%03d_sig%06.4f_dc%06.3f_der%07.4f.bin" % filen)

    #print "===> ", filen
    [run, sig, dc, der] = scanf("sim_run%d_sig%f_dc%f_der%f.bin",  filen)

    df_data.append([run, sig, der, dc, caP_mean[0], caP_mean[1], caP_mean[2], caP_mean[3], caP_mean[4], 
                               caP_ER_mean[0], caP_ER_mean[1], caP_ER_mean[2], caP_ER_mean[3], caP_ER_mean[4], 
                                caS_mean, caS_ER_mean])

df = pd.DataFrame(df_data, columns=columns)
df.to_csv("mean_data.csv", columns=columns)


