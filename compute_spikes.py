# Simulation parameters and library paths
# plot currents with respect to time
from input_data import *

import glob
import sys, os
os.environ['MPLCONFIGDIR'] = base
from scanf import scanf
import spikes
import pandas as pd

#from IPython import embed
import numpy as np
import tools


# compute spike indices for soma and processes for one very long run. 
# store the indices in npy files

# Store data in a pandas datastructure

#====== GENERATE DATA ===========

df_data = []

# changed name not to delete saved files that start with sim_
# Read the bin files and compute averages of calcium traces in process and soma as a function of sigma 
binFiles = glob.glob("diagnostic000_*.bin")
solFiles = glob.glob("sim_000*.bin")


def plotData(sig, dc, der):
    #binFiles = ["diagnostic000_sig%6.4f_dc00.050_der%07.4f.bin" % (sig,der)]
    simFiles = ["sim_run000_sig%6.4f_dc%06.3f_der%07.4f.bin" % (sig,dc,der)]

    #assert(len(binFiles) == 1)
    assert(len(simFiles) == 1)
    
    # Use data from only a single run. Use only a single process
    
    
    for filen in simFiles:
        t, stimuli, soma, randomInput,stim_er = tools.loadDataFast([filen], numProc)
        [run, sig, dc, der] = scanf("sim_run%d_sig%f_dc%f_der%f.bin",  filen)
        times = t

        
        caP = {}
        caP_ER = {}
    
        for i in range(numProc):
           caP[i]    = np.asarray(stimuli["caP%d" % (i+1)])
           caP_ER[i] = np.asarray(stim_er["caP_ER%d" % (i+1)])
    
        caS    = np.asarray(soma["caS"])
        caS_ER = np.asarray(soma["caS_ER"])
        
        process_spike_indices = {}
        for p in range(numProc):
            process_spike_indices[p] = spikes.findSpikes(t, caP[p], thresh=1.2)
            #print "proc %d, spikes: " % i, process_spike_indices[i][0:10]

        soma_spike_indices = spikes.findSpikes(t, caS, thresh=1.2)

        # Compute ISIs
        p_spikes = {}
        p_ISI = {}

        for p in range(numProc):
            p_spikes[p] = process_spike_indices[p]
            p_ISI[p] = times[-1] / len(p_spikes[p])
        s_ISI    = times[-1] / len(soma_spike_indices)

        data = [run, sig, der, dc]
        data.extend([p_ISI[p] for p in range(numProc)])
        data.append(s_ISI)
        df_data.append(data)
    
    columns = ["run", "sig", "der", "dc"]
    columns.extend(["ISI_P%d" % i for i in range(0,numProc)])
    columns.append("ISI_S")

    df = pd.DataFrame(df_data, columns=columns)
    df.to_csv("spike_data.csv", columns=columns)

    spike_indices = {}
    spike_indices['soma_spikes'] = soma_spike_indices

    for p in range(numProc):
        spike_indices['process%d_spikes' % p] = process_spike_indices[p]

    # Same format as simulation files
    np.save(spike_index_file % (run,sig,dc,der), spike_indices)
    np.save("times", t)   # no need to waste file space. Will also speed up loading in later stages
#-----------------------------------

der = [0.00, 0.0002, 0.0004, 0.0006, 0.0008]
der = [0.00, 0.0008, 0.0016, 0.0024, 0.0032, 0.0040]
der = [0.00, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020]
der = der_vals

for s, sig in enumerate(sig_vals):
    for e, de in enumerate(der): 
       for c, dc in enumerate(dc_vals): 
           print "sig= %f, e = %f" % (sig, de)
           plotData(sig, dc, de) 
    
