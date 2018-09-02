import sys, os
import subprocess
repo_dir = subprocess.Popen(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).communicate()[0].rstrip()

#base = os.path.join(repo_dir, "comp_astro_evan", "refactored")
#libs = os.path.join(repo_dir, base, "libs")

base = os.path.join(repo_dir)
libs = os.path.join(repo_dir, base, ".")

os.environ['MPLCONFIGDIR'] = base

import numpy as np
sys.path.extend([base,libs])


################################################################
# INPUT DATA
#
#src_file = "astro_2017.c"
# copy of astro-2007.c + output of all terms in equations for balance equations
src_file = "astro_2017_diagnostics.c"
src_path = os.path.join(libs, src_file)

#data_to = data directory with generated data
#data_from = where data is taken for plotting

volumeDiff = 1
time = 500000
time = 1000000
time = 1000000
time = 100000
nSteps = 2*time  # dt = time / nSteps
mu = 0.0
sig = .2
der = 0
corr = .2
addOn = 1.0
numProc = 5
sig_num =  40
der_num = 5

# .001 -> .05
dc = .05

seed = 123456789

nRuns = 1 # must be 1 to use with create_multiple_runs.py
# allows runs to be added to the same directory
min_run = 0
max_run = min_run + nRuns

tau = 8
eps = .04

#generate values for simulations
#der_vals = np.linspace(.0,.001*tau*eps/alpha,der_num)
der_num = 11
der_vals = np.linspace(.0,.004, der_num)
der_num = 6
der_vals = np.linspace(.0,.002, der_num)
der_vals = [0.0]
#print der_vals; quit()

dc_vals = [0.025, 0.050, 0.10]

sig_max = 0.3
sig_vals1 = np.linspace(0, 0.1, 11)[0:-1]   # max(sig) approx 0.3
sig_vals2 = np.linspace(0.1, sig_max, 5)[1:]   # max(sig) approx 0.3
sig_vals3 = np.linspace(0.08, .12, 5)   # max(sig) approx 0.3

sig_vals1 = np.linspace(0, 0.07, 8)   # max(sig) approx 0.3
sig_vals2 = np.linspace(0.08, .12, 9)   # max(sig) approx 0.3
sig_vals3 = np.linspace(0.15, sig_max, 4)   # max(sig) approx 0.3
sig_vals = np.hstack([sig_vals1, sig_vals2, sig_vals3])

sig_vals = [.1, .2, .3]
print sig_vals

# High resolution between 0.08 and 0.12
#sig_vals = np.linspace(0.08, 0.12, 21)  # max(sig) approx 0.3

#sig_vals = [.1,.2,.3] # TESTING
#der_vals = [.0002, .0004, .0006]  # TESTING
#der_num = len(der_vals)  # TESTING

sigNum = len(sig_vals)


output_file      = "sim_run%03d_sig%06.4f_dc%06.3f_der%07.4f.bin" 
diagnostic_file  = "diagnostic%03d_sig%06.4f_dc%06.3f_der%07.4f.bin" 
spike_index_file = "spike_index%03d_sig%06.4f_dc%06.3f_der%07.4f.npy" 

# GE testing
#sig_vals = [0.04, 0.06, 0.1]
#der_num = 1

print "len/sigvals: ", len(sig_vals), sig_vals


save_every = 1

#
# 
################################################################
#quit()

