import os
from input_data import *

runs = range(1)

for run in runs:
    print "*****  run %d ********" % run
    min_run = run
    max_run = run+1
    os.system("python generate_data.py")

	# creates balance.csv
    os.system("python compute_diagnostics.py")  

	# creates means.csv
    os.system("python compute_means.py")

	# creates spikes.csv
    os.system("python compute_spikes.py")

    run_dir = "run_data%03d" % run
    os.system("mkdir %s" % run_dir)
    os.system("mv balance_data.csv %s" % run_dir)
    os.system("mv mean_data.csv %s" % run_dir)
    os.system("mv spike_data.csv %s" % run_dir)
    os.system("mv *.npy %s" % run_dir)
    os.system("rm *.bin")

