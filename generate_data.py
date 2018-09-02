# Simulation parameters and library paths
from input_data import *

import glob
import sys, os
os.environ['MPLCONFIGDIR'] = base

#from IPython import embed
import numpy as np
import tools

import calendar
import time as timex  # time is use for other purposes


#====== GENERATE DATA ===========

if 1:
  count = 0

  for run in range(min_run, max_run):
    print "======= run ", run

    for dc_sim,dc in enumerate(dc_vals):
        for der_sim,der in enumerate(der_vals):
            print
            print '  ===== der = ',der
            print

            # changed name not to delete saved files that start with sim_
            binFiles = glob.glob("sim[0-9]*.bin")

            for binFile in binFiles:
                os.remove(binFile)

            cmd = "gcc %(src_path)s -I%(libs)s -o sim_sig -lm -lutils -L%(libs)s -DNPROC=%(numProc)s" % locals()
            returnVal = os.system(cmd)
            if returnVal != 0: sys.exit()
            
            #iterate through sig values
            for x,sig in enumerate(sig_vals):
                print "    ======= sigma ", sig
                count = count + 1
                print '\n\n Running sim  '+str(count)+' of '+ str(len(sig_vals)*len(der_vals)*nRuns) 
                # added an additional 10 when making a new run with different dER
                seedf = seed + 10*run + 100*x + 1000*der_sim + 10
                # Keep seed constant for this experiment (changing Der and examining peak structure as a function of time
                #seedf = seed 

                random_time = calendar.timegm(timex.gmtime())
                seedf = seedf + random_time

                # Gordon's version
                subprocess.call(["./sim_sig",str(time),str(nSteps),str(der),str(sig),str(dc),str(corr),str(addOn),"sim"+str(x)+".bin",str(seedf)])

                #if IC.bin exists remove it; we want individual simulations
                if os.path.exists("IC.bin"):
                    os.remove("IC.bin")

                #load data; we care about caS since we are doing ISI
                filen = 'sim'+str(x)+'.bin'
                file_diagnostic = 'concentration_balance.bin'

                # save full dataset every so often
                if run % save_every == 0: 
                    file_new = output_file % (run, sig, dc, der)
                    cmd = "cp %s %s" % (filen, file_new)
                    os.system(cmd)

                    file_new = diagnostic_file % (run, sig, dc, der)
                    cmd = "cp %s %s" % (file_diagnostic, file_new)
                    os.system(cmd)

