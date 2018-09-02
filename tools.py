import numpy as np
from pylab import setp
import sys, os
import STA
import numpy
import copy  as cpy
import struct
os.environ['MPLCONFIGDIR'] = '../'
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython import embed
from matplotlib.mlab import griddata
import matplotlib.colors as colors

import time as tim

from mpl_toolkits.mplot3d import Axes3D


def findSpikes(signal, thresh=1.2):

    # np.diff(a) : a[n+1] - a[n]
    return np.where(np.diff(np.sign(signal-thresh)) > 0)[0]


def globalSpikes(processes, soma, window=20, numProc=4):
    # ASSUMING THAT THERE ARE 5 PROCESSES TOTAL
    globalEvents = []
    somaSpikes = findSpikes(soma)
    #noSpikeFlag = False

    for somaSpike in somaSpikes:
        noSpikeCount = 0
        for process in processes:
            procSpikes = findSpikes(process[somaSpike-window:somaSpike+window])
            if procSpikes.size == 0:
                #noSpikeFlag = True
                noSpikeCount = noSpikeCount+1
                #break
        #if not noSpikeFlag:
        #globalEvents.append(somaSpike)

        if noSpikeCount <= 5 - numProc:
             globalEvents.append(somaSpike)

        #noSpikeFlag = False

    return globalEvents

def plotGlobalSpikes(processes, soma, window, numProc):
    # ASSUMING THAT THERE ARE 5 PROCESSES TOTAL
    gevents = globalSpikes(processes, soma, window,numProc)
    fig = plt.figure()
    numprocs = len(processes)

    for axnum, process in enumerate(processes):
        ax = fig.add_subplot(numprocs+1, 1, axnum+1)
        ax.plot(process)
        ax.set_xticks([])
        ax.set_yticks([0,.5,1.0,1.5])
        ax.legend()
        if(axnum==3):
                ax.set_ylabel('$Ca^{2+}$')
        for ge in gevents:
            #ax.plot([ge-window, ge-window], [0, 1.5], color='black')
            #ax.plot([ge+window, ge+window], [0, 1.5], color='black')
            ax.axvspan(ge-window,ge+window,color='red',alpha=.5)

    somaaxis = fig.add_subplot(numprocs+1, 1, numprocs+1)
    somaaxis.plot(soma)
    for ge in gevents:
        somaaxis.axvspan(ge-window,ge+window,color='red',alpha=.5)

    plt.xlabel("$T$")
    plt.savefig("globalSpikes.pdf")


def loadData(files = [],numproc = 5):
        """ loading data from .bin files
        args:
                files: (list)

        output:
                (t, stimuli, soma, randomInput) : list

        """
        print("In loadData: numproc = ", numproc, "\n")

        assert len(files) > 0, "ERROR: Driver requires data files!"



        struct_fmt = "=%dd" %(3*numproc+3)
        struct_len = struct.calcsize(struct_fmt)
        struct_unpack = struct.Struct(struct_fmt).unpack_from
        stimuli = {}
        stim_er = {}
        randomInput = {}
        for proc in np.arange(numproc):
                stimuli["caP%d" % (proc+1)] = []
                stim_er["caP_ER%d" % (proc+1)] = []
                randomInput["dr%d" % (proc+1)] = []
        soma = {}
        soma["caS"] = []
        soma["caS_ER"] = []
        t = []
        for x in files:
            #print len(t)
            #print x
            with open(x, "rb") as f:
                while True:
                    data = f.read(struct_len)
                    if not data: break
                    s = struct_unpack(data)
                    #print s
                    t.append(s[0])
                    for proc in np.arange(numproc):
                        stimuli["caP%d" % (proc+1)].append(s[2*proc+1])
                        stim_er["caP_ER%d" % (proc+1)].append(s[2*proc+2])
                    soma["caS"].append(s[2*numproc+1])
                    soma["caS_ER"].append(s[2*numproc+2])
                    for proc in np.arange(numproc):
                        randomInput["dr%d" % (proc+1)].append(s[2*numproc+3+proc])
                        if s[2*numproc+3+proc]>1e5:
                                print("oh shit")
                                sys.exit()

        return (t, stimuli, soma, randomInput,stim_er)

#----------------------------------------------------------------------
# GE, 3/1/18: 
def loadDataFast(files = [],numproc = 5):
        """ loading data from .bin files
        args:
                files: (list)

        output:
                (t, stimuli, soma, randomInput) : list
        """

        print "file: ", files[0]
        f = open(files[0], "rb")
        f.seek(0,2)
        nb_floats = f.tell() / 8  # 8 for double
        nb_vars = 3*numproc + 3
        
        t1 = tim.time()
        rdata = np.fromfile(files[0], dtype=float, count=nb_floats)
        # Order of variables 
        rdata = rdata.reshape(rdata.size / nb_vars, nb_vars).T

        stimuli = {}
        stim_er = {}
        randomInput = {}
        soma = {}
        
        for proc in range(numproc):
            stimuli["caP%d" % (proc+1)]    = rdata[2*proc+1]
            stim_er["caP_ER%d" % (proc+1)] = rdata[2*proc+2]
            randomInput["dr%d" % (proc+1)] = rdata[2*numproc+3+proc]
        soma["caS"]    = rdata[2*numproc+1]
        soma["caS_ER"] = rdata[2*numproc+2]
        #t2 = tim.time(); print "np.fromfile: ", t2-t1
        t = rdata[0]

        return (t, stimuli, soma, randomInput,stim_er)

#----------------------------------------------------------------------
# GE, 3/4/18: 
def loadDiagnosticsFast(files = [],numproc = 5):
        """ loading data from .bin files
        args:
                files: (list)

        output:
                (t, stimuli, soma, randomInput) : list
        """

        f = open(files[0], "rb")
        f.seek(0,2)
        nb_floats = f.tell() / 8
        nb_vars = 1 + 6*(numproc+1)
        # 0th var: time
        
        t1 = tim.time()
        rdata = np.fromfile(files[0], dtype=float, count=nb_floats)
        #print rdata.shape
        # Order of variables 
        rdata = rdata.reshape(rdata.size / nb_vars, nb_vars).T
        
        serca_current = {}
        cicr_current = {}
        linear_current = {}
        ca_removal = {}
        diffc = {}
        differ = {}


        t = rdata[0]

        for proc in range(numproc+1):
            serca_current["%d" % proc]   = rdata[proc+1]
            cicr_current["%d" % proc]    = rdata[(numproc+1)+proc+1]
            linear_current["%d" % proc]  = rdata[2*(numproc+1)+proc+1]
            ca_removal["%d" % proc]      = rdata[3*(numproc+1)+proc+1]
            diffc["%d" % proc]           = rdata[4*(numproc+1)+proc+1]
            differ["%d" % proc]          = rdata[5*(numproc+1)+proc+1]

        #t2 = tim.time(); print "np.fromfile: ", t2-t1
        t = rdata[0]

        return (t, serca_current, cicr_current, linear_current, ca_removal, diffc, differ)

#----------------------------------------------------------------------
def plotTrace(t , traces ,caS):
        """ Creating sta of subject with respect to stumli

        args:
                numprocs: (int)
                        the number of processes in the simulation

                traces: (ndarray of arrays)
                        calcium traces of each process. The number
                        of elements in stimuli is the number of
                        processes to calculate the sta for. Last
                        trace needs to be the soma

        """
        fig1 = plt.figure(figsize=(12,8))

        ax1 = fig1.add_subplot(111)
        ax1.set_xlabel('$t$ $(ms)$')
        ax1.set_ylabel('$Ca^{2+}$')
        set_pub(ax1)

        set_pub(ax1)
        for i in range(len(traces)):
                ax1.plot(t,traces[i],  label='Process '+str(i+1),linewidth=1)


        ax1.plot(t,caS, label='Soma',linewidth=2,marker='*',markevery=100)
        ax1.set_ylim([0,2])
        ax1.set_xlim([190000,191000])
        plt.legend()

        plt.savefig('doubleTap.pdf')
        plt.close()
        sys.exit()
        return

def createSTACurves(stimuli, subject,names):
        """ Creating sta of subject with respect to stumli

        args:
                stimuli: (ndarray of arrays)
                        calcium traces of each process. The number
                        of elements in stimuli is the number of
                        processes to calculate the sta for.

                subject: ndarray
                        The subject of the sta. We're calculating the
                        sta between each process in stimuli and this
                        subject. We use its spike as event around
                        which the averages are generated.
                names: list
                        list of strings which have names

        output:
                staCurve: list
                        This is a list of sta instances describing the
                        sta between each stimuli and the subject

        """

        staCurve = []
        for x,stim in enumerate(stimuli):
                #embed()
                stimulus = np.asarray(stim)
                subject = np.asarray(subject)
                staCurve.append(STA.STA(subject,stimulus,name = names[x]))
                staCurve[-1].LinearNonlinearPoisson()
                print('Found '+str(staCurve[-1].spikeCount)+' spikes')
                staCurve[-1].STAR()

        return staCurve

        #old sta procedure for reference 3/28
        """print "sta for process 5"
        stim = caP5
        stim = np.asarray(stim)
        subject = np.asarray(caS)

        sta5 = STA(subject,stim)
        sta5.LinearNonlinearPoisson()
        #sta5.STAR()
        procSpikes.append(sta5.preSpikeStim)
        print sta5.sta
        print len(sta5.sta)
        """



def plotSTACurves(staInstances,names = []):
        """ Plotting sta of subject with respect to stumli

        args:
                STAs: (ndarray of STA intances)
                        the STA curves to be compared.. note
                        they should have the same subject to
                        be validly compared

        """
        #initiate plotting for sta figures
        fig2 = plt.figure()

        ax2 = fig2.add_subplot(111)

        ax2.set_xlabel('t (s)')
        ax2.set_ylabel('Ca2+')

        set_pub(ax2)

        ax2.set_ylim([0,1.5])
        for staI in staInstances:
                ax2.plot(np.linspace(-staI.stim_size/1000,0,staI.stim_size),staI.sta,label=staI.name)
        ax2.plot(np.linspace(-staI.stim_size/1000,0,staInstances[0].stim_size),staInstances[0].staR,label='Random')
        leg = ax2.legend(fontsize=25,loc='upper left')
        plt.show()


def ISI(trace, thr=1.2, t_ref=104, inputs=None, time_full=0):
        """Calculating the ISI of a given trace and giving the
           average pre-spike spike count of given inputs
        args:
                trace: (ndarray)
                        trace of subject for ISI calculation
                inputs: (ndarray of ndarrays)(optional)
                        if given, inputs are used to find the number of spikes
                        leading to a subject trace spike
                t_ref: refractory time. 
        output:
                ISI: (float)
                        Inter Spike Interval of trace
                avg_pre_spikes: (float)
                        average number of spike leading to a subject spike
        """


        #find the trace ISI
        spikeCount = 0
        pre_spikes = []
        t_r = 0

        print("\ntime_full = ", time_full, "\n")
        for index, signal in enumerate(trace):
                t_r+=1
                # consider the calcium concentration of note and see if it is above threshold (spiking)
                if signal >=thr and t_r>t_ref :
                        spikeCount += 1
                        t_r = 0

                        #if inputs to subject trace are provided, count the number of spikes
                        #avg_prespikes no longer used. Should be debugged if used in the future
                        old_sig = 10
                        if (inputs):
                                count = 0
                                for inSignal in inputs:
                                        # GE:  What is 2000?
                                        for stim_sig in inSignal[index-(2000-1):index]:
                                                if stim_sig >= thr and old_sig < thr:
                                                        count += 1
                                                old_sig = cpy.copy(stim_sig)
                                pre_spikes.append(count)

#                        print 'spike count is ', spikeCount
        if(spikeCount==0):
            #isi = 5500
            isi = time_full + 1 # GE on 2/28/18
            average = 0
        else:
            isi = time_full/spikeCount
            average = np.mean(pre_spikes)


        return isi, average


"""
if 0 :
        #initiate plotting for calcium traces
        fig1 = plt.figure()

        ax1 = fig1.add_subplot(111)

        ax1.set_title("Astro PLot")
        ax1.set_xlabel('t')
        ax1.set_ylabel('Ca2+')

        ax1.plot(t,caP1, c='r', label='process 1',linewidth=.5)
        ax1.plot(t,caP2, c='y', label='process 2',linewidth=.5)
        ax1.plot(t,caP3, c='b', label='process 3',linewidth=.5)
        ax1.plot(t,caP4, c='m', label='process 4',linewidth=.5)
        ax1.plot(t,caP5, c='c', label='process 5',linewidth=.5)
        ax1.plot(t,actual_input, c='k', label='actual input',linewidth=2)
        ax1.plot(t,soma_output, label='soma out',linewidth=2)
        ax1.plot(t,caS,c='g', label='soma',linewidth=2)
        ax1.plot(t,caS_ER, label='soma ER',linewidth=2)

        leg = ax1.legend()




        #initialze vector to store process spikes counts preceding soma spikes
        #by 2 seconds
        procSpikes = []

"""


if 0:


        meanSTA = []
        for x,sig in enumerate(sta2.sta):
                meanSTA.append(np.mean([sta1.sta[x],sta2.sta[x],sta3.sta[x],sta4.sta[x],sta5.sta[x]]))


        print('This was done over..')

        print(len(sta1.spikeTimes))
        print(len(sta2.spikeTimes))
        print(len(sta3.spikeTimes))
        print(len(sta4.spikeTimes))
        print(len(sta5.spikeTimes))

        print('... spikes')

        print('Spike Counts in the processes in 2 s window preceding spikes')

        print(procSpikes)

        #"""

                #"""

        #-----------------------------------------------#
        #-----------------------------------------------#
        #----------------STA for Soma-------------------#
        print()
        print("sta for process soma")
        print()
        stim = caS
        stim = np.asarray(caS)
        subject = np.asarray(caS)

        staS = STA(subject,stim)
        staS.LinearNonlinearPoisson()
        #staS.STAR() # deleted 3/21 see line 254
        procSpikes.append(staS.preSpikeStim)
        print(staS.sta)
        print(len(staS.sta))

        #-----------------------------------------------#

        #-----------------------------------------------#
        #---------STA for Total Input----------------#

        actual_input = [i-.5 for i in actual_input ]
        print("sta for process for actual input")
        stim = actual_input
        subject = np.asarray(caS)

        staInput = STA(subject,stim)
        staInput.LinearNonlinearPoisson()

        #-----------------------------------------------#

        #-----------------------------------------------#
        #---------------------STA for ER----------------#

        print("sta for ER with respect to soma")
        stim = np.asarray(caS_ER)
        subject = np.asarray(caS)

        staER = STA(subject,stim)
        staER.LinearNonlinearPoisson()


        #-----------------------------------------------#

        #-----------------------------------------------#
        #----------------STA for Soma Output---------------------#

        print("sta for output from  soma")

        stim = np.asarray(soma_output)

        subject = np.asarray(caS)

        staOut = STA(subject,stim)
        staOut.LinearNonlinearPoisson()


        #-----------------------------------------------#


        #print staInput.sta
        #print len(staInput.sta)



        run_i = sys.argv[1]
        #comile a list of relative process spike times and bin and save
        rel_spikeTimes = sta1.rel_procSpikeTimes + sta2.rel_procSpikeTimes + sta3.rel_procSpikeTimes + sta4.rel_procSpikeTimes + sta5.rel_procSpikeTimes
        np.save('rel_spikeTimes'+str(run_i),rel_spikeTimes)
        np.save('staInput'+str(run_i),staInput)
        np.save('meanProc'+str(run_i),meanSTA)

        plt.hist(rel_spikeTimes, bins=50)  # plt.hist passes it's arguments to np.histogram
        plt.title('histogram')
        plt.show()

        #sta.publicationPlots(stim)
        #sta.plotData()

def set_pub(ax):
        """Function for setting plotting parameters uniformly for publication

        Input: ax: a axis instance

        Output: none
                function sets the parameters uniformly for plotting

        """


        plt.rc('font', weight='normal')    # bold fonts are easier to see
        plt.rc('lines', lw=10, color='k') # thicker black lines (no budget for color!)
        plt.rc('grid', c='0.5', ls='-', lw=2)  # solid gray grid lines
        plt.rc('savefig', dpi=300)       # higher res outputs
        plt.rc('legend',fontsize=20)


        #http://stackoverflow.com/questions/10947167/changing-font-size-of-legend-title-in-python-pylab-rose-polar-plot
        l = ax.legend()
        """
        #not working
        setp(l.get_texts(), fontsize=14)
        setp(l.get_title(), fontsize=14)
        #error
        File "plot_ISI_final.py", line 78, in <module>
                    tools.set_pub(ax1)
          File "../libs/tools.py", line 416, in set_pub
                    setp(l.get_texts(), fontsize=14)
        AttributeError: 'NoneType' object has no attribute 'get_texts'
        """

        # http://stackoverflow.com/questions/6390393/matplotlib-make-tick-labels-font-size-smaller
        for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(15)
        for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(15)
        print("here")
        ax.legend(loc='upper left')


