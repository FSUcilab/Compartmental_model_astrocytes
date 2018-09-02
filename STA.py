import numpy as np
import math
import copy  as cpy
import matplotlib.pyplot as plt
from collections import deque
from tqdm import tqdm
import random as rnd
from operator import add
from operator import sub
from IPython import embed

#--------------------------------------------------------#
#-----Parameters------#
# caS - main calcium trace of interest # caS - if one is interested in the multiple calcium traces this provides the second
# stim_size - number of elements in stimulus vector
# ^^ THIS DEPENDS ON DT!!!! ^^
# thresh - threshold for what is considered a spike
# t_ref - refractory time after a spike when (duration of a spike)


class STA:


    def __init__(self, caS, stim, stim_size=100, thresh=1.2, t_ref=104,name=""):
        self.caS = caS
        self.stimulus = stim
        self.stim_size = stim_size
        self.thresh = thresh
        self.t_ref = t_ref
        self.preSpike_procSpike = []
        self.spikeCount = 0
        self.name = name


    def LinearNonlinearPoisson(self):
        self.spikeTimes = []
        self.rel_procSpikeTimes = []
        self.sta = np.zeros(self.stim_size)
       # sigTrace = deque(maxlen = self.stim_size)
        old_signal = 2. * self.thresh

        # Only want complete ensembles (of size stim_size). So we populate the
        # dequese before processing. Minus one so we can check the
        # stim_size(th) element for a spike
        

        # for s, signal in enumerate(self.caS[:self.stim_size-1]):
        #     # keep track of calcium signal
        #     sigTrace.append(signal)
        #     #ensemble.append(signal)

        print("Finding Spikes...")
        #print self.caS
        #fig3 = plt.figure()        
        #ax3 = fig3.add_subplot(111)
        for s, signal in tqdm(enumerate(self.caS[self.stim_size-1:])):
            index = s + self.stim_size-1
            
            # keep track of calcium signal
            #sigTrace.append(signal)
            #ensemble.append(signal)

            # consider the calcium concentration of note and see if it is above threshold (spiking)
            # if signal >= self.thresh and t > self.t_ref : # signal goes above threshold from below

            if signal >= self.thresh and old_signal < self.thresh:
                   #print s, "s"
                #print signal,"signal"
                #print index, "index"
                #exit()
                self.spikeTimes.append(index)
                self.spikeCount += 1
        
                """
                print
                print 'sta length'
                print len(self.sta)                
                print 'total stim length'
                print len(self.stimulus)
                print 'stim size'
                print self.stim_size
                print 'index'
                print index
                """
                
                self.sta += self.stimulus[index-self.stim_size:index]
                #self.ensembles.append(cpy.copy(sigTrace))

                #FIND THE AVERAGE NUMBER OF SPIKES influencing the soma-----------------#
                #fill this average into a vector that can be used to caclulate an average
                #process spike count leading in the processes

                old_sig = 2. * self.thresh
                count = 0        
                for ss, stim_sig in enumerate(self.stimulus[index-(self.stim_size-1):index]):
                        if stim_sig >= self.thresh and old_sig < self.thresh:
                                count += 1
                                #self.rel_procSpikeTimes.append(ss)
                                #print ss
                        old_sig = cpy.copy(stim_sig)
                self.preSpike_procSpike.append(count)
                #print self.preSpike_procSpike
                #print self.preSpikeStim
                #print self.rel_procSpikeTimes
                #print len(np.linspace(0,2000,1999))
                #print len(self.stimulus[index-(self.stim_size-1):index])
                #ax3.plot(np.linspace(0,self.stim_size,self.stim_size-1),self.stimulus[index-(self.stim_size-1):index])
                #-----------------------------------------------------------------------#

            old_signal = cpy.copy(signal)
        #fig3.show()
        self.sta /= len(self.spikeTimes)

        print
        print("average number of spikes (if stimulus is a process) preceeding soma spike:")
        print
        print(np.mean(self.preSpike_procSpike)        )


    def STAR(self):
        self.spikeTimesR = []
        self.staR = np.zeros(self.stim_size)
        
        old_signal = 2. * self.thresh

        print("Picking random locations")
        randLoc = rnd.sample(xrange(len(self.caS)-self.stim_size),500)

        for index in randLoc:
                index += self.stim_size
                #actual_input = map(add,trace,actual_input)
                self.staR = map(add,self.staR, self.stimulus[index-self.stim_size:index])

        self.staR = [i/len(randLoc) for i in self.staR]

    #----------------------------
    def plotData(self):
        # Need 4 plots in single figure

        #plt.figure(1, figsize=(16,8))
        fig1, axes = plt.subplots(2, 2)

        nrows = 3
        ncols = 2
        plt.subplot(nrows, ncols,1)
        plt.plot(self.caS, 'r')

        ## ERROR: one red line appears to be incorrect.
        plt.subplot(nrows, ncols,2)
        for s in self.spikeTimes:
            plt.plot(range(0, self.stim_size), self.caS[s-self.stim_size:s])
        plt.plot(range(0, self.stim_size), self.sta,'go')

        # Compute histograms of P(s|spike)
        nbins = 50
        plt.subplot(nrows, ncols, 3)
        plt.hist(self.LFR, bins=nbins)

        #compute the linear filter response and use the corresponding histogram for P(s)
        plt.subplot(nrows, ncols,4)
        plt.hist(self.TFR, bins=nbins)

        # P(spike|s) = P(s|spike) / P(s)
        gen = np.zeros(nbins) # size 50
        histogram_LFR = np.histogram(self.LFR, bins=nbins)
        #print histogram_LFR[0]
        #print histogram_LFR[1]
        #print type(histogram_LFR)
        #print np.shape(histogram_LFR)
        histogram_TFR = np.histogram(self.TFR, bins=50)
        #print histogram_TFR[0]
        #print histogram_TFR[1]
        #print type(histogram_TFR)
        #print np.shape(histogram_TFR)
        lfr_freq = histogram_LFR[0]  # size 50
        tfr_freq = histogram_TFR[0]  # histogram_LFR[1]: size 50
        bins     = histogram_TFR[1]
        for i in range(nbins):
            #print "i, lfr_freq[i], tfr_freq[i] = ", i, lfr_freq[i], tfr_freq[i]
            if (tfr_freq[i] == 0):
                gen[i] = 0.
            else:
                gen[i] = lfr_freq[i] / float(tfr_freq[i])
            #print "gen[%d] = %f" % (i, gen[i])
            #print "gen[%d] = %f" % (i, gen[i])
        plt.subplot(nrows, ncols, 5)
        #print "gen: ", np.shape(gen)   # 50
        #print "histogram_TFR: ", np.shape(histogram_TFR)
        #print "histogram_TFR[0]: ", np.shape(histogram_TFR[0]) # 50
        #print "histogram_TFR[1]: ", np.shape(histogram_TFR[1]) # 50

        #plt.plot(histogram_TFR[1], gen, 'r')
        width = bins[1] - bins[0]
        plt.bar(bins[:-1], height=gen, width=width)

        plt.show()
        #-------------------------------------------------------------------------#


    def publicationPlots(self, stimulus=None):
        #plot the stimulus space
        if stimulus is None:
            stimulus = self.caS

        plt.figure(0)
        plt.subplot(211)
        sXpts = []
        sYpts = []
        xpts = []
        ypts = []

        spikeCount = 0
        for s in xrange(self.stim_size, stimulus.size):
            x = stimulus[s-self.stim_size]
            y = stimulus[s-self.stim_size/2]
            if spikeCount < len(self.spikeTimes) and s == self.spikeTimes[spikeCount]:
                sXpts.append(x)
                sYpts.append(y)
                spikeCount += 1
            else:
                xpts.append(x)
                ypts.append(y)

        plt.scatter(xpts, ypts, c='k')
        plt.scatter(sXpts, sYpts, c='g')

        #plot the stiumuls space with the average projected out.. 
        
        plt.subplot(212)
        sXpts = []
        sYpts = []
        xpts = []
        ypts = []

        spikeCount = 0
        for s in xrange(self.stim_size, stimulus.size):
            #using the first value of sta
            x = stimulus[s-self.stim_size] - self.sta[0]
#            print "sta length" , len(self.sta)
 #           print "val 1:", self.stim_size," before spike"
#            print "val 2:", self.stim_size/2," before spike"
#            print "s" , s
            #using the middle value of sta
            y = stimulus[s-self.stim_size/2] - self.sta[self.stim_size/2-1]
            if spikeCount < len(self.spikeTimes) and s == self.spikeTimes[spikeCount]:
                sXpts.append(x)
                sYpts.append(y)
                spikeCount += 1
            else:
                xpts.append(x)
                ypts.append(y)

        plt.scatter(xpts, ypts, c='k')
        plt.scatter(sXpts, sYpts, c='g')



        plt.figure(1)
        plt.hist(self.TFR, bins=500, histtype='step', normed=True, label="TFR")
        plt.hist(self.LFR, bins=500, histtype='step', normed=True, label="LFR")
        plt.legend()

        plt.figure(2)
        for s in self.spikeTimes:
            plt.plot(self.stimulus[s-self.stim_size:s])
        plt.plot(self.sta,color="r", lw=10)

        plt.figure(3)
        plt.plot(self.stimulus/(10*self.stimulus.max())+0.25, color="r")
        plt.plot(self.caS, color='b')
        plt.scatter(self.spikeTimes, self.caS[self.spikeTimes], color='k')
        


#        plt.show()

