# Routines to identify spike locations
import numpy as np

def findSpikes(time, signal, thresh):
    # Spike: signal > thresh and signal increasing (ndarray)
    # time array: real numbers: 0, dt, 2*dt, etc. (ndarray)
    # Return a list of global spike indexes into the "time" array

    # indexes of soma spike > threshold
    # where returns a n-tuple (results). Extract 0th element: results
    where = np.where(np.diff(np.sign(signal-thresh)) > 0)[0]
    return where



