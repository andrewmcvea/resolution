# Imported packages for histograms
from __future__ import division
import h5py
import numpy as np
from zmq_client import adc_to_voltage
import sys
from scipy.optimize import fmin

import argparse
import matplotlib.pyplot as plt

print "start program"

#finds time resolution (40% on peak rise time)
def find_res(v):
    t = np.empty(v.shape[0],dtype=float)
    for i in range(len(v)):
        if i % 100 == 0:
            print "\r%i/%i" % (i+1,len(v)),
            sys.stdout.flush()
        j = np.argmin(v[i])
        disc = v[i,j]*0.4
        while v[i,j] < disc and j > 0:
            j -= 1
        t[i] = j + (disc - v[i,j])/(v[i,j+1] - v[i,j])
    print()
    return t

if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="input filename")
    args = parser.parse_args()

    f = h5py.File(args.filename)
    dset = f['c1'][:100000]
    t1 = find_res(dset)
    t = t1.copy()
    t *= 0.5
    t -= np.mean(t)

    plt.hist(t, bins=np.arange(-50,50,0.5))
    plt.title("Histogram of Time Resolution")
    plt.xlabel("Time Resolution")

plt.show()
 

