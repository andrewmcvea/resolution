# Imported packages for histograms
from __future__ import division
import h5py
import numpy as np
from zmq_client import adc_to_voltage #may need this later
import sys 
from scipy.optimize import fmin
from scipy.stats import norm

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

#Finds the values to be used in the Gaussian curve
def gauss(v,bins):
    avg = np.mean(v)
    std = np.std(v)
    bincenters = (bins[1:] + bins[:-1])/2
    hist, _ = np.histogram(v,bins)
    hist_sigma = hist.copy()
    hist_sigma[hist == 0] = 1
#I should define foo more precisely later 
   def foo(args):
        mu, std, c = args
        pdf = c*norm.pdf(bincenters,mu,std)
        return np.sum((pdf-hist)**2/hist_sigma)

    return fmin(foo,[avg,std,1000])


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
    
    bins = np.arange(-50,50,0.5)
    x = np.linspace(-50,50,1000)
    result = gauss(t,bins)
    mu, std, c = result

    plt.hist(t, bins)
    plt.title("Histogram of Time Resolution")
    plt.xlabel("Time Resolution")
    plt.plot(x,c*norm.pdf(x,mu,std))


plt.show()
 

