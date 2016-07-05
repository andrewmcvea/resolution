#Imported packages               
from __future__ import division
import h5py
import numpy as np
from zmq_client import adc_to_voltage
import sys
from scipy.stats import norm
from scipy.optimize import fmin
import heapq

print "start program"

#finds time for each PMT (40% on peak rise time)
def find_time(v):
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

#creates a gaussian fit for the data
def gauss(v,bins):
    avg = np.mean(v)
    sig = np.std(v)
    bincenters = (bins[1:] + bins[:-1])/2
    hist, _ = np.histogram(v,bins)
    hist_sigma = hist.copy()
    hist_sigma[hist == 0] = 1
    def foo(args):
        mu, std, c = args
        pdf = c*norm.pdf(bincenters,mu,std)
        return np.sum((pdf-hist)**2/hist_sigma)
    return fmin(foo,[avg,sig,1000])

#finds the amplitude and filters out values below a certain value
def find_amp(v):
        amplitude = np.min(v,axis=1)
        filteramp = amplitude[amplitude < -200]
        return amplitude

if __name__ == '__main__':
    import argparse
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="input filename")
    args = parser.parse_args()

    f = h5py.File(args.filename)
    dset = f['c1'][:100000]
    amp = find_amp(dset)
    f_dset = dset[amp < -200]

    t_tr = find_time(f_dset)
    t = t_tr.copy()
    t *= 0.5 #ns conversion
    t -= np.mean(t)

    dset2 = f['c2'][:100000]
    amp2 = find_amp(dset2)
    size = len(f_dset)
    famp2 = heapq.nlargest(size, amp2) 
    f_dset2 = dset2[famp2]

    ts = find_time(f_dset2)
    t2 = ts.copy()
    t2 *= 0.5 #ns conversion
    t2 -= np.mean(t2)

    res = t2 - t

    bins = np.arange(-50,50,0.5)
    x = np.linspace(-50,50,100000)
    result = gauss(res,bins)
    mu, std, c = result

    plt.hist(res, bins)
    plt.title("Time Resolution Histogram")
    plt.xlabel("Time Resolution")
    plt.plot(x,c*norm.pdf(x,mu,std))

plt.show()
