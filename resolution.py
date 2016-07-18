#Imported packages
from __future__ import division
import h5py
import math
import numpy as np
from zmq_client import adc_to_voltage
from scipy.stats import norm
from scipy.optimize import fmin
from scipy.optimize import curve_fit

print "start program"

def get_times(y, fraction=0.4):
    """
    Returns pulse times in `y` by looking for the pulse
    to cross a constant fraction `fraction` of the pulse height in each
    waveform. `y` should be a 2 dimensional array with shape (N,M) where
    N is the number of waveforms, and M is the number of samples per
    waveform.
    """
    # samples below threshold
    mask1 = y > np.min(y,axis=-1)[:,np.newaxis]*fraction
    # samples before the minimum
    mask2 = np.arange(y.shape[1]) < np.argmin(y,axis=-1)[:,np.newaxis]

    # right side of threshold crossing
    r = y.shape[1] - np.argmax((mask1 & mask2)[:,::-1], axis=-1)
    r[r == 0] = 1
    r[r == y.shape[1]] = y.shape[1] - 1
    l = r - 1

    yl = y[np.arange(y.shape[0]),l]
    yr = y[np.arange(y.shape[0]),r]

    return (np.min(y,axis=-1)*fraction - yl)/(yr-yl) + l

def fft_filter(y, dt, cutoff=500e6):
    """
    Filter the array `y` by removing frequency components above
    `cutoff` in Hz.
    """
    out = np.fft.rfft(y)
    freq = np.fft.rfftfreq(y.shape[1], d=dt)
    out[:,freq > cutoff] = 0
    return np.fft.irfft(out)

#Fits the data with the given fittting function
def find_fit(x,amp, mu, sd):
    gaussian = amp*np.exp(-(((x-mu)**2)/(2*sd**2)))
    return gaussian

if __name__ == '__main__':
    import argparse
    import sys
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='output', default='out.txt')
    parser.add_argument('filenames', nargs='+', help='input files')
    parser.add_argument('-c', '--chunk', type=int, default=10000)
    args = parser.parse_args()

    t = []
    for filename in args.filenames:
        with h5py.File(filename) as f:
            for i in range(0, f['c1'].shape[0], args.chunk):
                y1 = adc_to_voltage(f['c1'][i:i+args.chunk])
                y2 = adc_to_voltage(f['c2'][i:i+args.chunk])
                # only accept c2 events below -10 mV
                mask = np.min(y2,axis=-1) < -10e-3
                #tony's new line
                mask &= np.min(y1, axis=-1) < -100e-3
                y1 = y1[mask]
                y2 = y2[mask]
                t1 = get_times(y1)*0.5 # sample -> ns
                t2 = get_times(y2)*0.5
                res = t2 - t1
                t.extend(res)

    g_t = [a for a in t if a > 0 and a < 25]
    l = len(g_t)
    bins = np.linspace(np.min(t),np.max(t),l)

    y, x = np.histogram(g_t, bins) #Outputs the values for a histogram
    center = (x[1:] + x[:-1])/2 #finds the center of the bins

    #Approxiamtions for the parameters
    guess_amp = np.amax(y,axis=0)
    guess_mu = np.average(g_t)
    guess_sd = np.std(g_t)

    popt, pcov = curve_fit(find_fit, center, y, p0 = [guess_amp,guess_mu,guess_sd])
    fit = find_fit(center, *popt)

    print 'Guesses:'
    print 'amp=', guess_amp
    print 'mu=', guess_mu
    print 'std=', guess_sd
    print popt

    plt.hist(t, bins)
    plt.xlabel("Time Resolution")
    plt.plot(center,fit)
    plt.title(r'$\sigma$ = %.2f' % guess_sd)
    plt.yscale('log')
    plt.ylim(ymin=.9)

plt.show()
