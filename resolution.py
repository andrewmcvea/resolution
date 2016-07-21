#Imported packages
from __future__ import division
import h5py
import math
import numpy as np
from zmq_client import adc_to_voltage
from scipy.stats import norm
from scipy.optimize import fmin
from scipy.optimize import curve_fit
from scipy.special import erf

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
def find_gauss(x, amp, mu, sd):
    gaussian = amp*np.exp(-(((x-mu)**2)/(2*sd**2))) #general gaussian functional form
    return gaussian

def find_fit(x, t1, mu, sigma, amp1, t2, amp2):
    lmb1 = 1/t1
    lmb2 = 1/t2
    exp_fit = (amp1*(lmb1/2)*(math.e**((lmb1/2)*(2*mu + lmb1*(sigma**2) - 2*x)))) * \
              (erf((mu + lmb1*(sigma**2) - x)/(math.sqrt(2)*sigma))) + \
              (amp2*(lmb2/2)*(math.e**((lmb2/2)*(2*mu + lmb2*(sigma**2) - 2*x)))) * \
              (erf((mu + lmb2*(sigma**2) - x)/(math.sqrt(2)*sigma)))
    return exp_fit

if __name__ == '__main__':
    import argparse
    import sys
    import matplotlib.pyplot as plt

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='output', default='out.txt')
    parser.add_argument('filenames', nargs='+', help='input files')
    parser.add_argument('-c', '--chunk', type=int, default=10000)
    args = parser.parse_args()

#Calculate the time resolution for each pulse and creates an array of these values
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
    """
    g_t = [a for a in t if a > 0 and a < 20] #filters out the tail for better gaussian fitting
    l = len(g_t)
    bins = np.linspace(np.min(t),np.max(t),l)

    y, x = np.histogram(g_t, bins) #Outputs the values for a histogram
    center = (x[1:] + x[:-1])/2 #finds the center of the bins

    #Approxiamtions for the parameters
    guess_amp = np.amax(y,axis=0)
    guess_mu = np.average(g_t)
    guess_sd = np.std(g_t)

    popt, pcov = curve_fit(find_gauss, center, y, p0 = [guess_amp,guess_mu,guess_sd]) #finds parameters for gaussian
    fit = find_gauss(center, *popt)

    print 'Guesses:'
    print 'amp=', guess_amp
    print 'mu=', guess_mu
    print 'std=', guess_sd
    print popt
    """

#    t = [a for a in t if a >= 0]
    bins = np.linspace(np.min(t),np.max(t),1000)

    y, x = np.histogram(t, bins) #Outputs the values for a histogram
    center = (x[1:] + x[:-1])/2 #finds the center of the bins

    #Approxiamtions for the parameters
    guess_mu = np.average(t)
    guess_sd = np.std(t)
    guess_amp1 = np.amax(y,axis=0)
    decay_amp = min(y, key=lambda x:abs(x-(guess_amp1/math.e)))
    guess_t1 = float(x[y == decay_amp])
    guess_amp2 = 9*10**4
    guess_t2 = 26

    #finds parameters
    popt, pcov = curve_fit(find_fit, center, y, p0 = [guess_t1,guess_mu,guess_sd,guess_amp1,guess_t2,guess_amp2])
    fit = find_fit(center, *popt)

    print pcov

    print 'Guesses:'
    print 't=', guess_t1
    print 'mu=', guess_mu
    print 'std=', guess_sd
    print 'amp=', guess_amp1
    print '[t1,mu,std,amp1,t2,amp2]:'
    print popt

    plt.hist(t, bins)
    plt.xlabel("Time Resolution")
    plt.plot(center,fit)
    plt.title(popt)
#    plt.yscale('log')
#    plt.ylim(ymin=.9)

plt.show()
