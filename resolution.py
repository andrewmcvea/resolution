#Imported packages
from __future__ import division
import math
import h5py
import numpy as np
from scipy.stats import norm
from scipy.optimize import fmin
from scipy.special import erf
from lmfit import minimize, Parameters, Parameter, report_fit

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

#Finds the best fit using the given equation
def find_fit(params, x, y):
    t1 = params['n0'].value
    mu = params['n1'].value
    sigma = params['n2'].value
    amp1 = params['n3'].value
    t2 = params['n4'].value
    amp2 = params['n5'].value
    t3 = params['n6'].value
    amp3 = params['n7'].value
    gamp = params['n8'].value
    gmu = params['n9'].value
    gsd = params['n10'].value

    lmb1 = 1/t1
    lmb2 = 1/t2
    lmb3 = 1/t3
    exp_fit = (amp1*(lmb1/2)*(math.e**((lmb1/2)*(2*mu + lmb1*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb1*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (amp2*(lmb2/2)*(math.e**((lmb2/2)*(2*mu + lmb2*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb2*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (amp3*(lmb3/2)*(math.e**((lmb3/2)*(2*mu + lmb3*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb3*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              gamp*np.exp(-(((x-gmu)**2)/(2*gsd**2)))
    return abs(exp_fit - y)

#Converts ADC counts to a voltage
def adc_to_voltage(a):
    return a/1900.77

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

    t = [a for a in t if a >= 0 and a <= 200]
    bins = np.linspace(np.min(t),np.max(t),1000)

    y, x = np.histogram(t, bins) #Outputs the values for a histogram
    center = (x[1:] + x[:-1])/2 #finds the center of the bins
    
    #Approxiamtions for the parameters
    guess_mu = np.average(t)
    guess_sd = np.std(t)
    guess_amp1 = np.amax(y,axis=0)
    decay_amp = min(y, key=lambda x:abs(x-(guess_amp1/math.e)))
    guess_t1 = float(x[y == decay_amp])

    print 'Guesses:'
    print 'amp=', guess_amp1
    print 'mu=', guess_mu
    print 'std=', guess_sd
    print 't1=', guess_t1
    print 'entries=', len(t)

    params = Parameters()
    params.add('n0', value= guess_t1, min=1, max=100)
    params.add('n1', value= guess_mu, min=(guess_mu*.5), max=(guess_t1*2))
    params.add('n2', value= 1.2, min=.8, max=1.5)
    params.add('n3', value= guess_amp1, min=(guess_amp1*.5), max=(guess_amp1*1.$
    params.add('n4', value= 20, min=1, max=1000)
    params.add('n5', value= 7, min=0, max=1000000)
    params.add('n6', value= 20, min=1, max=10000)
    params.add('n7', value= 7, min=0, max=1000)
    params.add('n8', value= 7, min=0, max=1000)
    params.add('n9', value= 38, min=30, max=50)
    params.add('n10', value= 7, min=0, max=20)
    params.add('n11', value= 10, min=1, max=20)

    result = minimize(find_fit, params, args=(center, y))

    t1 = result.params['n0']
    mu = result.params['n1']
    sigma = result.params['n2']
    amp1 = result.params['n3']
    t2 = result.params['n4']
    amp2 = result.params['n5']
    t3 = result.params['n6']
    amp3 = result.params['n7']
    gamp = result.params['n8']
    gmu = result.params['n9']
    gsd = result.params['n10']

    lmb1 = 1/t1
    lmb2 = 1/t2
    lmb3 = 1/t3
    x_fit = (amp1*(lmb1/2)*(math.e**((lmb1/2)*(2*mu + lmb1*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb1*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (amp2*(lmb2/2)*(math.e**((lmb2/2)*(2*mu + lmb2*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb2*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (1-(erf((mu + lmb3*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              gamp*np.exp(-(((x-gmu)**2)/(2*gsd**2)))

    print 'Values:'
    print 't1=', t1
    print 'mu=', mu
    print 'std=', sigma
    print 'amp1=', amp1
    print 't2=', t2
    print 'amp2=', amp2
    print 't3=', t3
    print 'amp3=', amp3
    print 'gamp=', gamp
    print 'gmu=', gmu
    print 'gsd=', gsd

    plt.hist(t, bins)
    plt.xlabel("Time Resolution")
    plt.plot(bins, x_fit)
    plt.yscale('log')
    plt.ylim(ymin=.9)

plt.show()
