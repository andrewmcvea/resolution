#Imported packages
from __future__ import division
import h5py
import numpy as np

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

#Converts ADC counts to a voltage
def adc_to_voltage(a):
    return a/1900.77

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='output', default='histo.txt')
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

    t = [a for a in t if a >= 0 and a <= 700]
    bins = np.linspace(np.min(t),np.max(t),1000)

    y, x = np.histogram(t, bins) #Outputs the values for a histogram
    
    data = zip(y, x)
    np.savetxt('outhist.txt', data, delimeter=',')