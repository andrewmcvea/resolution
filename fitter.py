#Imported packages
from __future__ import division
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from lmfit import minimize, Parameters, Parameter, report_fit

print "start program"

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
    c = params['n11'].value
    tau = params['n12'].value
    rise = params['n13'].value

    lmb1 = 1/t1
    lmb2 = 1/t2
    lmb3 = 1/t3
    # fitting function
    exp_fit = (amp1*(lmb1/2)*(math.e**((lmb1/2)*(2*mu + lmb1*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb1*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (amp2*(lmb2/2)*(math.e**((lmb2/2)*(2*mu + lmb2*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb2*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (amp3*(lmb3/2)*(math.e**((lmb3/2)*(2*mu + lmb3*(sigma**2) - 2*x)))) * \
              (1-(erf((mu + lmb3*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
              (gamp*np.exp(-(((x-gmu)**2)/(2*gsd**2))) + c) + \
              (math.e**(-x/tau) - math.e**(-x/rise))/(tau/rise)
    return exp_fit**2 - y**2

y, x = np.loadtxt('outhist.txt',unpack=True)

#Approxiamtions for the parameters
guess_mu = np.average(y)
guess_sd = np.std(y)
guess_amp1 = np.amax(y,axis=0)
decay_amp = min(y, key=lambda x:abs(x-(guess_amp1/math.e)))
guess_t1 = float(x[y == decay_amp])

print 'Guesses:'
print 'amp=', guess_amp1
print 'mu=', guess_mu
print 'std=', guess_sd
print 't1=', guess_t1
print 'entries=', len(y)

params = Parameters()
params.add('n0', value= 4.929, vary=True, min=1, max=5)
params.add('n1', value= 11.84, vary=True, min=0, max=(guess_mu*2))
params.add('n2', value= .86, vary=True, min=.1, max=1.8)
params.add('n3', value= 89168.95, vary=True, min=(guess_amp1/10), max=100000)
params.add('n4', value= 19, vary=True, min=5, max=50)
params.add('n5', value= 18000, vary=True, min=10000, max=20000)
params.add('n6', value= 122.85, vary=True, min=50, max=500)
params.add('n7', value= 2475.64, vary=True, min=1000, max=5000)
params.add('n8', value= 244.27, vary=True, min=1, max=300)
params.add('n9', value= 38.9, vary=True, min=35, max=40)
params.add('n10', value= 5.228, vary=True,  min=0, max=10)
params.add('n11', value= 1.75, vary=True, min=1, max=20)
params.add('n12', value= 1.3, vary=True, min=.01, max=100)
params.add('n13', value= 1.75, vary=True, min=.01, max=1000)

result = minimize(find_fit, params, method='nelder', args=(x, y))

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
c = result.params['n11']
tau = result.params['n12']
rise = result.params['n13']
chi = result.chisqr
dof = (len(x) - 14)
chi2 = chi / dof
lmb1 = 1/t1
lmb2 = 1/t2
lmb3 = 1/t3

x_fit = (amp1*(lmb1/2)*(math.e**((lmb1/2)*(2*mu + lmb1*(sigma**2) - 2*x)))) * \
          (1-(erf((mu + lmb1*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
          (amp2*(lmb2/2)*(math.e**((lmb2/2)*(2*mu + lmb2*(sigma**2) - 2*x)))) * \
          (1-(erf((mu + lmb2*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
          (amp3*(lmb3/2)*(math.e**((lmb3/2)*(2*mu + lmb3*(sigma**2) - 2*x)))) * \
          (1-(erf((mu + lmb3*(sigma**2) - x)/(math.sqrt(2)*sigma)))) + \
          (gamp*np.exp(-(((x-gmu)**2)/(2*gsd**2))) + c) + \
          (math.e**(-x/tau) - math.e**(-x/rise))/(tau/rise)

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
print 'c=', c
print 'tau=', tau
print 'Rise Time=', rise
print 'Chi Squared=', chi2

bins = np.linspace(np.min(x),np.max(x),len(x))

plt.xlabel("Time Resolution")
plt.step(x, y)
plt.plot(bins, x_fit)
plt.yscale('log')
plt.ylim(ymin=.9)

plt.show()

