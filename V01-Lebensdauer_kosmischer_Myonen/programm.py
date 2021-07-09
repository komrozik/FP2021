import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import scipy.constants as const
import os

if (not os.path.exists('plots')):
    os.mkdir('plots')
# Create a directory named path with numeric mode mode.

data1 = np.genfromtxt('data/data1.txt',)
delay = data1[:,0]
rate = data1[:,1]/10

plt.plot(delay,rate,label = 'Peak')
plt.ylabel('Rate')
plt.xlabel('Delay in ns')
plt.legend(loc = 'best')
plt.savefig('plots/Delay_plot.pdf')

data2 = np.genfromtxt('data/data2.txt',)
pulsabstand = data1[:,0]
bins = data1[:,1]