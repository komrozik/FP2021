import numpy as np
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import scipy.constants as const

x = np.linspace(0,100,1000)

plt.plot(x,x**2,label ='f(x)')
plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend(loc="best")
plt.savefig("plots/plot1.pdf")
plt.close()