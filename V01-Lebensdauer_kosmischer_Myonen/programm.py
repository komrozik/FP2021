import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import scipy.constants as const
import os
import pandas as pd

if (not os.path.exists('plots')):
    os.mkdir('plots')
# Create a directory named path with numeric mode mode.

def f(x,a,b):
    return a*x+b

def gauß(x,a,b,c):
    return a*np.exp(-(x-b)**2/(2*c**2))

df1 = pd.read_csv('data/data1.txt'
                 ,sep = ','
                 ,lineterminator='\n'
                 #,skiprows=(0)
                 ,header = (0)
                 )
df2 = pd.read_csv('data/data2.txt'
                 ,sep = ','
                 ,lineterminator='\n'
                 ,header = 0
                 ,skiprows = (0,1)
                 ,names = ['Puls','Bin']
                 )
df3 = pd.read_csv('data/data3.txt'
                 ,lineterminator='\n'
                 ,names = ['Count']
                 )


#Bins zuordnen
params1,cov = curve_fit(f,df2['Bin'],df2['Puls'])
errors = np.sqrt(np.diag(cov))
params1_err = unp.uarray(params1,errors)
err = params1_err[0]*df2['Bin']+params1_err[1]

plt.plot(df2['Bin'],f(df2['Bin'],*params1),alpha = 0.5)
plt.plot(df2['Bin'],df2['Puls'],'.')
plt.xlabel('Bin')
plt.ylabel('Pulslänge')
plt.savefig('plots/bins.pdf')

# E-FUNKTION FIT
def e(x,A,B,C):
    return A*np.exp(-B*x)+C

xdata = f(np.arange(512)[df3['Count'] != 0],*params1)
ydata = df3[df3['Count'] != 0]['Count'].to_numpy()

params2,cov = curve_fit(e
                        ,xdata
                        ,ydata
                       )
errors = np.sqrt(np.diag(cov))
params2_err = unp.uarray(params2,errors)
err = params2_err[0]*xdata+params2_err[1]

plot_errorbar = False

if plot_errorbar:
    plt.errorbar(f(np.arange(512)[df3['Count'] != 0],*params1)
                 ,df3[df3['Count'] != 0]['Count'].to_numpy()
                 ,yerr = np.sqrt(df3['Count'])
                 ,fmt = '.'
                 ,label = 'Messwerte'
                 ,alpha = 0.5
                 ,ecolor = 'k'
                )
else:
    plt.plot(xdata
             ,ydata
             ,'.'
             ,label = 'Messwerte'
            )
    
plt.plot(xdata
         ,e(xdata,*params2)
        )
plt.plot(f(np.arange(512)[df3['Count'] == 0],*params1)
         ,df3[df3['Count'] == 0]['Count'].to_numpy()
         ,'.'
         ,label = 'Nicht verwendete Daten (Count = 0)'
        )
plt.ylim((0,150))
# plt.xlim((0,10))
plt.xlabel(f'$\Delta$t/$\mu s$')
plt.ylabel(f'Counts')
plt.legend(loc = 'best')
plt.savefig('plots/e-fkt.pdf')