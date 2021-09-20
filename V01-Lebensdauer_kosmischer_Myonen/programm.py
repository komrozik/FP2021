import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib
font = {'size': 11.0}
matplotlib.rc('font', **font)
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

plt.figure(figsize=(6.4,3.96),dpi=300)
plt.plot(df2['Bin']
         ,df2['Puls']
         ,'.'
         ,label = 'Meßdaten'
        )
plt.plot(df2['Bin']
         ,f(df2['Bin'],*params1)
         ,alpha = 0.5
         ,label = 'Fit einer linearen Funktion'
        )
plt.xlabel('Bin')
plt.ylabel(f'Pulslänge $\Delta$t/$\mu s$')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/bins.pdf',bbox_inches = "tight")
plt.close()

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

plt.figure(figsize=(6.4,3.96),dpi=300)
plt.plot(xdata
         ,ydata
         ,'.'
         ,label = 'Messwerte'
        )
plt.plot(xdata
         ,e(xdata,*params2)
         ,label = 'Fit der e-Funktion'
        )
plt.plot(f(np.arange(512)[df3['Count'] == 0],*params1)
         ,df3[df3['Count'] == 0]['Count'].to_numpy()
         ,'.'
         ,label = 'Nicht verwendete Daten (Count = 0)'
        )
plot_errorbar = True
if plot_errorbar:
    plt.errorbar(f(np.arange(512)[df3['Count'] != 0],*params1)
                     ,df3[df3['Count'] != 0]['Count'].to_numpy()
                     ,yerr = np.sqrt(df3['Count'][df3['Count'] != 0])
                     ,fmt = '.'
                     ,label = 'Meßwerte'
                     ,ecolor = 'grey'
                     ,alpha = 0.3
                     ,markersize=0
                     ,capsize=3
                    )
plt.ylim((0,150))
# plt.xlim((0,10))
plt.xlabel(f'$\Delta$t/$\mu s$')
plt.ylabel(f'Counts')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/e-fkt.pdf',bbox_inches = "tight")
plt.close()

# Rate für 20 Hz finden 
def gauß(x,a,b,c):
    return a*np.exp(-(x-b)**2/(2*c**2))

params3,cov = curve_fit(gauß,df1['Delay'],df1['Rate'])
errors = np.sqrt(np.diag(cov))
params3_err = unp.uarray(params3,errors)
err = params1_err[0]*df1['Delay']+params1_err[1]

plt.figure(figsize=(6.4,3.96),dpi=300)
plt.errorbar(df1['Delay']
             ,df1['Rate']
             ,yerr = np.sqrt(df1['Rate'])
             ,fmt = '.'
             ,label = 'Fehler der Messwerte'
             ,ecolor = 'grey'
             ,capsize=3
            )
plt.plot(df1['Delay']
         ,gauß(df1['Delay'],*params3)
         ,label = 'Fit einer Gauß-Funktion'
        )
plt.hlines(y = 0.5*gauß(df1['Delay'],*params3)[8]
           ,xmin = -5.55
           ,xmax = 7.45
           ,label = 'Halbwertsbreite'
          )
plt.axvline(x = 1
            ,color = 'c'
            ,alpha = 0.4
            ,label = 'Maximum'
            )
plt.xlabel(f'$\Delta t_1 - \Delta t_2$')
plt.ylabel(f'Counts/10s')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/rate_20Hz.pdf',bbox_inches = "tight")
plt.close()
print('Durchgelaufen!')