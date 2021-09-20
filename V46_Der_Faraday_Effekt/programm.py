import numpy as np
import matplotlib
font = {'size': 11.0}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import uncertainties
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import scipy.constants

df1 = pd.read_csv('data/data1.txt'
                 ,sep = ','
                 ,lineterminator='\n'
                 ,skiprows=(3)
                 ,header = 0
                 ,names = ['lambda_0','theta_1','theta_1_min','theta_2','theta_2_min']
                 )
df2 = pd.read_csv('data/data2.txt'
                 ,sep = ','
                 ,lineterminator='\n'
                 ,skiprows=(3)
                 ,header = 0
                 ,names = ['lambda_0','theta_1','theta_1_min','theta_2','theta_2_min']
                 )
df3 = pd.read_csv('data/data3.txt'
                 ,sep = ','
                 ,lineterminator='\n'
                 ,skiprows=(3)
                 ,header = 0
                 ,names = ['lambda_0','theta_1','theta_1_min','theta_2','theta_2_min']
                 )
df4 = pd.read_csv('data/data4.txt'
                 ,sep = ','
                 ,lineterminator='\n'
                 ,skiprows=(0)
                 ,header = (0)
                 #,names = ['Abstand','B-Feld']
                 )

df1['theta1'] = (df1['theta_1'] + 1/60*df1['theta_1_min'])/180 * np.pi
df1['theta2'] = (df1['theta_2'] + 1/60*df1['theta_2_min'])/180 * np.pi
df1 = df1.drop(columns = ['theta_1','theta_1_min','theta_2','theta_2_min'])
df1['theta'] = np.abs(0.5*(df1['theta1']-df1['theta2']))

df2['theta1'] = (df2['theta_1'] + 1/60*df2['theta_1_min'])/180 * np.pi
df2['theta2'] = (df2['theta_2'] + 1/60*df2['theta_2_min'])/180 * np.pi
df2 = df2.drop(columns = ['theta_1','theta_1_min','theta_2','theta_2_min'])
df2['theta'] = np.abs(0.5*(df2['theta1']-df2['theta2']))

df3['theta1'] = (df3['theta_1'] + 1/60*df3['theta_1_min'])/180 * np.pi
df3['theta2'] = (df3['theta_2'] + 1/60*df3['theta_2_min'])/180 * np.pi
df3 = df3.drop(columns = ['theta_1','theta_1_min','theta_2','theta_2_min'])
df3['theta'] = np.abs(0.5*(df3['theta1']-df3['theta2']))

df1['theta_norm'] = df1['theta']/1.296
df2['theta_norm'] = df2['theta']/1.36
df3['theta_norm'] = df3['theta']/5.11

df1['theta_frei'] = np.abs(df1['theta_norm']-df3['theta_norm'])
df2['theta_frei'] = np.abs(df2['theta_norm']-df3['theta_norm'])

#------------------
#Magnetfeld Vermessen: 

def gauß(x,a,b,c):
    return a*np.exp(-((x-b)/(2*c))**(2))

params0,cov = curve_fit(gauß,df4['Abstand'],df4['B-Feld'])
errors0 = np.sqrt(np.diag(cov))
params0_err = unp.uarray(params0,errors0)
err0 = params0_err[0]*df1['lambda_0']**2+params0_err[1]

plt.figure(figsize=(6.4,3.96),dpi=300)
plt.plot(df4['Abstand']
         ,df4['B-Feld']
         ,marker = 'x'
         ,ls = ''
         ,label = 'Messdaten'
        )
plt.plot(np.linspace(-11,50,100)
         ,gauß(np.linspace(-11,50,100),*params0)
         ,ls = '--'
         ,alpha = 0.5
         ,label = 'Messdaten'
        )
plt.hlines(y = 0.5*gauß(np.linspace(-11,50,100),*params0)[44]
           ,xmin = 6.3
           ,xmax = 26
           ,label = 'Halbwertsbreite'
          )
plt.axvline(x = np.linspace(-11,50,100)[44]
            ,color = 'c'
            ,alpha = 0.4
            ,label = 'Maximum'
            )
plt.ylabel(f'$B$ in mT')
plt.xlabel('$d$ in mm')
plt.title('Magnetische Flussdichte')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/B_Feld.pdf',bbox_inches = "tight")
plt.close()

#-----------------
#Normierte Drehwinkel:
plt.figure(figsize=(6.4,3.96),dpi=300)
plt.plot(df1['lambda_0']**2
         ,np.abs(df1['theta_norm'])
         ,ls = '-'
         ,alpha = 0.2
         ,color = '#1f77b4'
        )
plt.plot(df1['lambda_0']**2
         ,np.abs(df1['theta_norm'])
         ,marker = 'x'
         ,ls = ''
         ,label = r'Messdaten N = $2,8 \cdot 10^{18} cm^{-3}$'
         ,color = '#1f77b4'
        )
plt.plot(df2['lambda_0']**2
         ,np.abs(df2['theta_norm'])
         ,ls = '-'
         ,alpha = 0.2
         ,color = '#ff7f0e'
        )
plt.plot(df2['lambda_0']**2
         ,np.abs(df2['theta_norm'])
         ,marker = 'x'
         ,ls = ''
         ,label = r'Messdaten N = $1,4 \cdot 10^{18} cm^{-3}$'
         ,color = '#ff7f0e'
        )
plt.plot(df3['lambda_0']**2
         ,np.abs(df3['theta_norm'])
         ,ls = '-'
         ,alpha = 0.2
         ,color = '#2ca02c'
        )
plt.plot(df3['lambda_0']**2
         ,np.abs(df3['theta_norm'])
         ,marker = 'x'
         ,ls = ''
         ,label = 'Messdaten reine Probe'
         ,color = '#2ca02c'
        )
plt.ylabel(r'$\theta_{norm}$ in mT')
plt.xlabel(r'$\lambda^2$ in $\mu m^2$')
plt.title('Normierte Drehwinkel')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/Drehwinkel.pdf',bbox_inches = "tight")
plt.close()

#-------------
#Drehung durch Dotierung

def f(x,m,b):
    return m*x+b

params1,cov = curve_fit(f,df1['lambda_0']**2,df1['theta_frei'])
errors1 = np.sqrt(np.diag(cov))
params1_err = unp.uarray(params1,errors1)
err1 = params1_err[0]*df1['lambda_0']**2+params1_err[1]


plt.figure(figsize=(6.4,3.96),dpi=300)
plt.plot(df1['lambda_0']**2
         ,df1['theta_frei']
         ,ls = '-'
         ,alpha = 0.2
         ,color = '#1f77b4'
        )
plt.plot(df1['lambda_0']**2
         ,df1['theta_frei']
         ,marker = 'x'
         ,ls = ''
         ,label = r'Messdaten'
         ,color = '#1f77b4'
        )
plt.plot(np.linspace(1,7,100)
         ,f(np.linspace(1,7,100),params1[0],params1[1])
         #,marker = 'x'
         ,ls = '--'
         ,label = r'Fit'
         #,color = '#1f77b4'
        )
plt.ylabel(r'$\theta$ in mT')
plt.xlabel(r'$\lambda^2$ in $\mu m^2$')
plt.title('Drehung durch Dotierung \n N = $2,8 \cdot 10^{18} cm^{-3}$')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/Dotierung1.pdf',bbox_inches = "tight")
plt.close()


params2,cov = curve_fit(f,df2['lambda_0']**2,df2['theta_frei'])
errors2 = np.sqrt(np.diag(cov))
params2_err = unp.uarray(params2,errors2)
err2 = params2_err[0]*df2['lambda_0']**2+params2_err[1]

params3,cov = curve_fit(f,df2.drop([1,7])['lambda_0']**2,df2.drop([1,7])['theta_frei'])
errors3 = np.sqrt(np.diag(cov))
params3_err = unp.uarray(params3,errors3)
err3 = params3_err[0]*df2['lambda_0']**2+params3_err[1]

plt.figure(figsize=(6.4,3.96),dpi=300)
plt.plot(df2['lambda_0']**2
         ,df2['theta_frei']
         ,ls = '-'
         ,alpha = 0.2
         ,color = '#ff7f0e'
        )
plt.plot(df2['lambda_0']**2
         ,df2['theta_frei']
         ,marker = 'x'
         ,ls = ''
         ,label = r'Messdaten'
         ,color = '#ff7f0e'
        )
plt.plot(np.linspace(1,7,100)
         ,f(np.linspace(1,7,100),params2[0],params2[1])
         #,marker = 'x'
         ,ls = '--'
         ,label = r'Fit'
         ,color = '#ff7f0e'
        )
plt.plot(np.linspace(1,7,100)
         ,f(np.linspace(1,7,100),params3[0],params3[1])
         #,marker = 'x'
         ,ls = ':'
         ,label = r'Fit 2 (Wert 2 und 8 entfernt)'
         ,color = '#ff7f0e'
        )
plt.ylabel(r'$\theta$ in mT')
plt.xlabel(r'$\lambda^2$ in $\mu m^2$')
plt.title('Drehung durch Dotierung \n N = $1,4 \cdot 10^{18} cm^{-3}$')
plt.legend(loc = 'best')
plt.tight_layout()
plt.savefig('plots/Dotierung2.pdf',bbox_inches = "tight")
plt.close()