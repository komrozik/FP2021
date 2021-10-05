import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
import os

path = './data/Wasserstoffatom/Ohne_Ring'
files = os.listdir(path)
D2=[4.896,4.210,3.480,2.803,0.771,0.712,1.789,3.024,5.671,7.205,10.182,9.851,15.075,17,18.952,17.812,21.649,19.413,21.426]#2,301kHz
D3=[0]*18#3,694kHz
D4=[0]*18#4,979kHz
D7=[0]*18#7,433kHz
a=np.arange(0,185,10)#winkel


#Polarplots
#-----------------------------------------
#Read in
path = './data/Wasserstoffatom/Ohne_Ring/Spectrums'
files = sorted(os.listdir(path))
D2=[] #[4.896,4.210,3.480,2.803,0.771,0.712,1.789,3.024,5.671,7.205,10.182,9.851,15.075,17,18.952,17.812,21.649,19.413,21.426]#2,301kHz
D3=[]#3,694kHz
D4=[]#4,979kHz
D7=[]#7,433kHz
a=np.arange(0,185,10)#winkel

#print(files)
#print(a)
for f in files:
    temp = np.genfromtxt(f"{path}/{f}",delimiter=" ")
    temp_f = temp[:,0]
    temp_D = temp[:,1]
    index2 = np.where(temp_f == 2300)
    index3 = np.where(temp_f == 3690)
    index4 = np.where(temp_f == 4970)
    index7 = np.where(temp_f == 7430)
    D2.append(temp_D[index2][0])
    D3.append(temp_D[index3][0])
    D4.append(temp_D[index4][0])
    D7.append(temp_D[index7][0])
arrays = [D2,D3,D4,D7]

#All together
fig,ax = plt.subplots(2,2,figsize=(6.4,6.4),dpi=300,subplot_kw={'projection': 'polar'})
count = 0
for i in (0,1):
    for j in (0,1):
        rho = np.concatenate([arrays[count],arrays[count][::-1],arrays[count],arrays[count][::-1]])
        r = np.linspace(1, 2, len(rho))
        theta = 2 * np.pi * r
        ax[i,j].set_rlim(0,35)
        ax[i,j].set_rticks([5, 15, 25])  # Less radial ticks
        ax[i,j].set_rlabel_position(-22.5)  # Move radial labels away from plotted line
        ax[i,j].plot(theta,rho)
        count = count+1
fig.tight_layout()
plt.close()

#Einzelne Plots
def cre_polar(data,name):
    rho = np.concatenate([data,data[::-1],data,data[::-1]])
    r = np.linspace(1, 2, len(rho))
    theta = 2 * np.pi * r
    a=np.arange(0,185,10)
    plt.polar(theta,rho)
    plt.savefig("plots/Hatom/polar_"+f"{name}.pdf")
    plt.close()
cre_polar(D2,"2300")
cre_polar(D3,"3960")
cre_polar(D4,"4970")
cre_polar(D7,"7430")

#Wasserstoffatom Zustandsaufspaltung 2095 2265
path = './data/Wasserstoffatom/Ring_9mm_winkelabh'
files = sorted(os.listdir(path))
D=[]
a=np.arange(0,185,10)#winkel

for f in files:
    temp = np.genfromtxt(f"{path}/{f}",delimiter=" ")
    temp_f = temp[:,0]
    temp_D = temp[:,1]
    index = np.where(temp_f == max(temp_f))
    print(index)
    #index = np.where(temp_f == 2095)
    D.append(temp_D[index][0])
cre_polar(D,"2095_9mm")


#Wasserstoff Molekül

#Read in
path = './data/Wasserstoffmolekuel/winkelverteilung_15mm/Spectrums'
files = sorted(os.listdir(path))
D2=[] #[4.896,4.210,3.480,2.803,0.771,0.712,1.789,3.024,5.671,7.205,10.182,9.851,15.075,17,18.952,17.812,21.649,19.413,21.426]#2,301kHz
D3=[]#3,694kHz
D4=[]#4,979kHz
D7=[]#7,433kHz
a=np.arange(0,185,10)#winkel

#print(files)
#print(a)
for f in files:
    temp = np.genfromtxt(f"{path}/{f}",delimiter=" ")
    temp_f = temp[:,0]
    temp_D = temp[:,1]
    index2 = np.where(temp_f == 2297)
    index3 = np.where(temp_f == 2304)
    index4 = np.where(temp_f == 2416)
    D2.append(temp_D[index2][0])
    D3.append(temp_D[index3][0])
    D4.append(temp_D[index4][0])
arrays = [D2,D3,D4,D7]
cre_polar(D2,"mol_2297")
cre_polar(D3,"mol_2304")
cre_polar(D4,"mol_2416")



#Zustandsaufspaltung
f=np.genfromtxt("data/Wasserstoffatom/Ring/Spektrum_Ring_3mm.txt",skip_header=1,usecols=(0))
D=np.genfromtxt("data/Wasserstoffatom/Ring/Spektrum_Ring_3mm.txt",skip_header=1,usecols=(1))
plt.title("H-Atom 180° 3mm")
plt.plot(f,D,label="H-Atom")
vlines(2292, 0, 23.5, colors='r', linestyles='--',alpha=0.6)
vlines(2229, 0, 13.5, colors='r', linestyles='--',alpha=0.6,label=r"$\Delta f=63\,$Hz")
plt.ylabel(r"$P\,/\,$will. Einheit")
plt.xlabel(r"$f\,/\,$Hz")
plt.legend(loc="best")
plt.savefig("plots/Hatom/zustandsaufspaltung.pdf")
plt.close()


f=np.genfromtxt("data/Wasserstoffatom/Ring/Spektrum_Ring_6mm.dat",skip_header=1,usecols=(0))
D=np.genfromtxt("data/Wasserstoffatom/Ring/Spektrum_Ring_6mm.dat",skip_header=1,usecols=(1))
plt.title("H-Atom 180° 6mm")
plt.plot(f,D,label="H-Atom")
vlines(2276, 0, 24.1, colors='r', linestyles='--',alpha=0.6)
vlines(2164, 0, 12.1, colors='r', linestyles='--',alpha=0.6,label=r"$\Delta f=112\,$Hz")
plt.ylabel(r"$P\,/\,$will. Einheit")
plt.xlabel(r"$f\,/\,$Hz")
plt.legend(loc="best")
plt.savefig("plots/Hatom/zustandsaufspaltung_6.pdf")
plt.close()

f=np.genfromtxt("data/Wasserstoffatom/Ring/Spektrum_Ring_9mm.dat",skip_header=1,usecols=(0))
D=np.genfromtxt("data/Wasserstoffatom/Ring/Spektrum_Ring_9mm.dat",skip_header=1,usecols=(1))
plt.title("H-Atom 180° 9mm")
plt.plot(f,D,label="H-Atom")
vlines(2266, 0, 12, colors='r', linestyles='--',alpha=0.6)
vlines(2095, 0, 15.5, colors='r', linestyles='--',alpha=0.6,label=r"$\Delta f=171\,$Hz")
plt.ylabel(r"$P\,/\,$will. Einheit")
plt.xlabel(r"$f\,/\,$Hz")
plt.legend(loc="best")
plt.savefig("plots/Hatom/zustandsaufspaltung_9.pdf")
plt.show()
plt.close()

plt.title("f-Aufspaltung")
plt.plot([3,6,9],[63,112,171],"o")
plt.plot([3,6,9],[63,112,171],"-")
plt.ylabel(r"$\Delta f\,/\,$Hz")
plt.xlabel(r"Ringdicke $d\,/\,$mm")
plt.savefig("plots/Hatom/faufspaltung.pdf")
plt.close()

#Wasserstoffmol
f=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_15mm.dat",skip_header=1,usecols=(0))
D=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_15mm.dat",skip_header=1,usecols=(1))
plt.title("H2-Mol. 180° 15mm")
plt.plot(f,D,label="H2-Mol.")
vlines(2297, 0, 5, colors='r', linestyles='--',alpha=0.6,label="2297 Hz")
vlines(2304, 0, 10.5, colors='r', linestyles='--',alpha=0.6,label="2304 Hz")
vlines(2416, 0, 6, colors='r', linestyles='--',alpha=0.6,label="2416 Hz")
plt.ylabel(r"$P\,/\,$will. Einheit")
plt.xlabel(r"$f\,/\,$Hz")
plt.legend(loc="best")
plt.savefig("plots/Hatom/zustandsaufspaltung_mol15.pdf")
#plt.show()
plt.close()

#resonanzfrequenzen plot
f5=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_5mm.dat",skip_header=1,usecols=(0))
D5=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_5mm.dat",skip_header=1,usecols=(1))
f10=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_10mm.dat",skip_header=1,usecols=(0))
D10=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_10mm.dat",skip_header=1,usecols=(1))
f15=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_15mm.dat",skip_header=1,usecols=(0))
D15=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_15mm.dat",skip_header=1,usecols=(1))
f20=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_20mm.dat",skip_header=1,usecols=(0))
D20=np.genfromtxt("data/Wasserstoffmolekuel/Spektrum_Ring_20mm.dat",skip_header=1,usecols=(1))
#plt.plot(f5,D5,label="5mm")
plt.plot(f10,D10,label="10mm")
plt.plot(f15,D15,label="15mm")
plt.plot(f20,D20,label="20mm")
plt.legend()
#plt.show()
plt.close()

d=[10,15,20]
R1=[2295,2295,2295]#1. Resonanz
R2=[2303,2304,2305]#2. Resonanz
R3=[2374,2417,2456]#2. Resonanz

plt.title("Resonanzfrequenz")
plt.plot(d,R1,"xg",label="1. Resonanz")
plt.plot(d,R1,"--g")
plt.plot(d,R2,"xk",label="2. Resonanz")
plt.plot(d,R2,"--k")
plt.plot(d,R3,"xb",label="3. Resonanz")
plt.plot(d,R3,"--b")
plt.ylabel(r"f$\,/\,$Hz")
plt.xlabel(r"Blendendurchmesser$\,/\,$mm")
plt.legend(loc="best")
plt.savefig("plots/Hatom/res_blende.pdf")
#plt.show()
plt.close()