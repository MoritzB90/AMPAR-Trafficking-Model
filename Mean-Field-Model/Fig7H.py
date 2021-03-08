# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 08:18:23 2020

@author: Moritz
"""


"""Fig 7H

This script reproduces the plots seen in Fig 7H of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.optimize`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("ticks")
from matplotlib.font_manager import FontProperties 
fontLgd = FontProperties()
fontLgd.set_size('x-small')
plt.rcParams['svg.fonttype'] = 'none'
from scipy.optimize import curve_fit

#%%

Matsuzaki2004_DVlong_x=np.loadtxt("Data\\Matsuzaki2004_DVlong_x.csv", delimiter=",") #This data from Matsuzaki et al. comprise the relationship between the initial size of a spine and the enlargment 70min. after sLTP induction.
Matsuzaki2004_DVlong_y=np.loadtxt("Data\\Matsuzaki2004_DVlong_y.csv", delimiter=",")

Yang2008_DVlong_x=np.loadtxt("Data\\Yang2008_DVlong_x.csv", delimiter=",") #This data from Yang et al. comprise the relationship between the initial size of a spine and the enlargment 45min. after sLTP induction.
Yang2008_DVlong_y=np.loadtxt("Data\\Yang2008_DVlong_y.csv", delimiter=",")
Yang2008_DVlong_y=(Yang2008_DVlong_y-Yang2008_DVlong_x)/Yang2008_DVlong_x*100

Tanaka2008_DVlong_x=np.loadtxt("Data\\Tanaka2008_DVlong_x.csv", delimiter=",") #This data from Tanaka et al. comprise the relationship between the initial size of a spine and the enlargment 40min. after sLTP induction.
Tanaka2008_DVlong_y=np.loadtxt("Data\\Tanaka2008_DVlong_y.csv", delimiter=",")

#%%

SaveFig=0

#The experimental data is fit by an exponential function.
def func(x,a,b,c):
    return a*np.exp(-x/b)+c


poptY, pcovY = curve_fit(func, Yang2008_DVlong_x, Yang2008_DVlong_y, p0=(400,0.05,0), maxfev = 1000)
print('Parameter (Y):', poptY)
print('Standard deviation error (Y):', np.sqrt(np.diag(pcovY)))
print('Relative error (Y):', np.sqrt(np.diag(pcovY))/poptY)
poptM, pcovM = curve_fit(func, Matsuzaki2004_DVlong_x, Matsuzaki2004_DVlong_y, p0=(400,0.05,0), maxfev = 1000)
print('Parameter (M):', poptM)
print('Standard deviation error (M):', np.sqrt(np.diag(pcovM)))
print('Relative error (M):', np.sqrt(np.diag(pcovM))/poptM)
poptT, pcovT = curve_fit(func, Tanaka2008_DVlong_x, Tanaka2008_DVlong_y, p0=(400,0.05,0), maxfev = 1000)
print('Parameter (T):', poptT)
print('Standard deviation error (T):', np.sqrt(np.diag(pcovT)))
print('Relative error (T):', np.sqrt(np.diag(pcovT))/poptT)

#%%

col=sns.color_palette("rocket")[0::]
Vinitial_List=np.linspace(0.00,0.6,1000)


plt.figure(figsize=(3.5,2), dpi=150)
plt.plot(Yang2008_DVlong_x,Yang2008_DVlong_y, marker='d',markersize=4, color=col[5], markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, linestyle='-', linewidth=0, label='Yang et al. 2008',zorder=0)
plt.xlim(0,0.6)
plt.ylim(-75,550)
plt.axhline(0, linestyle='--',linewidth=1, color='k')
plt.plot(Vinitial_List,func(Vinitial_List,*poptY), color=col[5], linestyle='--', linewidth=2,zorder=4)
plt.legend()
plt.xlabel('Initial spine volume $V_{spine}^{0} \, (\mu m^3)$')
plt.ylabel('Long-term spine vol. \n change $\Delta V_{spine}^{long} \, (\%)$')

plt.figure(figsize=(3.5,2), dpi=150)
plt.plot(Matsuzaki2004_DVlong_x,Matsuzaki2004_DVlong_y, marker='o',markersize=4, color=col[1], markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, linestyle='-', linewidth=0, label='Matsuzaki et al. 2004',zorder=0)
plt.xlim(0,0.6)
plt.ylim(-75,550)
plt.axhline(0, linestyle='--',linewidth=1, color='k')
plt.plot(Vinitial_List,func(Vinitial_List,*poptM), color=col[1], linestyle='--', linewidth=2,zorder=5)
plt.legend()
plt.xlabel('Initial spine volume $V_{spine}^{0} \, (\mu m^3)$')
plt.ylabel('Long-term spine vol. \n change $\Delta V_{spine}^{long} \, (\%)$')

plt.figure(figsize=(3.5,2), dpi=150)
plt.plot(Tanaka2008_DVlong_x,Tanaka2008_DVlong_y, marker='s',markersize=4, color=col[3], markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, linestyle='-', linewidth=0, label='Tanaka et al. 2008',zorder=0)
plt.xlim(0,0.6)
plt.ylim(-75,550)
plt.axhline(0, linestyle='--',linewidth=1, color='k')
plt.plot(Vinitial_List,func(Vinitial_List,*poptT), color=col[3], linestyle='--', linewidth=2,zorder=6)
plt.legend()
plt.xlabel('Initial spine volume $V_{spine}^{0} \, (\mu m^3)$')
plt.ylabel('Long-term spine vol. \n change $\Delta V_{spine}^{long} \, (\%)$')

#%%

popt=(poptY+poptM+poptT)/3


#%%

fig,axes=plt.subplots(figsize=(3.5,2), dpi=150)
plt.plot(Yang2008_DVlong_x,Yang2008_DVlong_y, marker='d',markersize=4, color=col[5], markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, linestyle='-', linewidth=0, label='Yang et al. 2008',zorder=0)
plt.plot(Matsuzaki2004_DVlong_x,Matsuzaki2004_DVlong_y, marker='o',markersize=4, color=col[1], markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, linestyle='-', linewidth=0, label='Matsuzaki et al. 2004',zorder=1)
plt.plot(Tanaka2008_DVlong_x,Tanaka2008_DVlong_y, marker='s',markersize=4, color=col[3], markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, linestyle='-', linewidth=0, label='Tanaka et al. 2008',zorder=2)


plt.plot(Vinitial_List,func(Vinitial_List,*popt), label=r'$\langle a \rangle={:1.0f}$'.format(popt[0])+r', $\langle b \rangle={:1.3f}$'.format(popt[1])+r', $\langle c \rangle={:1.1f}$'.format(popt[2]), color='k', linewidth=2,zorder=7)
plt.plot(Vinitial_List,func(Vinitial_List,*poptY), color=col[5], linestyle='--', linewidth=2,zorder=4)
plt.plot(Vinitial_List,func(Vinitial_List,*poptT), color=col[3], linestyle='--', linewidth=2,zorder=5)
plt.plot(Vinitial_List,func(Vinitial_List,*poptM), color=col[1], linestyle='--', linewidth=2,zorder=6)


plt.xlabel('Initial spine volume $V_{spine}^{0} \, (\mu m^3)$')
plt.ylabel('Long-term spine vol. \n change $\Delta V_{spine}^{long} \, (\%)$')
lgd=plt.legend(bbox_to_anchor=(0.1,0.3,0,0), loc="lower left", mode="normal", borderaxespad=0, ncol=1, prop=fontLgd).set_zorder(8)
plt.axhline(0, linestyle='--',linewidth=1, color='k')
plt.xlim(0,0.3)
plt.ylim(-75,550)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7H.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7H.svg', bbox_inches="tight", dpi=400)

