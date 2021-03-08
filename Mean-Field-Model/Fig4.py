# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:49:23 2020

@author: Moritz
"""


"""Fig 4

This script reproduces the plots seen in Fig 4 of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation". 
Note that depending on your machine this script can take 15-30 min to run.

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from ampartrafficking import parameter_sampling as ps
import seaborn as sns
sns.set_style("ticks")
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('xx-small')

# np.random.seed(1000)

print('Note that depending on your machine this script can take 15-30 min to run.')

#%%

# def main():

SaveFig=0
    
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[3]
col2=sns.dark_palette(sns.color_palette("colorblind", 10)[1],4)[3]

#%%

#Import experimental Data for comparison with the model. Taken from Penn et al 2017 and Barco et al. 2002 (see Paper for reference details):

Penn2017_x=np.loadtxt("Data\\Penn2017_x.csv", delimiter=",")
Penn2017_y=np.loadtxt("Data\\Penn2017_y.csv", delimiter=",")
Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")

Lee2009_x=np.loadtxt("Data\\Lee2009_x.csv", delimiter=",") 
Lee2009_y=np.loadtxt("Data\\Lee2009_y.csv", delimiter=",")

Patterson2010_x=np.loadtxt("Data\\Patterson2010_x.csv", delimiter=",")
Patterson2010_y=np.loadtxt("Data\\Patterson2010_y.csv", delimiter=",")

#%%


figB=plt.figure('B',figsize=(3.5,2), dpi=150)
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[1:4]
col2=sns.dark_palette(sns.color_palette("colorblind", 10)[1],4)[1:4]

figkexo=plt.figure('kexo',figsize=(3.5,2), dpi=150)
col_kexo=sns.dark_palette(sns.color_palette("colorblind", 10)[2],4)[1:4]

figkUB,axes=plt.subplots(figsize=(3.5,2), dpi=150, num='kUB')
col_kUB=sns.dark_palette(sns.color_palette("colorblind", 10)[3],4)[1:4]
ax2 = axes.twinx()
axes.zorder=1
axes.set_facecolor('xkcd:salmon')
axes.patch.set_alpha(0.0)
ax2.zorder=0

figU=plt.figure('U',figsize=(3.5,2), dpi=150)
col_U=sns.dark_palette(sns.color_palette("colorblind", 10)[8],4)[1:4]


plt.figure('B')
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, label='E-LTP', legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, label='LTP, $k_{exo}=0$', legend=False, zorder=0)

plt.figure('kexo')
sns.lineplot(Patterson2010_x/60,Patterson2010_y/60, color='k', marker="o", ms=5, label='Patterson et al. 2010', legend=False, zorder=0)

plt.figure('kUB')
sns.lineplot(Lee2009_x/60,Lee2009_y/max(Lee2009_y), ax=ax2, marker='o', color='k', label='Lee et al. 2009', legend=False, zorder=0)




#%%

#Basic AMPAR trafficking

Trials=200#2000# 10# 

P=70
kin=0.2
kout=0.018
kBU_max=1
TUB_max=500
kexo_max=0.01
aexo_max=20
Texo_max=4000

Init = [10,20,13]
t=np.arange(0,130*60)

#%%

sLTP=0
Cooperativity=0
aUB_max=50

B_Tr,B_ne_Tr, U_Tr, kUB_Tr,kexo_Tr, Aspine_Tr, kexo0_Tr,aexo_Tr,Texo_Tr=ps.parameterSampling(t,Init,sLTP,Cooperativity,P,kin,kout,Trials, kBU_max,aUB_max,TUB_max,kexo_max,aexo_max,Texo_max)

GoodMatch_index, GoodMatch_Value=ps.Matching(B_Tr,B_ne_Tr,Trials,t,Init[1])

plt.figure('B')
X=np.reshape(([t]*len(GoodMatch_index)), B_Tr[GoodMatch_index].size)/60
Y=np.reshape(B_Tr[GoodMatch_index]/Init[1]*100, B_Tr[GoodMatch_index].size)
X_nE=np.reshape(([t]*len(GoodMatch_index)), B_ne_Tr[GoodMatch_index].size)/60
Y_nE=np.reshape(B_ne_Tr[GoodMatch_index]/Init[1]*100, B_ne_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, linewidth=2, ci='sd', color=col1[0], label='Basic Model', legend=False)
sns.lineplot(X_nE,Y_nE, linewidth=2, ci='sd', color=col2[0], label='Basic Model'+',$k_{exo}=0$', legend=False)

plt.figure('kexo')
X=np.reshape(([t]*len(GoodMatch_index)), kexo_Tr[GoodMatch_index].size)/60
Y=np.reshape(kexo_Tr[GoodMatch_index], kexo_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', color=col_kexo[0], label='Basic Model', legend=False, zorder=0)
plt.axhline(np.mean(kexo_Tr[GoodMatch_index][::,0]), color='k', linestyle='--', linewidth=0.5)

plt.figure('kUB')
X=np.reshape(([t]*len(GoodMatch_index)), kUB_Tr[GoodMatch_index].size)/60
Y=np.reshape(kUB_Tr[GoodMatch_index]/kUB_Tr[GoodMatch_index][:,0][:,np.newaxis, ], kUB_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', ax=axes, color=col_kUB[0], label='Basic Model', legend=False, zorder=0)
axes.axhline(np.mean(kUB_Tr[GoodMatch_index][::,0]), color='k', linestyle='--', linewidth=0.5)

plt.figure('U')
V0=0.08
A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)
X=np.reshape(([t]*len(GoodMatch_index)), U_Tr[GoodMatch_index].size)/60
Y=np.reshape((U_Tr[GoodMatch_index]/Aspine_Tr)/(Init[0]/A0)*100, U_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', linewidth=2, color=col_U[0], label='Basic Model', legend=False)

#%%

sLTP=1
Cooperativity=0
aUB_max=50

B_Tr,B_ne_Tr, U_Tr, kUB_Tr,kexo_Tr, Aspine_Tr, kexo0_Tr,aexo_Tr,Texo_Tr=ps.parameterSampling(t,Init,sLTP,Cooperativity,P,kin,kout,Trials, kBU_max,aUB_max,TUB_max,kexo_max,aexo_max,Texo_max)

GoodMatch_index, GoodMatch_Value=ps.Matching(B_Tr,B_ne_Tr,Trials,t,Init[1])

plt.figure('B')
X=np.reshape(([t]*len(GoodMatch_index)), B_Tr[GoodMatch_index].size)/60
Y=np.reshape(B_Tr[GoodMatch_index]/Init[1]*100, B_Tr[GoodMatch_index].size)
X_nE=np.reshape(([t]*len(GoodMatch_index)), B_ne_Tr[GoodMatch_index].size)/60
Y_nE=np.reshape(B_ne_Tr[GoodMatch_index]/Init[1]*100, B_ne_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, linewidth=2, ci='sd', color=col1[1], label='sLTP Model', legend=False)
sns.lineplot(X_nE,Y_nE, linewidth=2, ci='sd', color=col2[1], label='sLTP Model'+',$k_{exo}=0$', legend=False)

plt.figure('kexo')
X=np.reshape(([t]*len(GoodMatch_index)), kexo_Tr[GoodMatch_index].size)/60
Y=np.reshape(kexo_Tr[GoodMatch_index], kexo_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', color=col_kexo[1], label='sLTP Model', legend=False, zorder=0)
plt.axhline(np.mean(kexo_Tr[GoodMatch_index][::,0]), color='k', linestyle='--', linewidth=0.5)

plt.figure('kUB')
X=np.reshape(([t]*len(GoodMatch_index)), kUB_Tr[GoodMatch_index].size)/60
Y=np.reshape(kUB_Tr[GoodMatch_index]/kUB_Tr[GoodMatch_index][:,0][:,np.newaxis, ], kUB_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', ax=axes, color=col_kUB[1], label='sLTP Model', legend=False, zorder=0)
axes.axhline(np.mean(kUB_Tr[GoodMatch_index][::,0]), color='k', linestyle='--', linewidth=0.5)

plt.figure('U')
V0=0.08
A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)
X=np.reshape(([t]*len(GoodMatch_index)), U_Tr[GoodMatch_index].size)/60
Y=np.reshape((U_Tr[GoodMatch_index]/Aspine_Tr)/(Init[0]/A0)*100, U_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', linewidth=2, color=col_U[1], label='sLTP Model', legend=False)

#%%

sLTP=0
Cooperativity=1
aUB_max=20

B_Tr,B_ne_Tr, U_Tr, kUB_Tr,kexo_Tr, Aspine_Tr, kexo0_Tr,aexo_Tr,Texo_Tr=ps.parameterSampling(t,Init,sLTP,Cooperativity,P,kin,kout,Trials, kBU_max,aUB_max,TUB_max,kexo_max,aexo_max,Texo_max)

GoodMatch_index, GoodMatch_Value=ps.Matching(B_Tr,B_ne_Tr,Trials,t,Init[1])

plt.figure('B')
X=np.reshape(([t]*len(GoodMatch_index)), B_Tr[GoodMatch_index].size)/60
Y=np.reshape(B_Tr[GoodMatch_index]/Init[1]*100, B_Tr[GoodMatch_index].size)
X_nE=np.reshape(([t]*len(GoodMatch_index)), B_ne_Tr[GoodMatch_index].size)/60
Y_nE=np.reshape(B_ne_Tr[GoodMatch_index]/Init[1]*100, B_ne_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, linewidth=2, ci='sd', color=col1[2], label='Coop. Model', legend=False)
sns.lineplot(X_nE,Y_nE, linewidth=2, ci='sd', color=col2[2], label='Coop. Model'+',$k_{exo}=0$', legend=False)

plt.figure('kexo')
X=np.reshape(([t]*len(GoodMatch_index)), kexo_Tr[GoodMatch_index].size)/60
Y=np.reshape(kexo_Tr[GoodMatch_index], kexo_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', color=col_kexo[2], label='Coop. Model', legend=False, zorder=0)
plt.axhline(np.mean(kexo_Tr[GoodMatch_index][::,0]), color='k', linestyle='--', linewidth=0.5)

plt.figure('kUB')
X=np.reshape(([t]*len(GoodMatch_index)), kUB_Tr[GoodMatch_index].size)/60
Y=np.reshape(kUB_Tr[GoodMatch_index]/kUB_Tr[GoodMatch_index][:,0][:,np.newaxis, ], kUB_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', ax=axes, color=col_kUB[2], label='Coop. Model', legend=False, zorder=0)
axes.axhline(np.mean(kUB_Tr[GoodMatch_index][::,0]), color='k', linestyle='--', linewidth=0.5)

plt.figure('U')
V0=0.08
A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)
X=np.reshape(([t]*len(GoodMatch_index)), U_Tr[GoodMatch_index].size)/60
Y=np.reshape((U_Tr[GoodMatch_index]/Aspine_Tr)/(Init[0]/A0)*100, U_Tr[GoodMatch_index].size)
sns.lineplot(X, Y, ci='sd', linewidth=2, color=col_U[2], label='Coop. Model', legend=False)

#%%

figB=plt.figure('B')
lgd=plt.legend(bbox_to_anchor=(0.0,0.4,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
figB.tight_layout()
if SaveFig==1:
    figB.savefig('Figures\\Fig4A.png', bbox_inches="tight", dpi=400)
    figB.savefig('Figures\\Fig4A.svg', bbox_inches="tight", dpi=400)

figkexo=plt.figure('kexo')
lgd=plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlabel('Time (min)')
plt.ylabel('Exocytosis event rate \n $k_{exo} \, (Events/s)$')
plt.xlim(-10/60,10)
sns.despine()
figkexo.tight_layout()
if SaveFig==1:
    figkexo.savefig('Figures\\Fig4C.png', bbox_inches="tight", dpi=400)
    figkexo.savefig('Figures\\Fig4C.svg', bbox_inches="tight", dpi=400)
    
figkUB=plt.figure('kUB')
axes.set_xlabel('Time (min)')
axes.set_ylabel('Binding rate $k_{UB}/k_{UB}^0$')
ax2.set_ylabel('CaMKII activity a.u.')
lgd=axes.figure.legend(bbox_to_anchor=(0.65,0.85,0,0), loc="upper right", mode="normal", borderaxespad=0, ncol=1, prop=fontLgd)
plt.xlim(-10/60,4)
# axes.set_xlim(0,4)
# ax2.set_xlim(0,4)
ax2.set_ylim(-0.4,1.5)     
sns.despine(right=False)
figkUB.tight_layout()
if SaveFig==1:
    figkUB.savefig('Figures\\Fig4B.png', bbox_inches="tight", dpi=400)
    figkUB.savefig('Figures\\Fig4B.svg', bbox_inches="tight", dpi=400)
    
figU=plt.figure('U')
lgd=plt.legend(bbox_to_anchor=(0.35,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Mobile AMPAR conc. \n $U/A_{spine} \, (\%)$')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
figU.tight_layout()
if SaveFig==1:
    figU.savefig('Figures\\Fig4D.png', bbox_inches="tight", dpi=400)
    figU.savefig('Figures\\Fig4D.svg', bbox_inches="tight", dpi=400)
    
    
#%%
