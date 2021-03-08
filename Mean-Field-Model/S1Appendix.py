# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 10:12:48 2020

@author: Moritz
"""


"""S2 Fig

This script reproduces the plots seen in Fig 1 in S1 Appendix of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ampartrafficking import rate_model as rm
from ampartrafficking import rate_model_suppl as rms
import seaborn as sns
sns.set_style("ticks")
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('xx-small')


#%%

# def main():

SaveFig=0

#Color palettes used in plots:
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[3]
col2=sns.dark_palette(sns.color_palette("colorblind", 10)[1],4)[3]
col3=sns.dark_palette(sns.color_palette("colorblind", 10)[8],4)[3]

col4=sns.dark_palette(sns.color_palette("colorblind", 10)[2],4)[3]
col5=sns.dark_palette(sns.color_palette("colorblind", 10)[3],4)[3]

col_I=['#EE6666', '#3388BB', '#9988DD',
                 '#EECC55', '#88BB44', '#FFBBBB']
#%%

#Import experimental Data for comparison with the model. Taken from Penn et al 2017 and Barco et al. 2002 (see Paper for reference details):

Penn2017_x=np.loadtxt("Data\\Penn2017_x.csv", delimiter=",")
Penn2017_y=np.loadtxt("Data\\Penn2017_y.csv", delimiter=",")
Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")

#%%

#Set parameter values and initial conditions:

Init2 = [5,10,13]

Init = [10,20,13,0,0,0]


V0=0.08
A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)

P=70
kin=0.2
kin2=0

kout=0.018
kout2=0.018

kexo0=0.0018
kexo02=0.0018#0#

kendo=rm.kendo_(Init,kexo0,kin,kout,A0)
kendo2=kendo

kBU=0.1
kBU2=0.025

t=np.arange(0,130*60)

#%%

#Basic AMPAR trafficking

sLTP=1
Cooperativity=0
kUB0=rms.kUB0_(Init,kBU,P,A0,Cooperativity)
kUB02=kUB0*10

k12=rm.Parameter(0.25)
k12.timecourse(rm.Stim_Resp,[0,1,5,180])
k21=0.0

#E-LTP with exocytosis:
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,30,5,60])
kUB2=rm.Parameter(kUB02)
kUB2.timecourse(rm.Stim_Resp,[1,30,5,60])

kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
kexo2=rm.Parameter(kexo02)
kexo2.timecourse(rm.Stim_Resp,[1,5,25,60])

Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])

kin_RE=0.1
kout_RE=((kin_RE)*V0)/Init[2]

Model=rms.Model_system_2()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kin2,kout,kout2,kexo,kexo2,kendo,kendo2,kUB,kUB2,kBU,kBU2,
                                                Cooperativity,sLTP,kin_RE,kout_RE,k12,k21),full_output=True)

#E-LTP without exocytosis:
kexo=rm.Parameter(0)
kexo.timecourse(rm.Stim_Resp,[1,0,25,60])
kexo2=rm.Parameter(0)
kexo2.timecourse(rm.Stim_Resp,[1,0,25,60])

Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,False])

Model=rms.Model_system_2()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kin2,kout,kout2,kexo,kexo2,kendo,kendo2,kUB,kUB2,kBU,kBU2,
                                                Cooperativity,sLTP,kin_RE,kout_RE,k12,k21),full_output=True)

#%%

fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(x=Barco2002_x,y=Barco2002_y, color='k', marker="s", ms=5, label='E-LTP', legend=False, zorder=0)
sns.lineplot(x=Penn2017_x,y=Penn2017_y, color='k', marker="v", ms=5, label='LTP, $k_{exo}=0$', legend=False, zorder=0)
plt.plot(t/60,(solve.T[1]+solve.T[4])/(Init[1]+Init[4])*100, color=col1, linewidth=2, label='Exoc. factor model')
plt.plot(t/60,(solve2.T[1]+solve2.T[4])/(Init[1]+Init[4])*100, color=col2, linewidth=2, label='Exoc. factor model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.05,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs \n $(B^{I}+B^{II})$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\S2FigA.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\S2FigA.svg', bbox_inches="tight", dpi=400)


#%%

fig=plt.figure(figsize=(3.5,2), dpi=150)
plt.plot(t/60,(solve.T[1]), color=col_I[0], linewidth=2, label='Type I')
plt.plot(t/60,(solve2.T[1]), color=col_I[1], linewidth=2, label='Type I, $k_{exo}=0$')
plt.plot(t/60,(solve.T[4]), color=col_I[2], linewidth=2, label='Type II')
plt.plot(t/60,(solve2.T[4]), color=col_I[4], linewidth=2, label='Type II, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.05,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (#)')
plt.xlim(-3,130)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\S2FigC.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\S2FigC.svg', bbox_inches="tight", dpi=400)
    
#%%

fig=plt.figure(figsize=(3.5,2), dpi=150)
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])
Vspine.update(t)
Aspine=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)
plt.plot(t/60,((solve.T[0]+solve.T[3])/Aspine)/((Init[0]+Init[3])/A0)*100, color=col3, linewidth=2)
# lgd=plt.legend(bbox_to_anchor=(0.05,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Mobile AMPAR conc. \n $(U^{I}+U^{II})/A_{spine}$ (%)')
plt.xlim(-3,130)
plt.ylim(80)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\S2FigB.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\S2FigB.svg', bbox_inches="tight", dpi=400)
    
#%%

fig, ax1 = plt.subplots(figsize=(3.5,2), dpi=150)
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])
Vspine.update(t)
Aspine=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)

ax1.plot(t/60,(solve.T[0]/Aspine), color=col_I[0], linewidth=2, label='Type I')
ax1.set_ylabel('Mobile AMPAR conc. \n type I ($\#/\mu m^2$)')
ax1.set_ylim(0,15)
ax2.plot(t/60,(solve.T[3]/Aspine), color=col_I[2], linewidth=2, label='Type II')

fig.legend(bbox_to_anchor=(0.35,0.75,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)

ax1.set_xlabel('Time (min)')
ax2.set_ylabel('Mobile AMPAR conc. \n type II ($\#/\mu m^2$)')
ax1.set_xlim(-3,130)
ax2.set_ylim(0,1)
sns.despine(right=False)
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\S2FigD.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\S2FigD.svg', bbox_inches="tight", dpi=400)
    
#%%

fig=plt.figure(figsize=(3.5,2), dpi=150)
plt.plot(t/60, solve.T[2], color='k', linewidth=2, label='Type I')
plt.plot(t/60, solve.T[5], color='grey', linewidth=2, linestyle='--', label='Type II')
lgd=plt.legend(bbox_to_anchor=(0.05,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.xlabel('Time (min)')
plt.ylabel('Endosomal AMPARs \n $S_{exo}$ (#)')
plt.xlim(-3,130)
plt.ylim(0)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\S2FigE.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\S2FigE.svg', bbox_inches="tight", dpi=400)
    
