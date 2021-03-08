# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 17:01:58 2020

@author: Moritz
"""

"""Figs 5 B-D

This script reproduces the plots seen in Figs 5 B-D of "The biophysical basis underlying 
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
import seaborn as sns
sns.set_style("ticks")
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('xx-small')

#%%

# def main():
    
SaveFig=0

#Import experimental Data for comparison with the model. Taken from Penn et al 2017 and Barco et al. 2002 (see Paper for reference details):
Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")

#%%

#Color palettes used in plots:
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[1:4]

#%%

#Set parameter values and initial conditions:
    
kin=0.2
kout=0.018
kexo0=0.0018
kendo=0.002058
kBU=0.1
kUB0=0.0005048
kin_RE=0.1
kout_RE=0.000615

Vspine_basal=0.08
Aspine_basal=4*np.pi*(3*Vspine_basal/(4*np.pi))**(2/3)
P_basal=70
Sexo_basal=13

kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,5,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
    
t=np.arange(0,120*60)

#%%

#Fig 5 B:

Vspine0=0.08
Aspine0=4*np.pi*(3*Vspine0/(4*np.pi))**(2/3)
Vspine=rm.Parameter(Vspine0)
Vspine.timecourse(rm.DV,[Vspine0,True])

pList=[40,100,200]

sLTP=0
Cooperativity=1
    
fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, label='E-LTP', zorder=0)
    
for i,P in zip(range(len(pList)),pList):

    SexoFP=rm.SexoFP_(kin_RE,kout_RE,Vspine0)
    UFP=rm.UFP_(kexo0,SexoFP,kendo,kin,kout,Aspine0)
    BFP=rm.BFP_(UFP,kUB0,kBU,P,Aspine0,Cooperativity)
    Init=[UFP,BFP,SexoFP]
    
    Model=rm.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
    
    
    plt.plot(t/60,solve.T[1]/Init[1]*100, color=col1[i], linewidth=2, label='$P={0:1.1f},\, B^*{1:1.1f}$'.format(P,Init[1]))
    
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
lgd=plt.legend(bbox_to_anchor=(0.15,0.6,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlim(-3,130)
plt.ylim(80)
sns.despine()
plt.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig5B.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig5B.svg', bbox_inches="tight", dpi=400)
        
#%%

#Fig 5 C:

vList=np.linspace(0.001,0.501,41, endpoint=True)#[0.039,0.08,0.376]

sLTP=1
Cooperativity=1
    
fig=plt.figure(figsize=(3.5,1.25), dpi=150)

DB1=[]
DB20=[]
i=-1
for j, Vspine0 in zip(range(len(vList)),vList):
    
    #Calculate the rate at which AMPARs from endosomes exit the spine such that Sexo=Sexo_basal=13 independent of the initial spine volume:
    kout_RE=kin_RE*Vspine0/Sexo_basal
    
    Aspine0=4*np.pi*(3*Vspine0/(4*np.pi))**(2/3)
    P=Aspine0/Aspine_basal*P_basal
    
    SexoFP=rm.SexoFP_(kin_RE,kout_RE,Vspine0)
    UFP=rm.UFP_(kexo0,SexoFP,kendo,kin,kout,Aspine0)
    BFP=rm.BFP_(UFP,kUB0,kBU,P,Aspine0,Cooperativity)
    Init=[UFP,BFP,SexoFP]
    
    
    Vspine=rm.Parameter(Vspine0)
    Vspine.timecourse(rm.DV,[Vspine0,True])
    
    Model=rm.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
    
    ii_20=np.where(t/60==20)[0][0]
    DB20.append(solve.T[1][ii_20]-Init[1])
    DB1.append(max(solve.T[1]-Init[1]))
    
    if np.round(Vspine0,3) in [0.039,0.089,0.376]:
        i+=1
        plt.plot(t/60,solve.T[1]-Init[1], color=col1[i], label='$V_{spine}^{0}=$'+'{:1.3f} $\mu m^3$'.format(Vspine0))
    
plt.xlabel('Time (min)')
plt.ylabel('$\Delta B$ (#)')
lgd=plt.legend(bbox_to_anchor=(0.57,-0.2,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(0, color='k', linestyle='--', linewidth=0.5)
sns.despine()
plt.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig5C_top.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig5C_top.svg', bbox_inches="tight", dpi=400)    

#%%

col2=sns.cubehelix_palette(3, start=.5, rot=-.75, dark=0.2, light=0.7)[0]
col3=sns.cubehelix_palette(3, start=.5, rot=-.75, dark=0.2, light=0.7)[1]

fig=plt.figure(figsize=(3.5,1.25),dpi=150)
sns.lineplot(vList, DB20, color=col2, label='after 20 min.', linewidth=2, legend=False)
sns.lineplot(vList, DB1, color=col3, label='after 1 min.', linewidth=2, legend=False)
plt.xlabel(r'$V_{spine}^{0}$ in $\mu m^3$')
plt.ylabel(r'$\Delta B$ in #')
lgd=plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlim(0,0.5)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig5C_bottom.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig5C_bottom.svg', bbox_inches="tight", dpi=400)
    
#%%

#Fig 5 D:

vList=np.linspace(0.001,0.501,41, endpoint=True)#[0.039,0.08,0.376]
DB1=[]
DB20=[]

kout_RE=0.000615
sLTP=1
Cooperativity=1
    
fig=plt.figure(figsize=(3.5,1.25), dpi=150)

i=-1
for j, Vspine0 in zip(range(len(vList)),vList):
    
    Aspine0=4*np.pi*(3*Vspine0/(4*np.pi))**(2/3)
    P=Aspine0/Aspine_basal*P_basal
    
    SexoFP=rm.SexoFP_(kin_RE,kout_RE,Vspine0)
    UFP=rm.UFP_(kexo0,SexoFP,kendo,kin,kout,Aspine0)
    BFP=rm.BFP_(UFP,kUB0,kBU,P,Aspine0,Cooperativity)
    Init=[UFP,BFP,SexoFP]
    
    
    Vspine=rm.Parameter(Vspine0)
    Vspine.timecourse(rm.DV,[Vspine0,True])
    
    #Calculate fixed point of the spine volume dependent exocytosis event size Sexo and update the intial state:
    Sexo_0=kin_RE/kout_RE*Vspine0
    Init[2]=Sexo_0
    
    Model=rm.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
   
    ii_20=np.where(t/60==20)[0][0]
    DB20.append(solve.T[1][ii_20]-Init[1])
    DB1.append(max(solve.T[1]-Init[1]))
    
    if np.round(Vspine0,3) in [0.039,0.089,0.376]:
        i+=1
        plt.plot(t/60,solve.T[1]-Init[1], color=col1[i], label='$V_{spine}^{0}=$'+'{:1.3f} $\mu m^3$'.format(Vspine0))
    
plt.xlabel('Time (min)')
plt.ylabel('$\Delta B$ (#)')
lgd=plt.legend(bbox_to_anchor=(0.57,-0.2,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(0, color='k', linestyle='--', linewidth=0.5)
sns.despine()
plt.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig5D_top.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig5D_top.svg', bbox_inches="tight", dpi=400)

#%%

col2=sns.cubehelix_palette(3, start=.5, rot=-.75, dark=0.2, light=0.7)[0]
col3=sns.cubehelix_palette(3, start=.5, rot=-.75, dark=0.2, light=0.7)[1]

fig=plt.figure(figsize=(3.5,1.25),dpi=150)
sns.lineplot(vList, DB20, color=col2, label='after 20 min.', linewidth=2, legend=False)
sns.lineplot(vList, DB1, color=col3, label='after 1 min.', linewidth=2, legend=False)
plt.xlabel(r'$V_{spine}^{0}$ in $\mu m^3$')
plt.ylabel(r'$\Delta B$ in #')
lgd=plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlim(0,0.5)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig5D_bottom.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig5D_bottom.svg', bbox_inches="tight", dpi=400)
    
#%%

# if __name__ == "__main__":
#     main()
    
    