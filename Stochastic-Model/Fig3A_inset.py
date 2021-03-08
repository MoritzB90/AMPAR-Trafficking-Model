# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 14:08:22 2020

@author: Moritz
"""

"""Fig 3A inset

This script reproduces the plots seen in Fig 3A inset of "The biophysical basis 
underlying the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from ampartrafficking import stochastic_model as sm
from ampartrafficking import rate_model as rm
import seaborn as sns
sns.set_style("ticks")
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('small')

from scipy.integrate import odeint

import sys
sys.path.append('../Mean-Field-Model')

#%%

col1=sns.color_palette("colorblind", 10)[0]
col2=sns.color_palette("colorblind", 10)[4]

Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")
    
#%%

SaveFig=0

duration=20000#Duration in s
Nr_Trials=40

kBU=0.1
kout=0.018
kin=0.02
kUB0=0.0005048
kBU=0.1
kexo0=0.0018
kendo=0.002058
A_spine_basal=0.898

beta=1#0#
alpha=16#0#

N=9
D_0=10

tStim=duration/2
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,5,5,60])
S_exo=13
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])


#%%


ID_basal=1

    
A_spine=A_spine_basal

UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0) 
dt=sm.calcTimeStep(UFP,A_spine,kUB.current_value*5,alpha,kBU,kout+kendo,kin*D_0+kexo.current_value*5*S_exo)

B_Tr=[]
U_Tr=[]

for Trial in range(0,Nr_Trials):
    
    if Trial%10==0:
        print('Trial:',Trial)
        print('P:',N**2, 'U:', UFP)
        
    
    PSD=np.zeros((N,N))
    U=UFP
    
    Time=[]
    B_t=[]
    U_t=[]
    
    for t in np.arange(0,duration+dt,dt):
        
        if t>tStim:
            kexo.update(t-tStim)
            kUB.update(t-tStim)
        
        NN=sm.nearestNeighbours(PSD)
        
        Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
        Mub=sm.kUBcoop(kUB.current_value*U/A_spine, NN, PSD, alpha)*dt

        PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
            
        pout=(kout*U/A_spine+kendo*U/A_spine)*dt
        pin=(kin*D_0+kexo.current_value*S_exo)*dt
        U=sm.update_mobilePool(U,pin,pout,dBoff,dBon)
        
        Time.append(t)
        B_t.append(np.sum(PSD==ID_basal))
        U_t.append(U)
        
    B_Tr.append(B_t)
    U_Tr.append(U_t)
            

#%%

sLTP=0
Cooperativity=1#0#

kin=0.02*D_0
kin_RE=0.1
kout_RE=0.000615
Vspine0=0.08

kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,5,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
Vspine=rm.Parameter(Vspine0)
Vspine.timecourse(rm.DV,[Vspine0,True])

t=np.arange(0,120*60)

   
initList = [[10,23,13]]
pList=[81]


fig=plt.figure(figsize=(3,2), dpi=150)

sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, label='E-LTP', legend=False, zorder=0)

for i,P,Init in zip(range(len(pList)),pList,initList):


    Model=rm.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
    
    
    plt.plot(np.array(Time)/60-duration/2/60,np.mean(B_Tr, axis=0)/np.mean(np.mean(B_Tr, axis=0)[10000:17000])*100, color=col2, label='Stochastic model')
    plt.plot(t/60,solve.T[1]/23*100, color=col1, linewidth=2, label='Rate model')
    
    plt.xlabel('Time (min)')
    plt.ylabel('Bound AMPARs $B$ (%)')
    lgd=plt.legend(bbox_to_anchor=(0.3,0.45,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
    plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
    plt.xlim(-3,130)
    plt.ylim(40)
    sns.despine()
    plt.tight_layout()
    if SaveFig==1:
        print('Save')
        fig.savefig('Figures\\Fig3A_inset.png', bbox_inches="tight", dpi=400)
        fig.savefig('Figures\\Fig3A_inset.svg', bbox_inches="tight", dpi=400)

#%%


# if __name__ == "__main__":
#     main()
    
    
