# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 20:18:45 2020

@author: Moritz
"""

"""Fig 5A

This script reproduces the plots seen in Fig 5A of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from ampartrafficking import rate_model as rm
import seaborn as sns
sns.set_style("ticks")
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('small')

    
#%%

SaveFig=0

Cooperativity=1

Init=[10,20,13]
Vspine=0.08
Aspine=4*np.pi*(3*Vspine/(4*np.pi))**(2/3)

P=70
kin=0.2
kout=0.018
kexo0=0.0018
kendo=rm.kendo_(Init,kexo0,kin,kout,Aspine) #0.002058
kBU=0.1
kUB0=rm.kUB0_(Init,kBU,P,Aspine,Cooperativity) #0.0005048
kin_RE=0.1
kout_RE=0.000615

Vspine_List=np.linspace(0.0001,0.2501,81, endpoint=True)

BFP_Coop=[]
UFP_List=[]
SexoFP_List=[]
for i,Vspine in enumerate(Vspine_List):
    
    Aspine=4*np.pi*(3*Vspine/(4*np.pi))**(2/3)
    
    SexoFP=rm.SexoFP_(kin_RE,kout_RE,Vspine)
    UFP=rm.UFP_(kexo0,SexoFP,kendo,kin,kout,Aspine)
    
    BFP=rm.BFP_(UFP,kUB0,kBU,P,Aspine,Cooperativity)
    
    BFP_Coop.append(BFP)
    UFP_List.append(UFP)
    SexoFP_List.append(SexoFP)
    

#%%


Cooperativity=0

Init=[10,20,13]
Vspine=0.08
Aspine=4*np.pi*(3*Vspine/(4*np.pi))**(2/3)

P=70
kin=0.2
kout=0.018
kexo0=0.0018
kendo=rm.kendo_(Init,kexo0,kin,kout,Aspine) #0.002058
kBU=0.1
kUB0=rm.kUB0_(Init,kBU,P,Aspine,Cooperativity) #0.00359
kin_RE=0.1
kout_RE=0.000615

Vspine_List=np.linspace(0.0001,0.2501,81, endpoint=True)

BFP_Basic=[]
for i,Vspine in enumerate(Vspine_List):
    
    Aspine=4*np.pi*(3*Vspine/(4*np.pi))**(2/3)
    
    SexoFP=rm.SexoFP_(kin_RE,kout_RE,Vspine)
    UFP=rm.UFP_(kexo0,SexoFP,kendo,kin,kout,Aspine)
    
    BFP=rm.BFP_(UFP,kUB0,kBU,P,Aspine,Cooperativity)
    
    BFP_Basic.append(BFP)
    

#%%

col0=sns.light_palette(sns.color_palette("colorblind", 10)[0],3)[1]
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[2]
col2=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[1]
col=[col0,col1,col2,col0,col1,col2,'k','k']


fig,axes=plt.subplots(figsize=(3.5,2.35), dpi=150)
ax2 = axes.twiny()
axes.plot(SexoFP_List,BFP_Coop, Linewidth=2, label='Coop. Model', color=col[0])
axes.plot(SexoFP_List,BFP_Basic, Linewidth=2, label='Basic Model', color='k')
ax2.plot(Vspine_List,BFP_Coop, Linewidth=2, label='Coop. Model', color=col[0])
ax2.plot(Vspine_List,BFP_Basic, Linewidth=2, label='Basic Model', color='k')
lgd = plt.legend(bbox_to_anchor=(0.02,0.80,0,0), loc="lower left", mode="normal", borderaxespad=0, ncol=6, prop=fontLgd)
plt.axhline(0, c='k', Linewidth=0.5)
ax2.axvline(0.08, c='k', Linewidth=0.5)
ax2.axvline(0.16, c='k', Linewidth=0.5)
# plt.xlabel('$k_{out}$')
ax2.set_xlabel('Spine volume $V_{spine}$ ($\mu m^3$)')
axes.set_xlabel('Exocytosis event size $S_{exo}$ (#)')
axes.set_ylabel('Bound AMPARs $B^*$ (#)')
ax2.set_xlim(0,0.25)
axes.set_xlim(0,40)
plt.ylim(10,40)
sns.despine(top=False)
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig5A.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig5A.svg', bbox_inches="tight", dpi=400)
