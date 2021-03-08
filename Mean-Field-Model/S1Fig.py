# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 18:40:26 2020

@author: Moritz
"""

"""S1 Fig

This script reproduces the plots seen in S1 Fig of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation". Note that depending on your machine this script can take 3-4 hours to run.

This script requires that `numpy`,`pandas`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in. 
Note that depending on your machine this script can take 2-3 hours to run.
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
fontLgd.set_size('x-small')
import pandas as pd
import time

# np.random.seed(1000)

start = time.time()
print('Note that depending on your machine this script can take 2-3 hours to run.')

#%%

Run=False#True#
SaveFig=0
saveData=0
loadData=1

Pars=[[0,0,'Basic Model'],[1,0,'sLTP Model'],[0,1,'Cooperative Model']]


if Run==True:
    kBU_max=1
    TUB_max=500
    kexo_max=0.01
    aexo_max=20
    Texo_max=4000
    V0=0.08
    A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)
    
    t=np.arange(0,130*60)
        
    Trials=2000#200#
    
    Init_List=[[5,10,13],[10,20,13],[20,40,13]]#
    dwellTime_List=[20,50,100]
    
    
    dwellTime_container=[]
    kexo0_container=[]
    kexoA_container=[]
    Texo_container=[]
    hue=[]
    mode=[]
    Goodness=[]
    GoodnessMean=[]

    for sLTP,Cooperativity,Label in Pars:
        for Init in Init_List:
            for dwellTime in dwellTime_List:
                
                P=3.5*Init[1]
                
                kout=A0/dwellTime
                kin=kout/A0*Init[0]
                
                if Cooperativity==0:
                    aUB_max=50
                else:
                    aUB_max=20
                    
                print('$dwell Time$=',dwellTime,' | $BFP$=',Init[1],' | $P$=',P,' | Label=',Label)
                
    
                B_Tr,B_ne_Tr, U_Tr, kUB_Tr,kexo_Tr, Aspine_Tr, kexo0_Tr,aexo_Tr,Texo_Tr=ps.parameterSampling(t,Init,sLTP,Cooperativity,P,kin,kout,Trials, kBU_max,aUB_max,TUB_max,kexo_max,aexo_max,Texo_max)
                
                GoodMatch_index, GoodMatch_Value=ps.Matching(B_Tr,B_ne_Tr,Trials,t,Init[1])
    
    
                dwellTime_container.extend([dwellTime]*len(GoodMatch_index))
                kexo0_container.extend(kexo0_Tr[GoodMatch_index])
                kexoA_container.extend(kexo0_Tr[GoodMatch_index]*(1+aexo_Tr[GoodMatch_index]))
                Texo_container.extend(Texo_Tr[GoodMatch_index]/60)
                hue.extend([str(Init[1])+','+str(np.round(Init[0]/A0,1))]*len(GoodMatch_index))
                mode.extend([Label]*len(GoodMatch_index))
                Goodness.extend(GoodMatch_Value)
                GoodnessMean.append(np.mean(GoodMatch_Value))
    
    
    
    d_kexo0={r'$\tau_{dwell}$ (s)': dwellTime_container, '$k_{exo}^{0}$ (Events/s)': kexo0_container, '$B^*,U^*/A_{spine}$':hue, '-':mode}
    df_kexo0 = pd.DataFrame(data=d_kexo0)
    
    d_kexoA={r'$\tau_{dwell}$ (s)': dwellTime_container, '$k_{exo}^{A}$ (Events/s)': kexoA_container, '$B^*,U^*/A_{spine}$':hue, '-':mode}
    df_kexoA = pd.DataFrame(data=d_kexoA)
    
    d_Texo={r'$\tau_{dwell}$ (s)': dwellTime_container, r'$\tau_{exo}$ (min)': Texo_container, '$B^*,U^*/A_{spine}$':hue, '-':mode}
    df_Texo = pd.DataFrame(data=d_Texo)
    #dError={r'$\tau_{dwell}$ in s': x, r'$Error in #$': np.array(Goodness)/BFP_V, '$B^*,U^*$':hue, 'Model':mode}
    #df_Error = pd.DataFrame(data=dError)
    
    
    if saveData==1:
        df_kexo0.to_csv('Data_PS\\df_kexo0.csv', index=False)
        df_kexoA.to_csv('Data_PS\\df_kexoA.csv', index=False)
        df_Texo.to_csv('Data_PS\\df_Texo.csv', index=False)
        #dfError.to_csv('Data_PS\\df_Error.csv', index=False)
        
    
    


#%%


if loadData==1:
    df_kexo0=pd.read_csv('Data_PS\\df_kexo0.csv')
    df_kexoA=pd.read_csv('Data_PS\\df_kexoA.csv')
    df_Texo=pd.read_csv('Data_PS\\df_Texo.csv')


colorpalette='muted'

fig=plt.figure(figsize=(10,8), dpi=250)
plot0=sns.catplot(x=r'$\tau_{dwell}$ (s)', y='$k_{exo}^{0}$ (Events/s)', hue='$B^*,U^*/A_{spine}$', col="-", data=df_kexo0, ci='sd', palette=colorpalette, dodge=True, linewidth=1, edgecolor='gray', kind="bar", height=2, aspect=1)
for i in range(len(Pars)):
    plot0.axes[0][i].fill_between([-0.6,2.6], [0.0,0.0], [0.0018,0.0018], color=sns.color_palette(colorpalette)[3], alpha=0.5, zorder=0)
    plot0.axes[0][i].set_xlim(-0.5,2.5)
    plot0.axes[0][i].set_ylim(0,0.015)
fig.tight_layout()
if SaveFig==1:
    plot0.savefig('Figures\\kexo0_Tdwell_SD.png', bbox_inches="tight", dpi=400)
    plot0.savefig('Figures\\kexo0_Tdwell_SD.svg', bbox_inches="tight", dpi=400)

fig=plt.figure(figsize=(10,8), dpi=250)
plotA=sns.catplot(x=r'$\tau_{dwell}$ (s)', y='$k_{exo}^{A}$ (Events/s)', hue='$B^*,U^*/A_{spine}$', col="-", data=df_kexoA, ci='sd', palette=colorpalette, dodge=True, linewidth=1, edgecolor='gray', kind="bar", height=2, aspect=1)
for i in range(len(Pars)):
    plotA.axes[0][i].fill_between([-0.6,2.6], [0.0,0.0], [0.011,0.011], color=sns.color_palette(colorpalette)[3], alpha=0.5, zorder=0)
    plotA.axes[0][i].set_xlim(-0.5,2.5)
    plotA.axes[0][i].set_ylim(0)
fig.tight_layout()
if SaveFig==1:
    plotA.savefig('Figures\\kexoA_Tdwell_SD.png', bbox_inches="tight", dpi=400)
    plotA.savefig('Figures\\kexoA_Tdwell_SD.svg', bbox_inches="tight", dpi=400)

fig=plt.figure(figsize=(10,8), dpi=250)
plotT=sns.catplot(x=r'$\tau_{dwell}$ (s)', y=r'$\tau_{exo}$ (min)', hue='$B^*,U^*/A_{spine}$', col="-", data=df_Texo, ci='sd', palette=colorpalette, dodge=True, linewidth=1, edgecolor='gray', kind="bar", height=2, aspect=1)     
for i in range(len(Pars)):
    plotT.axes[0][i].fill_between([-0.6,2.6], [0.0,0.0], [2,2], color=sns.color_palette(colorpalette)[3], alpha=0.5, zorder=0)
    plotT.axes[0][i].set_xlim(-0.5,2.5)
    plotT.axes[0][i].set_ylim(0)
fig.tight_layout()
if SaveFig==1:
    plotT.savefig('Figures\\Texo_Tdwell_SD.png', bbox_inches="tight", dpi=400)
    plotT.savefig('Figures\\Texo_Tdwell_SD.svg', bbox_inches="tight", dpi=400)
    

end = time.time()
print('Runtime:', end - start)

