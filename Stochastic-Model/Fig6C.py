# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 08:51:58 2021

@author: Moritz
"""



"""Fig 3B

This script reproduces the plots seen in Fig 6C of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`matplotlib` and `seaborn` are installed within the 
Python environment you are running this script in.
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
import copy
plt.rcParams['svg.fonttype'] = 'none'

#%%

def D_FP(kendo,kexo0,S_exo,kout,kin,Conc):
    return ((kendo+kout)*Conc-kexo0*S_exo)/kin

#%%

# def main():

SaveFig=0

A_spine_basal=0.898
    
duration=20000#7200# #Duration in s
kBU=0.1
kout=0.018
kin=0.02
kexo0=0.0018
kendo=0.002058
S_exo=13

Nr_Trials=10#20

beta=1#0#
alpha=16#0#
kUB=0.000505#*7
kBU=0.1


Conc_List=[10/A_spine_basal,20/A_spine_basal,30/A_spine_basal]
N_List=[7,10]


#%%


ID_basal=1

B_Mean=np.zeros((len(N_List),len(Conc_List)))
U_Mean=np.zeros((len(N_List),len(Conc_List)))
B_Std=np.zeros((len(N_List),len(Conc_List)))
U_Std=np.zeros((len(N_List),len(Conc_List)))

DurationOccupied_trace=[] #Duration a slot was in an occupied state
NrOccupied_trace=[] #Number of transitions of a slot between occupied and free state
for i,N in enumerate(N_List):
    
    A_spine=A_spine_basal
    
    NN_occ_0=[]
    NN_free_0=[]

    for j,Conc in enumerate(Conc_List):
        
        D_0=D_FP(kendo,kexo0,S_exo,kout,kin,Conc)
                
        UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0)
        dt=sm.calcTimeStep(UFP,A_spine,kUB,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo)
        
        
        DurationOccupied_trace_Temp=np.zeros((N,N))
        NrOccupied_trace_Temp=np.zeros((N,N))
        
        B_Trial=[]
        U_Trial=[]
        
        for Trial in range(0,Nr_Trials):
            
            if Trial%10==0:
                print('Trial:',Trial)
                print('P:',N**2, 'U:', UFP)
                
            
            PSD=np.zeros((N,N))
            U=UFP

            started_collecting=0
            for t in np.arange(0,duration+dt,dt):
                
                
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
                Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD, alpha)*dt
                
                PSD_past=copy.copy(PSD)
                PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
                
                pout=(kout*U/A_spine+kendo*U/A_spine)*dt
                pin=(kin*D_0+kexo0*S_exo)*dt
                U=sm.update_mobilePool(U,pin,pout,dBoff,dBon)
                
                if t>duration/2 and t%0.5==0:
                    if started_collecting==0:
                        NrOccupied_trace_Temp+=PSD_past #When we start collecting the data, the No. of times a slot has been in an occupied state is increased by the current state of the PSD. 
                        #Note that this introduces a bias as we do not no for how long an element has been in an occupied state already. However this bias gets smaller the longer the trials.
                        started_collecting=1
                    
                    dPSD=PSD-PSD_past*2
                    DurationOccupied_trace_Temp[dPSD==-1]+=dt #If a slot is occupied and was also occupied before, dS==-1 is True and the dwell time increases by dt
                    NrOccupied_trace_Temp[dPSD==1]+=1 # if a slot that was not occupied is now occupied, dS==1 is True and the number of times that a slot was in the occupied state is increased by 1.

                    B_Trial.append(np.sum(PSD))
                    U_Trial.append(U)

        B_Mean[i,j]=np.mean(B_Trial)
        U_Mean[i,j]=np.mean(U_Trial)
        B_Std[i,j]=np.std(B_Trial)
        U_Std[i,j]=np.std(U_Trial)  
        DurationOccupied_trace.append(DurationOccupied_trace_Temp)
        NrOccupied_trace.append(NrOccupied_trace_Temp)
        


#%%

plt.figure()
plt.imshow(PSD)
plt.colorbar()

#%%
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(3*3.3,2*2), dpi=150)

axes[0][0].yaxis.set_major_formatter(FormatStrFormatter('%.1d'))

e=-1
e2=-1
i0=0
j0=0
for i,n in enumerate(N_List):
    for j,Conc in enumerate(Conc_List):
        
        e+=1
        
        if n in [7,10] and Conc in [10/A_spine_basal,20/A_spine_basal,30/A_spine_basal]:
            e2+=1
            
            DwellTime=DurationOccupied_trace[e]/NrOccupied_trace[e]
            im=axes[i0,j0].imshow(DwellTime)
            if n==7:
                axes[i0,j0].set_title(r'$\langle U/A_{spine} \rangle=$'+'{0:1.1f} $\pm$ {1:1.1f}'
                                      .format(np.array(U_Mean)[i][j]/A_spine, np.array(U_Std)[i][j]/A_spine)+r' $\frac{\#}{\mu m^2}$'+
                                      '\n'+r'$\langle B \rangle=$'+'{0:1.1f} $\pm$ {1:1.1f}'
                                      .format(np.array(B_Mean)[i][j],np.array(B_Std)[i][j]), fontsize=10)
            # else:
                # axes[i0,j0].set_title('$\overline{B}=$'+'{:1.1f}'.format(np.array(B_Mean)[i][j]) )
            cbar=fig.colorbar(im, ax=axes[i0,j0])
            cbar.set_label('PSD dwell time \n'+r'$\langle \tau_{dwell}^{PSD} \rangle$ (s)', rotation=270, labelpad=27, y=0.55)
            
            tick_locator = ticker.MaxNLocator(nbins=3)
            cbar.locator = tick_locator
            cbar.update_ticks()
            
            axes[i0,j0].locator_params(axis='y', nbins=2)
            
            if j0<2:
                j0+=1
                
            else:
                j0=0
            if e2<2:
                i0=0
            else:
                i0=1

axes[0,0].set_xlabel('Position')
axes[0,0].set_ylabel('Position')
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig6C.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig6C.svg', bbox_inches="tight", dpi=400)
    

#%%


# if __name__ == "__main__":
#     main()
    
    