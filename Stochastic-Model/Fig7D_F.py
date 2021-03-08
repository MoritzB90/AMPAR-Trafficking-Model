# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:09:36 2021

@author: Moritz
"""


"""Fig 3B

This script reproduces the plots seen in Fig 7D and F of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `pandas`,`numpy`,`matplotlib` and `seaborn` 
are installed within the Python environment you are running this script in.
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
fontLgd.set_size('xx-small')
plt.rcParams['svg.fonttype'] = 'none'

def D_FP(kendo,kexo0,S_exo,kout,kin,Conc):
    return ((kendo+kout)*Conc-kexo0*S_exo)/kin

#%%

# def main():

SaveFig=0#1

A_spine=0.898
    
duration=20000#Duration in s
kBU=0.1
kout=0.018
kin=0.02
kexo0=0.0018
kendo=0.002058
S_exo=13

Nr_Trials=1

beta=1#0#
alpha=16#0#
kUB=0.0005#*7
kBU=0.1


Conc_List=[15/A_spine,30/A_spine]
N_List=[7,10]

#%%

PSD_SnapShot=[]
U_SnapShot=[]

B_N=[]
U_N=[]

ID_basal=1

for N in N_List:
    
    B_U=[]
    U_U=[]

    for Conc in Conc_List:
        
        D_0=D_FP(kendo,kexo0,S_exo,kout,kin,Conc)
                
        UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0)
        dt=sm.calcTimeStep(UFP,A_spine,kUB,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo)
        
        # print('U_FP:', UFP)
        # print('Conc.:', Conc)
        B_Tr=[]
        U_Tr=[]

        
        for Trial in range(0,Nr_Trials):
            
            # if Trial%10==0:
            print('Trial:',Trial)
            print('P:',N**2, 'U:', UFP)
                
            
            PSD=np.zeros((N,N))

            Time=[]
            B_t=[]
            U_t=[]

            U=UFP

            for t in np.arange(0,duration+dt,dt):
                    
                # if Conc==15/A_spine and N==7 and Trial==0 and abs(t-duration/2)<dt/2:
                if abs(t-duration/2)<dt/2:
                    PSD_SnapShot.append(PSD)                
                    U_SnapShot.append(U)
                    
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
                Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD, alpha)*dt

                PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
                
                pout=(kout*U/A_spine+kendo*U/A_spine)*dt
                pin=(kin*D_0+kexo0*S_exo)*dt
                U=sm.update_mobilePool(U,pin,pout,dBoff,dBon)
                
                if t%0.5==0:
                    Time.append(t)
                    B_t.append(np.sum(PSD==1))
                    U_t.append(U)
                
            B_Tr.append(B_t)
            U_Tr.append(U_t)

        B_U.append(B_Tr)
        U_U.append(U_Tr)

    B_N.append(B_U)
    U_N.append(U_U)
    
    
#%%

BMean_N=[]
BStd_N=[]
UMean_N=[]
UStd_N=[]
for i,N in enumerate(N_List):
    BMean_U=[]
    BStd_U=[]
    UMean_U=[]
    UStd_U=[]
    for j,Conc in enumerate(Conc_List):
        
        
        BMean=np.mean(np.mean(np.array(B_N)[i,j,0::,int(len(Time)/2)::], axis=1))
        BStd=np.mean(np.std(np.array(B_N)[i,j,0::,int(len(Time)/2)::], axis=1))
        UMean=np.mean(np.mean(np.array(U_N)[i,j,0::,int(len(Time)/2)::], axis=1))
        UStd=np.mean(np.std(np.array(U_N)[i,j,0::,int(len(Time)/2)::], axis=1))
        
        print(BMean)
        print(BStd)
        print(UMean)
        print(UStd)
        
        BMean_U.append(BMean)
        BStd_U.append(BStd)

        UMean_U.append(UMean)
        UStd_U.append(UStd)

    BMean_N.append(BMean_U)
    BStd_N.append(BStd_U)
    UMean_N.append(UMean_U)
    UStd_N.append(UStd_U)
    
BMean_N=np.array(BMean_N)
BStd_N=np.array(BStd_N)
UMean_N=np.array(UMean_N)
UStd_N=np.array(UStd_N)

#%%

i=0
plt.figure(figsize=(2,2), dpi=150)
plt.imshow(PSD_SnapShot[i])
plt.locator_params(axis='x', nbins=2)
plt.locator_params(axis='y', nbins=2)
plt.tight_layout()
plt.title('$U/A_{spine}=$'+'${:1.1f}$, '.format(U_SnapShot[i]/A_spine)+'$B={:1.0f}$'.format(np.sum(PSD_SnapShot[i])), font=fontLgd)
plt.xlabel('Position')
plt.ylabel('Position')
if SaveFig==1:
    plt.savefig('Figures\\7F.png', bbox_inches="tight", dpi=400)
    plt.savefig('Figures\\7F.svg', bbox_inches="tight", dpi=400)
    
#%%

col=sns.color_palette('colorblind')
                         
DT=int(duration/10)

fig=plt.figure(figsize=(3,2), dpi=150)

dn=10
k=-1
for i,n in enumerate(N_List):
    for j,conc in enumerate(Conc_List):
        # if conc in [15/A_spine,30/A_spine] and n in [7,10]:
        k+=1
        if k==0:
            plt.plot(np.array(Time)[0::dn]/60,B_N[i][j][0][0::dn],c=col[k], 
                      label=r'$\langle U/A_{spine} \rangle =$'+'${0:1.1f}\pm{1:1.1f}$'
                      .format(UMean_N[i][j]/A_spine, UStd_N[i][j]/A_spine))
                            
        else:
            plt.plot(np.array(Time)[0::dn]/60,B_N[i][j][0][0::dn],c=col[k], label='$={0:1.1f}\pm{1:1.1f}$'
                      .format(UMean_N[i][j]/A_spine, UStd_N[i][j]/A_spine))
                            

# plt.axhline(0, c='k', linewidth=1)
plt.text(0.1,0.3,'$k_{UB}='+'{:1.4f}$ '.format(kUB)+'$\mu m^2/(\# s)^{-1}$', transform=fig.axes[0].transAxes)
plt.text(0.2,0.4,'$P=49$', transform=fig.axes[0].transAxes)
plt.text(0.2,0.5,'$P=100$', transform=fig.axes[0].transAxes)
plt.legend(bbox_to_anchor=(-0.1,0.3,0,0), loc="upper left", mode="normal", ncol=2, prop=fontLgd)
plt.ylim(-35,100)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs \n $B$ (#)')
plt.xlim(0,60)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7D.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7D.svg', bbox_inches="tight", dpi=400)