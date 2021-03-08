# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 08:18:23 2020

@author: Moritz
"""

"""Fig 6A and 7E

This script reproduces the plots seen in Fig 6A and 7E of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `numpy`, `matplotlib` and `seaborn` are installed within the Python 
environment you are running this script in.
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
plt.rcParams['svg.fonttype'] = 'none'
from scipy import stats

#%%

# def main():

SaveFig=0

A_spine=0.898
    
duration=7200*4 #Duration in s
kout=0.01796
kin=0.02
kexo0=0.0018
kendo=0.002058
S_exo=13

Nr_Trials=2#10

beta=1#0
alpha=16#0
kUB0=0.0005048#0.0036
kBU=0.1

D_List=np.array([5,7,8.5,10,12,15,20,25,30])
N_List=[3,5,8,10,12,15]


#%%

B_N=[]
U_N=[]

ID_basal=1

for N in N_List:
    
    B_U=[]
    U_U=[]

    for D_0 in D_List:
        
        UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0)
        dt=sm.calcTimeStep(UFP,A_spine,kUB0,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo)
        
        
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
                
                
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
                Mub=sm.kUBcoop(kUB0*U/A_spine, NN, PSD, alpha)*dt

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

plt.figure()
plt.imshow(PSD)
plt.colorbar()

#%%

BMean_N=[]
BStd_N=[]
BSEM_N=[]
UMean_N=[]
UStd_N=[]
USEM_N=[]
for i,N in enumerate(N_List):
    BMean_U=[]
    BStd_U=[]
    BSEM_U=[]
    UMean_U=[]
    UStd_U=[]
    USEM_U=[]
    for j,D_0 in enumerate(D_List):
        
        
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
        if Nr_Trials>3: #Calculate the standard error across mean values derived from different Trials
            BSEM_U.append(stats.sem(np.mean(np.array(B_N)[i,j,0::,int(len(Time)/2)::], axis=1), axis=0))
        UMean_U.append(UMean)
        UStd_U.append(UStd)
        if Nr_Trials>3:
            USEM_U.append(stats.sem(np.mean(np.array(U_N)[i,j,0::,int(len(Time)/2)::], axis=1), axis=0))

    BMean_N.append(BMean_U)
    BStd_N.append(BStd_U)
    BSEM_N.append(BSEM_U)
    UMean_N.append(UMean_U)
    UStd_N.append(UStd_U)
    USEM_N.append(USEM_U)
    
BMean_N=np.array(BMean_N)
BStd_N=np.array(BStd_N)
BSEM_N=np.array(BSEM_N)
UMean_N=np.array(UMean_N)
UStd_N=np.array(UStd_N)
USEM_N=np.array(USEM_N)
    
#%%

Cooperativity=1

BFP_N=[]
UFP_N=[]

for N in N_List:
    
    P=N**2
    BFP_U=[]
    UFP_U=[]

    for D_0 in np.linspace(0,max(D_List),40):
    
        UFP=rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine)
        BFP=rm.BFP_(UFP,kUB0,kBU,P,A_spine,Cooperativity)
        
        BFP_U.append(BFP)
        UFP_U.append(UFP)

    BFP_N.append(BFP_U)
    UFP_N.append(UFP_U)

BFP_N=np.array(BFP_N)
UFP_N=np.array(UFP_N)

#%%

plt.figure()
plt.plot(np.arange(0,duration+0.5,0.5)/60,np.array(U_t)/A_spine, label='Trace from stoch. sim.')
plt.axhline(np.mean(np.mean(np.array(U_Tr)[::,10000::]/A_spine)), color='r', linewidth=6, label='Mean Stochastic Model')
plt.axhline(UFP_U[-1]/A_spine, color='g', linewidth=4, label='Mean Rate Model')
plt.axhline((kexo0*S_exo+kin*30)/(kendo+kout), color='k', linewidth=2, label='Mean from fixed point eq.')
plt.legend()
#plt.ylim(0,60)
plt.xlabel('Time (min)')
plt.ylabel('$U/A_{spine}$')


#%%

col=sns.color_palette('colorblind')
ul=np.array([0,3,5,7])
fig=plt.figure(figsize=(3.5,2), dpi=150)

for c,i in enumerate([1,3,5]):
    plt.errorbar(x=UMean_N[i,ul]/A_spine, y=BMean_N[i,ul], xerr=UStd_N[i,ul]/A_spine, yerr=BStd_N[i,ul], 
                 color=np.array(col[c])*0.75, linewidth=1.5, fmt='o', markersize=4, label='P={:1.0f}'
                 .format(N_List[i]**2),zorder=0)
    # plt.plot(UMean_N[i,ul]/A_spine, BMean_N[i,ul], color=np.array(col[c])*0.75, linewidth=0, marker='o', ms=5, zorder=0)          
    
    plt.plot(UFP_N[i]/A_spine,BFP_N[i], color=np.array(col[c])*1, linestyle='--', linewidth=2, alpha=0.75, zorder=1)

plt.legend(prop=fontLgd)
plt.xlim(5,30)
plt.xlabel('Mobile AMPAR conc. $U/A_{spine} \, (\#/\mu m^2)$')
plt.ylabel('Bound AMPARs $B^*$ (#)')
sns.despine()
fig.tight_layout()
    

#%%

col=sns.color_palette('colorblind')
ul=np.array([0,3,5,7])
fig=plt.figure(figsize=(3.5,2), dpi=150)

for c,i in enumerate([1,2,5]):
    plt.errorbar(x=UMean_N[i,ul]/A_spine, y=BMean_N[i,ul]/N_List[i]**2*100, xerr=UStd_N[i,ul]/A_spine, 
                 yerr=(BStd_N[i,ul]/N_List[i]**2*100), color=np.array(col[c])*0.75, linewidth=1, fmt='o', markersize=4, 
                 capthick=1, capsize=2, label='P={:1.0f}'.format(N_List[i]**2),zorder=0)         
    # plt.plot(UMean_N[i,ul]/A_spine, BMean_N[i,ul]/N_List[i]**2*100, color=np.array(col[c])*0.75, marker='o', ms=4, linewidth=0)
    
    plt.plot(UFP_N[i]/A_spine,BFP_N[i]/N_List[i]**2*100, color=np.array(col[c]), linestyle='--', linewidth=2, zorder=1)

plt.legend(prop=fontLgd)
plt.xlim(5,31.7)
plt.ylim(0,100)
# plt.axhline(0.772*100)
plt.xlabel('Mobile AMPAR conc. $U/A_{spine} \, (\#/\mu m^2)$')
plt.ylabel('Bound AMPARs \n $B^*/P$ (%)')
sns.despine()
fig.tight_layout()


#%%

col=sns.color_palette('colorblind')

ul=np.array(range(len(D_List)))
fig=plt.figure(figsize=(3.5,2), dpi=150)

for c,i in enumerate([1,2,5]):
    # plt.errorbar(x=UMean_N[i,ul]/A_spine, y=BMean_N[i,ul]/N_List[i]**2*100, xerr=USEM_N[i,ul]/A_spine, 
                 # yerr=(BSEM_N[i,ul]/N_List[i]**2*100), color=np.array(col[c])*0.75, linewidth=1, fmt='o', markersize=4, 
                 # capthick=1, capsize=2,zorder=0)        
    plt.plot(UMean_N[i,ul]/A_spine, BMean_N[i,ul]/N_List[i]**2*100, color=np.array(col[c])*0.75, marker='o', ms=4, linewidth=0.5, label='$P={:1.0f}$'.format(N_List[i]**2))
    
    plt.plot(UFP_N[i]/A_spine,BFP_N[i]/N_List[i]**2*100, color=np.array(col[c]), linestyle='--', linewidth=2, alpha=0.75, zorder=0)

plt.legend(prop=fontLgd)
plt.xlim(5,31.7)
plt.ylim(0,100)
# plt.axhline(0.772*100)
plt.xlabel('Mobile AMPAR conc. $U/A_{spine} \, (\#/\mu m^2)$')
plt.ylabel('Bound AMPARs \n $B^*/P$ (%)')
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig6A.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig6A.svg', bbox_inches="tight", dpi=400)


#%%



Cooperativity=1

BFP_N=[]
UFP_N=[]

for N in N_List:
    
    P=N**2
    BFP_U=[]
    UFP_U=[]

    for D_0 in [5.0,10.0,15.0,25.0]:
    
        UFP=rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine)
        BFP=rm.BFP_(UFP,kUB0,kBU,P,A_spine,Cooperativity)
        
        BFP_U.append(BFP)
        UFP_U.append(UFP)

    BFP_N.append(BFP_U)
    UFP_N.append(UFP_U)

BFP_N=np.array(BFP_N)
UFP_N=np.array(UFP_N)
col=sns.color_palette('colorblind')

#%%


fig=plt.figure(figsize=(3,2), dpi=150)

i=-1
for j in np.arange(0,len(D_List)):
    if D_List[j] in [5,10,15,25]:
        i+=1
        if i==0:
            plt.errorbar(x=np.array(N_List)**2, y=BMean_N[::,j], yerr=BStd_N[::,j], 
                         color=np.array(col[i])*0.5, linewidth=1.5, fmt='o', markersize=4, 
                         label=r'$\langle U/A_{spine} \rangle =$'+'${0:1.1f}\pm{1:1.1f}$'
                         .format(UMean_N[0][j]/A_spine, UStd_N[0][j]/A_spine)+'$\, \#/\mu m^2$',zorder=1)
        else:
            plt.errorbar(x=np.array(N_List)**2, y=BMean_N[::,j], yerr=BStd_N[::,j], 
                         color=np.array(col[i])*0.5, linewidth=1.5, fmt='o', markersize=4, label='$={0:1.1f}\pm{1:1.1f}$'
                         .format(UMean_N[0][j]/A_spine, UStd_N[0][j]/A_spine)+'$\, \#/\mu m^2$',zorder=1)
    print(D_List[j])


for i in range(4): 
    plt.plot(np.array(N_List)**2,BFP_N.T[i], color=np.array(col[i])*1, linestyle='--', alpha=0.75, linewidth=2, zorder=0)
    
plt.xlim(0,230)
plt.ylim(-10,225)
plt.legend(ncol=1,prop=fontLgd)
plt.xlabel('Slots $P$')
plt.ylabel('Bound AMPARs \n $B^*$ (#)')
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7E.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7E.svg', bbox_inches="tight", dpi=400)
    
#%%

