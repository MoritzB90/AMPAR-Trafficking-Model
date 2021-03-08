# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 13:01:29 2020

@author: Moritz
"""

"""Fig 3B

This script reproduces the plots seen in Fig 3B of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation", which shows the dependence between 
normalized mean cooperative binding (unbinding) rate and the number of bound receptors B 
as well as the PSD size P.

This script requires that `pandas`,`numpy`,`scipy.optimize`,`matplotlib` and `seaborn` 
are installed within the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from ampartrafficking import stochastic_model as sm
from ampartrafficking import rate_model as rm
import pandas as pd
from scipy.optimize import curve_fit
import seaborn as sns
sns.set_style("ticks")
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('small')
plt.rcParams['svg.fonttype'] = 'none'



#%%

# def main():

SaveFig=0#1

A_spine_basal=0.898
    
duration=10000#Duration in s
kBU=0.1
kout=0.018
kin=0.02
kexo0=0.0018
kendo=0.002058
S_exo=13

Nr_Trials=2

beta=1#0#
alpha=16#0#
kUB=0.0005#*7
kBU=0.1


D_List=[3,7,9,10,11,12,15,20,30,60]
N_List=[5,7,10,12]


#%%

Boff_N=[]
kBUcoop_N=[]
Bon_N=[]
kUBcoop_N=[]

NN_occ=[]
NN_free=[]


ID_basal=1

for N in N_List:
    
    A_spine=A_spine_basal
    
    Boff_0=[]
    kBUcoop_0=[]
    Bon_0=[]
    kUBcoop_0=[]
    
    NN_occ_0=[]
    NN_free_0=[]

    for D_0 in D_List:
        
        UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0)
        dt=sm.calcTimeStep(UFP,A_spine,kUB,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo)
        
        for Trial in range(0,Nr_Trials):
            
            if Trial%10==0:
                print('Trial:',Trial)
                print('P:',N**2, 'U:', UFP)
                
            
            PSD=np.zeros((N,N))
            U=UFP

            for t in np.arange(0,duration+dt,dt):
                
                
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
                Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD, alpha)*dt
                
                if t>duration/2 and t%0.5==0:
                    if PSD[PSD==ID_basal].size>0: #if statement necessary, otherwise 
                        kBUcoop_0.append(np.mean(Mbu[PSD==ID_basal])/dt)
                        Boff_0.append(np.sum(PSD))
                        NN_occ_0.append(np.mean(NN[PSD==ID_basal]))
                    if PSD[PSD==0].size>0 and U>0:
                        kUBcoop_0.append(np.mean(Mub[PSD==0])/(U/A_spine)/dt)
                        Bon_0.append(np.sum(PSD))
                        NN_free_0.append(np.mean(NN[PSD==0]))
                        
                        
                PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
                
                pout=(kout*U/A_spine+kendo*U/A_spine)*dt
                pin=(kin*D_0+kexo0*S_exo)*dt
                U=sm.update_mobilePool(U,pin,pout,dBoff,dBon)
                


    Boff_N.append(Boff_0)
    kBUcoop_N.append(kBUcoop_0)
    Bon_N.append(Bon_0)
    kUBcoop_N.append(kUBcoop_0)
    
    NN_occ.append(NN_occ_0)
    NN_free.append(NN_free_0)


#%%

plt.figure()
plt.imshow(PSD)
plt.colorbar()

#%%

#Fig 3 B, bottom


def kBUcoop_Fit(B, a, b):
    return a/(b+B)-0.5

lambda_plt=[]
beta_plt=[]
dn=1

col=sns.cubehelix_palette(len(N_List), start=.5, rot=-.75, dark=0.2, light=0.7)

fig=plt.figure(figsize=(4,3), dpi=150)
for Boff,kBUc,N,i in zip(Boff_N,kBUcoop_N,N_List,range(len(N_List))):
    
    d={'value1': Boff[0::dn], 'value2': np.array(kBUc[0::dn])/kBU}
    df = pd.DataFrame(data=d)
    
    popt_kBU, pcov_kBU = curve_fit(kBUcoop_Fit, Boff[0::dn], np.array(kBUc[0::dn])/kBU, p0=(1,0.9), bounds=([0,0],[np.inf,np.inf]), maxfev=10000)
     
    sns.lineplot(x="value1", y="value2", data=df, color=col[i], linewidth=3, label='$P={:1d}$'.format(N**2), ci=None)
    plt.plot(np.arange(1,max(Boff[0::dn])), kBUcoop_Fit(np.arange(1,max(Boff[0::dn])), *popt_kBU), color='k', linestyle='--')

    lambda_plt.append(popt_kBU[0])
    beta_plt.append(popt_kBU[1])

plt.text(0.5,0.5,'$k_{UB}$='+'{:1.4f} '.format(kUB)+'$(\# s)^{-1}$', transform=fig.axes[0].transAxes)
sns.despine()
plt.xlabel('Bound AMPARs $B$ (#)')
plt.ylabel(r'Cooperative unbinding'+'\n'+r' $\langle \hat{k}_{BU}^{coop} \rangle/k_{BU}$')
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlim(0,150)
plt.ylim(0,1.01)
Ratio=fig.axes[0].get_xlim()[1]/fig.axes[0].get_ylim()[1]
fig.axes[0].set_aspect(aspect=Ratio/2, adjustable='box')
#fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig3B_bottom.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig3B_bottom.svg', bbox_inches="tight", dpi=400)

#%%

#Fig 7 B

Marker_List=['o','s', 'v', 'p', 'd', '*','h', 'x', 'X']

def lambda_P(x, a, c):
    return a*x+c

popt_lambda, pcov_lambda = curve_fit(lambda_P, np.array(N_List)**2, lambda_plt, p0=(3000,10), maxfev=800)

col=sns.cubehelix_palette(1, start=.5, rot=-.75, dark=0.2, light=0.5)  

fig=plt.figure(figsize=(3,2), dpi=150)
plt.plot(np.array(N_List)**2, lambda_plt, color='k', marker=Marker_List[0])
plt.plot(np.linspace(min(np.array(N_List)**2),max(np.array(N_List)**2),100), lambda_P(np.linspace(min(np.array(N_List)**2),max(np.array(N_List)**2),100), *popt_lambda), color=col[0], label=r'$m_{\lambda}$'+'$={:1.2f}$,'.format(popt_lambda[0])+'\n'+r'$c_{\lambda}$'+'$={:1.2f}$'.format(popt_lambda[1]))
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=1)
plt.xlabel('Slots $P$')
plt.ylabel(r'$\lambda$')
sns.despine()
plt.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7B.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7B.svg', bbox_inches="tight", dpi=400)

#%%

#Fig 7 C

def beta_P(x, a, c):
    return a*x+c


popt_beta, pcov_beta = curve_fit(beta_P, np.array(N_List)**2, beta_plt, p0=(3000,10), maxfev=800)
    
fig=plt.figure(figsize=(3,2), dpi=150)
plt.plot(np.array(N_List)**2, beta_plt, color='k', marker=Marker_List[0])
plt.plot(np.linspace(min(np.array(N_List)**2),max(np.array(N_List)**2),100), beta_P(np.linspace(min(np.array(N_List)**2),max(np.array(N_List)**2),100), *popt_beta), color=col[0], label=r'$m_{\beta}$'+'$={:1.2f}$,'.format(popt_beta[0])+'\n'+r'$c_{\beta}$'+'$={:1.2f}$'.format(popt_beta[1]))
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=1)
plt.xlabel('Slots $P$')
plt.ylabel(r'$\beta$')
sns.despine()
plt.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7C.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7C.svg', bbox_inches="tight", dpi=400)
    
#%%

#Fig 3 B, top

def kUBcoop_Fit(B, a):
    return (a*B**0.8+1)

m_plt=[]

col=sns.cubehelix_palette(len(N_List), start=.5, rot=0.0, dark=0.2, light=0.7)

fig=plt.figure(figsize=(4,3), dpi=150)
for Bon,kUBc,N,i in zip(Bon_N,kUBcoop_N,N_List,range(len(N_List))):
    
    d={'value1': Bon[0::dn], 'value2': np.array(kUBc[0::dn])/kUB}
    df = pd.DataFrame(data=d)
    
    popt_kUB, pcov3_kUB = curve_fit(kUBcoop_Fit, Bon[0::dn], np.array(kUBc[0::dn])/kUB, p0=(0.1), maxfev=10000)
       
    sns.lineplot(x="value1", y="value2", data=df, color=col[i], linewidth=3, label='$P={:1d}$'.format(N**2), ci=None)
    plt.plot(np.arange(1,max(Bon[0::dn])), kUBcoop_Fit(np.arange(1,max(Bon[0::dn])), *popt_kUB), color='k', linestyle='--')
    
    m_plt.append(popt_kUB[0])

plt.text(0.5,0.35,'$k_{UB}$='+'{:1.4f} '.format(kUB)+'$(\# s)^{-1}$', transform=fig.axes[0].transAxes)
plt.ylim(0)
sns.despine()
plt.xlabel('Bound AMPARs $B$ (#)')
plt.ylabel(r'Cooperative binding'+'\n'+r'$\langle \hat{k}_{UB}^{coop} \rangle/k_{UB}$')
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlim(0,150)
plt.ylim(0,11)
Ratio=fig.axes[0].get_xlim()[1]/fig.axes[0].get_ylim()[1]
fig.axes[0].set_aspect(aspect=Ratio/2, adjustable='box')
plt.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig3B_top.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig3B_top.svg', bbox_inches="tight", dpi=400)
    
#%%

#Fig 7 A

def m_P(x, a, b):
    return a/(b+x)

popt_m, pcov_m = curve_fit(m_P, np.array(N_List)**2, m_plt, p0=(3000,100), maxfev=800)

col=sns.cubehelix_palette(1, start=.7, rot=0.0, dark=0.2, light=0.5)
  
fig=plt.figure(figsize=(3,2), dpi=150)
plt.plot(np.array(N_List)**2, m_plt, color='k', marker=Marker_List[0])
plt.plot(np.linspace(min(np.array(N_List)**2),max(np.array(N_List)**2),100), m_P(np.linspace(min(np.array(N_List)**2),max(np.array(N_List)**2),100), *popt_m), color=col[0], label='$a_m={:1.2f}$,'.format(popt_m[0])+'\n'+'$b_m={:1.2f}$'.format(popt_m[1]))
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=1)
plt.xlabel('Slots $P$')
plt.ylabel('m')
sns.despine()
plt.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7A.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7A.svg', bbox_inches="tight", dpi=400)


#%%


col=sns.cubehelix_palette(len(N_List), start=.5, rot=-.75, dark=0.2, light=0.7)

fig=plt.figure(figsize=(4,3), dpi=150)
for Boff,NNocc,N,i in zip(Boff_N,NN_occ,N_List,range(len(N_List))):
    
    d={'value1': Boff[0::dn], 'value2': np.array(NNocc[0::dn])/8}
    df = pd.DataFrame(data=d)
    
    sns.lineplot(x="value1", y="value2", data=df, color=col[i], linewidth=3, label='$P={:1d}$'.format(N**2), ci=None)


plt.text(0.5,0.5,'$k_{UB}$='+'{:1.4f} '.format(kUB)+'$(\# s)^{-1}$', transform=fig.axes[0].transAxes)
sns.despine()
plt.xlabel('Bound AMPARs $B$ (#)')
plt.ylabel(r'$\chi$ (occpudied slots)')
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlim(0,150)
plt.ylim(0)
Ratio=fig.axes[0].get_xlim()[1]/fig.axes[0].get_ylim()[1]
fig.axes[0].set_aspect(aspect=Ratio/2, adjustable='box')
#fig.tight_layout()


#%%

col=sns.cubehelix_palette(len(N_List), start=.5, rot=0.0, dark=0.2, light=0.7)

fig=plt.figure(figsize=(4,3), dpi=150)
for Bon,NNfree,N,i in zip(Bon_N,NN_free,N_List,range(len(N_List))):
    
    d={'value1': Bon[0::dn], 'value2': np.array(NNfree[0::dn])/8}
    df = pd.DataFrame(data=d)
    
    sns.lineplot(x="value1", y="value2", data=df, color=col[i], linewidth=3, label='$P={:1d}$'.format(N**2), ci=None)


plt.text(0.5,0.5,'$k_{UB}$='+'{:1.4f} '.format(kUB)+'$(\# s)^{-1}$', transform=fig.axes[0].transAxes)
sns.despine()
plt.xlabel('Bound AMPARs $B$ (#)')
plt.ylabel(r'$\chi$ (free slots)')
lgd = plt.legend(loc="best", mode="normal", borderaxespad=0, ncol=2, prop=fontLgd)
plt.xlim(0,150)
plt.ylim(0)
Ratio=fig.axes[0].get_xlim()[1]/fig.axes[0].get_ylim()[1]
fig.axes[0].set_aspect(aspect=Ratio/2, adjustable='box')
#fig.tight_layout()
    

#%%


# if __name__ == "__main__":
#     main()
    
    