# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:43:19 2020

@author: Moritz
"""

"""Fig 6B

This script reproduces the plots seen in Fig 6B of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation", which shows the coefficient of variation of the cooperative model for differenr values of bound receptors B and PSD sizes P.

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
plt.rcParams['svg.fonttype'] = 'none'

#%%

# def main():

SaveFig=0
    
A_spine_basal=0.898
P_basal=70
 
duration=100000#Duration in s
kBU=0.1
kout=0.018
kin=0.02
kexo0=0#0.0018
kendo=0.002058
S_exo=13

Nr_Trials=1

beta=1#0#
alpha=16#0#
kUB=0.0005#*7
kBU=0.1#/100


N_List=[6,8,10]
D_List=np.array([3,5,7,10,12,14,18,30])
        


#%%

B_N=[]
U_N=[]

ID_basal=1

for N in N_List:
    
    A_spine=N**2/P_basal*A_spine_basal
    
    B_U=[]
    U_U=[]

    for D_0 in D_List:
        
        UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0)
        dt=sm.calcTimeStep(UFP,A_spine,kUB,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo)
        
        
        B_Tr=[]
        U_Tr=[]

        
        for Trial in range(0,Nr_Trials):
            
            if Trial%10==0:
                print('Trial:',Trial)
                print('P:',N**2, 'U:', UFP)
                
            
            PSD=np.zeros((N,N))

            Time=[]
            B_t=[]
            U_t=[]

            U=UFP

            for t in np.arange(0,duration+dt,dt):
                
                
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal)*dt
                Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD)*dt

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
BCV_N=[]
for i,N in enumerate(N_List):
    BMean_U=[]
    BCV_U=[]
    for j,D_0 in enumerate(D_List):
        
        DeltaT=len(Time)-int(len(Time)/5)
        BMean=np.mean(np.mean(np.array(B_N)[i,j,0::,int(len(Time)/5)::], axis=1))
        BStd=np.mean(np.std(np.array(B_N)[i,j,0::,int(len(Time)/5)::], axis=1))
        
        print(BMean)
        print(BStd)
        print(BStd/BMean*100)
        
        BCV_U.append(BStd/BMean*100)
        BMean_U.append(BMean)

    BCV_N.append(BCV_U)
    BMean_N.append(BMean_U)


#%%

beta=1
alpha=0
kUB=0.0036

N_noC_List=[10]
D_noC_List=np.array([0.5,1,3,6,10,15,20,30,50])#np.array([0.5,1,3,6,15,30,70])
        

B_noC_N=[]
U_noC_N=[]

ID_basal=1

for N in N_noC_List:
    
    A_spine=N**2/P_basal*A_spine_basal
    
    B_noC_U=[]
    U_noC_U=[]

    for D_0 in D_noC_List:
        
        UFP=np.round(rm.UFP_(kexo0,S_exo,kendo,kin*D_0,kout,A_spine),0)
        dt=sm.calcTimeStep(UFP,A_spine,kUB,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo)
        
        
        B_noC_Tr=[]
        U_noC_Tr=[]

        
        for Trial in range(0,Nr_Trials):
            
            if Trial%10==0:
                print('Trial:',Trial)
                print('P:',N**2, 'U:', UFP)
                
            
            PSD=np.zeros((N,N))

            Time=[]
            B_noC_t=[]
            U_noC_t=[]

            U=UFP

            for t in np.arange(0,duration+dt,dt):
                
                
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta=beta)*dt
                Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD, alpha=alpha)*dt

                PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
                
                
                pout=(kout*U/A_spine+kendo*U/A_spine)*dt
                pin=(kin*D_0+kexo0*S_exo)*dt
                U=sm.update_mobilePool(U,pin,pout,dBoff,dBon)
                
                if t%0.5==0:
                    Time.append(t)
                    B_noC_t.append(np.sum(PSD==1))
                    U_noC_t.append(U)
                
            B_noC_Tr.append(B_noC_t)
            U_noC_Tr.append(U_noC_t)

        B_noC_U.append(B_noC_Tr)
        U_noC_U.append(U_noC_Tr)

    B_noC_N.append(B_noC_U)
    U_noC_N.append(U_noC_U)

#%%

BMean_noC_N=[]
BCV_noC_N=[]
for i,N in enumerate(N_noC_List):
    BMean_noC_U=[]
    BCV_noC_U=[]
    for j,D_0 in enumerate(D_noC_List):
        
        DeltaT=len(Time)-int(len(Time)/5)
        BMean=np.mean(np.mean(np.array(B_noC_N)[i,j,0::,int(len(Time)/5)::], axis=1))
        BStd=np.mean(np.std(np.array(B_noC_N)[i,j,0::,int(len(Time)/5)::], axis=1))
        
        print(BMean)
        print(BStd)
        print(BStd/BMean*100)
        
        BCV_noC_U.append(BStd/BMean*100)
        BMean_noC_U.append(BMean)

    BCV_noC_N.append(BCV_noC_U)
    BMean_noC_N.append(BMean_noC_U)
    
#%%


col=sns.color_palette('colorblind')

fig=plt.figure(figsize=(3.5,2), dpi=150)
for i,x,y in (zip(range(len(N_List)),BMean_N,BCV_N)):
    plt.plot(x,y, marker='o', markersize=2, linewidth=1.5, linestyle='-', color=col[i+1], label='$P={:1.0f}$'.format(N_List[i]**2))

for i,x,y in (zip(range(len(N_noC_List)),BMean_noC_N,BCV_noC_N)):
    plt.plot(x,y, linewidth=1.5, linestyle='-.', color='k', label='$P={:1.0f}$'.format(N_noC_List[i]**2)+', $k_{UB}^{coop}=const.$')

plt.xlim(0,70)
plt.ylim(0,100)
plt.xlabel('Bound AMPARs $B$ (#)')
plt.ylabel('Coefficient of variation \n CV (%)')
plt.legend(prop=fontLgd)
sns.despine(top=True)
fig.tight_layout()

if SaveFig==1:
    fig.savefig('Figures\\Fig6B.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig6B.svg', bbox_inches="tight", dpi=400)
    

#%%


# if __name__ == "__main__":
#     main()



#%%  

# B_N=np.array(B_N)
# U_N=np.array(U_N)

# col=sns.color_palette('colorblind')[3::]

# dn=1#25

# from scipy import stats

# for i,N in enumerate(N_List):
    
#     A_spine=N**2/P_basal*A_spine_basal
    
#     for j,D_0 in enumerate(D_List):
        
#         DeltaT=len(Time)-int(len(Time)/5)
#         BMean=np.mean(np.mean(np.array(B_N)[i,j,0::,int(len(Time)/5)::], axis=1))
#         UMean=np.mean(np.mean(np.array(U_N)[i,j,0::,int(len(Time)/5)::], axis=1))
        
#         print(BMean)
#         print(UMean)
        
#         B_err=stats.sem((B_N[i][j][0:,0::dn]+U_N[i][j][0:,0::dn])/(BMean+UMean)*100)
#         Y=(np.mean(B_N[i][j][0:,0::dn], axis=0)+np.mean(U_N[i][j][0:,0::dn], axis=0))/(BMean+UMean)*100
#         # plt.errorbar(np.array(Time[0::dn])-duration/2/60,Y, yerr=B_err, label=r'P={0:1.0f}, $\langle B_0 \rangle$={1:1.0f}'.format(N**2, BMean))
        
#         fig=plt.figure(figsize=(3.5,2), dpi=150)
#         plt.plot(np.array(Time[0::dn])/60,Y, color=col[i], linewidth=2, label=r'P={0:1.0f}, $\langle B \rangle$={1:1.0f}'.format(N**2, BMean))
        
#         plt.text(0.3,0.55,r'$\langle U/A_{spine} \rangle$='+'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine), transform=fig.axes[0].transAxes)
                
#         plt.axvline(int(len(Time)/5)*0.5/60)
#         plt.axhline(100, color='k')
#         # plt.xlim(0,35)
#         # plt.ylim(0,120)
#         plt.xlabel('Time (min)')
#         plt.ylabel('AMPAR \n Recovery (%)')
#         plt.legend(ncol=2, prop=fontLgd)
#         sns.despine()
#         fig.tight_layout()