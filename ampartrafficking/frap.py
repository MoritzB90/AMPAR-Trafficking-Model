# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:35:31 2020

@author: Moritz
"""

import sys
sys.path.append('../')

import numpy as np
import ampartrafficking.stochastic_model as sm
import ampartrafficking.rate_model as rm

#%%

def D_FP(kendo,kexo0,S_exo_0,kout,kin,Conc):
    return ((kendo+kout)*Conc-kexo0*S_exo_0)/kin

#%%



def FRAP(N_List, Conc_List, beta, alpha, kUB, kBU, duration, Nr_Trials):

    """Returns evolution over time of photobleached and not photobleached mobile and bound receptors (FRAP simulation).

    Parameters
    ----------
    N_List : array_like
        List of PSD sizes P
    Conc_List : array_like
        List of target fixed points for the mobile receptor concentration. Sets the influx of receptors into spine. 
    beta : float
        cooperativity factor for the unbinding. (Should be set to 1 or 0)
    alpha : float
        cooperativity factor for the binding
    kUB : float
        bidning rate
    kBU : float
        unbidning rate
    duration : float
        Duration of the simulation
    Nr_Trials : integer
        Number of trials.

    Returns
    -------
    B_N : array_like
        Time evolution of bound receptors for different conditions and number of Trials; shape(len(N_List), len(Conc_List), Trials, duration/0.5+1)
    U_N : array_like
        Time evolution of mobile receptors for different conditions and number of Trials; shape(len(N_List), len(Conc_List), Trials, duration/0.5+1)
    B_notBleached_N
        Time evolution of bound receptors not photobleached for different conditions and number of Trials; shape(len(N_List), len(Conc_List), Trials, duration/0.5+1)
    U_notBleached_N
        Time evolution of mobile receptors not photobleached for different conditions and number of Trials; shape(len(N_List), len(Conc_List), Trials, duration/0.5+1)
    PSD : array_like
        Matrix representing the PSD grid and its receptors at the end of the simulation (bleached, not bleached).
    Time : arrayl_like
        Time; shape(duration/0.5+1,)
    
    """
    
    A_spine_basal=0.898
    P_basal=70
    
    kout=0.018
    kin=0.02
    kexo0=0.0018
    kendo=0.002058
    S_exo_0=13
    
    t_bleaching=duration/2# 1 #
    
    
    #%%
    
    B_N=[]
    U_N=[]
    B_notBleached_N=[]
    U_notBleached_N=[]
    
    ID_basal=1
    ID_notBleached=2
    
    for N in N_List:
        
        A_spine=N**2/P_basal*A_spine_basal
        
        B_U=[]
        U_U=[]
        B_notBleached_U=[]
        U_notBleached_U=[]

        for Conc in Conc_List:
            
            D_0=D_FP(kendo,kexo0,S_exo_0,kout,kin,Conc)
                    
            UFP=np.round(rm.UFP_(kexo0,S_exo_0,kendo,kin*D_0,kout,A_spine),0)
            dt=sm.calcTimeStep(UFP,A_spine,kUB,alpha,kBU,kout+kendo, kin*D_0+kexo0*S_exo_0)
            
            
            B_Tr=[]
            U_Tr=[]
            B_notBleached_Tr=[]
            U_notBleached_Tr=[]
            
            for Trial in range(0,Nr_Trials):
                
                if Trial%10==0:
                    print('Trial:',Trial)
                    print('P:',N**2, 'U:', UFP)
                    
                
                PSD=np.zeros((N,N))
    
                D=D_0
                D_notBleached=0
                S_exo=S_exo_0
                S_exo_notBleached=0
                
                Time=[]
                B_t=[]
                U_t=[]
                B_notBleached_t=[]
                U_notBleached_t=[]
    
                U=UFP
                U_notBleached=0
    
                for t in np.arange(0,duration+dt,dt):
                    
                    if t>t_bleaching:
                        D=0
                        D_notBleached=D_0
                        S_exo=0
                        S_exo_notBleached=S_exo_0
                    
                    NN=sm.nearestNeighbours(PSD)
                    
                    Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
                    Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD, alpha)*dt
                    Mbu_notBleached=sm.kBUcoop(kBU, NN, PSD, ID_notBleached, beta)*dt
                    Mub_notBleached=sm.kUBcoop(kUB*U_notBleached/A_spine, NN, PSD, alpha)*dt
    
                    PSD,dBoff,dBon,dBoff_notBleached,dBon_notBleached=sm.probabilityEval(Mub,Mbu,PSD,ID_basal,Mub_notBleached,Mbu_notBleached,ID_notBleached)
                    
                    
                    pout=(kout*U/A_spine+kendo*U/A_spine)*dt
                    pout_notBleached=(kout*U_notBleached/A_spine+kendo*U_notBleached/A_spine)*dt
                    pin=(kin*D+kexo0*S_exo)*dt
                    pin_notBleached=(kin*D_notBleached+kexo0*S_exo_notBleached)*dt
                    U,U_notBleached=sm.update_mobilePool(U,pin,pout,dBoff,dBon, U_notBleached,pin_notBleached,pout_notBleached,dBoff_notBleached,dBon_notBleached)
                    
                    if t%0.5==0:
                        Time.append(t)
                        B_t.append(np.sum(PSD==ID_basal))
                        U_t.append(U)
                        B_notBleached_t.append(np.sum(PSD==ID_notBleached))
                        U_notBleached_t.append(U_notBleached)
                    
                B_Tr.append(B_t)
                U_Tr.append(U_t)
                B_notBleached_Tr.append(B_notBleached_t)
                U_notBleached_Tr.append(U_notBleached_t)
                
            B_U.append(B_Tr)
            U_U.append(U_Tr)
            B_notBleached_U.append(B_notBleached_Tr)
            U_notBleached_U.append(U_notBleached_Tr)
    
        B_N.append(B_U)
        U_N.append(U_U)
        B_notBleached_N.append(B_notBleached_U)
        U_notBleached_N.append(U_notBleached_U)
        
    return np.array(B_N), np.array(U_N), np.array(B_notBleached_N), np.array(U_notBleached_N), PSD, np.array(Time)/60