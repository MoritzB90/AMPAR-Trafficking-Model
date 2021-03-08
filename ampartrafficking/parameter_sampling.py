# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 17:12:34 2020

@author: Moritz
"""

import sys
sys.path.append('../')

import numpy as np
from scipy.integrate import odeint
import ampartrafficking.rate_model as rm

#%%

def parameterSampling(t,Init,sLTP,Cooperativity,P,kin,kout,Trials, kBU_max,aUB_max,TUB_max,kexo_max,aexo_max,Texo_max):

    """Returns evolution over time of bound AMPARs for normal E-LTP and LTP with exocytosis blockage for various Trials. Parameter values for exocytosis event, receptor binding and unbinding rate are drawn randomly for each Trial. Also, evolution over time of mobile receptors, binding rate and exocytosis event rate and spine area are returned as well as values for the baseline exocytosis event rate kexo0, the factor of change in exocytosis event rate amplitide during LTP-induction aexo and amplitude decay time constant Texo are returned for each Trial.
    
    Parameters
    ----------
    t : array_like
        Time over which the model is integrated.
    Init : array_like
        Initial state values for U, B and Sexo ([U(0),B(0),Sexo(0)]) .
    sLTP : 0, 1
        Specifies whether sLTP is accounted for (=1) or not (=0).
    Cooperativity : 0, 1
        Specifies whether cooperative recepotr binding is accounted for (=1) or not (=0).
    P : float
        Number of receptor binding sites
    kin : float
            Rate at which receptors hop from the dednritic membrane compartment onto the spine membrane compartment.
    kout : float 
            Rate at which receptors hop from the spine membrane compartment onto the dendritic membrane compartment.
    Trials : integer
        Number of trials.
    kBU_max : float
        Maximum of the receptor unbinding rate. Values are drawn randomly between 0 and kBU_max. kUB0 is calculated from the unbinding rate (kUB0_(kBU)).
    aUB_max : float
        Maximum of the factor by which the receptor binding rate increases during LTP-induction. Values are drawn randomly between 0 and aUB_max.
    TUB_max : float
        Maximum of the decay time by which the receptor binding rate decreases to the baseline after LTP-induction. Values are drawn randomly between 0 and TUB_max.
    kexo_max : float
        Maximum of the baseline receptor exocytosis event rate kexo0. Values are drawn randomly between 0 and kexo_max. kendo is calculated from the exocytosis event rate (kendo_(kexo0)).
    aexo_max : float
        Maximum of the factor by which the receptor exocytosis event rate increases during LTP-induction. Values are drawn randomly between 0 and aexo_max.
    Texo_max : float
        Maximum of the decay time by which the receptor exocytosis event rate decreases to the baseline after LTP-induction. Values are drawn randomly between 0 and Texo_max.
    
    Returns
    -------
    B_Tr : array_like
        Time evolution of bound receptors during normal E-LTP for different Trials; shape(Trials, len(t))
    B_ne_Tr : array_like
        Time evolution of bound receptors during E-LTP with exocytosis blockage for different Trials; shape(Trials, len(t))
    U_Tr : array_like
        Time evolution of mobile receptors for different Trials; shape(Trials, len(t))
    kUB_Tr : array_like
        Time evolution of the receptor binding rate for different Trials; shape(Trials, len(t))
    kexo_Tr : array_like
        Time evolution of the receptor exocytosis event rate for different Trials; shape(Trials, len(t))
    Aspine_Tr : array_like
        Time evolution of the spine surface area for different Trials; shape(Trials, len(t))
    kexo0_Tr : array_like
        Values of the receptor exocytosis event rate at baseline for different Trials; shape(Trials)
    aexo_Tr : array_like
        Values of the factor by which the receptor exocytosis event rate increases during LTP-induction for different Trials; shape(Trials)
    Texo_Tr : array_like
        Values of the decay time by which the receptor exocytosis event rate decreases to the baseline after LTP-induction for different Trials; shape(Trials)
    """
    
    kin_RE=0.1
    kout_RE=0.000615
    V0=0.08
    A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)

    kUB_Tr=[]
    kexo_Tr=[]
    B_Tr=[]
    U_Tr=[]
    
    B_ne_Tr=[]
    U_ne_Tr=[]
    
    kexo0_Tr=[]
    aexo_Tr=[]
    Texo_Tr=[]
    
    for n in range(Trials):
        
        if n%100==0:
            print('Trial:',n)
        
        #E-LTP with exocytosis:
        kBU=np.random.rand()*kBU_max
        kUB0=rm.kUB0_(Init,kBU,P,A0,Cooperativity)
            
        a_UB=np.random.rand()*aUB_max
        T_UB=np.random.rand()*TUB_max
        kUB=rm.Parameter(kUB0)
        kUB.timecourse(rm.Stim_Resp,[1,a_UB,5,T_UB])
        
        kexo0=np.random.rand()*kexo_max
        a_exo=np.random.rand()*aexo_max
        T_exo=np.random.rand()*Texo_max
        kexo=rm.Parameter(kexo0)
        kexo.timecourse(rm.Stim_Resp,[1,a_exo,25,T_exo])
        kendo=rm.kendo_(Init,kexo0,kin,kout,A0)
            
        Vspine=rm.Parameter(V0)
        Vspine.timecourse(rm.DV,[V0,True])
        
        Model=rm.Model_system()
        
        solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
        
        
        kUB.update(t)
        kUB_Tr.append(kUB.current_value)
        kexo.update(t)
        kexo_Tr.append(kexo.current_value)
        B_Tr.append(solve.T[1])
        U_Tr.append(solve.T[0])
        
        kexo0_Tr.append(kexo0)
        aexo_Tr.append(a_exo)
        Texo_Tr.append(T_exo)
    
    
        #E-LTP without exocytosis:
        kexo=rm.Parameter(0)
        kexo.timecourse(rm.Stim_Resp,[1,a_exo,25,T_exo])
        Vspine=rm.Parameter(V0)
        Vspine.timecourse(rm.DV,[V0,False])
        
        Model2=rm.Model_system()
        
        solve2,infodict2 = odeint(Model2.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
    
        B_ne_Tr.append(solve2.T[1])
        U_ne_Tr.append(solve2.T[0])
        
        
    B_Tr=np.array(B_Tr)
    B_ne_Tr=np.array(B_ne_Tr)
    U_Tr=np.array(U_Tr)
    kUB_Tr=np.array(kUB_Tr)
    kexo_Tr=np.array(kexo_Tr)
    
    kexo0_Tr=np.array(kexo0_Tr)
    aexo_Tr=np.array(aexo_Tr)
    Texo_Tr=np.array(Texo_Tr)
    
    if sLTP==1:
        Vspine=rm.Parameter(V0)
        Vspine.timecourse(rm.DV,[V0,True])
        Vspine.update(t)
        Aspine_Tr=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)
    else:
        Aspine_Tr=A0*np.ones(len(t))
    
    return B_Tr,B_ne_Tr, U_Tr, kUB_Tr,kexo_Tr,Aspine_Tr, kexo0_Tr,aexo_Tr,Texo_Tr


def ELTP(Time):
    """Returns time evolution of EPSPs during E-LTP in % of baseline. The function has been fitted to data from Barco et al. 2002.
    
    Parameters
    ----------
    Time : array_like
        Time in s.
    
    Returns
    -------
    float
        Time evolution of EPSPs during E-LTP in % of baseline.
    """
    return rm.Stim_Resp(Time,0.98,2.53,82.13,2.35*10**-7)+rm.Stim_Resp(Time,-0.02,0.97,3445.37,128.40)


def LTPnoExo(Time):
    """Returns time evolution of EPSPs during LTP with exocytosis blockage in % of baseline. The function has been fitted to data from Penn et al. 2017.
    
    Parameters
    ----------
    Time : array_like
        Time in s.
    
    Returns
    -------
    float
        Time evolution of EPSPs during LTP with exocytosis blockage in % of baseline.
    """
    return rm.Stim_Resp(Time,0.84,1.85,41.62,18.12)+rm.Stim_Resp(Time,-0.16,0.57,605.16,41.62)


def Matching(B_Basal,B_nE,Trials,Time,BFP,Percentage=0.005):
    """Returns the indices and distance measure for the 0.5% of Trials that best match with experimental data from Barco et al. 2002 and Penn et al. 2017.
    
    Parameters
    ----------
    B_Basal : array_like
        Time evolution of bound AMPARs during normal E-LTP for different Trials (shape(Trials,t)).
    B_nE : array_like
        Time evolution of bound AMPARs during LTP with exocytosis blockage for different Trials (shape(Trials,len(Time))).
    Trials : int
        Number of trials.
    Time : array_like
        Time in s.
    BFP : float
        Fixed point/Basline level of the number of bound receptors.
    Percentage : float, optional
        By default Percentage=0.005. Percentage of Trials that are returned that best match with experimental data.
    
    Returns
    -------
    GoodMatch_index : array_like
        Indices of the 0.5% of Trials that best match with experimental data.
    GoodMatch_Value : array_like
        Distance measure of the 0.5% of Trials that best match with experimental data.
    """
    
    MinOutput=int(np.round(Trials*Percentage, 0))#/20 #200# 100# 

    B_Basal_Ref=ELTP(Time)*BFP
    B_nE_Ref=LTPnoExo(Time)*BFP

    DeltaT=Time[1]-Time[0]

    Dist=np.sum(np.abs(B_Basal_Ref-B_Basal),axis=1)*DeltaT/(Time[-1]-Time[0])+np.sum(np.abs(B_nE_Ref-B_nE),axis=1)*DeltaT/(Time[-1]-Time[0])
    Indices=np.arange(0,len(Dist),dtype=int)
    
    C=zip(Dist,Indices)
    Dist,Indices=zip(*sorted(C))
    Dist=np.array(Dist)
    Indices=np.array(Indices)
    
    GoodMatch_index=Indices[0:MinOutput]
    GoodMatch_Value=Dist[0:MinOutput]
                
    return GoodMatch_index, GoodMatch_Value