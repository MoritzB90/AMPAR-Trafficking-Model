# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 10:14:16 2020

@author: Moritz
"""

import numpy as np
import matplotlib.pyplot as plt

class Parameter_2():
    """This class defines some basic properties and functionality for parameters that change during LTP induction.
    
    Attributes
    ----------
    base_value : float
        Base value of the correspconding model parameter, i.e. before LTP-induction.

    Methods
    -------
    timecourse(function,pars)
        Sets up time course of the parameter during LTP-induction. Requires a function and the parameters of that function to be passed.
    update(t)
        Updates the parameter to the current value according to the function passed to timecourse.
    """

    def __init__(self,base_value):
        """    
        Parameters
        ----------
        base_value : float
            Base value of the correspconding model parameter, i.e. before LTP-induction.
        """
        self.base_value=base_value
        self.current_value=base_value
    
    def timecourse(self,function,pars):
        """
        Sets up time course of the parameter during LTP-induction. Requires a function and the parameters of that function to be passed.
        
        Parameters
        ----------
        function : function
            A function that describes the time course of the parameter during LTP-induction
        pars : array_like
            List of parameters for the function.
            
    
        Returns
        -------
        out: 
        """
            
        self.pars=np.array(pars,dtype=float)
        self.course=function
        
    def update(self,t):
        """
        Sets up time course of the parameter during LTP-induction. Requires a function and the parameters of that function to be passed.
        
        Parameters
        ----------
        t : float
            Time passed since the LTP induction stimulus.
    
        Returns
        -------
        out: Updates current value of the partameter.
        """
        self.current_value=self.course(t,*self.pars)
        
class Model_system_2():
    """This class contains the set of differential equations that describe AMPAR trafficking at spines.
    
    Attributes
    ----------

    Methods
    -------
    odes(Init,t,Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE)
        Defines the ODEs that describe the synapse.
    """    
    def __init__(self):
        """    
        Parameters
        ----------
        container : array_like
            container that can be filled with any parameter one is interested in during numerical integration, e.g. t, kexo or Vspine etc..
        """
        self.container=[]
        
    def odes(self,Init,t,Vspine,P,kin,kin2,kout,kout2,kexo,kexo2,kendo,kendo2,
             kUB,kUB2,kBU,kBU2,Cooperativity,sLTP,kin_RE,kout_RE,k12,k21):
        """
        Defines the ODEs that describe the AMPAR dynamics at the spine.
        
        Parameters
        ----------
        Init : [U(0),B(0),Sexo(0)]
            Initial conditions for the three variables of the system.
        t : array_like
            Time.
        Vspine : float
            Describes the spine volume and spine volume change during E-LTP.
        P : float
            Number of binding site/slots at the PSD.
        kin : float
            Rate at which receptors hop from the dednritic membrane compartment onto the spine membrane compartment.
        kout : float 
            Rate at which receptors hop from the spine membrane compartment onto the dendritic membrane compartment.
        kexo : float
            Rate of exocytosis events occuring at the spine.
        kendo : float
            Rate at which receptors are endocytosed at the spine.
        kUB : float
            Rate at which AMPARs bind to PSD slots.
        kBU : float
            Rate at which AMPARs unbind from PSD slots.
        Cooperativity: 1,0
            States whether binding is cooperative (1) or not cooperative (0).
        sLTP : 1,0
            States whether slTP does occur (1) or does not occur (0).
        kin_RE : float
            Rate at which AMPAR containing endosomes enter the spine (dynamics of exocytosis event size Sexo).
        kout_RE : float
            Rate at which AMPAR containing endosomes leave the spine.
            
        Returns
        -------
        out: [dU,dB,dS_exo]
        """
    
        self.container.append([Vspine.current_value,kexo.current_value,kUB.current_value,t])
        
        U=Init[0]
        B=Init[1]
        S_exo=Init[2]
        U2=Init[3]
        B2=Init[4]
        S_exo2=Init[5]
            
        kexo.update(t)
        kUB.update(t)
        
        kexo2.update(t)
        kUB2.update(t)
        
        k12.update(t)
        
        if sLTP==1:
            Vspine.update(t)
        
        Aspine=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)

        dS_exo=kin_RE-kout_RE*S_exo/Vspine.base_value-k12.current_value*S_exo+k21*S_exo2#-kexo.current_value*S_exo
        
        dS_exo2=-kout_RE*S_exo2/Vspine.base_value+k12.current_value*S_exo-k21*S_exo2#-kexo2.current_value*S_exo2
        
        
        dU=kexo.current_value*S_exo+kin+kBU_(B+B2,kBU,Cooperativity,P)*B-kendo*U/Aspine-kout*U/Aspine-kUB_(B+B2,kUB.current_value,Cooperativity,P)*(P-(B+B2))*U/Aspine
        
        dU2=kexo2.current_value*S_exo2+kin2+kBU_(B+B2,kBU2,Cooperativity,P)*B2-kendo2*U2/Aspine-kout2*U2/Aspine-kUB_(B+B2,kUB2.current_value,Cooperativity,P)*(P-(B+B2))*U2/Aspine
        
        dB=kUB_(B+B2,kUB.current_value,Cooperativity,P)*(P-(B+B2))*U/Aspine-kBU_(B+B2,kBU,Cooperativity,P)*B
        
        dB2=kUB_(B+B2,kUB2.current_value,Cooperativity,P)*(P-(B+B2))*U2/Aspine-kBU_(B+B2,kBU2,Cooperativity,P)*B2
    
        return [dU,dB,dS_exo,dU2,dB2,dS_exo2]
    
    
def kUB_(B,kUB,Cooperativity,P):
    """Returns the receptor binding rate kUB for either the cooperative or the non-cooperative binding model.

    Parameters
    ----------
    B : float
        Number of bound receptors.
    kUB : float
        Rate at which AMPARs bind to PSD slots.
    Cooperativity : 0, 1
        Specifies whether cooperative recepotr binding is accounted for (=1) or not (=0).
    P : float
        Number of binding site/slots at the PSD.
        
    Returns
    -------
    float
        kUB.
    """
    if Cooperativity==1:       
        m = 24.6/(12.5 + P);
        return kUB*(m*B**0.8 + 1);
    else:
        return kUB

def kBU_(B,kBU,Cooperativity,P):
    """Returns the receptor unbinding rate kUB for either the cooperative or the non-cooperative binding model.

    Parameters
    ----------
    B : float
        Number of bound receptors.
    kBU : float
        Rate at which AMPARs unbind from PSD slots.
    Cooperativity : 0, 1
        Specifies whether cooperative recepotr binding is accounted for (=1) or not (=0).
    P : float
        Number of binding site/slots at the PSD.
        
    Returns
    -------
    float
        kUB.
    """
    if Cooperativity==1:
        Lambda = 0.9*P + 15.5;
        Beta = 0.6*P + 9.3;
        return kBU*(Lambda/(Beta + B)-0.5);
    else:
        return kBU


def kUB0_(Init,kBU,P,Aspine,Cooperativity):
    """Returns the receptor binding rate at baseline, which is calculated from the fixed point equation of the AMPAR trafficking model.

    Parameters
    ----------
    Init : array_like
        Initial state values for U, B and Sexo ([U(0),B(0),Sexo(0)]) .
    kBU : float
        Rate at which AMPARs unbind from PSD slots.
    P : float
        Number of binding site/slots at the PSD.
    Aspine : float
        Spine surface area.
    Cooperativity : 0, 1
        Specifies whether cooperative recepotr binding is accounted for (=1) or not (=0).

    Returns
    -------
    float
        Receptor binding rate at baseline kUB0.
    """
    
    if Cooperativity==0:
        kUB0=((Init[1]*kBU)/((P-(Init[1]+Init[4]))*Init[0]))*Aspine
    else:
        kUB0=(Aspine*Init[1]*kBU*((15.5+0.9*P)/((Init[1]+Init[4])+9.3+0.6*P)-0.5))/((P-(Init[1]+Init[4]))*(1+(24.6*(Init[1]+Init[4])**0.8)/(12.5+P))*Init[0])
                    
    return kUB0

def kendo_(Init,kexo0,kin,kout,Aspine):
    """Returns the receptor endocytosis rate, which is calculated from the fixed point equation of the AMPAR trafficking model.

    Parameters
    ----------
    Init : array_like
        Initial state values for U, B and Sexo ([U(0),B(0),Sexo(0)]) .
    kexo0 : float
        Rate of exocytosis events occuring at the spine under basal conditions.
    kin : float
        Rate at which receptors hop from the dednritic membrane compartment onto the spine membrane compartment.
    kout : float 
        Rate at which receptors hop from the spine membrane compartment onto the dendritic membrane compartment.
    Aspine : float
        Spine surface area.

    Returns
    -------
    float
        Receptor endocytosis rate kendo.
    """
    return (kexo0*Init[2]+kin)*Aspine/Init[0]-kout

