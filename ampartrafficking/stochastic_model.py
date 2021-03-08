# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:48:22 2020

@author: Moritz
"""


import numpy as np
import sys

def nearestNeighbours(PSD):
    """Returns the number of nearest neighbours on a grid.

    Parameters
    ----------
    PSD : array_like
        Grid with occupied and free elements, shape(N, N)

    Returns
    -------
    out: array_like
        Matrix containing the number of nearest neighbors for each element of the grid matrix.
        
    Examples
    --------
    
    Import libraries:
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import ampartrafficking.stochastic_model as sm
    
    Create and populate grid and calculate nearest neighbour matrix:

    >>> N=10
    >>> PSD=np.zeros((N,N))
    >>> PSD[np.random.randint(0,N,20),np.random.randint(0,N,20)]=1
    >>> 
    >>> NN=sm.nearestNeighbours(PSD)
    
    Plot:
        
    >>> plt.figure(figsize=(3,3), dpi=150)
    >>> plt.imshow(PSD)
    >>> plt.colorbar()
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')
    >>> 
    >>> plt.figure(figsize=(3,3), dpi=150)
    >>> plt.imshow(NN)
    >>> plt.colorbar()
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')

    Output:
        
    .. image:: images/example1_nearestNeighbours.png
        :width: 30%
    .. image:: images/example2_nearestNeighbours.png
        :width: 30%
    """
    
    n=np.shape(PSD)[0]
    m=np.shape(PSD)[1]
    
    PSD_border=np.zeros((n+2,m+2));#PSD_boarder is the matrix that contains PSD with the addition of one row or column to each border of PSD%
    PSD_border[1:n+1,1:m+1]=PSD!=0;#PSD also accounts for receptor types (bleached, not not bleached) such that occupied sites can have float values >0. Therefore it has to be checked where PSD does not equal zero to genererate the matrix of occupied sites PSD_boarder.
            
    PSD_up=PSD_border[0:n,1:m+1];# The matrix PSD_up represent the number of occupied nearest neighbors up to each site. 
    PSD_down=PSD_border[2:n+2,1:m+1];# The matrix PSD_down represent the number of occupied nearest neighbors down to each site.
    PSD_left=PSD_border[1:n+1,0:m];# The matrix PSD_left represent the number of occupied nearest neighbors on the left to each site.
    PSD_right=PSD_border[1:n+1,2:m+2];# The matrix PSD_right represent the number of occupied nearest neighbors on the right of each site.
    PSD_upright=PSD_border[0:n,2:m+2];# The matrix PSD_upright represent the number of occupied nearest neighbors up-right to each site.
    PSD_downleft=PSD_border[2:n+2,0:m];# The matrix PSD_downleft represent the number of occupied nearest neighbors down-left to each site.
    PSD_downright=PSD_border[2:n+2,2:m+2];# The matrix PSD_downright represent the number of occupied nearest neighbors down-right to each site.
    PSD_upleft=PSD_border[0:n,0:m];# The matrix PSD_upleft represent the number of occupied nearest neighbors up-left to each site.
    
    NN=PSD_up+PSD_down+PSD_left+PSD_right+PSD_upright+PSD_downright+PSD_downleft+PSD_upleft;#The matrix nearestNeighbors represents the total number of occupied nearest neighbors for each site.
    
    return NN.astype(int)





def kBUcoop(kBU, NN, PSD, typeID, beta=1.0):
    
    """Returns the cooperative unbinding rate per bound receptor.

    Parameters
    ----------
    kBU : float
        unbinding rate
    beta : float, optional
        By default beta=1.0. Factor by which the fraction of occupied nearest neighbours lowers the unbinding rate.
    NN : array_like
        Matrix that contains the number of nearest neighbours for each grid element.
    PSD : array_like
        Matrix representing the PSD grid
    typeID : float>0
        Receptor-type ID.
        

    Returns
    -------
    out: array_like
        Matrix containing the unbinding rates at each occupied grid element.
        
    Examples
    --------
    
    Import libraries:
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import ampartrafficking.stochastic_model as sm
    
    Create and populate grid and calculate nearest neighbour matrix:

    >>> kBU=0.1
    >>> typeID=2
    >>> N=10
    >>> PSD=np.zeros((N,N))
    >>> PSD[np.random.randint(0,N,20),np.random.randint(0,N,20)]=typeID
    >>> NN=sm.nearestNeighbours(PSD)
    
    Plot:
    
    >>> plt.figure(figsize=(4,3), dpi=150)
    >>> plt.plot(sm.kBUcoop(kBU, np.arange(0,9), np.array([typeID]*9), typeID))
    >>> plt.xlabel('number of nearest neighbours')
    >>> plt.ylabel('unbinding rate $k_{BU}^{coop}$')
    >>> 
    >>> fig=plt.figure(figsize=(3,2.25), dpi=150)
    >>> plt.imshow(PSD)
    >>> cbar=plt.colorbar()
    >>> cbar.set_label('occupied (type)', rotation=90, labelpad=10, y=0.5)
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')
    >>> 
    >>> fig=plt.figure(figsize=(3,2.25), dpi=150)
    >>> plt.imshow(NN)
    >>> cbar=plt.colorbar()
    >>> cbar.set_label('nearest neighbours', rotation=90, labelpad=10, y=0.5)
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')
    >>> 
    >>> fig=plt.figure(figsize=(3,2.25), dpi=150)
    >>> plt.imshow(sm.kBUcoop(kBU, NN, PSD, typeID))
    >>> cbar=plt.colorbar()
    >>> cbar.set_label('unbinding rate $k_{BU}^{coop}$', rotation=90, labelpad=10, y=0.5)
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')

    Output:
        
    .. image:: images/example1_kBUcoop.png
        :width: 45%
    .. image:: images/example2_kBUcoop.png
        :width: 45%
    .. image:: images/example3_kBUcoop.png
        :width: 45%
    .. image:: images/example4_kBUcoop.png
        :width: 45%
    """
    
    M=np.zeros(np.shape(PSD))
    occupied=np.where(PSD==typeID)
    Chi=np.arange(0,9)/8
    M[occupied]=((kBU*(1-Chi*beta))[NN])[occupied]
    
    return M





def kUBcoop(kUB, NN, PSD, alpha=16):
    
    """Returns the cooperative binding rate per mobile receptor.

    Parameters
    ----------
    kUB : float
        binding rate
    alpha : float, optional
        By default alpha=16. Factor by which the fraction of occupied nearest neighbours increases the binding rate.
    NN : array_like
        Matrix that contains the number of nearest neighbours for each grid element.
    PSD : array_like
        Matrix representing the PSD grid
        

    Returns
    -------
    out: array_like
        Matrix containing the binding rates at each unoccupied grid element.
        
    Examples
    --------
    
    Import libraries:
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import ampartrafficking.stochastic_model as sm
    
    Create and populate grid and calculate nearest neighbour matrix:

    >>> kUB=0.0005
    >>> N=10
    >>> PSD=np.zeros((N,N))
    >>> PSD[np.random.randint(0,N,20),np.random.randint(0,N,20)]=1
    >>> NN=sm.nearestNeighbours(PSD)
    
    Plot:
    
    >>> plt.figure(figsize=(4,3), dpi=150)
    >>> plt.plot(sm.kUBcoop(kUB, np.arange(0,9), np.array([0]*9)))
    >>> plt.xlabel('number of nearest neighbours')
    >>> plt.ylabel('binding rate $k_{UB}^{coop}$')
    >>> 
    >>> fig=plt.figure(figsize=(3,2.25), dpi=150)
    >>> plt.imshow(PSD)
    >>> cbar=plt.colorbar()
    >>> cbar.set_label('occupied (type)', rotation=90, labelpad=10, y=0.5)
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')
    >>> 
    >>> fig=plt.figure(figsize=(3,2.25), dpi=150)
    >>> plt.imshow(NN)
    >>> cbar=plt.colorbar()
    >>> cbar.set_label('nearest neighbours', rotation=90, labelpad=10, y=0.5)
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')
    >>> 
    >>> fig=plt.figure(figsize=(3,2.25), dpi=150)
    >>> plt.imshow(sm.kUBcoop(kUB, NN, PSD))
    >>> cbar=plt.colorbar()
    >>> cbar.set_label('binding rate $k_{UB}^{coop}$', rotation=90, labelpad=10, y=0.5)
    >>> plt.xlabel('position')
    >>> plt.ylabel('position')

    Output:
        
    .. image:: images/example1_kUBcoop.png
        :width: 45%
    .. image:: images/example2_kUBcoop.png
        :width: 45%
    .. image:: images/example3_kUBcoop.png
        :width: 45%
    .. image:: images/example4_kUBcoop.png
        :width: 45%
    """

    M=np.zeros(np.shape(PSD))
    free=np.where(PSD==0)
    Chi=np.arange(0,9)/8
    M[free]=((kUB*(alpha*Chi+1))[NN])[free]
    
    return M






def probabilityEval(Mub,Mbu,PSD,ID_basal,Mub_notBleached=None,Mbu_notBleached=None,ID_notBleached=None):
    
    """Returns the updated PSD Matrix and the corresponding number of receptors that got bound and unbound. To types, "basal" and "not bleached" can be considered, which is necessary when simulation FRAP.

    Parameters
    ----------
    Mub : array_like
        Matrix containing binding probabilities for the type "basal".
    Mbu : array_like
        Matrix containing unbinding probabilities for the type "basal".
    Mub_notBleached : array_like, optional
        By default None. Matrix containing binding probabilities for the type "not bleached".
    Mbu_notBleached : array_like, optional
        By default None. Matrix containing unbinding probabilities for the type "not bleached".
    PSD : array_like
        Matrix representing the PSD grid and its bound receptors.
    ID_basal : float
        Receptor ID of the basal pool.
    ID_notBleached: float
        Receptor ID of the not bleached pool.
        

    Returns
    -------
    out: float, float, float, float, array_like
        Number of receptors that got bound and unbound of the two types "basal" and "not bleached" and the updated PSD matrix.
        
    Examples
    --------
    
    Import libraries:
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import ampartrafficking.stochastic_model as sm
    
    Set parameters:

    >>> U=10
    >>> U_notBleached=10
    >>> kUB=0.005
    >>> kBU=1
    >>> N=10
    >>> ID_basal=1
    >>> ID_notBleached=2
    >>> dt=0.5
    
    Create and populate grid and calculate nearest neighbour matrix:
    
    >>> PSD=np.zeros((N,N))
    >>> while np.sum(PSD)<20*ID_basal:
    >>>     i=np.random.randint(0,N)
    >>>     j=np.random.randint(0,N)
    >>>     if PSD[i,j]==0:
    >>>         PSD[i,j]=ID_basal
    >>>         
    >>> while np.sum(PSD)<20*ID_basal+20*ID_notBleached:
    >>>     i=np.random.randint(0,N)
    >>>     j=np.random.randint(0,N)
    >>>     if PSD[i,j]==0:
    >>>         PSD[i,j]=ID_notBleached
    >>>         
    >>> NN=sm.nearestNeighbours(PSD)
    
    Plot PSD:
        
    >>> plt.figure()
    >>> plt.imshow(PSD)
    >>> plt.colorbar()
    
    Calculate probability Matrices and update the PSD Matrix:
        
    >>> Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal)*dt
    >>> Mub=sm.kUBcoop(kUB*U, NN, PSD)*dt
    >>> Mbu_notBleached=sm.kBUcoop(kBU, NN, PSD, ID_notBleached)*dt
    >>> Mub_notBleached=sm.kUBcoop(kUB*U_notBleached, NN, PSD)*dt
    >>>
    >>> PSD,dBoff,dBon,dBoff_notBleached,dBon_notBleached=sm.probabilityEval(Mub,Mbu,PSD,ID_basal,Mub_notBleached,Mbu_notBleached,ID_notBleached)

    Plot PSD:
        
    >>> plt.figure()
    >>> plt.imshow(PSD)
    >>> plt.colorbar()
    
    Output: (left: before, right: after)
        
    .. image:: images/example1_probabilityEval.png
        :width: 45%
    .. image:: images/example2_probabilityEval.png
        :width: 45%
    """    
    
    n=np.shape(PSD)[0]
    m=np.shape(PSD)[1]
    
    R=np.random.rand(n,m)
    Mask_ub=R<Mub
    Mask_bu=R<Mbu
    
    if Mub_notBleached is not None:
        R=np.random.rand(n,m)
        Mask_ub_notBleached=R<Mub_notBleached
        Mask_bu_notBleached=R<Mbu_notBleached

    
    if Mub_notBleached is not None:
        R2=np.random.rand(n,m)
        ii_basal=np.where((Mask_ub==True)&(Mask_ub_notBleached==True)&(R2<0.5))
        ii_notBleached=np.where((Mask_ub==True)&(Mask_ub_notBleached==True)&(R2>=0.5))
        Mask_ub[ii_basal]=False
        Mask_ub_notBleached[ii_notBleached]=False
    
    
    dBoff=np.sum(Mask_bu)
    dBon=np.sum(Mask_ub)
    
    if Mub_notBleached is not None:
        dBoff_notBleached=np.sum(Mask_bu_notBleached)
        dBon_notBleached=np.sum(Mask_ub_notBleached)

    
    PSD[Mask_ub]=ID_basal
    PSD[Mask_bu]=0
    
    if Mub_notBleached is not None:
        PSD[Mask_ub_notBleached]=ID_notBleached
        PSD[Mask_bu_notBleached]=0
    
    
    if Mub_notBleached is not None:
        return PSD,dBoff, dBon, dBoff_notBleached, dBon_notBleached
    else:
        return PSD,dBoff, dBon
    






def update_mobilePool(U,pin,pout,dBoff,dBon, U_notBleached=None,pin_notBleached=None,pout_notBleached=None,dBoff_notBleached=None,dBon_notBleached=None):

    """Updates the value for the mobile receptor pool. When simulating FRAP, a second type "not bleached" is considered.

    Parameters
    ----------
    U : float
        Mobile AMPAR pool.
    pin : float
        Probability of a receptor to enter the spine's mobile pool.
    pout : float
        Probability of a receptor to leave the spine's mobile pool.
    dBoff : float
        Number of receptors that got unbound from the PSD grid.
    dBon : float
        Number of receptors that got bound to the PSD grid.
    U_notBleached : float, optional
        Mobile AMPAR pool (bleached).
    pin_notBleached : float, optional
        Probability of a not bleached receptor to enter the spine's mobile pool.
    pout_notBleached : float, optional
        Probability of a not bleached receptor to enter the spine's mobile pool.
    dBoff_notBleached : float, optional
        Number of not bleached receptors that got unbound from the PSD grid.
    dBon_notBleached : float, optional
        Number of not bleached receptors that got bound to the PSD grid.
        

    Returns
    -------
    out: float, float
        U, U_notBleached.
    """
    
    if np.random.rand()<pin:
        U+=1
    
    if U_notBleached is not None:
        if np.random.rand()<pin_notBleached:
            U_notBleached+=1
        
    if np.random.rand()<pout:
        U-=1
    
    if U_notBleached is not None:
        if np.random.rand()<pout_notBleached:
            U_notBleached-=1

    DeltaB=dBoff-dBon
    U+=DeltaB

    if U_notBleached is not None:
        DeltaB_notBleached=dBoff_notBleached-dBon_notBleached
        U_notBleached+=DeltaB_notBleached

    if U<0:
        U=0
        
    if U_notBleached is not None:
        if U_notBleached<0:
            U_notBleached=0
        return U, U_notBleached
    else:
        return U
        
    
    
def calcTimeStep(UFP_0,A_spine,kUB,alpha,kBU,kout,kin):

    """Returns integration time step dt. By default simulations are carried out at dt=0.5s. If parameter choices require a smaller time step, dt is set to 0.25s. If dt is still too large, the simulation is cancelled and an error message is displayed. In this case dt needs to be set manually.

    Parameters
    ----------
    UFP_0 : float
        Mobile receptor pool fixed points. Sets the influx of receptors into spine. 
    A_spine : float
        Spine surface area. 
    kUB : float
        bidning rate
    alpha : float
        cooperativity factor for the binding
    kBU : float
        unbidning rate
    kout : float
        Total rate at which receptors exit the spine membrane (e.g. kout+kendo).
    kin : float
        Total rate at which receptors enter the spine membrane (e.g. kexo*Sexo+kin).
    

    Returns
    -------
    float
        Integration time step dt. 
    """
    
    dt=0.5
    thr=0.5 
    #Note that a threshold for the probability of 0.5 (thr) is rather high.
    #However, it turned out to be sufficient for the current study. 
    #Also, below, the probabilities are calculated for values of U twice 
    #the fixed point (2*UFP_0) to account for fluctuations in U that are above the fixed point value. 
    
    if 2*UFP_0/A_spine*kUB*(alpha+1)*dt>thr or kBU*dt>thr or kout*2*UFP_0/A_spine*dt>thr or kin*dt>thr:
        dt=0.25

    if 2*UFP_0/A_spine*kUB*(alpha+1)*dt>thr or kBU*dt>thr or kout*2*UFP_0/A_spine*dt>thr or kin*dt>thr:
        dt=0.1

    if 2*UFP_0/A_spine*kUB*(alpha+1)*dt>thr or kBU*dt>thr or kout*2*UFP_0/A_spine*dt>thr or kin*dt>thr:
        print("Errors, Integration time step dt is too large!")
        sys.exit("Errors, Integration time step dt is too large!")

    print('dt=',dt, 'pmax=', max([2*UFP_0/A_spine*kUB*(alpha+1)*dt, kBU*dt, kout*2*UFP_0/A_spine*dt, kin*dt]))
    print('pUB=',2*UFP_0/A_spine*kUB*(alpha+1)*dt)
    print('pBU=',kBU*dt)
    print('pin=',kin*dt)
    print('pout=',kout*2*UFP_0/A_spine*dt)
        
    return dt

