# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 16:37:38 2020

@author: Moritz
"""

"""Fig 2

This script reproduces the plots seen in Fig 2, Fig 3C and D and Fig 7G of "The biophysical basis underlying 
the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within 
the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from ampartrafficking import rate_model as rm
import seaborn as sns
sns.set_style("ticks")
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('x-small')


#%%

# def main():

SaveFig=1

#Color palettes used in plots:
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[3]
col2=sns.dark_palette(sns.color_palette("colorblind", 10)[1],4)[3]
col3=sns.dark_palette(sns.color_palette("colorblind", 10)[8],4)

#%%

#Import experimental Data for comparison with the model. Taken from Penn et al 2017 and Barco et al. 2002 (see Paper for reference details):

Penn2017_x=np.loadtxt("Data\\Penn2017_x.csv", delimiter=",")
Penn2017_y=np.loadtxt("Data\\Penn2017_y.csv", delimiter=",")
Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")

#%%

#Set parameter values and initial conditions:

Init = [10,20,13]

P=70
kin=0.2
kout=0.018
kexo0=0.0018
kBU=0.1
kin_RE=0.1
V0=0.08
A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)

kout_RE=kin_RE*V0/Init[2]#0.000615
kendo=rm.kendo_(Init,kexo0,kin,kout,A0)#0.002058

U_trace=[]
U_trace_label=[]
t=np.arange(0,130*60)

#%%

#Basic AMPAR trafficking

sLTP=0
Cooperativity=0
kUB0=rm.kUB0_(Init,kBU,P,A0,Cooperativity)#0.00359

#E-LTP with exocytosis:
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,30,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])

Model=rm.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=rm.Parameter(0)
kexo.timecourse(rm.Stim_Resp,[1,0,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,False])

Model=rm.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, label='E-LTP', legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, label='LTP, $k_{exo}=0$', legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/Init[1]*100, color=col1, linewidth=2, label='Basic Model')
plt.plot(t/60,solve2.T[1]/Init[1]*100, color=col2, linewidth=2, label='Basic Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.05,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2A.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2A.svg', bbox_inches="tight", dpi=400)
    

U_trace.append((solve.T[0]/A0)/(Init[0]/A0)*100)
U_trace_label.append('Basic Model')


#%%

#sLTP

sLTP=1
Cooperativity=0
kUB0=rm.kUB0_(Init,kBU,P,A0,Cooperativity)#0.00359

#E-LTP with exocytosis:
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,30,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])

Model=rm.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=rm.Parameter(0)
kexo.timecourse(rm.Stim_Resp,[1,0,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,False])

Model=rm.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/Init[1]*100, color=col1, linewidth=2, label='sLTP Model')
plt.plot(t/60,solve2.T[1]/Init[1]*100, color=col2, linewidth=2, label='sLTP Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.4,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2B.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2B.svg', bbox_inches="tight", dpi=400)


fig=plt.figure(figsize=(3.5,2), dpi=150)
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])
Vspine.update(t)
Aspine=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)

U_trace.append((solve.T[0]/Aspine)/(Init[0]/A0)*100)
U_trace_label.append('sLTP Model')


#%%

#Cooperative receptor binding

sLTP=0
Cooperativity=1
kUB0=rm.kUB0_(Init,kBU,P,A0,Cooperativity)#0.0005048

#E-LTP with exocytosis:
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,5,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])

Model=rm.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=rm.Parameter(0)
kexo.timecourse(rm.Stim_Resp,[1,0,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,False])

Model=rm.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/Init[1]*100, color=col1, linewidth=2, label='Coop. Model')
plt.plot(t/60,solve2.T[1]/Init[1]*100, color=col2, linewidth=2, label='Coop. Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.4,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2C.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2C.svg', bbox_inches="tight", dpi=400)


U_trace.append((solve.T[0]/A0)/(Init[0]/A0)*100)
U_trace_label.append('Coop. Model')


#%%

#sLTP + cooperative receptor binding

sLTP=1
Cooperativity=1
kUB0=rm.kUB0_(Init,kBU,P,A0,Cooperativity)#0.0005048

#E-LTP with exocytosis:
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,5,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])

Model=rm.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=rm.Parameter(0)
kexo.timecourse(rm.Stim_Resp,[1,0,25,60])
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,False])

Model=rm.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/Init[1]*100, color=col1, linewidth=2, label='sLTP+Coop. Model')
plt.plot(t/60,solve2.T[1]/Init[1]*100, color=col2, linewidth=2, label='sLTP+Coop. Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.3,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2D.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2D.svg', bbox_inches="tight", dpi=400)


fig=plt.figure(figsize=(3.5,2), dpi=150)
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])
Vspine.update(t)
Aspine=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)

U_trace.append((solve.T[0]/Aspine)/(Init[0]/A0)*100)
U_trace_label.append('sLTP+Coop. Model')


#%%

fig=plt.figure(figsize=(3.5,2), dpi=150)
for i,u,l in zip(np.arange(4),U_trace,U_trace_label):

    plt.plot(t/60,u, color=col3[i], linewidth=2, label=l)

lgd=plt.legend(bbox_to_anchor=(0.45,0.35,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Mobile AMPAR conc. \n $U/A_{spine}$ (%)')
plt.xlim(-3,130)
plt.ylim(80,159)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig3D.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig3D.svg', bbox_inches="tight", dpi=400)
    
#%%
########################################################################
#Fig3C

Init = [10,20,13]

sLTP=1
Cooperativity=0
kUB0=rm.kUB0_(Init,kBU,P,A0,Cooperativity)#0.00359
kUB=rm.Parameter(kUB0)
kUB.timecourse(rm.Stim_Resp,[1,30,5,60])
kexo=rm.Parameter(kexo0)
kexo.timecourse(rm.Stim_Resp,[1,5,25,60])

#%%V0=0.08
V0=0.08
A0=4*np.pi*(3*V0/(4*np.pi))**(2/3)
kout_RE=kin_RE*V0/Init[2]
Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])

Model=rm.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

Vspine=rm.Parameter(V0)
Vspine.timecourse(rm.DV,[V0,True])
Vspine.update(t)

#%%V0=0.15

V0b=0.15
A0=4*np.pi*(3*V0b/(4*np.pi))**(2/3)
kout_RE=kin_RE*V0b/Init[2]
Vspine2=rm.Parameter(V0b)
Vspine2.timecourse(rm.DV,[V0b,True])

Model=rm.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine2,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

Vspine2=rm.Parameter(V0b)
Vspine2.timecourse(rm.DV,[V0b,True])
Vspine2.update(t)

#%%V0=0.3

V0c=0.3
A0=4*np.pi*(3*V0c/(4*np.pi))**(2/3)
kout_RE=kin_RE*V0c/Init[2]
Vspine3=rm.Parameter(V0c)
Vspine3.timecourse(rm.DV,[V0c,True])

Model=rm.Model_system()

solve3,infodict3 = odeint(Model.odes,Init,t,args=(Vspine3,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

Vspine3=rm.Parameter(V0c)
Vspine3.timecourse(rm.DV,[V0c,True])
Vspine3.update(t)

#%%

col=sns.color_palette("rocket", 8)

fig=plt.figure(figsize=(3.5,2), dpi=150)


fig=plt.figure(figsize=(3.5,2), dpi=150)
plt.plot(t/60, solve.T[2]/Init[2]*100, color='k', alpha=1, linewidth=2, label='$S_{exo}$', zorder=1)
plt.plot(t/60, solve2.T[2]/Init[2]*100, color='k',alpha=0.75, linewidth=2, zorder=1)
plt.plot(t/60, solve3.T[2]/Init[2]*100, color='k', alpha=0.5, linewidth=2, zorder=1)
plt.plot(t/60, Vspine.current_value/V0*100, color=col[5], linewidth=2, label='$V_{spine}$', zorder=0)
plt.plot(t/60, Vspine2.current_value/V0b*100, color=col[6], linewidth=2, zorder=0)
plt.plot(t/60, Vspine3.current_value/V0c*100, color=col[7], linewidth=2, zorder=0)
sns.despine()
plt.legend()
plt.xlabel('Time (min)')
plt.ylabel('% of baseline')
#lgd=plt.figure.legend(bbox_to_anchor=(0.88,0.88,0,0), loc="upper right", mode="normal", borderaxespad=0, ncol=1, prop=fontLgd)

plt.axhline(100, color='k', linestyle='--', linewidth=0.5)

plt.text(0.2,0.9,'$V_{spine}^{0}=$'+'{} $\mu m^3$'.format(V0), transform=fig.axes[0].transAxes)
plt.text(0.2,0.7,'{} $\mu m^3$'.format(V0b), transform=fig.axes[0].transAxes)
plt.text(0.2,0.5,'{} $\mu m^3$'.format(V0c), transform=fig.axes[0].transAxes)

plt.xlim(-10,8000/60)    
sns.despine(right=True)

fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig3C.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig3C.svg', bbox_inches="tight", dpi=400)
    
#%%

#Fig7H

#%%V0=0.08
V0=0.11
Vspine4=rm.Parameter(V0)
Vspine4.timecourse(rm.DV,[V0,True])
Vspine4.update(t)

#%%Blocked exocytosis
V0=0.11
Vspine_noExo=rm.Parameter(V0)
Vspine_noExo.timecourse(rm.DV,[V0,False])
Vspine_noExo.update(t)

#%%

Matsuzaki2004_x=np.loadtxt("Data\\Matsuzaki2004_x.csv", delimiter=",") #Matsuzaki et al. that spines with an initial volume smaller 0.1 mu m^3 stably increase in size for hours after LTP induction.
Matsuzaki2004_y=np.loadtxt("Data\\Matsuzaki2004_y.csv", delimiter=",")

Yang2008_x=np.loadtxt("Data\\Yang2008_x.csv", delimiter=",") #Yang et al. show that the maintenance of spine enlargement after LTP induction depends on exocytosis. In the absence, the spine volume decays to its baseline within minutes.
Yang2008_y=np.loadtxt("Data\\Yang2008_y.csv", delimiter=",")

fig,axes=plt.subplots(figsize=(3.5,2), dpi=150)
col=sns.color_palette("rocket", 8)

plot0=sns.lineplot(Matsuzaki2004_x-40, Matsuzaki2004_y, color='k', marker="o", label='sLTP, Matsuzaki et al. 2004', legend=False, zorder=0)
plot1=sns.lineplot(Yang2008_x, Yang2008_y, color='k', marker="s", label='sLTP+Botox, Yang et al. 2008', legend=False, zorder=0)

plt.plot(t/60, Vspine4.current_value/V0*100, color=col[4], linewidth=2, label='sLTP (Early Phase); $V_{spine}^{0}=$'+str(V0)+r'$\mu m^3$', zorder=1)
plt.plot(t/60, Vspine_noExo.current_value/V0*100, color=col[6], linewidth=2, label='sLTP, $k_{exo}=0$', zorder=1)

axes.set_xlabel('Time (min)')
axes.set_ylabel('$V_{spine}$ (%)')

lgd=axes.figure.legend(bbox_to_anchor=(0.88,0.88,0,0), loc="upper right", mode="normal", borderaxespad=0, ncol=1, prop=fontLgd)

axes.axhline(100, color='k', linestyle='--', linewidth=0.5)

plt.xlim(-10,130)    
sns.despine(right=True)

fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7G.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7G.svg', bbox_inches="tight", dpi=400)
    
#%%


# if __name__ == "__main__":
#     main()
