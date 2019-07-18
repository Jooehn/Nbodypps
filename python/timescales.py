#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 10:07:49 2019

@author: jooehn
"""

import numpy as np
import matplotlib.pyplot as plt
from diskmodel import *
from astrounit import *
from labellines import labelLines
from plotset import *

msuntome = Msun/Mearth

M_s      = 1
mdot_gas = 1e-7*M_s**2
mdot_peb = 4e-4 #earthmasses/yr
L_s      = M_s**2
alpha_v  = 1e-3
alpha_d  = 1e-3
kap      = 1e-2
G        = 4*np.pi**2
opt_vis  = True
st       = 1e-1

a = np.linspace(0.1,100,50)

r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow  = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

fig,ax = plt.subplots(1,2,figsize=(12,6),sharey=True)

masslist = [0.1,1,10,30,50,100]

clist = ['tab:blue','tab:purple','tab:orange','tab:green','tab:red','tab:brown']

for m,c in zip(masslist,clist):
    
    tmig1vals = tmig_1(m,a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis=True)
    tpebvals  = tpeb(m,a,mdot_gas,mdot_peb,L_s,M_s,alpha_d,alpha_v,kap,st,opt_vis=True)
    ax[0].plot(a,tmig1vals,color=c,label=r'$'+'{}'.format(m)+'\ M_\oplus $')
    ax[1].plot(a,tpebvals,'-.',color=c,label=r'$'+'{}'.format(m)+'\ M_\oplus $')

for i in ax:
    i.set_xlabel('$a\ [\mathrm{AU}]$')
    i.set_xscale('log')
    i.set_yscale('log')
    i.set_xlim(0.1,100)
    i.set_ylim(1e3,1e8)

    i.axvline(r_trans,linewidth=lw2, color='c', linestyle='--')
    i.axvline(r_snow,linewidth=lw2, color='m', linestyle='--')

    i.text(1.04*r_trans,3e6,'$\\rm r_{\\rm trans}$', fontsize= fs2,color='c')
    i.text(0.6*r_snow,3e6,'$\\rm r_{\\rm ice}$', fontsize= fs2,color='m')

ax[0].set_ylabel(r'$\tau_I\ [\mathrm{yr}]$')
ax[1].set_ylabel(r'$\tau_\mathrm{PA}\ [\mathrm{yr}]$')
fig.subplots_adjust(wspace=0.1)


labelLines(ax[0].lines,fontsize=10)
add_date(fig)