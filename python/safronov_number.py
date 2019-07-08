#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 14:26:51 2019

@author: jooehn
"""

import numpy as np
import matplotlib.pyplot as plt
from diskmodel import *
from astrounit import *
from labellines import labelLines
from plotset import *
from m6_funcs import safronov_number

msuntome = Msun/Mearth

M_s      = 1
mdot_gas = 1e-7*M_s**2
mdot_peb = 4e-4 #earthmasses/yr
L_s      = M_s**2
alpha_v  = 1e-3
alpha_d  = 1e-3
hg       = 0.03
kap      = 1e-2
G        = 4*np.pi**2
opt_vis  = True
st       = 1e-1

a = np.linspace(0.1,100,50)

r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow  = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

fig,ax = plt.subplots(figsize=(8,6))

masslist = [1,10,30,50,100]

clist = ['tab:blue','tab:purple','tab:orange','tab:green','tab:red','tab:cyan']

for m,c in zip(masslist,clist):
    
    saf = safronov_number(m/msuntome,M_s,a)
    ax.plot(a,saf,color=c,label=r'$'+'{}'.format(m)+'\ M_\oplus $')
#    ax[1].plot(a,tpebvals,'-.',color=c,label=r'$'+'{}'.format(m)+'\ M_\oplus $')

ax.set_xlabel('$a\ [\mathrm{AU}]$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,100)
ax.set_ylim(0.1,10)

ax.axvline(r_trans,linewidth=lw2, color='c', linestyle='--')
ax.axvline(r_snow,linewidth=lw2, color='m', linestyle='--')
ax.axhline(1,linewidth=lw2,color='tab:gray',linestyle='--')

ax.text(1.04*r_trans,.5,'$\\rm r_{\\rm trans}$', fontsize= fs2,color='c')
ax.text(0.6*r_snow,.5,'$\\rm r_{\\rm ice}$', fontsize= fs2,color='m')

ax.set_ylabel(r'$\Theta$')
fig.subplots_adjust(wspace=0.1)

labelLines(ax.lines,fontsize=10)
add_date(fig)