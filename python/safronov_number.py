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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm

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

r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow  = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

fig,ax = plt.subplots(figsize=(10,6))

masses  = np.linspace(0.1,50,100)
avals   = np.linspace(0.1,100,100)

xx, yy = np.meshgrid(avals,masses)

safvals = safronov_number(yy,M_s,xx)

levels = np.logspace(np.log10(safvals.min()),1,100)

contax = ax.contourf(xx,yy,safvals,levels=levels,cmap='plasma',norm=LogNorm())
cont1 = ax.contour(xx,yy,safvals,levels = [1],colors=['w'],linestyles='--',norm=LogNorm())

divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.1)

cbar = plt.colorbar(contax,cax=cax)
cbar.set_label(r'$\Theta$')
cbar.set_ticks([1e-2,1e-1,1e0,1e1])

ax.set_xlabel('$a\ [\mathrm{AU}]$')
ax.set_ylabel(r'$M_p\ [M_\oplus]$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,100)
ax.set_ylim(0.1,50)

fmt = {}
c_label = [r'$\Theta = 1$']
for l, s in zip(cont1.levels, c_label):
    fmt[l] = s

ax.clabel(cont1, cont1.levels, inline=True, fmt=fmt, colors='w', fontsize=12,\
          manual=True)

ax.axvline(r_trans,linewidth=lw2, color='c', linestyle='--')
ax.axvline(r_snow,linewidth=lw2, color='k', linestyle='--')
ax.text(1.04*r_trans,.2,'$\\rm r_{\\rm trans}$', fontsize= fs2,color='c')
ax.text(1.04*r_snow,.2,'$\\rm r_{\\rm ice}$', fontsize= fs2,color='k')
add_date(fig)