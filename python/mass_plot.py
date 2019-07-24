#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:12:25 2019

@author: jooehn
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from diskmodel import *
from astrounit import *
from labellines import labelLines
from plotset import *
from m6_funcs import get_disk_params

#If an input file is provided, we use the given inputs from it
try:
    sys.argv[1]
except (IndexError,NameError):
    source  = 'pa+mig10'
else:
    source = sys.argv[1] #The vars.ini file

cdir = os.getcwd()#+'/' #current dir
sourcedir = cdir+'/../sim/'+source+'/'
pydir = cdir # python dir
workdir = cdir+'/../figure/'
## in source direction, generating the aei output file
os.chdir(sourcedir)

avals = np.linspace(0.1,100,100)

_,alpha_v,mdot_gas,L_s,_,mdot_peb,st,kap,M_s,opt_vis = get_disk_params()

r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

# change into figure direction
os.chdir(workdir)

fig, ax = plt.subplots(figsize=(8,6))

alpha_d = np.array([1e-3,1e-4])

cols = ['tab:green','tab:blue']
lines = []

for i in range(len(alpha_d)):
    
    m_gap = gap_mass(avals,mdot_gas,L_s,M_s,alpha_d[i],alpha_v,kap,opt_vis)
    m_iso = m_gap/2.3
    m_opt = opt_mass(avals,mdot_gas,L_s,M_s,alpha_d[i],alpha_v,kap,opt_vis)
    
    line, = ax.plot(avals,m_gap,linestyle='-',color=cols[i],\
            label=r'$\alpha_t='+'{:.0e}'.format(alpha_d[i])+'$')
    lines.append(line)
    ax.plot(avals,m_iso,linestyle='--',color=cols[i])
    ax.plot(avals,m_opt,linestyle='-.',color=cols[i])

gaphandle, = ax.plot([],[],'k-',label='$M_\mathrm{gap}$')
isohandle, = ax.plot([],[],'k--',label='$M_\mathrm{iso}$')
opthandle, = ax.plot([],[],'k-.',label='$M_\mathrm{opt}$')
    
ax.axvline(r_trans,linewidth=lw2, color='c', linestyle='--')
ax.axvline(r_snow,linewidth=lw2, color='m', linestyle='--')
ax.text(1.04*r_trans,4,'$\\rm r_{\\rm trans}$', fontsize= fs2,color='c')
ax.text(1.04*r_snow,4,'$\\rm r_{\\rm ice}$', fontsize= fs2,color='m')
add_date(fig)    
    
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.1,100)
ax.set_ylim(1,200)

labelLines(lines,fontsize=12,xvals=[1,1])

ax.set_xlabel('$a\ \mathrm{[AU]}$')
ax.set_ylabel(r'$M\ \mathrm{[M}_\oplus]$')

ax.legend(handles = [gaphandle,isohandle,opthandle],prop={'size':13})