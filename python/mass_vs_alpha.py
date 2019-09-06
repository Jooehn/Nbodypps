#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 14:00:18 2019

@author: John Wimarsson

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
    source  = 't2migtest4'
else:
    source = sys.argv[1] #The vars.ini file

cdir = os.getcwd()#+'/' #current dir
sourcedir = cdir+'/../sim/'+source+'/'
pydir = cdir # python dir
workdir = cdir+'/../figure/'
## in source direction, generating the aei output file
os.chdir(sourcedir)

_,alpha_v,mdot_gas,L_s,_,mdot_peb,st,kap,M_s,opt_vis = get_disk_params()

r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

# change into figure direction
os.chdir(workdir)

fig, ax = plt.subplots(figsize=(8,6))

#alpha_d = np.array([1e-3,1e-4])
alpha_d = np.logspace(-5,1,500)

cols = ['tab:green','tab:blue','tab:orange']
    
m_gap = gap_mass(r_trans,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis)
m_iso = m_gap/2.3
m_opt = opt_mass(r_trans,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis)
    
ax.plot(alpha_d,m_gap,linestyle='-',color=cols[0],label='$M_\mathrm{gap}$')
ax.plot(alpha_d,m_iso,linestyle='-',color=cols[1],label='$M_\mathrm{iso}$')
ax.plot(alpha_d,m_opt,linestyle='-',color=cols[2],label='$M_\mathrm{opt}$')

add_date(fig)    
    
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e-5,1)
ax.set_ylim(1,2e3)

labelLines(ax.get_lines(),fontsize=12,xvals = [1e-4,1e-3,1e-2])

ax.set_xlabel(r'$\alpha_t$')
ax.set_ylabel(r'$M\ \mathrm{[M}_\oplus]$')
ax.set_title('$\mathrm{Mass\ limits\ at\ }r_\mathrm{trans}$')