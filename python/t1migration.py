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
L_s      = M_s**2
alpha_v  = 1e-3
alpha_d  = 1e-3
hg       = 0.03
kap      = 1e-2
G        = 4*np.pi**2
opt_vis  = True

def tmig_1(mp,a):
    
    """Computes the type I migration time for a given planetary mass at a given
    semi-major axis.
    
    mp: mass given in earth masses
    a : semi-major axis in AU"""
    
    sigma_g = siggas(a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)*(AU**2/Msun)
    
    t = (M_s/(mp/msuntome))*(M_s/(sigma_g*a**2))*hg**2*(a**3/(G*M_s))**(0.5)
    
    return t

a = np.linspace(0.1,100,50)

r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow  = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

fig = plt.figure(figsize=(8,6))

masslist = [1,10,30,50,100]

for m in masslist:
    plt.plot(a,tmig_1(m,a),label=r'$'+'{}'.format(m)+'\ M_\oplus $')

plt.xlabel('$a\ [\mathrm{AU}]$')
plt.ylabel(r'$\tau_I\ [\mathrm{yr}]$')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1,100)
plt.ylim(1e2,2e6)

plt.axvline(r_trans,linewidth=lw2, color='c', linestyle='--')
plt.axvline(r_snow,linewidth=lw2, color='m', linestyle='--')
plt.text(1.04*r_trans,5e5,'$\\rm r_{\\rm trans}$', fontsize= fs2,color='c')
plt.text(0.68*r_snow,5e5,'$\\rm r_{\\rm ice}$', fontsize= fs2,color='m')

labelLines(plt.gca().get_lines(),fontsize=12)
add_date(fig)