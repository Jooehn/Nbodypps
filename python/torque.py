# -*- coding: utf-8 -*- 
#!/usr/bin/env python
# type I migration torque Paardekooper2011
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib
from  astrounit import *
from diskmodel  import *
from plotset  import *
import subprocess as sp
from  plotset import *

nint = 1e+2
k1 = 200
k2 = 200


M_s = 1
mdot_gas = 1e-7*M_s**2
L_s = M_s**2
alpha_v = 1e-3
alpha_d = 1e-3
kap = 1e-2
opt_vis = True

# calculate the transtion radius for two disk regions
r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
# calculate temperature
q = np.zeros(k2)
r = np.zeros(k1)
q_gap = np.zeros(k1)
Temp_vis = np.zeros(k1)
Temp_irr = np.zeros(k1)
Temp = np.zeros(k1)
Surf = np.zeros(k1)
torq = np.zeros((k2,k1))
f = np.zeros((k2,k1))
f_tot = np.zeros((k2,k1))
# get the temerature at two disk regions
for j in range(0,k1):
    r[j] = (j+1.0)/nint     # different location
    r[j] = 10**(np.log10(5e-2) - j*(np.log10(5e-2)-np.log10(50)) /(k1-1))     # different location
    Temp[j] = Tgas(r[j],mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
    Surf[j] = siggas(r[j],mdot_gas,L_s,M_s,alpha_v,kap, opt_vis)
for i in range(0,k2):
    q[i] = 10**(np.log10(3e-8) - i*(np.log10(3e-8)-np.log10(3e-4)) /(k2-1))  # different planetary mass
    for j in range(0,k1):
        torq[i][j], f[i][j] = torque(q[i],r[j],mdot_gas,L_s,M_s,alpha_v,alpha_d,kap,opt_vis)


r_t = r_trans*np.ones(k2)
r_s = r_snow*np.ones(k2)
plt.clf()  # clear image
plt.close('all') # delete figure
#################
### figure 3 ####
#################

levels = np.linspace(-4.5,4.5,300)
cmap='seismic'
fig, [ax1,ax2] = plt.subplots(2,sharex= True, sharey = False, gridspec_kw={'height_ratios':[2.,4]}, figsize=(6,5), num= 3)
fig.subplots_adjust(hspace=0.07)
yy = np.linspace(0.1,1e+4,len(r_t))
ax1.plot(r,Surf,linewidth=lw1, color='k')
ax1.plot(r_t,yy,linewidth=lw2, color='c', linestyle='solid')
ax1.plot(r_s,yy,linewidth=lw2, color='m', linestyle='solid')
ax1.semilogy()
ax1.semilogx()
ax1.set_yticks([1,100])
ax1.set_yticklabels(['$1$','$10^{2}$'],fontsize=fs1)
ax1.set_ylabel('$ \\rm \\Sigma_{g} \\ [gcm^{-2}]$',fontsize=fs1)
ax1.set_ylim(0.5,1.e+3)
axx = ax1.twiny()
axx.semilogx()
axx.set_xlim(ax1.get_xlim())
axx.set_xticks([0.1,1,5,10])
axx.set_xticklabels(['$0.1$','$1$','$5$','$10$'],fontsize=fs1)
axx.set_xlabel('$\\rm Radius \\ [AU]$',fontsize=fs1)
axx.set_xlim(0.1,40)

ax11 = ax1.twinx()
ax11.set_ylabel('$ \\rm T_{g} \\ [K]$',fontsize=fs1)
ax11.set_ylim(11,4e+3)
ax11.set_xlim(ax1.get_xlim()) 
ax11.set_yticks([100,1000])
ax11.set_yticklabels(['$100$','$1000$'],fontsize=fs1,color='b')
ax11.spines['right'].set_color('b')
ax11.yaxis.label.set_color('b')
ax11.semilogy()
ax11.semilogx()
#[r.set_color('b') for r in ax11.yaxis.get_ticklabels()]
ax11.plot(r,Temp,linewidth=lw1, color='b')

ax2.plot(r_t,q*M_s,linewidth=lw2, color='c', linestyle='solid')
ax2.text(1.04*r_trans,1.8e-5,'$\\rm r_{\\rm trans}$', fontsize= fs2,color='c')
ax2.plot(r_s,q*M_s,linewidth=lw2, color='m', linestyle='solid')
ax2.text(0.68*r_snow,1.8e-5,'$\\rm r_{\\rm ice}$', fontsize= fs2,color='m')
CS=ax2.contourf(r,q*M_s,f,levels,cmap=cmap,alpha=1)
CD=ax2.contour(r,q*M_s,f,0,colors=(0.3,0.3,0.3))
#ax2.set_xlabel('$\\rm Radius \\ [AU]$',fontsize=fs1)
ax2.set_ylabel('$\\rm Planet \\ mass \\ [M_{\\oplus}$]',fontsize=fs1)
ax2.set_ylim(3e-8,3e-4)
ax2.semilogy()
ax2.semilogx()
ax2.set_xlim(0.1,40)
ax2.tick_params(axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    labelbottom=False)
#ax2.set_xticks([0.2,1,5,10])
#ax2.set_xticklabels(['$0.2$','$1$','$5$','$10$'],fontsize=fs1)
ax2.set_yticks([3e-7,9e-7,3e-6,9e-6,3e-5,9e-5])
ax2.set_yticklabels(['$0.1$','$0.3$','$1$','$3$','$10$','$30$'],fontsize=fs1)
#fig.colorbar(CS,orientation='horizontal',aspect=40,fraction=0.1, ticks=[-2,0,2],pad=0.7) 
#fig.colorbar(CS,orientation='vertical',aspect=25,ticks=[-2,0,2],pad=0.1) 
cbax = fig.add_axes([0.15, 0.04, 0.75, 0.03]) # setup colorbar axes. 
cb1 = fig.colorbar(CS, cmap=cmap,cax=cbax, norm=levels,ticks=[-4,-2,0,2,4], orientation='horizontal') 

fig.savefig('migration_map.pdf',orientation='landscape', format='pdf',bbox_inches='tight', pad_inches=0.01)







sp.call(['mv', 'migration_map.pdf', '../figure/'])
