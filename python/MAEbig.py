#!/usr/bin/env python
#make all plots for planetesmial grwoth, like semi-time, ecc-time and mass-time figure
import numpy as np
import math
import string
import os
import csv
import pylab
import sys
import matplotlib
import subprocess as sp
import matplotlib.pyplot as plt
from astrounit import *
from diskmodel import rtrans,rsnow
from m6_funcs import *
from plotset import *


###################################################
############# read data files into array ##########
###################################################
def readdata(filename):
    data = np.loadtxt(filename)
    if data.ndim != 1:
        time = data[:,0]/365.25 # in unit of year
        count = len(time) # total line of time array; or np.size(time)
        semi = data[:,1]
        ecc = data[:,2]
        inc = data[:,3]*np.pi/180
        peri = data[:,4] # periastron 
        node = data[:,5] # ascending node
        mean = data[:,6] #mean anormal
        mass = data[:,7]*msuntome
        masswater = data[:,14] # mass in water
        fwater = data[:,14]/data[:,7] #  water fraction
    else :
        time = data[0]/365.25 # in unit of year
        count = 1 # total line of time array; or np.size(time)
        semi = data[1]
        ecc = data[2]
        inc = data[3]*np.pi/180
        peri = data[4] # periastron 
        node = data[5] # ascending node
        mean = data[6] #mean anormal
        mass = data[7]*msuntome
        masswater = data[14] # mass in water
        fwater = data[14]/data[7] #  water fraction
    # calculate the mass increase rate dmass
    #dm = np.diff(mass)
    #dt = np.diff(time)
    #dmdt = dm/dt
    return time, semi, ecc, inc, mass, fwater



##### defualt simulation infomation ####
#If an input file is provided, we use the given inputs from it
try:
    sys.argv[1]
except (IndexError,NameError):
    source  = 'diskdeptest3'
else:
    source = sys.argv[1] #The vars.ini file
#    vars_   = process_input(inputfile)
#    source  = vars_[8][0].rstrip('\n')

big = True # plot big planets
generate_newdata = False
print('data file', source)

##### defualt figure parameters ####
tmin = 1e+5 # year
tmax = 2e+6 # year
amin = .1 # 1
amax = 100 #8
emin = 1e-7
emax = 1
imin = 1e-7
imax = 1
mmin = 1e-1#1e-11
mmax = 100
dmmin = 1e-16
dmmax = 3e-10
# linewidth 
lw1 = 0.25
lw2 = 0.5
lw3 = 0.85
lw4 = 2.5
lw5 = 5.5
# fontsize 
fs1 = 12
fs2 = 16
fs3 = 17

cdir = os.getcwd()#+'/' #current dir
sourcedir = cdir+'/../sim/'+source+'/'
pydir = cdir # python dir
workdir = cdir+'/../figure/'
## in source direction, generating the aei output file
os.chdir(sourcedir)

#We get the number of bodies
nbig = len(get_names())
#As well as the collision info from the simulation
collinfo = detect_merger()
#And finally the disk parameters used
alpha_d, alpha_v, mdot_gas, L_s, _, _, st, kap, M_s, opt_vis = get_disk_params()
# calculate the transtion radius for two disk regions
r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
r_snow = rsnow(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)

# change into figure direction
os.chdir(workdir)

plt.clf()  # clear image
plt.close('all') # delete figure

#### We get colour for the different planets
ccycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
if nbig > len(ccycle):
    ncycle = int(np.ceil(nbig/len(ccycle)))
    ccycle = np.concatenate([ccycle]*ncycle)
else:
    ccycle = ccycle[:nbig]

#########################################
###### Main rountine: make figures ######
#########################################

fig_all, ax = plt.subplots(2,2,figsize=(12,8))
axlist = np.ravel(ax)

for k in range(0,4): 
    fig, axes = plt.subplots(figsize=(8,5))
## k is the output figure type  ##
    if k == 0:
        output = 'mass_time' 
        print (output)
    if k == 1: 
        output = 'semi_time'
        
        axes.axhline(r_trans,linewidth=lw2, color='c', linestyle='--')
        axes.axhline(r_snow,linewidth=lw2, color='m', linestyle='--')
        axes.text(1.5e6,0.68*r_trans,'$\\rm r_{\\rm trans}$',color='c')
        axes.text(1.5e6,0.68*r_snow,'$\\rm r_{\\rm ice}$',color='m')
        
        axlist[k].axhline(r_trans,linewidth=lw2, color='c', linestyle='--')
        axlist[k].axhline(r_snow,linewidth=lw2, color='m', linestyle='--')
        axlist[k].text(1.5e6,0.68*r_trans,'$\\rm r_{\\rm trans}$',color='c')
        axlist[k].text(1.5e6,0.68*r_snow,'$\\rm r_{\\rm ice}$',color='m')
        print (output)
    if k == 2: 
        output = 'ecc_time'
        print (output)
    if k == 3: 
        output = 'inc_time'
        print (output)
    if big == True and nbig >0:
        for i in range(1, nbig+1):
            filename = sourcedir+ 'P'+ str(i)+'.aei'
            time_big, semi_big, ecc_big, inc_big, mass_big, fwater_big = \
            readdata(filename)
            if output =='mass_time':
                #We check if the planet reaches its gap opening mass at any point
                #of its evolution. If so, we increase the thickness of its
                #plot line width
                gm_idx = check_gapmass(semi_big,mass_big,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis)
                if gm_idx is None:
                    lw = lw3
                else:
                    lw = lw4
                axes.plot(time_big,mass_big,linewidth=lw,color=ccycle[i-1],\
                      zorder=1)
                axes.set_xlabel(r'$\mathrm{ Time \ (yr)}$')
                axes.set_ylabel('$ {\\rm Mass \\ (M_{\\oplus})}$')
                axes.set_ylim(mmin,mmax)
                #We also plot the data in a separate figure
                axlist[k].plot(time_big,mass_big,linewidth=lw,color=ccycle[i-1],\
                      zorder=1)
                axlist[k].set_xlabel(r'$\mathrm{ Time \ (yr)}$')
                axlist[k].set_ylabel('$ {\\rm Mass \\ (M_{\\oplus})}$')
                axlist[k].set_ylim(mmin,mmax)
                #We check if the planet has undergone collisions and plot them
                #if that's the case
                if i in collinfo[0]:  
                    idx = np.where(i==collinfo[0])[0]
                    axes.plot(collinfo[2,idx],collinfo[1,idx],'X',markersize=5,markerfacecolor=ccycle[i-1],\
                          markeredgecolor='k',markeredgewidth=0.5,zorder=2)
                    axlist[k].plot(collinfo[2,idx],collinfo[1,idx],'X',markersize=5,markerfacecolor=ccycle[i-1],\
                          markeredgecolor='k',markeredgewidth=0.5,zorder=2)
            if output =='semi_time':
                #We check if the planet reaches its gap opening mass at any point
                #of its evolution. If so, we increase the thickness of its
                #plot line width
                gm_idx = check_gapmass(semi_big,mass_big,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis)
                if gm_idx is None:
                    lw = lw3
                else:
                    lw = lw4
                axes.plot(time_big,semi_big,linewidth=lw,color=ccycle[i-1],\
                          zorder=1)
                axes.set_xlabel('${\\rm Time \\ (yr)}$')
                axes.set_ylabel('$ {\\rm Semimajor \\ axis \\ (AU)}$')
                axes.set_ylim(amin,amax)
                
                axlist[k].plot(time_big,semi_big,linewidth=lw,color=ccycle[i-1],\
                      zorder=1)
                axlist[k].set_xlabel('${\\rm Time \\ (yr)}$')
                axlist[k].set_ylabel('$ {\\rm Semimajor \\ axis \\ (AU)}$')
                axlist[k].set_ylim(amin,amax)
                #We check if the planet has reached isolation mass
                im_idx = check_isomass(semi_big,mass_big,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis)
                #If it has we plot the corresponding time and semi-major axis
                if im_idx is not None:
                    axes.plot(time_big[im_idx],semi_big[im_idx],'p',markersize=5,markerfacecolor=ccycle[i-1],\
                          markeredgecolor='k',markeredgewidth=0.5,zorder=2)
                    axlist[k].plot(time_big[im_idx],semi_big[im_idx],'p',markersize=5,markerfacecolor=ccycle[i-1],\
                          markeredgecolor='k',markeredgewidth=0.5,zorder=2)
            if output =='ecc_time':
                axes.plot(time_big,ecc_big,linewidth=lw3,color=ccycle[i-1])
                axes.set_xlabel('${\\rm Time \\ (yr)}$')
                axes.set_ylabel('$ {\\rm Eccentricity}$')
                axes.set_ylim(emin,emax)
                
                axlist[k].plot(time_big,ecc_big,linewidth=lw3,color=ccycle[i-1])
                axlist[k].set_xlabel('${\\rm Time \\ (yr)}$')
                axlist[k].set_ylabel('$ {\\rm Eccentricity}$')
                axlist[k].set_ylim(emin,emax)
            if output =='inc_time':
                axes.plot(time_big,inc_big,linewidth=lw3,color=ccycle[i-1])
                axes.set_xlabel('${\\rm Time \\ (yr)}$')
                axes.set_ylabel('$ {\\rm Inclination (radian)}$')
                axes.set_ylim(imin,imax)
                
                axlist[k].plot(time_big,inc_big,linewidth=lw3,color=ccycle[i-1])
                axlist[k].set_xlabel('${\\rm Time \\ (yr)}$')
                axlist[k].set_ylabel('$ {\\rm Inclination (radian)}$')
                axlist[k].set_ylim(imin,imax)
            
        axes.set_xlim(tmin,tmax)
        axes.semilogx()
        axes.semilogy()
#        axes.set_xticks([1e+2,1e+3,1e+4,1e+5,1e+6],['$10^{2}$','$10^{3}$','$10^{4}$','$10^{5}$','$10^{6}$'])
        
        axlist[k].set_xlim(tmin,tmax)
        axlist[k].semilogx()
        axlist[k].semilogy()
#        axlist[k].set_xticks([1e+2,1e+3,1e+4,1e+5,1e+6],['$10^{2}$','$10^{3}$','$10^{4}$','$10^{5}$','$10^{6}$'])
    
    fig.suptitle(r'$\tau_s = '+'{:.2E}'.format(st)+'$')
    add_date(fig)            
    fig.savefig(source+'_'+output+'.pdf',orientation='landscape', format='pdf',bbox_inches='tight', pad_inches=0.1)

fig_all.suptitle(r'$\tau_s = '+'{:.2E}'.format(st)+'$')
add_date(fig_all)
fig_all.savefig(source+'_'+'all'+'.pdf',orientation='landscape', format='pdf',bbox_inches='tight', pad_inches=0.1)

os.chdir(pydir)
