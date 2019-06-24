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
        mass = data[:,7]
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
        mass = data[7]
        masswater = data[14] # mass in water
        fwater = data[14]/data[7] #  water fraction
    # calculate the mass increase rate dmass
    #dm = np.diff(mass)
    #dt = np.diff(time)
    #dmdt = dm/dt
    return time, semi, ecc, inc, mass, fwater






##### defualt simulation infomation ####
# number of big planets
nbig = 1
big = True # plot big planets
generate_newdata = False
source = 'run1'
#source = 'run_md_11'
print('data file', source)


##### defualt figure parameters ####
tmin = 1e+2 # year
tmax = 2e+6 # year
amin = 1.# 1
amax = 50 #8
emin = 1e-6
emax = 0.02
imin = 1e-6
imax = 0.02
mmin = 3e-8#1e-11
mmax = 1e-4
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

# change into figure direction
os.chdir(workdir)

plt.clf()  # clear image
plt.close('all') # delete figure




#########################################
###### Main rountine: make figures ######
#########################################
for k in range(0,4): 
## k is the output figure type  ##
    if k == 0:
        output = 'mass_time' 
        print (output)
    if k == 1: 
        output = 'semi_time'
        print (output)
    if k == 2: 
        output = 'ecc_time'
        print (output)
    if k == 3: 
        output = 'inc_time'
        print (output)
    plt.figure(num=k,figsize=(8,5))
    if big == True and nbig >0:
        for i in range(1, nbig+1):
            filename = sourcedir+ 'P'+ str(i)+'.aei'
            time_big, semi_big, ecc_big, inc_big, mass_big, fwater_big = \
            readdata(filename)
            if output =='mass_time':
                plt.plot(time_big,mass_big,'b',linewidth=lw3)
                plt.xlabel(r'$\mathrm{ Time \ (yr)}$',fontsize=fs2)
                plt.ylabel('$ {\\rm Mass \\ (M_{\\oplus})}$',fontsize=fs2 )
                plt.ylim(mmin,mmax)
                plt.yticks([3.e-8,3.e-7,3e-6,3e-5,3e-4],['$10^{-2}$','$10^{-1}$','$1$','$10$','$100$'],fontsize=fs2)
            if output =='semi_time':
                plt.plot(time_big,semi_big,'b',linewidth=lw3)
                plt.xlabel('${\\rm Time \\ (yr)}$',fontsize=fs2)
                plt.ylabel('$ {\\rm Semimajor \\ axis \\ (AU)}$',fontsize=fs2)
                plt.ylim(amin,amax)
            if output =='ecc_time':
                plt.plot(time_big,ecc_big,'b',linewidth=lw3)
                plt.xlabel('${\\rm Time \\ (yr)}$',fontsize=fs2)
                plt.ylabel('$ {\\rm Eccentricity}$',fontsize=fs2 )
                plt.ylim(emin,emax)
            if output =='inc_time':
                plt.plot(time_big,ecc_big,'b',linewidth=lw3)
                plt.xlabel('${\\rm Time \\ (yr)}$',fontsize=fs2)
                plt.ylabel('$ {\\rm Inclination (radian)}$',fontsize=fs2 )
                plt.ylim(imin,imax)
        plt.xlim(tmin,tmax)
        plt.semilogx()
        plt.semilogy()
        plt.xticks([1e+2,1e+3,1e+4,1e+5,1e+6],['$10^{2}$','$10^{3}$','$10^{4}$','$10^{5}$','$10^{6}$'],fontsize=fs2)
    plt.savefig(source+'_'+output+'.pdf',orientation='landscape', format='pdf',bbox_inches='tight', pad_inches=0.1)


os.chdir(pydir)
