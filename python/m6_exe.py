#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:33:31 2019

@author: jooehn
"""

import numpy as np
import os
from subprocess import call
from astrounit import *
from sys import argv
from m6_funcs import *

############### Input ###############

#If an input file is provided in the call, e.g. 'python m6_exe.py vars.ini'
#we read the input from there. Otherwise, we use the input given in this script

try:
    argv[1]
except (IndexError,NameError):
    source  = 'migtest' #Source directory
    
    N       = 1   #Number of planet embryos
    amin    = 3  #Min semi-major axis val
    astep   = 5   #Step in between planets
    mrange  = np.array([50,50])  #Planetary mass in earth masses
    
    T       = 2e5 #End time in yr
else:
    inputfile   = argv[1] #The vars.ini file
    vars_       = process_input(inputfile)
    N           = int(vars_[0][0])
    mrange      = np.array(vars_[1][0].split(',')[:],dtype=float)
    amin        = float(vars_[2][0])
    astep       = float(vars_[3][0])
    T           = float(vars_[4][0])
    st          = float(vars_[5][0])
    source      = vars_[6][0].rstrip('\n')
############### Setup ###############

#We first swap into the source directory
cdir = os.getcwd()#+'/' #current dir
sourcedir = cdir+'/../sim/'+source+'/'
pydir = cdir # python dir
workdir = cdir+'/../figure/'

os.chdir(sourcedir)

#We set up the big.in file with the number of planetary embryos we want

#First we generate the semi-major axis values given our separation and a min value

amax  = amin + astep*N
avals = np.arange(amin,amax,astep)

#We also generate the masses in the range provided
pmass = np.linspace(mrange[0],mrange[1],N)

#Next we set up the physical properties and the phase of the planets and
#store them in an array of size (N,10)
mp = pmass*(Mearth/Msun)
ep = np.random.rayleigh(1e-2,N) #We draw eccentricities from a Rayleigh distr.
ip = 0.5*ep
rp = 1.0
dp = 1.5
xp = 1.5e-9

props = np.array([0,rp,dp,xp,0,ep,ip,0,0,0])

bigdata    = np.zeros((N,10))
bigdata[:] = props

#We also insert the semi-major axis values and masses
bigdata[:,0] = mp
bigdata[:,4] = avals

#Finally we generate the names and write it all into our big.in file
bignames = ['P{}'.format(i+1) for i in range(N)]

big_input(bignames,bigdata,asteroidal=True)

#We set up the end time of the run which is divided into shorter integrations of 
#length dt
dt  = 1e4
setup_end_time(dt)
#Furthermore, we set the input Stokes number
set_stokes_number(st)

#We also remove old files
bad_ext = ['*.dmp','*.tmp','*.aei','*.clo','*.out']        
for j in bad_ext:
    call(['find',os.getcwd(),'-maxdepth','1','-type','f','-name',j,'-delete'])

############### Executing ###############
print('Initiating integration')
t = 0
while t < T:
    
    #We first perform the integration by calling MERCURY
    
    call(['./mercury'])    
        
    #Since planets residing at our cavity at 0.2 AU will slow the integration
    #down significantly, we check if there is such a planet and 
    #remove it if this is indeed the case
    all_in_cavity = check_cavity_planet()
    if all_in_cavity:
        print('All planets have entered cavity after {} yr'.format(t+dt))
        print('Stopping integration')
        break
    extend_stop_time(dt)
    t += dt

print('Creating .aei files')
if not all_in_cavity:
    call(['./element'])
#Finally, we create figures
os.chdir(pydir)
print('Creating figures')
call(['python','MAEbig.py','vars.ini'])
print('The deed is done')