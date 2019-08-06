#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:21:35 2019

@author: John Wimarsson

Some functions that handles MERCURY input files
"""

import numpy as np
import os
import math
from tempfile import mkstemp
from subprocess import call
from shutil import move
from astrounit import *
from diskmodel import iso_mass

def setup_end_time(T,T_start=0):
    
    """Small function that sets up the duration of our integration.
        
        T: the total time of the integration given in yr
        T_start: we can also specify the start time. If no start time
            is given, it is set to zero by default. Should aso be given in yr"""
        
    #The string we want to change
    start_str = ' start time (days) = '
    end_str   = ' stop time (days) = '
    
    #Makes temporary file
    fh, abs_path = mkstemp()
    
    with os.fdopen(fh,'w') as new_file:
        with open('param.in') as old_file:
            for line in old_file:
                if start_str in line:
                    old_sstr = line
    
#                    old_stime = float(old_sstr.strip(start_str))
                    new_stime = T_start
                
                    new_sstr = start_str+str(new_stime)+'\n'
                    new_file.write(line.replace(old_sstr, new_sstr))
                    
                elif end_str in line:
                    old_estr = line
    
                    etime = T*365.25
                    
                    new_estr = end_str+str(etime)+'\n'
                    new_file.write(line.replace(old_estr, new_estr))
                else:
                    new_file.write(line)
    #Remove original file and move new file
    os.remove('param.in')
    move(abs_path, 'param.in')
    
def extend_stop_time(T):
    
    """Small function that updates the stop time in param.dmp to allow for an
    extended integration in case we have no collisions. Updates the old time
    value by adding the value T."""
    
    #The string we want to change
    stime_str = '  stop time (days) =    '
    
    #Makes temporary file
    fh, abs_path = mkstemp()
    
    with os.fdopen(fh,'w') as new_file:
        with open('param.dmp') as old_file:
            for line in old_file:
                if stime_str in line:
                    old_str = line
    
                    old_time = float(old_str.strip(stime_str))
                    new_time = old_time+T*365.25
                    
                    rep_str = stime_str+str(old_time)
                    new_str = stime_str+str(new_time)
                    new_file.write(line.replace(rep_str, new_str))
                else:
                    new_file.write(line)
    #Remove original file and move new file
    os.remove('param.dmp')
    move(abs_path, 'param.dmp')
    
def big_input(names,bigdata,asteroidal=False,epoch=0):
    
    """Function that generates the big.in input file for MERCURY6 given an Nx10
    array of data in the following format:
        
        Columns:
            
            0: mass of the object given in solar masses
            1: radius of the object in Hill radii
            2: density of the object
            3: x for the planet
            4: semi-major axis in AU
            5: eccentricity
            6: inclination in degrees
            7: argument of pericentre in degrees
            8: longitude of the ascending node
            9: mean anomaly in degrees
    
    We can also pass the argument asteroidal as True if we want that coordinate
    system. Also the epoch can be specified, it should be given in years."""
    
    N = len(bigdata)       
    
    if asteroidal:
        style = 'Asteroidal'
    else:
        style = 'Cartesian'
    
    initlist = [')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n',\
        ") Lines beginning with `)' are ignored.\n",\
        ')---------------------------------------------------------------------\n',\
        ' style (Cartesian, Asteroidal, Cometary) = {}\n'.format(style),\
        ' epoch (in days) = {}\n'.format(epoch*365.25),\
        ')---------------------------------------------------------------------\n']
    
    with open('big.in','w+') as bigfile:
        
        for i in initlist:
            bigfile.write(i)
        
        for j in range(N):
            
            bigfile.write(' {0:11}m={1:.17E} r={2:.0f}.d0 d={3:.2f} x={4:.2e}\n'.format(names[j],*bigdata[j,0:4]))
            bigfile.write(' {0: .17E} {1: .17E} {2: .17E}\n'.format(*bigdata[j,4:7]))
            bigfile.write(' {0: .17E} {1: .17E} {2: .17E}\n'.format(*bigdata[j,7:]))
            bigfile.write('  0. 0. 0.\n')
            
def process_input(file_): #Function that reads the input file if any
    h_file = []
    input_ = open(file_, 'r')
    for line in input_:
        if line.find("#") != -1:
            continue
        elif line.find("\n") == 0:
            continue
        else:
            h_file.append(line.split('\t'))
    return h_file

def check_timestep():
    """Checks if the timestep has gone below values of 1 day. Returns either
    True or False depending on the timestep size."""
    
    with open('param.dmp') as param:
        for line in param:
            if 'timestep' in line:
                timestep = float(line.split()[-1])
                print('Current timestep: {} days'.format(timestep))
                if timestep < 1:
                    return True
                else:
                    return False
            
def get_names(dmp=False):
    """Returns a list with the names of the current objects active in our simulation."""
    
    if dmp:
        file = 'big.dmp'
    else:
        file = 'big.in'
    
    with open(file,'r') as bigfile:
        biglines = bigfile.readlines()
        
    nameid = np.arange(6,len(biglines),4)
        
    names = [biglines[i].split()[0] for i in nameid]
    
    return names

def check_cavity_planet():
    """Removes planet that has entered our cavity in order to keep the timestep
    from converging at low values and thereby slowing down our integration.
        names: The names of our planets
        data: The corresponding data input data for the planets
    Returns True if all planets have entered the cavity and the integration
    should be terminated."""
    
    #We generate output files
    call(['./element'])
    names = get_names(True)
    N = len(names)
    in_cavity = np.full(N,False)
    
    #We check the element file for the semi-major axis values of our planets
    with open('element.out') as element:
        elelines = element.readlines()
        for i in range(len(elelines[2:])):
            params = elelines[2:][i].split()
            x = float(params[9])
            y = float(params[10])
            z = float(params[11])
            a = np.sqrt(x**2+y**2+z**2)
            #If a is less than our cavity value which is 0.2 AU we register the
            #corresponding planet
            if a <= 0.2:
                in_cavity[i] = True
    #We extract the ids of the planets that are in the cavity and remove them
    if np.all(in_cavity):
        return True
    elif np.any(in_cavity):
        cavity_ids  = np.where(in_cavity)[0]
        old_names = []
        for i in range(len(cavity_ids)):
            print('We remove P{} from the integration\n'.format(names[cavity_ids[i]]))
            old_names.append(names[cavity_ids[i]])
            names.remove(names[cavity_ids[i]])
        #We loop through each line in big.dmp and remove the information of the
        #planets inside the cavity
        
        #Makes temporary file
        fh, abs_path = mkstemp()
    
        with os.fdopen(fh,'w') as new_file:
            with open('big.dmp') as bigfile:
                biglines = bigfile.readlines()
                i = 6
                while (len(biglines)>6) & (len(old_names)>0):
                    if old_names[0] in biglines[i]:
                        del biglines[i:i+4]
                        del old_names[0]
                    else:
                        i += 4
            new_file.writelines(biglines)
        #Remove original file and move new file
        call(['find',os.getcwd(),'-maxdepth','1','-type','f','-name','*.aei','-delete']) 
        os.remove('big.dmp')
        move(abs_path, 'big.dmp')    
        return False
    else:
        call(['find',os.getcwd(),'-maxdepth','1','-type','f','-name','*.aei','-delete'])
        return False
    
def set_stokes_number(st):
    
        #The string we want to change
    st_str = " pebble's stokes number =    "
    
    #Makes temporary file
    fh, abs_path = mkstemp()
    
    with os.fdopen(fh,'w') as new_file:
        with open('param.in') as old_file:
            for line in old_file:
                if st_str in line:
                    rep_str = line
                    new_str = st_str+'{:.0E}'.format(st)+'\n'
                    new_file.write(line.replace(rep_str, new_str))
                else:
                    new_file.write(line)
    #Remove original file and move new file
    os.remove('param.in')
    move(abs_path, 'param.in')
    
def set_turbulent_alpha(alpha):
    
        #The string we want to change
    turb_str = " turbulent strength alphad =     "
    
    #Makes temporary file
    fh, abs_path = mkstemp()
    
    with os.fdopen(fh,'w') as new_file:
        with open('param.in') as old_file:
            for line in old_file:
                if turb_str in line:
                    rep_str = line
                    new_str = turb_str+'{:.17E}'.format(alpha)+'\n'
                    new_file.write(line.replace(rep_str, new_str))
                else:
                    new_file.write(line)
    #Remove original file and move new file
    os.remove('param.in')
    move(abs_path, 'param.in')
    
def get_disk_params():
    """Finds the disk parameters given in the param.in file and returns them
    as an array in the following order: alpha_t, alpha_v, mdot_gas, L_s, R_disk,
    mdot_peb, st, kappa, M_s, opt_vis."""
    
    parlist = []
    with open('param.in','r') as paramfile:
        parlines = paramfile.readlines()
        for line in parlines[37:45]:
            parlist.append(float(line.split()[-1]))
        parlist.append(float(parlines[28].split()[-1]))
        if parlines[46].split()[-1] in ['yes','Yes','y']:
            parlist.append(True)
        else:
            parlist.append(False)
        
    return parlist
 
def insert_planet_dmp(bigdata,name):
    """Converts the asteroidal input information for our new planet to cartesian
    coordinates and then writes it into our new big.dmp file."""

    cwd = os.getcwd()
    os.chdir('../conversion/')
    
    bad_ext = ['*.dmp','*.tmp','*.aei','*.clo','*.out']        
    for j in bad_ext:
        call(['find',os.getcwd(),'-maxdepth','1','-type','f','-name',j,'-delete'])
    
    initlist = [')O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n',\
        ") Lines beginning with `)' are ignored.\n",\
        ')---------------------------------------------------------------------\n',\
        ' style (Cartesian, Asteroidal, Cometary) = Asteroidal\n',\
        ' epoch (in days) = 0\n',\
        ')---------------------------------------------------------------------\n']
    
    with open('big.in','w+') as bigfile:
        
        for i in initlist:
            bigfile.write(i)
        
        bigfile.write(' {0:11}m={1:.17E} r={2:.0f}.d0 d={3:.2f} x={4:.2e}\n'.format(name,*bigdata[0:4]))
        bigfile.write(' {0: .17E} {1: .17E} {2: .17E}\n'.format(*bigdata[4:7]))
        bigfile.write(' {0: .17E} {1: .17E} {2: .17E}\n'.format(*bigdata[7:]))
        bigfile.write('  0. 0. 0.\n')
        
    call(['./mercury'])
    with open('big.dmp','r') as bigdmp:
        bigdmplines = bigdmp.readlines()
    os.chdir(cwd)
    with open('big.dmp','a') as bigdmpnew:
        for line in bigdmplines[6:]:
            bigdmpnew.write(line)

def insert_planet_embryo(mp,ap,emb_id):
    """Inserts a new embryo into the integration at the ice-line"""
    
    ep = np.random.rayleigh(1e-2) #We draw eccentricities from a Rayleigh distr.
    ip = 0.5*ep
    rp = 1.0
    dp = 1.5
    xp = 1.5e-9
    
    #We insert our parameters into our data array
    bigdata = np.array([mp,rp,dp,xp,ap,ep,ip,0,0,0])
    
    #Finally we generate the names and write it all into our big.in file
    name = 'P{}'.format(emb_id)
    
    with open('big.in','a') as bigfile:
        bigfile.write(' {0:11}m={1:.17E} r={2:.0f}.d0 d={3:.2f} x={4:.2e}\n'.format(name,*bigdata[0:4]))
        bigfile.write(' {0: .17E} {1: .17E} {2: .17E}\n'.format(*bigdata[4:7]))
        bigfile.write(' {0: .17E} {1: .17E} {2: .17E}\n'.format(*bigdata[7:]))
        bigfile.write('  0. 0. 0.\n')
        
    #We then convert our aei coordinates to a cartesian frame using an external
    #mercury run for 0.1 days and write it into our new big.dmp file
    
    insert_planet_dmp(bigdata,name)
    
def safronov_number(mp,ms,ap):
    """Computes the Safronov number assuming a constant density of the planetary
    embryos."""
    #We assume a constant density
    mp  = mp*Mearth
    ms  = ms*Msun
    ap  = ap*AU
    rho = 1.5

    rp  = (3*mp/(4*np.pi*rho))**(1/3)
    
    saf = np.sqrt(mp*ap/(ms*rp))
    
    return saf

def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def detect_merger():
    """Finds which planet has grown through collisions with another embryo at
    which time and mass."""
    
    coll = 'was hit by'
    
    plist = []
    tlist = []
    mlist = []
    
    with open('info.out','r') as infofile:
        for line in infofile:
            if coll in line:
                cinfo = line.split()
                plist.append(cinfo[0])
                tlist.append(float(cinfo[6]))
    
    for i in range(len(plist)):
        pdata = np.loadtxt(plist[i]+'.aei',skiprows=4)
        t = pdata[:,0]
        timeid = find_nearest(t,tlist[i])-1
        mlist.append(pdata[timeid,7])
    
    plistn = [int(i.strip('P')) for i in plist]
    
    return np.asarray([plistn,np.asarray(mlist)*msuntome,np.asarray(tlist)/365.25])

def check_isomass(a,mp,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis):
    """Goes through the semi-major axis values and mass values at every point in
    time and compares the mass to the isolation mass at the corresponding a value.
    If we see a sign change in the difference, we note the index where it occured
    and return it."""
    
    m_iso = iso_mass(a,mdot_gas,L_s,M_s,alpha_d,alpha_v,kap,opt_vis)
    
    rdsign = np.sign(m_iso-mp)
    
    if np.any(rdsign==-1):
        sz = rdsign == 0
        while sz.any():
            rdsign[sz] = np.roll(rdsign, 1)[sz]
            sz = rdsign == 0
        schange = ((np.roll(rdsign, 1) - rdsign) != 0).astype(int)
        #We set the first element to zero due to the circular shift of
        #numpy.roll
        schange[0] = 0
        #Finally we check where the array is equal to one and extract the
        #corresponding index
        scidx = np.where(schange)[0][0]
        return scidx
    else:
        return None 