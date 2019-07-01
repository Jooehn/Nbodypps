#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:21:35 2019

@author: John Wimarsson

Some functions that handles MERCURY input files
"""

import numpy as np
import os
from tempfile import mkstemp
from subprocess import call
from shutil import move

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

def check_cavity_planet(namelist):
    """Removes planet that has entered our cavity in order to keep the timestep
    from converging at low values and thereby slowing down our integration.
        names: The names of our planets
        data: The corresponding data input data for the planets
    Returns True if all planets have entered the cavity and the integration
    should be terminated."""
    
    #We generate output files
#    call(['./element'])
    names = namelist.copy()
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
            print('We remove P{} from the integration\n'.format(cavity_ids[i]+1))
            old_names.append(namelist[cavity_ids[i]])
            names.remove(namelist[cavity_ids[i]])
        #We loop through each line in big.dmp and remove the information of the
        #planets inside the cavity
        
        #Makes temporary file
        fh, abs_path = mkstemp()
    
        with os.fdopen(fh,'w') as new_file:
            with open('big.dmp') as bigfile:
                biglines = bigfile.readlines()
                i = 6
                while len(old_names)>0:
                    if old_names[0] in biglines[i]:
                        del biglines[i:i+4]
                        del old_names[0]
                    else:
                        i += 4
            new_file.writelines(biglines)
        #Remove original file and move new file
        for i in range(len(names)):
            os.remove(names[i]+'.aei')    
        os.remove('big.dmp')
        move(abs_path, 'big.dmp')    
        return False
    else:
        call(['find',os.getcwd(),'-maxdepth','1','-type','f','-name','*.aei','-delete'])
        return False
os.chdir(os.getcwd()+'/../sim/migtest3')
check_cavity_planet(bignames)