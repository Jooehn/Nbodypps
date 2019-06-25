#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 10:21:35 2019

@author: John Wimarsson

Some functions that handles MERCURY input files
"""

import numpy as np
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

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
    
    with fdopen(fh,'w') as new_file:
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
    remove('param.in')
    move(abs_path, 'param.in')
    
def extend_stop_time(T):
    
    """Small function that updates the stop time in param.dmp to allow for an
    extended integration in case we have no collisions. Updates the old time
    value by adding the value T."""
    
    #The string we want to change
    stime_str = '  stop time (days) =    '
    
    #Makes temporary file
    fh, abs_path = mkstemp()
    
    with fdopen(fh,'w') as new_file:
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
    remove('param.dmp')
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