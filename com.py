#!/usr/bin/env python
#COM position

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

#Define arrays:
time = []
header = 0

#Import all the timesteps for the series:

ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0*/data0*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0014/data0014")
f = open('coms.csv','r+')
f.write("time, particles:x, y, z, gas:x, y, z, total:x, y, z\n")

#Define values:
dim = init_pf.domain_dimensions[0]
Msun = init_pf['MassUnits']       #grams/Msun
lu = init_pf['LengthUnits']     #length units in cm
tu = init_pf['TimeUnits']           #time units in sec
DensityUnits = init_pf['DensityUnits']
gl = lu/dim             #cm in grid length
Rsun = 69600000000.      # solar radii in  cm

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMass"][0]
pm2 = init_dd["ParticleMass"][1]
pm = pm1 + pm2

for pf in ts:
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()

    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]

    #set counters
    TM = pm
    x_top = pm1*pp1[0] + pm2*pp2[0]
    y_top = pm1*pp1[1] + pm2*pp2[1]
    z_top = pm1*pp1[2] + pm2*pp2[2]
    
    com_p = [(x_top/TM)/lu, (y_top/TM)/lu, (z_top/TM)/lu]
    
    x_top_g = 0
    y_top_g = 0
    z_top_g = 0
    TMg = 0
    
    #Get the grids in the data dump
    g = pf.h.grids[0]
    for x in range(dim):
        for y in range(dim):
            for z in range(dim):
                pos = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMass"][x, y, z]
                TM = TM + mass
                x_top = x_top + mass*pos[0]
                y_top = y_top + mass*pos[1]
                z_top = z_top + mass*pos[2]
                TMg = TMg + mass
                x_top_g = x_top_g + mass*pos[0]
                y_top_g = y_top_g + mass*pos[1]
                z_top_g = z_top_g + mass*pos[2]       
                
    x_pos = (x_top/TM)/lu
    y_pos = (y_top/TM)/lu
    z_pos = (z_top/TM)/lu
    com_t = [x_pos, y_pos, z_pos]
    
    x_pos_g = (x_top_g/TMg)/lu
    y_pos_g = (y_top_g/TMg)/lu
    z_pos_g = (z_top_g/TMg)/lu
    com_g = [x_pos_g, y_pos_g, z_pos_g]
    
    print 'com position:', com_t
                
    f.write(str(time_val) + "," + str(com_p[0]) + "," + str(com_p[1]) + "," + str(com_p[2]) + "," + str(com_g[0]) + "," + str(com_g[1]) + "," + str(com_g[2]) + "," + str(com_t[0]) + "," + str(com_t[1]) + "," + str(com_t[2]) + '\n')

