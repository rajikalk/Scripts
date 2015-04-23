#!/usr/bin/env python
# Calculate the z-components of angular momentum for the simulation.

from yt.mods import *
import matplotlib.pyplot as plt
import csv


def linear_momentum(mass, velocity):
    # Calculates linear momentum p=mv
    # Everything should be in cgs
    # mass in grams
    # velocity in cm/s
    # return in gcm/s
    px = mass * velocity[0]
    py = mass * velocity[1]
    pz = mass * velocity[2]
    p = [px, py, pz]
    return p
    
def rel_position(pos1, pos2):
    # Finds the relative position of the cell position to the CoM.
    # Output units are the same as input units.
    # This is just to work out the sign of the linear momentum.
    x = pos1[0] - pos2[0]
    y = pos1[1] - pos2[1]
    z = pos1[2] - pos2[2]
    rel = [x, y, z]
    return rel

def z_momentum(rel_position, momentum):
    # Works out the z-component of the angular momentum
    # Using the cross product of radius and linear momentum
    # i.e., L = r x p
    # Everything in cgs
    Lz = rel_position[0] * momentum[1] - rel_position[1] * momentum[0]
    return Lz

# Load data:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0*/data0*.hierarchy")
pf = load("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0014/data0014")

#Define values:
dim = pf.domain_dimensions[0]
Msun = pf['MassUnits']       #grams/Msun
lu = pf['LengthUnits']     #length units in cm
tu = pf['TimeUnits']           #time units in sec
DensityUnits = pf['DensityUnits']
gl = lu/dim             #cm in grid length

#Define arrays:
time = []	# Time in simulation
L_tot = []	# Total L_z momentum in box
L_p = []	# Total L_z momentum of particles
L_g = []	# Total L_z momentum of gas
coms = []	# Center of mass of system

header = 0
with open('coms.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            x = float(row[7]) * lu	# Multiplies CoM by length of grid in cm
            y = float(row[8]) * lu	# This is because CoM is given in grid units
            z = float(row[9]) * lu
            com_val = [x,y,z]
            coms.append(com_val)
        else:
            header = 1

#Import all the timesteps for the series:
f = open('angular_momentum_conservation_z.csv','r+')
f.write("time, L_tot, L_p, L_g\n")

#Find the mass of the particles
dd = pf.h.all_data()
pm1 = dd["ParticleMass"][0]
pm2 = dd["ParticleMass"][1]

it = 0

for pf in ts:
    L_z = 0.
    Lp = 0.
    Lg = 0.
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #Find velocity of the particles
    pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_position_z'][0]]
    pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_position_z'][1]]
    
    com = coms[it]
    #com_err = com_errs[it]
    print com
    
    #Calculate angular momentum of the particles:
    p1 = linear_momentum(pm1, pv1)
    p2 = linear_momentum(pm2, pv2)
    rel1 = rel_position(pp1, com)
    rel2 = rel_position(pp2, com)
    L1 = z_momentum(rel1, p1)
    L2 = z_momentum(rel2, p2)
    Lp = L1 + L2
    L_p.append(Lp)

    #Get the grids in the data dump
    g = pf.h.grids[0]

    #Define counters

    for x in range(dim):
        for y in range(dim):
            for z in range(dim):
                pos = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMass"][x, y, z]
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = (z_momentum(rpos, p))#/mass
                L_z = L_z + L #Adds up all the z-compnents to find the total angular momentum of gas

    Lg = L_z
    L_g.append(Lg)
    
    Lt = Lg + Lp

    L_tot.append(Lt)
    it = it + 1
    
    f.write(str(time_val) + ',' + str(Lt) + ',' + str(Lp) + ',' + str(Lg) + '\n')
    print "The total amount of angular momentum is:", Lt, "(gcm^2)/s"
    print "The total amount of orbital angular momentum is:", Lp, "(gcm^2)/s"
    print "The total amount of angular momentum in gas:", Lg, "(gcm^2)/s"

