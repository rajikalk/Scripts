# find relative velocity using angle between vectors

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv
import math

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def projection(v1, v2):
  return (dotproduct(v1, v2)) / (length(v2))

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

def angular_momentum(rel_position, momentum):
    #Everything in cgs
    Lx = rel_position[1]*momentum[2] - rel_position[2]*momentum[1]
    Ly = rel_position[2]*momentum[0] - rel_position[0]*momentum[2]
    Lz = rel_position[0]*momentum[1] - rel_position[1]*momentum[0]
    L = [Lx, Ly, Lz]
    return L

#Define arrays:
time = []
rL_1 = []
rL_2 = []
coms = []	# Center of mass of system

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
pf = load ("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")

#Define values:
lu = pf['LengthUnits']     #length units in cm

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

#Find the mass of the particles
dd = pf.h.all_data()
pm1 = dd["ParticleMass"][0]
pm2 = dd["ParticleMass"][1]

it = 0

for pf in ts:
	#get time
	t = pf.current_time
	time.append(t)

	#import all data
	dd = pf.h.all_data()
	g = pf.h.grids[0]

	#set center of mass
	com = coms[it]

	#get particle positions
	pp1 = [dd['particle_position_x'][0], dd['particle_position_y'][0], dd['particle_position_z'][0]]
	pp2 = [dd['particle_position_x'][1], dd['particle_position_y'][1], dd['particle_position_z'][1]]

	#get particle velocites
	pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_velocity_z'][0]]
	pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_velocity_z'][1]]
	
	#Calculate angular momentum of the particles
	p1 = linear_momentum(pm1, pv1)
    	p2 = linear_momentum(pm2, pv2)
    	rel1 = rel_position(pp1, com)
    	rel2 = rel_position(pp2, com)
    	Lp1 = angular_momentum(rel1, p1)
    	Lp2 = angular_momentum(rel2, p2)
	Lpm1 = (Lp1[0]**2. + Lp1[1]**2. + Lp1[2]**2.)**0.5
	Lpm2 = (Lp2[0]**2. + Lp2[1]**2. + Lp2[2]**2.)**0.5

	#make sphere around particles to get gas vector:
	sp1 = pf.h.sphere(pp1, 5.0/pf['rsun']) #sphere has radius of 2.5Rsun
	sp2 = pf.h.sphere(pp2, 5.0/pf['rsun'])

	#get gas vector
	bv1 = sp1.quantities["BulkVelocity"]()
	bv2 = sp2.quantities["BulkVelocity"]()
	
	#get mass of surrounding gas
	gm1 = sum(sp1['CellMass'])
	gm2 = sum(sp2['CellMass'])

	#find momentum of the gas
	pg1 = linear_momentum(gm1, bv1)
    	pg2 = linear_momentum(gm2, bv2)
	Lg1 = angular_momentum(rel1, pg1)
    	Lg2 = angular_momentum(rel2, pg2)
	Lgm1 = (Lg1[0]**2. + Lg1[1]**2. + Lg1[2]**2.)**0.5
	Lgm1 = (Lg2[0]**2. + Lg2[1]**2. + Lg2[2]**2.)**0.5

	#find projection of gas vector on the particle vector
	Lgm1 = projection(Lp1, Lg1)
	Lgm2 = projection(Lp1, Lg1)

	#find angular momentum of gas relative to particle
	rL1 = Lpm1 - Lgm1
	rL2 = Lpm2 - Lgm2
	rL_1.append(rL1)
	rL_2.append(rL2)

	#create plot
	plt.clf()
	plt.figure(figsize=(8, 6))
	plt.plot(time, rL_1, 'r', label='P1')
	plt.plot(time, rL_2, 'b', label='P2')
	plt.axhline(y=0.0, color='k')
	plt.axvline(x=0.01, color='k')
	plt.axvline(x=0.04, color='k')
	plt.axvline(x=0.10, color='k')
	plt.axvline(x=0.06, color='b', ls='--')
	plt.axvline(x=0.02, color='r', ls=':')
	plt.axvline(x=0.03, color='r', ls=':')
	plt.axvline(x=0.05, color='r', ls=':')
	plt.axvline(x=0.07, color='r', ls=':')
	plt.axvline(x=0.08, color='r', ls=':')
	plt.axvline(x=0.09, color='r', ls=':')
	plt.axvline(x=0.06, color='b', ls='--')
	plt.xlabel('Time ($Years$)')
	plt.ylabel('L$_{z,p}$-L$_{z,g}$ ($gcms^{-1}$)')
	plt.savefig('relative_total_angular_momentum.png')
	
	it = it + 1
