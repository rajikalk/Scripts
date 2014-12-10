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

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

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

#Define arrays:
time = []
rv_1 = []
rv_2 = []
coms = []	# Center of mass of system
it = 0

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
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
pf = load ("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")

dim = pf.domain_dimensions[0]
lu = pf['LengthUnits']     #length units in cm

for pf in ts:
	#get time
	t = pf.current_time
	time.append(t)

	#import all data
	dd = pf.h.all_data()
	g = pf.h.grids[0]

	#get particle positions
	pp1 = [dd['particle_position_x'][0], dd['particle_position_y'][0], dd['particle_position_z'][0]]
	pp2 = [dd['particle_position_x'][1], dd['particle_position_y'][1], dd['particle_position_z'][1]]

	#get particle velocites
	pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_velocity_z'][0]]
	pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_velocity_z'][1]]
	
	#get particle speed
	ps1 = (pv1[0]**2. + pv1[1]**2. + pv1[2]**2.)**0.5
	ps2 = (pv2[0]**2. + pv2[1]**2. + pv2[2]**2.)**0.5

	#make sphere around particles to get gas vector:
	sp1 = pf.h.sphere(pp1, 5.0/pf['rsun']) #sphere has radius of 2.5Rsun
	sp2 = pf.h.sphere(pp2, 5.0/pf['rsun'])

	#get gas vector
	bv1 = sp1.quantities["BulkVelocity"]()
	bv2 = sp2.quantities["BulkVelocity"]()
	
	#get mass of surrounding gas
	mg1 = sum(sp1['CellMass'])
	mg2 = sum(sp2['CellMass'])

	#get gas speed
	bs1 = (bv1[0]**2. + bv1[1]**2. + bv1[2]**2.)**0.5
	bs2 = (bv2[0]**2. + bv2[1]**2. + bv2[2]**2.)**0.5

	#find angle between particle and gas vector (angles are all in radians)
	ang1 = angle(pv1, bv1)
	ang2 = angle(pv2, bv2)

	#find projection of gas vector on the particle vector
	proj1 = bs1*math.cos(ang1)
	proj2 = bs2*math.cos(ang2)

	#find momentum of gas relative to particle
	rv1 = (proj1 - ps1)*mg1
	rv2 = (proj2 - ps2)*mg2
	rv_1.append(Lz1)
	rv_2.append(Lz2)

	#create plot
	plt.clf()
	plt.figure(figsize=(8, 6))
	plt.plot(time, rv_1, 'r', label='P1')
	plt.plot(time, rv_2, 'b', label='P2')
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
	plt.ylabel('momentum ($gcms^{-1}$)')
	plt.savefig('relative_momentum.png')
