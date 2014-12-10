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

#Define arrays:
time = []
rv_1 = []
rv_2 = []


#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD00*/CE00*.hierarchy")
pf = load ("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD0000/CE0000")

for pf in ts:
	#get time
	t = pf.current_time
	time.append(t)

	#import all data
	dd = pf.h.all_data()

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
	rv_1.append(rv1)
	rv_2.append(rv2)

#create plot
plt.clf()
plt.figure(figsize=(8, 6))
plt.plot(time, rv_1, 'r', label='P1')
plt.plot(time, rv_2, 'b', label='P2')
plt.axhline(y=0.0, color='k')
plt.xlabel('Time ($Years$)')
plt.ylabel('momentum ($gcms^{-1}$)')
plt.savefig('relative_momentum.png')
