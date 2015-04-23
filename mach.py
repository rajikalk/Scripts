#find the velocity and mach number ofthe particles, and local sound speed.

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

def soundSpeed(gamma, pressure, density):
    ss = ((gamma*pressure)/(density))**0.5
    return ss

#Define arrays:
time = []
r_1 = []
r_2 = []
v_1 = []
v_2 = []
den_1 = []
den_2 = []
mach_1 = []
mach_2 = []

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/reggie/Hot_fb_0.4k/DD00*/CE00*.hierarchy")
pf = load ("/disks/ceres/makemake/acomp/reggie/Hot_fb_0.4k/DD0000/CE0000")

#Define values:
dim = pf.domain_dimensions[0]
Msun = pf['MassUnits']       #grams/Msun
lu = pf['cm']     #length units in cm
tu = pf['TimeUnits']           #time units in sec
DensityUnits = pf['DensityUnits']
gl = lu/dim             #cm in grid length
G = 6.674e-8            #gravitational constant in cgs
Rsun = lu/pf['rsun']  #cm in a solar radius

for pf in ts:
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    g = pf.h.grids[0]
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0], dd['particle_position_y'][0], dd['particle_position_z'][0]]
    pp2 = [dd['particle_position_x'][1], dd['particle_position_y'][1], dd['particle_position_z'][1]]
    
    #Find velocity of the particles
    pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_velocity_z'][0]]
    pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_velocity_z'][1]]
    v1 = ((pv1[0])**2. + (pv1[1])**2. + (pv1[2])**2.)**0.5
    v2 = ((pv2[0])**2. + (pv2[1])**2. + (pv2[2])**2.)**0.5
    v_1.append(v1)
    v_2.append(v2)
    
    #Make spheres around the particles and find bulk velocity
    sp1 = pf.h.sphere(pp1, 2.5/pf['rsun'])
    sp2 = pf.h.sphere(pp2, 2.5/pf['rsun'])
    bv1 = sp1.quantities["BulkVelocity"]()
    bv2 = sp2.quantities["BulkVelocity"]()
    ss1 = sum(sp1['SoundSpeed'])/len(sp1['SoundSpeed'])
    ss2 = sum(sp2['SoundSpeed'])/len(sp2['SoundSpeed'])
    den1 = sum(sp1['Density'])/len(sp1['Density'])
    den2 = sum(sp2['Density'])/len(sp2['Density'])
    den_1.append(den1)
    den_2.append(den2)    

    #find mach number
    m1 = v1/ss1
    m2 = v2/ss2
    mach_1.append(m1)
    mach_2.append(m2)

    #find velocity contrast
    if np.sign(pv1[0])==np.sign(bv1[0]) and np.sign(pv1[1])==np.sign(bv1[1]):
	if abs(bv1[0])>abs(pv1[0]) or abs(bv1[1])>abs(pv1[1]):
    	    r1 = -1*(((pv1[0]-bv1[0])**2. + (pv1[1]-bv1[1])**2. + (pv1[2]-bv1[2])**2.)**0.5)
	else:
	    r1 = ((pv1[0]-bv1[0])**2. + (pv1[1]-bv1[1])**2. + (pv1[2]-bv1[2])**2.)**0.5
    else:
	r1 = ((pv1[0]-bv1[0])**2. + (pv1[1]-bv1[1])**2. + (pv1[2]-bv1[2])**2.)**0.5
    if np.sign(pv2[0])==np.sign(bv2[0]) and np.sign(pv2[1])==np.sign(bv2[1]):
	if abs(bv2[0])>abs(pv2[0]) or abs(bv2[1])>abs(pv2[1]):
    	    r2 = -1*(((pv2[0]-bv2[0])**2. + (pv2[1]-bv2[1])**2. + (pv2[2]-bv2[2])**2.)**0.5)
	else:
	    r2 = ((pv2[0]-bv2[0])**2. + (pv2[1]-bv2[1])**2. + (pv2[2]-bv2[2])**2.)**0.5
    else:
	r2 = ((pv2[0]-bv2[0])**2. + (pv2[1]-bv2[1])**2. + (pv2[2]-bv2[2])**2.)**0.5
    r_1.append(r1)
    r_2.append(r2)
    
    print "For particle 1, v/bv:", r1
    print "For particle 2, v/bv:", r2
    
    
plt.clf()
plt.figure(figsize=(8, 6))
top = plt.subplot(311)
top.plot(time, mach_1, 'r', label='0.6$M_\odot$')
top.plot(time, mach_2, 'b--', label='0.4$M_\odot$')
top.set_ylabel("Mach number")
mid = plt.subplot(312)
mid.plot(time, r_1, 'r', label='0.6$M_\odot$')
mid.plot(time, r_2, 'b--', label='0.4$M_\odot$')
plt.axhline(y=0.0, color='k')
mid.set_ylabel("Relative velocity ($v_p - v_g$)(cms^{-1})")
bot = plt.subplot(313)
bot.plot(time, den_1, 'r', label='0.6$M_\odot$')
bot.plot(time, den_2, 'b--', label = '0.4$M_\odot$')
bot.set_xlabel("Time (years)")
bot.set_ylabel("encountered density (gcm^{-3})")
plt.savefig('relative_velocity.png')
