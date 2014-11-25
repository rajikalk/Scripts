#find the velocity and mach number ofthe particles, and local sound speed.

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000.     #length units in cm
tu = 31558000.           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256.             #cm in grid length
G = 6.674e-8            #gravitational constant in cgs
Rsun = 6.96e10  #cm in a solar radius

def soundSpeed(gamma, pressure, density):
    ss = ((gamma*pressure)/(density))**0.5
    return ss

#Define arrays:
time = []
bv_1 = []
bv_2 = []
r_1 = []
r_2 = []
v_1 = []
v_2 = []

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/hot2/DD00*/CE00*.hierarchy")
f = open('velocity_contrast.csv','r+')
f.write("time, pv1, bv1, ratio1, pv2, bv2, ratio2 \n")

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
    sp1 = pf.h.sphere(pp1, 1/pf['rsun'])
    sp2 = pf.h.sphere(pp2, 1/pf['rsun'])
    bv1_vec = sp1.quantities['BulkVelocity']()
    bv2_vec = sp2.quantities['BulkVelocity']()
    bv1 = ((bv1_vec[0])**2. + (bv1_vec[1])**2. + (bv1_vec[2])**2.)**0.5
    bv2 = ((bv2_vec[0])**2. + (bv2_vec[1])**2. + (bv2_vec[2])**2.)**0.5
    bv_1.append(bv1)
    bv_2.append(bv2)
    
    #find velocity contrast
    r1 = v1/bv1
    r2 = v2/bv2
    
    r_1.append(r1)
    r_2.append(r2)
    
    print "For particle 1, v/bv:", r1
    print "For particle 2, v/bv:", r2
    
    f.write(str(time_val) + ',' + str(v1) + ',' + str(bv1) + ',' + str(r1) + ',' + str(v2) + ',' + str(bv2) + ',' + str(r2) + '\n')
    
plt.clf()
plt.figure(figsize=(8, 7))
top = plt.subplot(211)
#top.plot(time, m1, label='Mach Number')
top.plot(time, v_1, 'r:', label = 'speed of particle')
top.plot(time, bv_1, 'g--', label = 'bulk velocity of gas')
#top.set_title("Angular momentum loss per dump vs time")
#top.set_xlabel("Time (years)")
top.set_ylabel("Speed ($cms^{-1}$)")
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
topR = top.twinx()
topR.yaxis.tick_right()
topR.yaxis.set_label_position("right")
topR.plot(time, r_1, 'b', label='v/bv')
topR.set_ylabel("Ratio ($v/bv$)")
#topR.set_ylim([0.0, 8.0])
#top.legend(loc='upper right')
bot = plt.subplot(212)
#bot.plot(time, m2, label='Mach Number')
bot.plot(time, v_2, 'r:', label = 'speed of particle')
bot.plot(time, bv_2, 'g--', label = 'bulk velocity of gas')
#bot.set_title("Cumulative Angular Momentum loss")
bot.set_xlabel("Time (years)")
bot.set_ylabel("Speed ($cms^{-1}$)")
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
botR = bot.twinx()
botR.yaxis.tick_right()
botR.yaxis.set_label_position("right")
botR.plot(time, r_2, 'b', label='Mach Number')
botR.set_ylabel("Ratio ($v/bv$)")
#botR.set_ylim([0.0, 8.0])
plt.savefig('velocity_contrast.png')
