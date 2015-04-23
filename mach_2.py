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
grids_1 = []
grids_2 = []
mach_1 = []
mach_2 = []
ss_1 = []
ss_2 = []
v_1 = []
v_2 = []
pv_1 = []
pv_2 = []
den_1 = []
den_2 = []

#Import all the timesteps for the series:
'''
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot3/DD00*/CE00*.hierarchy")
pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot3/DD0000/CE0000")
'''
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/reggie/Hot_fb_0.4k/DD00*/CE00*.hierarchy")
pf = load ("/disks/ceres/makemake/acomp/reggie/Hot_fb_0.4k/DD0000/CE0000")

'''f = open('mach_corrected.txt','r+')
f.write("time, Mach:1, 2, ss:1, 2, v:1, 2 \n")'''

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
    pp1 = [int(round(dd['particle_position_x'][0]*dim)), int(round(dd['particle_position_y'][0]*dim)), int(round(dd['particle_position_z'][0]*dim))]
    pp2 = [int(round(dd['particle_position_x'][1]*dim)), int(round(dd['particle_position_y'][1]*dim)), int(round(dd['particle_position_z'][1]*dim))]
    
    #Find velocity of the particles
    pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_velocity_z'][0]]
    pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_velocity_z'][1]]
    ps1 = ((pv1[0])**2. + (pv1[1])**2. + (pv1[2])**2.)**0.5
    ps2 = ((pv2[0])**2. + (pv2[1])**2. + (pv2[2])**2.)**0.5
    pv_1.append(ps1)
    pv_2.append(ps2)
    #find center of box:
    #for particle 1:
    if pv1[0] > 0.0:
        cc1x = pp1[0] + 2
        cc1xn = pp1[0] + 1
    elif pv1[0] < 0.0:
        cc1x = pp1[0] - 2
        cc1xn = pp1[0] - 1
    else:
        cc1x = pp1[0]
        cc1xn = pp1[0]
    if pv1[1] > 0.0:
        cc1y = pp1[1] + 2
        cc1yn = pp1[1] + 1
    elif pv1[1] < 0.0:
        cc1y = pp1[1] - 2
        cc1yn = pp1[1] - 1
    else:
        cc1y = pp1[1]
        cc1yn = pp1[1]
    if pv1[2] > 0.0:
        cc1z = pp1[2] + 2
        cc1zn = pp1[2] + 1
    elif pv1[2] < 0.0:
        cc1z = pp1[2] - 2
        cc1zn = pp1[2] - 1
    else:
        cc1z = pp1[2]
        cc1zn = pp1[2]
    cc1 = [cc1x, cc1y, cc1z]
    cc1n = [cc1xn, cc1yn, cc1zn]
    #for particle 2:
    if pv2[0] > 0.0:
        cc2x = pp2[0] + 2
        cc2xn = pp2[0] + 1
    elif pv2[0] < 0.0:
        cc2x = pp2[0] - 2
        cc2xn = pp2[0] - 1
    else:
        cc2x = pp2[0]
        cc2xn = pp2[0]
    if pv2[1] > 0.0:
        cc2y = pp2[1] + 2
        cc2yn = pp2[0] + 1
    elif pv2[1] < 0.0:
        cc2y = pp2[1] - 2
        cc2yn = pp2[0] - 1
    else:
        cc2y = pp2[1]
        cc2yn = pp2[0]
    if pv2[2] > 0.0:
        cc2z = pp2[2] + 2
        cc2zn = pp2[2] + 1
    elif pv2[2] < 0.0:
        cc2z = pp2[2] - 2
        cc2zn = pp2[2] - 1
    else:
        cc2z = pp2[2]
        cc2zn = pp2[2]
    cc2 = [cc2x, cc2y, cc2z]
    cc2n = [cc2xn, cc2yn, cc2zn]
    
    #Get the relavent grids, and get the sound speed, and get bulk velocitys:
    ssg_1 = []
    ssg_2 = []
    dg_1 = []
    dg_2 = []
    bvx_1 = []
    bvy_1 = []
    bvz_1 = []
    bvx_2 = []
    bvy_2 = []
    bvz_2 = []
    #for particle 1:
    for x in range(cc1[0]-1, cc1[0]+2):
        for y in range(cc1[1]-1, cc1[1]+2):
            for z in range(cc1[2]-1, cc1[2]+2):
                if [x,y,z] != cc1n:
                    ssg1 = g['SoundSpeed'][x,y,z]
                    ssg_1.append(ssg1)
		    dg1 = g['Density'][x,y,z]
		    dg_1.append(dg1)
                    bvx = g['x-velocity'][x,y,z]
                    bvx_1.append(bvx)
                    bvy = g['y-velocity'][x,y,z]
                    bvy_1.append(bvy)
                    bvz = g['z-velocity'][x,y,z]
                    bvz_1.append(bvz)
    ss1 = sum(ssg_1)/len(ssg_1)
    ss_1.append(ss1)
    den1 = sum(dg_1)/len(dg_1)
    den_1.append(den1)
    bvx1 = sum(bvx_1)/len(bvx_1)
    bvy1 = sum(bvy_1)/len(bvy_1)
    bvz1 = sum(bvz_1)/len(bvz_1)
    bv1 = [bvx1, bvy1, bvz1]
    #for particle 2:
    for x in range(cc2[0]-1, cc2[0]+2):
        for y in range(cc2[1]-1, cc2[1]+2):
            for z in range(cc2[2]-1, cc2[2]+2):
                if [x,y,z] != cc2n:
                    ssg2 = g['SoundSpeed'][x,y,z]
                    ssg_2.append(ssg2)
		    dg2 = g['Density'][x,y,z]
		    dg_2.append(dg2)
                    bvx = g['x-velocity'][x,y,z]
                    bvx_2.append(bvx)
                    bvy = g['y-velocity'][x,y,z]
                    bvy_2.append(bvy)
                    bvz = g['z-velocity'][x,y,z]
                    bvz_2.append(bvz)
    ss2 = sum(ssg_2)/len(ssg_2)
    ss_2.append(ss2)
    den2 = sum(dg_2)/len(dg_2)
    den_2.append(den2)
    bvx2 = sum(bvx_2)/len(bvx_2)
    bvy2 = sum(bvy_2)/len(bvy_2)
    bvz2 = sum(bvz_2)/len(bvz_2)
    bv2 = [bvx2, bvy2, bvz2]
    
    #Correcting velocity of particles for bulk velocity of the gas
    if np.sign(pv1[0])==np.sign(bv1[0]) and np.sign(pv1[1])==np.sign(bv1[1]):
	if abs(bv1[0])>abs(pv1[0]) or abs(bv1[1])>abs(pv1[1]):
    	    v1 = -1*(((pv1[0]-bv1[0])**2. + (pv1[1]-bv1[1])**2. + (pv1[2]-bv1[2])**2.)**0.5)
    else:
	v1 = ((pv1[0]-bv1[0])**2. + (pv1[1]-bv1[1])**2. + (pv1[2]-bv1[2])**2.)**0.5
    if np.sign(pv2[0])==np.sign(bv2[0]) and np.sign(pv2[1])==np.sign(bv2[1]):
	if abs(bv2[0])>abs(pv2[0]) or abs(bv2[1])>abs(pv2[1]):
    	    v2 = -1*(((pv2[0]-bv2[0])**2. + (pv2[1]-bv2[1])**2. + (pv2[2]-bv2[2])**2.)**0.5)
    else:
	v2 = ((pv2[0]-bv2[0])**2. + (pv2[1]-bv2[1])**2. + (pv2[2]-bv2[2])**2.)**0.5
    v_1.append(v1)
    v_2.append(v2)
    
    m1 = ps1/ss1
    m2 = ps2/ss2
    mach_1.append(m1)
    mach_2.append(m2)
    
    print "For particle 1, v:", v1, "ss:", ss1, "Mach No.:", m1
    print "For particle 2, v:", v2, "ss:", ss2, "Mach No.:", m2
    '''f.write(str(time_val) + ',' + str(m1) + ',' + str(m2) + ',' + str(ss1) + ',' + str(ss2) + ',' + str(v1) + ',' + str(v2) + '\n')'''

plt.clf()
plt.figure(figsize=(4, 3.5))
top = plt.subplot(311)
top.plot(time, mach_1, 'r', label='0.6$M_\odot$')
top.plot(time, mach_2, 'b--', label='0.4$M_\odot$')
top.set_ylabel("Mach number")
mid = plt.subplot(312)
mid.plot(time, v_1, 'r', label='0.6$M_\odot$')
mid.plot(time, v_2, 'b--', label='0.4$M_\odot$')
mid.set_ylabel("Relative velocity ($v_p - v_g$)(cms^{-1})")
bot = plt.subplot(313)
bot.plot(time, den_1, 'r', label='0.6$M_\odot$')
bot.plot(time, den_2, 'b--', label = '0.4$M_\odot$')
bot.set_xlabel("Time (years)")
bot.set_ylabel("encountered density (gcm^{-3})")
plt.savefig('relative_velocity_2.png')
