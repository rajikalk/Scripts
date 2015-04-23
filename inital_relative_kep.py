#Keplerian velocity calculation:

from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np
import csv
    
#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000     #length units in cm
tu = 31558000           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256             #cm in grid length

#Define arrays:
time = []
pv1 = []
pv2 = []
'''
coms = []
comx = []
comy = []
comz = []
it = 0
header = 0
t2 = 0
x2 = 0
y2 = 0
z2 = 0

with open('com.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            x1 = x2
            x2 = float(row[11])
            y1 = y2
            y2 = float(row[12])
            z1 = z2
            z2 = float(row[13])
            t1 = t2
            t2 = float(row[0])*tu
            com_val = [x2,y2,z2]
            coms.append(com_val)
            if it == 0:
                comx_val = 0
                comy_val = 0
                comz_val = 0
            else:
                comx_val = (x2 - x1)/(t2 - t1)
                comy_val = (y2 - y1)/(t2 - t1)
                comz_val = (z2 - z1)/(t2 - t1)
                it = it + 1
            print comx_val, comy_val, comz_val
            comx.append(comx_val)
            comy.append(comy_val)
            comz.append(comz_val)
            it = it + 1
        else:
            header = 1
'''
#Load data dump:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
#f = open('particle_velocity.txt','r+')
#f.write("time, pv1, pv2\n")
it = 0

for pf in ts:
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Get the grids in the data dump
    g = pf.h.grids[0]
    '''
    com = [comx[it], comy[it], comz[it]]
    comv = (com[0]**2. + com[1]**2. + com[2]**2.)**(0.5)
    '''
    p1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_position_z'][0]]
    p2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_position_z'][1]]
    p1v = (p1[0]**2. + p1[1]**2. + p1[2]**2.)**(0.5)
    p2v = (p2[0]**2. + p2[1]**2. + p2[2]**2.)**(0.5)
    
    pv1_val = ((p1[0] - p2[0])**2. + (p1[1] - p2[1])**2. + (p1[2] - p2[2])**2.)**(0.5)
    #pv2_val = ((p2[0] - com[0])**2. + (p2[1] - com[1])**2. + (p2[2] - com[2])**2.)**(0.5)
    #print comv
    print 'p1:', p1v, 'p2:', p2v
    #print 'pv1rel:', comv/p1v, 'pv2rel:', comv/p2v
    
    pv1.append(pv1_val)
    #pv2.append(pv2_val)
   
    it = it + 1

plt.clf()
#fig = plt.figure(figsize=(18.5, 10.5))
plt.subplot(111)
plt.plot(time, pv1)
#plt.plot(time, pv2, label = 'velocity of particle 2')
#plt.legend(loc='lower right')
#plt.semilogy()
plt.title("Velocitiy of particles relative to the other particle vs time")
plt.xlabel("Time (years)")
plt.ylabel("Velocity ($cms^{-1}$)")
plt.savefig('velocity_timeseries_other_particle.png')

