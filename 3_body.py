#psuedo 3 body sim

from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np
import csv

#Define gravitational acceleration funtion
def g_accel(Mass, Radius):
    if Radius == 0.0:
        result = 0
    else:
        top = 6.67259e-8*Mass
        bot = Radius**2.0
        result = top/bot
    return result
    
def rel_position(pos1, pos2):
    x = abs(pos1[0] - pos2[0])
    y = abs(pos1[1] - pos2[1])
    z = abs(pos1[2] - pos2[2])
    rel = [x,y,z]
    return rel

def velocity(v_i, a, t):
    v_x = v_i[0] + a[0]*t
    v_y = v_i[1] + a[1]*t
    v_z = v_i[2] + a[2]*t
    v = [v_x, v_y, v_z]
    return v

def position_change(p_i, v_i, t, a):
    p_x = p_i[0] + (v_i[0]*t) + ((0.5)*a[0]*(t**2.))
    p_y = p_i[1] + (v_i[1]*t) + ((0.5)*a[1]*(t**2.))
    p_z = p_i[2] + (v_i[2]*t) + ((0.5)*a[2]*(t**2.))
    p = [p_x, p_y, p_z]
    return p
    
#Open initial data dump:
pf = load('/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000')
dd = pf.h.all_data()

#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000.     #length units in cm
tu = 31558000.           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256.             #cm in grid length
G = 6.674e-8            #gravitational constant in cgs

#Initiate counters
it = 0
t = 0.
dt = 60.
dump_cyc = 315576. #0.01yrs in seconds
header = 0
pm1 = dd["ParticleMassMsun"][0]*Msun
pm2 = dd["ParticleMassMsun"][1]*Msun
pm3s = []
pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
pp3s = []
pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_velocity_z'][0]]
pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_velocity_z'][1]]

with open('separation.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            pm = float(row[4])*Msun
            pm3s.append(pm)
            x = float(row[5])*lu
            y = float(row[6])*lu
            z = float(row[7])*lu
            pos = [x, y, z]
            pp3s.append(pos)
        else:
            header = 1

#print pm3s
#print pp3s

#open write out file:
f = open('3_body.txt','r+')
f.write('time, pp1:x, y, z, pp2:x, y, z \n')
f.write(str(t) + ',' + str(pp1[0]) + ',' + str(pp1[1]) + ',' + str(pp1[2]) + ',' + str(pp2[0]) + ',' + str(pp2[1]) + ',' + str(pp2[2]) + '\n')
print 'time:', t/tu
print 'pp1:[', pp1[0]/lu, pp1[1]/lu, pp1[2]/lu, ']'
print 'pp2:[', pp2[0]/lu, pp2[1]/lu, pp2[2]/lu, ']'
print 'pv1:', pv1
print 'pv2:', pv2
print 'pm1', pm1
print 'pm2', pm2

while t < 50*dump_cyc:
    pp3 = pp3s[it]
    pm3 = pm3s[it]
    
    #Find relative positions:
    rp12 = rel_position(pp1, pp2)
    rp13 = rel_position(pp1, pp3)
    rp23 = rel_position(pp2, pp3)
    print 'rel_pos:[', rp12[0]/lu, rp12[1]/lu, rp12[2]/lu, ']'
    
    #Acceleration on particle 1:
    A1x1 = g_accel(pm2, rp12[0])
    A1x2 = g_accel(pm3, rp13[0])
    if pp1[0] > pp2[0]:
        A1x1 = (-1)*A1x1
    if pp1[0] > pp3[0]:
        A1x2 = (-1)*A1x2
    A1x = A1x1 + A1x2
    A1y1 = g_accel(pm2, rp12[1])
    A1y2 = g_accel(pm3, rp13[1])
    if pp1[1] > pp2[1]:
        A1y1 = (-1)*A1y1
    if pp1[1] > pp3[1]:
        A1y2 = (-1)*A1y2
    A1y = A1y1 + A1y2
    A1z1 = g_accel(pm2, rp12[2])
    A1z2 = g_accel(pm3, rp13[2])
    if pp1[2] > pp2[2]:
        A1z1 = (-1)*A1z1
    if pp1[2] > pp3[2]:
        A1z2 = (-1)*A1z2
    A1z = A1z1 + A1z2
    A1 = [A1x, A1y, A1z]
        
    #Acceleration on particle 2:
    A2x1 = g_accel(pm1, rp12[0])
    A2x2 = g_accel(pm3, rp23[0])
    if pp2[0] > pp1[0]:
        A2x1 = (-1)*A2x1
    if pp2[0] > pp3[0]:
        A2x2 = (-1)*A2x2
    A2x = A2x1 + A2x2
    A2y1 = g_accel(pm1, rp12[1])
    A2y2 = g_accel(pm3, rp23[1])
    if pp2[1] > pp1[1]:
        A2y1 = (-1)*A2y1
    if pp2[1] > pp3[1]:
        A2y2 = (-1)*A2y2
    A2y = A2y1 + A2y2
    A2z1 = g_accel(pm1, rp12[2])
    A2z2 = g_accel(pm3, rp23[2])
    if pp2[2] > pp1[2]:
        A2z1 = (-1)*A2z1
    if pp2[2] > pp3[2]:
        A2z2 = (-1)*A2z2
    A2z = A2z1 + A2z2
    A2 = [A2x, A2y, A2z]
    
    if pp1[0] > lu or pp1[0] < 0 or pp1[1] > lu or pp1[1] < 0 or pp1[2] > lu or pp1[2] < 0:
        t =50*dump_cyc
    
    if pp2[0] > lu or pp2[0] < 0 or pp2[1] > lu or pp2[1] < 0 or pp2[2] > lu or pp2[2] < 0:
        t =50*dump_cyc
    
    
    print 'time:', t
    print 'pp1:[', pp1[0]/lu, pp1[1]/lu, pp1[2]/lu, '], pp2:[', pp2[0]/lu, pp2[1]/lu, pp2[2]/lu, ']'
    print 'pv1:', pv1, ', pv2:', pv2
    print 'A1:', A1, ', A2:', A2
    
    
    #Calculate new particle positions:
    pp1 = position_change(pp1, pv1, dt, A1)
    pp2 = position_change(pp2, pv2, dt, A2)
    
    #Calculate new particle velocities:
    pv1 = velocity(pv1, A1, dt)
    pv2 = velocity(pv2, A2, dt)
    
    #print "A1:", A1, ", A2:", A2
    
    t = t + dt
    if t%dump_cyc == 0:
        it = it + 1
        print 'time:', t/tu
        print 'pp1:', pp1
        print 'pp2:', pp2
        print 'pv1:', pv1
        print 'pv2:', pv2
        print 'A1:', A1
        print 'A2:', A2
    f.write(str(t/tu) + ',' + str(pp1[0]) + ',' + str(pp1[1]) + ',' + str(pp1[2]) + ',' + str(pp2[0]) + ',' + str(pp2[1]) + ',' + str(pp2[2]) + '\n')
    
        
    
    
    
    
    
    
    
    
    
