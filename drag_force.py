# Calculate drag force of the particles

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
r1 = 0.5*Rsun
r2 = 0.01*Rsun

def curly_F(Mass, Velocity, Density):
    top = 4.0*np.pi*((G*Mass)**2.0)*Density
    bot = Velocity**2.0
    result = top/bot
    return result

def I_sub(Mach):
    result = 0.5*np.log((1+Mach)/(1-Mach)) - Mach
    return result
    
def I_super(Mach, r_min, cs, t):
    result = 0.5*np.log((Mach+1)/(Mach-1)) + np.log((Mach-1)/(r_min/(cs*t)))    
    return result
    
#Define arrays:
time = []
mach_1 = []
mach_2 = []
ss_1 = []
ss_2 = []
v_1 = []
v_2 = []
DF_1 = []
DF_2 = []
header = 0
it = 0
t1 = 'off'
t2 = 'off'
dt1 = 0.
dt2 = 0.

with open('mach.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            m1 = float(row[1])
            mach_1.append(m1)
            m2 = float(row[2])
            mach_2.append(m2)
            ss1 = float(row[3])
            ss_1.append(ss1)
            ss2 = float(row[4])
            ss_2.append(ss2)
            v1 = float(row[5])
            v_1.append(v1)
            v2 = float(row[6])
            v_2.append(v2)
        else:
            header = 1
            
#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/media/DATA/Simulations/smallbox/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
f = open('drag_force.txt','r+')
f.write("time, drag:p1, p2 \n")

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun
r1 = 0.5*Rsun
r2 = 0.01*Rsun

for pf in ts:
    dt = 0.01*tu
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0], dd['particle_position_y'][0], dd['particle_position_z'][0]]
    pp2 = [dd['particle_position_x'][1], dd['particle_position_y'][1], dd['particle_position_z'][1]]
    
    v1 = v_1[it]
    v2 = v_2[it]

    #Make spheres around the particles
    sp1 = pf.h.sphere(pp1, 10/pf['rsun'])
    sp2 = pf.h.sphere(pp2, 10/pf['rsun'])
    ss1 = sp1['Density']
    ss2 = sp2['Density']
    d1 = (sum(ss1))/(len(ss1))
    d2 = (sum(ss2))/(len(ss2))
    
    #Calculate curly F:
    F1 = curly_F(pm1, v1, d1)
    F2 = curly_F(pm2, v2, d2)
    
    m1 = mach_1[it]
    m2 = mach_2[it]
    ss1 = ss_1[it]
    ss2 = ss_2[it]
    
    if m1 < 1:
        #subsonic case:
        I1 = I_sub(m1)
    else:
        if t1 == 'off':
            t1 = 'on'
            dt1 = dt
        else:
            dt1 = dt1 + dt
        I1 = I_super(m1, r1, ss1, dt1)
        
    if m2 < 1:
        #subsonic case:
        I2 = I_sub(m2)
    else:
        if t2 == 'off':
            t2 = 'on'
            dt2 = dt
        else:
            dt2 = dt2 + dt
        I2 = I_super(m2, r2, ss2, dt2)
    
    DF1 = -I1*F1
    DF2 = -I2*F2
    DF_1.append(DF1)
    DF_2.append(DF2)
    
    f.write(str(time_val) + ',' + str(DF1) + ',' + str(DF2) + '\n')

plt.clf()
plt.figure(figsize=(8, 7))
plt.plot(time, DF_1, label='Particle 1')
plt.plot(time, DF_2, label='Particle 2')
plt.axvline(x=0.01, color='k')
plt.axvline(x=0.04, color='k')
plt.axvline(x=0.1, color='k')
plt.axvline(x=0.06, color='b', ls='--')
plt.axvline(x=0.02, color='r', ls=':')
plt.axvline(x=0.03, color='r', ls=':')
plt.axvline(x=0.05, color='r', ls=':')
plt.axvline(x=0.07, color='r', ls=':')
plt.axvline(x=0.08, color='r', ls=':')
plt.axvline(x=0.09, color='r', ls=':')
plt.legend(loc='lower right')
plt.xlabel('Time ($years$)')
plt.ylabel('Drag Froce ($dyne$)')
#plt.ylim([-6.e35, 0.0])
plt.savefig('drag_force.png')
