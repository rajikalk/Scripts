#Momentum conservation wrt to center

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

def linear_momentum(mass, velocity):
    #Everything should be in cgs
    #mass in grams
    #velocity in cm/s
    #return in gcm/s
    px = mass*velocity[0]
    py = mass*velocity[1]
    pz = mass*velocity[2]
    p = [px, py, pz]
    return p
    
def angular_momentum(rel_position, momentum):
    #Everything in cgs
    Lx = rel_position[1]*momentum[2] - rel_position[2]*momentum[1]
    Ly = rel_position[2]*momentum[0] - rel_position[0]*momentum[2]
    Lz = rel_position[0]*momentum[1] - rel_position[1]*momentum[0]
    L = [Lx, Ly, Lz]
    return L

def center_of_mass(mass1, mass2, position1, position2):
    #Everything in cgs
    x_top = mass1*position1[0] + mass2*position2[0]
    y_top = mass1*position1[1] + mass2*position2[1]
    z_top = mass1*position1[2] + mass2*position2[2]
    bot = mass1 + mass2
    x = x_top/bot
    y = y_top/bot
    z = z_top/bot
    com = [x,y,z]
    return com

def rel_position(pos1, pos2):
    x = pos1[0] - pos2[0]
    y = pos1[1] - pos2[1]
    z = pos1[2] - pos2[2]
    rel = [x,y,z]
    return rel
    
def z_momentum(rel_position, momentum):
    #Everything in cgs
    Lz = rel_position[0]*momentum[1] - rel_position[1]*momentum[0]
    return Lz

def L_mag(L):
    L = (L[0]**2. + L[1]**2. + L[2]**2.)**0.5
    return L

#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000     #length units in cm
tu = 31558000           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256             #cm in grid length

#Define arrays:
time = []
com_p = []
com_g = []
com_t = []
vel_p = []
vel_g = []
vel_t = []
pm = 0
mass_g = []
mass_t = []
Ltot = []
Lp = []
Lg = []
header = 0

with open('coms.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            time_val = float(row[0])
            time.append(time_val)
            x = float(row[1])*lu
            y = float(row[2])*lu
            z = float(row[3])*lu
            p = [x,y,z]
            com_p.append(p)
            x = float(row[4])*lu
            y = float(row[5])*lu
            z = float(row[6])*lu
            g = [x, y, z]
            com_g.append(g)
            x = float(row[7])*lu
            y = float(row[8])*lu
            z = float(row[9])*lu
            t = [x, y, z]
            com_t.append(t)
            v_x = float(row[10])
            v_y = float(row[11])
            v_z = float(row[12])
            v_p = [v_x, v_y, v_z]
            vel_p.append(v_p)
            v_x = float(row[13])
            v_y = float(row[14])
            v_z = float(row[15])
            v_g = [v_x, v_y, v_z]
            vel_g.append(v_g)
            v_x = float(row[16])
            v_y = float(row[17])
            v_z = float(row[18])
            v_t = [v_x, v_y, v_z]
            vel_t.append(v_t)
            pm = float(row[19])
            m_g = float(row[20])
            mass_g.append(m_g)
            m_t = float(row[21])
            mass_t.append(m_t)
        else:
            header = 1

#Import all the timesteps for the series:
f = open('angular_momentum_conservation_center.txt','r+')
f.write("time, L_tot, L_p, L_g \n")

#Save directory:
save_dir = "~/YT_output/run1/e-6lessetot/longer/"

center = [0.5*lu, 0.5*lu, 0.5*lu]

for it in range(len(time)):
    #For the Particle CoM
    com = com_p[it]
    p = linear_momentum(pm, vel_p[it])
    rel = rel_position(com, center)
    L_p = angular_momentum(rel, p)
    L_p_mag = L_mag(L_p)
    Lp.append(L_p_mag)
    
    #For the Particle CoM
    com = com_g[it]
    p = linear_momentum(mass_g[it], vel_g[it])
    rel = rel_position(com, center)
    L_g = angular_momentum(rel, p)
    L_g_mag = L_mag(L_g)
    Lg.append(L_g_mag)
    
    #For the Particle CoM
    com = com_t[it]
    p = linear_momentum(mass_t[it], vel_t[it])
    rel = rel_position(com, center)
    L_t = angular_momentum(rel, p)
    L_t_mag = L_mag(L_t)
    Ltot.append(L_t_mag)
    
    f.write(str(time_val) + "," + str(L_t_mag) + ',' + str(L_p_mag) + ',' + str(L_g_mag) + '\n')
    print 'time:', time[it]

plt.clf()
#fig.set_size_inches(7,7)
#fig.tight_layout()
#fig.yscale("log")
plt.plot(time, Ltot, label='Total CoM')
plt.plot(time, Lp, label='Particle CoM')
plt.plot(time, Lg, label='Gas CoM')
plt.legend()
plt.xlabel('time ($years$)')
plt.ylabel('Angular momentum ($gcm^2/s$)')
plt.savefig("CoM_angular_momentum.png")

