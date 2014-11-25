#COM position

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

def distance(pos1, pos2):
    x = (pos1[0] - pos2[0])**2.
    y = (pos1[1] - pos2[1])**2.
    z = (pos1[2] - pos2[2])**2.
    rel = (x+y+z)**0.5
    return rel

#Define constants and lists:
lu = 10000000000000.     #length units in cm
Rsun = 69600000000.      # solar radii in cm
center = [0.5*lu, 0.5*lu, 0.5*lu]
header = 0
time = []
com_t = []
com_p = []
com_g = []
dis_t = []
dis_p = []
dis_g = []

with open('coms.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            t = float(row[0])
            time.append(t)
            x_t = float(row[7])*lu
            y_t = float(row[8])*lu
            z_t = float(row[9])*lu
            com_total = [x_t,y_t,z_t]
            com_t.append(com_total)
            x_p = float(row[1])*lu
            y_p = float(row[2])*lu
            z_p = float(row[3])*lu
            com_particles = [x_p,y_p,z_p]
            com_p.append(com_particles)
            x_g = float(row[4])*lu
            y_g = float(row[5])*lu
            z_g = float(row[6])*lu
            com_gas = [x_g,y_g,z_g]
            com_g.append(com_gas)
        else:
            header = 1

for it in range(len(time)):
    com1 = com_t[it]
    dis1 = distance(com1, center)/Rsun
    dis_t.append(dis1)
    com2 = com_p[it]
    dis2 = distance(com2, center)/Rsun
    dis_p.append(dis2)
    com3 = com_g[it]
    dis3 = distance(com3, center)/Rsun
    dis_g.append(dis3)

fig = plt.figure(1)
fig.clf()
#fig.set_size_inches(7,7)
#fig.tight_layout()
#fig.yscale("log")
plt.plot(time, dis_t, label='Total CoM')
plt.plot(time, dis_p, label='Particle CoM')
plt.plot(time, dis_g, label='Gas CoM')
plt.axvline(x=0.0099951271, color='k')
plt.axvline(x=0.0399955081, color='k')
plt.axvline(x=0.0999962685, color='k')
plt.axvline(x=0.0599950006, color='b')
plt.legend(loc = 'lower right')
plt.xlabel("Time ($years$)")
plt.ylabel("Distance from center of grid ($R_\odot$)")
#plt.xlim([0.0, 0.1])
#plt.ylim([0.0, 0.5])
plt.savefig("com.png")
