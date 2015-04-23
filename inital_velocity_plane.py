#relative keplering velocity on plane

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

def distance(point1, point2):
#Takes in the position as grid points and give the separation in grid points
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result
    
def KeplerianVelocity(mass, radius):
#All units are in CGS
    v = ((6.67259e-8*mass)/radius)**0.5
    return v
    
#load data cube:
pf = load ("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr/DD0000/CE0000")
g = pf.h.grids[0]

#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000     #length units in cm
tu = 31558000           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256             #cm in grid length

#Set x and y position you wish to fins vertical velocity distribution for.
#x_pos = 128
y_pos = 128
z_pos = 0


#Find size of domain
dim = pf.domain_dimensions[0]

#Find the mass of the particles
dd = pf.h.all_data()
pm1 = dd["ParticleMassMsun"][0]*Msun
pm2 = dd["ParticleMassMsun"][1]*Msun
mass = pm1 + pm2

#define arrays:
#z_pos = []  #z position relative to orbital plane
x_pos = []  #x position along orbital plane
kep_v= []   #keplerian velocity at that height above adn below orbital plane
vel = []    #Initial velocity of gas above and below the orbital plane
rel_v = []  #Relative velocity of gas compared to keplerian velocity.

for x in range(dim):
    x_pos.append(x-128) #append z position. Take away 128 to have it relative to orbital plane.
    radius = distance([x*gl,y_pos*gl,z_pos*gl],pf.domain_center*pf["cm"])
    v_kep = KeplerianVelocity(mass, radius)
    kep_v.append(v_kep)
    v_x = g["x-velocity"][x,y_pos,z_pos]
    v_y = g["y-velocity"][x,y_pos,z_pos]
    v_z = g["z-velocity"][x,y_pos,z_pos]
    v_mag = (v_x**2 + v_y**2 + v_z**2)**0.5
    vel.append(v_mag)
    v_rel = v_mag/v_kep
    rel_v.append(v_rel)
    print "radius:", radius, "cm, mass:", mass, "g, v_kep:", v_kep, "cm/s, v:", v_mag
    
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(7,7)
#fig.yscale("log")
top = fig.add_subplot(311)
top.plot(x_pos, kep_v)
#top.set_ylim([0.0, 0.15])
#top.set_xlim([0.0, 0.1])
top.set_title('Calculated Keplerian (y = ' + str(y_pos) + ', z = ' + str(z_pos) + ')')
top.set_ylabel('Keplerian Velocity ($cms^{-1}$)')
middle = fig.add_subplot(312, sharex=top)
middle.plot(x_pos, vel)
#middle.set_ylim([0.0, 0.25])
middle.set_title('Initial velocity (y = ' + str(y_pos) + ', z = ' + str(z_pos) + ')')
middle.set_ylabel('Velocity ($cms^{-1}$)')
bot = fig.add_subplot(313, sharex=top)
bot.plot(x_pos, rel_v)
#bot.set_ylim([0.0, 0.25])
bot.set_title('Relative velocity (y = ' + str(y_pos) + ', z = ' + str(z_pos) + ')')
bot.set_xlabel('x position (grid units)')
bot.set_ylabel('Relative velocity ($v/v_k$)')
#plt.tight_layout()

filename = 'initial_velocity_y_'+str(y_pos)+'_z_'+str(z_pos)+'.png'
plt.savefig(filename)
