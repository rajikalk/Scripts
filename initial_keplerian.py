#Find ratio of gas that is super keplerian

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

def keplerian_vel(mass, radius):
#All units are in CGS
    v = ((6.67259e-8*mass)/radius)**0.5
    return v
    
def distance(point1, point2):
#Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result

#Import all the timesteps for the series:
pf = load("/media/DATA/Simulations/smallbox/hot3/DD0000/CE0000")

#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000.     #length units in cm
tu = 31558000.           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256.             #cm in grid length

#Define iterators:
sub_kep = 0
super_kep = 0

#find size of grid:
dim = pf.domain_dimensions[0]

#Find the mass of the particles
dd = pf.h.all_data()
pm1 = dd["ParticleMassMsun"][0]*Msun
pm2 = dd["ParticleMassMsun"][1]*Msun
mass_tot = pm1 + pm2
    
g = pf.h.grids[0]

for x in range(dim):
    for y in range(dim):
        for z in range(dim):
            position = [x*gl,y*gl,z*gl]
            center = pf.domain_center*pf["cm"]
            radius = distance(position,center)
            ke = keplerian_vel(mass_tot, radius)
            v_x = g["x-velocity"][x,y,z]
            v_y = g["y-velocity"][x,y,z]
            v_z = g["z-velocity"][x,y,z]
            velocity = (v_x**2 + v_y**2 + v_z**2)**0.5
            if velocity > ke:
                super_kep = super_kep + 1
            else:
                sub_kep = sub_kep + 1

super_percent = (super_kep*100.)/(dim**3.)
sub_percent = (sub_kep*100.)/(dim**3.)
print "Initially:"
print "The amount of super keplerian gas is", super_percent, "%"
print "The amount of sub keplerian gas is", sub_percent, "%"
            
            
