#A Simple program that find the total mass in the grid, and returns data on the particles, such as their mass, separation and position.

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

#Load DataCube
pf = load("/media/DATA/Simulations/run1.e-6/DD0001/CE0001")

#Find Total Mass in the grid:
#Create region that covers entire grid:
reg = pf.h.region(pf.domain_center, pf.domain_left_edge, pf.domain_right_edge)

#Find domain size:
size = pf.domain_width* pf['au']

#Find total mass in region:
TotalMass = reg.quantities["TotalQuantity"]("TotalMassMsun")

#Find particle information:
dd = pf.h.all_data()
particles = dd["particle_type"] == 11
separation = ((((dd['particle_position_x'][particles][0]-dd['particle_position_x'][particles][1])**2) + ((dd['particle_position_y'][particles][0]-dd['particle_position_y'][particles][1])**2) + ((dd['particle_position_z'][particles][0]-dd['particle_position_z'][particles][1])**2))**(0.5))*pf['rsun']

print "+++++++++++++++++++++++++++++++++"
print "domain size in AU =", size
print "Current Time =",  pf.current_time, "Years"
print "Total Mass in the grid is =", TotalMass, "Msun"
print "Particle information:"
print "Primary:"
print "mass =", dd["ParticleMassMsun"][particles][0], "Msun"
print "position = (", dd['particle_position_x'][particles][0], ",", dd['particle_position_y'][particles][0], ",", dd['particle_position_z'][particles][0], ")"
print ""
print "Companion:"
print "mass =", dd["ParticleMassMsun"][particles][1], "Msun"
print "position = (", dd['particle_position_x'][particles][1], ",", dd['particle_position_y'][particles][1], ",", dd['particle_position_z'][particles][1], ")"
print " "
print "Separation =", separation, "rsun"
print "+++++++++++++++++++++++++++++++++"


