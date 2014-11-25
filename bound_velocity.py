#Calculates bounda dn unbound gas using velocities

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

#Define functions:
def escape_vel(mass, radius):
#All units are in CGS
    v = ((2.*6.67259e-8*mass)/radius)**0.5
    return v
    
def distance(point1, point2):
#Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result
    
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
    
#Define values:
Msun = 1.9891e+33        #grams/Msun
lu = 10000000000000.     #length units in cm
tu = 31558000.           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256.             #cm in grid length

#Define arrays:
time = []
bound = []
unbound = []
bm = []
ubm = []

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
f = open('velocity_bound_edge.txt','r+')
f.write('time, bound, unbound, bound_mass, unbound_mass\n')

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun
mass_tot = pm1 + pm2

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]
edge = [0.0, dim-1]

for pf in ts:
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    com = center_of_mass(pm1, pm2, pp1, pp2)

    #Get the grids in the data dump
    g = pf.h.grids[0]

    #Define counters
    it = 0
    bound_val = 0
    unbound_val = 0
    
    #Left and right edge:
    for x in edge:
        for y in range(dim):
            for z in range(dim):
                position = [x*gl,y*gl,z*gl]
                radius = distance(position,com)
                cellMass = g["CellMassMsun"][x,y,z]*Msun #in grams
                ve = escape_vel(mass_tot, radius)
                v_x = g["x-velocity"][x,y,z]
                v_y = g["y-velocity"][x,y,z]
                v_z = g["z-velocity"][x,y,z]
                velocity = (v_x**2 + v_y**2 + v_z**2)**0.5
                if velocity > ve:
                    unbound_val = unbound_val + cellMass
                    print "unbound"
                else:
                    bound_val = bound_val + cellMass
                    print "bound"
                it = it + 1
    
    #Front and back edge:
    it = 0
    for x in range(dim):
        for y in edge:
            for z in range(dim):
                position = [x*gl,y*gl,z*gl]
                radius = distance(position,com)
                cellMass = g["CellMassMsun"][x,y,z]*Msun #in grams
                ve = escape_vel(mass_tot, radius)
                v_x = g["x-velocity"][x,y,z]
                v_y = g["y-velocity"][x,y,z]
                v_z = g["z-velocity"][x,y,z]
                velocity = (v_x**2 + v_y**2 + v_z**2)**0.5
                if velocity > ve:
                    unbound_val = unbound_val + cellMass
                    print "unbound"
                else:
                    bound_val = bound_val + cellMass
                    print "bound"
                it = it + 1
    
    #Top and bottom edge:
    it = 0
    for x in range(dim):
        for y in range(dim):
            for z in edge:
                position = [x*gl,y*gl,z*gl]
                radius = distance(position,com)
                cellMass = g["CellMassMsun"][x,y,z]*Msun #in grams
                ve = escape_vel(mass_tot, radius)
                v_x = g["x-velocity"][x,y,z]
                v_y = g["y-velocity"][x,y,z]
                v_z = g["z-velocity"][x,y,z]
                velocity = (v_x**2 + v_y**2 + v_z**2)**0.5
                if velocity > ve:
                    unbound_val = unbound_val + cellMass
                    print "unbound"
                else:
                    bound_val = bound_val + cellMass
                    print "bound"
                it = it + 1
                
    TotalMass = unbound_val + bound_val

    unbound_percent = (unbound_val*100.)/(TotalMass)
    bound_percent = (bound_val*100.)/(TotalMass)
    unbound.append(unbound_percent)
    bound.append(bound_percent)
    bm.append(bound_val)
    ubm.append(unbound_val)
    
    f.write(str(time_val) + ',' + str(bound_percent) + ',' + str(unbound_percent) + ',' + str(bound_val) + ',' + str(unbound_val) + '\n')
    print "The amount of unbound gas at edge is", unbound_percent, "%"
    print "The amount of bound gas at edge is", bound_percent, "%"
        
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(7,7)
#fig.tight_layout()
#fig.yscale("log")
top = fig.add_subplot(111)
top.plot(time, unbound)
top.set_ylim([0.0, 100.0])
#top.set_xlim([0.0, 0.1])
top.set_title('Unbound gas (using velocities), at the edge')
top.set_xlabel('Time ($years$)')
top.set_ylabel('Unbound ratio (%)')
plt.savefig("velocity_edge_unbound.png")
