#Calculates fallback time of gas from edge of grid

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import pylab as P

#Defines function:
def g_accel(mass, radius):
#Takes mass and radius in cgs and spits out value in cgs
    g_top = (-6.67259e-8)*mass
    g_bot = (radius)**2.
    g = g_top/g_bot
    return g
    
def projectile_height(velocity, time_step, acceleration):
#Find the change in height for a projectile over a time step. All units given in CGS
    first_term = velocity*time_step
    second_term = 0.5*acceleration*(time_step**2.)
    h = first_term + second_term
    return h
    
def velocity(initial_velocity, acceleration, time_step):
#Calculates the velocity over a given time step. All units given in CGS.
    v = initial_velocity + acceleration*time_step
    return v

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

#load data:

ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
'''
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
'''
 
#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000.     #length units in cm
tu = 31558000.           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256.             #cm in grid length
dt = 60.*60*24.*7. #Time step
day = 60.*60.*24.*7. #seconds in a day
Rsun = 6.955e10 #cm in Rsun


#Save directory:
#save_dir = "~/YT_output/run1.e-6lessetot_Gcorr_0.75k/Fallback/"

#Find the mass of the particles, and the center of mass:
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun
mass = pm1 + pm2

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]
edge = [0.0, dim-1]
it = 0
region = 1.0*dim
lower = int(region)
upper = int(dim-region)
zedge = [lower, upper]

for pf in ts:
    #Get current time step:
    time_val = pf.current_time

    #Get the grids in the data dump:
    g = pf.h.grids[0]
    dd = pf.h.all_data()
    
    #Define Arrays:
    cmass = []
    time = []
    radius = []
    cmass_perp = []
    time_perp = []
    radius_perp = []
    
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    com = center_of_mass(pm1, pm2, pp1, pp2)

    #Left and right edge:
    print "For left and right edge"
    for x in edge:
        for y in range(dim):
            for z in range(dim):
                position = [x*gl,y*gl,z*gl]
		cell = g["CellMass"][x,y,z]
                r = distance(position, com)
                r_0 = r-1
                radius.append(r/Rsun)
                v_x = g["x-velocity"][x,y,z]
                v_y = g["y-velocity"][x,y,z]
                v_z = g["z-velocity"][x,y,z]
                v_rad = (v_x**2. + v_y**2. + v_z**2.)**0.5
                ve = escape_vel(mass, r)
                a = g_accel(mass, r)
                if v_rad < ve:
                    t = 0.
                    while r > r_0:
                        height = projectile_height(v_rad, dt, a)
                        r = r + height
                        v_rad = velocity(v_rad, a, dt)
                        a = g_accel(mass, r)
                        t = t + dt
                        #if t/day > 52.:
                            #r = 0.0
                    cell = g["CellMass"][x,y,z]
		    cmass.append(cell)
                    time.append(t/day)
                #print y, t/day, "weeks"
        
    #Front and back edge:
    print "For front and back edge"
    for x in range(dim):
        for y in edge:
            for z in range(dim):
                position = [x*gl,y*gl,z*gl]
                r = distance(position, com)
                r_0 = r-1
                radius.append(r/Rsun)
                v_x = g["x-velocity"][x,y,z]
                v_y = g["y-velocity"][x,y,z]
                v_z = g["z-velocity"][x,y,z]
                v_rad = (v_x**2. + v_y**2. + v_z**2.)**0.5
                ve = escape_vel(mass, r)
                a = g_accel(mass, r)
                if v_rad < ve:
                    t = 0.
                    while r > r_0:
                        height = projectile_height(v_rad, dt, a)
                        r = r + height
                        v_rad = velocity(v_rad, a, dt)
                        a = g_accel(mass, r)
                        t = t + dt
                        #if t/day > 52.:
                            #r = 0.0
                    cell = g["CellMass"][x,y,z]
		    cmass.append(cell)
                    time.append(t/day)
                #print x, t/day, "weeks"

    #Top and bottom edge:
    print "For top and bottom edge"
    for x in range(dim):
        for y in range(dim):
            for z in edge:
                position = [x*gl,y*gl,z*gl]
                r = distance(position, com)
                r_0 = r-1
                radius_perp.append(r/Rsun)
                v_x = g["x-velocity"][x,y,z]
                v_y = g["y-velocity"][x,y,z]
                v_z = g["z-velocity"][x,y,z]
                v_rad = (v_x**2. + v_y**2. + v_z**2.)**0.5
                ve = escape_vel(mass, r)
                a = g_accel(mass, r)
                if v_rad < ve:
                    t = 0.
                    while r > r_0:
                        height = projectile_height(v_rad, dt, a)
                        r = r + height
                        v_rad = velocity(v_rad, a, dt)
                        a = g_accel(mass, r)
                        t = t + dt
                        #if t/day > 52.:
                            #r = 0.0
                    cell = g["CellMass"][x,y,z]
		    cmass_perp.append(cell)
                    time_perp.append(t/day)
                #print x, t/day, "weeks"
    
    bins = range(40)

    plt.clf()
    '''
    plt.hist(time, bins, alpha=0.5, label='Orbital Plane')
    plt.hist(time_perp, bins, alpha=0.5, label="Perpendicular")
    plt.legend(loc='upper right')
    plt.xlabel("Fallback time ($weeks$)")
    plt.ylabel("count")
    '''
    plt.scatter(cmass, time, c='r', marker='o')
    plt.scatter(cmass_perp, time_perp, c='b', marker='^')
    #plt.legend(loc='upper right')
    plt.xlabel("Cell Mass ($g$)")
    plt.ylabel("Fallback time ($weeks$)")
    plt.xlim([0.0, 0.8e26])
    plt.ylim([0.0, 140])
    file_name = "box_frame_" + str(it) + "_time_" + str(time_val) + ".png"
    plt.savefig(file_name)
    it = it + 1


