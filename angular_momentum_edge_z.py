#Momentum conservation

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
    
def rel_position(pos1, pos2):
    x = pos1[0] - pos2[0]
    y = pos1[1] - pos2[1]
    z = pos1[2] - pos2[2]
    rel = [x,y,z]
    return rel
    
def angular_momentum(rel_position, momentum):
    #Everything in cgs
    Lx = rel_position[1]*momentum[2] - rel_position[2]*momentum[1]
    Ly = rel_position[2]*momentum[0] - rel_position[0]*momentum[2]
    Lz = rel_position[0]*momentum[1] - rel_position[1]*momentum[0]
    L = [Lx, Ly, Lz]
    return L
    
def z_momentum(rel_position, momentum):
    #Everything in cgs
    Lz = rel_position[0]*momentum[1] - rel_position[1]*momentum[0]
    return Lz

def L_mag(L):
    L = (L[0]**2. + L[1]**2. + L[2]**2.)**0.5
    return L

#Define arrays:
time = []
coms = []
header = 0

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0*/data0*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0014/data0014")
f = open('angular_momentum_edge_z.csv','r+')
f.write("time, Lpos: left, right, front, back, bot, top, Lneg: left, right, front, back, bot, top, M: left, right, front, back, bot, top, Mneg: left, right, front, back, bot, top \n")

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMass"][0]
pm2 = init_dd["ParticleMass"][1]

#Define values:
dim = init_pf.domain_dimensions[0]
Msun = init_pf['MassUnits']       #grams/Msun
lu = init_pf['LengthUnits']     #length units in cm
tu = init_pf['TimeUnits']           #time units in sec
DensityUnits = init_pf['DensityUnits']
gl = lu/dim             #cm in grid length

with open('coms.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            x = float(row[7])*lu
            y = float(row[8])*lu
            z = float(row[9])*lu
            com_val = [x,y,z]
            coms.append(com_val)
        else:
            header = 1

f = open('angular_momentum_edge_z.csv','r+')
f.write("time, Lpos: left, right, front, back, bot, top, Lneg: left, right, front, back, bot, top, M: left, right, front, back, bot, top, Mneg: left, right, front, back, bot, top \n")

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]
edge = [0, dim-1]
it = 0

for pf in ts:
    mass_tot = 0.
    m_left = 0.
    m_right = 0.
    m_front = 0.
    m_back = 0.
    m_top = 0.
    m_bot = 0.
    L_left = 0.
    L_right = 0.
    L_front = 0.
    L_back = 0.
    L_top = 0.
    L_bot = 0.
    neg_L_left = 0.
    neg_L_right = 0.
    neg_L_front = 0.
    neg_L_back = 0.
    neg_L_top = 0.
    neg_L_bot = 0.
    neg_left = 0.
    neg_right = 0.
    neg_front = 0.
    neg_back = 0.
    neg_top = 0.
    neg_bot = 0.
    left_count = 0.
    right_count = 0.
    front_count = 0.
    back_count = 0.
    bot_count = 0.
    top_count = 0.
    neg_left_count = 0.
    neg_right_count = 0.
    neg_front_count = 0.
    neg_back_count = 0.
    neg_bot_count = 0.
    neg_top_count = 0.
    time_val = pf.current_time
    time.append(time_val)
    
    #Reset angular momentum counter:
    L_total = 0.
    
    dd = pf.h.all_data()
    
    com = coms[it]

    #Get the grids in the data dump
    g = pf.h.grids[0]
    
    #Left and right edge:
    side = 0
    for x in edge:
        for y in range(dim):
            for z in range(dim):
                pos = [x*gl,y*gl,z*gl]
                mass = g["CellMass"][x, y, z]
                mass_tot = mass_tot + mass
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = z_momentum(rpos, p)
                mag = L/mass
                if side == 0:
                    m_left = m_left + mass
                    if L < 0:
                        neg_left = neg_left + mass
                        neg_L_left = neg_L_left + mag
                        neg_left_count = neg_left_count + 1
                    else:
                        L_left = L_left + mag
                        left_count = left_count + 1
                if side == 1:
                    m_right = m_right + mass
                    if L < 0:
                        neg_right = neg_right + mass
                        neg_L_right = neg_L_right + mag
                        neg_right_count = neg_right_count + 1
                    else:
                        L_right = L_right + mag
                        right_count = right_count + 1
        side = 1
    
    #Front and back edge:
    side = 0
    for y in edge:
        for x in range(dim):
            for z in range(dim):
                pos = [x*gl,y*gl,z*gl]
                mass = g["CellMass"][x, y, z]
                mass_tot = mass_tot + mass
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = z_momentum(rpos, p)
                mag = L/mass
                if side == 0:
                    m_front = m_front + mass
                    if L < 0:
                        neg_front = neg_front + mass
                        neg_L_front = neg_L_front + mag
                        neg_front_count = neg_front_count + 1
                    else:
                        L_front = L_front + mag
                        front_count = front_count + 1
                if side == 1:
                    m_back = m_back + mass
                    if L < 0:
                        neg_back = neg_back + mass
                        neg_L_back = neg_L_back + mag
                        neg_back_count = neg_back_count + 1
                    else:
                        L_back = L_back + mag
                        back_count = back_count + 1
        side = 1
    
    #Top and bottom edge:
    side = 0
    for z in edge:
        for y in range(dim):
            for x in range(dim):
                pos = [x*gl,y*gl,z*gl]
                mass = g["CellMass"][x, y, z]
                mass_tot = mass_tot + mass
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = z_momentum(rpos, p)
                mag = L/mass
                if side == 0:
                    m_bot = m_bot + mass
                    if L < 0:
                        neg_bot = neg_bot + mass
                        neg_L_bot = neg_L_bot + mag
                        neg_bot_count = neg_bot_count + 1
                    else:
                        L_bot = L_bot + mag
                        bot_count = bot_count + 1
                if side == 1:
                    m_top = m_top + mass
                    if L < 0:
                        neg_top = neg_top + mass
                        neg_L_top = neg_L_top + mag
                        neg_top_count = neg_top_count + 1
                    else:
                        L_top = L_top + mag
                        top_count = top_count + 1
        side = 1 
    
    #Calculate specific angular momentum per cell:
    if left_count > 0:
        L_left = L_left/left_count
    if neg_left_count > 0:
        neg_L_left = neg_L_left/neg_left_count
    if right_count > 0:
        L_right = L_right/right_count
    if neg_right_count > 0:
        neg_L_right = neg_L_right/neg_right_count
    if front_count > 0:
        L_front = L_front/front_count
    if neg_front_count > 0:
        neg_L_front = neg_L_front/neg_front_count
    if back_count > 0:
        L_back = L_back/back_count
    if neg_back_count > 0:
        neg_L_back = neg_L_back/neg_back_count
    if bot_count > 0:
        L_bot = L_bot/bot_count
    if neg_bot_count > 0:
        neg_L_bot = neg_L_bot/neg_bot_count
    if top_count > 0:
        L_top = L_top/top_count
    if neg_top_count > 0:
        neg_L_top = neg_L_top/neg_top_count
    
    #Calculate mass fractions:
    left = m_left/mass_tot
    right = m_right/mass_tot
    front = m_front/mass_tot
    back = m_back/mass_tot
    bot = m_bot/mass_tot
    top = m_top/mass_tot
    neg_left = neg_left/m_left
    neg_right = neg_right/m_right
    neg_front = neg_front/m_front
    neg_back = neg_back/m_back
    neg_bot = neg_bot/m_bot
    neg_top = neg_top/m_top
    
    it = it + 1
    
    print 'time:', time_val
    
    f.write(str(time_val) + ',' + str(L_left) + ',' + str(L_right) + ',' + str(L_front) + ',' + str(L_back) + ',' + str(L_bot) + ',' + str(L_top) + ',' + str(neg_L_left) + ',' + str(neg_L_right) + ',' + str(neg_L_front) + ',' + str(neg_L_back) + ',' + str(neg_L_bot) + ',' + str(neg_L_top) + ',' + str(left) + ',' + str(right) + ',' + str(front) + ',' + str(back) + ',' + str(bot) + ',' + str(top) + ',' + str(neg_left) + ',' + str(neg_right) + ',' + str(neg_front) + ',' + str(neg_back) + ',' + str(neg_bot) + ',' + str(neg_top) + '\n')
    
"""
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(7,7)
#fig.tight_layout()
#fig.yscale("log")
top = fig.add_subplot(111)
top.plot(time, L_tot, 'k-', time, L_p, 'b--', time, L_g, 'r-.')
#top.set_ylim([0.0, 0.15])
#top.set_xlim([0.0, 0.1])
top.set_title('Total Angular Momentum')
top.set_xlabel('Time ($years$)')
top.set_ylabel('Total Angular Momentum ($gcm^2/s$)')
plt.savefig("Total_Angular_Momentum.png")
"""
