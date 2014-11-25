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
    mag = (L[0]**2. + L[1]**2. + L[2]**2.)**0.5
    if L[2] < 0:
        mag = (-1)*mag
    return mag

#Define values:
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000.     #length units in cm
tu = 31558000.           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256.             #cm in grid length

#Define arrays:
time = []
coms = []
header = 0

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

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
f = open('angular_momentum_edge.csv','r+')
f.write("time, L:plane, perp, M_frac:plane, perp \n")

#Save directory:
save_dir = "~/YT_output/run1/e-6lessetot/longer/"

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]
edge = [0, dim-1]
it = 0

for pf in ts:
    mass_tot = 0.
    m_plane = 0.
    m_perp = 0.
    L_plane = 0.
    L_perp = 0.
    time_val = pf.current_time
    time.append(time_val)
    
    #Reset angular momentum counter:
    dd = pf.h.all_data()
    
    com = coms[it]

    #Get the grids in the data dump
    g = pf.h.grids[0]
    
    #Left and right edge:
    for x in edge:
        for y in range(dim):
            for z in range(dim):
                pos = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMassMsun"][x, y, z]*Msun
                mass_tot = mass_tot + mass
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = z_momentum(rpos, p)
                mag = L/mass
		m_plane = m_plane + mass
		L_plane = L_plane + mag
    
    #Front and back edge:
    for y in edge:
        for x in range(dim):
            for z in range(dim):
                pos = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMassMsun"][x, y, z]*Msun
                mass_tot = mass_tot + mass
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = z_momentum(rpos, p)
                mag = L/mass
		m_plane = m_plane + mass
		L_plane = L_plane + mag
    
    #Top and bottom edge:
    for z in edge:
        for y in range(dim):
            for x in range(dim):
                pos = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMassMsun"][x, y, z]*Msun
                mass_tot = mass_tot + mass
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = z_momentum(rpos, p)
                mag = L/mass
		m_perp = m_perp + mass
		L_perp = L_perp + mag
    
    #Calculate specific angular momentum per cell:
    Lsp_plane = L_plane/((dim**2.)*4.)
    Lsp_perp = L_perp/((dim**2.)*2.)
    
    #Calculate mass fractions:
    m_plane = m_plane/mass_tot
    m_perp = m_perp/mass_tot
    
    it = it + 1
    "time, L:plane, perp, M_frac:plane, perp"
    f.write(str(time_val) + ',' + str(L_plane) + ',' + str(L_perp) + ',' + str(m_plane) + ',' + str(m_perp) + '\n')
    
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
