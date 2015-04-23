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
Msun = 1.9891e+33       #grams/Msun
lu = 10000000000000     #length units in cm
tu = 31558000           #time units in sec
DensityUnits = 1.9891e-06
gl = lu/256             #cm in grid length

#Define arrays:
time = []
L_tot = []
L_p = []
L_g = []
coms = []
header = 0
Lzpos = []
Lpos_val = []
Lzneg = []
Lneg_val = []

with open('com.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            x = float(row[4])*lu
            y = float(row[5])*lu
            z = float(row[6])*lu
            com_val = [x,y,z]
            coms.append(com_val)
        else:
            header = 1

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")
f = open('angular_momentum_conservation_total.txt','r+')
f.write("time, L_tot:x, y, z, L_p:x, y, z, L_g:x, y, z \n")

#Save directory:
save_dir = "~/YT_output/run1/e-6lessetot/longer/"

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]
it = 0

for pf in ts:
    Lx = 0
    Ly = 0
    Lz = 0
    Lx_g = 0
    Ly_g = 0
    Lz_g = 0
    time_val = pf.current_time
    time.append(time_val)
    
    #Reset angular momentum counter:
    L_total = 0.
    L_tot_g = 0.
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #Find velocity of the particles
    pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_position_z'][0]]
    pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_position_z'][1]]
    
    com = coms[it]
    print com
    
    #Calculate angular momentum of the particles:
    p1 = linear_momentum(pm1, pv1)
    p2 = linear_momentum(pm2, pv2)
    rel1 = rel_position(pp1, com)
    rel2 = rel_position(pp2, com)
    L1 = angular_momentum(rel1, p1)
    L2 = angular_momentum(rel2, p2)
    Lp_x = L1[0] + L2[0]
    Lp_y = L1[1] + L2[1]
    Lp_z = L1[2] + L2[2]
    Lx = Lp_x
    Ly = Lp_y
    Lz = Lp_z
    '''
    if L1z > 0:
        Lpos = Lpos + 1
        Lpos_value = Lpos_value + L1z
    else:
        Lneg = Lneg + 1
        Lneg_value = Lneg_value + L1z
    L2z = z_momentum(rel2, p2)
    if L2z > 0:
        Lpos = Lpos + 1
        Lpos_value = Lpos_value + L2z
    else:
        Lneg = Lneg + 1
        Lneg_value = Lneg_value + L2z
    '''
    lp = (Lx**2. + Ly**2. + Lz**2.)**0.5
    L_p.append(lp)

    #Get the grids in the data dump
    g = pf.h.grids[0]

    #Define counters

    for x in range(dim):
        for y in range(dim):
            for z in range(dim):
                pos = [x*gl,y*gl,z*gl]
                mass = g["TotalMassMsun"][x, y, z]*Msun
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                L = (angular_momentum(rpos, p))#/mass
                '''
                if L > 0:
                    Lpos = Lpos + 1
                    Lpos_value = Lpos_value + L
                else:
                    Lneg = Lneg + 1
                    Lneg_value = Lneg_value + L
                '''
                Lx_g = Lx_g + L[0]
                Ly_g = Ly_g + L[1]
                Lz_g = Lz_g + L[2]
                
    lg = (Lx_g**2. + Ly_g**2. + Lz_g**2.)**0.5
    L_g.append(lg)
    Lx = Lx + Lx_g
    Ly = Ly + Ly_g
    Lz = Lz + Lz_g
    L = (Lx**2. + Ly**2. + Lz**2.)**0.5
    L_tot.append(L)
    it = it + 1
    '''
    Lzpos.append(Lpos)
    Lzneg.append(Lneg)
    Lpos_val.append(Lpos_value)
    Lneg_val.append(Lneg_value)
    L_tot.append(L_total)
    L_g.append(L_tot_g)
    '''
    f.write(str(time_val) + ',' + str(Lx) + ',' + str(Ly) + ',' + str(Lz) + ',' + str(Lp_x) + ',' + str(Lp_y) + ',' + str(Lp_z) + ',' + str(Lx_g) + ',' + str(Ly_g) + ',' + str(Lz_g) + '\n')
    print "The total amount of angular momentum is:", L, "(gcm^2)/s"
    print "The total amount of orbital angular momentum is:", lp, "(gcm^2)/s"
    print "The total amount of angular momentum in gas:", lg, "(gcm^2)/s"

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

