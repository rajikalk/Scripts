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

#Define distance function
def distance(point1, point2):
#Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result

def L_mag(L):
    L = (L[0]**2. + L[1]**2. + L[2]**2.)**0.5
    return L
#load initial
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD0000/CE0000")

#Define values:
dim = init_pf.domain_dimensions[0]
Msun = init_pf['MassUnits']       #grams/Msun
lu = init_pf['LengthUnits']     #length units in cm
tu = init_pf['TimeUnits']           #time units in sec
DensityUnits = init_pf['DensityUnits']
gl = lu/dim             #cm in grid length
Rsun = 69600000000.      # solar radii in  cm
G = 6.674e-8            #gravitational constant in cgs

#Define arrays:
time = []
L_tot = []
L_p = []
L_g = []
coms = []
com_errs = []
header = 0
Lzpos = []
Lpos_val = []
Lzneg = []
Lneg_val = []
Ecc = []
err_l = []
err_h = []

with open('coms.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if header != 0:
            x = float(row[7])*lu
            y = float(row[8])*lu
            z = float(row[9])*lu
            com_val = [x,y,z]
            coms.append(com_val)
            '''
            x_err = float(row[10])*lu
            y_err = float(row[11])*lu
            z_err = float(row[12])*lu
            com_err = [x_err, y_err, z_err]
            com_errs.append(com_err)
            ec = float(row[7])
            Ecc.append(ec)
            '''
        else:
            header = 1

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD00*/CE00*.hierarchy")
f = open('angular_momentum_conservation_z.txt','r+')
f.write("time, L_tot, L_p, L_g, Error_l, Error_h \n")

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun
pmr = (pm1*pm2)/(pm1 + pm2)

it = 0

for pf in ts:
    L_z = 0.
    #L_z_l = 0.
    #L_z_h = 0.
    Lp = 0.
    Lg = 0.
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #Find velocity of the particles
    pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_position_z'][0]]
    pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_position_z'][1]]
    
    '''
    #Calulate angular momentum of particles:
    separation = distance(pp1, pp2)
    alpha = (-1)*G*pm1*pm2
    KE_tot = (0.5)*(pm1*(pv1**2.) + pm2*(pv2**2.))
    GPE_tot = alpha/(separation)
    E_tot = KE_tot + GPE_tot
    e_val = Ecc[it]
    Lp = ((e_val**2. - 1)*((pmr*alpha**2.)/(2.*E_tot)))**0.5
    '''
    com = coms[it]
    #com_err = com_errs[it]
    print com
    
    #Calculate angular momentum of the particles:
    p1 = linear_momentum(pm1, pv1)
    p2 = linear_momentum(pm2, pv2)
    rel1 = rel_position(pp1, com)
    rel2 = rel_position(pp2, com)
    L1 = z_momentum(rel1, p1)
    L2 = z_momentum(rel2, p2)
    Lp = L1 + L2
    L_p.append(Lp)
    '''
    rel_l_err1 = rel_position(pp1, [com[0]-com_err[0], com[1]-com_err[1], com[2]-com_err[2]])
    rel_h_err1 = rel_position(pp1, [com[0]+com_err[0], com[1]+com_err[1], com[2]+com_err[2]])
    rel_l_err2 = rel_position(pp2, [com[0]-com_err[0], com[1]-com_err[1], com[2]-com_err[2]])
    rel_h_err2 = rel_position(pp2, [com[0]+com_err[0], com[1]+com_err[1], com[2]+com_err[2]])
    L1_l = z_momentum(rel_l_err1, p1)
    L1_h = z_momentum(rel_h_err1, p1)
    L2_l = z_momentum(rel_l_err2, p2)
    L2_h = z_momentum(rel_h_err2, p2)
    Lp_l = L1_l + L2_l
    Lp_h = L1_h + L2_h
    
    if Lp_l > Lp_h:
        temp = Lp_l
        Lp_l = Lp_h
        Lp_h = temp
    '''
    #Get the grids in the data dump
    g = pf.h.grids[0]

    #Define counters

    for x in range(dim):
        for y in range(dim):
            for z in range(dim):
                pos = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["TotalMassMsun"][x, y, z]*Msun
                velocity = [g["x-velocity"][x, y, z], g["y-velocity"][x, y, z], g["z-velocity"][x, y, z]]
                p = linear_momentum(mass, velocity)
                rpos = rel_position(pos, com)
                '''
                rpos_l = rel_position(pos, [com[0]-com_err[0], com[1]-com_err[1], com[2]-com_err[2]])
                rpos_h = rel_position(pos, [com[0]+com_err[0], com[1]+com_err[1], com[2]+com_err[2]])           
                '''
                L = (z_momentum(rpos, p))#/mass
                '''
                L_l = (z_momentum(rpos_l, p))
                L_h = (z_momentum(rpos_h, p))
                if L_l > L_h:
                    temp = L_l
                    L_l = L_h
                    L_h = temp
                '''
                L_z = L_z + L
                '''
                L_z_l = L_z_l + L_l
                L_z_h = L_z_h + L_h
                '''
    
    Lg = L_z
    L_g.append(L)
    
    Lt = L_z + Lp
    '''
    Lt_l = Lp_l + L_z_l
    Lt_l = Lt - Lt_l
    Lt_h = Lp_h + L_z_h
    Lt_h = Lt_h - Lt
    '''
    L_tot.append(Lt)
    it = it + 1
    
    f.write(str(time_val) + ',' + str(Lt) + ',' + str(Lp) + ',' + str(Lg) + '\n')
    print "The total amount of angular momentum is:", Lt, "(gcm^2)/s"
    print "The total amount of orbital angular momentum is:", Lp, "(gcm^2)/s"
    print "The total amount of angular momentum in gas:", Lg, "(gcm^2)/s"

