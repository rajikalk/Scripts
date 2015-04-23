#Calculating the bound gas

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt

#Define distance function
def distance(point1, point2):
#Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result

#Define GPE function
def GPE(Mass, cellMass, distance):
#Takes in the masses in g and the separation in cm
#gives energy in erg
    top = -6.67259*(10**(-8.))*Mass*cellMass
    bottom = distance
    result = top/bottom
    return result
    
#Define KE function
def Kinetic_Energy(mass, velocity):
#Takes in mass in g and the velocity in cm/s
#gives energy in ergs
    result = 0.5 * mass * (velocity**2.)
    return result

#Define arrays:
time = []
bound = []
unbound = []
KE_g = []
KE_p = []
KE = []
TE = []
GPE_g = []
GPE_pg = []
GPE_pp = []
GPE_tot = []
Total = []

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD0000/CE0000")
f = open('energy_conservation.txt','r+')
f.write('time, bound, unbound, Total_E, Total_TE, Total_KE, Total_GPE, Total_KE_g, Total_KE_p, Total_GPE_g, Total_GPE_pg, Total_GPE_pp\n')

#Save directory:
save_dir = "~/YT_output/run1.e-6lessetot/longer/"

#Find size of domain
dim = init_pf.domain_dimensions[0]

#Define values:
Msun = init_pf['MassUnits']       #grams/Msun
lu = init_pf['cm']     #length units in cm
tu = init_pf['TimeUnits']           #time units in sec
DensityUnits = init_pf['DensityUnits']
gl = lu/dim             #cm in grid length

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun
pm = pm1 + pm2

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]

for pf in ts:
    counter = 0
    time_val = pf.current_time
    time.append(time_val)
    
    #find total mass:
    reg = pf.h.region([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    TotalMass = (reg.quantities["TotalQuantity"]("TotalMassMsun")[0]*Msun)-pm
    print "TotalMass:", TotalMass
    
    #Define counters
    bound_val = 0.
    unbound_val = 0.
    TotalE = 0.
    TotalKEg = 0.
    TotalKEp = 0.
    TotalKE = 0.
    TotalTE = 0.
    TotalGPEg = 0.
    TotalGPEpg = 0.
    TotalGPEpp = 0.
    TotalGPE = 0.
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #Find velocity of the particles
    pv1 = (((dd['particle_velocity_x'][0])**2 + (dd['particle_velocity_y'][0])**2 + (dd['particle_position_z'][0])**2)**0.5)
    pv2 = (((dd['particle_velocity_x'][1])**2 + (dd['particle_velocity_y'][1])**2 + (dd['particle_position_z'][1])**2)**0.5)
    
    #Find kinetic energy of particles
    KEp1 = Kinetic_Energy(pm1, pv1)
    KEp2 = Kinetic_Energy(pm2, pv2)
    KEp = KEp1 + KEp2
    KE_p.append(KEp)
    TotalKE = TotalKE + KEp
    TotalKEp = KEp
    
    #Find gravitational potential energy of the particles on each other
    pp_sep = distance(pp1, pp2)
    GPEpp = GPE(pm1, pm2, pp_sep)
    GPE_pp.append(GPEpp)
    TotalGPE = TotalGPE + GPEpp
    TotalGPEpp = GPEpp

    #Get the grids in the data dump
    g = pf.h.grids[0]

    for x in range(dim):
        for y in range(dim):
            for z in range(dim):
                position = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMassMsun"][x,y,z]*Msun #in grams
                velocity = (((g["x-velocity"][x,y,z])**2.)+((g["y-velocity"][x,y,z])**2.)+((g["z-velocity"][x,y,z])**2.))**0.5
                
                #Calculate KE in grid:
                KineticEnergy = Kinetic_Energy(mass, velocity)
                TotalKE = TotalKE + KineticEnergy
                TotalKEg = TotalKEg + KineticEnergy
                
                #Calculate TE in grid:
                ThermalEnergy = g["ThermalEnergy"][x,y,z]*(mass)
                TotalTE = TotalTE +ThermalEnergy
                
                #Calculate GPE in grid:
                GPEg = g["Grav_Potential"][x,y,z]*(mass)
                TotalGPEg = TotalGPEg + GPEg
                radius1 = distance(position, pp1)
                radius2 = distance(position, pp2)
                GPE1 = GPE(pm1, mass, radius1)
                GPE2 = GPE(pm2, mass, radius2)
                TotalGPEpg = TotalGPEpg + GPE1 + GPE2
                
                total = KineticEnergy + GPEg + GPE1 + GPE2 + ThermalEnergy
                if total > 0.:
                    unbound_val = unbound_val + mass
                else:
                    bound_val = bound_val + mass
    
    #Add total energies to arrays
    KE.append(TotalKE)
    KE_g.append(TotalKEg)
    TE.append(TotalTE)
    TotalGPE = TotalGPE + TotalGPEpg + TotalGPEg/2
    GPE_tot.append(TotalGPE)
    GPE_g.append(TotalGPEg/2) #divide by two to avoid doubling up the potential measured between points.
    GPE_pg.append(TotalGPEpg)
    TotalE = TotalKE + TotalGPE + TotalTE
    Total.append(TotalE)
    

    #Calculate bound and unbound percentages
    unbound_percent = (unbound_val*100)/(TotalMass)
    bound_percent = (bound_val*100.)/(TotalMass)
    unbound.append(unbound_percent)
    bound.append(bound_percent)
    
    f.write(str(time_val) + ',' + str(bound_percent) + ',' + str(unbound_percent) + ',' + str(TotalE) + ',' + str(TotalTE) + ',' + str(TotalKE) + ',' + str(TotalGPE) + ',' + str(TotalKEg) + ',' + str(TotalKEp) + ',' + str(TotalGPEg) + ',' + str(TotalGPEpg) + ',' + str(TotalGPEpp) +  '\n')
    
    print "The amount of unbound gas is", unbound_percent, "%"
    print "The amount of bound gas is", bound_percent, "%"
    print "The total energy is", TotalE, "ergs"
    
'''
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(7,7)
#fig.tight_layout()
#fig.yscale("log")
top = fig.add_subplot(211)
top.plot(time, unbound, 'k-', time, bound, 'b--')
#top.set_ylim([0.0, 0.15])
#top.set_xlim([0.0, 0.1])
top.set_title('Unbound gas (using energy)')
top.set_xlabel('Time ($years$)')
top.set_ylabel('Bound and Unbound ratios (%)')
bot = fig.add_subplot(212)
bot.plot(time, Total, 'k-', time, KE, 'b--', time, GPE_tot, 'r-.', time, KE_g, 'c:', time, KE_p, 'ms', time, GPE_g, 'y^', time,GPE_pg , 'ko', time, GPE_pp, 'b+')
#top.set_ylim([0.0, 0.15])
#top.set_xlim([0.0, 0.1])
top.set_title('Total Energy')
top.set_xlabel('Time ($years$)')
top.set_ylabel('Total Energy (ergs)')
#fig.legend()
plt.savefig("unbound.png")
'''
