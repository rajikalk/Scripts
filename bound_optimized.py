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
def Grav_Pot(Mass, cellMass, distance):
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
GPE = []
Total_g = []
Total_p = []
Total = []

#Import all the timesteps for the series:
'''
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/highres_hot2/DD00*/CE00*.hierarchy")
init_pf = load("/media/DATA/Simulations/smallbox/highres_hot2/DD0000/CE0000.hierarchy")
'''
ts = TimeSeriesData.from_filenames("/home/science/staff/reggie/Simulation/Hot_fb_0.5k/DD00*/CE00*.hierarchy")
init_pf = load("/home/science/staff/reggie/Simulation/Hot_fb_0.5k/DD0000/CE0000")

f = open('energy_conservation.csv','r+')
f.write('time, bound, unbound, Total_E, Total_TE, Total_KE, Total_GPE, Total_KE_g, Total_KE_p, Total_GPE_g, Total_GPE_pg, Total_GPE_pp, Total_p, Total_g \n')

#Save directory:
save_dir = "~/YT_output/run1.e-6lessetot/longer/"

#Find size of domain
dim = init_pf.domain_dimensions[0]

#Define values:
Msun = init_pf['MassUnits']       #grams/Msun
lu = init_pf['LengthUnits']     #length units in cm
tu = init_pf['TimeUnits']           #time units in sec
gl = lu/dim             #cm in grid length

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm = init_dd["ParticleMass"]

for pf in ts:
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #Find velocity of the particles
    pv = (((dd['particle_velocity_x'])**2 + (dd['particle_velocity_y'])**2 + (dd['particle_position_z'])**2)**0.5)
    
    #Find kinetic energy of particles
    KEp1 = Kinetic_Energy(pm[0], pv[0])
    KEp2 = Kinetic_Energy(pm[1], pv[1])
    TotalKEp = KEp1 + KEp2
    KE_p.append(TotalKEp)
    
    #Find kinetic energy of gas
    KEg = dd['KineticEnergy']*dd['CellVolume']
    TotalKEg = sum(KEg)
    KE_g.append(TotalKEg)
    TotalKE = TotalKEp + TotalKEg
    KE.append(TotalKE)
    
    #Find gravitational potential energy of the particles on each other
    pp_sep = distance(pp1, pp2)
    TotalGPEpp = Grav_Pot(pm[0], pm[1], pp_sep)
    GPE_pp.append(TotalGPEpp)
    
    #Find gravitational potential energy of the gas
    GPEg = dd['Grav_Potential']*dd['CellMass']
    TotalGPEg = sum(GPEg)/2.0 #divide by two to avoid doubling up the potential measured between points.
    GPE_g.append(TotalGPEg)
    
    #Find thermal energy of gas
    TE_val = dd['ThermalEnergy']*dd['CellMass']
    TotalTE = sum(TE_val)
    TE.append(TotalTE)

    #Get the grids in the data dump
    g = pf.h.grids[0]
    
    #Calculate GPE of particles on gas
    GPEpg = []
    for x in range(dim):
        for y in range(dim):
            for z in range(dim):
                position = [(x+0.5)*gl,(y+0.5)*gl,(z+0.5)*gl]
                mass = g["CellMass"][x,y,z] #in grams
                
                #Calculate GPE in grid:
                radius1 = distance(position, pp1)
                radius2 = distance(position, pp2)
                GPE1 = Grav_Pot(pm[0], mass, radius1)
                GPE2 = Grav_Pot(pm[1], mass, radius2)
                GPEpg_val = GPE1 + GPE2
                GPEpg.append(GPEpg_val)
    
    GPEpg = np.array(GPEpg)
    TotalGPEpg = sum(GPEpg)
    GPE_pg.append(TotalGPEpg)
    TotalGPE = TotalGPEpp + TotalGPEpg + TotalGPEg
    GPE.append(TotalGPE)
    
    #Now calculate the total energy of the gas
    totalg = KEg + GPEg + GPEpg + TE_val
    tg_val = sum(totalg)
    Total_g.append(tg_val)
    CM = dd['CellMass']
    TM = sum(CM)
    
    bm = 0
    um = 0
    #Find unbound mass:
    for i in range(totalg.size):
        if totalg[i] < 0:
            bm = bm + CM[i]
        else:
            um = um + CM[i]
    bp = (bm/TM)*100.
    up = (um/TM)*100.
    bound.append(bp)
    unbound.append(up)
    
    #Calculate the total particle energy
    totalp = TotalKEp + TotalGPEpp
    Total_p.append(totalp)
    
    #Finally calculate the total energy
    TotalE = totalp + tg_val
    Total.append(TotalE)
    
    f.write(str(time_val) + ',' + str(bp) + ',' + str(up) + ',' + str(TotalE) + ',' + str(TotalTE) + ',' + str(TotalKE) + ',' + str(TotalGPE) + ',' + str(TotalKEg) + ',' + str(TotalKEp) + ',' + str(TotalGPEg) + ',' + str(TotalGPEpg) + ',' + str(TotalGPEpp) +  ',' + str(totalp) + ',' + str(tg_val) + '\n')
    
    print "The amount of unbound gas is", up, "%"
    print "The amount of bound gas is", bp, "%"
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
