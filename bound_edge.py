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
def Kinetic_Energy(Mass, velocity):
#Takes in mass in g and the velocity in cm/s
#gives energy in ergs
    result = 0.5 * mass * velocity**2.
    return result

#Define arrays:
time = []
bound = []
unbound = []
bm = []
ubm = []
KE = []
TE = []
GPE_g = []
GPE_pg = []
GPE_tot = []
Total = []
mass_ratio = []
E_loss = []
E_loss_cumul = []
M_tot = []

#Import all the timesteps for the series:
'''
ts = TimeSeriesData.from_filenames("/media/DATA/Simulations/smallbox/highres_hot2/DD00*/CE00*.hierarchy")
init_pf = load("/media/DATA/Simulations/smallbox/highres_hot2/DD0000/CE0000")
'''
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/rotation/run1.e-6lessetot_Gcorr_0.75k/DD0000/CE0000")

f = open('energy_bound_edge.csv','r+')
f.write('time, bound, unbound, Total_E, Total_TE, Total_KE, Total_GPE, Total_GPE_g, Total_GPE_pg, E_loss:Eq, Perpendicular, total\n')

#Find the mass of the particles
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMass"][0]
pm2 = init_dd["ParticleMass"][1]

#Find the grid dimensions
dim = init_pf.domain_dimensions[0]
edge = [0, dim-1]
it = 0

#Define values:
Msun = init_pf['MassUnits']       #grams/Msun
lu = init_pf['LengthUnits']     #length units in cm
tu = init_pf['TimeUnits']           #time units in sec
gl = lu/dim             #cm in grid length

for pf in ts:
    time_val = pf.current_time
    time.append(time_val)
    
    dd = pf.h.all_data()
    
    #Find the position of the particles
    pp1 = [dd['particle_position_x'][0]*lu, dd['particle_position_y'][0]*lu, dd['particle_position_z'][0]*lu]
    pp2 = [dd['particle_position_x'][1]*lu, dd['particle_position_y'][1]*lu, dd['particle_position_z'][1]*lu]
    
    #Get the grids in the data dump
    g = pf.h.grids[0]

    #Define counters
    bound_val = 0.
    unbound_val = 0.
    TotalE = 0.
    TotalKE = 0.
    TotalTE = 0.
    TotalGPEg = 0.
    TotalGPEpg = 0.
    TotalGPE = 0.
    mass_tot = 0.
    
    #Left and right edge:
    for x in edge:
        for y in range(dim):
            for z in range(dim):
                position = [x*gl,y*gl,z*gl]
                mass = g["CellMassMsun"][x,y,z]*Msun #in grams
                mass_tot = mass_tot + mass
                velocity = (((g["x-velocity"][x,y,z])**2.)+((g["y-velocity"][x,y,z])**2.)+((g["z-velocity"][x,y,z])**2.))**0.5
                KE_val = Kinetic_Energy(mass, velocity)/mass
                TotalKE = TotalKE + KE_val
                TE_val = g["ThermalEnergy"][x,y,z] #appends thermal energy in ergs/g.
                TotalTE = TotalTE + TE_val
                GPEg = g["Grav_Potential"][x,y,z] #appends thermal energy in ergs/g.
                TotalGPEg = TotalGPEg + GPEg
                radius1 = distance(position, pp1)
                radius2 = distance(position, pp2)
                GPE1 = GPE(pm1, mass, radius1)/mass
                GPE2 = GPE(pm2, mass, radius2)/mass
                TotalGPEpg = TotalGPEpg + GPE1 + GPE2
                total = KE_val + TE_val + GPEg + GPE1 + GPE2
                TotalE = TotalE + total
                if total > 0.:
                    unbound_val = unbound_val + mass
                else:
                    bound_val = bound_val + mass

    
    #Front and back edge:
    for x in range(dim):
        for y in edge:
            for z in range(dim):
                position = [x*gl,y*gl,z*gl]
                mass = g["CellMassMsun"][x,y,z]*Msun #in grams
                mass_tot = mass_tot + mass
                velocity = (((g["x-velocity"][x,y,z])**2.)+((g["y-velocity"][x,y,z])**2.)+((g["z-velocity"][x,y,z])**2.))**0.5
                KE_val = Kinetic_Energy(mass, velocity)/mass
                TotalKE = TotalKE + KE_val
                TE_val = g["ThermalEnergy"][x,y,z] #appends thermal energy in ergs/g.
                TotalTE = TotalTE + TE_val
                GPEg = g["Grav_Potential"][x,y,z] #appends thermal energy in ergs/g.
                TotalGPEg = TotalGPEg + GPEg
                radius1 = distance(position, pp1)
                radius2 = distance(position, pp2)
                GPE1 = GPE(pm1, mass, radius1)/mass
                GPE2 = GPE(pm2, mass, radius2)/mass
                TotalGPEpg = TotalGPEpg + GPE1 + GPE2
                total = KE_val + TE_val + GPEg + GPE1 + GPE2
                TotalE = TotalE + total
                if total > 0.:
                    unbound_val = unbound_val + mass
                else:
                    bound_val = bound_val + mass
    
    #Calculate average energy on equatorial plane:
    Ave_Eq = TotalE/((dim**2.)*4.)
    
    #Top and bottom edge:
    Totop = 0.
    mass_perp = 0.
    for x in range(dim):
        for y in range(dim):
            for z in edge:
                position = [x*gl,y*gl,z*gl]
                mass = g["CellMassMsun"][x,y,z]*Msun #in grams
                mass_tot = mass_tot + mass
                mass_perp = mass_perp + mass
                velocity = (((g["x-velocity"][x,y,z])**2.)+((g["y-velocity"][x,y,z])**2.)+((g["z-velocity"][x,y,z])**2.))**0.5
                KE_val = Kinetic_Energy(mass, velocity)/mass
                TotalKE = TotalKE + KE_val
                TE_val = g["ThermalEnergy"][x,y,z] #appends thermal energy in ergs/g.
                TotalTE = TotalTE + TE_val
                GPEg = g["Grav_Potential"][x,y,z] #appends thermal energy in ergs/g.
                TotalGPEg = TotalGPEg + GPEg
                radius1 = distance(position, pp1)
                radius2 = distance(position, pp2)
                GPE1 = GPE(pm1, mass, radius1)/mass
                GPE2 = GPE(pm2, mass, radius2)/mass
                TotalGPEpg = TotalGPEpg + GPE1 + GPE2
                total = KE_val + TE_val + GPEg + GPE1 + GPE2
                TotalE = TotalE + total
                Totop = Totop + total
                if total > 0.:
                    unbound_val = unbound_val + mass
                else:
                    bound_val = bound_val + mass
                
    #Calculate average energy and mass fraction on top and bot:
    Ave_Perp = Totop/((dim**2.)*2.)
    mass_frac = mass_perp/mass_tot
    
    #Find energy loss:
    mtot = sum(dd['CellMass'])
    if time_val == 0.0:
        mloss = 0
    else:
        mloss = M_tot[it-1] - mtot
    M_tot.append(mtot)
    #on the plane:
    E_pl = Ave_Eq*(1 - mass_frac)*mloss
    #perpendidular:
    E_perp = Ave_Perp*mass_frac*mloss
    #total loss
    E_loss_tot = E_pl + E_perp
    E_loss.append(E_loss_tot)
    if time_val == 0.0:
        cumul_val = 0.0
    else:
        cumul_val = E_loss_cumul[it-1] + E_loss_tot
    E_loss_cumul.append(cumul_val)
             
    KE.append(TotalKE)
    TE.append(TotalTE)
    GPE_g.append(TotalGPEg)
    GPE_pg.append(TotalGPEpg)
    TotalGPE = TotalGPEg + TotalGPEpg
    GPE_tot.append(TotalGPE)
    Total.append(TotalE)
    TotalMass = unbound_val + bound_val
    unbound_percent = (unbound_val*100.)/(TotalMass)
    bound_percent = (bound_val*100.)/(TotalMass)
    unbound.append(unbound_percent)
    bound.append(bound_percent)
    bm.append(bound_val)
    ubm.append(unbound_val)
    
    it = it + 1

    f.write(str(time_val) + ',' + str(bound_percent) + ',' + str(unbound_percent) + ',' + str(TotalE) + ',' + str(TotalTE) + ',' + str(TotalKE) + ',' + str(TotalGPE) + ',' + str(TotalGPEg) + ',' + str(TotalGPEpg) + ',' + str(E_pl) + ',' + str(E_perp) + ',' + str(E_loss_tot) + '\n')
    print "The amount of unbound gas at edge is", unbound_percent, "%"
    print "The amount of bound gas at edge is", bound_percent, "%"
    print "Total energy loss is", E_loss_tot, "ergs"
    print "Cumulative energy loss", cumul_val, "ergs"
'''
fig = plt.figure(1)
fig.clf()
fig.set_size_inches(7,7)
#fig.tight_layout()
#fig.yscale("log")
top = fig.add_subplot(111)
top.plot(time, unbound)
#top.set_ylim([90, 100])
#top.set_xlim([0.0, 0.1])
top.set_title('Unbound gas (using velocities), at the edge')
top.set_xlabel('Time ($years$)')
top.set_ylabel('Unbound ratio (%)')
plt.savefig("edge_unbound-zoom.png")
'''
