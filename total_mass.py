# This program plots the timeseries of the the total mass in the grid

from yt.mods import *
import matplotlib.pyplot as plt
import csv

#Import all the timesteps for the series:
#ts = DatasetSeries("DATA/Simulations/rotation/run1.e-6lessetot_Gcorr_0.75k/hot/DD00*/CE00*.hierarchy")
ts = TimeSeriesData.from_filenames('/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0*/data0*.hierarchy')
pf = load("/disks/ceres/makemake/acomp/riaconi/rxi552/Binary-JC_primary-0.6Msun_companion-3_stellar_radii-256grid/DD0014/data0014")
f = open('total_mass.csv','r+')
f.write('time, mass, massloss\n')

#Define arrays
times = []      #Stores the time
mass = []

#mass of particles
init_dd = pf.h.all_data()
Msun = pf['MassUnits'] #grams/Msun
Rsun = 6.96*(10**10.)   #cm/rsun
pm = (init_dd['ParticleMass'][0] + init_dd['ParticleMass'][1])/Msun
it = 0

for pf in ts:
    dd = pf.h.all_data()
    time = pf.current_time
    times.append(time)
    #Create region that covers entire grid:
    reg = pf.h.region([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    #Find total mass in region:
    TotalMass = sum(reg['CellMassMsun'])
    mass.append(TotalMass)
    if it == 0:
        massloss = 0
    else:
        massloss = (mass[it-1]-TotalMass)*Msun
    string = str(time) + ',' + str(TotalMass) + ',' + str(massloss) + '\n'
    print string
    f.write(string)
    it = it + 1

for i in range(len(mass)):
    print times[i], mass[i]
    
plt.clf()
p = plt.plot(times, mass)
plt.title("Total mass of gas within box")
plt.xlabel("Time (years)")
#plt.xlim([0.0, 0.14])
plt.ylabel("Total Mass ($M_\odot$)")
plt.savefig('/home/reggie/DATA/YT_output/hot/mass_timeseries.png')
