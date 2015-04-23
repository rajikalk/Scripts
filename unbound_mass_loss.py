#calculate unbound ratio of gas that left the grid.

from yt.mods import *
import numpy as np
import matplotlib.pyplot as plt
import csv

#Import all the timesteps for the series:
init_pf = load("/home/reggie/DATA/Simulations/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD0000/CE0000")
ts = DatasetSeries("/home/reggie/DATA/Simulations/rotation/run1.e-6lessetot_Gcorr_0.75k/hot2/DD00*/CE00*.hierarchy")

txt = open('mass_loss.csv', 'r+')
txt.write('time, unbound massloss \n')

#find particle mass
init_dd = init_pf.h.all_data()
Msun = init_pf['MassUnits'] #grams/Msun
pm = (init_dd['particle_mass'][0] + init_dd['particle_mass'][1])/Msun

#Define arrays
time = []   #time
TM = []   #total mass in the regions
dm_array = []
totalmassloss = []
init_tm = 0
it = 0
for pf in ts:
    #Get current time
    time_val = pf.current_time
    time.append(time_val)
    #Find total mass in region:
    reg = pf.h.region([0.5, 0.5, 0.5], [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    TotalMass = reg.quantities.total_quantity("cell_mass")
    TM.append(TotalMass)
    if it > 0:
        dm_val = TM[it-1] - TotalMass
        tml = init_tm - TotalMass
    else:
        dm_val = 0
        init_tm = TotalMass
        tml = 0
    totalmassloss.append(tml)
    dm_array.append(dm_val)
    it = it + 1
    
dm = np.average(dm_array[1:])
TotalMassLoss = TM[0] - TM[it-1]

UnboundLoss = 0
massloss = []
TotalUnboundLoss = []
line = 0
with open('energy_bound_edge.csv', 'r+') as f:
    reader = csv.reader(f)
    for row in reader:
        if line != 0:
            if it == 0 | it == len(time)-1:
                dm_first = dm/2.
                ratio = float(row[2])/100.
                loss = dm_first * ratio
                massloss.append(loss)
                UnboundLoss = UnboundLoss + loss
                TotalUnboundLoss.append(UnboundLoss)
            else:
                ratio = float(row[2])/100.
                loss = dm * ratio
                massloss.append(loss)
                UnboundLoss = UnboundLoss + loss
                TotalUnboundLoss.append(UnboundLoss)
            print "Unbound loss:", UnboundLoss, "Msun"
        else:
            line = 1
            
for i in range(len(time)):
    txt.write(str(float(time[i])) + ',' + str(float(TotalUnboundLoss[i]/Msun)) + '\n')
          
print "Total mass loss:", TotalMassLoss, "Msun"
print "Unbound mass:", UnboundLoss, "Msun"

plt.clf()
fig = plt.figure()
top = plt.subplot(211)
top.plot(time, totalmassloss, label = 'cumulative mass loss')
top.plot(time, TotalUnboundLoss, label = 'cumulative unbound mass loss')
'''
plt.axvline(x=0.01, color='k')
plt.axvline(x=0.04, color='k')
plt.axvline(x=0.10, color='k')
plt.axvline(x=0.06, color='b', ls='--')
plt.axvline(x=0.02, color='r', ls=':')
plt.axvline(x=0.03, color='r', ls=':')
plt.axvline(x=0.05, color='r', ls=':')
plt.axvline(x=0.07, color='r', ls=':')
plt.axvline(x=0.08, color='r', ls=':')
plt.axvline(x=0.09, color='r', ls=':')
'''
top.legend(loc = 2)
top.set_ylabel("Mass ($M_\odot$)")
bot = plt.subplot(212)
bot.plot(time, massloss, 'r-', label = 'unbound loss at dump')
bot.plot(time, dm_array, 'b-', label = 'total mass loss at dump')
'''
plt.axvline(x=0.01, color='k')
plt.axvline(x=0.04, color='k')
plt.axvline(x=0.10, color='k')
plt.axvline(x=0.06, color='b', ls='--')
plt.axvline(x=0.02, color='r', ls=':')
plt.axvline(x=0.03, color='r', ls=':')
plt.axvline(x=0.05, color='r', ls=':')
plt.axvline(x=0.07, color='r', ls=':')
plt.axvline(x=0.08, color='r', ls=':')
plt.axvline(x=0.09, color='r', ls=':')
'''
#plt.plot(time, Ve, 'b--', label = '$v/v_{escape}$, edge grids')
bot.legend()
#plt.semilogy()
#plt.title("Amount of unbound gas using different methods for $0.75v_{keplerian}$")
bot.set_xlabel("Time (years)")
bot.set_ylabel("Mass ($M_\odot$)")
fig.tight_layout()
plt.savefig('mass_loss_timeseries.png')

