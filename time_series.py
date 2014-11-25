# This program plots the timeseries of the particle separation and maximum density after time.

from yt.mods import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy as np

#Define distance function
def distance(point1, point2):
#Takes in the position as cm and give the separation in cm
    x_diff_sq = (point1[0] - point2[0])**2.
    y_diff_sq = (point1[1] - point2[1])**2.
    z_diff_sq = (point1[2] - point2[2])**2.
    result = (((x_diff_sq) + (y_diff_sq) + (z_diff_sq))**(0.5))
    return result

#Import all the timesteps for the series:
ts = TimeSeriesData.from_filenames("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/run2.e-4/DD00*/CE00*.hierarchy")
init_pf = load("/disks/ceres/makemake/acomp/jstaff/rajika/smallbox/run2.e-4/DD0000/CE0000")
f = open('separation.csv','r+')
f.write('time, separation\n')
#f.write('time, TE\n')

#Find the mass of the particles
Msun = init_pf['MassUnits']       #grams/Msun       #grams/Msun
init_dd = init_pf.h.all_data()
pm1 = init_dd["ParticleMassMsun"][0]*Msun
pm2 = init_dd["ParticleMassMsun"][1]*Msun

#Define arrays
times = []
sep = []
pp_1 = []
pp_2 = []
p1_dist = []
p2_dist = []
com = []
m = []
density = []
pv_1 = []
pv_2 = []
force_1 = []
force_2 = []
it = 0
lu = init_pf['LengthUnits']
Rsun = 6.96e10  #cm in a solar radius
TE = []

for pf in ts:
    time_val = pf.current_time * pf["years"]
    times.append(time_val)
    
    dd = pf.h.all_data()        
    particles = pf.h.find_particles_by_type(11)     #Passy's PM particle is type 11
    '''
    TE_full = dd['ThermalEnergy']*dd['CellMass']
    TE_val = sum(TE_full)
    TE.append(TE_val)
    
    '''
    #Then calculate separation using pythagoras
    pp1 = [dd['particle_position_x'][0], dd['particle_position_y'][0], dd['particle_position_z'][0]]
    pp2 = [dd['particle_position_x'][1], dd['particle_position_y'][1], dd['particle_position_z'][1]]
    pp_1.append(pp1)
    pp_2.append(pp2)
    separation = distance(pp1, pp2)
    sep.append((separation*lu)/Rsun)
    pc1 = distance(pp1, [0.5, 0.5, 0.5])
    p1_dist.append((pc1*lu)/Rsun)
    pc2 = distance(pp2, [0.5, 0.5, 0.5])/Rsun
    p2_dist.append((pc2*lu)/Rsun)
    
    pv1 = [dd['particle_velocity_x'][0], dd['particle_velocity_y'][0], dd['particle_position_z'][0]]
    pv2 = [dd['particle_velocity_x'][1], dd['particle_velocity_y'][1], dd['particle_position_z'][1]]
    pv_1.append(pv1)
    pv_2.append(pv2)
    
    #Make sphere between particles
    center = [(pp1[0] + pp2[0])/2., (pp1[1] + pp2[1])/2., (pp1[2] + pp2[2])/2.]
    sp = pf.h.sphere(center, (separation/2.))
    mass = sp.quantities['TotalQuantity']('CellMassMsun')[0]# - sp.quantities['TotalQuantity']('ParticleMassMsun')[0]
    m.append(mass)
    com_val = sp.quantities['CenterOfMass']()
    den = sp['Density']
    den_val = (sum(den))/(len(den))
    density.append(den_val)
    print "Initial:", sep[0], " separation:", separation
    '''
    if it ==0:
        force = (-6.67259e-8*(dd['ParticleMassMsun'][0] + dd['ParticleMassMsun'][1])*1.9891e+33)/((separation*lu)**2.)
        force_1.append(force)
        force_2.append(force)
    else:
        fx1 = (pv_1[it][0] - pv_1[it-1][0])/(times[it] - times[it-1])
        fy1 = (pv_1[it][1] - pv_1[it-1][1])/(times[it] - times[it-1])
        fz1 = (pv_1[it][2] - pv_1[it-1][2])/(times[it] - times[it-1])
        f1 = -((fx1**2 + fy1**2 + fz1**2)**0.5)*pm1
        force_1.append(f1)
        fx2 = (pv_2[it][0] - pv_2[it-1][0])/(times[it] - times[it-1])
        fy2 = (pv_2[it][1] - pv_2[it-1][1])/(times[it] - times[it-1])
        fz2 = (pv_2[it][2] - pv_2[it-1][2])/(times[it] - times[it-1])
        f2 = -((fx2**2 + fy2**2 + fz2**2)**0.5)*pm2
        force_2.append(f2)
    '''
    f.write(str(time_val) + ',' + str((separation*lu)/Rsun) + '\n')
    #f.write(str(time_val) + ',' + str(mass) + ',' + str(den_val) + '\n')

for i in range(len(sep)):
    print times[i], sep[i]

'''
# Plot the time arrays. Taken from an example on the internet: http://matplotlib.org/examples/axes_grid/demo_parasite_axes2.html)
plt.figure(figsize=(12,8), dpi=100)
host = host_subplot(111, axes_class=AA.Axes)
par1 = host.twinx()
host.set_xlabel("Time (years)")
host.set_ylabel("Separation ($R_\odot$)")
par1.set_ylabel("Mass ($M_\odot$)")
p1, = host.plot(times, sep, label="Separation")
#p2, = host.plot(times, sepave, label = "Average separation")
p3, = par1.plot(times, m, label="Mass within orbit")
par1.set_ylim(0.0, 0.00012)
#host.set_ylim([14, 21])
host.legend()
host.axis["left"].label.set_color(p1.get_color())
par1.axis["right"].label.set_color(p3.get_color())
'''
print 'building plot'
plt.clf()
plt.plot(times, sep, label='separation')
#plt.plot(times, force_1, label='0.6$M_\odot$')
#plt.plot(times, force_2, label='0.392$M_\odot$')
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
print 'labeling axis'
plt.xlabel("Time ($years$)")
plt.ticklabel_format(style='sci', axis='y')
plt.ylabel("Separation ($R_\odot$)")
#plt.legend()
#plt.tight_layout()
plt.savefig('separation.png')
