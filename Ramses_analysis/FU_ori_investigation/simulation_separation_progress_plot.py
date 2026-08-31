import numpy as np
from pyramses import rsink
import glob
import yt
import yt.units
from yt.units import g, s, cm, Lsun
import matplotlib.pyplot as plt
import os


#================================================================================
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override

last_n = int(sorted(glob.glob("data/output*"))[-1].split("_")[-1])
stars_output_file = 'data/output_'+("%05d" % last_n)+'/stars_output.dat'
path = os.getcwd()
sink_ind = int(path.split('/Sink_')[1].split('/')[0])
'''
loaded_sink_data_last = rsink(last_n, datadir='data/')
del last_n
sink_ind = np.argmin(loaded_sink_data_last['u'])
import pdb
pdb.set_trace()
del loaded_sink_data_last
'''
print("Sink_ind =", sink_ind)
loaded_sink_data = rsink(datadir='data/', all=True)
    
particle_data = {}
particle_data.update({'time':np.array([])})
particle_data.update({'mass':np.array([])})
particle_data.update({'mdot':np.array([])})
particle_data.update({'separation':np.array([])})
particle_data.update({'closest_sink':np.array([])})
particle_data.update({'closest_mass':np.array([])})
particle_data.update({'closest_mdot':np.array([])})
        
counter = 0
sink_form_time = 0

print("starting to update pickles, current counter=", counter)
while len(loaded_sink_data)>0:
    sink_data = loaded_sink_data[0]
    loaded_sink_data = loaded_sink_data[1:]
    counter = counter + 1
    if sink_form_time == 0:
        sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
    time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
    particle_data['time'] = np.append(particle_data['time'], time_val)
    
    particle_mass = yt.YTArray(sink_data['m'][np.array([sink_ind])]*units['mass_unit'].in_units('msun'), 'msun')
    particle_data['mass'] = np.append(particle_data['mass'], particle_mass)
    
    d_mass = sink_data['dm'][np.array([sink_ind])]*units['mass_unit'].in_units('msun')
    d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    if acc_val == 0:
        acc_val =1.e-12
    particle_data['mdot'] = np.append(particle_data['mdot'], acc_val)

    position = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
    dx = sink_data['x']*units['length_unit'].in_units('au') - position[0]
    dy = sink_data['y']*units['length_unit'].in_units('au') - position[1]
    dz = sink_data['z']*units['length_unit'].in_units('au') - position[2]
    sep = np.sqrt(dx**2 + dy**2 + dz**2)
    closest_ind = np.argsort(sep)[1]
    particle_data['closest_sink'] = np.append(particle_data['closest_sink'], closest_ind)
    separation = sep[closest_ind]
    particle_data['separation'] = np.append(particle_data['separation'], separation)
    
    other_mass = yt.YTArray(sink_data['m'][np.array([closest_ind])]*units['mass_unit'].in_units('msun'), 'msun')
    d_mass = sink_data['dm'][np.array([closest_ind])]*units['mass_unit'].in_units('msun')
    d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
    acc_val = d_mass/d_time
    if acc_val == 0:
        acc_val =1.e-12
    particle_data['closest_mass'] = np.append(particle_data['closest_mass'], other_mass)
    particle_data['closest_mdot'] = np.append(particle_data['closest_mdot'], acc_val)
    
plt.clf()
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
axs.flatten()[0].plot(particle_data['time'], particle_data['mass'])
axs.flatten()[0].plot(particle_data['time'], particle_data['closest_mass'])
axs.flatten()[0].set_ylabel('Mass (M$_\odot$)')
axs.flatten()[0].set_ylim(bottom=0)
axs.flatten()[1].semilogy(particle_data['time'], particle_data['mdot'])
axs.flatten()[1].semilogy(particle_data['time'], particle_data['closest_mdot'])
axs.flatten()[1].set_ylabel('Accretion rate (M$_\odot$/yr)')
axs.flatten()[1].set_ylim(bottom=1.e-8)
axs.flatten()[2].semilogy(particle_data['time'], particle_data['separation'])
axs.flatten()[2].set_ylabel('Separation (AU)')
axs.flatten()[2].set_xlabel('Time (yr)')
res_limit = 2*3.15
if 'Level_19' in path:
    res_limit = 2*1.57
if 'Level_20' in path:
    res_limit = 2*0.79
axs.flatten()[2].axhline(y=res_limit, ls=':')
axs.flatten()[2].set_xlim([particle_data['time'][0], particle_data['time'][-1]])
figname= "Sink_"+str(sink_ind)
if 'Level' in path:
    level_it = path.split('Level_')[-1].split('/')[0]
    figname = figname + "_Level_" + str(level_it)
figname = figname + '_progress_plot.png'
plt.savefig(figname)
