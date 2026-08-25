import numpy as np
from pyramses import rsink
import glob
import yt
import yt.units
from yt.units import g, s, cm, Lsun
import matplotlib.pyplot as plt

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
loaded_sink_data_last = rsink(last_n, datadir='data/')
del last_n
sink_ind = np.argmin(loaded_sink_data_last['u'])
del loaded_sink_data_last
print("Sink_ind =", sink_ind)
loaded_sink_data = rsink(datadir='data/', all=True)
    
particle_data = {}
particle_data.update({'time':np.array([])})
particle_data.update({'separation':np.array([])})
        
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

    position = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
    dx = sink_data['x']*units['length_unit'].in_units('au') - position[0]
    dy = sink_data['y']*units['length_unit'].in_units('au') - position[1]
    dz = sink_data['z']*units['length_unit'].in_units('au') - position[2]
    sep = np.sqrt(dx**2 + dy**2 + dz**2)
    particle_data['separation'] = np.append(particle_data['separation'], np.min(sep))
