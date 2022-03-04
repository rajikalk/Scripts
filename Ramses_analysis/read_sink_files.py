import csv
import numpy as np
import sys
import matplotlib.pyplot as plt
import yt
from pyramses import rsink
import pickle
import shutil
import os

#============================================================================================
def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-radius", "--proximity_radius", help="within what radius (in AU) do you want to save sink particle data?", type=float, default=10000.0)
    parser.add_argument("-def_sys", "--define_system", help="Is there a particular system that you would like to plot?", type=str, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-end_t", "--end_time", help="dow you want to add a y limit?", type=float, default=None)
    parser.add_argument("-plt_matches", "--plot_matched_times", help="do you want to plot to match times that you found?", default='False', type=str)
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

if args.simulation_density_id == 'G50':
    units_override.update({"mass_unit":(1500,"Msun")})
elif args.simulation_density_id == 'G200':
    units_override.update({"mass_unit":(6000,"Msun")})
elif args.simulation_density_id == 'G400':
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    units_override.update({"mass_unit":(2998,"Msun")})

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
momentum_unit = yt.YTQuantity(units_override['mass_unit'][0], units_override['mass_unit'][1]).in_units('g')*yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})


directory = sys.argv[1]
save_dir = sys.argv[2]
pickle_file = save_dir + 'particle_data.pkl'
sink_data = rsink(datadir=directory,all=True)
starting_sink_index = np.argmin(sink_data[-1]['u'])
n_particles_total = np.arange(len(sink_data[-1]['u']))[starting_sink_index:]

if os.path.exists(pickle_file):
    try:
        file_open = open(pickle_file, 'rb')
        particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
        file_open.close()
    except:
        shutil.copy(pickle_file.split('.pkl')[0]+'_tmp.pkl',pickle_file)
        file_open = open(pickle_file, 'rb')
        particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
        file_open.close()
else:
    particle_data = {}
    particle_data.update({'particle_tag':[]})
    particle_data.update({'time':[]})
    particle_data.update({'posx':[]})
    particle_data.update({'posy':[]})
    particle_data.update({'posz':[]})
    particle_data.update({'velx':[]})
    particle_data.update({'vely':[]})
    particle_data.update({'velz':[]})
    particle_data.update({'momx':[]})
    particle_data.update({'momy':[]})
    particle_data.update({'momz':[]})
    particle_data.update({'mass':[]})
    particle_data.update({'mdot':[]})
    particle_data.update({'speed':[]})
    particle_data.update({'age':[]})

    particle_data['particle_tag'] = n_particles_total
    sink_form_time = sink_data[-1]['tcreate'][starting_sink_index]*time_unit.in_units('yr').value
    init_line_counter = 0


start_it = init_line_counter
for sit in range(start_it,len(sink_data)):
    if len(sink_data[sit]['x'])>starting_sink_index:
        n_particles = np.arange(len(sink_data[-1]['u']))[starting_sink_index:]
        time_val = sink_data[sit]['snapshot_time']*time_unit.in_units('yr').value - sink_form_time
        particle_data['time'].append(time_val)
        particle_data['posx'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['posy'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['posz'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['velx'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['vely'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['velz'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['momx'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['momy'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['momz'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['mass'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['mdot'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['speed'].append(np.empty(np.shape(n_particles_total))*np.nan)
        particle_data['age'].append(np.empty(np.shape(n_particles_total))*np.nan)
        
        
        particle_data['posx'][-1][:len(n_particles)] = sink_data[sit]['x'][starting_sink_index:]*length_unit.in_units('AU').value
        particle_data['posy'][-1][:len(n_particles)] = sink_data[sit]['y'][starting_sink_index:]*length_unit.in_units('AU').value
        particle_data['posz'][-1][:len(n_particles)] = sink_data[sit]['z'][starting_sink_index:]*length_unit.in_units('AU').value
        particle_data['velx'][-1][:len(n_particles)] = sink_data[sit]['ux'][starting_sink_index:]*velocity_unit.in_units('cm/s').value
        particle_data['vely'][-1][:len(n_particles)] = sink_data[sit]['uy'][starting_sink_index:]*velocity_unit.in_units('cm/s').value
        particle_data['velz'][-1][:len(n_particles)] = sink_data[sit]['uz'][starting_sink_index:]*velocity_unit.in_units('cm/s').value
        particle_data['momx'][-1][:len(n_particles)] = sink_data[sit]['px'][starting_sink_index:]*momentum_unit.in_units('g*cm/s').value
        particle_data['momy'][-1][:len(n_particles)] = sink_data[sit]['py'][starting_sink_index:]*momentum_unit.in_units('g*cm/s').value
        particle_data['momz'][-1][:len(n_particles)] = sink_data[sit]['pz'][starting_sink_index:]*momentum_unit.in_units('g*cm/s').value
        particle_data['mass'][-1][:len(n_particles)] = sink_data[sit]['m'][starting_sink_index:]*mass_unit.in_units('Msun').value
        numerator = sink_data[sit]['dm'][starting_sink_index:]*mass_unit.in_units('Msun').value
        denominator = (sink_data[sit]['snapshot_time'] - sink_data[sit]['tflush'])*time_unit.in_units('yr').value
        particle_data['mdot'][-1][:len(n_particles)] = numerator/denominator
        particle_data['speed'][-1][:len(n_particles)] = sink_data[sit]['u'][starting_sink_index:]*velocity_unit.in_units('cm/s').value
        particle_data['age'][-1][:len(n_particles)] = (sink_data[sit]['snapshot_time']-sink_data[sit]['tcreate'][starting_sink_index:])*time_unit.in_units('yr').value
        
        if np.remainder(sit,10000) == 0:
            line_counter = sit
            file_open = open(pickle_file, 'wb')
            pickle.dump((particle_data, sink_form_time, line_counter),file_open)
            file_open.close()
            print("dumped pickle after line", line_counter)
            shutil.copy(pickle_file, pickle_file.split('.pkl')[0]+'_tmp.pkl')
            file_open = open(pickle_file, 'rb')
            particle_data, sink_form_time, line_counter = pickle.load(file_open)
            file_open.close()
        
        print("read data for snapshot", sit, "of", len(sink_data))
        
for key in list(particle_data.keys()):
    try:
        particle_data[key] = np.array(particle_data[key])
        particle_data[key] = particle_data[key][sorted_inds]
    except:
        continue
if len(particle_data['particle_tag']) == 2:
    particle_data.update({'separation':np.sqrt((particle_data['posx'].T[0] - particle_data['posx'].T[1])**2. + (particle_data['posy'].T[0] - particle_data['posy'].T[1])**2. + (particle_data['posz'].T[0] - particle_data['posz'].T[1])**2.)})
    zero_inds = np.where(particle_data['separation'] == 0)[0]
    particle_data['separation'][zero_inds] = np.nan
print("sorted data and save")
line_counter = sit
file_open = open(pickle_file, 'wb')
pickle.dump((particle_data, sink_form_time, line_counter),file_open)
file_open.close()
        
