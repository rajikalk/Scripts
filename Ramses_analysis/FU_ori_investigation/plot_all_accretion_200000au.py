import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from pyramses import rsink
import sys
import os
import yt

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================
args = parse_inputs()

path = sys.argv[1]
save_dir = sys.argv[2]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == "False":
    os.makedirs(save_dir)
    
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

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

print("Reading particle data")
loaded_sink_data = rsink(datadir=path, all=True)
updating = False

if args.sink_number == None:
    last_n = int(sorted(glob.glob(path+"output*"))[-1].split("_")[-1])
    stars_output_file = path + 'output_'+("%05d" % last_n)+'/stars_output.dat'
    while os.path.exists(stars_output_file) == False:
        last_n = last_n - 1
        stars_output_file = path + 'output_'+("%05d" % last_n)+'/stars_output.dat'
    loaded_sink_data_last = rsink(last_n, datadir=path)
    sink_ind = np.argmin(loaded_sink_data_last['u'])
else:
    sink_ind = args.sink_number
    
for sink_data in loaded_sink_data:
    if len(sink_data['u']) > sink_ind:
       target_sink_formation_location = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
       break
    
plotted_sinks = []
for sink_data in loaded_sink_data:
    if len(sink_data['u']) < sink_ind:
        import pdb
        pdb.set_trace()
    if len(sink_data['u']) > sink_ind:
        tags = np.arange(len(sink_data['u']))[sink_ind-1:sink_ind+1]
        for tag in tags:
            if tag not in nearby_accretion['particle_tag']:
                nearby_accretion['particle_tag'].append(tag)
        pos_prim = yt.YTArray(np.array([sink_data['x'][sink_ind-1], sink_data['y'][sink_ind-1], sink_data['z'][sink_ind-1]])*units['length_unit'].in_units('au'), 'au')
        pos_second = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
        separation = np.sqrt(np.sum((pos_second - pos_prim)**2))
        nearby_accretion['separation'].append(separation)
        if sink_form_time == 0:
            sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
        time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
        #if len(particle_tags) == 1:
        nearby_accretion['time'].append(time_val)
        nearby_accretion['mass'].append(yt.YTArray(sink_data['m'][sink_ind-1:sink_ind+1]*units['mass_unit'].in_units('msun'), 'msun'))
        
        d_mass = sink_data['dm'][sink_ind-1:sink_ind+1]*units['mass_unit'].in_units('msun')
        d_time = (sink_data['snapshot_time'] - sink_data['tflush'])*units['time_unit'].in_units('yr')
        acc_val = d_mass/d_time
        acc_val[np.where(acc_val == 0)[0]]=1.e-12
        nearby_accretion['mdot'].append(yt.YTArray(acc_val, 'msun/yr'))
#write lastest pickle
file = open(save_dir+'nearby_accretion.pkl', 'wb')
pickle.dump((nearby_accretion, counter, sink_ind, sink_form_time), file)
file.close()
os.system('cp '+save_dir+'nearby_accretion.pkl '+save_dir+'nearby_accretion_tmp.pkl ')
print('read', counter, 'snapshots of sink particle data, and saved pickle')
