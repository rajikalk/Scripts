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
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472

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
    if len(sink_data['u']) > sink_ind:
        target_sink_location = yt.YTArray(np.array([sink_data['x'][sink_ind], sink_data['y'][sink_ind], sink_data['z'][sink_ind]])*units['length_unit'].in_units('au'), 'au')
        dx = sink_data['x']*units['length_unit'].in_units('au') - target_sink_location[0]
        dy = sink_data['y']*units['length_unit'].in_units('au') - target_sink_location[1]
        dz = sink_data['z']*units['length_unit'].in_units('au') - target_sink_location[2]
        sep = np.sqrt(dx**2 + dy**2 + dz**2)
        close_sinks = np.where(sep<20000)[0]
    else:
        dx = sink_data['x']*units['length_unit'].in_units('au') - target_sink_formation_location[0]
        dy = sink_data['y']*units['length_unit'].in_units('au') - target_sink_formation_location[1]
        dz = sink_data['z']*units['length_unit'].in_units('au') - target_sink_formation_location[2]
        sep = np.sqrt(dx**2 + dy**2 + dz**2)
        close_sinks = np.where(sep<20000)[0]
    for close_sink in close_sinks:
        if close_sink not in plotted_sinks:
            time_arr = []
            acc_arr = []
            mass_arr = []
            for sink_data_acc in loaded_sink_data:
                if len(sink_data_acc['u']) > close_sink:
                    time_val = sink_data_acc['snapshot_time']*units['time_unit'].in_units('yr')
                    mass_val = sink_data_acc['m'][close_sink]*units['mass_unit'].in_units('msun')
                    d_mass = sink_data_acc['dm'][close_sink]*units['mass_unit'].in_units('msun')
                    d_time = (sink_data_acc['snapshot_time'] - sink_data_acc['tflush'])*units['time_unit'].in_units('yr')
                    acc_val = d_mass/d_time
                    time_arr.append(time_val)
                    acc_arr.append(acc_val)
                    mass_arr.append(mass_val)
            time_arr = yt.YTArray(time_arr) - time_arr[0]
            plt.clf()
            fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(single_col_width, single_col_width*1.5), sharex=True)
            axs[0].plot(time_arr, mass_arr)
            axs[1].semilogy(time_arr, acc_arr)
            axs[0].set_title("Sink "+str(close_sink))
            axs[1].set_xlabel("Simulation time")
            axs[0].set_ylabel("Mass (Msun)")
            axs[1].set_ylabel("Accretion rate (Msun/yr)")
            plt.savefig("Sink_"+str(close_sink)+".png", bbox_inches='tight')
            plotted_sinks.append(close_sink)
            print("plotted accretion history for sink", close_sink+1)
