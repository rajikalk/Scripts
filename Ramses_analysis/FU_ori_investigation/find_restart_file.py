import numpy as np
import pickle
import glob
from pyramses import rsink
import sys
import os
import yt
import yt.units
from yt.units import g, s, cm, Lsun
import csv

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("-sink", "--sink_id", type=int, default=None)
    parser.add_argument("-end_time", "--end_time_val", default=75000, type=float)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================
args = parse_inputs()
target_times = {"17": 83816.50575862321, "51": 48419.8594782825, "101": 105290.43936672136, "103": 100963.20291334442, "177": 149168.43773802178, "221": 63131.04522576899, "262": 36958.441396587965, "272": 45604.98441964608}

path = sys.argv[1]
save_dir = sys.argv[2]
info_files = sorted(glob.glob(path+"*/info*.txt"))
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
    
loaded_sink_data = rsink(datadir=path, all=True)
for tag in target_times.keys():
    print('Searching for restart for sink', tag)
    sink_ind = int(tag)
    sink_form_time = 0
    counter = -1
    for sink_data in loaded_sink_data:
        counter = counter + 1
        if len(sink_data['tcreate']) > sink_ind:
            if sink_form_time == 0:
                sink_form_time = sink_data['tcreate'][sink_ind]*units['time_unit'].in_units('yr')
                    
            time_val = sink_data['snapshot_time']*units['time_unit'].in_units('yr') - sink_form_time
            if time_val >= yt.YTQuantity(target_times[tag], 'yr'):
                trunc_ind = counter
                prev_flush_time = sink_data['tflush']
                
                #binary search for file
                low_ind = 0
                high_ind = len(info_files) - 1
                mid_ind = int(len(info_files)/2)
                
                found = False
                while found == False:
                    with open(info_files[mid_ind], "r") as f:
                        reader = csv.reader(f, delimiter=' ')
                        for row in reader:
                            if row[0] == 'time':
                                curr_time = eval(row[-1])
                                break
                    f.close()
                    
                    if curr_time > prev_flush_time:
                        #if current time higher than target, then get mid index between the lowest and current ind
                        high_ind = mid_ind
                        mid_ind = int(np.mean([low_ind, mid_ind]))
                    elif curr_time < prev_flush_time:
                        #if current time lower than target, then get mid index between the current and highest ind
                        low_ind = mid_ind
                        mid_ind = int(np.mean([mid_ind, high_ind]))
                    else:
                        import pdb
                        pdb.set_trace()
                    
                    if low_ind == mid_ind:
                        if curr_time > prev_flush_time:
                            restart_file = info_files[mid_ind]
                            loaded_sink_data = loaded_sink_data[counter:]
                            found = True
                    elif high_ind == mid_ind:
                        import pdb
                        pdb.set_trace()
                print("Restart file for sink", sink_ind, "is", restart_file)
                break
