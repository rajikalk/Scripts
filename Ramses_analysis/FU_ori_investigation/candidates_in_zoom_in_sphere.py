#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import my_ramses_fields_short as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.pyplot as plt
#plt.rcParams['figure.dpi'] = 300

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======
rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

#File files
usable_files = sorted(glob.glob(input_dir+"*/info*.txt"))
if size == 1:
    import pdb
    pdb.set_trace()

sys.stdout.flush()
CW.Barrier()

#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
mym.set_units(units_override)
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3").in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
if rank == 0:
    print("set units")

#find sink particle to center on and formation time
#del units_override['density_unit']
ds = yt.load(usable_files[-1], units_override=units_override)
if rank == 0:
    print("Doing initial ds.all_data() load")
dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
    
#Start iterating through files
time_arr = []
nearest_sinks = []
max_inds = 0
for fn in yt.parallel_objects(usable_files):
    ds = yt.load(fn, units_override=units_override)
    dd = ds.all_data()
    
    #calculate separation to every other sink
    if len(dd['sink_particle_posx'])>sink_id:
        file_time = ds.current_time.in_units('yr')
        time_arr.append(file_time)
        center_pos = yt.YTArray([dd['sink_particle_posx'][sink_id], dd['sink_particle_posy'][sink_id], dd['sink_particle_posz'][sink_id]])
        dx = dd['sink_particle_posx'] - center_pos[0]
        dy = dd['sink_particle_posy'] - center_pos[1]
        dz = dd['sink_particle_posz'] - center_pos[2]
        separations = np.sqrt(dx**2 + dy**2 + dz**2)
        close_sinks = np.where(separations.in_units('au')<10000)[0]
        nearest_sinks.append(close_sinks)
        if len(nearest_sinks) > max_inds:
            max_inds = len(nearest_sinks)

        file = open('Candidates_in_zoom_sphere_' + str(rank) + '.pkl', 'wb')
        pickle.dump((time_arr, nearest_sinks, max_inds), file)
        file.close()
        if len(close_sinks) > 1:
            print('Found ' + str(close_sinks) + ' within 10000au of ' + str(sink_id) + ' for ' + fn)
        else:
            print('No sinks found within 10000au of ' + str(sink_id) + ' for ' + fn)
        
#Compile pickles
Candidate_pickles = glob.glob('*_sphere_*.pkl')
time_arr = []
nearest_sinks = []
max_inds = 0
for pick_file in Candidate_pickles:
    file = open(pick_file, 'rb')
    nearest_data = pickle.load(file)
    file.close()
    
    time_arr = time_arr + nearest_data[0]
    nearest_sinks = nearest_sinks + nearest_data[1]
    #if nearest_data[2] > max_inds:
    #    max_inds = nearest_data[2]
 
#Square off the nearest sinks
for nearest_ind in range(len(nearest_sinks)):
    if len(nearest_sinks[nearest_ind]) < max_inds:
        add_length = max_inds - len(nearest_sinks[nearest_ind])
        nearest_sinks[nearest_ind] = np.append(nearest_sinks[nearest_ind], np.zeros(add_length)*np.nan)
    if len(nearest_sinks[nearest_ind]) == max_inds:
        print('All nearby sinks are', nearest_sinks[nearest_ind])

sorted_inds = np.argsort(time_arr)
time_arr = np.array(time_arr)[sorted_inds]
nearest_sinks = np.array(nearest_sinks)[sorted_inds]

file = open('Candidates_in_zoom_sphere.pkl', 'wb')
pickle.dump((time_arr, nearest_sinks, max_inds), file)
file.close()
