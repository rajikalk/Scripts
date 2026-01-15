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
from matplotlib.colors import LogNorm
cm = plt.cm.get_cmap('seismic')

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("-make_pickles", "--make_pickle_files", type=str, default="True")
    parser.add_argument("-make_plots", "--make_plot_figures", type=str, default="True")
    parser.add_argument("-sphere_radius", "--sphere_radius_cells", type=float)
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

if args.make_pickle_files == "True":
    #File files
    usable_files = sorted(glob.glob(input_dir+"*/info*.txt"))

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
    del units_override['density_unit']
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
    for fn in yt.parallel_objects(usable_files):
        ds = yt.load(fn, units_override=units_override)
        dd = ds.all_data()
        file_time = ds.current_time.in_units('yr')
        time_arr.append(file_time)
        
        #calculate separation to every other sink
        center_pos = yt.YTArray([dd['sink_particle_posx'][sink_id], dd['sink_particle_posy'][sink_id], dd['sink_particle_posz'][sink_id]])
        dx = dd['sink_particle_posx'] - center_pos[0]
        dy = dd['sink_particle_posy'] - center_pos[1]
        dz = dd['sink_particle_posz'] - center_pos[2]
        separations = np.sqrt(dx**2 + dy**2 + dz**2)
        close_sinks = np.where(separations.in_units('au')<10000)[0]
        if len(np.where(separations.in_units('au')<10000)[0]) > 1:
            
        import pdb
        pdb.set_trace()
