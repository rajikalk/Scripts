#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import my_ramses_fields as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("-make_pickles", "--make_pickle_files", help="do you want to update pickles?", default='True', type=str)
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

#Set some plot variables independant on data files
sys.stdout.flush()
CW.Barrier()

#File files
files = sorted(glob.glob(input_dir+"*/info*.txt"))

sys.stdout.flush()
CW.Barrier()

#Define units to override:
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
mym.set_units(units_override)


#find sink particle to center on and formation time
ds = yt.load(files[-1], units_override=units_override)
#try:
dd = ds.all_data()
if args.sink_number == None:
    sink_id = np.argmin(dd['sink_particle_speed'])
else:
    sink_id = args.sink_number
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_com_pos_use_gas(False)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
del dd

if args.make_pickle_files == 'True':

    sink_dict = {'time':[], 'mass':[], 'mdot':[], 'max_outflow_speed':[], 'mean_density':[]}

    file_int = -1
    for fn in yt.parallel_objects(files):
        ds = yt.load(fn, units_override=units_override)
        if ds.current_time.in_units('yr') > sink_form_time:
            dd = ds.all_data()
            
            sink_dict['time'].append(ds.current_time.in_units('yr'))
            sink_dict['mass'].append(dd['sink_particle_mass'][sink_id].in_units('msun'))
            sink_dict['mdot'].append(dd['sink_particle_accretion_rate'][sink_id].in_units('msun/yr'))
            #Define box:
            center_pos = yt.YTArray([dd['sink_particle_posx'][sink_id].in_units('au'), dd['sink_particle_posy'][sink_id].in_units('au'), dd['sink_particle_posz'][sink_id].in_units('au')])
            
            sph = ds.sphere(center_pos, (20, "kpc"))
            
            
            
        
