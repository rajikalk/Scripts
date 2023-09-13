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
myf.set_com_pos_use_part(True)
myf.set_com_vel_use_gas(False)
myf.set_com_vel_use_part(True)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
#find start file
usable_files = mym.find_files([0], files, sink_form_time, sink_id, verbatim=False)
start_ind = files.index(usable_files[0])
files = files[start_ind+1:]
del dd

if args.make_pickle_files == 'True':

    sink_dict = {'time':[], 'mass':[], 'mdot':[], 'max_outflow_speed':[], 'mean_density':[]}

    for fn in yt.parallel_objects(files):
        pickle_file = 'burst_analysys_sink_'+str(sink_id)+'_'+str(rank)+'.pkl'
        ds = yt.load(fn, units_override=units_override)
        dd = ds.all_data()
        
        sink_dict['time'].append(ds.current_time.in_units('yr') - sink_form_time)
        sink_dict['mass'].append(dd['sink_particle_mass'][sink_id].in_units('msun'))
        sink_dict['mdot'].append(dd['sink_particle_accretion_rate'][sink_id].in_units('msun/yr'))
        #Define box:
        #center_pos = yt.YTArray([dd['sink_particle_posx'][sink_id].in_units('au'), dd['sink_particle_posy'][sink_id].in_units('au'), dd['sink_particle_posz'][sink_id].in_units('au')])
        #center_vel = yt.YTArray([dd['sink_particle_velx'][sink_id].in_units('cm/s'), dd['sink_particle_vely'][sink_id].in_units('cm/s'), dd['sink_particle_velz'][sink_id].in_units('cm/s')])
        
        center_pos = dd['Center_Position'].in_units('au')
        center_vel = dd['Center_Velocity'].in_units('cm/s')
        
        sph = ds.sphere(center_pos, (20, "au"))
        L_vec = yt.YTArray([np.sum(sph['Angular_Momentum_x']), np.sum(sph['Angular_Momentum_y']), np.sum(sph['Angular_Momentum_z'])])
        L_mag = np.sqrt(np.sum(L_vec**2))
        L_norm = L_vec/L_mag
        
        disk = ds.disk(center_pos, L_norm, (20, "au"), (20, "au"))
        sep_vector = yt.YTArray([disk['x'].in_units('cm')-center_pos[0].in_units('cm'), disk['y'].in_units('cm')-center_pos[1].in_units('cm'), disk['z'].in_units('cm')-center_pos[2].in_units('cm')]).T
        gas_vel = yt.YTArray([disk['Corrected_velx'], disk['Corrected_vely'], disk['Corrected_velz']]).in_units('cm/s').T
        proj_vectors = myf.projected_vector(gas_vel, sep_vector)
        vel_dot = np.diag(np.dot(gas_vel, proj_vectors.T))
        outflow_inds = np.where(vel_dot>1)[0]
        max_outflow_vel = np.max(disk['Corrected_vel_mag'][outflow_inds])
        sink_dict['max_outflow_speed'].append(max_outflow_vel.in_units('km/s'))
        sink_dict['mean_density'].append(np.mean(disk['density']))
        
        file = open(pickle_file, 'wb')
        pickle.dump((sink_dict), file)
        file.close()
        
        print('updated', pickle_file, 'for file', files.index(fn), 'of', len(files))
    
