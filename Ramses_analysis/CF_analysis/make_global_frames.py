#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_fields_short as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False and rank==0:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

pickle_files = glob.glob(input_dir+'movie_frame*.pkl')








#Set some plot variables independant on data files
#File files
sink_files = sorted(glob.glob(input_dir+"output*/*.dat"))
files = sorted(glob.glob(input_dir+"*/info*.txt"))
rm_files = []
for info_name in files:
    sink_file = info_name.split('info')[0]+'stars_output.dat'
    if sink_file not in sink_files:
        rm_files.append(info_name)
for rm_file in rm_files:
    files.remove(rm_file)
del sink_files

print('Collected files on rank', rank)

sys.stdout.flush()
CW.Barrier()

#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = save_dir.split('G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})

print('Calculated the units on rank', rank)
sys.stdout.flush()
CW.Barrier()

usable_files = files
del files
sys.stdout.flush()
CW.Barrier()
'''
rit = -1
for fn_it in range(len(usable_files)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
'''
for fn_it in yt.parallel_objects(range(len(usable_files)), njobs=int(size/24)):
    pickle_file = save_dir + "movie_frame_" + ("%06d" % fn_it) + ".pkl"
    if os.path.exists(pickle_file) == False:
        fn = usable_files[fn_it]
        print('Making projection of file', fn, 'on rank', rank)
        ds = yt.load(fn, units_override=units_override)
        fn
        
        time_real = ds.current_time.in_units('yr')
        time_val = np.round(time_real.in_units('yr'))
        del time_real
        
        proj = yt.ProjectionPlot(ds, "z", ("ramses", "Density"), width=units['length_unit'], method='integrate')
        proj_array = np.array(proj.frb.data[("ramses", "Density")]/units['length_unit'].in_units('cm'))
        image = proj_array*units['density_unit'].in_units('g/cm**3')
        del proj
        del proj_array

        dd = ds.all_data()
        particle_positions = yt.YTArray(np.array([dd['sink_particle_posx'], dd['sink_particle_posy'], dd['sink_particle_posz']]).T, dd['sink_particle_posx'].units)
        particle_masses = dd['sink_particle_mass']
        del ds
        del dd
        

        if np.remainder(rank,24) == 0:
            file = open(pickle_file, 'wb')
            #pickle.dump((X, Y, proj_array, time_val, particle_positions, particle_masses), file)
            pickle.dump((image, time_val, particle_positions, particle_masses), file)
            file.close()
            print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
        del image
        del particle_positions
        del particle_masses
        del time_val
            
print('FINISHED MAKING YT PROJECTIONS ON RANK', rank)
