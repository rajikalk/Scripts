#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
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
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def has_sinks(ds):
    '''
    Checks particle file to see if particles exists, or tries the plot file.
    '''
    dd = ds.all_data()
    if len(dd['sink_particle_tag'][myf.get_centred_sink_id():].astype(int)) != 0:
        del dd
        return True
    else:
        del dd
        return False
        
def projected_vector(vector, proj_vector):
    """
    Calculates the projection of vecter projected onto vector
    """
    vector_units = vector.units
    if len(proj_vector)>3:
        #Calc vector.proj
        v_dot_pv = vector.T[0]*proj_vector.T[0] + vector.T[1]*proj_vector.T[1] + vector.T[2]*proj_vector.T[2]
        pv_dot_pv = proj_vector.T[0]**2 + proj_vector.T[1]**2 + proj_vector.T[2]**2
        proj_v_x = (v_dot_pv/pv_dot_pv)*proj_vector.T[0]
        proj_v_y = (v_dot_pv/pv_dot_pv)*proj_vector.T[1]
        proj_v_z = (v_dot_pv/pv_dot_pv)*proj_vector.T[2]
        proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T)
    else:
        proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
        proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
        proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
        proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
    return proj_v

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

sys.stdout.flush()
CW.Barrier()

#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}
units_override.update({"mass_unit":(2998,"Msun")})
units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm').value # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s').value         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3').value  # 2998 Msun / (4 pc)^3
mym.set_units(units_override)
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
if rank == 0:
    print("CENTERED SINK ID:", sink_id)
myf.set_centred_sink_id(sink_id)
sink_form_time = dd['sink_particle_form_time'][sink_id]
start_file = mym.find_files([0.0], usable_files, sink_form_time,sink_id)
usable_files = usable_files[usable_files.index(start_file[0]):]


dx_min = np.min(dd['dx'].in_units('au'))
sphere_radius = 4*dx_min
del dd
    
sys.stdout.flush()
CW.Barrier()

for fn in yt.parallel_objects(usable_files, njobs=int(size/6)):
    ds = yt.load(fn, units_override=units_override)
    dd = ds.all_data()
    #Get secondary position
    
    particle_position = yt.YTArray([dd['sink_particle_posx'][sink_id], dd['sink_particle_posy'][sink_id], dd['sink_particle_posz'][sink_id]])
    particle_velocity = yt.YTArray([dd['sink_particle_velx'][sink_id], dd['sink_particle_vely'][sink_id], dd['sink_particle_velz'][sink_id]])
    measuring_sphere = ds.sphere(particle_position.in_units('au'), sphere_radius)
    
    #Let's measure the angular momentum vector.
    sph_dx = measuring_sphere['x'].in_units('cm') - particle_position[0].in_units('cm')
    sph_dy = measuring_sphere['y'].in_units('cm') - particle_position[1].in_units('cm')
    sph_dz = measuring_sphere['z'].in_units('cm') - particle_position[2].in_units('cm')
    sph_radial_vector = yt.YTArray([sph_dx, sph_dy, sph_dz]).T
    
    sph_dvx = measuring_sphere['velocity_x'].in_units('cm/s') - particle_velocity[0].in_units('cm/s')
    sph_dvy = measuring_sphere['velocity_y'].in_units('cm/s') - particle_velocity[1].in_units('cm/s')
    sph_dvz = measuring_sphere['velocity_z'].in_units('cm/s') - particle_velocity[2].in_units('cm/s')
    sph_velocity_vector = yt.YTArray([sph_dvx, sph_dvy, sph_dvz]).T
    sph_specific_ang = yt.YTArray(np.cross(sph_radial_vector, sph_velocity_vector), 'cm**2/s')
    sph_ang = (measuring_sphere['mass'].in_units('g')*sph_specific_ang.T).T
    sph_total_ang = yt.YTArray([np.sum(sph_ang.T[0]), np.sum(sph_ang.T[1]), np.sum(sph_ang.T[2])])
    sph_total_ang_mag = np.sqrt(np.sum(sph_total_ang**2))
    sph_total_ang_unit = sph_total_ang/sph_total_ang_mag
    
    import pdb
    pdb.set_trace()
        
    del dd
