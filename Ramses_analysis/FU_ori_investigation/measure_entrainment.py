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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-sink_id", "--sink_number", help="Which sink do you want to measure around? default is the sink with lowest velocity", type=int, default=None)
    parser.add_argument("-make_pickles", "--make_pickle_files", type=str, default="True")
    parser.add_argument("-make_plots", "--make_plot_figures", type=str, default="True")
    parser.add_argument("-sphere_radius", "--sphere_radius_cells", type=float, default=4)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
        
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

if args.make_pickle_files == "True":
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
    dx_min = np.min(dd['dx'].in_units('au'))
    sphere_radius = args.sphere_radius_cells*dx_min
    r_acc = 4*dx_min
    if args.sink_number == None:
        sink_id = np.argmin(dd['sink_particle_speed'])
    else:
        sink_id = args.sink_number
    if rank == 0:
        print("CENTERED SINK ID:", sink_id)
    #myf.set_centred_sink_id(sink_id)
    sink_form_time = dd['sink_particle_form_time'][sink_id]
    del dd
        
    sys.stdout.flush()
    CW.Barrier()

    file_int = -1
    no_files = len(usable_files)
    for fn in yt.parallel_objects(usable_files, njobs=int(size)):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
            if usable_files[file_int] == usable_files[file_int-1]:
                os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        make_pickle = False
        pickle_file = save_dir + "movie_frame_" + ("%06d" % file_int + ".pkl")
        if os.path.isfile(pickle_file) == False:
            make_pickle = True
        if make_pickle:
            ds = yt.load(fn, units_override=units_override)
            time_val = ds.current_time.in_units('yr') - sink_form_time
            dd = ds.all_data()
            
            #Get secondary position
            particle_position = yt.YTArray([dd['sink_particle_posx'][sink_id], dd['sink_particle_posy'][sink_id], dd['sink_particle_posz'][sink_id]])
            particle_velocity = yt.YTArray([dd['sink_particle_velx'][sink_id], dd['sink_particle_vely'][sink_id], dd['sink_particle_velz'][sink_id]])
            measuring_sphere = ds.sphere(particle_position.in_units('au'), sphere_radius)
            print("Got particle position and velocity")
            
            #Let's measure the angular momentum vector.
            import pdb
            pdb.set_trace()

sys.stdout.flush()
CW.Barrier()

