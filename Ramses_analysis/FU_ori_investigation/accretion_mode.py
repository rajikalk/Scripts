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
    if args.sink_number == None:
        sink_id = np.argmin(dd['sink_particle_speed'])
    else:
        sink_id = args.sink_number
    if rank == 0:
        print("CENTERED SINK ID:", sink_id)
    #myf.set_centred_sink_id(sink_id)
    sink_form_time = dd['sink_particle_form_time'][sink_id]
    start_file = mym.find_files([0.0], usable_files, sink_form_time,sink_id)
    usable_files = usable_files[usable_files.index(start_file[0]):]

    dx_min = np.min(dd['dx'].in_units('au'))
    sphere_radius = 4*dx_min
    del dd
        
    sys.stdout.flush()
    CW.Barrier()

    file_int = -1
    no_files = len(usable_files)
    for fn in yt.parallel_objects(usable_files, njobs=int(size/6)):
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
            del dd
            
            #Let's measure the angular momentum vector.
            sph_dx = measuring_sphere['x'].in_units('cm') - particle_position[0].in_units('cm')
            sph_dy = measuring_sphere['y'].in_units('cm') - particle_position[1].in_units('cm')
            sph_dz = measuring_sphere['z'].in_units('cm') - particle_position[2].in_units('cm')
            sph_radial_vector = yt.YTArray([sph_dx, sph_dy, sph_dz]).T
            sph_radial_vector_mag = np.sqrt(np.sum(sph_radial_vector**2, axis=1)).value
            sph_radial_vector_unit = yt.YTArray([sph_dx/sph_radial_vector_mag, sph_dy/sph_radial_vector_mag, sph_dz/sph_radial_vector_mag]).T
            
            sph_dvx = measuring_sphere['velocity_x'].in_units('cm/s') - particle_velocity[0].in_units('cm/s')
            sph_dvy = measuring_sphere['velocity_y'].in_units('cm/s') - particle_velocity[1].in_units('cm/s')
            sph_dvz = measuring_sphere['velocity_z'].in_units('cm/s') - particle_velocity[2].in_units('cm/s')
            sph_velocity_vector = yt.YTArray([sph_dvx, sph_dvy, sph_dvz]).T
            
            #Calculate radial velocity
            shape = np.shape(measuring_sphere['x'])
            radial_vel_vec = yt.YTArray(projected_vector(sph_velocity_vector.in_units('km/s'), sph_radial_vector_unit.in_units('au')).value, 'km/s')
            radial_vel_mag = np.sqrt(np.sum(radial_vel_vec**2, axis=1))
            radial_vel_unit = (radial_vel_vec.T/radial_vel_mag).T
            sign = np.diag(np.dot(sph_radial_vector_unit.in_units('au'), radial_vel_unit.T))
            sign = np.sign(sign)
            
            rv_mag = radial_vel_mag*sign
            rv_mag = np.reshape(rv_mag, shape)
            if np.inf in rv_mag.value or np.nan in rv_mag.value:
                rv_mag = yt.YTArray(np.nan_to_num(rv_mag.value), 'km/s')
                
            #Calcualte radial momentum
            radial_momentum = rv_mag.in_units('cm/s') * measuring_sphere['mass'].in_units('g')
            
            file = open(pickle_file, 'wb')
            pickle.dump((time_val, measuring_sphere['density'], radial_momentum), file)
            file.close()
            print("wrote file", pickle_file, "for file_int", file_int, "of", no_files)



sys.stdout.flush()
CW.Barrier()

if args.make_plot_figures == "True":
    import matplotlib.pyplot as plt
    #plt.rcParams['figure.dpi'] = 300
    from matplotlib.colors import LogNorm

    #Plotting
    pickle_files = sorted(glob.glob(save_dir + "movie_frame_*.pkl"))
    xmin = np.nan
    xmax = np.nan
    lin_thresh = np.nan
    rit = -1
    fit = 0
    while fit < len(pickle_files):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            plot_pickle = pickle_files[fit]
            file_name = save_dir + plot_pickle[:-3]+'jpg'
            if os.path.isfile(file_name) == False:
                file = open(plot_pickle, 'rb')
                time_val, density, radial_momentum = pickle.load(file)
                file.close()
                
                if np.isnan(xmin):
                    xmin = np.min(density)
                elif np.min(density) < xmin:
                    xmin = np.min(density)
                
                if np.isnan(xmax):
                    xmax = np.max(density)
                elif np.max(density) > xmax:
                    xmax = np.max(density)
                    
                if np.isnan(lin_thresh):
                    lin_thresh = np.min(np.abs(radial_momentum))
                elif np.min(np.abs(radial_momentum)) < lin_thresh:
                    lin_thresh = np.min(np.abs(radial_momentum))
                    
                #Plot figure
                plt.clf()
                plt.xscale('log')
                plt.yscale('symlog', linthresh=lin_thresh)
                plt.scatter(density, radial_momentum)
                plt.xlim([xmin.value,xmax.value])
                plt.xlabel('density (g/cm$^3$)')
                plt.ylabel('radial momentum (cm$\,$g/s)')

                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                print("Plotted", file_name, "for pickle", fit, "of", len(pickle_files))
        fit = fit + 1
