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
    parser.add_argument("-sphere_radius", "--sphere_radius_cells", type=float, default=10)
    parser.add_argument("-cone_angle", "--entrainment_cone_angle", type=float, default=30)
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
    sphere_radius = yt.YTQuantity(args.sphere_radius_cells, 'au')
    cone_angle = yt.YTQuantity(args.entrainment_cone_angle, '')
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
            accretion_rate = dd['sink_particle_accretion_rate'][sink_id].in_units('msun/yr')
            del dd
            particle_speed = np.sqrt(np.sum(particle_velocity**2))
            particle_velocity_unit = particle_velocity/particle_speed
            del particle_speed
            measuring_sphere = ds.sphere(particle_position.in_units('au'), sphere_radius)
            print("Got particle position and velocity")
            
            #Let's calculate the entrainment vector.
            entrainment_vector = -1*particle_velocity_unit
            
            #Calcalate relative positions to particle
            dx = measuring_sphere['x'] - particle_position[0]
            dy = measuring_sphere['y'] - particle_position[1]
            dz = measuring_sphere['z'] - particle_position[2]
            gas_position =yt.YTArray([dx, dy, dz])
            del dx, dy, dz
            separations = np.sqrt(np.sum(gas_position**2, axis=0))
            gas_position_unit = (gas_position/separations).T
            del separations
            angles = np.rad2deg(np.arccos(np.dot(gas_position_unit, entrainment_vector)))
            cone_inds = np.where(angles<=cone_angle)[0]
            del angles
            
            #Save arrays
            cone_densities = measuring_sphere[('gas','Density')][cone_inds]
            
            #Relative_velocities
            dvx = measuring_sphere['x-velocity'][cone_inds].in_units('km/s') - particle_velocity[0]
            dvy = measuring_sphere['y-velocity'][cone_inds].in_units('km/s') - particle_velocity[1]
            dvz = measuring_sphere['z-velocity'][cone_inds].in_units('km/s') - particle_velocity[2]
            gas_velocity = yt.YTArray([dvx, dvy, dvz]).T
            del dvx, dvy, dvz
            radial_vel = projected_vector(gas_velocity, gas_position_unit[cone_inds])
            radial_sign = np.sign(np.diag(np.dot(gas_velocity, gas_position_unit[cone_inds].T)))
            radial_speed = np.sqrt(np.sum(radial_vel**2, axis=1)) * radial_sign
            del radial_vel, radial_sign
            radial_speed = yt.YTArray(radial_speed, 'km/s')
            radial_momentum = radial_speed * measuring_sphere['mass'][cone_inds].in_units('g')
            
            write_dict = {'time':time_val, 'mdot':accretion_rate, 'density':cone_densities, 'radial_speed':radial_speed, 'radial_momentum':radial_momentum}
            file = open(pickle_file, 'wb')
            pickle.dump((write_dict), file)
            file.close()
            print("wrote file", pickle_file, "for file_int", file_int, "of", no_files)

sys.stdout.flush()
CW.Barrier()

if args.make_plot_figures == "True":
    import matplotlib.pyplot as plt
    #plt.rcParams['figure.dpi'] = 300
    from matplotlib.colors import LogNorm
    cm = plt.cm.get_cmap('seismic')
    
    two_col_width = 7.20472 #inches
    single_col_width = 3.50394 #inches
    page_height = 10.62472 #inches
    font_size = 10
    
    sink_pickle = "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()

    #read pickles and make frames!
    pickle_files = sorted(glob.glob(save_dir + "movie_frame_*.pkl"))
    #Get start and end time
    file = open(pickle_files[0], 'rb')
    write_dict = pickle.load(file)
    file.close()
    start_time = write_dict['time']
    start_ind = np.argmin(abs(particle_data['time'] - start_time))
    
    file = open(pickle_files[-1], 'rb')
    write_dict = pickle.load(file)
    file.close()
    end_time = write_dict['time']
    end_ind = np.argmin(abs(particle_data['time'] - end_time))
    
    time_arr = []
    acc_arr = []
    mean_dens_array = []
    
    xmin = np.nan
    xmax = np.nan
    ymin = np.nan
    ymax = np.nan
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
                plt.clf()
                plt.figure(figsize=(12,8))
                #fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(single_col_width, 2.4*single_col_width))#, sharey=True)
                #plt.subplots_adjust(wspace=0.0)
                #plt.subplots_adjust(hspace=0.0)
            
                file = open(plot_pickle, 'rb')
                write_dict = pickle.load(file)
                file.close()
                
                curr_it = np.argmin(abs(particle_data['time'] - write_dict['time']))
                
                time_arr.append(write_dict['time'])
                acc_arr.append(write_dict['mdot'])
                mean_dens_array.append(np.mean(write_dict['density']))
                
                if np.isnan(xmin):
                    xmin = np.min(write_dict['density'].value)
                elif np.min(write_dict['density'].value) < xmin:
                    xmin = np.min(write_dict['density'].value)
                
                if np.isnan(xmax):
                    xmax = np.max(write_dict['density'].value)
                elif np.max(write_dict['density'].value) > xmax:
                    xmax = np.max(write_dict['density'].value)
                    
                if np.isnan(ymin):
                    ymin = np.nanmin(write_dict['radial_momentum'].value)
                elif np.nanmin(write_dict['radial_momentum'].value) < ymin:
                    ymin = np.nanmin(write_dict['radial_momentum'].value)
                    
                if np.isnan(ymax):
                    ymax = np.nanmax(write_dict['radial_momentum'].value)
                elif np.nanmax(write_dict['radial_momentum'].value) > ymax:
                    ymax = np.nanmax(write_dict['radial_momentum'].value)
                    
                if np.isnan(lin_thresh):
                    lin_thresh = np.nanmin(np.abs(write_dict['radial_momentum'].value))
                elif np.nanmin(np.abs(write_dict['radial_momentum'].value)) < lin_thresh:
                    lin_thresh = np.nanmin(np.abs(write_dict['radial_momentum'].value))


                plt.subplot(2,2,1)
                plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind])
                plt.scatter(particle_data['time'][curr_it], particle_data['separation'][curr_it])
                plt.xlim([particle_data['time'][start_ind], particle_data['time'][end_ind]])
                plt.ylim([np.min(particle_data['separation'][start_ind:end_ind]), np.max(particle_data['separation'][start_ind:end_ind])])
                plt.xlabel('Time (yr)')
                plt.ylabel('Separation (au)')
                
                plt.subplot(2,2,3)
                plt.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'][start_ind:end_ind])
                plt.scatter(particle_data['time'][curr_it], particle_data['mdot'][curr_it][0])
                plt.scatter(particle_data['time'][curr_it], particle_data['mdot'][curr_it][1])
                plt.xlim([particle_data['time'][start_ind], particle_data['time'][end_ind]])
                plt.ylim([np.min(particle_data['mdot'][start_ind:end_ind]), np.max(particle_data['mdot'][start_ind:end_ind])])
                plt.xlabel('Time (yr)')
                plt.ylabel('Separation (au)')
                
                plt.subplot(2,2,(2,4))
                plt.xscale('log')
                plt.yscale('symlog', linthresh=lin_thresh)
                plot = plt.scatter(write_dict['density'], write_dict['radial_momentum'])
                plt.xlim([xmin,xmax])
                plt.ylim([ymin,ymax])
                plt.xlabel('density (g/cm$^3$)')
                plt.ylabel('radial momentum (cm$\,$g/s)')
                
                '''
                axs.set_xscale('log')
                axs.set_yscale('symlog', linthresh=lin_thresh)
                plot = axs.scatter(write_dict['density'], write_dict['radial_momentum'])
                cbar = plt.colorbar(plot, pad=0.0)
                cbar.set_label(r"v$_{radial}$/v$_{magnitude}$", rotation=270, labelpad=14)
                axs.set_xlim([xmin,xmax])
                axs.set_ylim([ymin,ymax])
                axs.set_xlabel('density (g/cm$^3$)')
                axs.set_ylabel('radial momentum (cm$\,$g/s)')
                '''

                plt.savefig(file_name, bbox_inches='tight', dpi=300)
                print("Plotted", file_name, "for pickle", fit, "of", len(pickle_files))
        fit = fit + 1
