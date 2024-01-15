#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import pickle5 as pickle
import os
import my_flash_module as mym
import my_flash_fields as myf

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-make_pickles", "--make_movie_pickles", type=str, default='True')
    parser.add_argument("-make_frames", "--make_movie_frames", type=str, default='True')
    parser.add_argument("-width", "--plot_width", type=float, default=2000)
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="dens", type=str)
    #parser.add_argument("-cbar_lim", "-cbar_limits", type=str, default=[])
    
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", type=int, default = 25)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-start", "--start_time", help="What time do you want to the movie to finish at?", default=0, type=int)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    parser.add_argument("-no_quiv", "--quiver_arrows", default=31., type=float)
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=None)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=None)
    parser.add_argument("-weight", "--weight_field", help="set weight field?", type=str, default=None)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=None)
    parser.add_argument("-all_files", "--use_all_files", help="Do you want to make frames using all available files instead of at particular time steps?", type=str, default='False')
    parser.add_argument("-update_vel", "--update_velocity_field", help="update velocity field to be wrt to primary", type=str, default='False')
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def projected_vector(vector, proj_vector):
    """
    Calculates the position of vector projected onto proj_vector
    """
    vector_units = vector.units
    proj_v_x = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[0]
    proj_v_y = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[1]
    proj_v_z = (np.dot(vector, proj_vector)/np.dot(proj_vector,proj_vector))*proj_vector[2]
    proj_v = yt.YTArray(np.array([proj_v_x,proj_v_y,proj_v_z]).T, vector_units)
    return proj_v

#-------------------------------------------------------
#get input and save_dir directory and arguments
input_dir = sys.argv[1]
save_dir = sys.argv[2]
args = parse_inputs()
font_size = 10
mym.set_global_font_size(font_size)

if args.make_movie_pickles == 'True':
    
    files = sorted(glob.glob(input_dir + '*plt_cnt*'))
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=args.presink_frames, end_time=args.end_time, start_time=args.start_time)
    print('generated frame times')
    sys.stdout.flush()
    CW.Barrier()
    
    if args.use_all_files == 'False':
        no_frames = len(m_times)
        m_times = m_times[args.start_frame:]
        usable_files = mym.find_files(m_times, files)
        frames = list(range(args.start_frame, no_frames))
    elif args.use_all_files != 'False' and args.plot_time != None:
        usable_files = mym.find_files([args.plot_time], files)
        start_index = files.index(usable_files[0])
        args.plot_time = None
        end_file = mym.find_files([args.end_time], files)
        end_index = files.index(end_file[0])
        usable_files = files[start_index:end_index]
        frames = list(range(len(usable_files)))
        no_frames = len(usable_files)
    else:
        start_file = mym.find_files([0], files)
        usable_files = files[files.index(start_file[0]):]
        frames = list(range(len(usable_files)))
        no_frames = len(usable_files)
    print('found usable files for frames')

    #Get movie files
    #movie_files = sorted(glob.glob(input_dir + '*plt_cnt*'))
    #if rank == 1:
    #    print("Movie files=", movie_files)

    #Calculate image grid:
    fn = usable_files[-1]
    part_file = 'part'.join(fn.split('plt_cnt'))
    ds = yt.load(fn, particle_filename=part_file)
    x_image_min = yt.YTQuantity(-1*args.plot_width/2, 'au')
    x_image_max = yt.YTQuantity(args.plot_width/2, 'au')
    #x_image_min = -1*ds.domain_width.in_units('au')[0]/2
    #x_image_max = ds.domain_width.in_units('au')[0]/2
    x_range = np.linspace(x_image_min, x_image_max, 800)
    X_image, Y_image = np.meshgrid(x_range, x_range)
    annotate_space = (x_image_max - x_image_min)/args.quiver_arrows
    x_ind = []
    y_ind = []
    counter = 0
    while counter < args.quiver_arrows:
        val = annotate_space*counter + annotate_space/2. + x_image_min
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    X_image_vel, Y_image_vel = np.meshgrid(x_ind, y_ind)

    #Now let's iterate over the files and get the images we want to plot
    file_int = -1
    if size > 1:
        njobs = int(size/5)
    else:
        njobs = 1
    for fn in yt.parallel_objects(usable_files, njobs=njobs):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
            if usable_files[file_int] == usable_files[file_int-1]:
                os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        make_pickle = False
        if args.plot_time is None:
            pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl"
        else:
            pickle_file = save_dir + "time_" + str(args.plot_time) +".pkl"
        if os.path.isfile(pickle_file) == False:
            make_pickle = True
        elif os.path.isfile(pickle_file) == True:
            if os.stat(pickle_file).st_size == 0:
                make_pickle = True
        if usable_files[file_int] == usable_files[file_int-1]:
            os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        #if os.path.isfile(pickle_file) == True:
        #if len(glob.glob(pickle_file)) == 1:
        #    make_pickle = False
        if make_pickle:
            #print(fn, "is going to rank", rank)
            slice_root_rank = int(rank/5)*5
            part_file = 'part'.join(fn.split('plt_cnt'))
            ds = yt.load(fn, particle_filename=part_file)
            if args.plot_time == None:
                time_val = m_times[file_int]#ds.current_time.in_units('yr')
            else:
                dd = ds.all_data()
                time_val = int(yt.YTQuantity(ds.current_time.value - np.min(dd['particle_creation_time']).value, 's').in_units('yr').value)
            
                        #Get particle data:
            dd = ds.all_data()
            #Load fields
            #del test_fields
            if len([field for field in ds.field_list if 'particle_mass' in field[1]]) > 0:
                if len(dd['particle_mass']) > 1:
                    Primary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[0], dd['particle_posy'].in_units('au')[0], dd['particle_posz'].in_units('au')[0]])
                    Secondary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[1], dd['particle_posy'].in_units('au')[1], dd['particle_posz'].in_units('au')[1]])
                    d_pos = Secondary_pos - Primary_pos

                    Primary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[0], dd['particle_vely'].in_units('km/s')[0], dd['particle_velz'].in_units('km/s')[0]])
                    Secondary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[1], dd['particle_vely'].in_units('km/s')[1], dd['particle_velz'].in_units('km/s')[1]])
                    d_vel = Secondary_vel - Primary_vel
                    L_vec = np.cross(d_pos, d_vel).T
                    proj_vector_unit = L_vec/np.sqrt(np.sum(L_vec**2))

                    part_info = {'particle_mass':dd['particle_mass'][:2].in_units('msun'),
                                 'particle_position':yt.YTArray([Primary_pos, Secondary_pos]).T,
                                 'particle_velocities':yt.YTArray([Primary_vel, Secondary_vel]).T,
                                 'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
                                 'particle_tag':dd['particle_tag'][:2],
                                 'particle_form_time':dd['particle_creation_time'][:2]}
                    pos_array = yt.YTArray([Primary_pos, Secondary_pos])
                    east_unit_vector = [1, 0, 0]
                    north_unit = [0, 1, 0]
                    
                    center_pos = Primary_pos
                    center_vel = Primary_vel

                    projected_particle_posy = projected_vector(pos_array, north_unit)
                    slice_part_y_mag = np.sqrt(np.sum((projected_particle_posy**2), axis=1))
                    slice_part_y_unit = (projected_particle_posy.T/slice_part_y_mag).T
                    north_sign = np.dot(north_unit, slice_part_y_unit.T)
                    slice_part_y = slice_part_y_mag*north_sign
                    slice_part_y = np.nan_to_num(slice_part_y)

                    projected_particle_posx = projected_vector(pos_array, east_unit_vector)
                    slice_part_x_mag = np.sqrt(np.sum((projected_particle_posx**2), axis=1))
                    slice_part_x_unit = (projected_particle_posx.T/slice_part_x_mag).T
                    east_sign = np.dot(east_unit_vector, slice_part_x_unit.T)
                    slice_part_x = slice_part_x_mag*east_sign
                    slice_part_x = np.nan_to_num(slice_part_x)
                    
                    projected_particle_posz = projected_vector(pos_array, proj_vector_unit)
                    slice_part_z_mag = np.sqrt(np.sum((projected_particle_posz**2), axis=1))
                    slice_part_z_unit = (projected_particle_posz.T/slice_part_z_mag).T
                    slice_sign = np.dot(proj_vector_unit, slice_part_z_unit.T)
                    slice_part_z = slice_part_z_mag*slice_sign
                    slice_part_z = np.nan_to_num(slice_part_z)

                    part_info['particle_position'] = yt.YTArray([slice_part_x, slice_part_y])
                    
                    center_vel_proj_y = projected_vector(center_vel, north_unit)
                    center_vel_y = np.sqrt(center_vel_proj_y.T[0]**2 + center_vel_proj_y.T[1]**2 + center_vel_proj_y.T[2]**2).in_units('cm/s')
                    
                    center_vel_proj_x = projected_vector(center_vel, east_unit_vector)
                    center_vel_x = np.sqrt(center_vel_proj_x.T[0]**2 + center_vel_proj_x.T[1]**2 + center_vel_proj_x.T[2]**2).in_units('cm/s')
                    
                    center_vel_proj_rv = projected_vector(center_vel, proj_vector_unit)
                    center_vel_rv_mag = np.sqrt(np.sum(center_vel_proj_rv**2))
                    center_vel_rv_unit = center_vel_proj_rv/center_vel_rv_mag
                    rv_sign = np.dot(proj_vector_unit, center_vel_rv_unit)
                    center_vel_rv = center_vel_rv_mag*rv_sign
                    
                    center_vel_image = np.array([center_vel_x, center_vel_y])
                
                    #set vectors:
                    myf.set_normal(proj_vector_unit)
                    myf.set_east_vector(east_unit_vector)
                    myf.set_north_vector(north_unit)
                
                else:
                    Primary_pos = yt.YTArray([dd['particle_posx'].in_units('au')[0], dd['particle_posy'].in_units('au')[0], dd['particle_posz'].in_units('au')[0]])
                    Primary_vel = yt.YTArray([dd['particle_velx'].in_units('km/s')[0], dd['particle_vely'].in_units('km/s')[0], dd['particle_velz'].in_units('km/s')[0]])
                    
                    proj_vector_unit = [0, 0, 1]

                    part_info = {'particle_mass':dd['particle_mass'].in_units('msun'),
                                 'particle_position':yt.YTArray([Primary_pos]).T,
                                 'particle_velocities':yt.YTArray([Primary_vel]).T,
                                 'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
                                 'particle_tag':dd['particle_tag'],
                                 'particle_form_time':dd['particle_creation_time']}
                    pos_array = yt.YTArray([Primary_pos])
                    east_unit_vector = [1, 0, 0]
                    north_unit = [0, 1, 0]
                    
                    center_pos = Primary_pos
                    center_vel = Primary_vel
            else:
                import pdb
                pdb.set_trace()
                center_pos = yt.YTArray([0, 0, 0], 'cm')
                has_particles = False
                part_info = {}
                
            del dd
            
            myf.set_normal(proj_vector_unit)
            myf.set_east_vector(east_unit_vector)
            myf.set_north_vector(north_unit)
            
            #make list of sliceection fields: density, velocity, magnetic field
            if args.field == 'dens':
                slice_field_list = [('flash', 'dens')]
            else:
                try:
                    slice_field_list = [field for field in ds.derived_field_list if (args.field == field[1])]
                except:
                    slice_field_list = [field for field in ds.derived_field_list if (args.field in field[1])]
                    if len(slice_field_list) > 1:
                        slice_field_list = [slice_field_list[0]]
                
            slice_field_list = slice_field_list + [('gas', 'Proj_x_velocity'), ('gas', 'Proj_y_velocity')]
        
            #define sliceection region
            plot_width = yt.YTQuantity(args.plot_width, 'au')
            
            slice_dict = {}
            for sto, field in yt.parallel_objects(slice_field_list, storage=slice_dict):
                slice = yt.OffAxisSlicePlot(ds, proj_vector_unit, field, width=(plot_width, 'au'), center=Primary_pos, north_vector=[0, 1, 0])
                if args.weight_field == None:
                    #thickness = (slice.bounds[1] - slice.bounds[0]).in_cgs() #MIGHT HAVE TO UPDATE THIS LATER
                    slice_array = slice.frb.data[field].in_cgs()
                else:
                    slice_array = slice.frb.data[field].in_cgs()
                if args.axis == 'y':
                    slice_array = slice_array.T
                #print(field, "sliceection =", slice_array)
                sto.result_id = field[1]
                sto.result = slice_array
            
            if rank == slice_root_rank and size > 1:
                velx, vely, velz = mym.get_quiver_arrays(0, 0, X_image, slice_dict[list(slice_dict.keys())[1]], slice_dict[list(slice_dict.keys())[2]], no_of_quivers=args.quiver_arrows)
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, slice_dict[list(slice_dict.keys())[0]], slice_dict[list(slice_dict.keys())[3]], slice_dict[list(slice_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
                file.close()
                print("created pickle for frame", file_int, "on rank", rank)
            elif size == 1:
                velx, vely, velz = mym.get_quiver_arrays(0, 0, X_image, slice_dict[list(slice_dict.keys())[1]], slice_dict[list(slice_dict.keys())[2]], no_of_quivers=args.quiver_arrows)
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, slice_dict[list(slice_dict.keys())[0]], slice_dict[list(slice_dict.keys())[3]], slice_dict[list(slice_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
                file.close()
                print("created pickle for frame", file_int, "of", len(m_times))

    print("finished making movie frame pickles on rank", rank)
    
    sys.stdout.flush()
    CW.Barrier()
    
    if rank == 0:
        for sim_file in usable_files:
            file_int = file_int + 1
            try:
                if usable_files[file_int] == usable_files[file_int-1]:
                    os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
            except:
                continue
        print('Finished copying pickles that use the same file for the same frame')

sys.stdout.flush()
CW.Barrier()

if args.make_movie_frames == 'True':
    #if size == 1:
    #Make frames.
    import matplotlib as mpl
    #mpl.rcParams['pdf.fonttype'] = 42
    #mpl.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    #plt.rcParams['figure.dpi'] = 300
    from matplotlib.colors import LogNorm
    import matplotlib.patheffects as path_effects
    import my_flash_module as mym

    #Let's get the pickle files
    if args.plot_time != None:
        pickle_files = [save_dir + "time_" + str(args.plot_time) +".pkl"]
    else:
        pickle_files = sorted(glob.glob(save_dir+"movie_frame_*.pkl"))
    no_frames = len(pickle_files)

    rit = -1
    for pickle_file in pickle_files:
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            frame_no = int(pickle_file.split('_')[-1].split('.')[0])
            if args.plot_time != None:
                file_name = save_dir + "plot_time_" + ("%06d" % frame_no)
            else:
                file_name = save_dir + "movie_frame_" + ("%06d" % frame_no)
            if os.path.isfile(file_name+'.jpg') == False:
                print('making frame from', pickle_file, 'on rank', rank)
                file = open(pickle_file, 'rb')
                X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
                file.close()
                
                if args.update_velocity_field != 'False':
                    primary_ind = np.argmin(part_info['particle_form_time'])
                    part_velx = part_info['particle_velocities'][0][primary_ind]
                    part_vely = part_info['particle_velocities'][1][primary_ind]
                    velx = velx - part_velx.value
                    vely = vely - part_vely.value

                plt.clf()
                fig, ax = plt.subplots()
                ax.set_xlabel('AU', labelpad=-1, fontsize=10)
                ax.set_ylabel('AU', fontsize=10) #, labelpad=-20
                xlim = [np.min(X_image).value, np.max(X_image).value]
                ylim = [np.min(Y_image).value, np.max(Y_image).value]
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                
                
                if args.colourbar_min == None:
                    cmin = np.min(image)
                else:
                    cmin = args.colourbar_min
                    
                if args.colourbar_max == None:
                    cmax = np.max(image)
                else:
                    cmax = args.colourbar_max
                
                '''
                if args.colourbar_min == None:
                    if args.field == 'dens':
                        cmin = 1.e-16
                    else:
                        if 'spec' in args.field:
                            cmin = 1.e18
                        else:
                            cmin = 1.e45
                else:
                    cmin = args.colourbar_min
                
                if args.colourbar_max == None:
                    if args.field == 'dens':
                        cmax = 5.e-14
                    else:
                        if 'spec' in args.field:
                            cmax = 1.e20
                        else:
                            cmax = 1.e48
                else:
                    cmax = args.colourbar_max
                '''
                
                cbar_lims = [cmin, cmax]
                
                if args.standard_vel == None:
                    if args.field == 'dens':
                        if args.axis == 'z':
                            stdvel = 2
                        else:
                            stdvel = 5
                    else:
                        stdvel = 5
                else:
                    stdvel = args.standard_vel
                
                cmap=plt.cm.gist_heat
                if 'Relative_keplerian_velocity_wrt_primary' in args.field or cmin < 0:
                    plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.YlGn, vmin=cbar_lims[0], vmax=cbar_lims[1], rasterized=True, zorder=1)
                else:
                    if np.isnan(cbar_lims[0]):
                        plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(), rasterized=True, zorder=1)
                    else:
                        plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
                plt.gca().set_aspect('equal')

                if 'L_gas_wrt_primary_spec' in args.field:
                    #make contour plot
                    exp_min = np.log10(cbar_lims[0])
                    exp_max = np.log10(cbar_lims[1])
                    n_level = 11#(exp_max-exp_min)*2 + 1
                    contour_levels = np.logspace(exp_min, exp_max, int(n_level))
                    CS = ax.contour(X_image,Y_image,image, locator=plt.LogLocator(), linewidths=0.5, levels=contour_levels, colors='k')
                else:
                    if frame_no > 0 or time_val > -1.0:
                        plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
                    else:
                        plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, minlength=0.5, zorder=2)
                        #plt.streamplot(X_image, Y_image, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
                cbar = plt.colorbar(plot, pad=0.0)
                try:
                    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
                except:
                    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)

                if len(part_info.keys())>0:
                    if 'particle_form_time' in part_info.keys():
                        if len(part_info['particle_form_time']) > 1:
                            if np.min(part_info['particle_form_time'][1:] - part_info['particle_form_time'][:-1]) < 0:
                                sort_inds = np.argsort(part_info['particle_form_time'])
                                part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                                part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                                part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                                part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
                    else:
                        print("pickle doesn't have sink formation time")
                        os.remove(pickle_file)
                            
                    primary_ind = np.argmax(part_info['particle_mass'])
                    part_info['particle_position'][0] = part_info['particle_position'][0] - part_info['particle_position'][0][primary_ind]
                    part_info['particle_position'][1] = part_info['particle_position'][1] - part_info['particle_position'][1][primary_ind]                    #Update particle position
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=7)

                plt.tick_params(axis='both', which='major')# labelsize=16)
                for line in ax.xaxis.get_ticklines():
                    line.set_color('white')
                for line in ax.yaxis.get_ticklines():
                    line.set_color('white')
                    
                cbar.set_label(args.field + " (" + str(image.units)+")", rotation=270, labelpad=14, size=10)
                '''
                if args.field == 'dens':
                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
                elif 'Relative_keplerian_velocity_wrt_primary' in args.field:
                    cbar.set_label(r"Relative Keplerian Velocity", rotation=270, labelpad=14, size=10)
                elif 'spec' in args.field:
                    cbar.set_label(r"Specific angular momentum (g$\,$cm$^{2}/s$)", rotation=270, labelpad=14, size=10)
                else:
                    cbar.set_label(r"Angular momentum (g$\,$cm$^{2}/s$)", rotation=270, labelpad=14, size=10)
                '''
                time_string = "$t$="+str(int(time_val))+"yr"
                time_string_raw = r"{}".format(time_string)
                time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
                time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                if size > 1:
                    try:
                        plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
                        #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                        print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                    except:
                        print("couldn't save for the dviread.py problem. Make frame " + str(frame_no) + " on ipython")
                else:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
                    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    print('Created frame', (frame_no), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')

    print("Finished plotting frames on rank", rank)
