#!/usr/bin/env python
import sys
import yt
yt.enable_parallelism()
import glob
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
                slice = yt.OffAxisSlicePlot(ds, proj_vector_unit, field, width=plot_width, center=center_pos, north_vector=north_unit)
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
