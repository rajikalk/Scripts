#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
#import pickle5 as pickle
import pickle
import os
import my_flash_module as mym

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="z")
    parser.add_argument("-make_pickles", "--make_movie_pickles", type=str, default='True')
    parser.add_argument("-make_frames", "--make_movie_frames", type=str, default='True')
    parser.add_argument("-width", "--plot_width", type=float, default=2000)
    
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-pf", "--presink_frames", help="How many frames do you want before the formation of particles?", type=int, default = 25)
    parser.add_argument("-end", "--end_time", help="What time do you want to the movie to finish at?", default=None, type=int)
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default = 0, type=int)
    
    parser.add_argument("-all_files", "--use_all_files", help="Do you want to make frames using all available files instead of at particular time steps?", type=str, default='False')
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#-------------------------------------------------------
#get input and output directory and arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
args = parse_inputs()
center_pos = [0, 0, 0]

if args.make_movie_pickles == 'True':
    
    files = sorted(glob.glob(input_dir + '*plt_cnt*'))
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=args.presink_frames, end_time=args.end_time)
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
    annotate_space = (x_image_min - x_image_max)/32.
    x_ind = []
    y_ind = []
    counter = 0
    while counter < 32:
        val = annotate_space*counter + annotate_space/2. + x_image_min
        x_ind.append(int(val))
        y_ind.append(int(val))
        counter = counter + 1
    X_image_vel, Y_image_vel = np.meshgrid(x_ind, y_ind)

    #Now let's iterate over the files and get the images we want to plot
    file_int = -1
    for fn in yt.parallel_objects(usable_files, njobs=int(size/5)):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
            if file_int > 0
                if usable_files[file_int] == usable_files[file_int-1]:
                    os.system('cp '+ output_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + output_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        file_counter = usable_files.index(fn)
        pickle_file = output_dir+"movie_frame_"+("%06d" % file_counter)+".pkl"
        make_pickle = True
        #if os.path.isfile(pickle_file) == True:
        if len(glob.glob(pickle_file)) == 1:
            make_pickle = False
        if make_pickle:
            #print(fn, "is going to rank", rank)
            proj_root_rank = int(rank/5)*5
            part_file = 'part'.join(fn.split('plt_cnt'))
            ds = yt.load(fn, particle_filename=part_file)
            time_val = ds.current_time.in_units('yr')
            
                        #Get particle data:
            dd = ds.all_data()
            if len([field for field in ds.field_list if 'particle_mass' in field[1]]) > 0:
                has_particles = True
                part_mass = dd['particle_mass'].in_units('msun')
                part_pos_fields = [field for field in ds.field_list if ('particle_pos' in field[1])&(field[0]=='all')&(field[1]!='particle_pos'+args.axis)]
                part_pos_x = dd[part_pos_fields[0]].in_units('au')
                part_pos_y = dd[part_pos_fields[1]].in_units('au')
                positions = np.array([part_pos_x,part_pos_y])
                part_info = {'particle_mass':part_mass,
                         'particle_position':positions,
                         'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
                         'particle_tag':dd['particle_tag']}
            else:
                has_particles = False
                part_info = {}
            
            #make list of projection fields: density, velocity, magnetic field
            proj_field_list = [('flash', 'dens')] + \
                [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vel'+args.axis not in field[1])] + \
                [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('mag'+args.axis not in field[1])]
        
            #define projection region
            if args.axis == 'z':
                left_corner = yt.YTArray([center_pos[0]-((args.plot_width/2)+100), center_pos[1]-((args.plot_width/2)+100), center_pos[2]-(0.5*(args.plot_width/2))], 'AU')
                right_corner = yt.YTArray([center_pos[0]+((args.plot_width/2)+100), center_pos[1]+((args.plot_width/2)+100), center_pos[2]+(0.5*(args.plot_width/2))], 'AU')
                region = ds.box(left_corner, right_corner)
            else:
                left_corner = yt.YTArray([center_pos[0]-((args.plot_width/2)+100), center_pos[1]-(0.5*(args.plot_width/2)), center_pos[2]-((args.plot_width/2)+100)], 'AU')
                right_corner = yt.YTArray([center_pos[0]+((args.plot_width/2)+100), center_pos[1]+(0.5*(args.plot_width/2)), center_pos[2]+((args.plot_width/2)+100)], 'AU')
                region = ds.box(left_corner, right_corner)
            
            #Make projections of each field
            proj_dict = {}
            for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
                #print("Projecting field", field, "on rank", rank)
                proj = yt.ProjectionPlot(ds, args.axis, field, method='integrate', data_source=region, width=(args.plot_width,'au'), weight_field=None, center=(center_pos, 'AU'))
                thickness = (proj.bounds[1] - proj.bounds[0]).in_cgs() #MIGHT HAVE TO UPDATE THIS LATER
                proj_array = proj.frb.data[field].in_cgs()/thickness
                if args.axis == 'y':
                    proj_array = proj_array.T
                #print(field, "projection =", proj_array)
                sto.result_id = field[1]
                sto.result = proj_array
                #if rank == proj_root_rank:
                #    proj_dict[field[1]] = proj_array
                #else:
                #    file = open(pickle_file.split('.pkl')[0] + '_proj_data_' + str(proj_root_rank)+ str(proj_field_list.index(field)) + '.pkl', 'wb')
                #    pickle.dump((field[1], proj_array), file)
                #    file.close()
            
            if rank == proj_root_rank and size > 1:
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], part_info, time_val), file)
                file.close()
                print("created pickle for frame", file_counter, "on rank", rank)
            elif size == 1:
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], part_info, time_val), file)
                file.close()
                print("created pickle for frame", file_counter)

    print("finished making movie frame pickles on rank", rank)

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
    pickle_files = sorted(glob.glob(output_dir+"movie_frame_*.pkl"))
    no_frames = len(pickle_files)

    rit = -1
    for pickle_file in pickle_files:
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            frame_no = int(pickle_file.split('_')[-1].split('.')[0])
            file_name = output_dir + "movie_frame_" + ("%06d" % frame_no)
            if os.path.isfile(file_name+'.jpg') == False:
                print('making frame from', pickle_file, 'on rank', rank)
                file = open(pickle_file, 'rb')
                X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
                file.close()

                plt.clf()
                fig, ax = plt.subplots()
                ax.set_xlabel('AU', labelpad=-1, fontsize=10)
                ax.set_ylabel('AU', fontsize=10) #, labelpad=-20
                xlim = [np.min(X_image).value, np.max(X_image).value]
                ylim = [np.min(Y_image).value, np.max(Y_image).value]
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                
                cmap=plt.cm.gist_heat
                plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(), rasterized=True, zorder=1)
                plt.gca().set_aspect('equal')

                if frame_no > 0 or time_val > -1.0:
                    plt.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
                else:
                    plt.streamplot(X_image, Y_image, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
                cbar = plt.colorbar(plot, pad=0.0)
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=5)

                if len(part_info.keys())>0:
                    mym.set_global_font_size(8)
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7)

                cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
                plt.tick_params(axis='both', which='major')# labelsize=16)
                for line in ax.xaxis.get_ticklines():
                    line.set_color('white')
                for line in ax.yaxis.get_ticklines():
                    line.set_color('white')
                    
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
