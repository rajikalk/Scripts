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
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="z")
    parser.add_argument("-make_pickles", "--make_movie_pickles", type=str, default='True')
    parser.add_argument("-make_frames", "--make_movie_frames", type=str, default='True')
    parser.add_argument("-width", "--plot_width", type=float, default=2000)
    parser.add_argument("-thickness", "--proj_thickness", type=float, default=None)
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
            proj_root_rank = int(rank/5)*5
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
                has_particles = True
                part_mass = dd['particle_mass'].in_units('msun')
                part_pos_fields = [field for field in ds.field_list if ('particle_pos' in field[1])&(field[0]=='all')&(field[1]!='particle_pos'+args.axis)]
                part_pos_x = dd[part_pos_fields[0]].in_units('au')
                part_pos_y = dd[part_pos_fields[1]].in_units('au')
                positions = yt.YTArray([part_pos_x,part_pos_y])
                part_vel_fields = [field for field in ds.field_list if ('particle_vel' in field[1])&(field[0]=='all')&(field[1]!='particle_vel'+args.axis)]
                part_vel_x = dd[part_vel_fields[0]].in_units('cm/s')
                part_vel_y = dd[part_vel_fields[1]].in_units('cm/s')
                velocities = yt.YTArray([part_vel_x,part_vel_y])
                part_info = {'particle_mass':part_mass,
                         'particle_position':positions,
                         'particle_velocities':velocities,
                         'accretion_rad':2.5*np.min(dd['dx'].in_units('au')),
                         'particle_tag':dd['particle_tag'],
                         'particle_form_time':dd['particle_creation_time']}
                primary_ind = np.argmin(dd['all', 'particle_creation_time'])
                center_pos = yt.YTArray([dd['all', 'particle_posx'][primary_ind].in_units('cm').value, dd['all', 'particle_posy'][primary_ind].in_units('cm').value, dd['all', 'particle_posz'][primary_ind].in_units('cm').value], 'cm')
            else:
                center_pos = yt.YTArray([0, 0, 0], 'cm')
                has_particles = False
                part_info = {}
                
            del dd
            
            #make list of projection fields: density, velocity, magnetic field
            if args.field == 'dens':
                proj_field_list = [('flash', 'dens')]
            else:
                try:
                    proj_field_list = [field for field in ds.derived_field_list if (args.field == field[1])]
                except:
                    proj_field_list = [field for field in ds.derived_field_list if (args.field in field[1])]
                    if len(proj_field_list) > 1:
                        proj_field_list = [proj_field_list[0]]
                
            proj_field_list = [field for field in ds.field_list if ('vel'in field[1])&(field[0]=='flash')&('vel'+args.axis not in field[1])] + [field for field in ds.field_list if ('mag'in field[1])&(field[0]=='flash')&('mag'+args.axis not in field[1])] + proj_field_list
        
            #define projection region
            plot_width = yt.YTQuantity(args.plot_width, 'au')
            if args.axis == 'z':
                if args.proj_thickness == None:
                    thickness = yt.YTQuantity(args.plot_width/2, 'au')
                else:
                    thickness = yt.YTQuantity(args.proj_thickness, 'au')
                
                left_corner = yt.YTArray([center_pos[0].in_units('au')-((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au')), center_pos[1].in_units('au')-((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au')), center_pos[2].in_units('au')-(thickness.in_units('au')/2)], 'AU')
                right_corner = yt.YTArray([center_pos[0].in_units('au')+((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au')), center_pos[1].in_units('au')+((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au')), center_pos[2].in_units('au')+(thickness.in_units('au')/2)], 'AU')
                region = ds.box(left_corner, right_corner)
            else:
                left_corner = yt.YTArray([center_pos[0].in_units('au')-((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au')), center_pos[1].in_units('au')-(plot_width.in_units('au')/2), center_pos[2].in_units('au')-((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au'))], 'AU')
                right_corner = yt.YTArray([center_pos[0].in_units('au')+((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au')), center_pos[1].in_units('au')+(plot_width.in_units('au')/2), center_pos[2].in_units('au')+((plot_width.in_units('au')/2)+yt.YTQuantity(100, 'au'))], 'AU')
                region = ds.box(left_corner, right_corner)
                
            test_fields = region['x'], region['y'], region['z'], region['velx'], region['vely'], region['velz'], region['mass']
            del test_fields
            test_fields = region['nearest_particle_index']
            del test_fields
            test_fields = region['L_gas_wrt_nearest_sink']
            del test_fields
            test_fields = region['Radial_velocity_wrt_primary']
            
            #Make projections of each field
            #proj_depth = yt.ProjectionPlot(ds, args.axis, [('flash', 'z'), ('gas', 'Neg_z'), ('flash', 'dz'), ('gas', 'Neg_dz')], width=(args.plot_width,'au'), weight_field=None, data_source=region, method='mip', center=(center_pos, 'AU'))
            #thickness = ((proj_depth.frb.data[('gas', 'Neg_z')].in_units('cm') + proj_depth.frb.data[('gas', 'Neg_dz')].in_units('cm')/2.) + (proj_depth.frb.data[('flash', 'z')].in_units('cm') + proj_depth.frb.data[('flash', 'dz')].in_units('cm')/2.))
            

            thickness_proj = yt.ProjectionPlot(ds, args.axis, ('gas', 'N_cells'), method='integrate', data_source=region, width=plot_width, weight_field=None, center=center_pos)
            thickness_arr = (thickness_proj.frb.data[('gas', 'N_cells')].in_cgs())
            fix_thickness = np.ones(np.shape(thickness_arr))*((plot_width/np.shape(thickness_arr)[0])**2*thickness).in_units('cm**3')
            
            proj_dict = {}
            for sto, field in yt.parallel_objects(proj_field_list, storage=proj_dict):
                #print("Projecting field", field, "on rank", rank)
                proj = yt.ProjectionPlot(ds, args.axis, field, method='integrate', data_source=region, width=plot_width, weight_field=args.weight_field, center=center_pos)
                if args.weight_field == None:
                    #thickness = (proj.bounds[1] - proj.bounds[0]).in_cgs() #MIGHT HAVE TO UPDATE THIS LATER
                    if field[1] == 'L_gas_wrt_primary_density':
                        proj_array = proj.frb.data[field].in_cgs()*fix_thickness/thickness_arr
                    else:
                        proj_array = proj.frb.data[field].in_cgs()/thickness_arr
                else:
                    proj_array = proj.frb.data[field].in_cgs()
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
                velx, vely, velz = mym.get_quiver_arrays(0, 0, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=args.quiver_arrows)
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
                file.close()
                print("created pickle for frame", file_int, "on rank", rank)
            elif size == 1:
                velx, vely, velz = mym.get_quiver_arrays(0, 0, X_image, proj_dict[list(proj_dict.keys())[1]], proj_dict[list(proj_dict.keys())[2]], no_of_quivers=args.quiver_arrows)
                file = open(pickle_file, 'wb')
                pickle.dump((X_image, Y_image, proj_dict[list(proj_dict.keys())[0]], proj_dict[list(proj_dict.keys())[3]], proj_dict[list(proj_dict.keys())[4]], X_image_vel, Y_image_vel, velx, vely, part_info, time_val), file)
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
    
    del usable_files
    del frames

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
