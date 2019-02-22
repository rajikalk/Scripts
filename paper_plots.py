import yt
import my_module as mym
import matplotlib.pyplot as plt
import my_fields as myf
from subprocess import call
import sys
import os
import glob
import numpy as np
import h5py
import pickle
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import csv
import matplotlib.patheffects as path_effects
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sp", "--slice_plot", help="Did you want to plot create a slice plot?", type=str, default="False")
    parser.add_argument("-pp", "--profile_plot", help="Did you want to plot a profile plot?", type=str, default="False")
    parser.add_argument("-col_mean", "--calc_column_mean", help="Did you want calculate the column mean? eg. for column density", type=str, default="False")
    parser.add_argument("-fc", "--force_comp", help="Did you want to create a plot comparing pressure?", type=str, default="False")
    parser.add_argument("-fp", "--force_on_particles", help="Did you want to create a plot of the force on the particles?", type=str, default="False")
    parser.add_argument("-bp", "--b_mag", help="Did you want to create a plot of where the magnetic fiedl is 30 degrees", type=str, default="False")
    parser.add_argument("-op", "--outflow_pickle", help="Do you want to measure the outflows?", type=str, default="False")
    parser.add_argument("-po", "--plot_outflows", help="Do you want to plot the outflows now that you've measured them?", type=str, default="False")
    parser.add_argument("-sep", "--separation", help="Do you want to plot the separation of the particles?", type=str, default="False")
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default=4.1, type=float)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int, default=0)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-c", "--center", help="What center do you want to set for everything?, if 3 it combines all centers", type=int, default=0)
    parser.add_argument("-ic", "--image_center", help="Where would you like to center the image?", type=int, default=0)
    parser.add_argument("-ppm", "--profile_plot_multi", help="Did you want to plot a profile plot with multiple lines?", type=str, default="False")
    
    #movie plot args
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="True")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-16)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    
    #slice plot args
    parser.add_argument("-f", "--field", help="What is the field you want to have a sliceplot of?", type=str, default="Relative_Keplerian_Velocity")
    parser.add_argument("-res", "--resolution", help="What resolution do you want the slices to have?", type=int, default=1024)
    parser.add_argument("-af", "--annotate_field", help="What field of the particles do you want to annotate", type=str, default="particle_mass")
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=None)
    
    #profile plot args
    parser.add_argument("-rmax", "--r_max", help="radius of measuing volume", type=float, default=500.)
    parser.add_argument("-dt", "--disk_thickness", help="How far above and below the midplane do you want your profile to go?", type=float, default=100.)
    parser.add_argument("-xf", "--x_field", help="x axis of the profile plot?", type=str, default="Distance_from_Center")
    parser.add_argument("-wf", "--weight_field", help="any weight field?", type=str, default=None)
    parser.add_argument("-log", "--logscale", help="Want to use a log scale?", type=str, default="False")
    parser.add_argument("-pb", "--profile_bins", help="how many bins do you want for the profile?", type=int, default=None)
    parser.add_argument("-zb", "--z_bins", help="how many z bins do you want when sampling points?", type=int, default=5.)
    parser.add_argument("-nsp", "--no_sampled_points", help="how many random points do you want to randomly sample?", type=int, default=2500)
    parser.add_argument("-xu", "--x_units", help="x units for profile plot", type=str, default="AU")
    parser.add_argument("-yu", "--y_units", help="y units for profile plot", type=str, default=None)
    parser.add_argument("-yl", "--y_label", help="what is the y_label you will use", type=str, default=None)
    parser.add_argument("-sn", "--save_name", help="If not defined, it's the same as the field", type=str, default=None)
    parser.add_argument("-cumul", "--cumulative", help="do you want to calculate the cumulative bin values?", type=str, default='False')
    
    #pressure plot args
    parser.add_argument("-ts", "--time_step", help="would you like to plot multiple times?", type=int, default=500)
    
    parser.add_argument("-pd", "--pickle_dump", help="do you want to pickle data?", type=str, default="False")
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=12)
    parser.add_argument("-cu", "--c_units", help="What units do you want the colorbar in")
    parser.add_argument("-et", "--endtime", default=2001, type=int)
    parser.add_argument("-stdv", "--standard_vel", default=5, type=float)
    
    #yt_slice
    parser.add_argument("-ys", "--yt_slice", help="did you want to make a yt-slice?", type=str, default="False")
    parser.add_argument("-st", "--slice_thickness", help="How thick do you want the slice to be?", default=300.0, type=float)
    
    parser.add_argument("-mov", "--produce_movie", help="Do you want to make a series of plots for a movie?", type=str, default="False")
    parser.add_argument("-sf", "--start_frame", help="if you don't want to start at frame 0, what frame do you want to start at?", type=int, default=0)
    parser.add_argument("-ef", "--end_frame", help="do you only want to plot up to t a certain number of frames?", default=None, type=int)
    
    parser.add_argument("-totm", "--total_mass_of_system", help="do you want to plot the total mass of your system?", type=str, default="False")
    parser.add_argument("-ni", "--non_ideal", help="do you want to plot the non-ideal regimes we should think about", type=str, default="False")
    parser.add_argument("-proj_or", "--projection_orientation", help="Do you want to set the projection orientation? give as angle (in degrees) from positive y-axis", default=None, type=float)
    parser.add_argument("-proj_ax", "--projection_axis", help="defaults to 'xz' and takes into account any projection orientation. Or can be set to xy", default="xz", type=str)
    parser.add_argument("-amb", "--angular_momentum_budget", help="do you want plot the angular momenutum budget", type=str, default="False")
    
    parser.add_argument("-read_part_file", "--read_particle_file", help="dow you want to read in the particle file adn save as a pickle?", type=str, default="False")
    parser.add_argument("-phasefolded", "--phasefolded_accretion", help="do you want to plot the phasefolded accretion", type=str, default="False")
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#======================================================================================================

args = parse_inputs()
mym.set_global_font_size(args.text_font)

path = sys.argv[1]
save_dir = sys.argv[2]
if save_dir[-1] != '/':
    save_dir = save_dir + '/'
if os.path.exists(save_dir) == "False":
    os.makedirs(save_dir)

myf.set_center(args.center)

movie_files = sorted(glob.glob(path + 'WIND_slice_*'))
sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt_*'))
if len(movie_files) > 0:
    sink_form_time = mym.find_sink_formation_time(movie_files)

    if args.produce_movie != "False":
        m_times = mym.generate_frame_times(sim_files, args.time_step, presink_frames=0)
    else:
        m_times = [args.plot_time]
    no_of_frames = len(m_times)
    if args.end_frame != None:
        final_frame = args.end_frame + 1
        m_times = m_times[args.start_frame:final_frame]
    else:
        m_times = m_times[args.start_frame:]

    if args.slice_plot != "False":
        usable_movie_files = mym.find_files(m_times, movie_files)
    usable_files = mym.find_files(m_times, sim_files)
    del movie_files
    del sim_files

    #Now iterate over the files
    if args.end_frame == None:
        final_frame = len(usable_files)
    else:
        final_frame = args.end_frame + 1
    for fit in range(len(usable_files)):
        if fit == 0 or usable_files[fit] != usable_files[fit-1]:
            print "CREATING FRAME", args.start_frame + fit, "WITH FILE", usable_files[fit]
            if args.slice_plot != "False":
                movie_file = usable_movie_files[fit]
            file = usable_files[fit]
            part_file = file[:-12] + 'part' + file[-5:]
            ds = yt.load(file, particle_filename=part_file)
            dd = ds.all_data()
        else:
            print "CREATING FRAME", args.start_frame + fit, "WITH FILE", usable_files[fit-1]

        if args.non_ideal == 'True':
            inds = np.where((dd['temperature']<281)&(dd['temperature']>279))
            times_str = str(m_times[fit])
            plt.clf()
            plt.plot(np.log10(dd['H_nuclei_density'][inds].value), np.log10(dd['magnetic_field_strength'][inds].value), 'o')
            plt.plot([12.2, 14.2], [-0.5, 1.5], 'k-')
            plt.plot([13.8, 15], [-0.5, 0.7], 'k-')
            #plt.xlim([12, 15])
            #plt.ylim([-0.5, 1.5])
            plt.xlabel('log n$_\mathrm{H}$ (cm$^3$)')
            plt.ylabel('log B (G)')
            plt.title('At time '+ times_str + 'yr')
            save_name = save_dir + 'konigl_salmeron_time_'+times_str+'.eps'
            plt.savefig(save_name)
            print "created file:", save_name

        if args.slice_plot == 'True':
            print "MAKING SLICE PLOT WIHT FILE:", file
            n_bins = np.ceil(np.sqrt((args.resolution/2.)**2. + (args.resolution/2.)**2.))
            #myf.set_n_bins(n_bins)
            save_image_name = save_dir + "Slice_Plot_time_" + str(args.plot_time) + ".eps"
            X, Y, X_vel, Y_vel, cl = mym.initialise_grid(movie_file, zoom_times=args.zoom_times)
            part_info = mym.get_particle_data(file, axis='xy')
            center_vel = [0.0, 0.0, 0.0]
            if args.image_center != 0:
                center_vel = [dd['particle_velx'][args.image_center-1].value, dd['particle_vely'][args.image_center-1].value, dd['particle_velz'][args.image_center-1].value]
                x_pos = np.round(part_info['particle_position'][0][args.image_center-1]/cl)*cl
                y_pos = np.round(part_info['particle_position'][1][args.image_center-1]/cl)*cl
                X = X + x_pos
                Y = Y + y_pos
                X_vel = X_vel + x_pos
                Y_vel = Y_vel + y_pos
            #myf.set_coordinate_system('cylindrical')
            myf.set_coordinate_system('spherical')
            dummy1, dummy2, dummy3, magx_grid, dummy4 = mym.sliceplot(ds, X, Y, 'magx', resolution=args.resolution, center=args.center, units=args.c_units)
            dummy1, dummy2, dummy3, magy_grid, dummy4 = mym.sliceplot(ds, X, Y, 'magy', resolution=args.resolution, center=args.center, units=args.c_units)
            del dummy1
            del dummy2
            del dummy3
            if args.field == 'Relative_Keplerian_Velocity':
                fig, ax, xy, field_grid, weight_field = mym.sliceplot(ds, X, Y, args.field, resolution=args.resolution, center=args.center, units=args.c_units, weight=args.weight_field)
            else:
                
                fig, ax, xy, field_grid, weight_field = mym.sliceplot(ds, X, Y, args.field, resolution=args.resolution, center=args.center, units=args.c_units, weight=args.weight_field, log=True, cbar_label='Density (gcm$^{-3}$)', cmap=plt.cm.gist_heat, vmin=args.colourbar_min, vmax=args.colourbar_max)
            f = h5py.File(movie_file, 'r')
            time_stamp = (ds.current_time.value - np.min(dd['particle_creation_time'].value))/yt.units.yr.in_units('s').value
            time_string = '$t=$'+str(int(np.round(time_stamp/1000.)*1000.)) + 'yr'
            dim = np.shape(f['dens_slice_xy'])[0]
            xmin_full = f['minmax_xyz'][0][0]/yt.units.AU.in_units('cm').value
            xmax = f['minmax_xyz'][0][1]/yt.units.AU.in_units('cm').value
            cl = (xmax-xmin_full)/dim
            x_pos_min = int(np.round(np.min(X) - xmin_full))/cl
            y_pos_min = int(np.round(np.min(Y) - xmin_full))/cl
            if f['velx_slice_xy'].shape == (2048,2048):
                velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X, f['velx_slice_xy'], f['vely_slice_xy'], center_vel=center_vel[:2])
            else:
                velx, vely = mym.get_quiver_arrays(y_pos_min, x_pos_min, X, f['velx_slice_xy'][:,:,0], f['vely_slice_xy'][:,:,0], center_vel=center_vel[:2])
            f.close()

            if args.ax_lim != None:
                if args.image_center == 0:
                    xlim = [-1*args.ax_lim, args.ax_lim]
                    ylim = [-1*args.ax_lim, args.ax_lim]
                else:
                    xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center-1], args.ax_lim + part_info['particle_position'][0][args.image_center-1]]
                    ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center-1], args.ax_lim + part_info['particle_position'][1][args.image_center-1]]
            else:
                xlim = [np.min(X), np.max(X)]
                ylim = [np.min(Y), np.max(Y)]
            limits = [xlim, ylim]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            if args.pickle_dump == "False":
                if fit == 0 and len(usable_files) != 1:
                    plt.streamplot(xy[0], xy[1], magx_grid.value, magy_grid.value, density=4, linewidth=0.25, minlength=0.5)
                else:
                    plt.streamplot(xy[0], xy[1], magx_grid.value, magy_grid.value, density=4, linewidth=0.25, minlength=0.5, arrowstyle='-')
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, standard_vel=args.standard_vel)
                if args.annotate_field == 'particle_mass':
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits, annotate_field=part_info['particle_mass'])
                elif 'angular' in args.annotate_field:
                    annotate_field = dd[args.annotate_field].value
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits, annotate_field=annotate_field, field_symbol='L', units=args.c_units)
                else:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits)
                for line in ax.xaxis.get_ticklines():
                    line.set_color('white')
                for line in ax.yaxis.get_ticklines():
                    line.set_color('white')

                fig.savefig(save_image_name)
                print "created slice plot:", save_image_name
            else:
                print "CREATING PICKLE"
                pickle_file = save_dir + 'slice_pickle.pkl'
                file = open(pickle_file, 'w+')
                pickle.dump((X_vel, Y_vel, xy, field_grid, weight_field, velx, vely, magx_grid, magy_grid, part_info, limits, time_string),file)
                file.close()
                print "Created Pickle"

        if args.profile_plot == 'True':
            if args.calc_column_mean != "False":
                time_stamp = (ds.current_time.value - np.min(dd['particle_creation_time'].value))/yt.units.yr.in_units('s').value
                time_string = '$t=$'+str(int(np.round(time_stamp/1000.)*1000.)) + 'yr'
                pickle_file = path + args.field + '_movie_time_' + str(m_times[fit]) + '.0_unweighted.pkl'
                file = open(pickle_file, 'r')
                X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo = pickle.load(file)
                file.close()
                #calc center of mass:
                image = image.in_units('g/cm**2')
                X = yt.YTArray(X, 'au')
                Y = yt.YTArray(Y, 'au')
                x_dim = (np.max(X) - np.min(X)).in_units('cm')
                y_dim = (np.max(Y) - np.min(Y)).in_units('cm')
                cube = ((np.max(X[0][-1]) - np.min(X[0][-2])).in_units('cm'))*((np.max(X[0][-1]) - np.min(X[0][-2])).in_units('cm'))
                TM = np.sum(image*(x_dim*y_dim))
                x_top = yt.YTArray(0.0, 'cm*g')
                y_top = yt.YTArray(0.0, 'cm*g')
                part_mass = yt.YTArray(part_info['particle_mass'], 'msun')
                part_pos_x = yt.YTArray(part_info['particle_position'][0], 'AU')
                part_pos_y = yt.YTArray(part_info['particle_position'][1], 'AU')
                TM = TM + np.sum(part_mass.in_units('g'))
                for part in range(len(part_info['particle_mass'])):
                    x_top = x_top + part_mass[part].in_units('g')*part_pos_x[part].in_units('cm')
                    y_top = y_top + part_mass[part].in_units('g')*part_pos_y[part].in_units('cm')
                x_top = x_top + np.sum(image*cube*X.in_units('cm'))
                y_top = y_top + np.sum(image*cube*Y.in_units('cm'))
                com = yt.YTArray([(x_top.value/TM.value), (y_top.value/TM.value), 0.0], 'cm')
                x_arr = yt.YTArray(np.sqrt((X-com[0].in_units('au'))**2. + (Y-com[1].in_units('au'))**2.), 'au')
                y_arr = image/(2.8*yt.units.mass_hydrogen.in_units('g'))
                w_arr = None
                z_arr = None
            else:
                myf.set_center(args.center)
                myf.set_coordinate_system('cylindrical')
                center_pos = dd['Center_Position']
                center_vel = dd['Center_Velocity'].value
                part_pos = dd['All_Particle_Positions']
                part_mass = dd['All_Particle_Masses']
                print "np.mean(dd['Particle_Potential'])", np.mean(dd['Particle_Potential'])
                print "Doing Profile Plot"
                save_image_name = save_dir + "Profile_Plot_time_" + str(args.plot_time) + ".pdf"
                L = [0.0, 0.0, 1.0]
                if myf.get_center() == 0:
                    measuring_volume = ds.disk(center_pos.value, L, (args.r_max*2., 'au'), (args.disk_thickness*2., 'au'))
                    tot_vec = [np.sum(measuring_volume['Angular_Momentum_x']).value, np.sum(measuring_volume['Angular_Momentum_y']).value, np.sum(measuring_volume['Angular_Momentum_z']).value]
                else:
                    tot_vec = [dd['particle_x_ang'][args.center-1].value, dd['particle_y_ang'][args.center-1].value, dd['particle_z_ang'][args.center-1].value]
                tot_mag = np.sqrt(tot_vec[0]**2. + tot_vec[1]**2. + tot_vec[2]**2.)
                L = tot_vec/tot_mag
                measuring_volume = ds.disk(center_pos.value, L, (args.r_max, 'au'), (args.disk_thickness, 'au'))
                if m_times[fit] == 5000:
                    disk_radii = [30, 35, 65]
                    sim_int = int(10*float(ds.directory.split('Mach_')[-1]))
                    disk_mass_measuring_volume = ds.disk(center_pos.value, L, (disk_radii[sim_int], 'au'), (args.disk_thickness, 'au'))
                    disk_mass = np.sum(disk_mass_measuring_volume['cell_mass'].in_units('msun'))
                    print "DISK MASS FOR MACH", sim_int/10., "IS", disk_mass
                if args.weight_field != None:
                    w_arr = measuring_volume[args.weight_field]
                else:
                    w_arr = args.weight_field
                x_arr = measuring_volume[args.x_field].in_units(args.x_units)
                if args.y_units == None:
                    y_arr = measuring_volume[args.field]
                    args.y_units = str(y_arr.units)
                else:
                    y_arr = measuring_volume[args.field].in_units(args.y_units)
                z_arr = measuring_volume['z'].in_units('AU')

            if args.field == "Tangential_Velocity":
                cal_vd_bool = True
                print "CALCULATING VELOCITY DISPERSION"
            else:
                cal_vd_bool = False

            if args.cumulative == 'False':
                cumulative = False
            else:
                cumulative = True

            prof_x, prof_y = mym.profile_plot(x_arr, y_arr, weight=w_arr, log=args.logscale, n_bins=args.profile_bins, bin_min=0.0, bin_max=args.r_max, calc_vel_dispersion=cal_vd_bool, cumulative=cumulative)
            if z_arr is not None:
                sampled_points = mym.sample_points(x_arr, y_arr, z_arr, bin_no=args.z_bins, no_of_points=args.no_sampled_points, weight_arr=w_arr)
            else:
                sampled_points = None
            
            if args.pickle_dump == "False":
                plt.clf()
                cm = plt.cm.get_cmap('RdYlBu')
                plot = plt.scatter(sampled_points[1], sampled_points[2], c=sampled_points[0], alpha=(1-w_arr/np.max(w_arr)), cmap=cm, edgecolors='none')
                if args.logscale != "False":
                    plt.semilogy(prof_x, prof_y, 'k-', linewidth=2.)
                else:
                    plt.plot(prof_x, prof_y, 'k-', linewidth=2.)
                cbar = plt.colorbar(plot, pad=0.0)
                cbar.set_label('|z position| (AU)', rotation=270, labelpad=13, size=14)
                plt.xlabel('Cyclindral Radius (AU)', labelpad=-1)
                plt.ylabel('Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)', labelpad=-3)
                plt.ylabel('$\theta$', labelpad=-3)
                plt.xlim([0.0, args.r_max])
                plt.ylim([0.0, 2.0])
                #plt.axhline(y=60.0, color='k', linestyle='--')
                #plt.axes().set_aspect(1./plt.axes().get_data_ratio())
                #plt.axes().set_aspect((args.r_max)/(2.0))
                plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
                print "created profile plot:", save_image_name
            else:
                if args.calc_column_mean != "False":
                    pickle_file = path + 'column_density_profile_pickle_'+str(args.plot_time)+'.pkl'
                else:
                    pickle_file = path + args.field + '_profile_pickle_'+str(args.plot_time)+'.pkl'
                file = open(pickle_file, 'w+')
                pickle.dump((prof_x, prof_y, sampled_points),file)
                file.close()
                print "created profile pickle:", pickle_file

        if args.profile_plot_multi == 'True':
            files = sorted(glob.glob(path + args.field + '_profile_pickle_*.pkl'))
            times = []
            plt.clf()
            colors = ['k', 'b', 'c', 'g', 'r', 'm']
            dash_list =  [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
            for f_it, file in enumerate(files):
                time_val = int(file.split('_profile_pickle_')[-1].split('.pkl')[0])
                times.append(time_val)
                open_file = open(file, 'r')
                prof_x, prof_y, sampled_points = pickle.load(open_file)
                open_file.close()
                plt.plot(prof_x, prof_y, c=colors[-len(files) + f_it], dashes=dash_list[-len(files) + f_it], label=str(time_val)+"yr")
            plt.legend(loc='best')
            plt.xlim([np.min(prof_x), np.max(prof_x)])
            plt.xlabel('Radius (AU)')
            plt.ylabel('velocity dispersion (km/s)')
            plt.savefig('velocity_dispersion.eps', bbox_inches='tight', pad_inches = 0.02)

        if args.b_mag == 'True':
            files = sorted(glob.glob(path + '*_plt_cnt*'))
            times = mym.generate_frame_times(files, args.time_step, start_time=0, end_time=args.endtime)
            plot_files = mym.find_files(times, files)
            if args.save_name == None:
                save_image_name = "angle_scatter.pdf"
            else:
                save_image_name = args.save_name + ".pdf"
            myf.set_center(0)
            plt.clf()
            x = []
            y = []
            fit = 0
            for file in plot_files:
                time = times[fit]
                part_file = file[:-12] + 'part' + file[-5:]
                ds = yt.load(file, particle_filename=part_file)
                dd = ds.all_data()
                measuring_volume = ds.disk(dd['Center_Position'], [0.0, 0.0, 1.0], (args.r_max+100, 'au'), (args.disk_thickness, 'au'))
                inds = np.where((measuring_volume['B_angle']>59.5)&(measuring_volume['B_angle']<60.5))[0]
                x.append(measuring_volume['Distance_from_Center'].in_units('AU')[inds])
                y.append(measuring_volume['dz_from_Center'].in_units('AU')[inds])
                if args.pickle_dump == "False":
                    plt.clf()
                    for time in range(len(x)):
                        plt.scatter(x[time], y[time], alpha=0.4, edgecolors='none')
                        plt.xlim([0., args.r_max])
                        plt.ylim([0.,50.])
                
                    plt.xlabel('Cylindrical distance (AU)')
                    plt.ylabel('Z-distance (AU)')
                    plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
                    print "created force comparison plot:", save_image_name
                fit = fit + 1
            if args.pickle_dump == "False":
                plt.clf()
                colors = ['k', 'b', 'c', 'g', 'r', 'm']
                dash_list =  [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
                #dash_list = [(None, None), (None, None), (None, None), (None, None), (None, None), (None, None)]
                for time in range(len(x)):
                    plt.scatter(x[time], y[time], c=colors[-len(x) + time], label=str(times[time])+"yr", alpha=0.4, edgecolors='none')
                    plt.xlim([0., args.r_max])
                    plt.ylim([0.,50.])

                plt.legend(loc='best')
                plt.xlabel('Cylindrical distance (AU)')
                plt.ylabel('Z-distance (AU)')
                plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
                print "created force comparison plot:", save_image_name
            else:
                pickle_file = save_dir + 'force_comp_pickle.pkl'
                file = open(pickle_file, 'w+')
                pickle.dump((x, y, times, args.y_label),file)
                file.close()
                print "created force comp pickle:", pickle_file

        if args.outflow_pickle == 'True':
            pickle_file = save_dir + "outflow_quantities.pkl"
            mass = []
            maximum_speed = []
            momentum = []
            ang_momentum = []
            height=500.0
            radius = 500.0
            thickness = 500.0
            files = sorted(glob.glob(path + '*_plt_cnt*'))
            times = mym.generate_frame_times(files, args.time_step, start_time=0)
            plot_files = mym.find_files(times, files)
            for file in plot_files:
                part_file = file[:-12] + 'part' + file[-5:]
                ds = yt.load(file, particle_filename=part_file)
                tube_1 = ds.disk([0.0, 0.0, height*yt.units.au.in_units('cm')], [0.0, 0.0, 1.0], (radius, 'au'), (thickness/2., 'au'))
                tube_2 = ds.disk([0.0, 0.0, -height*yt.units.au.in_units('cm')], [0.0, 0.0, 1.0], (radius, 'au'), (thickness/2., 'au'))
                outflow_pos_1 = np.where(tube_1['velocity_z'].value > 0.0)
                outflow_pos_2 = np.where(tube_2['velocity_z'].value < 0.0)

                mass_1 = tube_1['cell_mass'][outflow_pos_1].in_units('Msun')
                mass_2 = tube_2['cell_mass'][outflow_pos_2].in_units('Msun')
                outflow_mass = np.sum(mass_1) + np.sum(mass_2)
                
                if outflow_mass == 0.0:
                    outflow_mass = 0.0000001
                    max_speed = 0.0000001
                    mom = 0.0000001
                    L = 0.0000001
                else:
                    speed_1 = tube_1['velocity_magnitude'][outflow_pos_1].in_units("km/s")
                    mom_1 = speed_1*mass_1
                    speed_2 = tube_2['velocity_magnitude'][outflow_pos_2].in_units("km/s")
                    mom_2 = speed_2*mass_2

                    if len(speed_1) == 0:
                        max_speed_1 = yt.YTArray(0.0, 'km/s')
                    else:
                        max_speed_1 = np.max(speed_1)
                    if len(speed_2) == 0:
                        max_speed_2 = yt.YTArray(0.0, 'km/s')
                    else:
                        max_speed_2 = np.max(speed_2)
                    if max_speed_1 > max_speed_2:
                        max_speed = max_speed_1
                    else:
                        max_speed = max_speed_2
                    mom = np.sum(mom_1) + np.sum(mom_2)

                    L_x = np.sum(tube_1['Angular_Momentum_x'][outflow_pos_1].in_units("Msun * km**2 / s")) + np.sum(tube_2['Angular_Momentum_x'][outflow_pos_2].in_units("Msun * km**2 / s"))
                    L_y = np.sum(tube_1['Angular_Momentum_y'][outflow_pos_1].in_units("Msun * km**2 / s"))+ np.sum(tube_2['Angular_Momentum_y'][outflow_pos_2].in_units("Msun * km**2 / s"))
                    L_z = np.sum(tube_1['Angular_Momentum_z'][outflow_pos_1].in_units("Msun * km**2 / s")) + np.sum(tube_2['Angular_Momentum_z'][outflow_pos_2].in_units("Msun * km**2 / s"))
                    L = np.sqrt((L_x)**2. + (L_y)**2. + (L_z)**2.)

                mass.append(outflow_mass)
                maximum_speed.append(max_speed)
                momentum.append(mom)
                ang_momentum.append(L)
                print "outflow_mass", outflow_mass
                print "max_speed", max_speed
                print "mom", mom
                print "L", L
                print "DONE FILE", file
            file = open(pickle_file, 'w+')
            pickle.dump((times, mass, maximum_speed, momentum, ang_momentum),file)
            file.close()

        if args.yt_slice == 'True':
            title_parts = args.title.split('_')
            title = ''
            for part in title_parts:
                if part != title_parts[-1]:
                    title = title + part + ' '
                else:
                    title = title + part

            if args.projection_axis == 'xy':
                part_info = mym.get_particle_data(file, axis='xy')
            elif ('all', u'particle_mass') in ds.field_list:
                part_info = mym.get_particle_data(file)
                part_plane_position = np.array([dd['particle_posx'].in_units('AU'), dd['particle_posy'].in_units('AU')])
                part_info['particle_position'][0] = np.sign(part_plane_position[0])*np.sqrt((part_plane_position[0])**2. + (part_plane_position[1])**2.)
            print "MAKING YT SLICE FOR TIME =", m_times[fit]
            if args.image_center == 0 or ('all', u'particle_mass') not in ds.field_list:
                c = np.array([0.0, 0.0, 0.0])
            else:
                c = np.array([dd['particle_posx'][args.image_center-1].in_units('AU'), dd['particle_posy'][args.image_center-1].in_units('AU'), dd['particle_posz'][args.image_center-1].in_units('AU')])
            if ('all', u'particle_mass') not in ds.field_list:
                L = [0.0, 1.0, 0.0]
            elif len(dd['particle_posx']) == 1:
                L = [0.0, 1.0, 0.0]
            elif args.projection_orientation != None:
                y_val = 1./np.tan(np.deg2rad(args.projection_orientation))
                L = [1, y_val, 0]
                print "SET PROJECTION ORIENTATION L=", L
            elif args.projection_axis == 'xy':
                L = [0.0, 0.0, 1.0]
            else:
                pos_vec = [np.diff(dd['particle_posx'].value)[0], np.diff(dd['particle_posy'].value)[0]]
                L = [-1*pos_vec[-1], pos_vec[0]]
                L.append(0.0)
                if L[0] > 0:
                    L = [-1*L[0], -1*L[1], 0.0]
            myf.set_normal(L)
            if args.ax_lim != None:
                if args.image_center == 0 or ('all', u'particle_mass') not in ds.field_list:
                    xlim = [-1*args.ax_lim, args.ax_lim]
                    ylim = [-1*args.ax_lim, args.ax_lim]
                else:
                    xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center-1], args.ax_lim + part_info['particle_position'][0][args.image_center-1]]
                    ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center-1], args.ax_lim + part_info['particle_position'][1][args.image_center-1]]
            else:
                xlim = [-1000, 1000]
                ylim = [-1000, 1000]
            x_width = (xlim[1] -xlim[0])
            y_width = (ylim[1] -ylim[0])
            thickness = yt.YTArray(args.slice_thickness, 'AU')
            field = args.field
            for field_ind in ds.derived_field_list:
                if field_ind[1] == args.field:
                    full_field = field_ind
            temp = dd[full_field]
            temp = dd['velx']
            temp = dd['vely']
            temp = dd['velz']
            if ('all', u'particle_mass') in ds.field_list:
                temp = dd['particle_posx']
                temp = dd['particle_posy']
            temp = dd['velocity_magnitude']
            del temp

            if fit == 0 or usable_files[fit] != usable_files[fit-1]:
                L = np.array(L)
                if args.projection_axis == 'xz':
                    proj = yt.OffAxisProjectionPlot(ds, L, [field, 'Projected_Velocity_mw', 'velz_mw', 'Projected_Magnetic_Field_mw', 'magz_mw', 'cell_mass'], center=(c, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))
                    image_data = (proj.frb.data[full_field]/thickness.in_units('cm'))#.T#[:][::-1]
                    velx_full = (proj.frb.data[('gas', 'Projected_Velocity_mw')].in_units('g*cm**2/s')/thickness.in_units('cm'))#.T#[:][::-1]
                    vely_full = (proj.frb.data[('gas', 'velz_mw')].in_units('g*cm**2/s')/thickness.in_units('cm'))#.T#[:][::-1]
                    magx = (proj.frb.data[('gas', 'Projected_Magnetic_Field_mw')].in_units('g*gauss*cm')/thickness.in_units('cm'))#.T#[:][::-1]
                    magy = (proj.frb.data[('gas', 'magz_mw')].in_units('g*gauss*cm')/thickness.in_units('cm'))#.T#[:][::-1]
                    mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm'))#.T#[:][::-1]
                    velx_full = velx_full/mass
                    vely_full = vely_full/mass
                    magx = magx/mass
                    magy = magy/mass
                else:
                    proj = yt.OffAxisProjectionPlot(ds, L, [field, 'velx', 'vely', 'magx', 'magy', 'cell_mass'], center=(c, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))
                    image_data = (proj.frb.data[full_field]/thickness.in_units('cm'))#.T#[:][::-1]
                    velx_full = (proj.frb.data[('flash', u'velx')].in_units('cm**2/s')/thickness.in_units('cm'))#.T#[:][::-1]
                    vely_full = (proj.frb.data[('flash', u'vely')].in_units('cm**2/s')/thickness.in_units('cm'))#.T#[:][::-1]
                    magx = (proj.frb.data[('flash', u'magx')].in_units('gauss*cm')/thickness.in_units('cm'))#.T#[:][::-1]
                    magy = (proj.frb.data[('flash', u'magy')].in_units('gauss*cm')/thickness.in_units('cm'))#.T#[:][::-1]
                    mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm'))#.T#[:][::-1]
                '''
                if np.median(image_data[200:600,400]) > 1.e-15:
                    image_data = image_data.T
                    velx_full = velx_full.T
                    vely_full = vely_full.T
                    magx = magx.T
                    magy = magy.T
                    mass = mass.T
                '''
            del proj

            res = np.shape(image_data)[0]

            x = np.linspace(xlim[0], xlim[1], res)
            y = np.linspace(ylim[0], ylim[1], res)
            X, Y = np.meshgrid(x, y)

            cl = float((xlim[1] - xlim[0]))/float(res)
            x = np.arange(xlim[0]+cl/2., xlim[1]+cl/2., cl)
            annotate_freq = res/31.
            x_ind = []
            y_ind = []
            counter = 0
            while counter < 31:
                val = annotate_freq*counter + annotate_freq/2.
                x_ind.append(int(val))
                y_ind.append(int(val))
                counter = counter + 1
            x_vel = []
            y_vel = []
            for x_counter in x_ind:
                x_val = xlim[0]+cl/2. + x_counter*cl
                y_val = ylim[0]+cl/2. + x_counter*cl
                x_vel.append(x_val)
                y_vel.append(y_val)
            x_vel = np.array(x_vel)
            y_vel = np.array(y_vel)
            X_vel, Y_vel = np.meshgrid(x_vel, y_vel)

            velx, vely = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full)

            if args.pickle_dump == "False":
                fig, ax = plt.subplots()
                plot = ax.pcolormesh(X, Y, image_data.value, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=args.colourbar_min, vmax=args.colourbar_max), rasterized=True)
                plt.gca().set_aspect('equal')
                cbar = plt.colorbar(plot, pad=0.0)
                cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)

                plt.streamplot(X, Y, magx.value, magy.value, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)
                if ('all', u'particle_mass') in ds.field_list:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'])
                time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(m_times[fit]))+'yr', va="center", ha="left", color='w', fontsize=args.text_font)
                time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                #ax.set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
                ax.set_xlabel('Distance from center (AU)', labelpad=-1, fontsize=args.text_font)
                ax.set_ylabel('$z$ (AU)', labelpad=-20, fontsize=args.text_font)
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                if args.save_name == None:
                    save_image_name = save_dir + "yt_slice"
                else:
                    save_image_name = args.save_name
                plt.savefig(save_image_name+'.eps', format='eps', bbox_inches='tight')
                plt.savefig(save_image_name+'.pdf', format='pdf', bbox_inches='tight')
                print "SAVED PLOT AS", save_image_name
            else:
                args_dict = {}
                if args.annotate_time == "True":
                    args_dict.update({'annotate_time': '$t$='+str(int(args.plot_time))+'yr'})
                args_dict.update({'annotate_velocity': args.plot_velocity_legend})
                args_dict.update({'cbar_min': args.colourbar_min})
                args_dict.update({'cbar_max': args.colourbar_max})
                args_dict.update({'title': title})
                args_dict.update({'axlim':args.ax_lim})
                args_dict.update({'xlim':xlim})
                args_dict.update({'ylim':ylim})
                args_dict.update({'standard_vel':args.standard_vel})
                args_dict.update({'yabel':'$z$ (AU)'})

                pickle_file = save_dir + 'yt_proj_pickle.pkl'
                file = open(pickle_file, 'w+')
                pickle.dump((X, Y, X_vel, Y_vel, image_data.value, velx, vely, magx.value, magy.value, part_info, args_dict),file)
                file.close()
                print "created force comp pickle:", pickle_file


        #saving the figure
        if args.produce_movie != "False":
            save_image_name = save_dir + "movie_frame_" + ("%06d" % (args.start_frame + fit))
            plt.savefig(save_image_name + ".eps", format='eps', bbox_inches='tight')
            call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', save_image_name+'.eps', save_image_name+'.jpg'])
            os.remove(save_image_name + '.eps')
            print "CREATED MOVIE FRAME No.", args.start_frame + fit, "OF", no_of_frames, "saved as:", save_image_name

#=============================================================================
#These plots don't need to iterate over multiple files

if args.angular_momentum_budget == 'True':
    files = sorted(glob.glob(path + 'angular_momentum_budget_*.csv'))
    plt.clf()
    for file in files:
        times = []
        ang_labels = []
        ang_arrays = []
        header = True
        with open(file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if header == False:
                    times.append(float(row[0]))
                    for col_it in range(len(row[1:])):
                        ang_arrays[col_it].append(float(row[col_it+1]))
                if header == True:
                    for col in row[1:]:
                        ang_labels.append(col)
                        ang_arrays.append([])
                    header = False
        for arr in range(len(ang_arrays)):
            plt.semilogy(times, ang_arrays[arr], label=ang_labels[arr])
    plt.xlabel('Time since first sink particle formation (years)')
    plt.ylabel('Total specific angular momentum')
    plt.legend(loc='best')
    plt.savefig('Angular_momentum_budget.eps', bbox_inches='tight', pad_inches = 0.02)

if args.force_comp  == 'True':
    field = args.field
    files = sorted(glob.glob(path + '*_plt_cnt*'))
    #times = [0.0, 500.0, 1000.0, 1500.0]
    times = [2000.0, 3000.0, 4000.0, 5000.0]
    plot_files = mym.find_files(times, files)
    if args.save_name == None:
        save_image_name = save_dir + field + "_abs.eps"
    else:
        save_image_name = args.save_name + ".eps"
    myf.set_coordinate_system('sph')
    myf.set_center(args.center)
    plt.clf()
    x = []
    y = []
    fit = 0
    radius = 100
    height = 1000
    for file in plot_files:
        time = times[fit]
        part_file = file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        column = ds.disk(dd['Center_Position'], [0.0, 0.0, 1.0], (radius, 'au'), (height, 'au'))
        bin_data = yt.YTArray(np.linspace(0, 1000, 1001), 'AU')
        x_field = column['dz_from_Center'].in_units('AU')
        dummy = column['magx']
        dummy = column['magy']
        dummy = column['magz']
        if args.y_units == None:
            y_field = column[field]
            args.y_units = str(y_field.units)
        else:
            y_field = column[field].in_units(args.y_units)
        if np.max(y_field) < 0:
            y_field = np.abs(y_field)
        if args.weight_field != None:
            w_field = column[args.weight_field]
        else:
            w_field = None
        prof_x, prof_y= mym.profile_plot(x_field, y_field, weight=w_field, bin_data=bin_data, log=False, n_bins=args.profile_bins, bin_min=0.1)
        prof_x = np.array(prof_x)
        prof_y = np.array(prof_y)
        x.append(prof_x)
        y.append(prof_y)
        if args.pickle_dump == "False":
            plt.clf()
            plt.axhline(y=1.0, color='k', linestyle='--')
            for time in range(len(x)):
                if args.logscale == True:
                    plt.loglog(x[time], y[time], linewidth=1.5, alpha=0.75)
                    plt.xlim([1., height])
                else:
                    plt.plot(x[time], y[time], linewidth=1.5)
                    plt.xlim([0., height])
        
            plt.xlabel('Z-distance (AU)')
            if args.y_label == None:
                plt.ylabel(field+" ("+args.y_units+")")
        else:
            plt.ylabel(args.y_label, fontsize=14)
            plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
            print "created force comparison plot:", save_image_name
        fit = fit + 1
            
    if args.pickle_dump == "False":
        plt.clf()
        colors = ['k', 'b', 'c', 'g', 'r', 'm']
        dash_list = [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
        dash_list = [(None, None), (None, None), (None, None), (None, None), (None, None), (None, None)]
        for time in range(len(x)):
            plt.axhline(y=1.0, color='k', linestyle='--')
            if args.logscale == True:
                plt.loglog(x[time], y[time], c=colors[-len(x) + time], dashes=dash_list[-len(x) + time], label=str(times[time])+"yr", linewidth=1.5, alpha=0.75)
                plt.xlim([1., height])
            else:
                plt.plot(x[time], y[time], c=colors[-len(x) + time], dashes=dash_list[-len(x) + time], label=str(times[time])+"yr", linewidth=1.5)
                plt.xlim([0., height])
        
        plt.legend(loc='best')
        plt.xlabel('Z-distance (AU)')
        if args.y_label == None:
            plt.ylabel(field+" ("+args.y_units+")")
        else:
            plt.ylabel(args.y_label, fontsize=14)
        plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
        print "created force comparison plot:", save_image_name
    else:
        pickle_file = path + 'force_comp_pickle.pkl'
        file = open(pickle_file, 'w+')
        pickle.dump((x, y, times, args.y_label),file)
        file.close()
        print "created force comp pickle:", pickle_file

if args.separation == 'True':
    #image_name = save_dir + "separation"
    image_name = save_dir + "binary_system_time_evolution"
    files = ["/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.0/sinks_evol_copy.dat", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.1/sinks_evol_copy.dat", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.2/sinks_evol.dat"]
    csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
    line_style = ['k-', 'b-', 'r-']
    labels=["NT", "T1", "T2"]
    lit = 0
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(6, 8.)
    gs = gridspec.GridSpec(2, 1)
    #gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=0.0)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    #ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
    sim_times = []
    sim_total_mass = []
    for file in files:
        print "reading file", file
        sink_form_time = 0
        particle_tag = []
        times = [[],[]]
        x_pos = []
        y_pos = []
        z_pos = []
        mass = []
        with open(file, 'r') as f:
            reader = csv.reader(f, dialect='dat')
            counter = 0
            for row in reader:
                counter = counter
                if row[0][0] != '[':
                    if sink_form_time == 0:
                        sink_form_time = float(row[-1])
                    part_tag = int(row[0])
                    if part_tag == 65583:
                        part_tag = 65582
                    if part_tag not in particle_tag:
                        particle_tag.append(part_tag)
                        x_pos.append([])
                        y_pos.append([])
                        z_pos.append([])
                        mass.append([])
                        pit = len(particle_tag) - 1
                    elif particle_tag[0] == part_tag:
                        pit = 0
                    else:
                        pit = 1
                    x = float(row[2])/yt.units.AU.in_units('cm').value
                    y = float(row[3])/yt.units.AU.in_units('cm').value
                    z = float(row[4])/yt.units.AU.in_units('cm').value
                    m = float(row[14])/yt.units.msun.in_units('g').value
                    time = (float(row[1]) - sink_form_time)/yt.units.yr.in_units('s').value
                    #if time not in times[pit]:
                    times[pit].append(time)
                    x_pos[pit].append(x)
                    y_pos[pit].append(y)
                    z_pos[pit].append(z)
                    mass[pit].append(m)
                    '''
                    else:
                        ind = times[pit].index(time)
                        x_pos[pit][ind] = x
                        y_pos[pit][ind] = y
                        z_pos[pit][ind] = z
                        mass[pit][ind] = (m)
                    '''
        times = np.array(times)
        sorted_inds_1 = np.argsort(times[0])
        sorted_inds_2 = np.argsort(times[1])
        first_inds_1 = np.where((np.array(times[0])[sorted_inds_1][1:] - np.array(times[0])[sorted_inds_1][:-1])>0)[0]
        first_inds_2 = np.where((np.array(times[1])[sorted_inds_2][1:] - np.array(times[1])[sorted_inds_2][:-1])>0)[0]
        
        refined_time = []
        refined_time.append(np.array(times[0])[sorted_inds_1][first_inds_1])
        refined_time.append(np.array(times[1])[sorted_inds_2][first_inds_2])
        refined_x = []
        refined_x.append(np.array(x_pos[0])[sorted_inds_1][first_inds_1])
        refined_x.append(np.array(x_pos[1])[sorted_inds_2][first_inds_2])
        refined_y = []
        refined_y.append(np.array(y_pos[0])[sorted_inds_1][first_inds_1])
        refined_y.append(np.array(y_pos[1])[sorted_inds_2][first_inds_2])
        refined_z = []
        refined_z.append(np.array(z_pos[0])[sorted_inds_1][first_inds_1])
        refined_z.append(np.array(z_pos[1])[sorted_inds_2][first_inds_2])
        refined_mass = []
        refined_mass.append(np.array(mass[0])[sorted_inds_1][first_inds_1])
        refined_mass.append(np.array(mass[1])[sorted_inds_2][first_inds_2])
        dx = refined_x[0][-len(refined_x[1]):] - refined_x[1]
        dy = refined_y[0][-len(refined_y[1]):] - refined_y[1]
        dz = refined_z[0][-len(refined_z[1]):] - refined_z[1]
        sep = np.sqrt(dx**2. + dy**2. + dz**2.)
        total_mass = refined_mass[0].tolist()[:-len(refined_mass[1])] + (np.array(refined_mass[0][-len(refined_mass[1]):])+np.array(refined_mass[1])).tolist()
        mass_ratio = (np.nan*np.zeros(np.shape(refined_mass[0].tolist()[:-len(refined_mass[1])]))).tolist() +  ((np.array(refined_mass[1]))/(np.array(refined_mass[0][-len(refined_mass[1]):]))).tolist()
        dt = 25
        inds = [0]
        times_sort = [refined_time[0][0]]
        for t_it, time in enumerate(refined_time[0]):
            if time - times_sort[-1] > dt:
                times_sort.append(time)
                inds.append(t_it)
        times_sort = np.array(times_sort)
        total_mass_sort = np.array(total_mass)[inds]
        m_dot = (total_mass_sort[1:] - total_mass_sort[:-1])/(times_sort[1:]-times_sort[:-1])
        time_m_dot = (times_sort[:-1] + times_sort[1:])/2.
        
        ax1.semilogy(refined_time[1], sep, line_style[lit], label=labels[lit])
        #ax2.plot(refined_time[0], mass_ratio, line_style[lit], label=labels[lit])
        ax2.plot(time_m_dot, m_dot, line_style[lit], label=labels[lit])
        #ax3.plot(refined_time[0], refined_mass[0], line_style[lit], linewidth=1, alpha=0.5)
        #ax3.plot(refined_time[1], refined_mass[1], line_style[lit], linewidth=1, alpha=0.5)
        #ax3.plot(refined_time[0], total_mass, line_style[lit], linewidth=2, label=labels[lit])
        sim_times.append(times[1])
        sim_total_mass.append(total_mass)
        lit = lit + 1
    ax1.set_ylim([1e0,5e2])
    ax1.set_xlim([0.0, 5000.0])
    ax1.axhline(y=4.89593796548, linestyle='--', color='k', alpha=0.5)
    #ax1.legend(loc='best')
    #plt.xlabel("Time since formaton of first protostar (yr)", fontsize=14)
    ax1.set_ylabel("Separation (AU)", fontsize=args.text_font)
    ax1.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    ax1.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    ax1.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
    ax1.yaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    #ax3.yaxis.set_ticks_position('both')
    plt.setp([ax1.get_xticklabels() for ax1 in fig.axes[:-1]], visible=False)
    ax2.set_xlabel("Time since first protostar formation (yr)", fontsize=args.text_font)
    #ax2.set_ylabel("Mass ratio ($M_s/M_p$)", fontsize=args.text_font)
    #ax2.set_ylabel("Accreted Mass (M$_\odot$)", fontsize=args.text_font)
    ax2.set_ylabel("Accretion Rate (M$_\odot$/yr)", fontsize=args.text_font)
    ax2.legend(loc='best', fontsize=args.text_font)
    ax2.set_xlim([0, 5000])
    #ax2.set_ylim([1.e-10, 1.e-20])
    ax2.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    #ax3.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    #ax3.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    ax2.set_ylim([0.0, 0.0006])
    #ax3.set_ylim(bottom=0.0)
    plt.setp([ax2.get_yticklabels()[-1]], visible=False)
    plt.savefig(image_name + ".eps", bbox_inches='tight')
    plt.savefig(image_name + ".pdf", bbox_inches='tight')
    print "Created image", image_name

if args.force_on_particles == 'True':
    image_name = save_dir + "force_on_sinks"
    files = sorted(glob.glob(path+'*.csv'))
    time = []
    F_rad_tot = []
    F_tan_tot = []
    for file in files:
        time.append([])
        F_rad_tot.append([])
        F_tan_tot.append([])
        header = True
        with open(file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if header == False:
                    time_val = float(row[0])
                    Part_1_f_rad = float(row[1])
                    Part_1_f_tan = float(row[2])
                    if len(row) > 3:
                        Part_2_f_rad = float(row[3])
                        Part_2_f_tan = float(row[4])
                    else:
                        Part_2_f_rad = 0.0
                        Part_2_f_tan = 0.0
                    F_rad_val = Part_1_f_rad + Part_2_f_rad
                    F_tan_val = Part_1_f_tan + Part_2_f_tan
                    time[-1].append(time_val)
                    F_rad_tot[-1].append(F_rad_val)
                    F_tan_tot[-1].append(F_tan_val)
                if header == True:
                    header = False

    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(6, 7.)
    gs = gridspec.GridSpec(2, 1)
    gs.update(hspace=0.0)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    labels = ['Mach 0.0', 'Mach 0.1', 'Mach 0.2']
    for sim in range(len(time)):
        ax1.plot(time[sim], F_rad_tot[sim], label = labels[sim])
        ax2.plot(time[sim], F_tan_tot[sim], label = labels[sim])
    plt.xlabel('time since formation of first sink (yr)')
    ax1.set_ylabel('Radial force on particles')
    ax2.set_ylabel('Tangential force on particles')
    ax1.legend(loc='best')
    plt.savefig(image_name+'.eps')

if args.read_particle_file == 'True':
    file = path + 'sinks_evol.dat'
    csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
    pickle_file = save_dir + 'particle_data.pkl'
    if os.path.exists(pickle_file):
        file_open = open(pickle_file, 'r')
        particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
        file_open.close()
    else:
        init_line_counter = 0
        particle_data = {}
        particle_data.update({'particle_tag':[]})
        particle_data.update({'time':[]})
        particle_data.update({'posx':[]})
        particle_data.update({'posy':[]})
        particle_data.update({'posz':[]})
        particle_data.update({'velx':[]})
        particle_data.update({'vely':[]})
        particle_data.update({'velz':[]})
        particle_data.update({'accelx':[]})
        particle_data.update({'accely':[]})
        particle_data.update({'accelz':[]})
        particle_data.update({'anglx':[]})
        particle_data.update({'angly':[]})
        particle_data.update({'anglz':[]})
        particle_data.update({'mass':[]})
        particle_data.update({'mdot':[]})
        sink_form_time = 0
    line_counter = 0
    with open(file, 'r') as f:
        reader = csv.reader(f, dialect='dat')
        for row in reader:
            if line_counter >= init_line_counter:
                if row[0][0] != '[':
                    if sink_form_time == 0:
                        sink_form_time = float(row[-1])
                    time_val = (float(row[1])-sink_form_time)/yt.units.yr.in_units('s').value
                    part_tag = int(row[0])
                    if part_tag not in particle_data['particle_tag']:
                        print "adding particle", part_tag
                        particle_data['particle_tag'].append(part_tag)
                    if len(particle_data['particle_tag']) == 1:
                        part_ind = 0
                    else:
                        part_ind = particle_data['particle_tag'].index(part_tag)
                    if time_val not in particle_data['time']:
                        particle_data['time'].append(time_val)
                        particle_data['posx'].append([np.nan, np.nan])
                        particle_data['posy'].append([np.nan, np.nan])
                        particle_data['posz'].append([np.nan, np.nan])
                        particle_data['velx'].append([np.nan, np.nan])
                        particle_data['vely'].append([np.nan, np.nan])
                        particle_data['velz'].append([np.nan, np.nan])
                        particle_data['accelx'].append([np.nan, np.nan])
                        particle_data['accely'].append([np.nan, np.nan])
                        particle_data['accelz'].append([np.nan, np.nan])
                        particle_data['anglx'].append([np.nan, np.nan])
                        particle_data['angly'].append([np.nan, np.nan])
                        particle_data['anglz'].append([np.nan, np.nan])
                        particle_data['mass'].append([np.nan, np.nan])
                        particle_data['mdot'].append([np.nan, np.nan])
                    time_ind = particle_data['time'].index(time_val)
                    if np.isnan(particle_data['posx'][time_ind][part_ind]) == False:
                        print "replacing data at time", time_valparti
                    particle_data['posx'][time_ind][part_ind] = float(row[2])/yt.units.au.in_units('cm').value
                    particle_data['posy'][time_ind][part_ind] = float(row[3])/yt.units.au.in_units('cm').value
                    particle_data['posz'][time_ind][part_ind] = float(row[4])/yt.units.au.in_units('cm').value
                    particle_data['velx'][time_ind][part_ind] = float(row[5])
                    particle_data['vely'][time_ind][part_ind] = float(row[6])
                    particle_data['velz'][time_ind][part_ind] = float(row[7])
                    particle_data['accelx'][time_ind][part_ind] = float(row[8])
                    particle_data['accely'][time_ind][part_ind] = float(row[9])
                    particle_data['accelz'][time_ind][part_ind] = float(row[10])
                    particle_data['anglx'][time_ind][part_ind] = float(row[11])
                    particle_data['angly'][time_ind][part_ind] = float(row[12])
                    particle_data['anglz'][time_ind][part_ind] = float(row[13])
                    particle_data['mass'][time_ind][part_ind] = float(row[14])/yt.units.msun.in_units('g').value
                    particle_data['mdot'][time_ind][part_ind] = float(row[15])/(yt.units.msun.in_units('g').value/yt.units.yr.in_units('s').value)
                if np.remainder(line_counter,10000) == 0:
                    pickle_file = save_dir + 'particle_data.pkl'
                    file_open = open(pickle_file, 'w+')
                    pickle.dump((particle_data, sink_form_time, line_counter),file_open)
                    file_open.close()
                    print "dumped pickle after line", line_counter
            line_counter = line_counter + 1
        pickle_file = save_dir + 'particle_data.pkl'
        file_open = open(pickle_file, 'w+')
        pickle.dump((particle_data, sink_form_time, line_counter),file_open)
        file_open.close()
        print "dumped pickle after line", line_counter

if args.phasefolded_accretion == 'True':
    pickle_file = save_dir + 'particle_data.pkl'
    file_open = open(pickle_file, 'r')
    particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
    file_open.close()
    sorted_inds = np.argsort(particle_data['time'])
    particle_data['time'] = np.array(particle_data['time'])[sorted_inds]
    for key in particle_data.keys():
        if key != 'particle_tag' or key != 'time':
            particle_data[key] = np.array(particle_data[key]).T
            particle_data[key][0] = particle_data[key][0][sorted_inds]
            particle_data[key][1] = particle_data[key][1][sorted_inds]
