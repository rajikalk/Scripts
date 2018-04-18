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
    parser.add_argument("-sp", "--slice_plot", help="Did you want to plot create a slice plot?", default=False)
    parser.add_argument("-pp", "--profile_plot", help="Did you want to plot a profile plot?", default=False)
    parser.add_argument("-fc", "--force_comp", help="Did you want to create a plot comparing pressure?", default=False)
    parser.add_argument("-bp", "--b_mag", help="Did you want to create a plot of where the magnetic fiedl is 30 degrees", default=False)
    parser.add_argument("-op", "--outflow_pickle", help="Do you want to measure the outflows?", default=False)
    parser.add_argument("-po", "--plot_outflows", help="Do you want to plot the outflows now that you've measured them?", default=False)
    parser.add_argument("-sep", "--separation", help="Do you want to plot the separation of the particles?", default=False)
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default=4.1, type=float)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int, default=0)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-c", "--center", help="What center do you want to set for everything?, if 3 it combines all centers", type=int, default=0)
    parser.add_argument("-ic", "--image_center", help="Where would you like to center the image?", type=int, default=0)
    
    #movie plot args
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", default=True)
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
    parser.add_argument("-log", "--logscale", help="Want to use a log scale?", type=bool, default=False)
    parser.add_argument("-pb", "--profile_bins", help="how many bins do you want for the profile?", type=int, default=None)
    parser.add_argument("-zb", "--z_bins", help="how many z bins do you want when sampling points?", type=int, default=2.)
    parser.add_argument("-nsp", "--no_sampled_points", help="how many random points do you want to randomly sample?", type=int, default=2000)
    parser.add_argument("-xu", "--x_units", help="x units for profile plot", type=str, default="AU")
    parser.add_argument("-yu", "--y_units", help="y units for profile plot", type=str, default=None)
    parser.add_argument("-yl", "--y_label", help="what is the y_label you will use", type=str, default=None)
    parser.add_argument("-sn", "--save_name", help="If not defined, it's the same as the field", type=str, default=None)
    
    #pressure plot args
    parser.add_argument("-ts", "--time_step", help="would you like to plot multiple times?", type=int, default=500)
    
    parser.add_argument("-pd", "--pickle_dump", help="do you want to pickle data?", default=False)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-cu", "--c_units", help="What units do you want the colorbar in")
    parser.add_argument("-et", "--endtime", default=2001, type=int)
    parser.add_argument("-stdv", "--standard_vel", default=5, type=float)
    
    #yt_slice
    parser.add_argument("-ys", "--yt_slice", help="did you want to make a yt-slice?", default=False)
    parser.add_argument("-st", "--slice_thickness", help="How thick do you want the slice to be?", default=300.0, type=float)
    
    parser.add_argument("-mov", "--produce_movie", help="Do you want to make a series of plots for a movie?", default=False, type=bool)
    parser.add_argument("-sf", "--start_frame", help="if you don't want to start at frame 0, what frame do you want to start at?", type=int, default=0)
    parser.add_argument("-ef", "--end_frame", help="do you only want to plot up to t a certain number of frames?", default=None, type=int)
    
    parser.add_argument("-totm", "--total_mass_of_system", help="do you want to plot the total mass of your system?", default=False, type=bool)
    parser.add_argument("-ni", "--non_ideal", help="do you want to plot the non-ideal regimes we shoudl think about", default=False, type=bool)
    parser.add_argument("-proj_or", "--projection_orientation", help="Do you want to set the projection orientation? give as angle (in degrees) from positive y-axis", default=None, type=float)
    
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
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

myf.set_center(args.center)

movie_files = sorted(glob.glob(path + 'WIND_slice_*'))
sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt_*'))
file = sim_files[-1]
part_file = file[:-12] + 'part' + file[-5:]
ds = yt.load(file, particle_filename=part_file)
dd = ds.all_data()
sink_form_time = np.min(dd['particle_creation_time'].value/yt.units.yr.in_units('s').value)

if args.produce_movie != False:
    m_times = mym.generate_frame_times(sim_files, args.time_step, presink_frames=0)
else:
    m_times = [args.plot_time]
no_of_frames = len(m_times)
if args.end_frame != None:
    final_frame = args.end_frame + 1
    m_times = m_times[args.start_frame:final_frame]
else:
    m_times = m_times[args.start_frame:]

if args.slice_plot != False:
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
        if args.slice_plot != False:
            movie_file = usable_movie_files[fit]
        file = usable_files[fit]
        part_file = file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
    else:
        print "CREATING FRAME", args.start_frame + fit, "WITH FILE", usable_files[fit-1]

    if args.non_ideal == 'True':
        args.non_ideal = true
    if args.non_ideal:
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
        args.slice_plot = True
    if args.slice_plot:
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
        fig, ax, xy, field_grid, weight_field = mym.sliceplot(ds, X, Y, args.field, resolution=args.resolution, center=args.center, units=args.c_units, weight=args.weight_field)
        f = h5py.File(movie_file, 'r')
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
        if args.pickle_dump == False:
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
            pickle.dump((X_vel, Y_vel, xy, field_grid, weight_field, velx, vely, part_info, limits),file)
            file.close()
            print "Created Pickle"

    if args.profile_plot == 'True':
        args.profile_plot = True
    if args.profile_plot:
        myf.set_center(args.center)
        myf.set_coordinate_system('spherical')
        print "Doing Profile Plot"
        save_image_name = save_dir + "Profile_Plot_time_" + str(args.plot_time) + ".pdf"
        if myf.get_center() == 0:
             tot_vec = [np.sum(dd['specific_angular_momentum_x']).value, np.sum(dd['specific_angular_momentum_y']).value, np.sum(dd['specific_angular_momentum_z']).value]
        else:
            tot_vec = [dd['particle_x_ang'][args.center-1].value, dd['particle_y_ang'][args.center-1].value, dd['particle_z_ang'][args.center-1].value]
        tot_mag = np.sqrt(tot_vec[0]**2. + tot_vec[1]**2. + tot_vec[2]**2.)
        L = tot_vec/tot_mag
        measuring_volume = ds.disk(dd['Center_Position'], L, (args.r_max, 'au'), (args.disk_thickness, 'au'))
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

        prof_x, prof_y = mym.profile_plot(x_arr, y_arr, weight=w_arr, log=args.logscale, n_bins=args.profile_bins, bin_min=0.1)
        sampled_points = mym.sample_points(x_arr, y_arr, z_arr, bin_no=args.z_bins, no_of_points=args.no_sampled_points, weight_arr=w_arr)
        
        if args.pickle_dump == False:
            plt.clf()
            cm = plt.cm.get_cmap('RdYlBu')
            plot = plt.scatter(sampled_points[1], sampled_points[2], c=sampled_points[0], alpha=(1-w_arr/np.max(w_arr)), cmap=cm, edgecolors='none')
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
            pickle_file = save_dir + 'profile_pickle.pkl'
            file = open(pickle_file, 'w+')
            pickle.dump((prof_x, prof_y, sampled_points),file)
            file.close()
            print "created profile pickle:", pickle_file

    if args.b_mag == 'True':
        args.b_mag = True
    if args.b_mag:
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
            if args.pickle_dump == False:
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
        if args.pickle_dump == False:
            plt.clf()
            colors = ['k', 'b', 'c', 'g', 'r', 'm']
            dash_list = [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
            dash_list = [(None, None), (None, None), (None, None), (None, None), (None, None), (None, None)]
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
        args.outflow_pickle = True
    if args.outflow_pickle:
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
        args.yt_slice = True
    if args.yt_slice:
        title_parts = args.title.split('_')
        title = ''
        for part in title_parts:
            if part != title_parts[-1]:
                title = title + part + ' '
            else:
                title = title + part
        
        if ('all', u'particle_mass') in ds.field_list:
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
        else:
            pos_vec = [np.diff(dd['particle_posx'].value)[0], np.diff(dd['particle_posy'].value)[0]]
            L = [-1*pos_vec[-1], pos_vec[0]]
            L.append(0.0)
            if L[0] > 0:
                L = [-1*L[0], -1*L[1], 0.0]
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
        if args.field == "Relative_Keplerian_Velocity":
            field = "dens"
        else:
            field = args.field
        temp = dd['velx']
        temp = dd['vely']
        temp = dd['velz']
        if ('all', u'particle_mass') in ds.field_list:
            temp = dd['particle_posx']
            temp = dd['particle_posy']
        temp = dd['velocity_magnitude']

        if fit == 0 or usable_files[fit] != usable_files[fit-1]:
            L = np.array(L)
            proj = yt.OffAxisProjectionPlot(ds, L, [field, 'Projected_Velocity_mw', 'velz_mw', 'Projected_Magnetic_Field_mw', 'magz_mw', 'cell_mass'], center=(c, 'AU'), width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))
            image_data = (proj.frb.data[('flash', 'dens')]/thickness.in_units('cm')).T#[:][::-1]
            velx_full = (proj.frb.data[('gas', 'Projected_Velocity_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).T#[:][::-1]
            vely_full = (proj.frb.data[('gas', 'velz_mw')].in_units('g*cm**2/s')/thickness.in_units('cm')).T#[:][::-1]
            magx = (proj.frb.data[('gas', 'Projected_Magnetic_Field_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).T#[:][::-1]
            magy = (proj.frb.data[('gas', 'magz_mw')].in_units('g*gauss*cm')/thickness.in_units('cm')).T#[:][::-1]
            mass = (proj.frb.data[('gas', 'cell_mass')].in_units('cm*g')/thickness.in_units('cm')).T#[:][::-1]
            if np.median(image_data[200:600,400]) > 1.e-15:
                image_data = image_data.T
                velx_full = velx_full.T
                vely_full = vely_full.T
                magx = magx.T
                magy = magy.T
                mass = mass.T

        velx_full = velx_full/mass
        vely_full = vely_full/mass
        magx = magx/mass
        magy = magy/mass

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

        if args.pickle_dump == False:
            fig, ax = plt.subplots()
            plot = ax.pcolormesh(X, Y, image_data.value, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=1.e-16, vmax=1.e-14), rasterized=True)
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
                save_image_name = save_dir + "yt_slice.eps"
            else:
                save_image_name = args.save_name + ".eps"
            plt.savefig(save_image_name, format='eps', bbox_inches='tight')
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
    if args.produce_movie != False:
        save_image_name = save_dir + "movie_frame_" + ("%06d" % (args.start_frame + fit))
        plt.savefig(save_image_name + ".eps", format='eps', bbox_inches='tight')
        call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', save_image_name+'.eps', save_image_name+'.jpg'])
        os.remove(save_image_name + '.eps')
        print "CREATED MOVIE FRAME No.", args.start_frame + fit, "OF", no_of_frames, "saved as:", save_image_name

#=============================================================================
#These plots don't need to iterate over multiple files
'''
if args.measure_disks == 'True':
    args.measure_disks = True
if args.measure_disks:
    field = args.field
    files = sorted(glob.glob(path + '*_plt_cnt*'))
    times = [0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0]
    plot_files = mym.find_files(times, files)
    save_image_name = save_dir + field + "_profile_plot.eps"
    myf.set_center(args.center)
    height = 40
    radius = 20
    for file in plot_files:
        time = times[fit]
        part_file = file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
'''

if args.force_comp == 'True':
    args.force_comp = True
if args.force_comp:
    field = args.field
    files = sorted(glob.glob(path + '*_plt_cnt*'))
    times = [1000.0, 2000.0, 3000.0, 4000.0, 5000.0]
    plot_files = mym.find_files(times, files)
    if args.save_name == None:
        save_image_name = save_dir + field + "_abs.eps"
    else:
        save_image_name = args.save_name + ".eps"
    myf.set_coordinate_system('sph')
    #myf.set_center(1)
    plt.clf()
    x = []
    y = []
    fit = 0
    height = 1000
    for file in plot_files:
        time = times[fit]
        part_file = file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        column = ds.disk(dd['Center_Position'], [0.0, 0.0, 1.0], (50, 'au'), (height+100, 'au'))
        bin_data = column['dz_from_Center'].in_units('AU') - column['dz'].in_units('AU')
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
        prof_x, prof_y= mym.profile_plot(x_field, y_field, weight=w_field, log=args.logscale, n_bins=args.profile_bins, bin_data=bin_data, bin_min=0.1)
        prof_x = np.array(prof_x)
        prof_y = np.array(prof_y)
        x.append(prof_x)
        y.append(prof_y)
        if args.pickle_dump == False:
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
            
    if args.pickle_dump == False:
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
        pickle_file = save_dir + 'force_comp_pickle.pkl'
        file = open(pickle_file, 'w+')
        pickle.dump((x, y, times, args.y_label),file)
        file.close()
        print "created force comp pickle:", pickle_file

if args.plot_outflows == 'True':
    args.plot_outflows = True
if args.plot_outflows:
    files = glob.glob('/home/100/rlk100/Paper_plots/Outflows/*/outflow*')
    fig = plt.figure()
    gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=0.0)
    pit = 0
    times = []
    mass = []
    maximum_speed = []
    momentum = []
    ang_momentum = []
    max_time = 3000.
    for file in files:
        legend_label =  file.split('Outflows/')[-1].split('/')[0]
        f = open(file, 'r')
        t, m, ms, mom, ang_mom = pickle.load(f)
        times.append(t)
        mass.append(m)
        maximum_speed.append(ms)
        momentum.append(mom)
        ang_momentum.append(ang_mom)
        f.close()
    linestyles = ['k:', 'b-.', 'g--', 'r-', 'm+', 'cx']
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
    for file_no in range(len(times)):
        inds = np.where(mass[file_no] == 0.0)
        if len(inds) == 0:
            mass[file_no][inds] = 0.0000001
        ax1.semilogy(times[file_no], mass[file_no], linestyles[file_no])
    ax1.set_ylabel('Outflow Mass (M$_\odot$)')
    ax1.set_xlim([0, max_time])
    ax1.set_ylim(bottom=1.e-5)
    for file_no in range(len(times)):
        inds = np.where(momentum[file_no] == 0.0)
        if len(inds) == 0:
            momentum[file_no][inds] = 0.0000001
        ax2.semilogy(times[file_no], momentum[file_no], linestyles[file_no])
    ax2.set_ylabel('Outflow Momentum (M$_\odot$kms$^{-1}$)')
    ax2.set_xlim([0, max_time])
    ax2.set_ylim(bottom=1.e-5)
    for file_no in range(len(times)):
        inds = np.where(ang_momentum[file_no] == 0.0)
        if len(inds) == 0:
            ang_momentum[file_no][inds] = 0.0000001
        ax3.semilogy(times[file_no], ang_momentum[file_no], linestyles[file_no])
    ax3.set_xlabel('Time since protostar formation (yr)')
    ax3.set_ylabel('Outflow Ang. Mom. (M$_\odot$km$^2$s$^{-1}$)')
    ax3.set_xlim([0, max_time])
    ax3.set_ylim(bottom=1.e5)
    plt.tight_layout()
    image_name = save_dir + "outflow_plot.eps"
    plt.savefig(image_name, bbox_inches='tight')
    print "Created image", image_name

if args.separation == 'True':
    args.separation = True
if args.separation:
    image_name = save_dir + "separation"
    files = ["/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/CircumbinaryOutFlow_0.50_lref_10/sinks_evol.dat", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.1/sinks_evol.dat", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.2/sinks_evol.dat"]
    csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
    line_style = ['k-', 'b-', 'r-']
    labels=["Mach 0.0", "Mach 0.1", "Mach 0.2"]
    lit = 0
    plt.clf()
    for file in files:
        sink_form_time = 0
        particle_tag = []
        times = [[],[]]
        x_pos = []
        y_pos = []
        z_pos = []
        with open(file, 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if row[0][0] != '[':
                    if sink_form_time == 0:
                        sink_form_time = float(row[-1])
                    part_tag = int(row[0])
                    if part_tag not in particle_tag:
                        particle_tag.append(part_tag)
                        x_pos.append([])
                        y_pos.append([])
                        z_pos.append([])
                        pit = len(particle_tag) - 1
                    elif particle_tag[0] == part_tag:
                        pit = 0
                    else:
                        pit = 1
                    time = (float(row[1]) - sink_form_time)/yt.units.yr.in_units('s').value
                    times[pit].append(time)
                    x = float(row[2])/yt.units.AU.in_units('cm').value
                    y = float(row[3])/yt.units.AU.in_units('cm').value
                    z = float(row[4])/yt.units.AU.in_units('cm').value
                    x_pos[pit].append(x)
                    y_pos[pit].append(y)
                    z_pos[pit].append(z)
        times = np.array(times)
        sorted_inds_1 = np.argsort(times[0])
        sorted_inds_2 = np.argsort(times[1])
        times[0] = np.array(times[0])[sorted_inds_1]
        times[1] = np.array(times[1])[sorted_inds_2]
        
        x_pos = np.array(x_pos)
        y_pos = np.array(y_pos)
        z_pos = np.array(z_pos)
        x_pos[0] = np.array(x_pos[0])[sorted_inds_1]
        x_pos[1] = np.array(x_pos[1])[sorted_inds_2]
        y_pos[0] = np.array(y_pos[0])[sorted_inds_1]
        y_pos[1] = np.array(y_pos[1])[sorted_inds_2]
        z_pos[0] = np.array(z_pos[0])[sorted_inds_1]
        z_pos[1] = np.array(z_pos[1])[sorted_inds_2]
        x_pos[0] = x_pos[0][-len(x_pos[1]):]
        y_pos[0] = y_pos[0][-len(y_pos[1]):]
        z_pos[0] = z_pos[0][-len(z_pos[1]):]
        dx = x_pos[0] - x_pos[1]
        dy = y_pos[0] - y_pos[1]
        dz = z_pos[0] - z_pos[1]
        sep = np.sqrt(dx**2. + dy**2. + dz**2.)
        '''
        if file == "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/CircumbinaryOutFlow_0.50_lref_10/sinks_evol.dat":
            times_new = np.array([0])
            sep_new = np.array([500])
            times = np.append(times_new, times)
            sep = np.append(sep_new, sep)
        '''
        plt.semilogy(times[1], sep, line_style[lit], label=labels[lit])
        lit = lit + 1
    plt.xlim([0.0, 5000.0])
    plt.legend(loc='best')
    plt.xlabel("Time since formaton of first protostar (yr)", fontsize=14)
    plt.ylabel("Separation (AU)", fontsize=14)
    plt.tick_params(axis='x', which='major', labelsize=14)
    plt.tick_params(axis='y', which='major', labelsize=14)
    plt.savefig(image_name + '.eps', bbox_inches='tight')
    plt.savefig(image_name + '.pdf', bbox_inches='tight')
    print "Created image", image_name

if args.total_mass_of_system == 'True':
    args.total_mass_of_system = True
if args.total_mass_of_system:
    legend_labels = ['Mach 0.0', 'Mach 0.1', 'Mach 0.2']
    linestyles = ['k-', 'b--', 'r-.']
    image_name = save_dir + "total_mass"
    dirs = ["/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/CircumbinaryOutFlow_0.50_lref_10", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.1", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/Turbulent_sims/CircumbinaryOutFlow_0.50_lref_10/Mach_0.2"]
    case_it = 0
    plt.clf()
    for dir in dirs:
        files = sorted(glob.glob(dir + '/WIND_proj*'))
        m_times = mym.generate_frame_times(files, 100, presink_frames=0, end_time=None)
        usable_files = mym.find_files(m_times, files)
        sink_form_time = mym.find_sink_formation_time(files)
        times = []
        mass = []
        for file in usable_files:
            f = h5py.File(file, 'r')
            file_time = (f['time'][0]/yt.units.yr.in_units('s').value)-sink_form_time
            part_mass = np.sum(np.array(f['particlemasses'])/yt.units.msun.in_units('g').value)
            times.append(file_time)
            mass.append(part_mass)
            f.close()
        plt.plot(times, mass, linestyles[case_it], label=legend_labels[case_it])
        case_it = case_it + 1
    plt.xlabel("Time since first protostar formation (yr)", fontsize=args.text_font)
    plt.ylabel("Total accreted mass (M$_\odot$)", fontsize=args.text_font)
    plt.legend(loc='best')
    plt.xlim([0, 5000])
    plt.tick_params(axis='y', which='major', labelsize=args.text_font)
    plt.tick_params(axis='x', which='major', labelsize=args.text_font)
    plt.savefig(image_name + ".eps", bbox_inches='tight')
    plt.savefig(image_name + ".pdf", bbox_inches='tight')
    print "Created image", image_name


