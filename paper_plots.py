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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sp", "--slice_plot", help="Did you want to plot create a slice plot?", default=False)
    parser.add_argument("-pp", "--profile_plot", help="Did you want to plot a profile plot?", default=False)
    parser.add_argument("-fc", "--force_comp", help="Did you want to create a plot comparing pressure?", default=False)
    parser.add_argument("-bp", "--b_mag", help="Did you want to create a plot of where the magnetic fiedl is 30 degrees", default=False)
    parser.add_argument("-op", "--outflow_pickle", help="Do you want to measure the outflows?", default=False)
    parser.add_argument("-po", "--plot_outflows", help="Do you want to plot the outflows nwo that you've measured them?", default=False)
    parser.add_argument("-sep", "--separation", help="Do you want to plot the separation of the particles?", default=False)
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default=4.1, type=float)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int, default=0)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="True")
    parser.add_argument("-c", "--center", help="What center do you want to set for everything?, if 3 it combines all centers", type=int, default=0)
    parser.add_argument("-ic", "--image_center", help="Where would you like to center the image?", type=int, default=0)
    
    #movie plot args
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", default=True)
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-15)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-13)
    
    #slice plot args
    parser.add_argument("-f", "--field", help="What is the field you want to have a sliceplot of?", type=str, default="Relative_Keplerian_Velocity")
    parser.add_argument("-res", "--resolution", help="What resolution do you want the slices to have?", type=int, default=1024)
    parser.add_argument("-af", "--annotate_field", help="What field of the particles do you want to annotate", type=str, default="particle_mass")
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=None)
    
    #profile plot args
    parser.add_argument("-rmax", "--r_max", help="radius of measuign volume", type=float, default=500.)
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
    parser.add_argument("-st", "--slice_thickness", help="How thick do you want the slice to be?", default=100.0, type=float)
    
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
movie_file = mym.find_files([args.plot_time], movie_files)[0]
del movie_files
file_no = int(movie_file.split('WIND')[1][-6:])*4
file = movie_file.split('WIND')[0]+'WIND_hdf5_plt_cnt_' + ("%04d" % file_no)
if os.path.isfile(file) == False:
    sim_files = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt_*'))
    file = mym.find_files([args.plot_time], sim_files)[0]
part_file = file[:-12] + 'part' + file[-5:]
ds = yt.load(file, particle_filename=part_file)
dd = ds.all_data()
sink_form_time = np.min(dd['particle_creation_time'].value/yt.units.yr.in_units('s').value)

if args.slice_plot == 'True':
    args.slice_plot = True
if args.slice_plot:
    '''
    movie_files = sorted(glob.glob(path + 'WIND_slice_*'))
    movie_file = mym.find_files([args.plot_time], movie_files)[0]
    del movie_files
    file_no = int(movie_file.split('WIND')[1][-6:])*4
    file = movie_file.split('WIND')[0]+'WIND_hdf5_plt_cnt_' + ("%04d" % file_no)
    part_file = file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    dd = ds.all_data()
    '''
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
    myf.set_coordinate_system('cylindrical')
    dummy1, dummy2, dummy3, magx_grid = mym.sliceplot(ds, X, Y, 'magx', resolution=args.resolution, center=args.center, units=args.c_units)
    dummy1, dummy2, dummy3, magy_grid = mym.sliceplot(ds, X, Y, 'magy', resolution=args.resolution, center=args.center, units=args.c_units)
    del dummy1
    del dummy2
    del dummy3
    fig, ax, xy, field_grid = mym.sliceplot(ds, X, Y, args.field, resolution=args.resolution, center=args.center, units=args.c_units)
    f = h5py.File(movie_file, 'r')
    dim = np.shape(f['dens_slice_xy'])[0]
    xmin_full = f['minmax_xyz'][0][0]/yt.units.AU.in_units('cm').value
    xmax = f['minmax_xyz'][0][1]/yt.units.AU.in_units('cm').value
    cl = (xmax-xmin_full)/dim
    x_pos_min = int(np.round(np.min(X) - xmin_full))/cl
    y_pos_min = int(np.round(np.min(Y) - xmin_full))/cl
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
        pickle.dump((X_vel, Y_vel, xy, field_grid, velx, vely, part_info, limits),file)
        file.close()
        print "Created Pickle"

if args.profile_plot == 'True':
    args.profile_plot = True
if args.profile_plot:
    myf.set_center(args.center)
    print "Doing Profile Plot"
    save_image_name = save_dir + "Profile_Plot_time_" + str(args.plot_time) + ".pdf"
    measuring_volume = ds.disk(dd['Center_Position'], [0.0, 0.0, 1.0], (args.r_max+100, 'au'), (args.disk_thickness, 'au'))
    '''
    if args.profile_bins == None:
        args.profile_bins = args.resolution/2.
    '''
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
    sampled_points = mym.sample_points(x_arr, y_arr, z_arr, bin_no=args.z_bins, no_of_points=args.no_sampled_points)
    
    if args.pickle_dump == False:
        plt.clf()
        cm = plt.cm.get_cmap('RdYlBu')
        plot = plt.scatter(sampled_points[1], sampled_points[2], c=sampled_points[0], alpha=0.4, cmap=cm)
        plt.plot(prof_x, prof_y, 'k-', linewidth=2.)
        cbar = plt.colorbar(plot, pad=0.0)
        cbar.set_label('|z position| (AU)', rotation=270, labelpad=13, size=14)
        plt.xlabel('Cyclindral Radius (AU)', labelpad=-1)
        #plt.ylabel('Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)', labelpad=-3)
        plt.ylabel('$\theta$', labelpad=-3)
        plt.xlim([0.0, args.r_max])
        plt.ylim([0.0, 90.0])
        plt.axhline(y=60.0, color='k', linestyle='--')
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


if args.force_comp == 'True':
    args.force_comp = True
if args.force_comp:
    field = args.field
    files = sorted(glob.glob(path + '*_plt_cnt*'))
    times = [0.0, 500.0, 1000.0, 2000.0]
    plot_files = mym.find_files(times, files)
    if args.save_name == None:
        save_image_name = save_dir + field + "_abs.eps"
    else:
        save_image_name = args.save_name + ".eps"
    myf.set_coordinate_system('sph')
    myf.set_center(1)
    plt.clf()
    x = []
    y = []
    fit = 0
    for file in plot_files:
        time = times[fit]
        part_file = file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        height = 1000
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
        #prof_x = x_field
        #prof_y = y_field
        #prof_y = (prof_y + prof_y[::-1])/2.
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
            #plt.ylim([np.min(np.min(y)), np.max(np.max(y))])
            #plt.axes().set_aspect((1000.)/(plt.ylim()[-1]))
            #plt.axhline(y=1.0, color='k', linestyle='--')
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
        #plt.ylim([np.min(np.min(y)), np.max(np.max(y))])
        #plt.axes().set_aspect((1000.)/(plt.ylim()[-1]))
        plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
        print "created force comparison plot:", save_image_name
    else:
        pickle_file = save_dir + 'force_comp_pickle.pkl'
        file = open(pickle_file, 'w+')
        pickle.dump((x, y, times, args.y_label),file)
        file.close()
        print "created force comp pickle:", pickle_file

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
                plt.scatter(x[time], y[time], alpha=0.4)
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
            plt.scatter(x[time], y[time], c=colors[-len(x) + time], label=str(times[time])+"yr", alpha=0.4)
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
    image_name = save_dir + "separation.eps"
    files = ["/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.25/CircumbinaryOutFlow_0.25_lref_10/sinks_evol.dat", "/short/ek9/rlk100/Output/omega_t_ff_0.20/CircumbinaryOutFlow_0.50/CircumbinaryOutFlow_0.50_lref_10/sinks_evol.dat"]
    csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)
    line_style = ['b-', 'r-']
    labels=["Tight Binary", "Wide Binary"]
    lit = 0
    plt.clf()
    for file in files:
        particle_tag = []
        times = []
        x_pos = []
        y_pos = []
        z_pos = []
        with open(file, 'r') as f:
            reader = csv.reader(f, dialect='dat')
            for row in reader:
                if row[0][0] != '[':
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
                    time = (float(row[1]) - float(row[-1]))/yt.units.yr.in_units('s').value
                    times.append(time)
                    x = float(row[2])/yt.units.AU.in_units('cm').value
                    y = float(row[3])/yt.units.AU.in_units('cm').value
                    z = float(row[4])/yt.units.AU.in_units('cm').value
                    x_pos[pit].append(x)
                    y_pos[pit].append(y)
                    z_pos[pit].append(z)
        times = np.array(times)
        times = times[::2]
        x_pos = np.array(x_pos)
        y_pos = np.array(y_pos)
        z_pos = np.array(z_pos)
        dx = x_pos[0] - x_pos[1]
        dy = y_pos[0] - y_pos[1]
        dz = z_pos[0] - z_pos[1]
        sep = np.sqrt(dx**2. + dy**2. + dz**2.)
        plt.semilogy(times, sep, line_style[lit], label=labels[lit])
        lit = lit + 1
    plt.xlim([0.0, 3000.0])
    plt.legend(loc='best')
    plt.xlabel("Time (yr)", fontsize=14)
    plt.ylabel("Separation (AU)", fontsize=14)
    plt.tick_params(axis='x', which='major', labelsize=14)
    plt.tick_params(axis='y', which='major', labelsize=14)
    plt.savefig(image_name, bbox_inches='tight')
    print "Created image", image_name

if args.yt_slice == 'True':
    args.yt_slice = True
if args.yt_slice:
    print "MAKING YT SLICE"
    if args.image_center == 0:
        c = np.array([0.5, 0.5, 0.5])
    else:
        c = np.array([(dd['particle_posx'][args.image_center-1]-ds.domain_left_edge[0])/ds.domain_width[0], (dd['particle_posy'][args.image_center-1]-ds.domain_left_edge[1])/ds.domain_width[1], (dd['particle_posz'][args.image_center-1]-ds.domain_left_edge[2])/ds.domain_width[2]])
    if len(dd['particle_posx']) == 1:
        L = [0.0, 1.0, 0.0]
    else:
        pos_vec = [np.diff(dd['particle_posx'].value)[0], np.diff(dd['particle_posy'].value)[0]]
        L = [-1*pos_vec[-1], pos_vec[0]]
        L.append(0.0)
        if L[0] > 0:
            L = [-1*L[0], -1*L[1], 0.0]
    if args.ax_lim != None:
        if args.image_center == 0:
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
    temp = dd['particle_posx']
    temp = dd['particle_posy']
    temp = dd['velocity_magnitude']

    proj = yt.OffAxisProjectionPlot(ds, L, [field,'Projected_Velocity', 'velz'], center=c, width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'))
    image_data = (proj.frb.data[('flash', 'dens')]/thickness.in_units('cm')).T
    velx_full = (proj.frb.data[('gas', 'Projected_Velocity')].in_units('cm**2/s')/thickness.in_units('cm')).T
    vely_full = (proj.frb.data[('flash', 'velz')].in_units('cm**2/s')/thickness.in_units('cm')).T


    #weighted_proj = yt.OffAxisProjectionPlot(ds, L, ['Projected_Velocity', 'velz', 'magx', 'magz'], center=c, width=(x_width, 'AU'), depth=(args.slice_thickness, 'AU'), weight_field='density')

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
        val = xlim[0]+cl/2. + x_counter*cl
        x_vel.append(val)
        y_vel.append(val)
    x_vel = np.array(x_vel)
    y_vel = np.array(y_vel)
    X_vel, Y_vel = np.meshgrid(x_vel, y_vel)

    velx, vely = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full)

    part_info = mym.get_particle_data(file)
    part_plane_position = np.array([dd['particle_posx'].in_units('AU'), dd['particle_posy'].in_units('AU')])
    cen_AU = c*ds.domain_width.in_units('AU')[0].value - ds.domain_width.in_units('AU')[0].value/2.
    part_info['particle_position'][0] = np.sign(part_plane_position[0])*np.sqrt((part_plane_position[0] - cen_AU[0])**2. + (part_plane_position[1] - cen_AU[1])**2.)

    fig, ax = plt.subplots()
    plot = ax.pcolormesh(X, Y, image_data.value, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=1.e-16, vmax=1.e-14), rasterized=True)
    plt.gca().set_aspect('equal')
    cbar = plt.colorbar(plot, pad=0.0)
    cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)

    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)

    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'])
    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(args.plot_time))+'yr', va="center", ha="left", color='w', fontsize=12)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

    ax.set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
    ax.set_ylabel('$z$ (AU)', labelpad=-20, fontsize=args.text_font)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.savefig("yt_slice.eps", format='eps', bbox_inches='tight')
    print "SAVED PLOT AS yt_slice.eps"

