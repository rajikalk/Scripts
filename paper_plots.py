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
import csv

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-mp", "--movie_plot", help="did you want to plot a slice of the movie?", default=False)
    parser.add_argument("-sp", "--slice_plot", help="Did you want to plot create a slice plot?", default=False)
    parser.add_argument("-pp", "--profile_plot", help="Did you want to plot a profile plot?", default=False)
    parser.add_argument("-fc", "--force_comp", help="Did you want to create a plot comparing pressure?", default=False)
    parser.add_argument("-op", "--outflow_pickle", help="Do you want to measure the outflows?", default=False)
    parser.add_argument("-po", "--plot_outflows", help="Do you want to plot the outflows nwo that you've measured them?", default=False)
    parser.add_argument("-sep", "--separation", help="Do you want to plot the separation of the particles?", default=False)
    parser.add_argument("-zt", "--zoom_times", help="4x is default zoom", default=4.1, type=float)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=int, default=0)
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="True")
    parser.add_argument("-c", "--center", help="What center do you want to set for everything?", type=int, default=0)
    parser.add_argument("-ic", "--image_center", help="Where would you like to center the image?", type=int, default=0)
    
    #movie plot args
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", default=True)
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=float, default=1.e-15)
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-13)
    
    #slice plot args
    parser.add_argument("-sf", "--slice_field", help="What is the field you want to have a sliceplot of?", type=str, default="Relative_Keplerian_Velocity")
    parser.add_argument("-res", "--resolution", help="What resolution do you want the slices to have?", type=int, default=1024)
    parser.add_argument("-af", "--annotate_field", help="What field of the particles do you want to annotate", type=str, default="particle_mass")
    parser.add_argument("-cc", "--combine_centers", help="Did you want to take the average of all centers for the slice plot?", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=None)
    
    #profile plot args
    parser.add_argument("-rmax", "--r_max", help="radius of measuign volume", type=float, default=500.)
    parser.add_argument("-dt", "--disk_thickness", help="How far above and below the midplane do you want your profile to go?", type=float, default=100.)
    parser.add_argument("-xf", "--x_field", help="x axis of the profile plot?", type=str, default="Distance_from_Center")
    parser.add_argument("-yf", "--y_field", help="y axis of the profile plot?", type=str, default="Relative_Keplerian_Velocity")
    parser.add_argument("-wf", "--weight_field", help="any weight field?", type=str, default="cell_mass")
    parser.add_argument("-log", "--logscale", help="Want to use a log scale?", type=bool, default=False)
    parser.add_argument("-pb", "--profile_bins", help="how many bins do you want for the profile?", type=int, default=None)
    parser.add_argument("-zb", "--z_bins", help="how many z bins do you want when sampling points?", type=int, default=2.)
    parser.add_argument("-nsp", "--no_sampled_points", help="how many random points do you want to randomly sample?", type=int, default=2000)
    parser.add_argument("-xy", "--x_units", help="x units for profile plot", type=str, default="AU")
    
    #pressure plot args
    parser.add_argument("-ts", "--time_step", help="would you like to plot multiple times?", type=int, default=500)
    
    parser.add_argument("-pd", "--pickle_dump", help="do you want to pickle data?", default=False)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-cu", "--c_units", help="What units do you want the olorbar in")
    
    #outflow plots
    
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

if args.movie_plot == 'True':
    args.movie_plot = True
if args.movie_plot:
    call(['mpirun', '-np', '16', 'python', '/home/100/rlk100/Scripts/movie_script_mod.py', path, save_dir, '-pt', str(args.plot_time), '-t', args.title, '-z', 'True', '-zt', str(args.zoom_times), '-cmin', str(args.colourbar_min), '-cmax', str(args.colourbar_max), '-at', str(args.annotate_time), '-ax', 'xy', '-wr', str(args.working_rank), '-pvl', str(args.plot_velocity_legend), '-ic', str(args.image_center)])

files = sorted(glob.glob(path + '*_plt_cnt*'))
file = mym.find_files([args.plot_time], files)[0]
print "FILE IS:", file
del files
part_file = file[:-12] + 'part' + file[-5:]
ds = yt.load(file, particle_filename=part_file)
dd = ds.all_data()

if args.slice_plot == 'True':
    args.slice_plot = True
if args.slice_plot:
    movie_files = sorted(glob.glob(path + 'WIND_slice_*'))
    movie_file = mym.find_files([args.plot_time], movie_files)[0]
    del movie_files
    n_bins = np.ceil(np.sqrt((args.resolution/2.)**2. + (args.resolution/2.)**2.))
    myf.set_n_bins(n_bins)
    save_image_name = save_dir + "Slice_Plot_time_" + str(args.plot_time) + ".eps"
    X, Y, X_vel, Y_vel = mym.initialise_grid(movie_file, zoom_times=args.zoom_times, center=args.image_center)
    if args.combine_centers == 'True':
        args.combine_centers = True
    myf.set_coordinate_system('cylindrical')
    fig, ax, xy, field_grid = mym.sliceplot(ds, X, Y, args.slice_field, resolution=args.resolution, comb=args.combine_centers, units=args.c_units)
    f = h5py.File(movie_file, 'r')
    velx, vely = mym.get_quiver_arrays(f['velx_slice_xy'][:,:,0], f['vely_slice_xy'][:,:,0])
    f.close()
    part_info = mym.get_particle_data(movie_file, axis='xy')
    if args.ax_lim == None:
        ax.set_xlim([np.min(X), np.max(X)])
        ax.set_ylim([np.min(Y), np.max(Y)])
    else:
        ax.set_xlim([-1*args.ax_lim, args.ax_lim])
        ax.set_ylim([-1*args.ax_lim, args.ax_lim])
    if args.pickle_dump == False:
        mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend)
        if args.annotate_field == 'particle_mass':
            mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], np.max(X), annotate_field=part_info['particle_mass'])
        elif 'angular' in args.annotate_field:
            annotate_field = dd[args.annotate_field].value
            mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], np.max(X), annotate_field=annotate_field, field_symbol='L', units=args.c_units)
        else:
            mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], np.max(X))

        fig.savefig(save_image_name)
        print "created slice plot:", save_image_name
    else:
        print "CREATING PICKLE"
        pickle_file = save_dir + 'slice_pickle.pkl'
        file = open(pickle_file, 'w+')
        pickle.dump((X_vel, Y_vel, xy, field_grid, velx, vely, part_info),file)
        file.close()
        print "Created Pickle"

if args.profile_plot == 'True':
    args.profile_plot = True
if args.profile_plot:
    print "Doing Profile Plot"
    save_image_name = save_dir + "Profile_Plot_time_" + str(args.plot_time) + ".pdf"
    measuring_volume = ds.disk(dd['Center_Position'], [0.0, 0.0, 1.0], (args.r_max, 'au'), (args.disk_thickness, 'au'))
    if args.profile_bins == None:
        args.profile_bins = args.resolution/2.
    prof_x, prof_y = mym.profile_plot(measuring_volume, args.x_field, [args.y_field], weight_field=args.weight_field, log=args.logscale, n_bins=args.profile_bins, x_units=args.x_units)
    sampled_points = mym.sample_points(measuring_volume, args.x_field, args.y_field, bin_no=args.z_bins, no_of_points=args.no_sampled_points, x_units=args.x_units)
    
    if args.pickle_dump == False:
        plt.clf()
        cm = plt.cm.get_cmap('RdYlBu')
        plot = plt.scatter(sampled_points[1], sampled_points[2], c=sampled_points[0], alpha=0.4, cmap=cm)
        plt.plot(prof_x, prof_y[args.y_field], 'k-', linewidth=2.)
        cbar = plt.colorbar(plot, pad=0.0)
        cbar.set_label('|z position| (AU)', rotation=270, labelpad=13, size=14)
        plt.xlabel('Cyclindral Radius (AU)', labelpad=-1)
        plt.ylabel('Relative Keplerian Velocity ($v_\phi$/$v_\mathrm{kep}$)', labelpad=-3)
        plt.xlim([0.0, args.r_max])
        plt.xlim([0.0, 2.0])
        plt.axhline(y=1.0, color='k', linestyle='--')
        plt.axes().set_aspect(1./plt.axes().get_data_ratio())
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
    field = 'Force_Ratio'
    files = sorted(glob.glob(path + '*_plt_cnt*'))
    times = mym.generate_frame_times(files, args.time_step, start_time=0)
    plot_files = mym.find_files(times, files)
    save_image_name = save_dir + field + "_abs.eps"
    myf.set_coordinate_system('sph')
    plt.clf()
    x = []
    y = []
    fit = 0
    for file in plot_files:
        time = times[fit]
        fit = fit + 1
        if time < 2500:
            part_file = file[:-12] + 'part' + file[-5:]
            ds = yt.load(file, particle_filename=part_file)
            dd = ds.all_data()
            column = ds.disk(dd['Center_Position'], [0.0, 0.0, 1.0], (25, 'au'), (1000, 'au'))
            prof_x, prof_y = mym.profile_plot(column, 'dz_from_Center', [field], weight_field=args.weight_field, log=args.logscale, n_bins=args.profile_bins, x_units=args.x_units)
            x.append(prof_x)
            prof_y = (np.array(prof_y[field][::-1]) + np.array(prof_y[field]))/2.
            y.append(prof_y)
    colors = ['k', 'b', 'c', 'g', 'r', 'm']
    dash_list =  [[1,3], [5,3,1,3,1,3,1,3], [5,3,1,3,1,3], [5,3,1,3], [5,5], (None, None)]
    for time in range(len(x)):
        plt.plot(x[time], y[time], c=colors[-len(x) + time], dashes=dash_list[-len(x) + time], label=str(times[time])+"yr")
        #plt.semilogy(x[time], y[time], c=colors[len(dash_list) - len(times) + time], dashes=dash_list[len(dash_list) - len(times) + time], label=str(times[time])+"yr")
    plt.legend(loc='best')
    plt.xlabel('Z-distance (AU)', labelpad=-1)
    plt.ylabel(field, labelpad=-3)
    plt.xlim([0.0, 1000.0])
    plt.axes().set_aspect(1./plt.axes().get_data_ratio())
    plt.savefig(save_image_name, bbox_inches='tight', pad_inches = 0.02)
    print "created profile plot:", save_image_name

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
    plt.xlabel("Time (yr)", labelpad=-1, fontsize=14)
    plt.ylabel("Separation (AU)", labelpad=-2, fontsize=14)
    plt.tick_params(axis='x', which='major', labelsize=14)
    plt.tick_params(axis='y', which='major', labelsize=14)
    plt.savefig(image_name, bbox_inches='tight')
    print "Created image", image_name
