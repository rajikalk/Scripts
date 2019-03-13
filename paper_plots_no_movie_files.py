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
import shutil
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-fc", "--force_comp", help="Did you want to create a plot comparing pressure?", type=str, default="False")
    parser.add_argument("-fp", "--force_on_particles", help="Did you want to create a plot of the force on the particles?", type=str, default="False")
    parser.add_argument("-bp", "--b_mag", help="Did you want to create a plot of where the magnetic fiedl is 30 degrees", type=str, default="False")
    parser.add_argument("-op", "--outflow_pickle", help="Do you want to measure the outflows?", type=str, default="False")
    parser.add_argument("-sep", "--separation", help="Do you want to plot the separation of the particles?", type=str, default="False")
    parser.add_argument("-c", "--center", help="What center do you want to set for everything?, if 3 it combines all centers", type=int, default=0)
    parser.add_argument("-ppm", "--profile_plot_multi", help="Did you want to plot a profile plot with multiple lines?", type=str, default="False")
    
    #slice plot args
    parser.add_argument("-f", "--field", help="What is the field you want to have a sliceplot of?", type=str, default="Relative_Keplerian_Velocity")
    
    #profile plot args
    parser.add_argument("-rmax", "--r_max", help="radius of measuing volume", type=float, default=500.)
    parser.add_argument("-dt", "--disk_thickness", help="How far above and below the midplane do you want your profile to go?", type=float, default=100.)
    parser.add_argument("-wf", "--weight_field", help="any weight field?", type=str, default=None)
    parser.add_argument("-log", "--logscale", help="Want to use a log scale?", type=str, default="False")
    parser.add_argument("-pb", "--profile_bins", help="how many bins do you want for the profile?", type=int, default=None)
    parser.add_argument("-yu", "--y_units", help="y units for profile plot", type=str, default=None)
    parser.add_argument("-yl", "--y_label", help="what is the y_label you will use", type=str, default=None)
    parser.add_argument("-sn", "--save_name", help="If not defined, it's the same as the field", type=str, default=None)
    
    #pressure plot args
    parser.add_argument("-ts", "--time_step", help="would you like to plot multiple times?", type=int, default=500)
    
    parser.add_argument("-pd", "--pickle_dump", help="do you want to pickle data?", type=str, default="False")
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=12)
    parser.add_argument("-et", "--endtime", default=2001, type=int)
    
    parser.add_argument("-mov", "--produce_movie", help="Do you want to make a series of plots for a movie?", type=str, default="False")
    parser.add_argument("-sf", "--start_frame", help="if you don't want to start at frame 0, what frame do you want to start at?", type=int, default=0)
    
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
    pickle_file = path + 'particle_data.pkl'
    if os.path.exists(pickle_file):
        try:
            file_open = open(pickle_file, 'r')
            particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
            file_open.close()
        except:
            shutil.copy(pickle_file.split('.pkl')[0]+'_tmp.pkl',pickle_file)
            file_open = open(pickle_file, 'r')
            particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
            file_open.close()
        if 'velx' in particle_data.keys():
            particle_data.pop('velx')
            particle_data.pop('vely')
            particle_data.pop('velz')
            particle_data.pop('accelx')
            particle_data.pop('accely')
            particle_data.pop('accelz')
            particle_data.pop('anglx')
            particle_data.pop('angly')
            particle_data.pop('anglz')
    else:
        init_line_counter = 0
        particle_data = {}
        particle_data.update({'particle_tag':[]})
        particle_data.update({'time':[]})
        particle_data.update({'posx':[]})
        particle_data.update({'posy':[]})
        particle_data.update({'posz':[]})
        #particle_data.update({'velx':[]})
        #particle_data.update({'vely':[]})
        #particle_data.update({'velz':[]})
        #particle_data.update({'accelx':[]})
        #particle_data.update({'accely':[]})
        #particle_data.update({'accelz':[]})
        #particle_data.update({'anglx':[]})
        #particle_data.update({'angly':[]})
        #particle_data.update({'anglz':[]})
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
                        #particle_data['velx'].append([np.nan, np.nan])
                        #particle_data['vely'].append([np.nan, np.nan])
                        #particle_data['velz'].append([np.nan, np.nan])
                        #particle_data['accelx'].append([np.nan, np.nan])
                        #particle_data['accely'].append([np.nan, np.nan])
                        #particle_data['accelz'].append([np.nan, np.nan])
                        #particle_data['anglx'].append([np.nan, np.nan])
                        #particle_data['angly'].append([np.nan, np.nan])
                        #particle_data['anglz'].append([np.nan, np.nan])
                        particle_data['mass'].append([np.nan, np.nan])
                        particle_data['mdot'].append([np.nan, np.nan])
                    time_ind = particle_data['time'].index(time_val)
                    if np.isnan(particle_data['posx'][time_ind][part_ind]) == False:
                        print "replacing data at time", time_val
                    particle_data['posx'][time_ind][part_ind] = float(row[2])/yt.units.au.in_units('cm').value
                    particle_data['posy'][time_ind][part_ind] = float(row[3])/yt.units.au.in_units('cm').value
                    particle_data['posz'][time_ind][part_ind] = float(row[4])/yt.units.au.in_units('cm').value
                    #particle_data['velx'][time_ind][part_ind] = float(row[5])
                    #particle_data['vely'][time_ind][part_ind] = float(row[6])
                    #particle_data['velz'][time_ind][part_ind] = float(row[7])
                    #particle_data['accelx'][time_ind][part_ind] = float(row[8])
                    #particle_data['accely'][time_ind][part_ind] = float(row[9])
                    #particle_data['accelz'][time_ind][part_ind] = float(row[10])
                    #particle_data['anglx'][time_ind][part_ind] = float(row[11])
                    #particle_data['angly'][time_ind][part_ind] = float(row[12])
                    #particle_data['anglz'][time_ind][part_ind] = float(row[13])
                    particle_data['mass'][time_ind][part_ind] = float(row[14])/yt.units.msun.in_units('g').value
                    particle_data['mdot'][time_ind][part_ind] = float(row[15])/(yt.units.msun.in_units('g').value/yt.units.yr.in_units('s').value)
                if np.remainder(line_counter,10000) == 0:
                    file_open = open(pickle_file, 'w+')
                    pickle.dump((particle_data, sink_form_time, line_counter),file_open)
                    file_open.close()
                    print "dumped pickle after line", line_counter
                    shutil.copy(pickle_file, pickle_file.split('.pkl')[0]+'_tmp.pkl')
                    del particle_data
                    del sink_form_time
                    del line_counter
                    file_open = open(pickle_file, 'r')
                    particle_data, sink_form_time, line_counter = pickle.load(file_open)
                    file_open.close()
            line_counter = line_counter + 1
        file_open = open(pickle_file, 'w+')
        pickle.dump((particle_data, sink_form_time, line_counter),file_open)
        file_open.close()
        print "dumped pickle after line", line_counter
        shutil.copy(pickle_file, pickle_file.split('.pkl')[0]+'_tmp.pkl')
    print "sorting data"
    sorted_inds = np.argsort(particle_data['time'])
    particle_data['time'] = np.array(particle_data['time'])[sorted_inds]
    for key in particle_data.keys():
        if key != 'particle_tag' and key != 'time':
            particle_data[key] = np.array(particle_data[key]).T
            try:
                particle_data[key] = particle_data[key][:,sorted_inds]
            except:
                particle_data[key] = particle_data[key][sorted_inds]
    particle_data.update({'separation':np.sqrt((particle_data['posx'][0] - particle_data['posx'][1])**2. + (particle_data['posy'][0] - particle_data['posy'][1])**2. + (particle_data['posz'][0] - particle_data['posz'][1])**2.)})
    usable_inds = np.where(np.isnan(particle_data['separation']) == False)[0]
    for key in particle_data.keys():
        if key != 'particle_tag':
            try:
                particle_data[key] = particle_data[key][:,usable_inds]
            except:
                particle_data[key] = particle_data[key][usable_inds]
    print "sorted data and save"
    file_open = open(pickle_file, 'w+')
    pickle.dump((particle_data, sink_form_time, line_counter),file_open)
    file_open.close()

if args.phasefolded_accretion == 'True':
    pickle_file = path + 'particle_data.pkl'
    file_open = open(pickle_file, 'r')
    particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
    file_open.close()
    
    apsis_pickle = path + 'apsis_data.pkl'
    if os.path.exists(apsis_pickle):
        file_open = open(apsis_pickle, 'r')
        periastron_inds, apastron_inds = pickle.load(file_open)
        file_open.close()
    else:
        window = 10
        moving_index = window
        moving_average_time = []
        moving_average_sep = []
        moving_average_accretion = [[],[]]
        while moving_index < len(particle_data['time']):
            moving_average_time.append(np.mean(particle_data['time'][moving_index-window: moving_index]))
            moving_average_sep.append(np.mean(particle_data['separation'][moving_index-window: moving_index]))
            moving_average_accretion[0].append(np.mean(particle_data['mdot'][0][moving_index-window: moving_index]))
            moving_average_accretion[1].append(np.mean(particle_data['mdot'][1][moving_index-window: moving_index]))
            moving_index = moving_index + 1
        particle_data['time'] = np.array(moving_average_time)
        particle_data['separation'] = np.array(moving_average_sep)
        particle_data['mdot'] = np.array(moving_average_accretion)
        periastron_inds = [np.argmin(particle_data['separation'])]
        apastron_inds = []
        index_interval = 500
        index = periastron_inds[0] + index_interval
        passed_apastron = False
        while index < len(particle_data['time']):
            if np.min(particle_data['separation'][index-index_interval:index]) != np.min(particle_data['separation'][np.array([index-index_interval, index-index_interval+1, index-2, index-1])]) and passed_apastron == True:
                periastron_ind = np.argmin(particle_data['separation'][index-index_interval:index]) + index-index_interval
                periastron_inds.append(periastron_ind)
                print "found periastron at time", particle_data['time'][periastron_ind], "of separation", particle_data['separation'][periastron_ind]
                passed_apastron = False
            elif np.max(particle_data['separation'][index-index_interval:index]) != np.max(particle_data['separation'][np.array([index-index_interval, index-index_interval+1, index-2, index-1])]) and passed_apastron == False:
                apastron_ind = np.argmax(particle_data['separation'][index-index_interval:index]) + index-index_interval
                apastron_inds.append(apastron_ind)
                print "found apastron at time", particle_data['time'][apastron_ind], "of separation", particle_data['separation'][apastron_ind]
                passed_apastron = True
            index = index + index_interval/2
        plt.clf()
        plot_start_ind = 0
        plot_end_ind = len(periastron_inds)-1
        plt.plot(particle_data['time'][periastron_inds[plot_start_ind]:periastron_inds[plot_end_ind]], particle_data['separation'][periastron_inds[plot_start_ind]:periastron_inds[plot_end_ind]], 'k-')
        plt.xlabel('time (yr)')
        plt.ylabel('separation (au)')
        for periastron in periastron_inds[plot_start_ind:plot_end_ind]:
            plt.axvline(x=particle_data['time'][periastron], alpha=0.5)
        for apastron in apastron_inds[plot_start_ind:plot_end_ind]:
            plt.axvline(x=particle_data['time'][apastron], color='r', alpha=0.5)
        plt.savefig(save_dir + 'separation.eps')
        plt.clf()

        file_open = open(apsis_pickle, 'w+')
        pickle.dump((periastron_inds, apastron_inds),file_open)
        file_open.close()
        print "Created separation evolution plot"


    #periastron_inds = periastron_inds[10:]

    hist_ind = 1
    ap_ind = 0
    phase = np.linspace(0, 1, 20)
    #FIND BINNED DATA
    plt.clf()
    averaged_binned_accretion = [[],[]]
    averaged_total_accretion = []
    fig, axs = plt.subplots(2, 1, sharex=True)
    while hist_ind < len(periastron_inds):
        binned_accretion = [[],[]]
        total_accretion = []
        time_bins_1 = np.linspace(particle_data['time'][periastron_inds[hist_ind-1]], particle_data['time'][apastron_inds[ap_ind]],11)
        time_bins_2 = np.linspace(particle_data['time'][apastron_inds[ap_ind]], particle_data['time'][periastron_inds[hist_ind]],11)
        bin_ind = 1
        while bin_ind < len(time_bins_1):
            time_bin_inds_1 = np.where((particle_data['time'] > time_bins_1[bin_ind-1]) & (particle_data['time'] < time_bins_1[bin_ind]))[0]
            binned_accretion[0].append(np.sum(particle_data['mdot'][0][time_bin_inds_1]))
            binned_accretion[1].append(np.sum(particle_data['mdot'][1][time_bin_inds_1]))
            total_accretion.append(np.sum(particle_data['mdot'][:,time_bin_inds_1]))
            #binned_accretion[0].append(np.median(particle_data['mdot'][0][time_bin_inds]))
            #binned_accretion[1].append(np.median(particle_data['mdot'][1][time_bin_inds]))
            #total_accretion.append(np.median(particle_data['mdot'][:,time_bin_inds]))
            bin_ind = bin_ind + 1
        bin_ind = 1
        while bin_ind < len(time_bins_2):
            time_bin_inds_2 = np.where((particle_data['time'] > time_bins_2[bin_ind-1]) & (particle_data['time'] < time_bins_2[bin_ind]))[0]
            binned_accretion[0].append(np.sum(particle_data['mdot'][0][time_bin_inds_2]))
            binned_accretion[1].append(np.sum(particle_data['mdot'][1][time_bin_inds_2]))
            total_accretion.append(np.sum(particle_data['mdot'][:,time_bin_inds_2]))
            bin_ind = bin_ind + 1
        axs[0].plot(phase, binned_accretion[0], 'k-')
        axs[1].plot(phase, binned_accretion[1], 'k-')
        hist_ind = hist_ind + 1
        ap_ind = ap_ind + 1
        averaged_binned_accretion[0].append(binned_accretion[0])
        averaged_binned_accretion[1].append(binned_accretion[1])
        averaged_total_accretion.append(total_accretion)
    plt.xlabel("Phase")
    plt.xlim([0.0, 1.0])
    plt.savefig(save_dir + 'accretion.pdf')
    plt.clf()
    phase_2 = np.linspace(1, 2, 20)
    phase_2 = phase.tolist() + phase_2[1:].tolist()
    start_ind = 0
    ymax = None
    while start_ind+15 < len(averaged_binned_accretion[0]):
        plt.clf()
        median_accretion = []
        standard_deviation = []
        median_accretion.append(np.median(averaged_binned_accretion[0][start_ind:start_ind+15], axis=0))
        median_accretion.append(np.median(averaged_binned_accretion[1][start_ind:start_ind+15], axis=0))
        standard_deviation.append(np.std(averaged_binned_accretion[0][start_ind:start_ind+15], axis=0))
        standard_deviation.append(np.std(averaged_binned_accretion[1][start_ind:start_ind+15], axis=0))
        median_total = np.median(averaged_total_accretion[start_ind:start_ind+15], axis=0)
        standard_deviation_total = np.std(averaged_total_accretion[start_ind:start_ind+15], axis=0)
        #plt.plot(phase_2, median_accretion[0].tolist()+median_accretion[0][1:].tolist(), ls='steps', alpha=0.5, label='Primary')
        #plt.plot(phase_2, median_accretion[1].tolist()+median_accretion[1][1:].tolist(), ls='steps', alpha=0.5, label='Secondary')
        #plt.plot(phase_2, median_total.tolist() + median_total[1:].tolist(),ls='steps', label='Total')
        plt.errorbar(phase_2, median_accretion[0].tolist()+median_accretion[0][1:].tolist(), yerr=standard_deviation[0].tolist()+standard_deviation[0][1:].tolist(), ls='steps-mid', alpha=0.5, label='Primary')
        plt.errorbar(phase_2, median_accretion[1].tolist()+median_accretion[1][1:].tolist(), yerr=standard_deviation[1].tolist()+standard_deviation[1][1:].tolist(), ls='steps-mid', alpha=0.5, label='Secondary')
        plt.errorbar(phase_2, median_total.tolist() + median_total[1:].tolist(), yerr=standard_deviation_total.tolist()+standard_deviation_total[1:].tolist(), ls='steps-mid', label='Total')
        
        plt.legend(loc='best')
        plt.xlabel("Orbital Phase")
        plt.ylabel("Normalised accretion")
        plt.title("Start orbit:" + str(start_ind))
        plt.xlim([0.0, 1.3])
        if ymax == None:
            plt.ylim(bottom=0.0)
            ymax = plt.ylim()[1]
        else:
            plt.ylim([0.0, ymax])
        file_name = save_dir + 'accretion_median_start_orbit_'+ ("%02d" % (start_ind))
        plt.savefig(file_name +'.eps')
        call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
        start_ind = start_ind + 1
