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
from scipy.stats import skewnorm
import scipy.stats as stats
from scipy import optimize
from scipy.integrate import simps
from scipy.optimize import curve_fit
import math as math
import scipy.special as sp

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


Munoz = [0.6746769309013558, 0.7145197218701833,0.5569433848856828, 0.38752803014541026, 0.2541207203458393, 0.1926647710197189, 0.16721686663430102, 0.15719305545923135, 0.15231060868761137, 0.14225905633770175, 0.13740435074091906, 0.13252190396929908, 0.13278082160112703, 0.13819035069468377, 0.21555124026169104, 0.7145225050280866, 2.1646554314908557, 2.519659245902396, 1.8566914025475634, 0.9674850313244088, 0.6746769309013558, 0.7212266222807857, 0.5620755022308517, 0.39267864160713817, 0.25926208474928814, 0.1978153824814477, 0.17236747809602893, 0.15719305545923135]


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-fc", "--force_comp", help="Did you want to create a plot comparing pressure?", type=str, default="False")
    parser.add_argument("-fp", "--force_on_particles", help="Did you want to create a plot of the force on the particles?", type=str, default="False")
    parser.add_argument("-bp", "--b_mag", help="Did you want to create a plot of where the magnetic fiedl is 30 degrees", type=str, default="False")
    parser.add_argument("-op", "--outflow_pickle", help="Do you want to measure the outflows?", type=str, default="False")
    parser.add_argument("-sep", "--separation", help="Do you want to plot the separation of the particles?", type=str, default="False")
    parser.add_argument("-plot_e", "--plot_eccentricity", help="Do you want to plot the eccentricity with the separation?", type=str, default="False")
    parser.add_argument("-c", "--center", help="What center do you want to set for everything?, if 3 it combines all centers", type=int, default=0)
    parser.add_argument("-ppm", "--profile_plot_multi", help="Did you want to plot a profile plot with multiple lines?", type=str, default="False")
    
    #slice plot args
    parser.add_argument("-f", "--field", help="What is the field you want to have a sliceplot of?", type=str, default="Relative_Keplerian_Velocity")
    
    #profile plot args
    parser.add_argument("-rmax", "--r_max", help="radius of measuring volume", type=float, default=500.)
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
    parser.add_argument("-phase_multi", "--phasefolded_multi", help="do you want to plot the phasefolded accretion", type=str, default="False")
    parser.add_argument("-window", "--smoothing_window", help="How wide in years do you want the smoothing window to be?", default=1.0, type=float)
    parser.add_argument("-rep", "--repeats", help="how many phasefoles do you want to plot on one image?", type=int, default=1)
    parser.add_argument("-plot_reps", "--plot_repetitions", help="Do you want to plot multiple phasefolds next to each other?", type=str, default="False")
    parser.add_argument("-start_orb", "--starting_orbit", help="which orbit do you want to start folding at", type=int, default=1)
    parser.add_argument("-n_orb", "--n_orbits", help="How many orbits do you want to fold over?", type=int, default=5)
    parser.add_argument("-calc_a", "--calculate_apsis", help="do you want to calculate the apsis times?", type=str, default="False")
    parser.add_argument("-calc_e", "--calculate_eccentricity", help="do you want to calculate eccentricity", type=str, default="False")
    parser.add_argument("-overlap", "--overlapping_folds", help="do you want cosnsecutive folds to overlap in when orbits are used?", type=str, default="False")
    parser.add_argument("-skip", "--skip_iterator", help="If you do want to to overlap, but not overlap all orbits how many do you want to skip?", type=int, default=0)
    parser.add_argument("-p_beta", "--plot_beta", help="Do you want to plot the ratio of the max accretion against the base accretion?", type=str, default="False")
    parser.add_argument("-p_best", "--plot_best_fit", help="Do you want to plot the best fitting phasefolded curve", type=str, default="False")
    parser.add_argument("-method", "--beta_method", help="How do you want to calculate beta?", type=int, default=1)
    parser.add_argument("-res_study", "--resolution_study", help="plot resolution study?", type=str, default="False")
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def func(x, sigmag, mu, alpha, c, amp):
    #normal distribution
    normpdf = (1/(sigmag*np.sqrt(2*math.pi)))*np.exp(-(np.power((x-mu),2)/(2*np.power(sigmag,2))))
    normcdf = (0.5*(1+sp.erf((alpha*((x-mu)/sigmag))/(np.sqrt(2)))))
    return 2*amp*normpdf*normcdf + c

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
            print("created force comparison plot:", save_image_name)
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
        print("created force comparison plot:", save_image_name)
    else:
        pickle_file = save_dir + 'force_comp_pickle.pkl'
        file = open(pickle_file, 'w+')
        pickle.dump((x, y, times, args.y_label),file)
        file.close()
        print("created force comp pickle:", pickle_file)

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
        print("outflow_mass", outflow_mass)
        print("max_speed", max_speed)
        print("mom", mom)
        print("L", L)
        print("DONE FILE", file)
    file = open(pickle_file, 'w+')
    pickle.dump((times, mass, maximum_speed, momentum, ang_momentum),file)
    file.close()

#saving the figure
if args.produce_movie != "False":
    save_image_name = save_dir + "movie_frame_" + ("%06d" % (args.start_frame + fit))
    plt.savefig(save_image_name + ".eps", format='eps', bbox_inches='tight')
    call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', save_image_name+'.eps', save_image_name+'.jpg'])
    os.remove(save_image_name + '.eps')
    print("CREATED MOVIE FRAME No.", args.start_frame + fit, "OF", no_of_frames, "saved as:", save_image_name)

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
            print("created force comparison plot:", save_image_name)
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
        print("created force comparison plot:", save_image_name)
    else:
        pickle_file = path + 'force_comp_pickle.pkl'
        file = open(pickle_file, 'w+')
        pickle.dump((x, y, times, args.y_label),file)
        file.close()
        print("created force comp pickle:", pickle_file)

if args.separation == 'True':
    #image_name = save_dir + "separation"
    image_name = save_dir + "binary_system_time_evolution_Mach_0.2"
    line_style = ['b:', 'r-.', 'g--', 'c-']
    #labels=["T1", "T2"]
    labels=["$L_\mathrm{ref}$ = 11", "$L_\mathrm{ref}$ = 12", "$L_\mathrm{ref}$ = 13", "$L_\mathrm{ref}$ = 14"]
    lit = 0
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(6, 8.)
    if args.plot_eccentricity == 'True':
        #files = ["Mach_0.1/Lref_09/particle_data.pkl", "Mach_0.2/Lref_09/particle_data.pkl"]
        files = ["Mach_0.2/Lref_09/particle_data.pkl", "Mach_0.2/Lref_10/particle_data.pkl", "Mach_0.2/Lref_11/particle_data.pkl", "Mach_0.2/Lref_12/particle_data.pkl"]
        gs = gridspec.GridSpec(3, 1)
        #gs = gridspec.GridSpec(4, 1)
        gs.update(hspace=0.0)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
        ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
        #ax4 = fig.add_subplot(gs[3,0], sharex=ax1)
        for file in files:
            file_open = open(file, 'rb')
            particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
            file_open.close()
            ax1.semilogy(particle_data['time'][1:], particle_data['separation'][1:], line_style[lit], label=labels[lit])
            total_mass = np.array(particle_data['mass']).T[0] + np.array(particle_data['mass']).T[1]
            dt = 10
            prev_time = particle_data['time'][1]
            time_short = [particle_data['time'][1]]
            total_mass_short= [total_mass[1]]
            for time_it in range(len(particle_data['time'][1:])):
                if particle_data['time'][1:][time_it] - prev_time > dt:
                    time_short.append(particle_data['time'][1:][time_it])
                    total_mass_short.append(total_mass[time_it])
                    prev_time = particle_data['time'][1:][time_it]
            time_short = np.array(time_short)
            total_mass_short = np.array(total_mass_short)
            m_dot = (total_mass_short[1:] - total_mass_short[:-1])/(time_short[1:]-time_short[:-1])
            time_m_dot = (time_short[:-1] + time_short[1:])/2.
            
            ax2.semilogy(time_m_dot, m_dot, line_style[lit], label=labels[lit])
            ax3.semilogy(particle_data['time'][1:], particle_data['eccentricity'][1:], line_style[lit], label=labels[lit])
            #ax4.plot(particle_data['time'][1:], np.array(particle_data['mass']).T[1][1:]/np.array(particle_data['mass']).T[0][1:], line_style[lit], label=labels[lit])
            lit = lit + 1
        ax3.yaxis.set_ticks_position('both')
        ax3.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
        ax3.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
        ax3.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
        #ax4.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
        #ax4.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
        #ax4.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
        #ax1.set_xlim([0, 6500])
        ax3.set_ylim(top=1.3)
        ax2.set_ylim([1.e-6, 6.e-4])
        ax2.set_ylabel("Accretion Rate (M$_\odot$/yr)", fontsize=args.text_font)
        ax3.set_ylabel("Eccentricity", fontsize=args.text_font)
        #ax4.set_ylabel("Mass ratio ($q=M_\mathrm{p}/M_\mathrm{s}$)", fontsize=args.text_font)
        ax3.set_xlabel("Time since first protostar formation (yr)", fontsize=args.text_font)
        #ax4.set_xlabel("Time since first protostar formation (yr)", fontsize=args.text_font)
        plt.setp([ax1.get_xticklabels() for ax2 in fig.axes[:-1]], visible=False)
        #plt.setp([ax3.get_yticklabels()[-1]], visible=False)
        #plt.setp([ax3.get_yticklabels()[-1]], visible=False)
    else:
        gs = gridspec.GridSpec(2, 1)
        gs.update(hspace=0.0)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
        files = [path+'sinks_evol.dat']
        #files = ["/groups/astro/rlk/rlk/Simulations/Turbulent/Mach_0.1/Lref_10/sinks_evol.dat", "/groups/astro/rlk/rlk/Simulations/Turbulent/Mach_0.2/Lref_10/sinks_evol.dat"]
        csv.register_dialect('dat', delimiter=' ', skipinitialspace=True)

        #ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
        sim_times = []
        sim_total_mass = []
        for file in files:
            print("reading file", file)
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
            ax2.set_xlabel("Time since first protostar formation (yr)", fontsize=args.text_font)
            #ax3.plot(refined_time[0], refined_mass[0], line_style[lit], linewidth=1, alpha=0.5)
            #ax3.plot(refined_time[1], refined_mass[1], line_style[lit], linewidth=1, alpha=0.5)
            #ax3.plot(refined_time[0], total_mass, line_style[lit], linewidth=2, label=labels[lit])
            sim_times.append(times[1])
            sim_total_mass.append(total_mass)
            lit = lit + 1
    ax1.set_ylim(top=5e2)
    ax1.set_xlim(left=0.0)
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
    #ax2.set_ylabel("Mass ratio ($M_s/M_p$)", fontsize=args.text_font)
    #ax2.set_ylabel("Accreted Mass (M$_\odot$)", fontsize=args.text_font)
    ax2.set_ylabel("Accretion Rate (M$_\odot$/yr)", fontsize=args.text_font)
    ax1.legend(loc='best', fontsize=args.text_font)
    #ax2.set_xlim([0, 5000])
    #ax2.set_ylim([1.e-10, 1.e-20])
    ax2.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    #ax3.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    #ax3.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    #ax2.set_ylim([0.0, 0.006])
    #ax3.set_ylim(bottom=0.0)
    plt.setp([ax2.get_yticklabels()[-1]], visible=False)
    plt.savefig(image_name + ".eps", bbox_inches='tight', pad_inches=0.02)
    plt.savefig(image_name + ".pdf", bbox_inches='tight', pad_inches=0.02)
    print("Created image", image_name)

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
            file_open = open(pickle_file, 'rb')
            particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
            file_open.close()
        except:
            shutil.copy(pickle_file.split('.pkl')[0]+'_tmp.pkl',pickle_file)
            file_open = open(pickle_file, 'rb')
            particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
            file_open.close()
        if 'eccentricity' in list(particle_data.keys()):
            particle_data.pop('eccentricity')
        if 'accelx' in list(particle_data.keys()):
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
        particle_data.update({'velx':[]})
        particle_data.update({'vely':[]})
        particle_data.update({'velz':[]})
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
                        print("adding particle", part_tag)
                        particle_data['particle_tag'].append(part_tag)
                    if len(particle_data['particle_tag']) == 1:
                        part_ind = 0
                    else:
                        part_ind = particle_data['particle_tag'].index(part_tag)
                    if time_val not in particle_data['time']:
                        try:
                            particle_data['time'].append(time_val)
                            particle_data['posx'].append([np.nan, np.nan])
                            particle_data['posy'].append([np.nan, np.nan])
                            particle_data['posz'].append([np.nan, np.nan])
                            particle_data['velx'].append([np.nan, np.nan])
                            particle_data['vely'].append([np.nan, np.nan])
                            particle_data['velz'].append([np.nan, np.nan])
                            particle_data['mass'].append([np.nan, np.nan])
                            particle_data['mdot'].append([np.nan, np.nan])
                        except:
                            particle_data['time'] = np.append(particle_data['time'], time_val)
                            particle_data['posx'] = np.vstack((particle_data['posx'].T, [np.nan, np.nan])).T
                            particle_data['posy'] = np.vstack((particle_data['posy'].T, [np.nan, np.nan])).T
                            particle_data['posz'] = np.vstack((particle_data['posz'].T, [np.nan, np.nan])).T
                            particle_data['velx'] = np.vstack((particle_data['velx'].T, [np.nan, np.nan])).T
                            particle_data['vely'] = np.vstack((particle_data['vely'].T, [np.nan, np.nan])).T
                            particle_data['velz'] = np.vstack((particle_data['velz'].T, [np.nan, np.nan])).T
                            #particle_data['accelx'].append([np.nan, np.nan])
                            #particle_data['accely'].append([np.nan, np.nan])
                            #particle_data['accelz'].append([np.nan, np.nan])
                            #particle_data['anglx'].append([np.nan, np.nan])
                            #particle_data['angly'].append([np.nan, np.nan])
                            #particle_data['anglz'].append([np.nan, np.nan])
                            particle_data['mass'] = np.vstack((particle_data['mass'].T, [np.nan, np.nan])).T
                            particle_data['mdot'] = np.vstack((particle_data['mdot'].T, [np.nan, np.nan])).T
                    try:
                        time_ind = particle_data['time'].index(time_val)
                    except:
                        time_ind = particle_data['time'].tolist().index(time_val)
                    if np.isnan(particle_data['posx'][time_ind][part_ind]) == False:
                        print("replacing data at time", time_val)
                    particle_data['posx'][time_ind][part_ind] = float(row[2])/yt.units.au.in_units('cm').value
                    particle_data['posy'][time_ind][part_ind] = float(row[3])/yt.units.au.in_units('cm').value
                    particle_data['posz'][time_ind][part_ind] = float(row[4])/yt.units.au.in_units('cm').value
                    particle_data['velx'][time_ind][part_ind] = float(row[5])
                    particle_data['vely'][time_ind][part_ind] = float(row[6])
                    particle_data['velz'][time_ind][part_ind] = float(row[7])
                    #particle_data['accelx'][time_ind][part_ind] = float(row[8])
                    #particle_data['accely'][time_ind][part_ind] = float(row[9])
                    #particle_data['accelz'][time_ind][part_ind] = float(row[10])
                    #particle_data['anglx'][time_ind][part_ind] = float(row[11])
                    #particle_data['angly'][time_ind][part_ind] = float(row[12])
                    #particle_data['anglz'][time_ind][part_ind] = float(row[13])
                    particle_data['mass'][time_ind][part_ind] = float(row[14])/yt.units.msun.in_units('g').value
                    particle_data['mdot'][time_ind][part_ind] = float(row[15])/(yt.units.msun.in_units('g').value/yt.units.yr.in_units('s').value)
                if np.remainder(line_counter,10000) == 0:
                    file_open = open(pickle_file, 'wb')
                    pickle.dump((particle_data, sink_form_time, line_counter),file_open)
                    file_open.close()
                    print("dumped pickle after line", line_counter)
                    shutil.copy(pickle_file, pickle_file.split('.pkl')[0]+'_tmp.pkl')
                    file_open = open(pickle_file, 'rb')
                    particle_data, sink_form_time, line_counter = pickle.load(file_open)
                    file_open.close()
            line_counter = line_counter + 1
        file_open = open(pickle_file, 'wb')
        pickle.dump((particle_data, sink_form_time, line_counter),file_open)
        file_open.close()
        print("dumped pickle after line", line_counter)
        shutil.copy(pickle_file, pickle_file.split('.pkl')[0]+'_tmp.pkl')
    print("sorting data")
    sorted_inds = np.argsort(particle_data['time'])
    particle_data['time'] = np.array(particle_data['time'])[sorted_inds]
    for key in list(particle_data.keys()):
        if key != 'particle_tag' and key != 'time' and key != 'separation':
            particle_data[key] = np.array(particle_data[key])
            particle_data[key] = particle_data[key][sorted_inds]
    particle_data.update({'separation':np.sqrt((particle_data['posx'].T[0] - particle_data['posx'].T[1])**2. + (particle_data['posy'].T[0] - particle_data['posy'].T[1])**2. + (particle_data['posz'].T[0] - particle_data['posz'].T[1])**2.)})
    usable_inds = np.where(np.isnan(particle_data['separation']) == False)[0]
    for key in list(particle_data.keys()):
        if key != 'particle_tag':
            particle_data[key] = particle_data[key][usable_inds].tolist()
    print("sorted data and save")
    file_open = open(pickle_file, 'wb')
    pickle.dump((particle_data, sink_form_time, line_counter),file_open)
    file_open.close()

if args.calculate_eccentricity == 'True':
    pickle_file = path + 'particle_data.pkl'
    file_open = open(pickle_file, 'rb')
    particle_data, sink_form_time, line_counter = pickle.load(file_open)
    file_open.close()
    
    if 'eccentricity' in list(particle_data.keys()):
        e = particle_data['eccentricity']
        if len(e) == 0:
            calculate_e = True
        else:
            calculate_e = False
    else:
        calculate_e = True
    if calculate_e == True:
        Mass = yt.YTArray(np.array(particle_data['mass']).T, 'msun')
        velx = yt.YTArray(np.array(particle_data['velx']).T, 'cm/s')
        vely = yt.YTArray(np.array(particle_data['vely']).T, 'cm/s')
        velz = yt.YTArray(np.array(particle_data['velz']).T, 'cm/s')
        posx = yt.YTArray(np.array(particle_data['posx']).T, 'au')
        posy = yt.YTArray(np.array(particle_data['posy']).T, 'au')
        posz = yt.YTArray(np.array(particle_data['posz']).T, 'au')
        comx = (Mass[0].in_units('g')*posx[0].in_units('cm') + Mass[1].in_units('g')*posx[1].in_units('cm'))/(Mass[0].in_units('g')+Mass[1].in_units('g'))
        comy = (Mass[0].in_units('g')*posy[0].in_units('cm') + Mass[1].in_units('g')*posy[1].in_units('cm'))/(Mass[0].in_units('g')+Mass[1].in_units('g'))
        comz = (Mass[0].in_units('g')*posz[0].in_units('cm') + Mass[1].in_units('g')*posz[1].in_units('cm'))/(Mass[0].in_units('g')+Mass[1].in_units('g'))
        dist_from_com = np.sqrt((posx - comx)**2. + (posy - comy)**2. + (posz - comz)**2.)
        period = np.sqrt((4*np.pi*(dist_from_com[0]**3.)*((Mass[0]+Mass[1])**2.))/(yt.units.G.in_units('au**3/(msun*yr**2)')*(Mass[1]**3)))
        mu = yt.units.G*(Mass[0].in_units('g') + Mass[1].in_units('g'))
        dvx = velx[0] - velx[1]
        dvy = vely[0] - vely[1]
        dvz = velz[0] - velz[1]
        dv = yt.YTArray(np.array([dvx, dvy, dvz]).T, 'cm/s')
        drx = (posx[0] - posx[1]).in_units('cm')
        dry = (posy[0] - posy[1]).in_units('cm')
        drz = (posz[0] - posz[1]).in_units('cm')
        dr = yt.YTArray(np.array([drx, dry, drz]).T, 'cm/s')
        v = np.sqrt(dvx**2. + dvy**2. + dvz**2.)
        r_tot = np.sqrt(drx**2. + dry**2. + drz**2.)
        r_1 = [posx[0].in_units('cm')-comx.in_units('cm'), posy[0].in_units('cm')-comy.in_units('cm'), posz[0].in_units('cm')-comz.in_units('cm')]
        r_2 = [posx[1].in_units('cm')-comx.in_units('cm'), posy[1].in_units('cm')-comy.in_units('cm'), posz[1].in_units('cm')-comz.in_units('cm')]
        v_1 = [velx[0], vely[0], velz[0]]
        v_2 = [velx[1], vely[1], velz[1]]
        E_pot = -1*(yt.units.G*Mass[0].in_units('g')*Mass[1].in_units('g'))/r_tot
        E_kin = 0.5*Mass[0].in_units('g')*(velx[0]**2.+vely[0]**2.+velz[0]**2.) + 0.5*Mass[1].in_units('g')*(velx[1]**2.+vely[1]**2.+velz[1]**2.)
        reduced_mass = (Mass[0].in_units('g')*Mass[1].in_units('g'))/(Mass[0].in_units('g') + Mass[1].in_units('g'))
        epsilon = (E_pot + E_kin)/reduced_mass.in_units('g')
        #epsilon = (v**2.)/2. + mu/r
        m_x_r_1 = yt.YTArray(np.cross(np.array(r_1).T, np.array(v_1).T), 'cm**2/s').T
        m_x_r_2 = yt.YTArray(np.cross(np.array(r_2).T, np.array(v_2).T), 'cm**2/s').T
        L_1 = [Mass[0].in_units('g')*m_x_r_1[0], Mass[0].in_units('g')*m_x_r_1[1], Mass[0].in_units('g')*m_x_r_1[2]]
        L_2 = [Mass[1].in_units('g')*m_x_r_2[0], Mass[1].in_units('g')*m_x_r_2[1], Mass[1].in_units('g')*m_x_r_2[2]]
        L_tot = yt.YTArray(np.array(L_1), 'cm**2*g/s') + yt.YTArray(np.array(L_2), 'cm**2*g/s')
        h = np.sqrt(L_tot[0]**2. + L_tot[1]**2. + L_tot[2]**2.)/reduced_mass
        e = np.sqrt(1 + (2.*epsilon*h**2.)/(mu**2.))
        try:
            particle_data['eccentricity'] = e.tolist()
        except:
            particle_data.update({'eccentricity': e.tolist()})
        particle_data.update({'period': period.tolist()})
        file_open = open(pickle_file, 'wb')
        pickle.dump((particle_data, sink_form_time, line_counter),file_open)
        file_open.close()
    
    #plot separation and eccentricity
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(6, 8.)
    gs = gridspec.GridSpec(3, 1)
    #gs = gridspec.GridSpec(3, 1)
    gs.update(hspace=0.0)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax3 = fig.add_subplot(gs[2,0], sharex=ax1)
    ax1.semilogy(particle_data['time'], particle_data['separation'])
    ax2.semilogy(particle_data['time'], e)
    
    #Calculate smoothed accretion
    
    window = 500
    moving_index = window
    moving_average_time = []
    moving_average_accretion = [[],[]]
    particle_data['mdot'] = np.array(particle_data['mdot']).T
    while moving_index < len(particle_data['time']):
        moving_average_time.append(np.mean(particle_data['time'][moving_index-window: moving_index]))
        moving_average_accretion[0].append(np.mean(particle_data['mdot'][0][moving_index-window: moving_index]))
        moving_average_accretion[1].append(np.mean(particle_data['mdot'][1][moving_index-window: moving_index]))
        moving_index = moving_index + 1
    moving_average_time = np.array(moving_average_time)
    moving_average_accretion = np.array(moving_average_accretion)
    ax3.semilogy(moving_average_time, (moving_average_accretion[0]+moving_average_accretion[1]))
    '''
    smoothed_time = (particle_data['time'][::100][1:]+particle_data['time'][::100][:-1])/2.
    smoothed_accretion_1 = (particle_data['mass'][0,::100][1:]-particle_data['mass'][0,::100][:-1])/(particle_data['time'][::100][1:]-particle_data['time'][::100][:-1])
    smoothed_accretion_2 = (particle_data['mass'][1,::100][1:]-particle_data['mass'][1,::100][:-1])/(particle_data['time'][::100][1:]-particle_data['time'][::100][:-1])
    ax3.semilogy(smoothed_time, (smoothed_accretion_1+smoothed_accretion_2))
    '''
    ax3.set_xlabel("Time since first protostar formation (yr)", fontsize=args.text_font)
    ax1.set_ylabel("Separation (AU)", fontsize=args.text_font)
    ax2.set_ylabel("Eccentricity", fontsize=args.text_font)
    ax3.set_ylabel("Accretion (M$_\odot$/yr)", fontsize=args.text_font)
    ax1.set_ylim([1e0,5e2])
    ax2.set_ylim(top=2.0)
    ax3.set_ylim([1e-6,1e-3])
    ax1.set_xlim(left=0.0)
    #ax2.set_xlim([0, 5000])
    ax1.axhline(y=4.89593796548, linestyle='--', color='k', alpha=0.5)
    ax1.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    ax1.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    ax1.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
    ax2.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    ax3.tick_params(axis='y', which='major', labelsize=args.text_font, direction="in")
    ax3.tick_params(axis='y', which='minor', labelsize=args.text_font, direction="in")
    ax3.tick_params(axis='x', which='major', labelsize=args.text_font, direction="in")
    ax1.yaxis.set_ticks_position('both')
    ax2.yaxis.set_ticks_position('both')
    ax3.yaxis.set_ticks_position('both')
    plt.setp([ax1.get_xticklabels() for ax1 in fig.axes[:-1]], visible=False)
    plt.setp([ax2.get_xticklabels() for ax1 in fig.axes[:-1]], visible=False)
    plt.setp([ax3.get_yticklabels()[-1]], visible=False)
    
    apsis_pickle = path + 'apsis_data.pkl'
    file_open = open(apsis_pickle, 'rb')
    periastron_inds, apastron_inds = pickle.load(file_open)
    file_open.close()
    
    #e_max_1 = np.max(e[periastron_inds[0]:periastron_inds[6]])
    #e_min_1 = np.min(e[periastron_inds[0]:periastron_inds[6]])
    #ax2.fill_between([0, particle_data['time'][-1]], e_min_1, e_max_1, facecolor='grey', alpha=0.5)
    
    #e_max_2 = np.max(e[periastron_inds[-16]:periastron_inds[-1]])
    #e_min_2 = np.min(e[periastron_inds[-16]:periastron_inds[-1]])
    #ax2.fill_between([0, particle_data['time'][-1]], e_min_2, e_max_2, facecolor='grey', alpha=0.5)
    
    plt.savefig(save_dir+'system_evolution.eps', bbox_inches='tight')
    plt.savefig(save_dir+'system_evolution.pdf', bbox_inches='tight')

if args.calculate_apsis == 'True':
    pickle_file = path + 'particle_data.pkl'
    file_open = open(pickle_file, 'rb')
    particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
    file_open.close()
    
    apsis_pickle = path + 'apsis_data.pkl'
    if os.path.exists(apsis_pickle):
        file_open = open(apsis_pickle, 'rb')
        periastron_inds, apastron_inds = pickle.load(file_open)
        file_open.close()
    else:
        window = args.smoothing_window
        moving_index = int(window)
        moving_average_time = []
        moving_average_sep = []
        particle_data['time'] = np.array(particle_data['time'])
        particle_data['separation'] = np.array(particle_data['separation'])
        '''
        for window_int in range(window):
            particle_data['time'] = np.array(particle_data['time'])
            moving_average_time.append(particle_data['time'][moving_index-int(window/2):-int(window)])
            moving_average_sep.append(particle_data['separation'][moving_index-int(window/2):-int(window)])
        '''
        while moving_index < len(particle_data['time']):
            window_inds = np.where((particle_data['time'] < particle_data['time'][moving_index])&(particle_data['time'] > particle_data['time'][moving_index]-window))[0]
            moving_average_time.append(np.mean(particle_data['time'][window_inds]))
            moving_average_sep.append(np.mean(particle_data['separation'][window_inds]))
            moving_index = moving_index + 1
        particle_data['time'] = np.array(moving_average_time)
        particle_data['separation'] = np.array(moving_average_sep)
        periastron_inds = [np.argmin(particle_data['separation'][:np.argmin(abs(particle_data['time']-1550))])]
        apastron_inds = []
        index_interval = 1000
        index = periastron_inds[0] + index_interval
        passed_apastron = False
        while index < len(particle_data['time']):
            if np.min(particle_data['separation'][index-index_interval:index]) != np.min(particle_data['separation'][np.array([index-index_interval, index-1])]) and passed_apastron == True:
                periastron_ind = np.argmin(particle_data['separation'][index-index_interval:index]) + index-index_interval
                periastron_inds.append(periastron_ind)
                print("found periastron at time", particle_data['time'][periastron_ind], "of separation", particle_data['separation'][periastron_ind])
                passed_apastron = False
            elif np.max(particle_data['separation'][index-index_interval:index]) != np.max(particle_data['separation'][np.array([index-index_interval, index-1])]) and passed_apastron == False:
                apastron_ind = np.argmax(particle_data['separation'][index-index_interval:index]) + index-index_interval
                apastron_inds.append(apastron_ind)
                print("found apastron at time", particle_data['time'][apastron_ind], "of separation", particle_data['separation'][apastron_ind])
                passed_apastron = True
            index = index + int(index_interval/2)
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
        #plt.xlim(left=4000)
        plt.savefig(save_dir + 'separation.eps')
        plt.clf()

        file_open = open(apsis_pickle, 'wb')
        pickle.dump((periastron_inds, apastron_inds),file_open)
        file_open.close()
        print("Created separation evolution plot")

if args.phasefolded_accretion == 'True':
    use_e_bins = True
    e_bins = [1.1, 0.6]#, 0.4, 0.2, 0.0]
    #e_bins = [1.1, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]
    #e_bins = [1.1, 0.7, 0.5, 0.3, 0.1, 0.0]

    pickle_file = path + 'particle_data.pkl'
    file_open = open(pickle_file, 'rb')
    particle_data, sink_form_time, init_line_counter = pickle.load(file_open)
    file_open.close()
    
    apsis_pickle = path + 'apsis_data.pkl'
    file_open = open(apsis_pickle, 'rb')
    periastron_inds, apastron_inds = pickle.load(file_open)
    if periastron_inds[0] == 0:
        periastron_inds = periastron_inds[1:]
    if apastron_inds[0] == 0:
        apastron_inds = apastron_inds[1:]
    file_open.close()
    n_folded_orbits = args.n_orbits
    
    file_name = save_dir + 'multiple_folds_over_' + str(n_folded_orbits) + '_orbits'
    folded_pickle = path + file_name + '.pkl'
    repeats = args.repeats
    
    if use_e_bins:
        particle_data['time'] = np.array(particle_data['time'])
        particle_data['mdot'] = np.array(particle_data['mdot']).T
    
        e_bin_it = 1
        plot_it = 0
        phase = np.linspace(0.0, 1.0, 21)
        phase_centers = (phase[1:] + phase[:-1])/2.
        phase_pre = np.linspace(-1.0, 0.0, 21)
        phase_2 = np.linspace(1.0, 2.0, 21)
        phase_2 = phase_pre.tolist()[:-1] + phase.tolist() + phase_2[1:].tolist()
        phase_centers = (np.array(phase_2[1:]) + np.array(phase_2[:-1]))/2.
        
        median_eccentricity = []
        std_eccentricity = []
        
        multiple_folds = []
        accretion_err = []
        y_fits = []
        
        plt.clf()
        fig = plt.figure()
        columns = 4
        rows = 2#int((len(e_bins)-1)/4)
        fig.set_size_inches(3.25*columns, 3.25*rows)
        gs = gridspec.GridSpec(rows, columns)
        
        gs.update(wspace=0.0, hspace=0.0)
        axes_dict = {}
        
        while e_bin_it < len(e_bins):
            file_name = save_dir + 'accretion_median_start_orbit_from_' + str(e_bins[e_bin_it-1]) + '_' + str(e_bins[e_bin_it])
            usable_periastrons = np.where((np.array(particle_data['eccentricity'])[periastron_inds]<e_bins[e_bin_it-1])&(np.array(particle_data['eccentricity'])[periastron_inds]>e_bins[e_bin_it]))[0]
            median_e = np.median(particle_data['eccentricity'][periastron_inds[usable_periastrons[0]]:periastron_inds[usable_periastrons[-1]]])
            std_e = np.std(particle_data['eccentricity'][periastron_inds[usable_periastrons[0]]:periastron_inds[usable_periastrons[-1]]])
            mean_e = np.mean(particle_data['eccentricity'][periastron_inds[usable_periastrons[0]]:periastron_inds[usable_periastrons[-1]]])
            median_eccentricity.append(median_e)
            std_eccentricity.append([median_e-(mean_e-std_e), (mean_e+std_e)-median_e])
            
            rit = 1
            averaged_binned_accretion = [[],[]]
            averaged_total_accretion = []
            print("e=["+str(e_bins[e_bin_it-1])+','+str(e_bins[e_bin_it])+"] --> "+str(len(usable_periastrons)) + " orbits")
            while rit < len(usable_periastrons):
                binned_accretion = [[],[]]
                total_accretion = []

                time_bins_1 = sorted(np.linspace(particle_data['time'][periastron_inds[usable_periastrons[rit-1]]], particle_data['time'][apastron_inds[usable_periastrons[rit-1]]],11))
                time_bins_2 = sorted(np.linspace(particle_data['time'][apastron_inds[usable_periastrons[rit-1]]], particle_data['time'][periastron_inds[usable_periastrons[rit]]],11))
                bin_ind = 1
                while bin_ind < len(time_bins_1):
                    time_bin_inds_1 = np.where((particle_data['time'] > time_bins_1[bin_ind-1]) & (particle_data['time'] < time_bins_1[bin_ind]))[0]
                    intergrated_values = np.trapz(particle_data['mdot'][:,time_bin_inds_1], particle_data['time'][time_bin_inds_1])/(particle_data['time'][time_bin_inds_1[-1]]-particle_data['time'][time_bin_inds_1[0]])
                    binned_accretion[0].append(intergrated_values[0])
                    binned_accretion[1].append(intergrated_values[1])
                    total_accretion.append(np.sum(intergrated_values))
                    #binned_accretion[0].append(np.median(particle_data['mdot'][0][time_bin_inds]))
                    #binned_accretion[1].append(np.median(particle_data['mdot'][1][time_bin_inds]))
                    #total_accretion.append(np.median(particle_data['mdot'][:,time_bin_inds]))
                    bin_ind = bin_ind + 1
                bin_ind = 1
                while bin_ind < len(time_bins_2):
                    time_bin_inds_2 = np.where((particle_data['time'] > time_bins_2[bin_ind-1]) & (particle_data['time'] < time_bins_2[bin_ind]))[0]
                    intergrated_values = np.trapz(particle_data['mdot'][:,time_bin_inds_2], particle_data['time'][time_bin_inds_2])/(particle_data['time'][time_bin_inds_2[-1]]-particle_data['time'][time_bin_inds_2[0]])
                    binned_accretion[0].append(intergrated_values[0])
                    binned_accretion[1].append(intergrated_values[1])
                    total_accretion.append(np.sum(intergrated_values))
                    bin_ind = bin_ind + 1
                averaged_binned_accretion[0].append(binned_accretion[0])
                averaged_binned_accretion[1].append(binned_accretion[1])
                averaged_total_accretion.append(total_accretion)
                rit = rit + 1
            
            ax_label = 'ax' + str(plot_it)
            if plot_it == 0:
                axes_dict.update({ax_label:fig.add_subplot(gs[0,plot_it])})
                axes_dict[ax_label].tick_params(axis="x",direction="in")
                axes_dict[ax_label].set_xlim([0.0, 1.3])
                axes_dict[ax_label].set_ylim([0.0, 7.0])
                axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                yticklabels = axes_dict[ax_label].get_yticklabels()
                plt.setp(yticklabels[0], visible=False)
                xticklabels = axes_dict[ax_label].get_xticklabels()
                plt.setp(xticklabels[-2], visible=False)
                plt.setp(xticklabels, visible=False)
            else:
                axes_dict.update({ax_label:fig.add_subplot(gs[int(plot_it/columns),np.remainder(plot_it,columns)], sharex=axes_dict['ax0'], sharey=axes_dict['ax0'])})
                if plot_it < 1:
                    xticklabels = axes_dict[ax_label].get_xticklabels()
                    plt.setp(xticklabels, visible=False)
                elif plot_it > 0:
                    if plot_it != len(e_bins)-2:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels[0], visible=False)
                elif plot_it > columns-1:
                    if plot_it != len(e_bins)-2:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels[-2], visible=False)
                if np.remainder(plot_it,columns) != 0:
                    yticklabels = axes_dict[ax_label].get_yticklabels()
                    plt.setp(yticklabels, visible=False)
                    yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
                    plt.setp(yticklabels, visible=False)
                else:
                    axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                axes_dict[ax_label].set_xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
                axes_dict[ax_label].tick_params(axis="x",direction="in")
                
            median_accretion = []
            median_accretion_total = []
            standard_deviation = []
            standard_deviation_total = []
            median_acc = np.median(averaged_binned_accretion[0], axis=0)
            mean_acc = np.mean(averaged_binned_accretion[0], axis=0)
            std_acc = np.std(averaged_binned_accretion[0], axis=0)
            median_accretion.append(median_acc)
            standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
            median_acc = np.median(averaged_binned_accretion[1], axis=0)
            mean_acc = np.mean(averaged_binned_accretion[1], axis=0)
            std_acc = np.std(averaged_binned_accretion[1], axis=0)
            median_accretion.append(median_acc)
            standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
            
            median_acc_tot = np.median(averaged_total_accretion, axis=0)
            mean_acc_tot = np.mean(averaged_total_accretion, axis=0)
            std_acc_tot = np.std(averaged_total_accretion, axis=0)
            median_accretion_total.append(median_acc_tot)
            standard_deviation_total.append([median_acc_tot-(mean_acc_tot-std_acc_tot), (mean_acc_tot+std_acc_tot)-median_acc_tot])
            #plt.plot(phase_2, median_accretion[0].tolist()+median_accretion[0][1:].tolist(), ls='steps', alpha=0.5, label='Primary')
            #plt.plot(phase_2, median_accretion[1].tolist()+median_accretion[1][1:].tolist(), ls='steps', alpha=0.5, label='Secondary')
            #plt.plot(phase_2, median_total.tolist() + median_total[1:].tolist(),ls='steps', label='Total')
            long_median_accretion_1 = median_accretion[0].tolist()+median_accretion[0].tolist()+median_accretion[0].tolist()
            long_median_accretion_2 = median_accretion[1].tolist()+median_accretion[1].tolist()+median_accretion[1].tolist()
            long_median_total = median_accretion_total[0].tolist() + median_accretion_total[0].tolist() + median_accretion_total[0].tolist()
            yerr_1 = [standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist(), standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()]
            yerr_2 = [standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist(), standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()]
            yerr = [standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist(), standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()]
            accretion_err.append(np.array(yerr)/np.array(long_median_total))
            axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_1)*(1.e4), yerr=np.array(yerr_1)*(1.e4), ls='steps-mid', alpha=0.5, label='Primary')
            axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_2)*(1.e4), yerr=np.array(yerr_2)*(1.e4), ls='steps-mid', alpha=0.5, label='Secondary')
            axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_total)*(1.e4), yerr=np.array(yerr)*(1.e4), ls='steps-mid', label='Total')
            time_text = plt.text(0.1, axes_dict[ax_label].get_ylim()[1]-0.3, '$e$=['+str(e_bins[e_bin_it-1])+','+str(e_bins[e_bin_it])+'] using '+str(len(usable_periastrons))+' orbits', va="center", ha="left", color='k', fontsize=args.text_font)
            #popt, pcov = optimize.curve_fit(skew_norm_pdf, phase_centers[22:43], (np.array(long_median_total[22:43])*(1.e4)))
            #x_fit = np.linspace(phase_centers[22], phase_centers[43], 1000)
            #y_fit = skew_norm_pdf(x_fit, *popt)
            #plt.plot(x_fit, y_fit)
            multiple_folds.append(np.array(long_median_total)*(1.e4))
            #Calculate fits
            x_data = phase_centers[23:-15]
            x = np.linspace(np.min(x_data),np.max(x_data),100)
            if os.path.exists(file_name+'.pkl'):
                file_open = open(file_name+'.pkl', 'rb')
                phase_centers, long_median_accretion, yerr_tot, popt = pickle.load(file_open)
                file_open.close()
                y_fit= func(x, *popt)
            else:
                y_data = np.array(long_median_total)[23:-15]
                results = []
                pulsed_likelihoods = []
                no_pulsed_likelihoods = []
                for tries in range(1000):
                    sigma = np.random.random()*2*0.2
                    amp = np.random.random()*np.max(y_data)
                    p = np.array([sigma, x_data[12:20][np.argmax(y_data[12:20])],-5,np.min(y_data),amp])
                    try:
                        popt, pcov = curve_fit(func, x_data, y_data, p)
                        err = np.sum(np.abs(func(x_data, *popt) - y_data))
                        results.append((err, popt))
                        y_fit_data = func(x_data, *popt)
                        chi_squared_pulsed = np.sum((long_median_total[23:-15] - y_fit_data)**2./err**2.)
                        maximum_likelihood_pulsed = np.exp(-chi_squared_pulsed/2.)
                        y_no_pulsed = np.ones(np.shape(y_data))*np.random.normal(np.mean(long_median_total), np.std(long_median_total))
                        chi_squared_no_pulsed = np.sum((long_median_total[23:-15] - y_no_pulsed)**2./err**2.)
                        maximum_likelihood_no_pulsed = np.exp(-chi_squared_no_pulsed/2.)
                        pulsed_likelihoods.append(maximum_likelihood_pulsed)
                        no_pulsed_likelihoods.append(maximum_likelihood_no_pulsed)
                    except:
                        pass
                bayes_factor = np.mean(pulsed_likelihoods)/np.mean(no_pulsed_likelihoods)
                print("bayes_factor="+str(bayes_factor)+" for $e$=["+str(e_bins[e_bin_it-1])+","+str(e_bins[e_bin_it])+"]")
                err, popt = min(results, key=lambda x:x[0])
                y_fit= func(x, *popt)
                #import pdb
                #pdb.set_trace()
            y_fits.append(popt)
            #axes_dict[ax_label].plot(x, y_fit*(1.e4), 'k--')
            
            if e_bin_it == len(e_bins) - 1:
                #axes_dict[ax_label].legend(loc='center right', fontsize=args.text_font)#(loc='center left', bbox_to_anchor=(0.985, 0.5), fontsize=args.text_font)
                axes_dict[ax_label].legend(loc='center left', bbox_to_anchor=(0.985, 0.5), fontsize=args.text_font)
            if plot_it == 2:
                #fig.text(0.5, 0.07, "Orbital Phase ($\phi$)", va='center', ha='center')
                plt.xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
            plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
            plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)
            #call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
            e_bin_it = e_bin_it + 1
            plot_it = plot_it + 1
            
            if args.pickle_dump != 'False':
                long_median_accretion = [np.array(long_median_accretion_1)*(1.e4), np.array(long_median_accretion_2)*(1.e4), np.array(long_median_total)*(1.e4)]
                yerr_tot = [np.array(yerr_1)*(1.e4), np.array(yerr_2)*(1.e4), np.array(yerr)*(1.e4)]
                component_accretion_pickle = file_name + '.pkl'
                file = open(component_accretion_pickle, 'wb')
                pickle.dump((phase_centers, long_median_accretion, yerr_tot, popt),file)
                file.close()
    
        n_lines = len(multiple_folds)
        c_index = np.linspace(0.0, 0.95, n_lines)
        alpha_values = np.linspace(1.0, 0.2, n_lines)
        e_int = 1
        plt.clf()
        for i, shift in zip(c_index, multiple_folds):
            plt.plot(phase_centers, shift, color=plt.cm.gnuplot(i), label='$e$=['+str(e_bins[e_int-1])+','+str(e_bins[e_int])+']')#, ls='steps-mid')#, alpha=alpha_values[e_int])
            e_int = e_int + 1
        plt.legend(loc='best')
        plt.xlabel("Orbital Phase ($\phi$)")
        plt.ylabel("Accretion Rate ($10^{-4}$M$_\odot$/yr)")
        plt.xlim([0.0, 1.3])
        plt.ylim([0.0, 5.0])
        if args.pickle_dump != 'False':
            folded_pickle = path + 'using_e_bins.pkl'
            file = open(folded_pickle, 'wb')
            pickle.dump((multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits),file)
            file.close()
        file_name = save_dir + 'multiple_folds_using_e_bins'
        plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
        plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)
    
    elif os.path.exists(folded_pickle):
        file_open = open(folded_pickle, 'rb')
        multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits = pickle.load(file_open)
        file_open.close()
    else:
        #periastron_inds = periastron_inds[10:]
        rit = 0
        plot_it = 0
        n_folded_orbits = args.n_orbits
        
        plt.clf()
        fig = plt.figure()
        if repeats < 6:
            rows = 1
            fig.set_size_inches(3.25*repeats, 3.25)
            gs = gridspec.GridSpec(1, repeats)
        else:
            rows = int(repeats/5)
            fig.set_size_inches(3.25*5, 3.25*rows)
            gs = gridspec.GridSpec(rows, 5)
        
        gs.update(wspace=0.0, hspace=0.0)
        axes_dict = {}
        multiple_folds = []
        median_eccentricity = []
        std_eccentricity = []
        accretion_err = []
        y_fits = []
        particle_data['time'] = np.array(particle_data['time'])
        particle_data['mdot'] = np.array(particle_data['mdot']).T
        phase = np.linspace(0.0, 1.0, 21)
        phase_centers = (phase[1:] + phase[:-1])/2.
        phase_pre = np.linspace(-1.0, 0.0, 21)
        phase_2 = np.linspace(1.0, 2.0, 21)
        phase_2 = phase_pre.tolist()[:-1] + phase.tolist() + phase_2[1:].tolist()
        phase_centers = (np.array(phase_2[1:]) + np.array(phase_2[:-1]))/2.
        while rit < repeats:
            if args.overlapping_folds == 'False':
                hist_ind = args.starting_orbit + rit*n_folded_orbits
            elif args.skip_iterator >0:
                hist_ind = args.starting_orbit + rit*args.skip_iterator
            else:
                hist_ind = args.starting_orbit + rit#*args.skip_iterator
            ap_ind = hist_ind - 1
            #FIND BINNED DATA
            averaged_binned_accretion = [[],[]]
            averaged_total_accretion = []
            median_e = np.median(particle_data['eccentricity'][periastron_inds[hist_ind-1]:periastron_inds[hist_ind+n_folded_orbits-1]])
            if rit == 0:
                fold_bool = True
            elif (median_eccentricity[-1] - median_e) < 0.04:
                fold_bool = False
                repeats = repeats + 1
            else:
                fold_bool = True
            if fold_bool:
                ax_label = 'ax' + str(plot_it)
                if plot_it == 0:
                    axes_dict.update({ax_label:fig.add_subplot(gs[0,plot_it])})
                    axes_dict[ax_label].tick_params(axis="x",direction="in")
                    axes_dict[ax_label].set_xlim([0.0, 1.3])
                    axes_dict[ax_label].set_ylim(bottom=0.0)
                    axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                    yticklabels = axes_dict[ax_label].get_yticklabels()
                    plt.setp(yticklabels[0], visible=False)
                    xticklabels = axes_dict[ax_label].get_xticklabels()
                    plt.setp(xticklabels, visible=False)
                else:
                    if repeats < 6:
                        axes_dict.update({ax_label:fig.add_subplot(gs[0,plot_it], sharex=axes_dict['ax0'], sharey=axes_dict['ax0'])})
                    else:
                        axes_dict.update({ax_label:fig.add_subplot(gs[int(plot_it/5),np.remainder(plot_it,5)], sharex=axes_dict['ax0'], sharey=axes_dict['ax0'])})
                    if plot_it < 5:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels, visible=False)
                    elif plot_it != 9:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels[-2], visible=False)
                    if np.remainder(plot_it,5) != 0:
                        yticklabels = axes_dict[ax_label].get_yticklabels()
                        plt.setp(yticklabels, visible=False)
                        yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
                        plt.setp(yticklabels, visible=False)
                    else:
                        axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
                    if int(plot_it/5) == rows-1 and np.remainder(plot_it,5) == 2:
                        axes_dict[ax_label].set_xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
                    axes_dict[ax_label].tick_params(axis="x",direction="in")
                plot_it = plot_it + 1
                median_eccentricity.append(median_e)
                std_e = np.std(particle_data['eccentricity'][periastron_inds[hist_ind-1]:periastron_inds[hist_ind+n_folded_orbits-1]])
                mean_e = np.mean(particle_data['eccentricity'][periastron_inds[hist_ind-1]:periastron_inds[hist_ind+n_folded_orbits-1]])
                std_eccentricity.append([median_e-(mean_e-std_e), (mean_e+std_e)-median_e])
                if args.overlapping_folds == 'False':
                    end_ind = args.starting_orbit+n_folded_orbits + rit*n_folded_orbits
                elif args.skip_iterator >0:
                    end_ind = args.starting_orbit+n_folded_orbits + rit*args.skip_iterator
                else:
                    end_ind = args.starting_orbit+n_folded_orbits + rit
                while hist_ind < len(periastron_inds[:end_ind]):
                    binned_accretion = [[],[]]
                    total_accretion = []
                    time_bins_1 = sorted(np.linspace(particle_data['time'][periastron_inds[hist_ind-1]], particle_data['time'][apastron_inds[ap_ind]],11))
                    time_bins_2 = sorted(np.linspace(particle_data['time'][apastron_inds[ap_ind]], particle_data['time'][periastron_inds[hist_ind]],11))
                    bin_ind = 1
                    while bin_ind < len(time_bins_1):
                        time_bin_inds_1 = np.where((particle_data['time'] > time_bins_1[bin_ind-1]) & (particle_data['time'] < time_bins_1[bin_ind]))[0]
                        intergrated_values = np.trapz(particle_data['mdot'][:,time_bin_inds_1], particle_data['time'][time_bin_inds_1])/(particle_data['time'][time_bin_inds_1[-1]]-particle_data['time'][time_bin_inds_1[0]])
                        binned_accretion[0].append(intergrated_values[0])
                        binned_accretion[1].append(intergrated_values[1])
                        total_accretion.append(np.sum(intergrated_values))
                        #binned_accretion[0].append(np.median(particle_data['mdot'][0][time_bin_inds]))
                        #binned_accretion[1].append(np.median(particle_data['mdot'][1][time_bin_inds]))
                        #total_accretion.append(np.median(particle_data['mdot'][:,time_bin_inds]))
                        bin_ind = bin_ind + 1
                    bin_ind = 1
                    while bin_ind < len(time_bins_2):
                        time_bin_inds_2 = np.where((particle_data['time'] > time_bins_2[bin_ind-1]) & (particle_data['time'] < time_bins_2[bin_ind]))[0]
                        intergrated_values = np.trapz(particle_data['mdot'][:,time_bin_inds_2], particle_data['time'][time_bin_inds_2])/(particle_data['time'][time_bin_inds_2[-1]]-particle_data['time'][time_bin_inds_2[0]])
                        binned_accretion[0].append(intergrated_values[0])
                        binned_accretion[1].append(intergrated_values[1])
                        total_accretion.append(np.sum(intergrated_values))
                        bin_ind = bin_ind + 1
                    hist_ind = hist_ind + 1
                    ap_ind = ap_ind + 1
                    averaged_binned_accretion[0].append(binned_accretion[0])
                    averaged_binned_accretion[1].append(binned_accretion[1])
                    averaged_total_accretion.append(total_accretion)
                median_accretion = []
                median_accretion_total = []
                standard_deviation = []
                standard_deviation_total = []
                median_acc = np.median(averaged_binned_accretion[0], axis=0)
                mean_acc = np.mean(averaged_binned_accretion[0], axis=0)
                std_acc = np.std(averaged_binned_accretion[0], axis=0)
                median_accretion.append(median_acc)
                standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
                median_acc = np.median(averaged_binned_accretion[1], axis=0)
                mean_acc = np.mean(averaged_binned_accretion[1], axis=0)
                std_acc = np.std(averaged_binned_accretion[1], axis=0)
                median_accretion.append(median_acc)
                standard_deviation.append([median_acc-(mean_acc-std_acc), (mean_acc+std_acc)-median_acc])
                
                median_acc_tot = np.median(averaged_total_accretion, axis=0)
                mean_acc_tot = np.mean(averaged_total_accretion, axis=0)
                std_acc_tot = np.std(averaged_total_accretion, axis=0)
                median_accretion_total.append(median_acc_tot)
                standard_deviation_total.append([median_acc_tot-(mean_acc_tot-std_acc_tot), (mean_acc_tot+std_acc_tot)-median_acc_tot])
                #plt.plot(phase_2, median_accretion[0].tolist()+median_accretion[0][1:].tolist(), ls='steps', alpha=0.5, label='Primary')
                #plt.plot(phase_2, median_accretion[1].tolist()+median_accretion[1][1:].tolist(), ls='steps', alpha=0.5, label='Secondary')
                #plt.plot(phase_2, median_total.tolist() + median_total[1:].tolist(),ls='steps', label='Total')
                long_median_accretion_1 = median_accretion[0].tolist()+median_accretion[0].tolist()+median_accretion[0].tolist()
                long_median_accretion_2 = median_accretion[1].tolist()+median_accretion[1].tolist()+median_accretion[1].tolist()
                long_median_total = median_accretion_total[0].tolist() + median_accretion_total[0].tolist() + median_accretion_total[0].tolist()
                yerr_1 = [standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist()+standard_deviation[0][0].tolist(), standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()+standard_deviation[0][1].tolist()]
                yerr_2 = [standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist()+standard_deviation[1][0].tolist(), standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()+standard_deviation[1][1].tolist()]
                yerr = [standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist()+standard_deviation_total[0][0].tolist(), standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()+standard_deviation_total[0][1].tolist()]
                accretion_err.append(np.array(yerr)/np.array(long_median_total))
                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_1)*(1.e4), yerr=np.array(yerr_1)*(1.e4), ls='steps-mid', alpha=0.5, label='Primary')
                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_accretion_2)*(1.e4), yerr=np.array(yerr_2)*(1.e4), ls='steps-mid', alpha=0.5, label='Secondary')
                axes_dict[ax_label].errorbar(phase_centers, np.array(long_median_total)*(1.e4), yerr=np.array(yerr)*(1.e4), ls='steps-mid', label='Total')
                time_text = plt.text(0.1, 3.3, '$e$='+str(np.round(median_eccentricity[-1]*100)/100), va="center", ha="left", color='k', fontsize=args.text_font)
                #popt, pcov = optimize.curve_fit(skew_norm_pdf, phase_centers[22:43], (np.array(long_median_total[22:43])*(1.e4)))
                #x_fit = np.linspace(phase_centers[22], phase_centers[43], 1000)
                #y_fit = skew_norm_pdf(x_fit, *popt)
                #plt.plot(x_fit, y_fit)
                multiple_folds.append(np.array(long_median_total)*(1.e4))
                
                #Calculate fits
                file_name = save_dir + 'accretion_median_start_orbit_'+ ("%02d" % (hist_ind-n_folded_orbits)) + '_' + str(n_folded_orbits) + '_folded_orbits_e_'+str(np.round(median_eccentricity[-1]*100)/100)
                x_data = phase_centers[23:-15]
                x = np.linspace(np.min(x_data),np.max(x_data),100)
                if os.path.exists(file_name+'.pkl'):
                    file_open = open(file_name+'.pkl', 'rb')
                    phase_centers, long_median_accretion, yerr_tot, popt = pickle.load(file_open)
                    file_open.close()
                    y_fit= func(x, *popt)
                else:
                    y_data = np.array(long_median_total)[23:-15]
                    results = []
                    for tries in range(50):
                        sigma = np.random.random()*2*0.15
                        amp = np.random.random()*2*np.max(y_data)
                        p = np.array([sigma, x_data[12:20][np.argmax(y_data[12:20])], -5,np.min(y_data),amp])
                        try:
                            popt, pcov = curve_fit(func, x_data, y_data, p)
                            err = np.sum(np.abs(func(x_data, *popt) - y_data))
                            results.append((err, popt))
                            if err < 0.1:
                                break
                        except:
                            pass
                    err, popt = min(results, key=lambda x:x[0])
                    y_fit= func(x, *popt)
                y_fits.append(popt)
                axes_dict[ax_label].plot(x, y_fit*(1.e4), 'k--')
                
                if plot_it == 5:
                    axes_dict[ax_label].legend(loc='best', fontsize=args.text_font)
                plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
                plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)
            #call(['convert', '-antialias', '-quality', '100', '-density', '200', '-resize', '100%', '-flatten', file_name+'.eps', file_name+'.jpg'])
            rit = rit + 1
            if args.pickle_dump != 'False':
                long_median_accretion = [np.array(long_median_accretion_1)*(1.e4), np.array(long_median_accretion_2)*(1.e4), np.array(long_median_total)*(1.e4)]
                yerr_tot = [np.array(yerr_1)*(1.e4), np.array(yerr_2)*(1.e4), np.array(yerr)*(1.e4)]
                component_accretion_pickle = file_name + '.pkl'
                file = open(component_accretion_pickle, 'wb')
                pickle.dump((phase_centers, long_median_accretion, yerr_tot, popt),file)
                file.close()
    #multiple_folds = multiple_folds[::2]
    #mean_eccentricity = mean_eccentricity[::2]
    median_eccentricity = np.round(np.array(median_eccentricity)*100)/100.
    n_lines = len(multiple_folds)
    c_index = np.linspace(0.0, 0.95, n_lines)
    alpha_values = np.linspace(1.0, 0.2, n_lines)
    e_int = 0
    plt.clf()
    for i, shift in zip(c_index, multiple_folds):
        plt.plot(phase_centers, shift, color=plt.cm.gnuplot(i), label='e='+str(median_eccentricity[e_int]), ls='steps-mid', alpha=alpha_values[e_int])
        e_int = e_int + 1
    plt.legend(loc='best')
    plt.xlabel("Orbital Phase ($\phi$)")
    plt.ylabel("Accretion Rate ($10^{-4}$M$_\odot$/yr)")
    plt.xlim([0.0, 1.3])
    plt.ylim(bottom=0.0)
    if args.pickle_dump != 'False':
        file = open(folded_pickle, 'wb')
        pickle.dump((multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits),file)
        file.close()
    file_name = save_dir + 'multiple_folds_over_' + str(n_folded_orbits) + '_orbits_repeats_' + str(repeats) + '_overlap_' + args.overlapping_folds
    plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
    plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)

if args.phasefolded_multi == 'True':
    file_name = save_dir + 'multiple_folds_all'
    plt.clf()
    fig = plt.figure()
    fig.set_size_inches(4.0, 6.0)
    #files = ["Mach_0.1/multiple_folds_over_"+str(args.n_orbits)+"_orbits.pkl", "Mach_0.2/multiple_folds_over_"+str(args.n_orbits)+"_orbits.pkl"]
    #files = ["Mach_0.1/using_e_bins.pkl", "Mach_0.2/using_e_bins.pkl"]
    #files = ["Mach_0.2/Lref_10/using_e_bins.pkl", "Mach_0.2/Lref_11/using_e_bins.pkl"]
    files = ["Mach_0.2/Lref_11/using_e_bins.pkl", "Mach_0.2/Lref_12/using_e_bins.pkl"]
    #plot_eccentricities = [0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    gs = gridspec.GridSpec(2, 1)
    gs.update(hspace=0.0)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1, sharey=ax1)
    
    file_open = open(files[0], 'rb')
    #multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits = pickle.load(file_open)
    stuff = pickle.load(file_open)
    multiple_folds = stuff[0]
    phase_centers = stuff[1]
    median_eccentricity = stuff[2]
    std_eccentricity = stuff[3]
    accretion_err = stuff[4]
    n_lines = stuff[5]
    y_fits = stuff[6]
    #if len(stuff) > 6:
    #    multiple_folds_normalised = stuff[6]
    file_open.close()
    e_bins = [1.1, 0.6, 0.4, 0.2, 0.0]
    #e_bins = [1.1, 0.7, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0]

    x_data = phase_centers[23:-15]
    x = np.linspace(np.min(x_data),np.max(x_data),100)
    
    '''
    plotted_int = '0'
    n_lines = 0
    plot_es = []
    for e in mean_eccentricity:
        print(str(e)[2], plotted_int)
        if str(e)[2] != plotted_int:
            plot_es.append(e)
            n_lines = n_lines + 1
            plotted_int = str(e)[2]
    '''
    
    c_index = np.linspace(0.0, 0.95, n_lines)
    alpha_values = np.linspace(1.0, 0.2, n_lines)
    
    multiple_folds_normalised = []
    e_int = 0
    TWA_3A_e = 0.6280
    DQ_TAU_e = 0.568
    TWA_3A = [1.34735, 0.51477,0.47301,0.39734,0.41151,0.50567,0.45831,0.49353,0.38145,0.46946,0.55446,0.45741,0.67232,0.54981,0.55741,0.587,0.68482,0.84026,1.12274,1.49516]
    TWA_3A = TWA_3A + TWA_3A + TWA_3A
    TWA_3A_Area = simps(TWA_3A, dx=0.05)/(phase_centers[-1] - phase_centers[0])
    TWA_3A = TWA_3A/TWA_3A_Area
    DQ_TAU = [83.01204,46.33979,36.36598,23.6284,30.29718,24.75743,27.10763,25.15253,20.86345,21.43372,19.10391,15.56902,16.20591,20.36333,24.30141,20.54966,23.34257,46.89138,45.99525,73.52196]
    DQ_TAU = DQ_TAU + DQ_TAU + DQ_TAU
    DQ_TAU_Area = simps(DQ_TAU, dx=0.05)/(phase_centers[-1] - phase_centers[0])
    DQ_TAU = DQ_TAU/DQ_TAU_Area
    DQ_TAU_sim = []
    TWA_3A_sim = []
    DQ_TAU_err = []
    TWA_3A_err = []
    DQ_TAU_sim_e = []
    TWA_3A_sim_e = []
    TWA_3A_chi_T1 = []
    DQ_TAU_chi_T1 = []
    plot_bool = True
    
    for shift in multiple_folds:
        mean_folded_accretion = simps(shift, dx=0.05)/(phase_centers[-1] - phase_centers[0])
        print('mean_folded_accretion='+str(mean_folded_accretion))
        chi_squared = np.sum(((TWA_3A[20:46]-(shift/mean_folded_accretion)[20:46]**2)/((shift/mean_folded_accretion)[20:46])))
        TWA_3A_chi_T1.append(chi_squared)
        chi_squared = np.sum(((DQ_TAU[20:46]-(shift/mean_folded_accretion)[20:46]**2)/((shift/mean_folded_accretion)[20:46])))
        DQ_TAU_chi_T1.append(chi_squared)
        
        y = func(x, *y_fits[e_int])
        y_area = simps(y, dx=x[1]-x[0])/(x[-1] - x[0])
        #ax1.plot(x, y/y_area, color=plt.cm.magma(c_index[e_int]), label='e=['+str(e_bins[e_int])+','+str(e_bins[e_int+1]))#, alpha=0.5)
        if plot_bool:
            ax1.plot(phase_centers, shift/mean_folded_accretion, color=plt.cm.magma(c_index[e_int]), label='e=['+str(e_bins[e_int])+','+str(e_bins[e_int+1])+']', ls='steps-mid', linewidth=1)#, alpha=0.5)
            #ax1.plot(phase_centers, shift/mean_folded_accretion, color=plt.cm.magma(c_index[e_int]), label='e=['+str(median_eccentricity[e_int]), ls='steps-mid', linewidth=1)
            plot_bool = True
        else:
            plot_bool = True
        #multiple_folds_normalised.append(y_fits[e_int])
        multiple_folds_normalised.append(shift/mean_folded_accretion)
        if args.plot_best_fit != "False" and e_int == 0:
            TWA_3A_sim.append(shift/mean_folded_accretion)
            TWA_3A_err.append(accretion_err[e_int])
            TWA_3A_sim_e.append(median_eccentricity[e_int])
        elif e_int == np.argmin(np.abs(median_eccentricity - TWA_3A_e)) and args.plot_best_fit == "False":
            TWA_3A_sim.append(shift/mean_folded_accretion)
            TWA_3A_err.append(accretion_err[e_int])
            TWA_3A_sim_e.append(median_eccentricity[e_int])
        if (args.plot_best_fit != "False" and e_int == 0):
            DQ_TAU_sim.append(shift/mean_folded_accretion)
            DQ_TAU_err.append(accretion_err[e_int])
            DQ_TAU_sim_e.append(median_eccentricity[e_int])
        elif e_int == np.argmin(np.abs(median_eccentricity - DQ_TAU_e))and args.plot_best_fit == "False":
            DQ_TAU_sim.append(shift/mean_folded_accretion)
            DQ_TAU_err.append(accretion_err[e_int])
            DQ_TAU_sim_e.append(median_eccentricity[e_int])
        e_int = e_int + 1
    #ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.legend(loc='center left', bbox_to_anchor=(0.985, 0.5))
    #ax1.set_xlabel("Orbital Phase ($\phi$)")
    #ax1.xaxis.set_label_coords(1.16, -0.045)
    ax1.set_ylabel("Normalised Accretion")
    ax1.set_xlim([0.0, 1.3])
    #ax1.plot(phase_centers, TWA_3A, color='green', label='TWA 3A', ls='steps-mid')
    #ax1.plot(phase_centers, DQ_TAU, color='blue', label='DQ Tau', ls='steps-mid')
    xticklabels = ax1.get_xticklabels()
    plt.setp(xticklabels, visible=False)
    ax1.tick_params(axis='x', which='major', direction="in")
    print("TWA 3A e=" + str(median_eccentricity[np.argsort(TWA_3A_chi_T1)]))
    print("DQ Tau e=" + str(median_eccentricity[np.argsort(DQ_TAU_chi_T1)]))
    
    '''
    file = open(component_accretion_pickle, 'wb')
    pickle.dump((phase_centers, long_median_accretion, yerr),file)
    file.close()
    '''
    file_open = open(files[0], 'wb')
    pickle.dump((multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits, multiple_folds_normalised),file_open)
    file_open.close()
    '''
    file_open = open(files[1], 'rb')
    try:
        multiple_folds, phase_centers, mean_eccentricity, std_eccentricity, accretion_err, n_lines = pickle.load(file_open)
    except:
        multiple_folds, phase_centers, mean_eccentricity, std_eccentricity, accretion_err, n_lines, multiple_folds_normalised = pickle.load(file_open)
    file_open.close()
    '''
    file_open = open(files[1], 'rb')
    #multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits = pickle.load(file_open)
    stuff = pickle.load(file_open)
    multiple_folds = stuff[0]
    phase_centers = stuff[1]
    median_eccentricity = stuff[2]
    std_eccentricity = stuff[3]
    accretion_err = stuff[4]
    n_lines = stuff[5]
    y_fits = stuff[6]
    #if len(stuff) > 6:
    #    multiple_folds_normalised = stuff[6]
    file_open.close()
    '''
    plotted_int = '0'
    n_lines = 0
    plot_es = []
    for e in mean_eccentricity:
        print(str(e)[2], plotted_int)
        if str(e)[2] != plotted_int:
            plot_es.append(e)
            n_lines = n_lines + 1
            plotted_int = str(e)[2]
    '''
    #e_bins = [1.1, 0.7, 0.5, 0.3, 0.1, 0.0]
    e_bins = [1.1, 0.6, 0.4, 0.2, 0.0]
    c_index = np.linspace(0.0, 0.95, n_lines)

    multiple_folds_normalised = []
    e_int = 0
    TWA_3A_chi_T2 = []
    DQ_TAU_chi_T2 = []
    colour_int = 0
    plot_bool = True
    
    for shift in multiple_folds:
        mean_folded_accretion = simps(shift, dx=0.05)/(phase_centers[-1] - phase_centers[0])
        print('mean_folded_accretion='+str(mean_folded_accretion))
        inds = np.where((shift/mean_folded_accretion)[20:46]>0.0)[0]
        chi_squared = np.sum((TWA_3A[20:46][inds]-(shift/mean_folded_accretion)[20:46][inds]**2)/((shift/mean_folded_accretion)[20:46][inds]))
        TWA_3A_chi_T2.append(chi_squared)
        chi_squared = np.sum((DQ_TAU[20:46][inds]-(shift/mean_folded_accretion)[20:46][inds]**2)/((shift/mean_folded_accretion)[20:46][inds]))
        DQ_TAU_chi_T2.append(chi_squared)

        y = func(x, *y_fits[e_int])
        y_area = simps(y, dx=x[1]-x[0])/(x[-1] - x[0])
        #ax2.plot(x, y/y_area, color=plt.cm.magma(c_index[e_int]), label='e=['+str(e_bins[e_int])+','+str(e_bins[e_int+1]))#, ls='steps-mid')#, alpha=0.5)
        if plot_bool:
            ax2.plot(phase_centers, shift/mean_folded_accretion, color=plt.cm.magma(c_index[e_int]), label='e=['+str(e_bins[e_int])+','+str(e_bins[e_int+1])+']', ls='steps-mid', linewidth=1)#, alpha=0.5)
            #ax2.plot(phase_centers, shift/mean_folded_accretion, color=plt.cm.magma(c_index[e_int]), label='e=['+str(median_eccentricity[e_int]), ls='steps-mid', linewidth=1)#, alpha=0.5)
            plot_bool = True
        else:
            plot_bool = True
        #multiple_folds_normalised.append(y_fits[e_int])
        multiple_folds_normalised.append(shift/mean_folded_accretion)
        if args.plot_best_fit != "False" and e_int == 0:
            TWA_3A_sim.append(shift/mean_folded_accretion)
            TWA_3A_err.append(accretion_err[e_int])
            TWA_3A_sim_e.append(median_eccentricity[e_int])
        elif e_int == np.argmin(np.abs(median_eccentricity - TWA_3A_e)) and args.plot_best_fit == "False":
            TWA_3A_sim.append(shift/mean_folded_accretion)
            TWA_3A_err.append(accretion_err[e_int])
            TWA_3A_sim_e.append(median_eccentricity[e_int])
        if (args.plot_best_fit != "False" and e_int == 0):
            DQ_TAU_sim.append(shift/mean_folded_accretion)
            DQ_TAU_err.append(accretion_err[e_int])
            DQ_TAU_sim_e.append(median_eccentricity[e_int])
        elif e_int == np.argmin(np.abs(median_eccentricity - DQ_TAU_e))and args.plot_best_fit == "False":
            DQ_TAU_sim.append(shift/mean_folded_accretion)
            DQ_TAU_err.append(accretion_err[e_int])
            DQ_TAU_sim_e.append(median_eccentricity[e_int])
        e_int = e_int + 1
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_xlabel("Orbital Phase ($\phi$)")
    ax2.set_ylabel("Normalised Accretion")
    ax2.text(0.1, ax1.get_ylim()[1]*0.9, 'T2', va="center", ha="left", color='k', fontsize=args.text_font)
    ax1.text(0.1, ax1.get_ylim()[1]*0.9, 'T1', va="center", ha="left", color='k', fontsize=args.text_font)
    ax2.set_ylim(bottom=0.0)
    #ax2.plot(phase_centers, TWA_3A, color='green', ls='steps-mid')
    #ax2.plot(phase_centers, DQ_TAU, color='blue', ls='steps-mid')
    ax2.tick_params(axis='x', which='major', direction="in")
    yticklabels = ax2.get_yticklabels()
    plt.setp(yticklabels[-2], visible=False)

    print("TWA 3A e=" + str(median_eccentricity[np.argsort(TWA_3A_chi_T2)]))
    print("DQ Tau e=" + str(median_eccentricity[np.argsort(DQ_TAU_chi_T2)]))
    
    file_open = open(files[1], 'wb')
    pickle.dump((multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits, multiple_folds_normalised),file_open)
    file_open.close()

    plt.savefig(file_name +'_'+str(args.n_orbits)+'.eps', bbox_inches='tight', pad_inches = 0.02)
    plt.savefig(file_name +'_'+str(args.n_orbits)+'.pdf', bbox_inches='tight', pad_inches = 0.02)

    if args.plot_best_fit == "False":
        plt.clf()
        fig = plt.figure()
        fig.set_size_inches(5.0, 4.0)
        
        Munoz_area = simps(Munoz, dx=phase_centers[19:47][1]-phase_centers[19:47][0])/(phase_centers[19:47][-1] - phase_centers[19:47][0])
        plt.errorbar(phase_centers, TWA_3A_sim[0], label='T1 (e=[0.7, 0.5])', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=TWA_3A_err[0]*TWA_3A_sim[0])
        plt.errorbar(phase_centers, TWA_3A_sim[1], label='T2 (e=[0.7, 0.5])', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=TWA_3A_err[1]*TWA_3A_sim[1])
        plt.errorbar(phase_centers[19:47], Munoz/Munoz_area, label=r'Mu$\tilde{\mathrm{n}}$oz & Lai (2016)', drawstyle='steps-mid', ls='-.', color='black', alpha=0.5, linewidth=1)
        plt.errorbar(phase_centers, TWA_3A, label='TWA 3A', linewidth=1, drawstyle='steps-mid')
        plt.errorbar(phase_centers, DQ_TAU, label='DQ Tau', linewidth=1, drawstyle='steps-mid', color='r')
        plt.legend(loc='upper left')
        plt.xlim([0.0, 1.3])
        plt.ylim([0.0, np.max(TWA_3A_sim[1])+1])
        plt.ylabel("Normalised Accretion")
        plt.xlabel("Orbital Phase ($\phi$)")
        plt.savefig('obs_comparason_'+str(args.n_orbits)+'.eps', bbox_inches='tight', pad_inches = 0.02)
        plt.savefig('obs_comparason_'+str(args.n_orbits)+'.pdf', bbox_inches='tight', pad_inches = 0.02)
        '''
        plt.clf()
        fig = plt.figure()
        fig.set_size_inches(4.0, 6.0)
        gs = gridspec.GridSpec(2, 1)
        gs.update(hspace=0.0)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[1,0], sharex=ax1, sharey=ax1)
        ax1.errorbar(phase_centers, TWA_3A_sim[0], label='T1 (e='+str(TWA_3A_sim_e[0])+')', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=TWA_3A_err[0]*TWA_3A_sim[0])
        ax1.errorbar(phase_centers, TWA_3A_sim[1], label='T2 (e='+str(TWA_3A_sim_e[1])+')', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=TWA_3A_err[1]*TWA_3A_sim[1])
        ax1.errorbar(phase_centers, TWA_3A, label='Observed', drawstyle='steps-mid')
        ax1.legend(loc='upper left')
        ax1.set_ylabel("Normalised Accretion")
        ax1.set_xlim([0.0, 1.3])
        xticklabels = ax1.get_xticklabels()
        plt.setp(xticklabels, visible=False)
        ax1.tick_params(axis='x', which='major', direction="in")
        ax2.errorbar(phase_centers, DQ_TAU_sim[0], label='T1 (e='+str(DQ_TAU_sim_e[0])+')', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=DQ_TAU_err[0]*DQ_TAU_sim[0])
        ax2.errorbar(phase_centers, DQ_TAU_sim[1], label='T2 (e='+str(DQ_TAU_sim_e[1])+')', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=DQ_TAU_err[1]*DQ_TAU_sim[1])
        ax2.errorbar(phase_centers, DQ_TAU, label='Observed', drawstyle='steps-mid', color='r')
        ax2.legend(loc='upper left')
        ax2.set_xlabel("Orbital Phase ($\phi$)")
        ax2.set_ylabel("Normalised Accretion")
        ax2.text(1.01, ax1.get_ylim()[1]-0.2, 'DQ Tau', va="center", ha="left", color='k', fontsize=args.text_font)
        ax1.text(1.01, ax1.get_ylim()[1]-0.2, 'TWA 3A', va="center", ha="left", color='k', fontsize=args.text_font)
        ax1.set_ylim(bottom=0.0)
        yticklabels = ax2.get_yticklabels()
        plt.setp(yticklabels[-1], visible=False)
        plt.savefig('obs_comparason_'+str(args.n_orbits)+'.eps', bbox_inches='tight', pad_inches = 0.02)
        plt.savefig('obs_comparason_'+str(args.n_orbits)+'.pdf', bbox_inches='tight', pad_inches = 0.02)
        #ax1.set_ylim([0.0, 2.5])
        '''
    else:
        plt.clf()
        fig = plt.figure()
        fig.set_size_inches(5.0, 4.0)
        plt.errorbar(phase_centers, TWA_3A_sim[0], label='T1 (e='+str(TWA_3A_sim_e[0])+')', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=TWA_3A_err[0]*TWA_3A_sim[0], linewidth=1)
        plt.errorbar(phase_centers, TWA_3A_sim[1], label='T2 (e='+str(TWA_3A_sim_e[1])+')', drawstyle='steps-mid', alpha=0.5, ls='--', yerr=TWA_3A_err[1]*TWA_3A_sim[1], linewidth=1)
        plt.errorbar(phase_centers, TWA_3A, label='TWA 3A', drawstyle='steps-mid')
        plt.errorbar(phase_centers, DQ_TAU, label='DQ TAU', drawstyle='steps-mid')
        plt.legend(loc='upper left')
        plt.xlim([0.0, 1.3])
        plt.ylim(bottom=0.0)
        plt.ylabel("Normalised Accretion")
        plt.xlabel("Orbital Phase ($\phi$)")
        plt.savefig('best_fit_'+str(args.n_orbits)+'.eps', bbox_inches='tight', pad_inches = 0.02)
        plt.savefig('best_fit_'+str(args.n_orbits)+'.pdf', bbox_inches='tight', pad_inches = 0.02)

if args.plot_beta == "True":
    #files = ["Mach_0.1/multiple_folds_over_"+str(args.n_orbits)+"_orbits.pkl", "Mach_0.2/multiple_folds_over_"+str(args.n_orbits)+"_orbits.pkl"]
    #files = ["Mach_0.1/using_e_bins.pkl", "Mach_0.2/using_e_bins.pkl"]
    files = ["Mach_0.2/Lref_09/using_e_bins.pkl", "Mach_0.2/Lref_10/using_e_bins.pkl", "Mach_0.2/Lref_11/using_e_bins.pkl", "Mach_0.2/Lref_12/using_e_bins.pkl"]
    file_name = save_dir + 'beta_vs_e_'+str(args.n_orbits)
    top_bins = 3
    use_accretion_err = True
    if args.beta_method == 1:
        file_name = file_name + '_mean'
    elif args.beta_method == 2:
        file_name = file_name + '_median'
    elif args.beta_method == 3:
        file_name = file_name + '_max_bin'
    elif args.beta_method == 4:
        file_name = file_name + '_median_top_'+str(top_bins)
    elif args.beta_method == 5:
        file_name = file_name + '_mean_top_'+str(top_bins)
    elif args.beta_method == 6:
        file_name = file_name + '_gaussian_fit'
    #file_name = 'beta_takes_mean_of_accretion_err_bins'
    plt.clf()
    fig = plt.figure()
    label_it = 0
    periastron_inds = [0.8, 1.1]
    quiescent_ind = [0.2, 0.75]
    #markers = ['o', '^']
    #labels = ['T1', 'T2']
    markers = ['o', '^', 's', '+']
    labels = ['$L_\mathrm{ref}=11$', '$L_\mathrm{ref}=12$', '$L_\mathrm{ref}=13$', '$L_\mathrm{ref}=14$']
    for file in files:
        file_open = open(file, 'rb')
        multiple_folds, phase_centers, median_eccentricity, std_eccentricity, accretion_err, n_lines, y_fits, multiple_folds_normalised = pickle.load(file_open)
        file_open.close()
        
        accretion_err = np.nan_to_num(accretion_err)
        beta_total = []
        beta_total_err = []
        max_acc = []
        max_err_arr = []
        base_acc = []
        base_err_arr = []
        max_ind_1 = np.where((phase_centers - 0.025)==periastron_inds[0])[0][0]
        max_ind_2 = np.where(abs(np.around(((phase_centers - 0.025)-periastron_inds[1]), decimals=2))==0)[0][0]
        quiescent_ind_1 = np.where((phase_centers - 0.025)==quiescent_ind[0])[0][0]
        quiescent_ind_2 = np.where((phase_centers - 0.025)==quiescent_ind[1])[0][0]
        confidence_int = 1
        for orbit in range(len(multiple_folds_normalised)):
            if args.beta_method == 1:
                usable_inds_max = np.argsort(multiple_folds_normalised[orbit][max_ind_1:max_ind_2])[::-1]
                max_bin_values = multiple_folds_normalised[orbit][max_ind_1:max_ind_2][usable_inds_max]
                max_mean = np.mean(max_bin_values)
                max_err = [np.sqrt(np.sum(np.square(accretion_err[orbit][0][max_ind_1:max_ind_2][usable_inds_max]*max_bin_values)))/len(max_bin_values)/max_mean, np.sqrt(np.sum(np.square(accretion_err[orbit][1][max_ind_1:max_ind_2][usable_inds_max]*max_bin_values)))/len(max_bin_values)/max_mean]
                
                base_bin_values = multiple_folds_normalised[orbit][quiescent_ind_1:quiescent_ind_2]
                base_mean = np.mean(base_bin_values)
                #base_err = [np.mean(accretion_err[orbit][0][quiescent_ind_1:quiescent_ind_2]*base_bin_values)/base_median, np.mean(accretion_err[orbit][1][quiescent_ind_1:quiescent_ind_2]*base_bin_values)/base_median]
                base_err = [np.sqrt(np.sum(np.square(accretion_err[orbit][0][quiescent_ind_1:quiescent_ind_2]*base_bin_values)))/len(base_bin_values)/base_mean, np.sqrt(np.sum(np.square(accretion_err[orbit][1][quiescent_ind_1:quiescent_ind_2]*base_bin_values)))/len(base_bin_values)/base_mean]
                
                beta = max_mean/base_mean
                beta_err_rel = np.sqrt(np.array(base_err)**2. + np.array(max_err)**2.)
                beta_err = beta*beta_err_rel
                beta_total.append(beta)
                beta_total_err.append(beta_err)
                print("For e =", median_eccentricity[orbit], ", beta =", beta, "+/-", beta_err)
            
            elif args.beta_method == 2:
                usable_inds_max = np.argsort(multiple_folds_normalised[orbit][max_ind_1:max_ind_2])[::-1]
                max_bin_values = multiple_folds_normalised[orbit][max_ind_1:max_ind_2][usable_inds_max]
                max_median = np.median(max_bin_values)
                max_err = [np.sqrt(np.sum(np.square(accretion_err[orbit][0][max_ind_1:max_ind_2][usable_inds_max]*max_bin_values)))/len(max_bin_values)/max_median, np.sqrt(np.sum(np.square(accretion_err[orbit][1][max_ind_1:max_ind_2][usable_inds_max]*max_bin_values)))/len(max_bin_values)/max_median]
                
                base_bin_values = multiple_folds_normalised[orbit][quiescent_ind_1:quiescent_ind_2]
                base_median = np.median(base_bin_values)
                #base_err = [np.mean(accretion_err[orbit][0][quiescent_ind_1:quiescent_ind_2]*base_bin_values)/base_median, np.mean(accretion_err[orbit][1][quiescent_ind_1:quiescent_ind_2]*base_bin_values)/base_median]
                base_err = [np.sqrt(np.sum(np.square(accretion_err[orbit][0][quiescent_ind_1:quiescent_ind_2]*base_bin_values)))/len(base_bin_values)/base_median, np.sqrt(np.sum(np.square(accretion_err[orbit][1][quiescent_ind_1:quiescent_ind_2]*base_bin_values)))/len(base_bin_values)/base_median]
                
                beta = max_median/base_median
                beta_err_rel = np.sqrt(np.array(base_err)**2. + np.array(max_err)**2.)
                beta_err = beta*beta_err_rel
                beta_total.append(beta)
                beta_total_err.append(beta_err)
                print("For e =", median_eccentricity[orbit], ", beta =", beta, "+/-", beta_err)
            
            elif args.beta_method == 4:
                usable_inds_max = np.argsort(multiple_folds_normalised[orbit][max_ind_1:max_ind_2])[::-1][:top_bins]
                max_bin_values = multiple_folds_normalised[orbit][max_ind_1:max_ind_2][usable_inds_max]
                max_median = np.median(max_bin_values)
                max_err = [np.sqrt(np.sum(np.square(accretion_err[orbit][0][max_ind_1:max_ind_2][usable_inds_max]*max_bin_values)))/len(max_bin_values)/max_median, np.sqrt(np.sum(np.square(accretion_err[orbit][1][max_ind_1:max_ind_2][usable_inds_max]*max_bin_values)))/len(max_bin_values)/max_median]
                
                
                base_bin_values = multiple_folds_normalised[orbit][quiescent_ind_1:quiescent_ind_2]
                base_median = np.median(base_bin_values)
                #base_err = [np.mean(accretion_err[orbit][0][quiescent_ind_1:quiescent_ind_2]*base_bin_values)/base_median, np.mean(accretion_err[orbit][1][quiescent_ind_1:quiescent_ind_2]*base_bin_values)/base_median]
                base_err = [np.sqrt(np.sum(np.square(accretion_err[orbit][0][quiescent_ind_1:quiescent_ind_2]*base_bin_values)))/len(base_bin_values)/base_median, np.sqrt(np.sum(np.square(accretion_err[orbit][1][quiescent_ind_1:quiescent_ind_2]*base_bin_values)))/len(base_bin_values)/base_median]
                
                beta = max_median/base_median
                beta_err_rel = np.sqrt(np.array(base_err)**2. + np.array(max_err)**2.)
                beta_err = beta*beta_err_rel
                beta_total.append(beta)
                beta_total_err.append(beta_err)
                print("For e =", median_eccentricity[orbit], ", beta =", beta, "+/-", beta_err)
            
        
        plt.errorbar(median_eccentricity, beta_total, xerr=np.array(std_eccentricity).T, label=labels[label_it], fmt=markers[label_it], yerr=np.array(beta_total_err).T)
        label_it = label_it + 1
        print("beta="+str(beta_total))

    plt.xlabel('eccentricity')
    
    plt.legend(loc='best')
    plt.axhline(y=1.0, ls='--', color='k')
    #plt.axvline(x=3.5, ls='--')
    plt.xlim(left=0)
    plt.ylim([0.0, 55])
    plt.ylabel('$\\beta$')
    ymax = np.max(beta_total) + 1
    plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
    plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)

if args.resolution_study == 'True':
    file_name = 'resolution_study'
    e_bins = [1.1, 0.6, 0.4, 0.2, 0.0]
    dirs = ['Mach_0.2/Lref_09/', 'Mach_0.2/Lref_10/', 'Mach_0.2/Lref_11/', 'Mach_0.2/Lref_12/']
    
    plt.clf()
    fig = plt.figure()
    columns = 4
    rows = 2#int((len(e_bins)-1)/4)
    fig.set_size_inches(3.25*columns, 3.25*rows)
    gs = gridspec.GridSpec(rows, columns)
    
    gs.update(wspace=0.0, hspace=0.0)
    axes_dict = {}
    
    e_bin_it = 1
    plot_it=0
    
    while e_bin_it < len(e_bins):
        ax_label = 'ax' + str(plot_it)
        if plot_it == 0:
            axes_dict.update({ax_label:fig.add_subplot(gs[0,plot_it])})
            axes_dict[ax_label].tick_params(axis="x",direction="in")
            axes_dict[ax_label].set_xlim([0.0, 1.3])
            axes_dict[ax_label].set_ylim([0.0, 7.0])
            axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
            yticklabels = axes_dict[ax_label].get_yticklabels()
            #plt.setp(yticklabels[0], visible=False)
            xticklabels = axes_dict[ax_label].get_xticklabels()
            #plt.setp(xticklabels[-2], visible=False)
            #plt.setp(xticklabels, visible=False)
            axes_dict[ax_label].set_xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
        else:
            axes_dict.update({ax_label:fig.add_subplot(gs[int(plot_it/columns),np.remainder(plot_it,columns)], sharex=axes_dict['ax0'], sharey=axes_dict['ax0'])})
            if plot_it < 1:
                xticklabels = axes_dict[ax_label].get_xticklabels()
                plt.setp(xticklabels, visible=False)
            elif plot_it > 0:
                if plot_it != len(e_bins)-2:
                    xticklabels = axes_dict[ax_label].get_xticklabels()
                    plt.setp(xticklabels[0], visible=False)
            elif plot_it > columns-1:
                if plot_it != len(e_bins)-2:
                    xticklabels = axes_dict[ax_label].get_xticklabels()
                    plt.setp(xticklabels[-2], visible=False)
            if np.remainder(plot_it,columns) != 0:
                yticklabels = axes_dict[ax_label].get_yticklabels()
                plt.setp(yticklabels, visible=False)
                yticklabels = axes_dict[ax_label].get_yticklabels(minor=True)
                plt.setp(yticklabels, visible=False)
            else:
                axes_dict[ax_label].set_ylabel("Accretion Rate ($10^{-4}$ M$_\odot$/yr)", fontsize=args.text_font)
            axes_dict[ax_label].set_xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
            axes_dict[ax_label].tick_params(axis="x",direction="in")
            
        refinement_label = ["$L_\mathrm{ref}$ = 11", "$L_\mathrm{ref}$ = 12", "$L_\mathrm{ref}$ = 13", "$L_\mathrm{ref}$ = 14"]
        linestyles = ['steps-mid:', 'steps-mid-.', 'steps-mid--', 'steps-mid']
        color = ['b', 'r', 'g', 'c']
        for dir_it in range(len(dirs)):
            pickle_file = dirs[dir_it] + "accretion_median_start_orbit_from_" + str(e_bins[e_bin_it-1]) + "_" + str(e_bins[e_bin_it]) + ".pkl"
            
            if os.path.exists(pickle_file):
                file = open(pickle_file, 'rb')
                phase_centers, long_median_accretion, yerr_tot, popt = pickle.load(file)
                file.close()
                axes_dict[ax_label].errorbar(phase_centers, long_median_accretion[2], yerr=yerr_tot[2], ls=linestyles[dir_it], label=refinement_label[dir_it], color=color[dir_it])
            
        time_text = plt.text(0.1, axes_dict[ax_label].get_ylim()[1]-0.3, '$e$=['+str(e_bins[e_bin_it-1])+','+str(e_bins[e_bin_it])+']', va="center", ha="left", color='k', fontsize=args.text_font)
        if e_bin_it == 1:
        #axes_dict[ax_label].legend(loc='center right', fontsize=args.text_font)#(loc='center left', bbox_to_anchor=(0.985, 0.5), fontsize=args.text_font)
            axes_dict[ax_label].legend(loc='center left', fontsize=args.text_font)
        if plot_it == 2:
            #fig.text(0.5, 0.07, "Orbital Phase ($\phi$)", va='center', ha='center')
            plt.xlabel("Orbital Phase ($\phi$)", fontsize=args.text_font)
        
        plot_it = plot_it + 1
        e_bin_it = e_bin_it + 1
        plt.savefig(file_name +'.eps', bbox_inches='tight', pad_inches = 0.02)
        plt.savefig(file_name +'.pdf', bbox_inches='tight', pad_inches = 0.02)
    
            
    
