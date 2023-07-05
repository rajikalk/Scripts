import numpy as np
import glob
import os
import pickle
#import matplotlib.patches
import collections
import matplotlib
#import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_ramses_module as mym
import gc
import yt
import sys
from mpi4py.MPI import COMM_WORLD as CW
import csv
from pyramses import rsink
#from mpi4py.MPI import COMM_WORLD as CW

rank = CW.Get_rank()
size = CW.Get_size()

#matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
#matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
#matplotlib.rcParams['mathtext.rm'] = 'Arial'
#matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
#matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
#matplotlib.rcParams['mathtext.rm'] = 'Arial'
#matplotlib.rcParams['mathtext.sf'] = 'Arial'
#matplotlib.rcParams['mathtext.default'] = 'regular'
#matplotlib.rcParams['font.sans-serif'] = 'Arial'
#matplotlib.rcParams['font.family'] = 'sans-serif'
#matplotlib.rcParams['text.latex.preamble'] = [
#       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
#       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
#       r'\usepackage{helvet}',    # set the normal font here
#       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
#       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
#]


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifetime threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("-ts", "--timescale", help="Do you want to plot in terms of the of the actual simulation time, or the the time since the formation of the first sink?", type=str, default="first_sink")
    parser.add_argument("-fig_suffix", "--figure_suffix", help="Do you want add a suffix to the figures?", type=str, default="")
    parser.add_argument("-add_hist", "--add_histograms", help="Do you want to add the histograms at the end of the super plots?", type=str, default="False")
    parser.add_argument("-all_sep_evol", "--plot_all_separation_evolution", help="do you want to plot all separation evolution?", type=str, default="True")
    parser.add_argument("-x_field", "--x_field", help="Default for x-axis in the multiplot is time", type=str, default="Time")
    parser.add_argument("-tf", "--text_font", help="what font do you want the text to have?", type=int, default=10)
    parser.add_argument("-plt_key", "--plot_key", help="What dictionary key from superplot_dict do you want to plot?", type=str, default='System_seps')
    parser.add_argument("-smooth", "--smooth_bool", help="Do you want to smooth what you are plotting?", type=str, default='False')
    parser.add_argument("-smooth_window", "--smooth_window_val", help="What big (in yrs) do you want the smoothing window to be?", type=float, default=1000)
    parser.add_argument("-sim_id", "--simulation_id", type=str, default='G50')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

#rank = CW.Get_rank()
#size = CW.Get_size()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
plot_booleans = [[False, True], [False, False], [True, True], [True, False]]

args = parse_inputs()

sim_label = args.simulation_id

#Server pickles
pickle_files = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/"+sim_label+"/Max_iter_100/means_superplot.pkl"]

birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/"+sim_label+"/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl"] #"/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl"]

low_cadence_birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/"+sim_label+"/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl"]
#Local pickles

if sim_label != 'G100':
    files = sorted(glob.glob("/lustre/astro/troels/IMF_"+sim_label+"/data/*/info*.txt"))
    sink_files = sorted(glob.glob("/lustre/astro/troels/IMF_"+sim_label+"/data/output*/*.dat"))
else:
    files = sorted(glob.glob("/lustre/astro/troels/IMF_256_fixed_dt/data/*/info*.txt"))
    sink_files = sorted(glob.glob("/lustre/astro/troels/IMF_256_fixed_dt/data/output*/*.dat"))
global_pickle = '/groups/astro/rlk/Analysis_plots/Ramses/Global/'+sim_label+'/stars_imf_'+sim_label+'.pkl'

plt.clf()
'''
if plot_booleans[rank][0] == True:
    fig, axs = plt.subplots(len(pickle_files), 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, len(pickle_files)*3), sharex='col')
    iter_range = range(0, len(pickle_files)*2, 2)
    args.figure_suffix = args.figure_suffix + '_hist'
else:
'''

'''
#fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(two_col_width,page_height))
iter_range = range(0, len(pickle_files)*2)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.07)

#if plot_booleans[rank][1] == True:
#    args.figure_suffix = args.figure_suffix + '_thres'

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
GMC_mass_arr = [12000]
Shrinkage_hist = []
Shrinkage_bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
Shrinks = []
Lifetime_bins = [0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
Lifetime_hist = []
S_bins = np.logspace(1,4,13)
CF_hist = np.zeros((len(pickle_files),12)).tolist()
Times_all = []
M_tot_all = []
M_tot_multi_all = []
N_stars_all = []
N_multi_stars_all = []
SFE_all = []
t_max = np.nan
G50_t_max = 0
'''

#sink_files = sorted(glob.glob("/lustre/astro/troels/IMF_G400/data/output*/*.dat"))
#files = sorted(glob.glob("/lustre/astro/troels/IMF_G400/data/*/info*.txt"))

rm_files = []
for info_name in files:
    sink_file = info_name.split('info')[0]+'stars_output.dat'
    if sink_file not in sink_files:
        rm_files.append(info_name)
for rm_file in rm_files:
    files.remove(rm_file)
del sink_files
gc.collect()


#Define units to override:
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = global_pickle.split('/G')[-1].split('/')[0] #save_dir.split('G')[-1].split('/')[0]

if simulation_density_id == '50':
    Grho=50
    units_override.update({"mass_unit":(1500,"Msun")})
elif simulation_density_id == '100':
    Grho=100
    units_override.update({"mass_unit":(3000,"Msun")})
elif simulation_density_id == '125':
    Grho=125
    units_override.update({"mass_unit":(3750,"Msun")})
elif simulation_density_id == '150':
    Grho=150
    units_override.update({"mass_unit":(4500,"Msun")})
elif simulation_density_id == '200':
    Grho=200
    units_override.update({"mass_unit":(6000,"Msun")})
elif simulation_density_id == '400':
    Grho=400
    units_override.update({"mass_unit":(12000,"Msun")})
else:
    print("MASS UNIT NOT SET")
    import pdb
    pdb.set_trace()
    
del simulation_density_id
gc.collect()

units_override.update({"density_unit":(units_override['mass_unit'][0]/units_override['length_unit'][0]**3, "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
    
file_times = []
if rank == 0:
    for file in files:
        with open(file, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if row[0].split(' ')[0] == 'time':
                    file_time = eval(row[0].split(' ')[-1])
                    file_times.append(file_time)
                    
        f.close()
    file_times = np.array(file_times)*units['time_unit'].in_units('yr')
file_times = CW.bcast(file_times, root=0)

#print('Calculated the units on rank', rank)
sys.stdout.flush()
CW.Barrier()

#================================================================================================
#plot fractions of core fragmentation and dynamical capture

for pick_it in range(len(pickle_files)):
    #try:
    file = open(pickle_files[pick_it], 'rb')
    superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()
    
    file = open(birth_con_pickles[pick_it], 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()
    
    for sys_key in superplot_dict['System_times'].keys():
        if Lifetimes_sys[sys_key] > 1000:
            if len(flatten(eval(sys_key))) == 3 and str(np.max(flatten(eval(sys_key)))) in Sink_birth_all.keys():
                youngest_birth_con = Sink_birth_all[str(np.max(flatten(eval(sys_key))))]
                if superplot_dict['System_seps'][sys_key][0][0] < 1000 and superplot_dict['System_seps'][sys_key][0][0] > 100 and superplot_dict['System_seps'][sys_key][0][1] > 300 and superplot_dict['System_seps'][sys_key][0][1] < 3000:
                    if youngest_birth_con[0] == True or youngest_birth_con[1]==youngest_birth_con[2]:
                        print("Found Triple candidate:", sys_key)
                        #make frame!
                        #Find which global frame I need
                        if youngest_birth_con[0] == True:
                            goal_time = superplot_dict['System_times'][sys_key][0]
                            closest_time_ind = np.argmin(abs(file_times.value - goal_time))
                            if file_times[closest_time_ind] < goal_time:
                                pre_form_ind = closest_time_ind
                                post_form_ind = closest_time_ind + 1
                            else:
                                post_form_ind = closest_time_ind
                                if closest_time_ind != 0:
                                    pre_form_ind = closest_time_ind - 1
                                else:
                                    pre_form_ind = closest_time_ind
                            usable_files = [files[post_form_ind], files[pre_form_ind]]
                        else:
                            form_time = youngest_birth_con[-1]*units['time_unit'].in_units('yr')
                            sys_form_time = superplot_dict['System_times'][sys_key][0]
                            closest_form_ind = np.argmin(abs(file_times.value - form_time.value))
                            if file_times[closest_form_ind] < form_time:
                                pre_sink_ind = closest_form_ind
                                post_sink_ind = closest_form_ind + 1
                            else:
                                if closest_form_ind != 0:
                                    pre_sink_ind = closest_form_ind - 1
                                else:
                                    pre_sink_ind = closest_form_ind
                                post_sink_ind = closest_form_ind
                            closes_sys_form_ind = np.argmin(abs(file_times.value - sys_form_time))
                            if file_times[closes_sys_form_ind] < sys_form_time:
                                pre_sys_form_ind = closes_sys_form_ind
                                post_sys_form_ind = closes_sys_form_ind + 1
                            else:
                                if closes_sys_form_ind != 0:
                                    pre_sys_form_ind = closes_sys_form_ind - 1
                                else:
                                    pre_sys_form_ind = closes_sys_form_ind
                                post_sys_form_ind = closes_sys_form_ind
                            usable_files = [files[post_sys_form_ind], files[pre_sys_form_ind], files[post_sink_ind], files[pre_sink_ind]]
                        
                        sys.stdout.flush()
                        CW.Barrier()

                        
                        #usable_files = files[pre_form_ind:post_form_ind+1][::-1]
                        pickle_file_preffix = 'triple_'+str(flatten(eval(sys_key)))[1:-1]+'_'
                        pickle_file_preffix = pickle_file_preffix.replace(', ', '_')
                        if "'" in pickle_file_preffix:
                            pickle_file_preffix = pickle_file_preffix.replace("'", "")
                            
                        center_sink = np.max(flatten(eval(sys_key)))
                        center_positions = []
                        max_seps = []
                        for fn in usable_files: #range(len(usable_files)):
                            pit = len(usable_files) - usable_files.index(fn)
                            pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
                            if os.path.exists(pickle_file) == False:
                                print('Getting sink positions from', fn, 'on rank', rank)
                                #fn = usable_files[fn_it]
                                file_no = int(fn.split('output_')[-1].split('/')[0])
                                datadir = fn.split('output_')[0]
                                loaded_sink_data = rsink(file_no, datadir=datadir)
                                if fn == usable_files[1] and size == 1:
                                    import pdb
                                    pdb.set_trace()
                                try:
                                    if np.isnan(center_sink):
                                        import pdb
                                        pdb.set_trace()
                                    else:
                                        center_pos = yt.YTArray([loaded_sink_data['x'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['y'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['z'][center_sink]*units['length_unit'].in_units('au')])
                                        center_vel = yt.YTArray([loaded_sink_data['ux'][center_sink]*units['velocity_unit'].in_units('au/yr'), loaded_sink_data['uy'][center_sink]*units['velocity_unit'].in_units('au/yr'), loaded_sink_data['uz'][center_sink]*units['velocity_unit'].in_units('au/yr')])
                                        sink_creation_time_pick = loaded_sink_data['tcreate'][center_sink]*units['time_unit'].in_units('yr')
                                        t_snap = loaded_sink_data['snapshot_time']*units['time_unit'].in_units('yr')
                                    center_positions.append(center_pos)
                                except:
                                    curr_time = loaded_sink_data['snapshot_time']*units['time_unit'].in_units('yr')
                                    dt = curr_time - t_snap
                                    prev_center_pos = center_positions[-1]
                                    center_pos = prev_center_pos + center_vel*dt
                                    for new_pos_it in range(len(center_pos)):
                                        if center_pos[new_pos_it] < 0:
                                            center_pos[new_pos_it] = units['length_unit'].in_units('au') + center_pos[new_pos_it]
                                        elif center_pos[new_pos_it] > units['length_unit'].in_units('au'):
                                            center_pos[new_pos_it] = center_pos[new_pos_it] - units['length_unit'].in_units('au')
                                    center_positions.append(center_pos)
                                    sink_creation_time_pick = np.nan
                                existing_sinks = list(set(flatten(eval(sys_key))).intersection(np.arange(len(loaded_sink_data['m']))))
                                if len(existing_sinks)>0:
                                    particle_masses = loaded_sink_data['m'][existing_sinks]*units['mass_unit'].in_units('Msun')
                                    particle_x_pos = loaded_sink_data['x'][existing_sinks]*units['length_unit'].in_units('au')
                                    particle_y_pos = loaded_sink_data['y'][existing_sinks]*units['length_unit'].in_units('au')
                                else:
                                    particle_masses = yt.YTArray([], 'Msun')
                                    particle_x_pos = yt.YTArray([], 'au')
                                    particle_y_pos = yt.YTArray([], 'au')
                                try:
                                    youngest_sink = np.argmax(existing_sinks)
                                    dx = np.max(abs(particle_x_pos-particle_x_pos[youngest_sink]))
                                    dy = np.max(abs(particle_y_pos-particle_y_pos[youngest_sink]))
                                    if dx > dy:
                                        max_seps.append(dx)
                                    else:
                                        max_seps.append(dy)
                                except:
                                    pass
                                gc.collect()
                                #particle_masses = dd['sink_particle_mass']
                                
                                if np.isnan(sink_creation_time_pick) == False:
                                    sink_creation_time = sink_creation_time_pick

                                if np.remainder(rank, 3) == 0:
                                    #if np.remainder(rank,48) == 0:
                                    file = open(pickle_file, 'wb')
                                    #pickle.dump((image, time_val, particle_positions, particle_masses), file)
                                    pickle.dump((particle_x_pos, particle_y_pos, particle_masses, max_seps[-1], sink_creation_time_pick, center_pos), file)
                                    file.close()
                                    print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
                                #del x_lim
                                #del y_lim
                                #del z_lim
                                gc.collect()
                            else:
                                file = open(pickle_file, 'rb')
                                particle_x_pos, particle_y_pos, particle_masses, max_sep, sink_creation_time_pick, center_pos = pickle.load(file)
                                file.close()
                                max_seps.append(max_sep)
                                center_positions.append(center_pos)
                                if np.isnan(sink_creation_time_pick) == False:
                                    sink_creation_time = sink_creation_time_pick
                        
                        max_sep = np.max(max_seps)
                            
                        thickness = yt.YTQuantity(np.ceil(max_sep/100)*100+500, 'au')

                        #del units
                        gc.collect()

                        sys.stdout.flush()
                        CW.Barrier()
                        for usable in yt.parallel_objects(usable_files):
                            #for usable in usable_files:
                            pit = len(usable_files) - usable_files.index(usable)
                            pickle_file = pickle_file_preffix + str(pit) + '.pkl'
                            if os.path.exists(pickle_file) == False:
                                print('making projection of', usable, 'on rank', rank)
                                cit = usable_files.index(usable)
                                ds = yt.load(usable, units_override=units_override)
                                #dd = ds.all_data()
                                
                                center_pos = center_positions[cit]
                                time_val = ds.current_time.in_units('yr') - sink_creation_time
                                
                                left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
                                right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
                                region = ds.box(left_corner, right_corner)
                                del left_corner
                                del right_corner
                                gc.collect()
                                
                                axis_ind = 2
                                try:
                                    proj = yt.ProjectionPlot(ds, axis_ind, ("ramses", "Density"), width=thickness, data_source=region, method='integrate', center=(center_pos, 'AU'))
                                    image = (proj.frb.data[("ramses", "Density")]/thickness.in_units('cm')).value*units['density_unit'].in_units('g/cm**3')
                                    del proj
                                        
                                    gc.collect()
                                    
                                    file = open(pickle_file, 'wb')
                                    pickle.dump((image, time_val), file)
                                    file.close()
                                    print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
                                except:
                                    print("Density field doesn't exist")
                                    file = open(pickle_file, 'wb')
                                    pickle.dump((np.zeros((800,800)), time_val), file)
                                    file.close()

                        sys.stdout.flush()
                        CW.Barrier()
                            
                        #"""
                        #pickle_files = sorted(glob.glob(pickle_file_preffix + '*_part.pkl'))
                        if len(usable_files) == 2:
                            pickle_files = [pickle_file_preffix + '1_part.pkl', pickle_file_preffix + '2_part.pkl']
                        else:
                            pickle_files = [pickle_file_preffix + '1_part.pkl', pickle_file_preffix + '2_part.pkl', pickle_file_preffix + '3_part.pkl', pickle_file_preffix + '4_part.pkl']
                        #cit = 0
                        #for pickle_file in pickle_files:
                        for pickle_file in yt.parallel_objects(pickle_files, njobs=2):
                            pit = pickle_files.index(pickle_file)
                            file_name = pickle_file_preffix + ("%06d" % pit)
                            proj_pickle = "".join(pickle_file.split('_part'))
                            file = open(proj_pickle, 'rb')
                            image, time_val = pickle.load(file)
                                #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                            file.close()
                            if os.path.exists(file_name+ ".png") == False and np.nanmean(image) != 0:
                                print('making frame for pickle', pickle_file, 'on rank', rank)
                                file = open(pickle_file, 'rb')
                                particle_x_pos, particle_y_pos, particle_masses, max_seps, sink_creation_time_pick, center_pos = pickle.load(file)
                                #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                                file.close()
                                
                                center_pos = center_positions[::-1][pit]
                                
                                file = open(proj_pickle, 'rb')
                                image, time_val = pickle.load(file)
                                #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                                file.close()
                                
                                xlim = [-1*thickness, thickness]
                                ylim = [-1*thickness, thickness]
                                x = np.linspace(xlim[0], xlim[1], 800)
                                y = np.linspace(ylim[0], ylim[1], 800)
                                X, Y = np.meshgrid(x, y)
                                
                                #X = X + center_pos[0]
                                #Y = Y + center_pos[1]
                                particle_x_pos = particle_x_pos - center_pos[0]
                                particle_y_pos = particle_y_pos - center_pos[1]
                                
                                xlim = [np.min(X), np.max(X)]
                                ylim = [np.min(Y), np.max(Y)]
                                
                                annotate_space = (xlim[1] - xlim[0])/31
                                x_ind = []
                                y_ind = []
                                counter = 0
                                while counter < 31:
                                    val = annotate_space*counter + annotate_space/2. + xlim[0]
                                    x_ind.append(int(val))
                                    y_ind.append(int(val))
                                    counter = counter + 1
                                X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
                                      
                                has_particles = True
                                xabel = "X (AU)"
                                yabel = "Y (AU)"
                                plt.clf()
                                fig, ax = plt.subplots()
                                ax.set_xlabel(xabel, labelpad=-1, fontsize=10)
                                ax.set_ylabel(yabel, fontsize=10) #, labelpad=-20
                                ax.set_xlim(xlim)
                                ax.set_ylim(ylim)

                                
                                cmin = 10**(np.log10(np.mean(image))-1.5)
                                if np.isnan(cmin):
                                    cmin = 1.e-19
                                cmax = 10**(np.log10(np.mean(image))+1.5)
                                if np.isnan(cmax):
                                    cmax = 1.e-16
                                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cmin, vmax=cmax), rasterized=True)
                                plt.gca().set_aspect('equal')
                                plt.gca().set_aspect('equal')
                                #plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                                cbar = plt.colorbar(plot, pad=0.0)
                                plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)
                                #ax.scatter(particle_x_pos, particle_y_pos, color='c', s=1)
                                '''
                                try:
                                    mym.annotate_particles(ax, np.array([particle_x_pos, particle_y_pos]), 200, limits=[xlim, ylim], annotate_field=particle_masses)
                                except:
                                    try:
                                        mym.annotate_particles(ax, np.array([[particle_x_pos], [particle_y_pos]]), 200, limits=[xlim, ylim], annotate_field=[particle_masses])
                                    except:
                                        pass
                                '''
                                
                                cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
                                plt.savefig(file_name + ".png", format='png', bbox_inches='tight')

                                plt.tick_params(axis='both', which='major')# labelsize=16)
                                for line in ax.xaxis.get_ticklines():
                                    line.set_color('white')
                                for line in ax.yaxis.get_ticklines():
                                    line.set_color('white')

                                    
                                #part_color = ['cyan','magenta','r','b','y','w','k']
                                for part_it in range(len(particle_x_pos)):
                                    ax.scatter(particle_x_pos[part_it], particle_y_pos[part_it], marker='*', s=100, linewidth=1.5, edgecolor="k", zorder=11)
                                plt.savefig(file_name + ".png", format='png', bbox_inches='tight')

                                
                                try:
                                    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                    time_string = "$t$="+str(int(time_val))+"yr"
                                    time_string_raw = r"{}".format(time_string)
                                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
                                    try:
                                        plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                        time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                                    except:
                                        print("Couldn't outline time string")
                                except:
                                    print("Couldn't plot time string")
                                
                                
                                file = open(pickle_file.split('_part.pkl')[0]+'_all.pkl', 'wb')
                                pickle.dump((sys_key, particle_x_pos, particle_y_pos, particle_masses, center_pos, thickness, X, Y, image, time_val), file)
                                file.close()
                                    
                                try:
                                    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                except:
                                    pass
                                #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                                print('Created frame ' + file_name + '.png, on rank', rank)
                            #"""
            if len(flatten(eval(sys_key))) == 4 and str(np.max(flatten(eval(sys_key)))) in Sink_birth_all.keys():

                youngest_birth_con = Sink_birth_all[str(np.max(flatten(eval(sys_key))))]
                sys_brackets = ''
                for char in sys_key:
                    if char == '[' or char == ']':
                        sys_brackets = sys_brackets + char
                if sys_brackets == '[[][]]':
                    if np.max(superplot_dict['System_seps'][sys_key][0][:2]) < 600 and np.min(superplot_dict['System_seps'][sys_key][0][0]) > 100 and superplot_dict['System_seps'][sys_key][0][2] > 300 and superplot_dict['System_seps'][sys_key][0][2] < 3000:
                        if youngest_birth_con[0] == True or youngest_birth_con[1]==youngest_birth_con[2]:
                            print("Found Quadruple candidate:", sys_key)
                            if youngest_birth_con[0] == True:
                                goal_time = superplot_dict['System_times'][sys_key][0]
                                closest_time_ind = np.argmin(abs(file_times.value - goal_time))
                                if file_times[closest_time_ind] < goal_time:
                                    pre_form_ind = closest_time_ind
                                    post_form_ind = closest_time_ind + 1
                                else:
                                    post_form_ind = closest_time_ind
                                    pre_form_ind = closest_time_ind - 1
                                usable_files = [files[post_form_ind], files[pre_form_ind]]
                            else:
                                form_time = youngest_birth_con[-1]*units['time_unit'].in_units('yr')
                                sys_form_time = superplot_dict['System_times'][sys_key][0]
                                closest_form_ind = np.argmin(abs(file_times.value - form_time.value))
                                if file_times[closest_form_ind] < form_time:
                                    pre_sink_ind = closest_form_ind
                                    post_sink_ind = closest_form_ind + 1
                                else:
                                    pre_sink_ind = closest_form_ind - 1
                                    post_sink_ind = closest_form_ind
                                closes_sys_form_ind = np.argmin(abs(file_times.value - sys_form_time))
                                if file_times[closes_sys_form_ind] < sys_form_time:
                                    pre_sys_form_ind = closes_sys_form_ind
                                    post_sys_form_ind = closes_sys_form_ind + 1
                                else:
                                    pre_sys_form_ind = closes_sys_form_ind - 1
                                    post_sys_form_ind = closes_sys_form_ind
                                usable_files = [files[post_sys_form_ind], files[pre_sys_form_ind], files[post_sink_ind], files[pre_sink_ind]]
                            
                            
                            sys.stdout.flush()
                            CW.Barrier()

                            #usable_files = files[pre_form_ind:post_form_ind+1][::-1]
                            pickle_file_preffix = 'quad_'+str(flatten(eval(sys_key)))[1:-1]+'_'
                            pickle_file_preffix = pickle_file_preffix.replace(', ', '_')
                            if "'" in pickle_file_preffix:
                                pickle_file_preffix = pickle_file_preffix.replace("'", "")
                                
                            center_sink = np.max(flatten(eval(sys_key)))
                            center_positions = []
                            max_seps = []
                            for fn in yt.parallel_objects(usable_files, njobs=int(2)): #range(len(usable_files)):
                                pit = len(usable_files) - usable_files.index(fn)
                                pickle_file = pickle_file_preffix + str(pit) + '_part.pkl'
                                if os.path.exists(pickle_file) == False:
                                    print('Getting sink positions from', fn, 'on rank', rank)
                                    #fn = usable_files[fn_it]
                                    file_no = int(fn.split('output_')[-1].split('/')[0])
                                    datadir = fn.split('output_')[0]
                                    loaded_sink_data = rsink(file_no, datadir=datadir)
                                    try:
                                        if np.isnan(center_sink):
                                            import pdb
                                            pdb.set_trace()
                                        else:
                                            center_pos = yt.YTArray([loaded_sink_data['x'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['y'][center_sink]*units['length_unit'].in_units('au'), loaded_sink_data['z'][center_sink]*units['length_unit'].in_units('au')])
                                            center_vel = yt.YTArray([loaded_sink_data['ux'][center_sink]*units['velocity_unit'].in_units('au/yr'), loaded_sink_data['uy'][center_sink]*units['velocity_unit'].in_units('au/yr'), loaded_sink_data['uz'][center_sink]*units['velocity_unit'].in_units('au/yr')])
                                            sink_creation_time_pick = loaded_sink_data['tcreate'][center_sink]*units['time_unit'].in_units('yr')
                                            t_snap = loaded_sink_data['snapshot_time']*units['time_unit'].in_units('yr')
                                        center_positions.append(center_pos)
                                    except:
                                        curr_time = loaded_sink_data['snapshot_time']*units['time_unit'].in_units('yr')
                                        dt = curr_time - t_snap
                                        prev_center_pos = center_positions[-1]
                                        center_pos = prev_center_pos + center_vel*dt
                                        for new_pos_it in range(len(center_pos)):
                                            if center_pos[new_pos_it] < 0:
                                                center_pos[new_pos_it] = units['length_unit'].in_units('au') + center_pos[new_pos_it]
                                            elif center_pos[new_pos_it] > units['length_unit'].in_units('au'):
                                                center_pos[new_pos_it] = center_pos[new_pos_it] - units['length_unit'].in_units('au')
                                        center_positions.append(center_pos)
                                        sink_creation_time_pick = np.nan
                                    existing_sinks = list(set(flatten(eval(sys_key))).intersection(np.arange(len(loaded_sink_data['m']))))
                                    if len(existing_sinks)>0:
                                        particle_masses = loaded_sink_data['m'][existing_sinks]*units['mass_unit'].in_units('Msun')
                                        particle_x_pos = loaded_sink_data['x'][existing_sinks]*units['length_unit'].in_units('au')
                                        particle_y_pos = loaded_sink_data['y'][existing_sinks]*units['length_unit'].in_units('au')
                                    else:
                                        particle_masses = yt.YTArray([], 'Msun')
                                        particle_x_pos = yt.YTArray([], 'au')
                                        particle_y_pos = yt.YTArray([], 'au')
                                    try:
                                        youngest_sink = np.argmax(existing_sinks)
                                        dx = np.max(abs(particle_x_pos-particle_x_pos[youngest_sink]))
                                        dy = np.max(abs(particle_y_pos-particle_y_pos[youngest_sink]))
                                        if dx > dy:
                                            max_seps.append(dx)
                                        else:
                                            max_seps.append(dy)
                                    except:
                                        pass
                                    gc.collect()
                                    #particle_masses = dd['sink_particle_mass']
                                    
                                    if np.isnan(sink_creation_time_pick) == False:
                                        sink_creation_time = sink_creation_time_pick

                                    if np.remainder(rank, 3) == 0:
                                        #if np.remainder(rank,48) == 0:
                                        file = open(pickle_file, 'wb')
                                        #pickle.dump((image, time_val, particle_positions, particle_masses), file)
                                        pickle.dump((particle_x_pos, particle_y_pos, particle_masses, max_seps[-1], sink_creation_time_pick, center_pos), file)
                                        file.close()
                                        print("Created Pickle:", pickle_file, "for  file:", fn, "on rank", rank)
                                    #del x_lim
                                    #del y_lim
                                    #del z_lim
                                    gc.collect()
                                else:
                                    file = open(pickle_file, 'rb')
                                    particle_x_pos, particle_y_pos, particle_masses, max_sep, sink_creation_time_pick, center_pos = pickle.load(file)
                                    file.close()
                                    max_seps.append(max_sep)
                                    center_positions.append(center_pos)
                                    if np.isnan(sink_creation_time_pick) == False:
                                        sink_creation_time = sink_creation_time_pick
                                        
                                    
                            max_sep = np.max(max_seps)
                                
                            thickness = yt.YTQuantity(np.ceil(max_sep/100)*100+500, 'au')

                            #del units
                            gc.collect()

                            sys.stdout.flush()
                            CW.Barrier()
                            for usable in yt.parallel_objects(usable_files):
                                #for usable in usable_files:
                                pit = len(usable_files) - usable_files.index(usable)
                                pickle_file = pickle_file_preffix + str(pit) + '.pkl'
                                if os.path.exists(pickle_file) == False:
                                    print('making projection of', usable, 'on rank', rank)
                                    cit = usable_files.index(usable)
                                    ds = yt.load(usable, units_override=units_override)
                                    #dd = ds.all_data()

                                    center_pos = center_positions[cit]
                                    time_val = ds.current_time.in_units('yr') - sink_creation_time
                                    
                                    left_corner = yt.YTArray([center_pos[0]-(0.75*thickness), center_pos[1]-(0.75*thickness), center_pos[2]-(0.5*thickness)], 'AU')
                                    right_corner = yt.YTArray([center_pos[0]+(0.75*thickness), center_pos[1]+(0.75*thickness), center_pos[2]+(0.5*thickness)], 'AU')
                                    region = ds.box(left_corner, right_corner)
                                    del left_corner
                                    del right_corner
                                    gc.collect()
                                    
                                    axis_ind = 2
                                    try:
                                        proj = yt.ProjectionPlot(ds, axis_ind, ("ramses", "Density"), width=thickness, data_source=region, method='integrate', center=(center_pos, 'AU'))
                                        image = (proj.frb.data[("ramses", "Density")]/thickness.in_units('cm')).value*units['density_unit'].in_units('g/cm**3')
                                        del proj
                                        
                                        gc.collect()
                                    
                                        file = open(pickle_file, 'wb')
                                        pickle.dump((image, time_val), file)
                                        file.close()
                                        print("Created Pickle:", pickle_file, "for  file:", str(ds), "on rank", rank)
                                    except:
                                        print("Density field doesn't exist")
                                        file = open(pickle_file, 'wb')
                                        pickle.dump((np.zeros((800,800)), time_val), file)
                                        file.close()
                                        

                            sys.stdout.flush()
                            CW.Barrier()
                                
                                
                            #pickle_files = sorted(glob.glob(pickle_file_preffix + '*_part.pkl'))
                            if len(usable_files) == 2:
                                pickle_files = [pickle_file_preffix + '1_part.pkl', pickle_file_preffix + '2_part.pkl']
                            else:
                                pickle_files = [pickle_file_preffix + '1_part.pkl', pickle_file_preffix + '2_part.pkl', pickle_file_preffix + '3_part.pkl', pickle_file_preffix + '4_part.pkl']
                            #cit = 0
                            #for pickle_file in pickle_files:
                            for pickle_file in yt.parallel_objects(pickle_files, njobs=2):
                                pit = pickle_files.index(pickle_file)
                                file_name = pickle_file_preffix + ("%06d" % pit)
                                proj_pickle = "".join(pickle_file.split('_part'))
                                file = open(proj_pickle, 'rb')
                                image, time_val = pickle.load(file)
                                #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel =   pickle.load(file)
                                file.close()
                                if os.path.exists(file_name+ ".png") == False and np.nanmean(image) != 0:
                                    print('making frame for pickle', pickle_file, 'on rank', rank)
                                    file = open(pickle_file, 'rb')
                                    particle_x_pos, particle_y_pos, particle_masses, max_seps, sink_creation_time_pick, center_pos = pickle.load(file)
                                    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                                    file.close()
                                    
                                    center_pos = center_positions[::-1][pit]
                                    
                                    file = open(proj_pickle, 'rb')
                                    image, time_val = pickle.load(file)
                                    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
                                    file.close()
                                    
                                    xlim = [-1*thickness, thickness]
                                    ylim = [-1*thickness, thickness]
                                    x = np.linspace(xlim[0], xlim[1], 800)
                                    y = np.linspace(ylim[0], ylim[1], 800)
                                    X, Y = np.meshgrid(x, y)
                                    
                                    #X = X + center_pos[0]
                                    #Y = Y + center_pos[1]
                                    particle_x_pos = particle_x_pos - center_pos[0]
                                    particle_y_pos = particle_y_pos - center_pos[1]
                                    
                                    xlim = [np.min(X), np.max(X)]
                                    ylim = [np.min(Y), np.max(Y)]
                                    
                                    annotate_space = (xlim[1] - xlim[0])/31
                                    x_ind = []
                                    y_ind = []
                                    counter = 0
                                    while counter < 31:
                                        val = annotate_space*counter + annotate_space/2. + xlim[0]
                                        x_ind.append(int(val))
                                        y_ind.append(int(val))
                                        counter = counter + 1
                                    X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
                                          
                                    has_particles = True
                                    xabel = "X (AU)"
                                    yabel = "Y (AU)"
                                    plt.clf()
                                    fig, ax = plt.subplots()
                                    ax.set_xlabel(xabel, labelpad=-1, fontsize=10)
                                    ax.set_ylabel(yabel, fontsize=10) #, labelpad=-20
                                    ax.set_xlim(xlim)
                                    ax.set_ylim(ylim)

                                    
                                    cmin = 10**(np.log10(np.mean(image))-1.5)
                                    if np.isnan(cmin):
                                        cmin = 1.e-19
                                    cmax = 10**(np.log10(np.mean(image))+1.5)
                                    if np.isnan(cmax):
                                        cmax = 1.e-16
                                    plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cmin, vmax=cmax), rasterized=True)
                                    plt.gca().set_aspect('equal')
                                    plt.gca().set_aspect('equal')
                                    #plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                                    cbar = plt.colorbar(plot, pad=0.0)
                                    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                    #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=velz)
                                    #ax.scatter(particle_x_pos, particle_y_pos, color='c', s=1)
                                    '''
                                    try:
                                        mym.annotate_particles(ax, np.array([particle_x_pos, particle_y_pos]), 200, limits=[xlim, ylim], annotate_field=particle_masses)
                                    except:
                                        try:
                                            mym.annotate_particles(ax, np.array([[particle_x_pos], [particle_y_pos]]), 200, limits=[xlim, ylim], annotate_field=[particle_masses])
                                        except:
                                            pass
                                    '''
                                    
                                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
                                    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')

                                    plt.tick_params(axis='both', which='major')# labelsize=16)
                                    for line in ax.xaxis.get_ticklines():
                                        line.set_color('white')
                                    for line in ax.yaxis.get_ticklines():
                                        line.set_color('white')
                                        
                                    #part_color = ['cyan','magenta','r','b','y','w','k']
                                    for part_it in range(len(particle_x_pos)):
                                        ax.scatter(particle_x_pos[part_it], particle_y_pos[part_it], marker='*', s=100, linewidth=1.5, edgecolor="k", zorder=11)
                                    plt.savefig(file_name + ".png", format='png', bbox_inches='tight')

                                    try:
                                        plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                        time_string = "$t$="+str(int(time_val))+"yr"
                                        time_string_raw = r"{}".format(time_string)
                                        time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
                                        try:
                                            plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                                        except:
                                            print("Couldn't outline time string")
                                    except:
                                        print("Couldn't plot time string")
                                        
                                    
                                    file = open(pickle_file.split('_part.pkl')[0]+'_all.pkl', 'wb')
                                    pickle.dump((sys_key, particle_x_pos, particle_y_pos, particle_masses, center_pos, thickness, X, Y, image, time_val), file)
                                    file.close()
                                        

                                    try:
                                        plt.savefig(file_name + ".png", format='png', bbox_inches='tight')
                                    except:
                                        pass
                                    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                                    print('Created frame ' + file_name + '.png, on rank', rank)
                                

