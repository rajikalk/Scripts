import numpy as np
import matplotlib.pyplot as plt
import pyramses as pr
from pyramses import rsink
import multiplicity as m
import yt
import glob
from mpi4py.MPI import COMM_WORLD as CW
import sys
import collections
import os
import pickle
import gc

f_acc= 0.5

def losi(i, res):
    if (res['n'][i]==1) or (res['n'][i]==0):
        return i
    else:
        i1 = losi(res['index1'][i],res)
        i2 = losi(res['index2'][i],res)
        return [i1,i2]

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def accretion(sink_inds, time_ind):
    """
    Calculates the accretion of the given indeexes
    """
    global Accretion_array
    M_dot = Accretion_array[time_ind, sink_inds]
    return M_dot
    
def luminosity(global_data, sink_inds, global_ind):
    """
    Calculates the luminosity of the given indexes
    """
    global f_acc
    radius = yt.YTQuantity(2.0, 'rsun')
    M_dot = accretion(sink_inds, global_ind)
    M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
    L_acc = f_acc * (yt.units.G * M.in_units('g') * M_dot.in_units('g/s'))/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    return L_tot

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-global_data", "--global_data_pickle_file", help="Where is the directory of the global pickle data?", default='/groups/astro/rlk/Analysis_plots/Ramses/Global/G100/512/stars_red_512.pkl', type=str)
    parser.add_argument("-pickle", "--pickled_file", help="Define if you want to read this instead", type=str)
    parser.add_argument("-update", "--update_pickles", help="Do you want to remake the pickles?", type=str, default='True')
    parser.add_argument("-acc_lim", "--accretion_limit", help="What do you want to set the accretion limit to?", type=float, default=1.e-7)
    parser.add_argument("-upper_L", "--upper_L_limit", help="What is the upper Luminosity limit?", type=float, default=55.29)
    parser.add_argument("-lower_L", "--lower_L_limit", help="What is the upper Luminosity limit?", type=float, default=0.09)
    parser.add_argument("-bound", "--bound_check", help="Do you actually want to analyse bound systems?", type=str, default='True')
    parser.add_argument("-max_iter", "--max_iterations", help="How many iterations for multiplicity analysis", type=int, default=20)
    parser.add_argument("-start_time_it", "--start_time_index", help="What time index do you want to start at? mostly for diagnostic reasons.", type=int, default=0)
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifeimteimte threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("-plot_super", "--plot_superplot", help="do you want to plot the superplot?", type=str, default='True')
    parser.add_argument("-use_com_sep", "--use_separation_of_the_coms", help="Do you want to same the separation as the separation of the center of masses instead of relative to a particular sink?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

args = parse_inputs()

rank = CW.Get_rank()
size = CW.Get_size()


#Set units
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s")}

simulation_density_id = args.global_data_pickle_file.split('/')[-1].split('_')[0][1:]
try:
    simulation_density_id_int = int(simulation_density_id)
except:
    simulation_density_id = args.global_data_pickle_file.split('_G')[-1].split('.pkl')[0]

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
    

units_override.update({"density_unit":(units_override['mass_unit'][0]/(units_override['length_unit'][0]**3), "Msun/pc**3")})
    
scale_l = yt.YTQuantity(units_override['length_unit'][0], units_override['length_unit'][1]).in_units('cm') # 4 pc
scale_v = yt.YTQuantity(units_override['velocity_unit'][0], units_override['velocity_unit'][1]).in_units('cm/s')         # 0.18 km/s == sound speed
scale_t = scale_l/scale_v # 4 pc / 0.18 km/s
scale_d = yt.YTQuantity(units_override['density_unit'][0], units_override['density_unit'][1]).in_units('g/cm**3')  # 2998 Msun / (4 pc)^3

units={}
for key in units_override.keys():
    units.update({key:yt.YTQuantity(units_override[key][0], units_override[key][1])})
del units_override
gc.collect()

sys.stdout.flush()
CW.Barrier()

#Define arrays of quantities to be saved.
Times = []
SFE = []
N_vis_stars_UL = []
N_vis_stars_NUL = []

sys.stdout.flush()
CW.Barrier()

#loading global data
file_open = open(args.global_data_pickle_file, 'rb')
try:
    global_data = pickle.load(file_open,encoding="latin1")
except:
    file_open.close()
    import pickle5 as pickle
    file_open = open(args.global_data_pickle_file, 'rb')
    global_data = pickle.load(file_open,encoding="latin1")
file_open.close()
dm = global_data['dm']*units['mass_unit'].in_units('Msun')
dt = (global_data['time'] - global_data['tflush'])*units['time_unit'].in_units('yr')
Accretion_array = dm/dt
print('loaded global data')

sys.stdout.flush()
CW.Barrier()

#Get birth conditions
file_open = open(args.global_data_pickle_file, 'rb')
#Sink_bound_birth.append([born_bound, most_bound_sink_id, rel_sep])

sys.stdout.flush()
CW.Barrier()

if args.bound_check == 'True':
    bound_check = True
else:
    bound_check = False

sys.stdout.flush()
CW.Barrier()

#Find sink formation times
Sink_formation_times = []
time_it = 0
sink_it = 0
while time_it < np.shape(global_data['m'])[0]:
     while sink_it < np.shape(global_data['m'])[1]:
         if global_data['m'][time_it][sink_it]>0:
             Sink_formation_times.append(global_data['time'][time_it][sink_it])
             sink_it = sink_it + 1
         elif global_data['m'][time_it][sink_it] == 0:
             time_it = time_it + 1
             break
     if sink_it == np.shape(global_data['m'])[1]:
         break
Sink_formation_times = (Sink_formation_times*scale_t).in_units('yr')

sys.stdout.flush()
CW.Barrier()

#Define pickle name
pickle_file = args.pickled_file
            
sys.stdout.flush()
CW.Barrier()

start_time_it = args.start_time_index

if args.update_pickles == 'True':
    rit = -1
    for time_it in range(start_time_it, len(global_data['time'].T[0])):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            #preamble
            n_stars = np.where(global_data['m'][time_it]>0)[0]
            if len(n_stars) > 1:
                abspos = np.array([global_data['x'][time_it][n_stars], global_data['y'][time_it][n_stars], global_data['z'][time_it][n_stars]]).T#*scale_l
                absvel = np.array([global_data['ux'][time_it][n_stars], global_data['uy'][time_it][n_stars], global_data['uz'][time_it][n_stars]]).T#*scale_v
                mass = np.array(global_data['m'][time_it][n_stars])
                time = global_data['time'][time_it][n_stars][0]
                time_yr = yt.YTQuantity(scale_t*time, 's').in_units('yr').value
                Times.append(time_yr)
                SFE.append(np.sum(mass))

                L_tot = luminosity(global_data, n_stars, time_it)
                M_dot = accretion(n_stars, time_it)
                vis_upper_limit = np.where((L_tot>=0.09)&(M_dot>=1.e-7)&(L_tot<=55.29))[0]#Remember to update if you adjust criteria
                vis_no_upper_limit = np.where((L_tot>=0.09)&(M_dot>=1.e-7))[0]
                N_vis_stars.append(len(vis_no_upper_limit), )
            
    sys.stdout.flush()
    CW.Barrier()

#compile all the pickles:
if rank == 0:
    pickle_files = sorted(glob.glob(pickle_file.split('.pkl')[0] + "_*.pkl"))
    if len(pickle_files) > 0:
        file = open(pickle_files[0], 'rb')
        superplot_dict = pickle.load(file)
        file.close()
        
        full_dict = {}
        for key in superplot_dict.keys():
            if 'System' in key:
                full_dict.update({key+'_full':{}})
            else:
                full_dict.update({key+'_full':[]})

        for pick_file in pickle_files:
            file = open(pick_file, 'rb')
            superplot_dict = pickle.load(file)
            file.close()
            
            try:
                for time_key in superplot_dict['System_times'].keys():
                    if time_key not in full_dict['System_times_full'].keys():
                        for full_key in full_dict.keys():
                            if 'System' in full_key:
                                full_dict[full_key].update({time_key:superplot_dict[full_key.split('_full')[0]][time_key]})
                    else:
                        for full_key in full_dict.keys():
                            if 'System' in full_key:
                                full_dict[full_key][time_key] = full_dict[full_key][time_key] + superplot_dict[full_key.split('_full')[0]][time_key]
                for full_key in full_dict.keys():
                    if 'System' not in full_key:
                        full_dict[full_key] = full_dict[full_key] + superplot_dict[full_key.split('_full')[0]]
            except:
                for time_key in superplot_dict[0]['System_times'].keys():
                    if time_key not in full_dict['System_times_full'].keys():
                        for full_key in full_dict.keys():
                            if 'System' in full_key:
                                full_dict[full_key].update({time_key:superplot_dict[0][full_key.split('_full')[0]][time_key]})
                    else:
                        for full_key in full_dict.keys():
                            if 'System' in full_key:
                                full_dict[full_key][time_key] = full_dict[full_key][time_key] + superplot_dict[0][full_key.split('_full')[0]][time_key]
                for full_key in full_dict.keys():
                    if 'System' not in full_key:
                        full_dict[full_key] = full_dict[full_key] + superplot_dict[0][full_key.split('_full')[0]]
            os.remove(pick_file)
        
        #Let's sort the data
        for time_key in full_dict['System_times_full'].keys():
            for full_key in full_dict.keys():
                if 'System' in full_key:
                    dict_sort = np.argsort(np.array(full_dict['System_times_full'][time_key]).T)
                    sorted_array = np.array(full_dict[full_key][time_key])[dict_sort]
                    try:
                        superplot_dict[full_key.split('_full')[0]][time_key] = sorted_array.tolist()
                    except:
                        superplot_dict[0][full_key.split('_full')[0]][time_key] = sorted_array.tolist()
            
        sorted_inds = np.argsort(full_dict['Times_full'])
        try:
            for super_key in superplot_dict.keys():
                if 'System' not in super_key:
                    try:
                        superplot_dict[super_key] = np.array(full_dict[super_key+'_full'])[sorted_inds].tolist()
                    except:
                        print(super_key, 'does not exist')
        except:
            for super_key in superplot_dict[0].keys():
                if 'System' not in super_key:
                    try:
                        superplot_dict[0][super_key] = np.array(full_dict[super_key+'_full'])[sorted_inds].tolist()
                    except:
                        print(super_key, 'does not exist')
        
        #sort keys into chronological order of system formation
        Start_times = []
        Sort_keys = []
        try:
            for time_key in superplot_dict['System_times'].keys():
                Start_times.append(superplot_dict['System_times'][time_key][0])
                Sort_keys.append(time_key)
        except:
            for time_key in superplot_dict[0]['System_times'].keys():
                Start_times.append(superplot_dict[0]['System_times'][time_key][0])
                Sort_keys.append(time_key)
        sorted_inds = np.argsort(Start_times)
        Sorted_keys = np.array(Sort_keys)[sorted_inds]
        del Sort_keys
        del Start_times
        
        System_seps = {}
        System_midpoint_seps = {}
        System_semimajor = {}
        System_times = {}
        System_ecc = {}
        System_energies = {}
        for sorted_key in Sorted_keys:
            try:
                System_seps.update({sorted_key:superplot_dict['System_seps'][sorted_key]})
                System_midpoint_seps.update({sorted_key:superplot_dict['System_midpoint_seps'][sorted_key]})
                System_semimajor.update({sorted_key:superplot_dict['System_semimajor'][sorted_key]})
                System_times.update({sorted_key:superplot_dict['System_times'][sorted_key]})
                System_ecc.update({sorted_key:superplot_dict['System_ecc'][sorted_key]})
                System_energies.update({sorted_key:superplot_dict['System_energies'][sorted_key]})
            except:
                System_seps.update({sorted_key:superplot_dict[0]['System_seps'][sorted_key]})
                System_midpoint_seps.update({sorted_key:superplot_dict[0]['System_midpoint_seps'][sorted_key]})
                System_semimajor.update({sorted_key:superplot_dict[0]['System_semimajor'][sorted_key]})
                System_times.update({sorted_key:superplot_dict[0]['System_times'][sorted_key]})
                System_ecc.update({sorted_key:superplot_dict[0]['System_ecc'][sorted_key]})
                System_energies.update({sorted_key:superplot_dict[0]['System_energies'][sorted_key]})
        del Sorted_keys
        try:
            superplot_dict['System_seps'] = System_seps
            superplot_dict['System_midpoint_seps'] = System_midpoint_seps
            superplot_dict['System_semimajor'] = System_semimajor
            superplot_dict['System_times'] = System_times
            superplot_dict['System_ecc'] = System_ecc
            superplot_dict['System_energies'] = System_energies
            file = open(pickle_file+'.pkl', 'wb')
            pickle.dump((superplot_dict, Sink_bound_birth, Sink_formation_times),file)
            file.close()
        except:
            superplot_dict[0]['System_seps'] = System_seps
            superplot_dict[0]['System_midpoint_seps'] = System_midpoint_seps
            superplot_dict[0]['System_semimajor'] = System_semimajor
            superplot_dict[0]['System_times'] = System_times
            superplot_dict[0]['System_ecc'] = System_ecc
            superplot_dict[0]['System_energies'] = System_energies
            file = open(pickle_file+'.pkl', 'wb')
            pickle.dump((superplot_dict[0], Sink_bound_birth, Sink_formation_times),file)
            file.close()
        

del System_seps
del System_midpoint_seps
del System_semimajor
del System_ecc
del System_energies

sys.stdout.flush()
CW.Barrier()

print("gathered and sorted pickles and saved to", pickle_file+'.pkl')

#calculate means:
print("calculating means")
rit = -1
while rit < size:
    rit = rit + 1
    if rank == rit:
        file = open(pickle_file+'.pkl', 'rb')
        superplot_dict, Sink_bound_birth, Sink_formation_times = pickle.load(file)
        file.close()

sys.stdout.flush()
CW.Barrier()

#if args.update_pickles == 'True':
means_dict = {}
for super_key in superplot_dict.keys():
    if 'System' in super_key:
        means_dict.update({'System_mean_'+super_key.split('System_')[-1]:{}})
Lifetimes_sys = {}
Sep_maxs = []
Sep_mins = []
Initial_Seps = []
Final_seps = []
window = 50000#year
rit = -1
for time_key in superplot_dict['System_times'].keys():
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        mean_sys_dict = {}
        for super_key in superplot_dict.keys():
            if 'System' in super_key:
                mean_sys_dict.update({'mean_'+super_key:[]})
        for time in superplot_dict['System_times'][time_key]:
            start_time = time - window/2.
            end_time = time + window/2.
            start_ind = np.argmin(abs(np.array(superplot_dict['System_times'][time_key]) - start_time))
            end_ind = np.argmin(abs(np.array(superplot_dict['System_times'][time_key]) - end_time))
            if end_ind != start_ind:
                non_nan_inds = np.argwhere(np.isnan(np.array(superplot_dict['System_seps'][time_key][start_ind:end_ind])[:,0])==False)
                for super_key in superplot_dict.keys():
                    if 'System' in super_key:
                        if 'times' in super_key:
                            mean_t = np.mean(np.array((superplot_dict['System_times'][time_key]- np.array(superplot_dict['Times'][0]))[start_ind:end_ind])[non_nan_inds])
                            mean_sys_dict['mean_'+super_key].append(mean_t)
                        else:
                            mean_val = np.mean(np.array(superplot_dict[super_key][time_key][start_ind:end_ind])[non_nan_inds], axis=0)
                            mean_sys_dict['mean_'+super_key].append(mean_val[0])
        non_nan_inds = np.argwhere(np.isnan(np.array(superplot_dict['System_seps'][time_key]).T[0])==False).T[0]
        sys_life_time = np.array(superplot_dict['System_times'][time_key])[non_nan_inds][-1] - np.array(superplot_dict['System_times'][time_key])[non_nan_inds][0]
        Lifetimes_sys.update({time_key:sys_life_time})
        if sys_life_time > args.sys_lifetime_threshold:
            end_sep_time = np.array(superplot_dict['System_times'][time_key] - np.array(superplot_dict['Times'][0]))[non_nan_inds][0] + sys_life_time - 1000
            end_time_it = np.argmin(abs(np.array(superplot_dict['System_times'][time_key]- np.array(superplot_dict['Times'][0]))[non_nan_inds] - end_sep_time))
            for sys_it in range(np.shape(np.array(mean_sys_dict['mean_System_seps']).T)[0]):
                sep_max = np.nanmax(np.array(superplot_dict['System_seps'][time_key]).T[sys_it])
                sep_min = np.nanmin(np.array(superplot_dict['System_seps'][time_key]).T[sys_it])
                sep_init = np.array(superplot_dict['System_seps'][time_key]).T[sys_it][0]
                sep_final = np.mean(np.array(superplot_dict['System_seps'][time_key]).T[sys_it][non_nan_inds][end_time_it:])
                Sep_maxs.append(sep_max)
                Sep_mins.append(sep_min)
                Initial_Seps.append(sep_init)
                Final_seps.append(sep_final)
        for mean_key in means_dict.keys():
            means_dict[mean_key].update({time_key:mean_sys_dict['mean_'+'_'.join(mean_key.split('_mean_'))]})
        print("calculated means for", time_key)

sys.stdout.flush()
CW.Barrier()

#save pickles and compile means
pickle_file_rank = "means_" + pickle_file.split('.pkl')[0] + "_" + ("%03d" % rank) + ".pkl"
file = open(pickle_file_rank, 'wb')
pickle.dump((means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps), file)
file.close()
print("saved_means")

sys.stdout.flush()
CW.Barrier()

if rank == 0:
    pickle_files = sorted(glob.glob("means_" + pickle_file.split('.pkl')[0] + "_*.pkl"))
    if len(pickle_files) > 0:
        file = open(pickle_files[0], 'rb')
        means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
        file.close()
        means_full = {}
        for mean_key in means_dict.keys():
            means_full.update({mean_key+'_full':{}})
        Lifetimes_sys_full = {}
        Sep_maxs_full = []
        Sep_mins_full = []
        Initial_Seps_full = []
        Final_seps_full = []
        for pick_file in pickle_files:
            file = open(pick_file, 'rb')
            means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
            file.close()
            for time_key in means_dict['System_mean_times'].keys():
                for means_key in means_full.keys():
                    means_full[means_key].update({time_key:means_dict[means_key.split('_full')[0]][time_key]})
            for life_key in Lifetimes_sys.keys():
                Lifetimes_sys_full.update({life_key:Lifetimes_sys[life_key]})
            Sep_maxs_full = Sep_maxs_full + Sep_maxs
            Sep_mins_full = Sep_mins_full + Sep_mins
            Initial_Seps_full = Initial_Seps_full + Initial_Seps
            Final_seps_full = Final_seps_full + Final_seps
            os.remove(pick_file)
        
        Sep_maxs = Sep_maxs_full
        Sep_mins = Sep_mins_full
        Initial_Seps = Initial_Seps_full
        Final_seps = Final_seps_full
        Lifetimes_sys = Lifetimes_sys_full
        means_dict = {}
        for means_key in means_full.keys():
            means_dict.update({means_key.split('_full')[0]:means_full[means_key]})
        
        file = open('means_' + pickle_file + '.pkl', 'wb')
        pickle.dump((superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps),file)
        file.close()

print("gathered means")
sys.stdout.flush()
CW.Barrier()

#compile all the pickles
if rank == 0:
    file = open('means_' + pickle_file + '.pkl', 'rb')
    superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()

    #Create plot
    if args.plot_superplot == 'True':
        print("making superplot")
        plt.clf()
        fig, axs = plt.subplots(ncols=1, nrows=6, figsize=(8, 8))
        plt.subplots_adjust(hspace=0.0)
        gs = axs[3].get_gridspec()
        for ax in axs[3:]:
            ax.remove()
        axbig = fig.add_subplot(gs[3:])
        
        Times_adjusted = superplot_dict['Times'] - np.array(superplot_dict['Times'][0])
        
        axs[0].plot(Times_adjusted, superplot_dict['SFE'], label='SFE')
        axs[0].plot(Times_adjusted, superplot_dict['SFE_n'], label='SFE$_n$')
        axs[0].legend(loc='best')
        axs[0].set_ylabel('SFE')
        axs[0].set_ylim(bottom=0)
        axs[0].set_xlim([Times_adjusted[0], Times_adjusted[-1]])
        xticklabels = axs[0].get_xticklabels()
        plt.setp(xticklabels, visible=False)
        print("plotted SFE")
        
        axs[1].plot(Times_adjusted, superplot_dict['M_tot'], label='Total Mass')
        axs[1].plot(Times_adjusted, superplot_dict['M_tot_multi'], label='Mass in Multiple systems')
        axs[1].plot(Times_adjusted, superplot_dict['M_tot_vis'], label='Visible stars (based on Tobin)')
        axs[1].set_ylabel('M$_{tot}$ (M$_\odot$)')
        axs[1].set_ylim(bottom=0)
        axs[1].legend(loc='best')
        axs[1].set_xlim([Times_adjusted[0], Times_adjusted[-1]])
        xticklabels = axs[1].get_xticklabels()
        plt.setp(xticklabels, visible=False)
        print("plotted M_tot")
        
        axs[2].plot(Times_adjusted, superplot_dict['N_stars'], label='Total stars')
        axs[2].plot(Times_adjusted, superplot_dict['N_multi_stars'], label='Stars in Multiple systems')
        axs[2].plot(Times_adjusted, superplot_dict['N_vis_stars'], label='Visible stars (based on Tobin)')
        axs[2].set_ylabel('\# stars')
        axs[2].set_ylim(bottom=0)
        axs[2].legend(loc='best')
        axs[2].set_xlim([Times_adjusted[0], Times_adjusted[-1]])
        xticklabels = axs[2].get_xticklabels()
        plt.setp(xticklabels, visible=False)
        yticklabels = axs[2].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        plt.savefig('superplot.jpg', format='jpg', bbox_inches='tight')
        #plt.savefig('superplot.pdf', format='pdf', bbox_inches='tight')
        print("plotted N_star")
        
        window = 50000#year
        S_bins = np.logspace(1,4,13)
        for time_key in superplot_dict['System_times'].keys():
            print("plotting", time_key)
            axbig.semilogy((superplot_dict['System_times'][time_key]- np.array(superplot_dict['Times'][0])), superplot_dict['System_seps'][time_key], alpha=0.1, color='k')
            
            print("plotted system:", time_key)
        for bin_bound in S_bins:
            axbig.axhline(y=bin_bound, linewidth=0.5)
        axbig.set_ylabel('Separation (AU)')
        axbig.set_xlabel('Time ($yr$)')
        axbig.set_ylim([10, 10000])
        axbig.set_xlim([Times_adjusted[0], Times_adjusted[-1]])
        
        print("saving superplot")
        plt.savefig('superplot.jpg', format='jpg', bbox_inches='tight')
        print('Created superplot.jpg')
        #plt.savefig('superplot.pdf', format='pdf', bbox_inches='tight')
        
    d_sep = np.array(Initial_Seps) - np.array(Final_seps)
    shrinkage = 100*(d_sep/np.array(Initial_Seps))
    bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    shrink_hist, bins = np.histogram(shrinkage, bins=bins)
    shrink_hist = shrink_hist/np.sum(shrink_hist)
    bin_centers = np.linspace(-5, 95, 11)
    plt.clf()
    plt.bar(bin_centers, shrink_hist, width=10, edgecolor='black', alpha=0.5, label="Simulation")
    plt.xlabel('Orbital Shrinkage (\% of inital separation)')
    plt.ylabel('Fraction')
    plt.ylim(bottom=0.0)
    plt.savefig('shrinkage.jpg', format='jpg', bbox_inches='tight')
    print('created shrinkage.jpg')

sys.stdout.flush()
CW.Barrier()

print("Finished on rank", rank)
