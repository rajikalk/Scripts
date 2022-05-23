import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pickle
import matplotlib.patches
import collections
import matplotlib
import matplotlib.ticker
#from mpi4py.MPI import COMM_WORLD as CW

matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

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

pickle_files = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/means_superplot.pkl"]

#pickle_files = ["/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G50/G50_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G100/G100_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G125/G125_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G150/G150_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G200/G200_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G400/G400_pathway.pkl"]

birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/sink_birth_all.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/sink_birth_all.pkl"]

plt.clf()
'''
if plot_booleans[rank][0] == True:
    fig, axs = plt.subplots(len(pickle_files), 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, len(pickle_files)*3), sharex='col')
    iter_range = range(0, len(pickle_files)*2, 2)
    args.figure_suffix = args.figure_suffix + '_hist'
else:
'''
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(two_col_width,page_height))
iter_range = range(0, len(pickle_files)*2)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.07)

#if plot_booleans[rank][1] == True:
#    args.figure_suffix = args.figure_suffix + '_thres'

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
GMC_mass_arr = [1500, 3000, 3750, 4500, 6000, 12000]
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

#================================================================================================
#plot fractions of core fragmentation and dynamical capture
plot_frag_capt_frac = False
if plot_frag_capt_frac:
    plt.clf()
    fig, axs = plt.subplots(1, 1, figsize=(6, 4))
    frag_fracs = []
    capt_fracs = []
    shrink_bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    for pick_it in range(len(pickle_files)):
        #try:
        file = open(pickle_files[pick_it], 'rb')
        superplot_dict, Sink_bound_birth, Sink_E_tot, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
        file.close()
        
        frag_fraction = np.sum(Sink_bound_birth)/len(Sink_bound_birth)
        capt_fraction = 1 - frag_fraction
        
        frag_fracs.append(frag_fraction)
        capt_fracs.append(capt_fraction)
        
        SFE_5_ind = np.argmin(abs(np.array(superplot_dict['SFE'])-0.05))
        SFE_5_time = superplot_dict['Times'][SFE_5_ind]
        valid_sink_ind = np.where(Sink_formation_times.value < SFE_5_time)[0]
        
        used_sink_inds = []
        Frag_shrinkage = []
        Capt_shrinkage = []
        for sys_id in superplot_dict[args.plot_key].keys():
            sink_inds = flatten(eval(sys_id))
            new_inds = list(set(sink_inds) - set(used_sink_inds))
            if len(new_inds) > 0 and Lifetimes_sys[sys_id] > args.sys_lifetime_threshold:
                newest_ind = np.max(new_inds)
                if len(sink_inds) == 2:
                    sep_init = superplot_dict[args.plot_key][sys_id][0][0]
                    non_nan_inds = np.argwhere(np.isnan(np.array(superplot_dict[args.plot_key][sys_id]).T[0])==False).T[0]
                    sys_life_time = np.array(superplot_dict['System_times'][sys_id])[non_nan_inds][-1] - np.array(superplot_dict['System_times'][sys_id])[non_nan_inds][0]
                    end_sep_time = np.array(superplot_dict['System_times'][sys_id] - np.array(superplot_dict['Times'][0]))[non_nan_inds][0] + sys_life_time - 10000
                    end_time_it = np.argmin(abs(np.array(superplot_dict['System_times'][sys_id]- np.array(superplot_dict['Times'][0]))[non_nan_inds] - end_sep_time))
                    sep_final = np.mean(np.array(superplot_dict[args.plot_key][sys_id]).T[0][non_nan_inds][end_time_it:])
                    shrinkage = 100 * (sep_init-sep_final)/sep_init
                    if Sink_bound_birth[newest_ind] == True:
                        Frag_shrinkage.append(shrinkage)
                    else:
                        Capt_shrinkage.append(shrinkage)
                else:
                    sys_comps = sys_id
                    new_sys_id = 1
                    sep_ind = 0
                    reduced = False
                    while reduced == False:
                        open_braket_ind = []
                        for char_it in range(len(sys_comps)):
                            if sys_comps[char_it] == '[':
                                open_braket_ind.append(char_it)
                            if sys_comps[char_it] == ']':
                                open_ind = open_braket_ind.pop()
                                sub_sys = eval(sys_comps[open_ind:char_it+1])
                                if newest_ind in sub_sys:
                                    sep_init = np.array(superplot_dict[args.plot_key][sys_id]).T[sep_ind][0]
                                    non_nan_inds = np.where(np.isnan(np.array(superplot_dict[args.plot_key][sys_id]).T[sep_ind]))[0]
                                    sys_life_time = np.array(superplot_dict['System_times'][sys_id])[non_nan_inds][-1] - np.array(superplot_dict['System_times'][sys_id])[non_nan_inds][0]
                                    end_sep_time = np.array(superplot_dict['System_times'][sys_id] - np.array(superplot_dict['Times'][0]))[non_nan_inds][0] + sys_life_time - 10000
                                    end_time_it = np.argmin(abs(np.array(superplot_dict['System_times'][sys_id]- np.array(superplot_dict['Times'][0]))[non_nan_inds] - end_sep_time))
                                    sep_final = np.mean(np.array(superplot_dict[args.plot_key][sys_id]).T[sep_ind][non_nan_inds][end_time_it:])
                                    shrinkage = 100 * (sep_init-sep_final)/sep_init
                                    if Sink_bound_birth[newest_ind] == True:
                                        Frag_shrinkage.append(shrinkage)
                                    else:
                                        Capt_shrinkage.append(shrinkage)
                                    import pdb
                                    pdb.set_trace()
                                else:
                                    sep_ind = sep_ind + 1
                                    replace_string = str(new_sys_id)
                                    new_sys_id = new_sys_id + 1
                                    sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                                    if '[' not in sys_comps:
                                        reduced = True
                                    break
            used_sink_inds = used_sink_inds + new_inds
            
        #except:
        #    print("pickle doesn't exist")
        
    plt.plot(GMC_mass_arr[:len(frag_fracs)], frag_fracs, label='Core Fragmention')
    plt.plot(GMC_mass_arr[:len(capt_fracs)], capt_fracs, label='Dynamical capture')
    plt.legend(loc='best')
    plt.xlabel('GMC Mass (M$_\odot$})')
    plt.ylabel('Fraction')
    plt.xlim([np.min(GMC_mass_arr), np.max(GMC_mass_arr)])
    plt.ylim([0,1])
    plt.savefig('Frag_capt_fractions.jpg')
    plt.savefig('Frag_capt_fractions.pdf')
    print('saved Frag_capt_fractions')

#Plot simulation SFE eovlution
plot_truncated_super_mult = True
if plot_truncated_super_mult == True:
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(two_col_width, page_height), sharex=True)#, sharey=True)
    iter_range = range(0, len(pickle_files))
    plt.subplots_adjust(wspace=0.0)
    plt.subplots_adjust(hspace=0.02)

    CF_hist = np.zeros((len(pickle_files),12)).tolist()
    
    Formation_pathway = []
    Not_plotted_sinks = [[],[],[],[],[],[]]
    if args.plot_key == 'System_seps':
        alpha = 0.025
    else:
        alpha = 0.1
    for pick_it in iter_range:
        pathway_counters = [0,0,0,0]
        core_frag_marker_pos = []
        delayed_core_frag_marker_pos = []
        dynamical_capture_marker_pos = []
    
        file_it = pick_it
        file = open(pickle_files[file_it], 'rb')
        superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
        file.close()
        
        file = open(birth_con_pickles[file_it], 'rb')
        Sink_birth_all = pickle.load(file)
        file.close()
        
        Initial_Seps = [[],[],[],[]]
        Initial_Seps_100000 = [[],[],[],[]]
        
        Start_times = []
        Sort_keys = []
        for time_key in superplot_dict['System_times'].keys():
            Start_times.append(superplot_dict['System_times'][time_key][0])
            Sort_keys.append(time_key)
        if np.min(np.array(Start_times)[1:] - np.array(Start_times)[:-1]) < 0:
            import pdb
            pdb.set_trace()

        SFE_5_ind = np.argmin(abs(np.array(superplot_dict['SFE'])-0.05))
        SFE_5_time = superplot_dict['Times'][SFE_5_ind]
        if pick_it == 0:
            if args.x_field == 'Time':
                G50_t_max = superplot_dict['Times'][SFE_5_ind] - superplot_dict['System_times'][list(superplot_dict['System_times'].keys())[0]][0]
            elif args.x_field == 'SFE':
                G50_t_max = 0.05
        new_sys_id = superplot_dict['N_stars'][-1]

        sim_start_time = np.nan
        plotted_sinks = []
        sub_sys_dict = {}
        counter = 0
        for time_key in superplot_dict['System_times'].keys():
            counter = counter + 1
            if counter > 7:
                import pdb
                pdb.set_trace()
            for sit in range(len(S_bins[1:])):
                bin_inst = np.argwhere((superplot_dict[args.plot_key][time_key]>S_bins[sit])&(superplot_dict[args.plot_key][time_key]<S_bins[sit+1]))
                CF_hist[file_it][sit] = CF_hist[file_it][sit] + np.shape(bin_inst)[0]
            key_inds = flatten(eval(time_key))
            if args.x_field == 'Time':
                if args.timescale == "first_sink":
                    first_sink = np.argmin(Sink_formation_times[key_inds])
                    Time_adjusted_formation = superplot_dict['System_times'][time_key] - Sink_formation_times[key_inds[first_sink]].value
                    SFE_5_time = superplot_dict['Times'][SFE_5_ind] - Sink_formation_times[key_inds[first_sink]].value
                else:
                    if np.isnan(sim_start_time):
                        sim_start_time = superplot_dict['System_times'][time_key][0]
                    Time_adjusted_formation = superplot_dict['System_times'][time_key] - np.array(sim_start_time)
                    SFE_5_time = superplot_dict['Times'][SFE_5_ind] - np.array(sim_start_time)
                t_max_ind = np.argmin(abs(Time_adjusted_formation - SFE_5_time))
                if t_max_ind > 0:
                    axs.flatten()[pick_it].semilogy(Time_adjusted_formation[:t_max_ind], superplot_dict[args.plot_key][time_key][:t_max_ind], alpha=alpha, color='k', rasterized=True)
            elif args.x_field == 'SFE':
                if superplot_dict['System_times'][time_key][0] < SFE_5_time:
                    sep_end_ind = np.argmin(abs(np.array(superplot_dict['System_times'][time_key]) - SFE_5_time))
                    Time_arr = superplot_dict['System_times'][time_key][:sep_end_ind+1]
                    SFE_arr = []
                    for time_val in Time_arr:
                        SFE_arr.append(superplot_dict['SFE'][superplot_dict['Times'].index(time_val)])
                    Sep_arr = np.array(superplot_dict[args.plot_key][time_key][:sep_end_ind+1])
                    if np.nanmin(Sep_arr) == 0.0:
                        zero_inds = np.where(Sep_arr == 0.0)
                        Sep_arr[zero_inds] = np.nan
                    if args.smooth_bool == 'True':
                        window = args.smooth_window_val #yr
                        smooth_SFE = [SFE_arr[0]]
                        smooth_Sep = [Sep_arr[0]]
                        for sit in range(1, len(SFE_arr)):
                            curr_time = Time_arr[sit]
                            start_time = curr_time - window
                            if start_time < Time_arr[0]:
                                start_time = Time_arr[0]
                            start_ind = np.argmin(abs(np.array(Time_arr) - start_time))
                            non_nan_inds = np.where(np.isnan(Sep_arr[start_ind:sit+1].T[0]) == False)[0]
                            SFE_smooth_val = np.mean(np.array(SFE_arr[start_ind:sit+1])[non_nan_inds])
                            sep_smooth_val = np.mean(np.array(Sep_arr[start_ind:sit+1])[non_nan_inds], axis=0)
                            smooth_SFE.append(SFE_smooth_val)
                            smooth_Sep.append(sep_smooth_val)
                        SFE_arr = smooth_SFE
                        Sep_arr = smooth_Sep
                    axs.flatten()[pick_it].semilogy(SFE_arr, Sep_arr, alpha=alpha, color='k', rasterized=True)
            if superplot_dict['System_times'][time_key][0] < SFE_5_time:
                sep_ind = 0
                if set(key_inds).issubset(set(plotted_sinks)):
                    pathway_counters[3] = pathway_counters[3] + 1
                    Initial_Seps[3].append(Sep_arr[0][sep_ind])
                    if Lifetimes_sys[time_key]>100000:
                        Initial_Seps_100000[3].append(Sep_arr[0][sep_ind])
                sys_comps = time_key
                reduced = False
                while reduced == False:
                    open_braket_ind = []
                    for char_it in range(len(sys_comps)):
                        if sys_comps[char_it] == '[':
                            open_braket_ind.append(char_it)
                        if sys_comps[char_it] == ']':
                            open_ind = open_braket_ind.pop()
                            sub_sys = eval(sys_comps[open_ind:char_it+1])
                            if np.mean(np.array(sub_sys)<superplot_dict['N_stars'][-1]) == 1:
                                if str(np.max(sub_sys)) in Sink_birth_all.keys():
                                    other_sys = np.min(sub_sys)
                                    if Sink_birth_all[str(np.max(sub_sys))][0] == True and str(other_sys) == str(Sink_birth_all[str(np.max(sub_sys))][1]):
                                        marker_color = 'b'
                                        marker_shape = 's'
                                        #elif Sink_birth_all[str(np.max(sub_sys))][1] in flatten(eval(Sink_birth_all[str(np.max(sub_sys))][2])):
                                        #elif other_sys in flatten(eval(str(Sink_birth_all[str(np.max(sub_sys))][1]))):
                                    elif str(other_sys) == str(Sink_birth_all[str(np.max(sub_sys))][1]) and np.sum(np.array(flatten(eval(Sink_birth_all[str(np.max(sub_sys))][2])))>np.max(sub_sys))==0:
                                        marker_color = 'm'
                                        marker_shape = '^'
                                        #elif str(other_sys) != str(Sink_birth_all[str(np.max(sub_sys))][1]):
                                    else:
                                        marker_color = 'r'
                                        marker_shape = 'o'
                                    #else:
                                    #    import pdb
                                    #    pdb.set_trace()
                                    #elif str(other_sys) == str(Sink_birth_all[str(np.max(sub_sys))][1]):
                                    #    marker_color = 'm'
                                    #    marker_shape = '^'
                                    
                                else:
                                    #print("sink", np.max(sub_sys), "Not found in birth conditions")
                                    Not_plotted_sinks[pick_it].append(np.max(sub_sys))
                                    marker_color = 'k'
                                    marker_shape = 'x'
                                if set(sub_sys).issubset(set(plotted_sinks)) == False:
                                    plotted_sinks = plotted_sinks + sub_sys
                                    #print('plotted sinks', sub_sys)
                                    if args.x_field == 'Time':
                                        axs.flatten()[pick_it].scatter(Time_adjusted_formation[:t_max][0], superplot_dict[args.plot_key][time_key][:t_max][0][sep_ind], color=marker_color, marker=marker_shape)
                                    elif args.x_field == 'SFE':
                                        #axs.flatten()[pick_it].scatter(SFE_arr[0], superplot_dict[args.plot_key][time_key][:sep_end_ind+1][0][sep_ind], color=marker_color, marker=marker_shape)
                                        if marker_color == 'b':
                                            print("Core_frag | The birth conditions for", np.max(sub_sys), "is", Sink_birth_all[str(np.max(sub_sys))], "| full system:", time_key, "sub_sys:", sub_sys)
                                            print("-------------------------------------------------------")
                                            core_frag_marker_pos.append([Sink_birth_all[str(np.max(sub_sys))][-1], Sink_birth_all[str(np.max(sub_sys))][3]])
                                            Initial_Seps[0].append(Sink_birth_all[str(np.max(sub_sys))][3])
                                            if Lifetimes_sys[time_key]>100000:
                                                Initial_Seps_100000[0].append(Sep_arr[0][sep_ind])
                                            pathway_counters[0] = pathway_counters[0] + 1
                                        elif marker_color == 'm':
                                            print("Delayed_core_frag | The birth conditions for", np.max(sub_sys), "is", Sink_birth_all[str(np.max(sub_sys))], "| full system:", time_key, "sub_sys:", sub_sys)
                                            print("-------------------------------------------------------")
                                            delayed_core_frag_marker_pos.append([Sink_birth_all[str(np.max(sub_sys))][-1], Sink_birth_all[str(np.max(sub_sys))][3]])
                                            Initial_Seps[1].append(Sink_birth_all[str(np.max(sub_sys))][3])
                                            if Lifetimes_sys[time_key]>100000:
                                                Initial_Seps_100000[1].append(Sep_arr[0][sep_ind])
                                            pathway_counters[1] = pathway_counters[1] + 1
                                        elif marker_color == 'r':
                                            print("Dynamical_capt | The birth conditions for", np.max(sub_sys), "is", Sink_birth_all[str(np.max(sub_sys))], "| full system:", time_key, "sub_sys:", sub_sys)
                                            print("-------------------------------------------------------")
                                            dynamical_capture_marker_pos.append([Sink_birth_all[str(np.max(sub_sys))][-1], Sink_birth_all[str(np.max(sub_sys))][3]])
                                            Initial_Seps[2].append(Sink_birth_all[str(np.max(sub_sys))][3])
                                            if Lifetimes_sys[time_key]>100000:
                                                Initial_Seps_100000[2].append(Sep_arr[0][sep_ind])
                                            pathway_counters[2] = pathway_counters[2] + 1
                            elif np.mean(np.array(sub_sys)<superplot_dict['N_stars'][-1]) < 1:
                                real_sink_inds = np.where(np.array(sub_sys)<superplot_dict['N_stars'][-1])[0]
                                real_sinks = np.array(sub_sys)[real_sink_inds]
                                if len(real_sinks) > 0:
                                    if str(np.max(real_sinks)) in Sink_birth_all.keys():
                                        other_sys = sub_sys_dict[str(list(set(sub_sys).difference(set(real_sinks)))[0])]
                                        other_sys_str = str(other_sys)
                                        while True in (np.array(flatten(other_sys)) > superplot_dict['N_stars'][-1]):
                                            greater_than_N_inds = np.argwhere(np.array(flatten(other_sys)) > superplot_dict['N_stars'][-1])[0]
                                            for greater_ind in greater_than_N_inds:
                                                other_split = other_sys_str.split(str(flatten(other_sys)[greater_ind]))
                                                insert_str = str(sub_sys_dict[str(flatten(other_sys)[greater_ind])])
                                                other_sys_str = other_split[0] + insert_str + other_split[1]
                                                other_sys = eval(other_sys_str)
                                        if Sink_birth_all[str(np.max(real_sinks))][0] == True and other_sys_str == str(Sink_birth_all[str(np.max(real_sinks))][2]):
                                            marker_color = 'b'
                                            marker_shape = 's'
                                        elif Sink_birth_all[str(np.max(real_sinks))][1] not in flatten(eval(Sink_birth_all[str(np.max(real_sinks))][2])):
                                            marker_color = 'r'
                                            marker_shape = 'o'
                                        else:
                                            if str(other_sys) == Sink_birth_all[str(np.max(real_sinks))][2]:
                                                marker_color = 'm'
                                                marker_shape = '^'
                                            else:
                                                marker_color = 'r'
                                                marker_shape = 'o'
                                        '''
                                        elif Sink_birth_all[str(np.max(real_sinks))][1] in flatten(eval(Sink_birth_all[str(np.max(real_sinks))][2])) and np.sum(np.array(flatten(eval(Sink_birth_all[str(np.max(real_sinks))][2])))>np.max(real_sinks))==0 and str(other_sys) == str(Sink_birth_all[str(np.max(real_sinks))][2]):
                                            marker_color = 'm'
                                            marker_shape = '^'
                                            #elif Sink_birth_all[str(np.max(sub_sys))][1] in flatten(eval(Sink_birth_all[str(np.max(sub_sys))][2])):
                                            #elif other_sys in flatten(eval(str(Sink_birth_all[str(np.max(sub_sys))][1]))):
                                            #elif other_sys_str != str(Sink_birth_all[str(np.max(real_sinks))][1]):
                                            #marker_color = 'r'
                                            #marker_shape = 'o'
                                        else:
                                            marker_color = 'r'
                                            marker_shape = 'o'
                                        '''
                                        #elif other_sys_str == str(Sink_birth_all[str(np.max(real_sinks))][1]):
                                        #    marker_color = 'm'
                                        #    marker_shape = '^'
                                    else:
                                        #print("sink", np.max(real_sinks), "Not found in birth conditions")
                                        Not_plotted_sinks[pick_it].append(np.max(real_sinks))
                                        marker_color = 'k'
                                        marker_shape = 'x'
                                    if set(real_sinks).issubset(set(plotted_sinks)) == False:
                                        plotted_sinks = plotted_sinks + real_sinks.tolist()
                                        #print('plotted sinks', real_sinks.tolist())
                                        if args.x_field == 'Time':
                                            axs.flatten()[pick_it].scatter(Time_adjusted_formation[:t_max][0], superplot_dict[args.plot_key][time_key][:t_max][0][sep_ind], color=marker_color, marker=marker_shape)
                                        elif args.x_field == 'SFE':
                                            #axs.flatten()[pick_it].scatter(SFE_arr[0], superplot_dict[args.plot_key][time_key][:sep_end_ind+1][0][sep_ind], color=marker_color, marker=marker_shape)
                                            if marker_color == 'b':
                                                print("Core_frag | The birth conditions for", np.max(real_sinks), "is", Sink_birth_all[str(np.max(real_sinks))], "| full system:", time_key, "sub_sys:", sub_sys)
                                                print("-------------------------------------------------------")
                                                core_frag_marker_pos.append([Sink_birth_all[str(np.max(real_sinks))][-1], Sink_birth_all[str(np.max(real_sinks))][3]])
                                                Initial_Seps[0].append(Sink_birth_all[str(np.max(real_sinks))][3])
                                                if Lifetimes_sys[time_key]>100000:
                                                    Initial_Seps_100000[0].append(Sep_arr[0][sep_ind])
                                                pathway_counters[0] = pathway_counters[0] + 1
                                            elif marker_color == 'm':
                                                print("Delayed_core_frag | The birth conditions for", np.max(real_sinks), "is", Sink_birth_all[str(np.max(real_sinks))], "| full system:", time_key, "sub_sys:", sub_sys)
                                                print("-------------------------------------------------------")
                                                delayed_core_frag_marker_pos.append([Sink_birth_all[str(np.max(real_sinks))][-1], Sink_birth_all[str(np.max(real_sinks))][3]])
                                                Initial_Seps[1].append(Sink_birth_all[str(np.max(real_sinks))][3])
                                                if Lifetimes_sys[time_key]>100000:
                                                    Initial_Seps_100000[1].append(Sep_arr[0][sep_ind])
                                                pathway_counters[1] = pathway_counters[1] + 1
                                            elif marker_color == 'r':
                                                print("Dynamical_capt | The birth conditions for", np.max(real_sinks), "is", Sink_birth_all[str(np.max(real_sinks))], "| full system:", time_key, "sub_sys:", sub_sys)
                                                print("-------------------------------------------------------")
                                                dynamical_capture_marker_pos.append([Sink_birth_all[str(np.max(real_sinks))][-1], Sink_birth_all[str(np.max(real_sinks))][3]])
                                                Initial_Seps[2].append(Sink_birth_all[str(np.max(real_sinks))][3])
                                                if Lifetimes_sys[time_key]>100000:
                                                    Initial_Seps_100000[2].append(Sep_arr[0][sep_ind])
                                                pathway_counters[2] = pathway_counters[2] + 1
                            sep_ind = sep_ind + 1
                            replace_string = str(new_sys_id)
                            sub_sys_dict.update({str(new_sys_id):sub_sys})
                            new_sys_id = new_sys_id + 1
                            sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                            if '[' not in sys_comps:
                                reduced = True
                            break
        
        axs.flatten()[pick_it].scatter(np.array(core_frag_marker_pos).T[0], np.array(core_frag_marker_pos).T[1], color='b', marker='s', s=7, zorder=10)
        axs.flatten()[pick_it].scatter(np.array(delayed_core_frag_marker_pos).T[0], np.array(delayed_core_frag_marker_pos).T[1], color='m', marker='^', s=7, zorder=10)
        axs.flatten()[pick_it].scatter(np.array(dynamical_capture_marker_pos).T[0], np.array(dynamical_capture_marker_pos).T[1], color='r', marker='o', s=7, zorder=10)

        print('Finished plotting separation evolution')
        Formation_pathway.append(pathway_counters)
        
        if args.plot_key == 'System_seps':
            for bin_bound in S_bins:
                axs.flatten()[pick_it].axhline(y=bin_bound, linewidth=0.5)

            axs.flatten()[pick_it].set_ylabel('Separation (AU)', size=args.text_font, labelpad=-1)
            axs.flatten()[pick_it].set_ylim([10, 10000])
            plt.minorticks_on()
            if pick_it == 0:
                axs.flatten()[pick_it].set_yticks([10, 100, 1000, 10000])
            else:
                axs.flatten()[pick_it].set_yticks([10, 100, 1000])
            y_minor = matplotlib.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
            axs.flatten()[pick_it].yaxis.set_minor_locator(y_minor)
            axs.flatten()[pick_it].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        else:
            axs.flatten()[pick_it].set_ylabel('Semi-major axis (AU)', size=args.text_font)
            axs.flatten()[pick_it].set_ylim([1.e1, 1.e9])
        axs.flatten()[pick_it].set_xlim([0, 0.05])
        axs.flatten()[pick_it].tick_params(which='both', direction='in')
        axs.flatten()[pick_it].tick_params(axis='both', which='major', labelsize=font_size)
        axs.flatten()[pick_it].tick_params(axis='both', which='minor', labelsize=font_size)
        
        #axs.flatten()[pick_it].add_patch(matplotlib.patches.Rectangle((0.01*G50_t_max, 12), 0.004, 14, facecolor="white", edgecolor="black", zorder=10))
        axs.flatten()[pick_it].add_patch(matplotlib.patches.Rectangle((0.01*G50_t_max, 11), 0.006, 18, facecolor="white", edgecolor="black", zorder=10))
        axs.flatten()[pick_it].text((0.0134*G50_t_max), 14, subplot_titles[pick_it], zorder=11, size=args.text_font)
        #axs.flatten()[pick_it].text((0.015*G50_t_max), 23, "N$_{sys}$="+str(np.sum(np.array(list(Lifetimes_sys.values())[:SFE_5_ind])>args.sys_lifetime_threshold)), zorder=12, size=args.text_font)
        
        if args.add_histograms == "True":
            CF_norm = np.array(CF_hist[file_it])/np.sum(CF_hist[file_it])
            axs.flatten()[pick_it+1].step([CF_norm[0]] + CF_norm.tolist(), np.log10(S_bins), linewidth=2)
            axs.flatten()[pick_it+1].set_ylim([1,4])
            yticklabels =axs.flatten()[pick_it+1].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        '''
        if pick_it > 0:
            yticklabels = axs.flatten()[pick_it].get_yticklabels()
            if str(yticklabels[-1]) == "Text(0, 10000, '$\\\\mathdefault{10^{4}}$')"
                plt.setp(yticklabels[-1], visible=False)
        '''
        if file_it != len(pickle_files)-1:
            xticklabels =axs.flatten()[pick_it].get_xticklabels()
            plt.setp(xticklabels, visible=False)
        
        plt.savefig('superplot_multi_truncated'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight', pad_inches=0.02, dpi=300)
        print('plotted separations for pickle', pickle_files[file_it])
        
        file = open('formation_pathway_'+str(GMC_mass_arr[pick_it])+'.pkl', 'wb')
        pickle.dump((pathway_counters, Initial_Seps, Initial_Seps_100000), file)
        file.close()
        
        #plt.savefig('superplot_multi'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight')
    
    if args.timescale == "first_sink":
        axs.flatten()[pick_it].set_xlabel('Star formation efficiency ($M_\star/M_{gas}$)', size=args.text_font)
    else:
        axs.flatten()[pick_it].set_xlabel('Star formation efficiency ($M_\star/M_{gas}$)', size=args.text_font)
    if args.add_histograms == "True":
        axs.flatten()[pick_it+1].set_xlabel('% of systems', size=args.text_font)
    plt.savefig('superplot_multi_truncated'+args.figure_suffix+'.jpg', format='jpg', bbox_inches='tight', pad_inches=0.02, dpi=300)
    #plt.savefig('superplot_multi_truncated'+args.figure_suffix+'.pdf', format='pdf', bbox_inches='tight', pad_inches=0.02)
    print('Created superplot_multi_truncated'+args.figure_suffix+'.jpg')

#plot quantities VS cloud mass at SFE=0.05
plt.clf()
fig, axs = plt.subplots(1, 1, figsize=(4, 3), sharex=True)
plt.plot(GMC_mass_arr, Bound_fractions, label='All stars')
plt.scatter(GMC_mass_arr, Bound_fractions)
plt.plot(GMC_mass_arr, Bound_fraction_Multi, label='Stars in multiples')
plt.scatter(GMC_mass_arr, Bound_fraction_Multi)
plt.xlabel('GMC Mass (M$_\odot$)')
plt.ylabel('Fraction of bound sinks at birth')
plt.xlim([GMC_mass_arr[0], GMC_mass_arr[-1]])
plt.legend(loc='best')
plt.savefig('bound_fraction_vs_GMC_mass.jpg', format='jpg', bbox_inches='tight')
plt.savefig('bound_fraction_vs_GMC_mass.pdf', format='pdf', bbox_inches='tight')
print('Saved bound_fraction_vs_GMC_mass.jpg')

plt.clf()
fig, axs = plt.subplots(1, 1, figsize=(4, 3), sharex=True)
plt.plot(GMC_mass_arr, N_Stars_5)
plt.scatter(GMC_mass_arr, N_Stars_5)
plt.xlabel('GMC Mass (M$_\odot$)')
plt.ylabel('\# stars')
plt.xlim([GMC_mass_arr[0], GMC_mass_arr[-1]])
plt.savefig('N_stars_vs_GMC_mass.jpg', format='jpg', bbox_inches='tight')
plt.savefig('N_stars_vs_GMC_mass.pdf', format='pdf', bbox_inches='tight')
print('Saved N_stars_vs_GMC_mass.jpg')
