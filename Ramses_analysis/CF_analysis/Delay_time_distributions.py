import pickle
import collections
import numpy as np
import matplotlib.pyplot as plt
import gc
import sys

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10
def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl"]


full_global_data_pickles = ['/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/global_reduced.pkl']

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(birth_con_pickles), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(birth_con_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
markers = ['o', '^', 's', 'p', 'h']
marker_labels = ['Binary', 'Trinary', 'Quadruple', 'Quintuple', 'Sextuple']

for pickle_it in range(len(birth_con_pickles)):
    scale_t_yr = 21728716.033625457

    file_open = open(full_global_data_pickles[pickle_it], 'rb')
    global_data = pickle.load(file_open)
    file_open.close()
    del file_open
    gc.collect()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    print("Finding formation inds", flush=True)
    sys.stdout.flush()
    ##print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    del global_data['x']
    del global_data['y']
    del global_data['z']
    del global_data['ux']
    del global_data['uy']
    del global_data['uz']
    gc.collect()
    
    sink_ids = np.arange(np.shape(global_data['m'].T)[0])
    
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    formation_inds = []
    SFE_arr = np.sum(global_data['m'], axis=1)
    
    for sink_id in sink_ids:
        try:
            new_ind = np.argwhere(global_data['m'].T[sink_id]>0)[0][0]
            global_data['m'] = global_data['m'][new_ind:]
            if len(formation_inds) == 0:
                formation_ind = new_ind
            else:
                formation_ind = formation_inds[-1]+new_ind
            formation_inds.append(formation_ind)
        except:
            sink_ids = sink_ids[:np.argwhere(sink_ids == sink_id)[0][0]]
            break
    gc.collect()

    print("Found formation inds", flush=True)
    sys.stdout.flush()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    formation_inds = np.array(formation_inds)
    formation_times = global_data['time'][formation_inds]*scale_t_yr
    
    del global_data['m']
    gc.collect()

    SFE = [[], [], [], [], []]
    T_delay = [[], [], [], [], []]
    Initial_seps = [[], [], [], [], []]
    
    file = open(birth_con_pickles[pickle_it], 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()
    
    for sink_key in Sink_birth_all.keys():
        if Sink_birth_all[sink_key][0] == False and Sink_birth_all[sink_key][1] == Sink_birth_all[sink_key][2]:
            n_stars = len(flatten(eval(Sink_birth_all[sink_key][2])))
            SFE[n_stars-1].append(Sink_birth_all[sink_key][-1])
            sfe_it = np.argmin(abs(SFE_arr - Sink_birth_all[sink_key][-1]))
            time_sys_form = global_data['time'][sfe_it]*scale_t_yr
            sink_form = formation_times[eval(sink_key)]
            delay = time_sys_form - sink_form
            T_delay[n_stars-1].append(delay)
            Initial_seps[n_stars-1].append(Sink_birth_all[sink_key][-4])
            
    for T_del_it in range(len(T_delay)):
        axs[pickle_it].scatter(np.array(SFE[T_del_it])*100, T_delay[T_del_it], marker=markers[T_del_it], label=marker_labels[T_del_it])
    if pickle_it == 0:
        axs[pickle_it].legend(ncol=3, loc='lower left')
    axs[pickle_it].set_yscale('log')
    
    if pickle_it == len(birth_con_pickles)-1:
        axs[pickle_it].set_xlabel('SFE (%)', fontsize=font_size)
    axs[pickle_it].set_ylabel('T$_{delay}$ (yr)', fontsize=font_size)
    axs[pickle_it].set_xlim([0, 5])

    axs[pickle_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs[pickle_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pickle_it].tick_params(axis='x', direction='in')
    axs[pickle_it].tick_params(axis='y', direction='in')
    
    axs[pickle_it].text((4.4), 2, subplot_titles[pickle_it], zorder=11, fontsize=font_size)
    
    plt.savefig('delay_vs_SFE.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(birth_con_pickles), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(birth_con_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
markers = ['o', '^', 's', 'p', 'h']
marker_labels = ['Binary', 'Trinary', 'Quadruple', 'Quintuple', 'Sextuple']

for pickle_it in range(len(birth_con_pickles)):
    scale_t_yr = 21728716.033625457

    file_open = open(full_global_data_pickles[pickle_it], 'rb')
    global_data = pickle.load(file_open)
    file_open.close()
    del file_open
    gc.collect()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    print("Finding formation inds", flush=True)
    sys.stdout.flush()
    ##print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    del global_data['x']
    del global_data['y']
    del global_data['z']
    del global_data['ux']
    del global_data['uy']
    del global_data['uz']
    gc.collect()
    
    sink_ids = np.arange(np.shape(global_data['m'].T)[0])
    
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)
    formation_inds = []
    SFE_arr = np.sum(global_data['m'], axis=1)
    
    for sink_id in sink_ids:
        try:
            new_ind = np.argwhere(global_data['m'].T[sink_id]>0)[0][0]
            global_data['m'] = global_data['m'][new_ind:]
            if len(formation_inds) == 0:
                formation_ind = new_ind
            else:
                formation_ind = formation_inds[-1]+new_ind
            formation_inds.append(formation_ind)
        except:
            sink_ids = sink_ids[:np.argwhere(sink_ids == sink_id)[0][0]]
            break
    gc.collect()

    print("Found formation inds", flush=True)
    sys.stdout.flush()
    #print("Memory_useage:", virtual_memory().percent, "on line", getframeinfo(currentframe()).lineno)

    formation_inds = np.array(formation_inds)
    formation_times = global_data['time'][formation_inds]*scale_t_yr
    
    del global_data['m']
    gc.collect()

    SFE = [[], [], [], [], []]
    T_delay = [[], [], [], [], []]
    Initial_seps = [[], [], [], [], []]
    
    file = open(birth_con_pickles[pickle_it], 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()
    
    for sink_key in Sink_birth_all.keys():
        if Sink_birth_all[sink_key][0] == False and Sink_birth_all[sink_key][1] == Sink_birth_all[sink_key][2]:
            n_stars = len(flatten(eval(Sink_birth_all[sink_key][2])))
            SFE[n_stars-1].append(Sink_birth_all[sink_key][-1])
            sfe_it = np.argmin(abs(SFE_arr - Sink_birth_all[sink_key][-1]))
            time_sys_form = global_data['time'][sfe_it]*scale_t_yr
            sink_form = formation_times[eval(sink_key)]
            delay = time_sys_form - sink_form
            T_delay[n_stars-1].append(delay)
            Initial_seps[n_stars-1].append(Sink_birth_all[sink_key][-4])
            
    for T_del_it in range(len(T_delay)):
        axs[pickle_it].scatter(Initial_seps[T_del_it], T_delay[T_del_it], marker=markers[T_del_it], label=marker_labels[T_del_it])
    if pickle_it == 0:
        axs[pickle_it].legend(ncol=3, loc='lower center')
    axs[pickle_it].set_yscale('log')
    
    if pickle_it == len(birth_con_pickles)-1:
        axs[pickle_it].set_xlabel('Initial Separation (au)', fontsize=font_size)
    axs[pickle_it].set_ylabel('T$_{delay}$ (yr)', fontsize=font_size)
    axs[pickle_it].set_xlim([0, 10000])

    axs[pickle_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs[pickle_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pickle_it].tick_params(axis='x', direction='in')
    axs[pickle_it].tick_params(axis='y', direction='in')
    
    axs[pickle_it].text((9000), 2, subplot_titles[pickle_it], zorder=11, fontsize=font_size)
    
plt.savefig('delay_vs_sep.pdf', bbox_inches='tight', pad_inches=0.02)
