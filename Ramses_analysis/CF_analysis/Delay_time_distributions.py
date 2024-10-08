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

birth_con_pickles_low_cadence = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl"]

birth_con_pickles_high_cadence = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl"]


full_global_data_pickles = ['/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/global_reduced.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/global_reduced.pkl']

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(birth_con_pickles_low_cadence), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(birth_con_pickles_low_cadence))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
markers = ['o', '^', 's', 'p', 'h']
marker_labels = ['Binary', 'Trinary', 'Quadruple', 'Quintuple', 'Sextuple']

CF_delays = []
DCF_delays = []
DC_delays = []
CF_delays_diff = []
DCF_delays_diff = []
DC_delays_diff = []
Formation_pathway_high_cad = []
Formation_pathway_low_cad = []

for pickle_it in range(len(birth_con_pickles_low_cadence)):
    Formation_pathway_high_cad.append([0, 0, 0])
    Formation_pathway_low_cad.append([0, 0, 0])
    CF_delays.append([])
    DCF_delays.append([])
    DC_delays.append([])
    CF_delays_diff.append([])
    DCF_delays_diff.append([])
    DC_delays_diff.append([])
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
    
    file = open(birth_con_pickles_low_cadence[pickle_it], 'rb')
    Sink_birth_all_low_cad = pickle.load(file)
    file.close()
    '''
    for sink_key in Sink_birth_all_low_cad.keys():
        if Sink_birth_all_low_cad[sink_key][0] == False and Sink_birth_all_low_cad[sink_key][1] == Sink_birth_all_low_cad[sink_key][2]:
            n_stars = len(flatten(eval(Sink_birth_all_low_cad[sink_key][2])))
            SFE[n_stars-1].append(Sink_birth_all_low_cad[sink_key][-1])
            T_delay[n_stars-1].append(Sink_birth_all_low_cad[sink_key][-2])
            Initial_seps[n_stars-1].append(Sink_birth_all_low_cad[sink_key][-4])
    '''
    file = open(birth_con_pickles_high_cadence[pickle_it], 'rb')
    Sink_birth_all_high_cad = pickle.load(file)
    file.close()
    
    for sink_key in Sink_birth_all_high_cad.keys():
        if Sink_birth_all_high_cad[sink_key][0] == True:
            Formation_pathway_high_cad[-1][0] = Formation_pathway_high_cad[-1][0] + 1
        elif Sink_birth_all_high_cad[sink_key][0] == False and Sink_birth_all_high_cad[sink_key][1] == Sink_birth_all_high_cad[sink_key][2]:
            Formation_pathway_high_cad[-1][1] = Formation_pathway_high_cad[-1][1] + 1
            n_stars = len(flatten(eval(Sink_birth_all_high_cad[sink_key][2])))
            SFE[n_stars-1].append(Sink_birth_all_high_cad[sink_key][-1])
            sfe_it = np.argmin(abs(SFE_arr - Sink_birth_all_high_cad[sink_key][-1]))
            time_sys_form = global_data['time'][sfe_it]*scale_t_yr
            sink_form = formation_times[eval(sink_key)]
            delay = time_sys_form - sink_form
            T_delay[n_stars-1].append(delay)
            if Sink_birth_all_low_cad[sink_key][0] == True:
                CF_delays[-1].append(delay)
                CF_delays_diff.append(Sink_birth_all_low_cad[sink_key][-2] - delay)
            elif Sink_birth_all_low_cad[sink_key][0] == False and Sink_birth_all_low_cad[sink_key][1] == Sink_birth_all_low_cad[sink_key][2]:
                DCF_delays[-1].append(delay)
                DCF_delays_diff.append(Sink_birth_all_low_cad[sink_key][-2] - delay)
            else:
                DC_delays[-1].append(delay)
                DC_delays_diff.append(Sink_birth_all_low_cad[sink_key][-2] - delay)
            #T_delay[n_stars-1].append(Sink_birth_all[sink_key][-2])
            Initial_seps[n_stars-1].append(Sink_birth_all_high_cad[sink_key][-4])
        else:
            Formation_pathway_high_cad[-1][2] = Formation_pathway_high_cad[-1][2] + 1
    
        try:
            if Sink_birth_all_low_cad[sink_key][0] == True:
                Formation_pathway_low_cad[-1][0] = Formation_pathway_low_cad[-1][0] + 1
            elif Sink_birth_all_low_cad[sink_key][0] == False and Sink_birth_all_low_cad[sink_key][1] == Sink_birth_all_low_cad[sink_key][2]:
                Formation_pathway_low_cad[-1][1] = Formation_pathway_low_cad[-1][1] + 1
            else:
                Formation_pathway_low_cad[-1][2] = Formation_pathway_low_cad[-1][2] + 1
        except:
            import pdb
            pdb.set_trace()
            print("Sink", sink_key, "not found in low_cadence_data")
    
    for T_del_it in range(len(T_delay)):
        axs[pickle_it].scatter(np.array(SFE[T_del_it])*100, T_delay[T_del_it], marker=markers[T_del_it], label=marker_labels[T_del_it])
    if pickle_it == 0:
        axs[pickle_it].legend(ncol=3, loc='lower left')
    axs[pickle_it].set_yscale('log')
    
    if pickle_it == len(birth_con_pickles_low_cadence)-1:
        axs[pickle_it].set_xlabel('SFE (%)', fontsize=font_size)
    axs[pickle_it].set_ylabel('T$_{delay}$ (yr)', fontsize=font_size)
    axs[pickle_it].set_xlim([0, 5])

    axs[pickle_it].axhline(y=100)
    axs[pickle_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs[pickle_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pickle_it].tick_params(axis='x', direction='in')
    axs[pickle_it].tick_params(axis='y', direction='in')
    
    axs[pickle_it].text((4.4), 2, subplot_titles[pickle_it], zorder=11, fontsize=font_size)
    
    plt.savefig('delay_vs_SFE.pdf', bbox_inches='tight', pad_inches=0.02)


Core_frag_fracs_low = []
Delayed_core_frag_fracs_low = []
Dynamical_capt_fracs_low = []
Core_frag_fracs_high = []
Delayed_core_frag_fracs_high = []
Dynamical_capt_fracs_high = []
for sim in range(len(Formation_pathway_low_cad)):
    Total_sys_no = np.sum(Formation_pathway_low_cad[sim])
    Core_frag_frac = Formation_pathway_low_cad[sim][0]/Total_sys_no
    Delayed_core_frag_frac = Formation_pathway_low_cad[sim][1]/Total_sys_no
    Dynamical_capt_frac = Formation_pathway_low_cad[sim][2]/Total_sys_no
    import pdb
    pdb.set_trace()
    
    Core_frag_fracs_low.append(Core_frag_frac)
    Delayed_core_frag_fracs_low.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs_low.append(Dynamical_capt_frac)
    
    Total_sys_no = np.sum(Formation_pathway_high_cad[sim])
    Core_frag_frac = Formation_pathway_high_cad[sim][0]/Total_sys_no
    Delayed_core_frag_frac = Formation_pathway_high_cad[sim][1]/Total_sys_no
    Dynamical_capt_frac = Formation_pathway_high_cad[sim][2]/Total_sys_no
    
    Core_frag_fracs_high.append(Core_frag_frac)
    Delayed_core_frag_fracs_high.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs_high.append(Dynamical_capt_frac)

ind_low = np.arange(len(subplot_titles)*2)[::2]
ind_high = np.arange(len(subplot_titles)*2)[::2]+1

fig, ax = plt.subplots(1, 1, figsize=(5, 4))

p1 = plt.bar(ind_low, Core_frag_fracs_low, 0.95, color='b', linewidth=1, edgecolor='k')#, hatch='+'
p2 = plt.bar(ind_low, Delayed_core_frag_fracs_low, 0.95, bottom=Core_frag_fracs_low, color='m', linewidth=1, edgecolor='k')#, hatch='x'
p3 = plt.bar(ind_low, Dynamical_capt_fracs_low, 0.95, bottom=(np.array(Delayed_core_frag_fracs_low)+np.array(Core_frag_fracs_low)), color='r', linewidth=1, edgecolor='k')#, hatch='O'

p4 = plt.bar(ind_high, Core_frag_fracs_high, 0.95, color='b', linewidth=1, edgecolor='k')#, hatch='+'
p5 = plt.bar(ind_high, Delayed_core_frag_fracs_high, 0.95, bottom=Core_frag_fracs_high, color='m', linewidth=1, edgecolor='k')#, hatch='x'
p6 = plt.bar(ind_high, Dynamical_capt_fracs_high, 0.95, bottom=(np.array(Delayed_core_frag_fracs_high)+np.array(Core_frag_fracs_high)), color='r', linewidth=1, edgecolor='k')

#plt.xlim([-0.6, 5.6])
plt.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=font_size, right=True)
ax.tick_params(axis='both', which='minor', labelsize=font_size, left=True, right=True, top=False, bottom=False)
plt.xticks(ind_low, ("1500", "3000", "3750", "4500", "6000", "12000"))
ax.tick_params(which='both', direction='in')
plt.xlabel('Initial Gas Mass (M$_\odot$)', fontsize=font_size, labelpad=-0.5)

plt.legend((p3[0], p2[0], p1[0]), ('Dynamical capture', 'Delayed core frag.', 'Core fragmentation'), loc='upper right', fontsize=font_size)
plt.ylabel('Fraction', fontsize=font_size, labelpad=-0.5)
plt.ylim([0,1])
plt.savefig('formation_pathway_comp.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

import pdb
pdb.set_trace()

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(birth_con_pickles_low_cadence), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(birth_con_pickles_low_cadence))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
markers = ['o', '^', 's', 'p', 'h']
marker_labels = ['Binary', 'Trinary', 'Quadruple', 'Quintuple', 'Sextuple']

for pickle_it in range(len(birth_con_pickles_low_cadence)):
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
    
    file = open(birth_con_pickles_low_cadence[pickle_it], 'rb')
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
    
    if pickle_it == len(birth_con_pickles_low_cadence)-1:
        axs[pickle_it].set_xlabel('Initial Separation (au)', fontsize=font_size)
    axs[pickle_it].set_ylabel('T$_{delay}$ (yr)', fontsize=font_size)
    axs[pickle_it].set_xlim([0, 10000])

    axs[pickle_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs[pickle_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pickle_it].tick_params(axis='x', direction='in')
    axs[pickle_it].tick_params(axis='y', direction='in')
    
    axs[pickle_it].text((9000), 2, subplot_titles[pickle_it], zorder=11, fontsize=font_size)
    
    plt.savefig('delay_vs_sep.pdf', bbox_inches='tight', pad_inches=0.02)
