import pickle
import collections
import numpy as np
import matplotlib.pyplot as plt

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


plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(birth_con_pickles), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey='row')
iter_range = range(0, len(birth_con_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
markers = ['o', '^', 's', 'p', 'h']
marker_labels = ['Binary', 'Trinary', 'Quadruple', 'Quintuple', 'Sextuple']

for pickle_it in range(len(birth_con_pickles)):
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
            T_delay[n_stars-1].append(Sink_birth_all[sink_key][-2])
            Initial_seps[n_stars-1].append(Sink_birth_all[sink_key][-4])
            
    for T_del_it in range(len(T_delay)):
        axs[pickle_it].scatter(np.array(SFE[T_del_it])*100, T_delay[T_del_it], marker=markers[T_del_it], label=marker_labels[T_del_it])
    if pickle_it == 0:
        axs[pickle_it].legend()
    axs[pickle_it].set_yscale('log')
    
    if pickle_it == len(birth_con_pickles)-1:
        axs[pickle_it].set_xlabel('SFE (%)', fontsize=font_size)
    axs[pickle_it].set_ylabel('T$_{delay}$ (yr)', fontsize=font_size)
    axs[pickle_it].set_xlim([0, 10000])

    axs[pickle_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs[pickle_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pickle_it].tick_params(axis='x', direction='in')
    axs[pickle_it].tick_params(axis='y', direction='in')
    
plt.savefig('delay_vs_SFE.pdf', bbox_inches='tight', pad_inches=0.02)

'''
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.8*single_col_width), sharex=True)
ax = plt.gca()

for sim in range(len(subplot_titles)):
    plt.scatter(Initial_seps[sim], T_delay[sim], label=subplot_titles[sim], marker='.')

ax.set_yscale('log')
plt.legend()

plt.xlabel('Initial Separation (au)', fontsize=font_size)
plt.ylabel('T$_{delay}$ (yr)', fontsize=font_size)
plt.xlim([0, 10000])

axs.tick_params(axis='both', which='major', labelsize=font_size, right=True)
axs.tick_params(axis='both', which='minor', labelsize=font_size)
axs.tick_params(axis='x', direction='in')
axs.tick_params(axis='y', direction='in')

plt.savefig('delay_vs_sep.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.8*single_col_width), sharex=True)
ax = plt.gca()

for sim in range(len(subplot_titles)):
    plt.scatter(np.array(SFE[sim])*100, T_delay[sim], label=subplot_titles[sim], marker='.')
ax.set_yscale('log')
plt.legend()

plt.xlabel('SFE (%)', fontsize=font_size)
plt.ylabel('T$_{delay}$ (yr)', fontsize=font_size)
plt.xlim([0, 5])

axs.tick_params(axis='both', which='major', labelsize=font_size, right=True)
axs.tick_params(axis='both', which='minor', labelsize=font_size)
axs.tick_params(axis='x', direction='in')
axs.tick_params(axis='y', direction='in')

plt.savefig('delay_vs_SFE.pdf', bbox_inches='tight', pad_inches=0.02)
'''

            
