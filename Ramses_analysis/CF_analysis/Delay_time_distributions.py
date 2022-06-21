import pickle

birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Full_sink_data/Fast_analysis/sink_birth_all_delayed_core_frag_cleaned.pkl"]


SFE = []
T_delay = []
Initial_seps = []

for birth_con_pickle in birth_con_pickles:
    SFE.append([])
    T_delay.append([])
    Initial_seps.append([])
    
    file = open(birth_con_pickle, 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()
    
    for sink_key in Sink_birth_all.keys():
        if Sink_birth_all[sink_key][0] == False and Sink_birth_all[sink_key][1] == Sink_birth_all[sink_key][2]:
            SFE[-1].append(Sink_birth_all[sink_key][-1])
            T_delay[-1].append(Sink_birth_all[sink_key][-2])
            Initial_seps[-1].append(Sink_birth_all[sink_key][-4])

import matplotlib.pyplot as plt
subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.8*single_col_width), sharex=True)

for sim in range(len(subplot_titles)):
    plt.scatter(Initial_seps[sim], T_delay[sim], label=subplot_titles[sim])

plt.legend()

plt.xlabel('Initial Separation (au)', fontsize=font_size)
plt.ylabel('T$_{delay}$ (yr)', fontsize=font_size)

axs.tick_params(axis='both', which='major', labelsize=font_size, right=True)
axs.tick_params(axis='both', which='minor', labelsize=font_size)
axs.tick_params(axis='x', direction='in')
axs.tick_params(axis='y', direction='in')

plt.savefig('delay_vs_sep.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.8*single_col_width), sharex=True)

for sim in range(len(subplot_titles)):
    plt.scatter(SFE[sim], T_delay[sim], label=subplot_titles[sim])

plt.legend()
plt.xlabel('Initial Separation (au)', fontsize=font_size)
plt.ylabel('SFE', fontsize=font_size)

axs.tick_params(axis='both', which='major', labelsize=font_size, right=True)
axs.tick_params(axis='both', which='minor', labelsize=font_size)
axs.tick_params(axis='x', direction='in')
axs.tick_params(axis='y', direction='in')

plt.savefig('delay_vs_SFE.pdf', bbox_inches='tight', pad_inches=0.02)


            
