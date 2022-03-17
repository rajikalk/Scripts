import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pickle
import matplotlib.patches
import collections
import matplotlib
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


#rank = CW.Get_rank()
#size = CW.Get_size()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

pickle_files = ["/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G50/G50_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G100/G100_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G125/G125_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G150/G150_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G200/G200_pathway.pkl", "/Users/reggie/Documents/Simulation_analysis/Pathway_evolution/G400/G400_pathway.pkl"]


subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(single_col_width, single_col_width*1.5), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()
    
    axs.flatten()[0].plot(np.array(superplot_dict['Times'])/1.e6, superplot_dict['M_tot'], label=subplot_titles[pick_it])
    axs.flatten()[1].plot(np.array(superplot_dict['Times'])/1.e6, superplot_dict['N_stars'])
    axs.flatten()[2].plot(np.array(superplot_dict['Times'])/1.e6, superplot_dict['SFE'])
    
    SFE_5_ind = np.argmin(abs(np.array(superplot_dict['SFE'])-0.05))
    
    axs.flatten()[0].scatter(superplot_dict['Times'][SFE_5_ind]/1.e6, superplot_dict['M_tot'][SFE_5_ind], marker='o')
    axs.flatten()[1].scatter(superplot_dict['Times'][SFE_5_ind]/1.e6, superplot_dict['N_stars'][SFE_5_ind], marker='o')
    axs.flatten()[2].scatter(superplot_dict['Times'][SFE_5_ind]/1.e6, superplot_dict['SFE'][SFE_5_ind], marker='o')
    
    print('plotted data from pickle', pickle_files[file_it])

axs.flatten()[2].set_xlabel('Time in Simulation (Myr)', fontsize=font_size)
axs.flatten()[0].set_ylabel('$M_{stars}$ (M$_\odot$)', labelpad=1, fontsize=font_size)
axs.flatten()[0].legend(loc='upper right', fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.3)
axs.flatten()[1].set_ylabel('$N_{stars}$', labelpad=1, fontsize=font_size)
axs.flatten()[2].set_ylabel('$SFE$ (%)', labelpad=1, fontsize=font_size)

axs.flatten()[0].tick_params(axis='both', which='major', labelsize=font_size)
axs.flatten()[0].tick_params(axis='both', which='minor', labelsize=font_size)
axs.flatten()[0].tick_params(axis='x', direction='in')
axs.flatten()[0].tick_params(axis='y', direction='in')
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in')
axs.flatten()[0].set_ylim(bottom=0)

axs.flatten()[1].tick_params(axis='both', which='major', labelsize=font_size)
axs.flatten()[1].tick_params(axis='both', which='minor', labelsize=font_size)
axs.flatten()[1].tick_params(axis='x', direction='in')
axs.flatten()[1].tick_params(axis='y', direction='in')
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in')
axs.flatten()[1].set_ylim(bottom=0)

axs.flatten()[2].tick_params(axis='both', which='major', labelsize=font_size)
axs.flatten()[2].tick_params(axis='both', which='minor', labelsize=font_size)
axs.flatten()[2].tick_params(axis='x', direction='in')
axs.flatten()[2].tick_params(axis='y', direction='in')
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in')
axs.flatten()[2].set_ylim(bottom=0)

plt.savefig('simulation_evolution.pdf', bbox_inches='tight', pad_inches=0.02)
    
