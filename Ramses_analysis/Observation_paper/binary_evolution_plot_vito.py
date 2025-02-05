import csv
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import pickle

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

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-smooth_acc", "--smooth_accretion", help="do you want to smooth the accrtion", type=str, default='True')
    parser.add_argument("-window", "--smoothing_window", help="how big do you want the smoothing window to be, in term of phase", type=float, default=0.1)
    parser.add_argument("-savename", "--save_image_name", help="What would you like the plot to be saved as?", type=str, default="binary_evolution")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

end_peri = 13
args = parse_inputs()
labels = ['B1*', 'B1', 'B2', 'B3']
panel_tag = ['a)', 'b)', 'c)', 'd)']
linestyles = [':', '-', '-', '-']
colors = ['g', 'g', 'b', 'magenta']

pickle_files = ['/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_91/High_resolution/Remade_pickles/particle_data_neat.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_91/Remade_pickles/particle_data_neat.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_49/Remade_pickles/particle_data_neat.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Ramses/Sink_164/Remade_pickles/particle_data_neat.pkl']

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(single_col_width, two_col_width), sharex=True, sharey=False)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pit in range(len(pickle_files)):
    file_open = open(pickle_files[pit], 'rb')
    particle_data, line_counter, tag = pickle.load(file_open)
    file_open.close()
    
    if labels[pit] == 'B1':
        end_time_ind = np.argmin(abs(particle_data['time'].value - 110000))
    else:
        end_time_ind = len(particle_data['time']) - 1
    
    axs.flatten()[0].semilogy(particle_data['time'][:end_time_ind].in_units('kyr'), particle_data['separation'][:end_time_ind], label=labels[pit], linestyle=linestyles[pit], color=colors[pit], linewidth=1)
    if labels[pit] != 'B2':
        mass_ratio = particle_data['mass'][1]/particle_data['mass'][0]
    else:
        mass_ratio = particle_data['mass'][0]/particle_data['mass'][1]
        T_end = particle_data['time'][-1].in_units('kyr')
    total_mass = np.nansum(particle_data['mass'][:2], axis=0)
    axs.flatten()[1].plot(particle_data['time'][:end_time_ind].in_units('kyr'), total_mass[:end_time_ind], linestyle=linestyles[pit], color=colors[pit], label=labels[pit], linewidth=1)
    axs.flatten()[2].plot(particle_data['time'][:end_time_ind].in_units('kyr'), mass_ratio[:end_time_ind], linestyle=linestyles[pit], color=colors[pit], label=labels[pit], linewidth=1)
    axs.flatten()[3].plot(particle_data['time'][:end_time_ind].in_units('kyr'), particle_data['eccentricity'][:end_time_ind], linestyle=linestyles[pit], color=colors[pit], label=labels[pit], linewidth=1)
    
    #if labels[pit] == 'B1':
    #    import pdb
    #    pdb.set_trace()


axs.flatten()[0].set_xlim([0, T_end])
axs.flatten()[0].set_ylabel('Separation [AU]')
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[1].set_ylabel('Total Mass [M$_\odot$]')
axs.flatten()[1].set_ylim(bottom=0)
axs.flatten()[1].legend(loc='best')
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[2].set_ylabel('q [M$_s$/M$_p$]')
axs.flatten()[2].set_ylim([0, 1.1])
axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[3].set_ylabel('Eccentricity')
axs.flatten()[3].set_ylim([0.2, 1.1])
axs.flatten()[3].set_xlabel('Age [kyr]')
axs.flatten()[3].tick_params(axis='x', direction='in', top=True)
axs.flatten()[3].tick_params(axis='y', direction='in', right=True)
axs.flatten()[3].minorticks_on()
axs.flatten()[3].tick_params(which='both', direction='in', axis='both', right=True, top=True)
plt.savefig("system_evolution.pdf", bbox_inches='tight', pad_inches=0.02)
