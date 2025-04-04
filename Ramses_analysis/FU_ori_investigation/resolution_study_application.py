import numpy as np
import pickle
import glob
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import sys
import os
import yt
import yt.units
from yt.units import g, s, cm, Lsun

#matplotlib.rcParams['pdf.fonttype'] = 42
#matplotlib.rcParams['ps.fonttype'] = 42

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
matplotlib.rcParams['text.latex.preamble'] = r"\usepackage{siunitx}" "\sisetup{detect-all}" r"\usepackage{helvet}" r"\usepackage{sansmath}" "\sansmath"               # <- tricky! -- gotta actually tell tex to use!

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("-pc", "--personal_computer", default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================
args = parse_inputs()
if args.personal_computer == 'True':
    pickle_files = ["/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L18.pkl", "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L19.pkl", "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl", "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L21.pkl"]
else:
    pickle_files = ["/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L18.pkl", "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L19.pkl", "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl", "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L21.pkl"]

length_unit = yt.YTQuantity(4.0,"pc")
r_acc = [np.round(length_unit.in_units('au')/(2**18)*4, decimals=2), np.round(length_unit.in_units('au')/(2**19)*4, decimals=2), np.round(length_unit.in_units('au')/(2**20)*4, decimals=2), np.round(length_unit.in_units('au')/(2**21)*4, decimals=2)]

label = ["$\Delta x=3.15$AU", "$\Delta x=1.57$AU", "$\Delta x=0.79$AU", "$\Delta x=0.39$AU"]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(single_col_width, single_col_width*1.7), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

cit = -1
t_end_yr = 9000 # None
for pick_file in pickle_files:
    print("read pickle", pick_file)
    cit = cit + 1
    file_open = open(pick_file, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")

    #find time bounds:
    t_start_yr = 3000
    if t_end_yr == None:
        t_end_yr = (np.array(particle_data['time']) - t_start_yr)[-1]
    t_start = np.argmin(abs(np.array(particle_data['time']) - t_start_yr))
    t_end = np.argmin(abs(np.array(particle_data['time']) - t_end_yr))

    f_acc = 0.5
    radius = yt.YTQuantity(2.0, 'rsun')
    #M_dot = accretion(sink_inds, global_ind)
    #M = yt.YTArray(global_data['m'][global_ind,sink_inds]*units['mass_unit'].in_units('msun'), 'Msun')
    m_dot = yt.YTArray(particle_data['mdot']).in_units('g/s')
    mass = yt.YTArray(particle_data['mass']).in_units('g')
    L_acc = f_acc * (mass * m_dot * yt.units.gravitational_constant_cgs)/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')
    
    for part in range(len(L_tot[t_start:t_end].T)):
        if part == 0:
            axs.flatten()[0].plot(particle_data['time'][t_start:t_end], np.array(particle_data['mass'][t_start:t_end]).T[part], label=label[pickle_files.index(pick_file)], color=proj_colours[cit], ls="--")
        else:
            axs.flatten()[0].plot(particle_data['time'][t_start:t_end], np.array(particle_data['mass'][t_start:t_end]).T[part], color=proj_colours[cit], ls="-")
    axs.flatten()[0].set_ylabel('Mass (M$_\odot$)', size=font_size)
    axs.flatten()[0].legend(loc=(0.5, 0.13), fontsize=font_size)
    axs.flatten()[0].set_ylim(bottom=0.03)
    axs.flatten()[0].tick_params(axis='both', direction='in', top=True, right=True)
    
    for part in range(len(L_tot[t_start:t_end].T)):
        if part == 0:
            axs.flatten()[1].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[part], label=label[pickle_files.index(pick_file)], color=proj_colours[cit], ls="--")
        else:
            axs.flatten()[1].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[part], color=proj_colours[cit], ls="-")
    axs.flatten()[1].set_ylabel('Accretion rate (M$_\odot$/yr)', size=font_size)
    #axs.flatten()[0].set_title('Sink no ' + str(sink_ind))
    axs.flatten()[1].set_ylim([1.e-9, 1.e-4])
    axs.flatten()[1].tick_params(axis='both', direction='in', top=True, right=True)
    print("plotted Accretion rate")
    
    '''
    for part in range(len(L_tot[t_start:t_end].T)):
        axs.flatten()[1].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[part], color=proj_colours[cit])
    axs.flatten()[1].set_ylabel('Accretion rate (Msun/yr)')
    axs.flatten()[1].set_ylim(bottom=1.e-9)
    '''
    axs.flatten()[2].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['separation'][t_start:t_end]), color=proj_colours[cit])
    axs.flatten()[2].set_ylabel('Separation (AU)', size=font_size)
    #axs.flatten()[1].axhline(y=r_acc[cit],c=proj_colours[cit], alpha=0.5)
    axs.flatten()[2].axhline(y=r_acc[cit]*2,c=proj_colours[cit], alpha=0.75, linestyle=':')
    axs.flatten()[2].set_xlabel('Time (yr)', size=font_size)
    axs.flatten()[2].set_xlim([t_start_yr,t_end_yr])
    axs.flatten()[2].tick_params(axis='both', direction='in', top=True, right=True)
    #axs.flatten()[1].set_ylim()
    print("plotted separation")
    
    '''
    for part in range(len(L_tot[t_start:t_end].T)):
        axs.flatten()[2].plot(particle_data['time'][t_start:t_end], np.array(particle_data['mass'][t_start:t_end]).T[part], color=proj_colours[cit])
    axs.flatten()[2].set_xlabel('Time (yr)')
    axs.flatten()[2].set_xlim([t_start_yr,t_end_yr])
    axs.flatten()[2].set_ylim(bottom=0)
    axs.flatten()[2].set_ylabel('Mass (Msun)')
    '''
    print("plotted", pick_file)
    plt.savefig('resolution_study_sink_'+str(sink_ind)+'.pdf', bbox_inches='tight', pad_inches=0.02)
