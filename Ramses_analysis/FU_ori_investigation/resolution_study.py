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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sink", "--sink_number", help="do you want to specific which sink to center on?", type=int, default=None)
    parser.add_argument("-update", "--update_pickle", help="Do you want to update the pickle?", type=str, default='True')
    parser.add_argument("-sim_dens_id", "--simulation_density_id", help="G50, G100, G200 or G400?", type=str, default="G100")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args
    
    
#================================================================================

pickle_files = ["/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L18.pkl", "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L19.pkl", "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl"]#, "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L21.pkl"]#, "/groups/astro/rlk/rlk/FU_ori_investigation/Accretion_evolution/Sink_45/Level_19/Restart/Level_20/Level_21/particle_data.pkl"]

label = ["Lvl=18", "Lvl=19", "Lvl=20", "Lvl=21"]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(two_col_width, single_col_width*2), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

cit = -1
t_end_yr = None #9000 # None
for pick_file in pickle_files:
    cit = cit + 1
    file_open = open(pick_file, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()

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
    L_acc = f_acc * (yt.units.G * mass * m_dot)/radius.in_units('cm')
    L_tot = L_acc.in_units('Lsun')

    for part in range(len(L_tot[t_start:t_end].T)):
        if part == 0:
            axs.flatten()[0].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[part], label=label[pickle_files.index(pick_file)], color=proj_colours[cit])
        else:
            axs.flatten()[0].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[part], color=proj_colours[cit])
    axs.flatten()[0].set_ylabel('Accretion rate (Msun/yr)')
    axs.flatten()[0].set_title('Sink no ' + str(sink_ind))
    axs.flatten()[0].set_ylim(bottom=1.e-9)
    axs.flatten()[0].legend()
    
    '''
    for part in range(len(L_tot[t_start:t_end].T)):
        axs.flatten()[1].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[part], color=proj_colours[cit])
    axs.flatten()[1].set_ylabel('Accretion rate (Msun/yr)')
    axs.flatten()[1].set_ylim(bottom=1.e-9)
    '''
    axs.flatten()[1].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['separation'][t_start:t_end]), color=proj_colours[cit])
    axs.flatten()[1].set_ylabel('Separation (au)')
    #axs.flatten()[1].set_ylim()
    

    for part in range(len(L_tot[t_start:t_end].T)):
        axs.flatten()[2].plot(particle_data['time'][t_start:t_end], np.array(particle_data['mass'][t_start:t_end]).T[part], color=proj_colours[cit])
    axs.flatten()[2].set_xlabel('Time (yr)')
    axs.flatten()[2].set_xlim([t_start_yr,t_end_yr])
    axs.flatten()[2].set_ylim(bottom=0)
    axs.flatten()[2].set_ylabel('Mass (Msun)')
    plt.savefig('resolution_study_sink_'+str(sink_ind)+'.pdf', bbox_inches='tight', pad_inches=0.02)
