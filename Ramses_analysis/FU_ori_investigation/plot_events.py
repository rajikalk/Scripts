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
    pickle_file = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
else:
    pickle_file = "/lustre/astro/rlk/FU_ori_investigation/Sink_pickles/particle_data_L20.pkl"

length_unit = yt.YTQuantity(4.0,"pc")
r_acc = [np.round(length_unit.in_units('au')/(2**18)*4, decimals=2), np.round(length_unit.in_units('au')/(2**19)*4, decimals=2), np.round(length_unit.in_units('au')/(2**20)*4, decimals=2), np.round(length_unit.in_units('au')/(2**21)*4, decimals=2)]

label = ["3.15AU", "1.57AU", "0.79AU", "0.39AU"]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(single_col_width, single_col_width*2.5))#, sharey=True)
plt.subplots_adjust(wspace=0.0)
#plt.subplots_adjust(hspace=0.0)

time_bounds = [[3770, 4950],[5575, 5700], [6570, 6720], [7290, 7365], [7850, 7900]]


file_open = open(pickle_file, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()
print("finished reading in pickle")

e_it = -1
for t_bound in time_bounds:
    e_it = e_it + 1
    #find time bounds:
    t_start_yr = t_bound[0]
    t_end_yr = t_bound[1]
    t_start = np.argmin(abs(np.array(particle_data['time']) - t_start_yr))
    t_end = np.argmin(abs(np.array(particle_data['time']) - t_end_yr))
    
    #axs.flatten()[e_it].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[0])
    axs.flatten()[e_it].semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['mdot'][t_start:t_end]).T[1])
    axs.flatten()[e_it].set_ylabel('$\dot \mathrm{M}$ (M$_\odot$/yr)', size=font_size)
    axs.flatten()[e_it].set_xlim([t_start_yr,t_end_yr])
    axs.flatten()[e_it].tick_params(axis='both', direction='in')
    
    ax2 = axs.flatten()[e_it].twinx()
    ax2.semilogy(particle_data['time'][t_start:t_end], np.array(particle_data['separation'][t_start:t_end]), color='k', ls="--")
    ax2.set_ylabel('Separation (AU)')
    
axs.flatten()[e_it].set_ylabel('Time (yr)', size=font_size)
    
plt.savefig('suppression_events'+str(sink_ind)+'.pdf', bbox_inches='tight', pad_inches=0.02)
print("plot suppressino events")
