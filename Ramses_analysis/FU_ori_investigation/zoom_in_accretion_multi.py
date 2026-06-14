#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import matplotlib.colors as mcolors
import scipy.interpolate as interp
import pickle
import yt
import sys


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

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

args = parse_inputs()
save_dir = sys.argv[1]
'''
try:
    if sink_id == None:
        global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles_from_global_sim/particle_data_global.pkl"
    else:
        global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles_from_global_sim/particle_data_"+str(sink_id)+".pkl"
except:
    if sink_id == None:
        global_pickle = "/home/100/rlk100/rlk/RAMSES/Analysis/Sink_pickles/Low_res_pickles/particle_data_global.pkl"
    else:
'''
#global_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Long_term_evolution_pickles/particle_data_"+str(sink_id)+".pkl"

plot_sinks = [17, 45, 71, 75, 102]
Cand_labels = ['2', '3', '8', '11', '15+16']
particle_pickles = [['particle_data_17.pkl', 'particle_data_17_19.pkl'], ['array_particle_data_L18.pkl', 'array_particle_data_L19.pkl', 'array_particle_data_L20.pkl', 'array_particle_data_L21.pkl'], ['particle_data_71.pkl', 'particle_data_71_19.pkl', 'particle_data_71_20.pkl', 'particle_data_71_21.pkl'], ['particle_data_75.pkl', 'particle_data_75_19.pkl'], ['particle_data_102.pkl', 'particle_data_102_19.pkl', 'particle_data_102_20.pkl']]
length_unit = yt.YTQuantity(4.0,"pc")
r_acc = [np.round(length_unit.in_units('au')/(2**18)*4, decimals=2), np.round(length_unit.in_units('au')/(2**19)*4, decimals=2), np.round(length_unit.in_units('au')/(2**20)*4, decimals=2), np.round(length_unit.in_units('au')/(2**21)*4, decimals=2)]

res_label = ["$\Delta x=3.15$AU", "$\Delta x=1.57$AU", "$\Delta x=0.79$AU", "$\Delta x=0.39$AU"]
proj_colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(two_col_width, 1.5*two_col_width), sharey=True)#, sharey=True)

sink_it = -1
for plot_sink in plot_sinks:
    sink_it = sink_it + 1
    plot_pickles = particle_pickles[sink_it]
    Cand_label = 'Candidate '+ Cand_labels[sink_it]
    axs.flatten()[sink_it].set_title(Cand_label)
    pickle_it = -1
    for plot_pickle in plot_pickles:
        pickle_it = pickle_it+1
        if plot_sink != 45:
            file = open(plot_pickle, 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file)
            file.close()
            
            lns = []
            ln = axs.flatten()[sink_it].semilogy(particle_data['time'], particle_data['mdot'], label=res_label[pickle_it], color=proj_colours[pickle_it], ls='-')
            lns.append(ln)
            axs.flatten()[sink_it].semilogy(particle_data['time'], particle_data['closest_mdot'], color=proj_colours[pickle_it], ls=':')
            ax0 = axs.flatten()[sink_it].twinx()
            ln = ax0.semilogy(particle_data['time'], particle_data['separation'], color=proj_colours[pickle_it], ls="--", alpha=0.25, label="Separation")
            lns.append(ln)
            ax0.axhline(y=r_acc[pickle_it], color=proj_colours[pickle_it], ls='-.')
            lns.append(ln)
        else:
            file = open(plot_pickle, 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file)
            file.close()

            Cand_it = np.where(particle_data['particle_tag']==plot_sink)[0]
            Other_it = np.where(particle_data['particle_tag']!=plot_sink)[0]
            lns = []
            ln = axs.flatten()[sink_it].semilogy(particle_data['time'], particle_data['mdot'].T[Cand_it][0], label=res_label[pickle_it], color=proj_colours[pickle_it], ls='-')
            lns.append(ln)
            axs.flatten()[sink_it].semilogy(particle_data['time'], particle_data['mdot'].T[Other_it][0], color=proj_colours[pickle_it], ls=':')
            ax0 = axs.flatten()[sink_it].twinx()
            ln = ax0.semilogy(particle_data['time'], particle_data['separation'], color=proj_colours[pickle_it], ls="--", alpha=0.25, label="Separation")
            lns.append(ln)
            ax0.axhline(y=r_acc[pickle_it], color=proj_colours[pickle_it], ls='-.')
            lns.append(ln)

        
    #axs.flatten()[sink_it].set_xlim([time_bounds[0], time_bounds[1]])
    axs.flatten()[sink_it].set_ylabel("Accretion rate (M$_\odot$/yr)")
    axs.flatten()[sink_it].tick_params(axis='both', direction='in', top=True)
    ln_lab = lns[0]
    for ln_it in lns[1:]:
        ln_lab = ln_lab + ln_it
    labs = [l.get_label() for l in ln_lab]
    axs.flatten()[sink_it].legend(ln_lab, labs, loc='upper right', ncols=2, framealpha=0.9)
    ax0.set_ylabel('Separation (AU)')
    #ax0.set_ylim([5,1000])
    ax0.tick_params(axis='both', direction='in', top=True)
    print('plotted time panel', plot_sink)
    plt.savefig('zoom_in_multi.pdf', bbox_inches='tight', pad_inches=0.02)

 
