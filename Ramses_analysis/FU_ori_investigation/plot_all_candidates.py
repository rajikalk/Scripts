import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt

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

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global_all_candidates.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_form_time = pickle.load(file_open)
file_open.close()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

for tag in particle_data['tag']:
    tag_it = particle_data['tag'].index(tag)
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(two_col_width, single_col_width))
    axs.semilogy(particle_data['age'][tag_it], particle_data['mdot'][tag_it])
    axs.set_xlim([np.min(particle_data['age'][tag_it]), np.max(particle_data['age'][tag_it])])
    axs.set_ylim([np.min(particle_data['mdot'][tag_it]), np.max(particle_data['mdot'][tag_it])])
    axs.set_xlabel('Age (yr)')
    axs.set_ylabel('Accretion rate (M$_\odot$/yr)')

    plt.savefig('Sink_'+str(tag)+'_accretion_evol.pdf', bbox_inches='tight', pad_inches=0.02)

