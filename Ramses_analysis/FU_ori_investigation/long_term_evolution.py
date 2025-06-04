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

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(two_col_width, two_col_width), sharey=True)#, sharey=True)
#plt.subplots_adjust(wspace=0.0)
#plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[0].set_xlim([0, 25000])
axs.flatten()[0].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[1].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[1].set_xlim([25000, 50000])
axs.flatten()[1].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[2].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[2].set_xlim([50000, 75000])
axs.flatten()[2].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[2].set_ylim(bottom=1.e-3)

plt.savefig('long_term_evolution.pdf', bbox_inches='tight', pad_inches=0.02)

