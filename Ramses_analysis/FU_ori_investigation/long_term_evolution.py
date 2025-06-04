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
fig, axs = plt.subplots(ncols=1, nrows=5, figsize=(two_col_width, 1.5*two_col_width), sharey=True)#, sharey=True)
#plt.subplots_adjust(wspace=0.0)
#plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[0].set_xlim([0, 15000])
axs.flatten()[0].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[0].tick_params(axis='both', direction='in', top=True)
ax0 = axs.flatten()[0].twinx()
ax0.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax0.set_ylabel('Separation (AU)')
ax0.set_ylim([5,1000])
ax0.tick_params(axis='both', direction='in', top=True)

axs.flatten()[1].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[1].set_xlim([15000, 30000])
axs.flatten()[1].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[1].tick_params(axis='both', direction='in', top=True)
ax1 = axs.flatten()[1].twinx()
ax1.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax1.set_ylabel('Separation (AU)')
ax1.set_ylim([5,1000])
ax1.tick_params(axis='both', direction='in', top=True)

axs.flatten()[2].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[2].set_xlim([30000, 45000])
axs.flatten()[2].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[2].tick_params(axis='both', direction='in', top=True)
ax2 = axs.flatten()[2].twinx()
ax2.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax2.set_ylabel('Separation (AU)')
ax2.set_ylim([5,1000])
ax2.tick_params(axis='both', direction='in', top=True)

axs.flatten()[3].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[3].set_xlim([45000, 60000])
axs.flatten()[3].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[3].tick_params(axis='both', direction='in', top=True)
ax3 = axs.flatten()[3].twinx()
ax3.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax3.set_ylabel('Separation (AU)')
ax3.set_ylim([5,1000])
ax3.tick_params(axis='both', direction='in', top=True)

axs.flatten()[4].semilogy(particle_data['time'], particle_data['lacc'])
axs.flatten()[4].set_xlim([60000, 75000])
axs.flatten()[4].set_ylabel("L$_{acc}$ (M$_\odot$/yr)")
axs.flatten()[4].tick_params(axis='both', direction='in', top=True)
ax4 = axs.flatten()[4].twinx()
ax4.semilogy(particle_data['time'], particle_data['separation'], color='k', ls="--", alpha=0.25)
ax4.set_ylabel('Separation (AU)')
ax4.set_ylim([5,1000])
ax4.tick_params(axis='both', direction='in', top=True)
axs.flatten()[4].set_xlabel("Time (yr)")
axs.flatten()[4].set_ylim([5.e-2, 2.e1])

plt.savefig('long_term_evolution.pdf', bbox_inches='tight', pad_inches=0.02)

