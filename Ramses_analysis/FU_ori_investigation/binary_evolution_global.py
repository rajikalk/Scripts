import matplotlib.pyplot as plt
import numpy as np
import pickle
import yt
import matplotlib

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

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()

age = yt.YTArray(particle_data['time']).in_units('kyr')
mass_ratio = yt.YTArray(particle_data['mass']).T[1]/yt.YTArray(particle_data['mass']).T[0]

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(two_col_width, two_col_width), sharex=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

axs.flatten()[0].plot(age, yt.YTArray(particle_data['mass']).T[0], label='Primary')
axs.flatten()[0].plot(age, yt.YTArray(particle_data['mass']).T[1], label='Secondary')
axs.flatten()[0].set_xlim([age[0], age[-1]])
axs.flatten()[0].set_ylim([0, np.max(yt.YTArray(particle_data['mass']).T[0])])
axs.flatten()[0].legend()
axs.flatten()[0].set_ylabel('Mass (M$_\odot$)')
axs.flatten()[0].tick_params(axis='both', direction='in', top=True, right=True)

axs.flatten()[1].plot(age, mass_ratio)
axs.flatten()[1].set_ylim([0, 1])
axs.flatten()[1].set_ylabel('Mass ratio ($q$)')
axs.flatten()[1].tick_params(axis='both', direction='in', top=True, right=True)

axs.flatten()[2].semilogy(age, particle_data['separation'])
axs.flatten()[2].set_ylim([np.min(particle_data['separation']), np.max(particle_data['separation'])])
axs.flatten()[2].set_ylabel('Separation (AU)')
axs.flatten()[2].tick_params(axis='both', direction='in', top=True, right=True)

axs.flatten()[3].plot(age, particle_data['eccentricity'])
axs.flatten()[3].set_ylim([0, np.max(particle_data['eccentricity'])])
axs.flatten()[3].set_ylabel('Eccentricity')
axs.flatten()[3].set_xlabel('Time (kyr)')
axs.flatten()[3].tick_params(axis='both', direction='in', top=True, right=True)
plt.savefig('long_term_binary_evoluion.pdf', bbox_inches='tight', pad_inches=0.02)



