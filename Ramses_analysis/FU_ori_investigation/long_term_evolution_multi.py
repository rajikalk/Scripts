import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
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

sink_inds = [17, 45, 51, 71, 75, 85, 101, 103, 176, 177, 258, 272, 292]
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(two_col_width, single_col_width), sharex=True)
plt.subplots_adjust(hspace=0.0)

for sink_ind in sink_inds:
    pickle_file = '/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_'+str(sink_ind)+'_tmp.pkl'
    if os.path.isfile(pickle_file):
        print('reading ', pickle_file)
        file_open = open(pickle_file, 'rb')
        particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
        file_open.close()
        
        mass_ratio = yt.YTArray(particle_data['mass']).T[0]/yt.YTArray(particle_data['mass']).T[1]
        axs.flatten()[0].plot(particle_data['time'], mass_ratio)
        axs.flatten()[1].plot(particle_data['time'], particle_data['eccentricity'])

    axs.flatten()[0].set_xlim(left=0)
    axs.flatten()[0].set_ylabel('$q$ w.r.t closest sink')
    axs.flatten()[1].set_ylabel('Eccentricity')
    axs.flatten()[1].set_xlabel('Time (yr)')
    plt.savefig("q_and_e_evol_all_candidates.pdf", bbox_inches='tight', pad_inches=0.02)
    print('updated figure')
    
