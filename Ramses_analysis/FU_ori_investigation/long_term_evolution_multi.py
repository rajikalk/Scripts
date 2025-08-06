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
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(two_col_width, 1.5*single_col_width), sharex=True)
plt.subplots_adjust(hspace=0.0)
smoothing_window = yt.YTArray(100, 'yr')

for sink_ind in sink_inds:
    pickle_file = '/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_'+str(sink_ind)+'.pkl'
    if os.path.isfile(pickle_file):
        if os.path.isfile('/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl'):
            file_open = open('/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl', 'rb')
            smooth_t, smooth_q, smooth_e = pickle.load(file_open)
            file_open.close()
            
            axs.flatten()[0].plot(smooth_t, smooth_q, alpha=0.25, label=str(sink_ind))
            axs.flatten()[1].plot(smooth_t, smooth_e, alpha=0.25)
        else:
            print('reading ', pickle_file)
            file_open = open(pickle_file, 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            print('Finished reading pickle, calculating smoothed quantities')
            
            mass_ratio = yt.YTArray(particle_data['mass']).T[0]/yt.YTArray(particle_data['mass']).T[1]
            smooth_t = []
            smooth_q = []
            smooth_e = []
            for it in range(len(particle_data['time'])):
                curr_time = particle_data['time'][it]
                start_time = particle_data['time'][it] - smoothing_window/2
                if start_time < 0:
                    start_time = 0
                start_it = np.argmin(abs(yt.YTArray(particle_data['time']) - start_time))
                
                end_time = particle_data['time'][it] + smoothing_window/2
                if end_time > particle_data['time'][-1]:
                    end_time = particle_data['time'][-1]
                end_it = np.argmin(abs(yt.YTArray(particle_data['time']) - end_time))
                
                mean_t = np.mean(particle_data['time'][start_it:end_it])
                mean_q = np.mean(mass_ratio[start_it:end_it])
                mean_e = np.mean(particle_data['eccentricity'][start_it:end_it])
                smooth_t.append(mean_t)
                smooth_q.append(mean_q)
                smooth_e.append(mean_e)
            
            axs.flatten()[0].plot(smooth_t, smooth_q, alpha=0.25, label=str(sink_ind))
            axs.flatten()[1].plot(smooth_t, smooth_e, alpha=0.25)
            
            print('updating pickle')
            file = open('/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl', 'wb')
            pickle.dump((smooth_t, smooth_q, smooth_e), file)
            file.close()
            print('Finished updating pickle', '/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl')

    axs.flatten()[0].set_xlim([0, 75000])
    axs.flatten()[0].set_ylabel('$q$ w.r.t closest sink')
    axs.flatten()[0].set_ylim([0, 1.5])
    axs.flatten()[0].tick_params(axis='both', direction='in', top=True, right=True)
    axs.flatten()[1].set_ylabel('Eccentricity')
    axs.flatten()[1].set_ylim([0, 1.1])
    axs.flatten()[1].set_xlabel('Time (yr)')
    axs.flatten()[1].tick_params(axis='both', direction='in', top=True, right=True)
    axs.flatten()[0].legend(loc='center right', bbox_to_anchor=(1.2, 1.0),
          ncol=2)
    plt.savefig("q_and_e_evol_all_candidates.pdf", bbox_inches='tight', pad_inches=0.02)
    print('updated figure')
    
