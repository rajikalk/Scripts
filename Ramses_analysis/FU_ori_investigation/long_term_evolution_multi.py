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

sink_inds = [10, 17, 45, 48, 51, 54, 56, 71, 73, 75, 85, 93, 103, 109, 118, 141, 150, 151, 154, 159, 168, 176, 177, 195, 221, 239, 258, 275, 262, 292]

#sink_inds = [10, 17, 45, 51, 56, 71, 73, 75, 85, 93, 103, 109, 118, 141, 150, 151, 154, 159, 168, 176, 177, 221, 239, 258, 262, 292]

labels = ['1', '2', '3', '4*', '5', '6*', '7', '8', '9+10', '11^', '12+13', '14', '17', '18', '19', '20^', '21', '22', '23', '24', '25', '26+27', '28', '29*', '30', '31', '32^', '33', '35*', '36']
#labels = ['1', '2', '3', '5', '7', '8', '9+10', '11^', '12+13', '14', '17', '18', '19', '20^', '21', '22', '23', '24', '25', '26+27', '28', '30', '31', '32^', '33', '36']

#
#plot_window = {'17' : [[19000, 55000], [56000, 75000]], '45' : [[8500, 27300], [30600, 75000]], '51' : [[16000, 30500], [30900, 75000]], '71' : [[7500, 38000]], '75' : [[5900, 6100], [14000, 17000], [18250, 18500], [19000, 24500], [27000, 75000]], '85' : [[2250, 75000]], '101' : [[3000, 75000]], '103' : [[1000, 2000], [21900, 75000]], '176' : [[36500, 39000], [45000, 46000], [48250, 49000]], '177' : [[32000, 75000]], '258' : [[6500, 12500], [13900, 75000]], '272' : [[10100, 29750], [41000, 75000]], '292' : [[3000, 4000], [5900, 75000]]}

plot_window = {'10':[[238983, 278330]],
               '17':[[31323, 88478]],
               '45':[[27994, 76169]],
               '48':[[1037941, 1078831]],
               '51':[[26015, 81016]],
               '54':[[112554, 112569]],
               '56':[[28998, 47214]],
               '71':[[13236, 13245]],
               '73':[[65139, 65471]],
               '75':[[30918, 52999]],
               '85':[[47797, 98532]],
               '93':[[158528, 875146]],
               '103':[[76647, 100965]],
               '109':[[861997, 1034221]],
               '118':[[163707, 186713]],
               '141':[[601800, 604687]],
               '150':[[679968, 748192]],
               '151':[[193236, 246112]],
               '154':[[326996, 520258]],
               '159':[[157623, 199899]],
               '168':[[72453, 88947]],
               '176':[[73943, 186419]],
               '177':[[149168, 197064]],
               '195':[[421958, 557305]],
               '221':[[37296, 199944]],
               '239':[[100755, 123564]],
               '258':[[14885, 114601]],
               '262':[[34266, 36958]],
               '275':[[23027, 23033]],
               '292':[[19204, 95581]],
}
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(two_col_width, 2*single_col_width), sharex=True)
plt.subplots_adjust(hspace=0.1)
smoothing_window = 1000
do_smoothing = False
plot_whole_binary = True
plot_windows = False

#axs.flatten()[2].set_xlim([0, 600000])
#axs.flatten()[2].set_xlim(left=0)
axs.flatten()[2].set_ylabel('$M_{cand.}/M_{clos.}$')
axs.flatten()[2].set_ylim([0, 1.25])
axs.flatten()[2].tick_params(axis='both', direction='in', top=True, right=True)
#axs.flatten()[0].axhline(y=1.6, color='k', ls='--', label="2r$_{\mathrm{acc}}$")
axs.flatten()[1].set_ylabel('Eccentricity')
axs.flatten()[1].set_ylim([0, 1.5])
axs.flatten()[1].tick_params(axis='both', direction='in', top=True, right=True)
axs.flatten()[0].set_ylabel('Semimajor Axis (AU)')
#axs.flatten()[2].set_xlabel('Time since candidate formation (yr)')
axs.flatten()[2].set_xlabel('Simulation Time (yr)')
axs.flatten()[0].set_ylim([5, 1e3])
axs.flatten()[0].axhline(y=16.6, color='k', ls='--', label="r$_{soft}$")
#axs.flatten()[0].axhline(y=16.6, color='k', ls='--', lw=0.5)
axs.flatten()[0].tick_params(axis='both', direction='in', top=True, right=True)

label_it = -1
for sink_ind in sink_inds:
    label_it = label_it+1
    pickle_file = '/scratch/ek9/rlk100/RAMSES/Analysis/Long_term_evolution_pickles/particle_data_'+str(sink_ind)+'.pkl'
    if os.path.isfile(pickle_file):
        if do_smoothing == True:
            if os.path.isfile('/scratch/ek9/rlk100/RAMSES/Analysis/Long_term_evolution_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl'):
                file_open = open('/scratch/ek9/rlk100/RAMSES/Analysis/Long_term_evolution_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl', 'rb')
                smooth_t, smooth_q, smooth_e, smooth_sep = pickle.load(file_open)
                file_open.close()
                
                if plot_windows == True:
                    plot_colour = None
                    for time_window in plot_window[str(sink_ind)]:
                        start_ind = np.argmin(abs(np.array(smooth_t)-time_window[0]))
                        end_ind = np.argmin(abs(np.array(smooth_t)-time_window[1]))
                    
                        if sink_ind == 45:
                            alpha_val = 1.0
                        else:
                            alpha_val = 0.25
                        
                        label = "Cand. " + labels[label_it]
                        if '*' in labels[label_it]:
                            linestyle = ':'
                        elif '^' in labels[label_it]:
                            linestyle='--'
                        else:
                            linestyle='-'
                        
                        if plot_colour == None:
                            p = axs.flatten()[2].plot(np.array(smooth_t)[start_ind:end_ind], np.array(smooth_q)[start_ind:end_ind], alpha=alpha_val, label=label, ls=line_style)
                            axs.flatten()[1].plot(np.array(smooth_t)[start_ind:end_ind], np.array(smooth_e)[start_ind:end_ind], alpha=alpha_val, ls=line_style)
                            axs.flatten()[0].semilogy(np.array(smooth_t)[start_ind:end_ind], np.array(smooth_sep)[start_ind:end_ind], alpha=alpha_val, ls=line_style)
                            plot_colour = p[-1].get_color()
                        else:
                            axs.flatten()[2].plot(np.array(smooth_t)[start_ind:end_ind], np.array(smooth_q)[start_ind:end_ind], alpha=alpha_val, color=plot_colour, ls=line_style)
                            axs.flatten()[1].plot(np.array(smooth_t)[start_ind:end_ind], np.array(smooth_e)[start_ind:end_ind], alpha=alpha_val, color=plot_colour, ls=line_style)
                            axs.flatten()[0].semilogy(np.array(smooth_t)[start_ind:end_ind], np.array(smooth_sep)[start_ind:end_ind], alpha=alpha_val, color=plot_colour, ls=line_style)
                else:
                    import pdb
                    pdb.set_trace()
                    if sink_ind == 45:
                        alpha_val = 1.0
                    else:
                        alpha_val = 0.25
                    
                    label = "Cand. " + labels[label_it]
                    if '*' in labels[label_it]:
                        line_style = ':'
                    elif '^' in labels[label_it]:
                        line_style='--'
                    else:
                        line_style='-'
                    
                    p = axs.flatten()[2].plot(np.array(smooth_t), np.array(smooth_q), alpha=alpha_val, ls=line_style)
                    axs.flatten()[1].plot(np.array(smooth_t), np.array(smooth_e), alpha=alpha_val, ls=line_style)
                    axs.flatten()[0].semilogy(np.array(smooth_t), np.array(smooth_sep), alpha=alpha_val, ls=line_style, label=label)
            else:
                print('reading ', pickle_file)
                file_open = open(pickle_file, 'rb')
                particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
                file_open.close()
                print('Finished reading pickle, calculating smoothed quantities')
                
                mass_ratio = yt.YTArray(particle_data['mass'])/yt.YTArray(particle_data['closest_mass'])
                smooth_t = []
                smooth_q = []
                smooth_e = []
                smooth_sep = []
                closes_inds = np.unique(particle_data['closest_sink'][:end_it], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'][:end_it], return_index=True)[1])]
                import pdb
                pdb.set_trace()
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
                    mean_sep = np.mean(particle_data['semimajor_axis'][start_it:end_it])
                    smooth_t.append(mean_t)
                    smooth_q.append(mean_q)
                    smooth_e.append(mean_e)
                    smooth_sep.append(mean_sep)
                
                label = "Cand. " + labels[label_it]
                if '*' in labels[label_it]:
                    linestyle = ':'
                elif '^' in labels[label_it]:
                    linestyle='--'
                else:
                    linestyle='-'
                axs.flatten()[2].plot(smooth_t, smooth_q, alpha=0.25, ls=linestyle)
                axs.flatten()[1].plot(smooth_t, smooth_e, alpha=0.25, ls=linestyle)
                axs.flatten()[0].semilogy(smooth_t, smooth_sep, alpha=0.25, label=label, ls=linestyle)
                
                print('updating pickle')
                file = open('/scratch/ek9/rlk100/RAMSES/Analysis/Long_term_evolution_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl', 'wb')
                pickle.dump((smooth_t, smooth_q, smooth_e, smooth_sep), file)
                file.close()
                print('Finished updating pickle', '/scratch/ek9/rlk100/RAMSES/Analysis/Long_term_evolution_pickles/smoothed_particle_data_'+str(sink_ind)+'.pkl')
        else:
            print('reading ', pickle_file)
            file_open = open(pickle_file, 'rb')
            particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
            file_open.close()
            print('Finished reading pickle')
            
            if plot_whole_binary == False:
                
            
                for time_window in plot_window[str(sink_ind)]:
                    start_t = time_window[0]
                    end_t = time_window[-1]
                    start_it = np.argmin(abs(particle_data['time'] - start_t))
                    end_it = np.argmin(abs(particle_data['time'] - end_t))
                    closes_inds = np.unique(particle_data['closest_sink'][start_it:end_it], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'][start_it:end_it], return_index=True)[1])]
                    
                    plot_colour = None
                    for comp_ind in closes_inds:
                        if comp_ind == closes_inds[0]:
                            label = "Cand. " + labels[label_it]
                        else:
                            label = None
                        if '*' in labels[label_it]:
                            linestyle = ':'
                        elif '^' in labels[label_it]:
                            linestyle='--'
                        else:
                            linestyle='-'
                            
                        curr_inds = np.argwhere(np.array(particle_data['closest_sink'][start_it:end_it]) == comp_ind).T[0]
                        diff_inds = np.setdiff1d(np.arange(len(particle_data['time'][start_it:end_it])), curr_inds)
                        smooth_t = np.copy(particle_data['time'][start_it:end_it]+sink_form_time.value)
                        smooth_t[diff_inds] = np.nan
                        mass_ratio = yt.YTArray(particle_data['mass'][start_it:end_it])/yt.YTArray(particle_data['closest_mass'][start_it:end_it])
                        smooth_q = np.copy(mass_ratio)
                        smooth_q[diff_inds] = np.nan
                        smooth_e = np.copy(particle_data['eccentricity'][start_it:end_it])
                        smooth_e[diff_inds] = np.nan
                        smooth_sep = np.copy(particle_data['semimajor_axis'][start_it:end_it])
                        smooth_sep[diff_inds] = np.nan
                        
                        if plot_colour == None:
                            p = axs.flatten()[2].plot(smooth_t, smooth_q, alpha=0.25, ls=linestyle)
                            axs.flatten()[1].plot(smooth_t, smooth_e, alpha=0.25, ls=linestyle)
                            axs.flatten()[0].semilogy(smooth_t, smooth_sep, alpha=0.25, label=label, ls=linestyle)
                            plot_colour = p[-1].get_color()
                        else:
                            axs.flatten()[2].plot(smooth_t, smooth_q, alpha=0.25, ls=linestyle, color=plot_colour)
                            axs.flatten()[1].plot(smooth_t, smooth_e, alpha=0.25, ls=linestyle, color=plot_colour)
                            axs.flatten()[0].semilogy(smooth_t, smooth_sep, alpha=0.25, label=label, ls=linestyle, color=plot_colour)
            else:
                closes_inds = np.unique(particle_data['closest_sink'], return_index=True)[0][np.argsort(np.unique(particle_data['closest_sink'], return_index=True)[1])]
                
                plot_colour = None
                for comp_ind in closes_inds:
                    if comp_ind == closes_inds[0]:
                        label = "Cand. " + labels[label_it]
                    else:
                        label = None
                    if '*' in labels[label_it]:
                        linestyle = ':'
                    elif '^' in labels[label_it]:
                        linestyle='--'
                    else:
                        linestyle='-'
                        
                    curr_inds = np.argwhere(np.array(particle_data['closest_sink']) == comp_ind).T[0]
                    diff_inds = np.setdiff1d(np.arange(len(particle_data['time'])), curr_inds)
                    smooth_t = np.copy(particle_data['time']+sink_form_time.value)
                    smooth_t[diff_inds] = np.nan
                    mass_ratio = yt.YTArray(particle_data['mass'])/yt.YTArray(particle_data['closest_mass'])
                    smooth_q = np.copy(mass_ratio)
                    smooth_q[diff_inds] = np.nan
                    smooth_e = np.copy(particle_data['eccentricity'])
                    smooth_e[diff_inds] = np.nan
                    smooth_sep = np.copy(particle_data['semimajor_axis'])
                    smooth_sep[diff_inds] = np.nan
                    alpha_arr = np.ones(np.shape(smooth_t))*0.25
                    time_window = plot_window[str(sink_ind)][0]
                    start_t = time_window[0]
                    end_t = time_window[-1]
                    start_it = np.argmin(abs(particle_data['time'] - start_t))
                    end_it = np.argmin(abs(particle_data['time'] - end_t))
                    alpha_arr[start_t:end_t] = 1
                    
                    if plot_colour == None:
                        p = axs.flatten()[2].plot(smooth_t, smooth_q, alpha=alpha_arr, ls=linestyle)
                        axs.flatten()[1].plot(smooth_t, smooth_e, alpha=alpha_arr, ls=linestyle)
                        axs.flatten()[0].semilogy(smooth_t, smooth_sep, alpha=alpha_arr, label=label, ls=linestyle)
                        plot_colour = p[-1].get_color()
                    else:
                        axs.flatten()[2].plot(smooth_t, smooth_q, alpha=alpha_arr, ls=linestyle, color=plot_colour)
                        axs.flatten()[1].plot(smooth_t, smooth_e, alpha=alpha_arr, ls=linestyle, color=plot_colour)
                        axs.flatten()[0].semilogy(smooth_t, smooth_sep, alpha=alpha_arr, label=label, ls=linestyle, color=plot_colour)
                
        axs.flatten()[0].legend(loc='upper center', bbox_to_anchor=(0.5, 2.1), ncol=4)
        plt.savefig("q_and_e_evol_all_candidates.pdf", bbox_inches='tight', pad_inches=0.02)

'''
    axs.flatten()[0].set_xlim([0, 75000])
    axs.flatten()[0].set_ylabel('$M_{cand}/M_{clos}$')
    axs.flatten()[0].set_ylim([0, 1.5])
    axs.flatten()[0].tick_params(axis='both', direction='in', top=True, right=True)
    axs.flatten()[1].set_ylabel('Eccentricity')
    axs.flatten()[1].set_ylim([0, 1.1])
    axs.flatten()[1].tick_params(axis='both', direction='in', top=True, right=True)
    axs.flatten()[2].set_ylabel('Separation (AU)')
    axs.flatten()[2].set_xlabel('Time (yr)')
    axs.flatten()[2].tick_params(axis='both', direction='in', top=True, right=True)
    axs.flatten()[0].legend(loc='center right', bbox_to_anchor=(1.3, 0.5), ncol=1)
    plt.savefig("q_and_e_evol_all_candidates.pdf", bbox_inches='tight', pad_inches=0.02)
    print('updated figure')
'''
