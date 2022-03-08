import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pickle
import matplotlib.patches
import collections
import sys
from scipy import stats
#from mpi4py.MPI import COMM_WORLD as CW

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

pickle_file = sys.argv[1]
plot_gradient = True
#plot_key = sys.argv[2]
plot_keys = ['System_semimajor']#, 'System_ecc', 'System_energies']

for plot_key in plot_keys:
    savename = 'pathway_evolution_'+plot_key
    if plot_gradient:
        savename = savename + '_grad'
    plt.clf()
    fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(12, 9), sharey=True, sharex=True)
    if 'energies' in plot_key:
        axs.flatten()[0].set_yscale('symlog')
        axs.flatten()[1].set_yscale('symlog')
        axs.flatten()[2].set_yscale('symlog')
        axs.flatten()[3].set_yscale('symlog')
    plt.subplots_adjust(wspace=0.0)
    plt.subplots_adjust(hspace=0.07)

    file = open(pickle_file, 'rb')
    superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()

    SFE_5_ind = np.argmin(abs(np.array(superplot_dict['SFE'])-0.05))
    SFE_5_time = superplot_dict['Times'][SFE_5_ind]
     
    new_sys_id = superplot_dict['N_stars'][-1]

    sim_start_time = np.nan
    plotted_sinks = []
    set_start_time = True
    Initial_gradients_1000 = [[],[],[],[]]
    Initial_gradients = [[],[],[],[]]
    for time_key in superplot_dict['System_times'].keys():
        key_inds = flatten(eval(time_key))
        if set_start_time == True:
            start_time = superplot_dict['System_times'][time_key][0]
            set_start_time = False
        if superplot_dict['System_times'][time_key][0] < SFE_5_time:
            sep_end_ind = np.argmin(abs(np.array(superplot_dict['System_times'][time_key]) -SFE_5_time))
            Time_arr = superplot_dict['System_times'][time_key][:sep_end_ind+1]
            Time_arr = np.array(Time_arr) - Time_arr[0]

            sys_comps = time_key
            reduced = False
            sep_ind = 0
            while reduced == False:
                open_braket_ind = []
                for char_it in range(len(sys_comps)):
                    if sys_comps[char_it] == '[':
                        open_braket_ind.append(char_it)
                    if sys_comps[char_it] == ']':
                        open_ind = open_braket_ind.pop()
                        sub_sys = eval(sys_comps[open_ind:char_it+1])
                        real_sink_inds = np.where(np.array(sub_sys)<superplot_dict['N_stars'][-1])[0]
                        real_sinks = np.array(sub_sys)[real_sink_inds]
                        not_plotted_sinks = list(set(real_sinks).difference(set(plotted_sinks)))
                        if len(not_plotted_sinks) > 0:
                            birth_conditions = Sink_bound_birth[np.max(not_plotted_sinks)]
                            if birth_conditions[0] == True and birth_conditions[1] in key_inds:
                                axis_ind = 0
                                line_style = '-'
                                color = 'b'
                            elif birth_conditions[0] == False and birth_conditions[1] in key_inds:
                                axis_ind = 1
                                line_style = ':'
                                color='r'
                            else:
                                axis_ind = 2
                                line_style = '-'
                                color='k'
                        else:
                            axis_ind = 3
                            line_style = '-'
                            color='k'
                        Sep_arr = np.array(superplot_dict[plot_key][time_key]).T[sep_ind][:sep_end_ind+1]
                        dsep = Sep_arr[1:] - Sep_arr[:-1]
                        dtime = Time_arr[1:] - Time_arr[:-1]
                        grad = dsep/dtime
                        mean_x = (Time_arr[1:] + Time_arr[:-1])/2
                        try:
                            time_1000_ind = np.where(mean_x<1000)[0][-1] + 1
                            Initial_gradients_1000[axis_ind].append(grad[:time_1000_ind])
                        except:
                            print('system has not mean times < 1000yr')
                        if len(grad) > 0:
                            Initial_gradients[axis_ind].append(grad[0])
                        #lets trying smoothing over a 1000year window, but not centred

                        if plot_gradient == False:
                            if 'ecc' in plot_key:
                                axs.flatten()[axis_ind].plot(np.array(Time_arr), Sep_arr, alpha=0.2, color=color, rasterized=True, ls=line_style)
                            else:
                                axs.flatten()[axis_ind].semilogy(np.array(Time_arr), Sep_arr, alpha=0.2, color=color, rasterized=True, ls=line_style)
                            if 'energies' in plot_key:
                                axs.flatten()[0].set_yscale('symlog')
                                axs.flatten()[1].set_yscale('symlog')
                                axs.flatten()[2].set_yscale('symlog')
                        else:
                            axs.flatten()[axis_ind].loglog(mean_x, grad, alpha=0.2, color=color, rasterized=True, ls=line_style)
                            axs.flatten()[0].set_yscale('symlog')
                            axs.flatten()[1].set_yscale('symlog')
                            axs.flatten()[2].set_yscale('symlog')
                        plotted_sinks = plotted_sinks + not_plotted_sinks
                        #print('plotted sinks', not_plotted_sinks)
                        replace_string = str(new_sys_id)
                        new_sys_id = new_sys_id + 1
                        sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                        if '[' not in sys_comps:
                            reduced = True
                        break
            
        if plot_gradient == False:
            if 'semimajor' in plot_key or 'seps' in plot_key:
                axs.flatten()[0].set_ylabel('Separation (AU)')
                axs.flatten()[0].set_ylim([10, 10000])
                axs.flatten()[1].set_ylabel('Separation (AU)')
                axs.flatten()[1].set_ylim([10, 10000])
                axs.flatten()[2].set_ylabel('Separation (AU)')
                axs.flatten()[2].set_ylim([10, 10000])
                axs.flatten()[0].set_xlim(left=200)
        else:
            #axs.flatten()[0].set_ylim([-10, 10])
            axs.flatten()[0].set_ylabel('Gradient')
            #axs.flatten()[1].set_ylim([-10, 10])
            axs.flatten()[1].set_ylabel('Gradient')
            #axs.flatten()[2].set_ylim([-10, 10])
            axs.flatten()[2].set_ylabel('Gradient')
        axs.flatten()[1].set_xlabel('Time in simulation (yr)')
        plt.savefig(savename+'.png', format='png', bbox_inches='tight')
        
        #import pdb
        #pdb.set_trace()

    plt.savefig(savename+'.pdf', format='pdf', bbox_inches='tight')
    print('Created '+ savename+'.png')

grad_pickle = 'grad_pickle.pkl'
file = open(grad_pickle, 'wb')
pickle.dump((Initial_gradients, Initial_gradients_1000), file)
file.close()
print('saved_gradients')

grad_bins = np.concatenate((-1*np.logspace(5,-3,9), np.array([0, 1.e10])))
#grad_bins = np.concatenate((-1*np.logspace(4,-3,15), np.array([0, 1.e10])))
grad_hist_core, grad_bins = np.histogram(Initial_gradients[0], bins=grad_bins)
grad_hist_core_delayed, grad_bins = np.histogram(Initial_gradients[1], bins=grad_bins)
grad_hist_capt, grad_bins = np.histogram(Initial_gradients[2], bins=grad_bins)
grad_hist_misc, grad_bins = np.histogram(Initial_gradients[3], bins=grad_bins)

#calculate means
core_mean = []
for grit in range(len(Initial_gradients_1000[0])):
    core_mean.append(np.mean(Initial_gradients_1000[0][grit]))
grad_hist_core_mean, grad_bins = np.histogram(core_mean, bins=grad_bins)

core_delayed_mean = []
for grit in range(len(Initial_gradients_1000[1])):
    core_delayed_mean.append(np.mean(Initial_gradients_1000[1][grit]))
grad_hist_core_delayed_mean, grad_bins = np.histogram(core_delayed_mean, bins=grad_bins)

capt_mean = []
for grit in range(len(Initial_gradients_1000[2])):
    capt_mean.append(np.mean(Initial_gradients_1000[2][grit]))
grad_hist_capt_mean, grad_bins = np.histogram(capt_mean, bins=grad_bins)

misc_mean = []
for grit in range(len(Initial_gradients_1000[3])):
    misc_mean.append(np.mean(Initial_gradients_1000[3][grit]))
grad_hist_misc_mean, grad_bins = np.histogram(misc_mean, bins=grad_bins)

grad_bin_centers = (grad_bins[1:] + grad_bins[:-1])/2

bin_widths = grad_bins[1:] - grad_bins[:-1]
ticklabels = []
for bin_val in grad_bins:
    if int(np.sign(bin_val)) < 0:
        tick_str = str(int(np.sign(bin_val)))+'0$^{'+str(int(np.log10(abs(bin_val))))+'}$'
    elif int(np.sign(bin_val)) > 0:
        tick_str = '>0'
    else:
        tick_str = '0'
    ticklabels.append(r'{}'.format(tick_str))
ticklabels.append("")

grad_hist_core_norm = grad_hist_core/np.sum(grad_hist_core)
grad_hist_core_delayed_norm = grad_hist_core_delayed/np.sum(grad_hist_core_delayed)
grad_hist_capt_norm = grad_hist_capt/np.sum(grad_hist_capt)
grad_hist_misc_norm = grad_hist_misc/np.sum(grad_hist_misc)

grad_hist_core_norm = np.concatenate((grad_hist_core_norm, np.array([grad_hist_core_norm[-1]])))
grad_hist_core_delayed_norm = np.concatenate((grad_hist_core_delayed_norm, np.array([grad_hist_core_delayed_norm[-1]])))
grad_hist_capt_norm = np.concatenate((grad_hist_capt_norm, np.array([grad_hist_capt_norm[-1]])))
grad_hist_misc_norm = np.concatenate((grad_hist_misc_norm, np.array([grad_hist_misc_norm[-1]])))

x_range = np.arange(len(grad_hist_core)+1)

#Let's calculate a fit!
x_bins = np.log10(abs(grad_bin_centers[:-1]))[::-1]
y_core = grad_hist_core_norm[:-2][::-1]
plt.clf()
plt.semilogy(10**x_bins, y_core)
plt.xlabel('abs grad')
plt.ylabel('fraction')
plt.savefig('fit_test.png')
#samples_fit_log = scipy.stats.lognorm.pdf(x_vals, shape, loc=loc, scale=scale )

plt.clf()
fig, ax = plt.subplots(figsize=(6,6))
ax.step(x_range, grad_hist_core_norm, where='post', label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
ax.step(x_range, grad_hist_core_delayed_norm, where='post', label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
ax.step(x_range, grad_hist_capt_norm, where='post', label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
ax.step(x_range, grad_hist_misc_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
'''
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core/np.sum(grad_hist_core), label="Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='b')
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_delayed/np.sum(grad_hist_core_delayed), label="Delayed Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='purple')
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_capt/np.sum(grad_hist_capt), label="Dynamical Capture", width=1, color='None', linewidth=2, edgecolor='red')
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_misc/np.sum(grad_hist_misc), label="Other", width=1, color='None', linewidth=2, edgecolor='orange')
'''
ax.set_xlim([x_range[0], x_range[-1]])
ax.set_xticklabels(ticklabels)
ax.set_xlabel('Inspiral rate (au/yr)')
ax.set_ylabel('#')
ax.set_ylim(bottom=0)
ax.legend(loc='best')
plt.savefig('Initial_grad_hist.png')

grad_hist_core_mean_norm = grad_hist_core_mean/np.sum(grad_hist_core_mean)
grad_hist_core_delayed_mean_norm = grad_hist_core_delayed_mean/np.sum(grad_hist_core_delayed_mean)
grad_hist_capt_mean_norm = grad_hist_capt_mean/np.sum(grad_hist_capt_mean)
grad_hist_misc_mean_norm = grad_hist_misc_mean/np.sum(grad_hist_misc_mean)

grad_hist_core_mean_norm = np.concatenate((grad_hist_core_mean_norm, np.array([grad_hist_core_mean_norm[-1]])))
grad_hist_core_delayed_mean_norm = np.concatenate((grad_hist_core_delayed_mean_norm, np.array([grad_hist_core_delayed_mean_norm[-1]])))
grad_hist_capt_mean_norm = np.concatenate((grad_hist_capt_mean_norm, np.array([grad_hist_capt_mean_norm[-1]])))
grad_hist_misc_mean_norm = np.concatenate((grad_hist_misc_mean_norm, np.array([grad_hist_misc_mean_norm[-1]])))

plt.clf()
fig, ax = plt.subplots()
ax.step(x_range, grad_hist_core_mean_norm, where='post', label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
ax.step(x_range, grad_hist_core_delayed_mean_norm, where='post', label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
ax.step(x_range, grad_hist_capt_mean_norm, where='post', label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
ax.step(x_range, grad_hist_misc_mean_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
'''
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_mean/np.sum(grad_hist_core_mean), label="Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='b')
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_delayed_mean/np.sum(grad_hist_core_delayed_mean), label="Delayed Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='purple')
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_capt_mean/np.sum(grad_hist_capt_mean), label="Dynamical Capture", width=1, color='None', linewidth=2, edgecolor='red')
ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_misc_mean/np.sum(grad_hist_misc_mean), label="Other", width=1, color='None', linewidth=2, edgecolor='orange')
'''
ax.set_xlim([x_range[0], x_range[-1]])
ax.set_xticklabels(ticklabels)
ax.set_xlabel('Inspiral rate (au/yr)')
ax.set_ylabel('#')
ax.set_ylim(bottom=0)
ax.legend(loc='best')
plt.savefig('Mean_grad_hist.png')
