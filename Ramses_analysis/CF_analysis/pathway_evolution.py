import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pickle
import matplotlib.patches
import collections
import sys
from scipy import stats
from mpi4py.MPI import COMM_WORLD as CW

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

rank = CW.Get_rank()
size = CW.Get_size()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10

sys.stdout.flush()
CW.Barrier()

pickle_file = sys.argv[1]
true_birth_con_pickle = sys.argv[2]
plot_gradient = False
read_pickle = bool(sys.argv[3])
#plot_key = sys.argv[2]
plot_keys = ['System_semimajor']#, 'System_ecc', 'System_energies']

sys.stdout.flush()
CW.Barrier()

if read_pickle == True:
    for plot_key in plot_keys:
        savename = 'pathway_evolution_'+plot_key
        if plot_gradient:
            savename = savename + '_grad'
        '''
        plt.clf()
        fig, axs = plt.subplots(ncols=1, nrows=4, figsize=(12, 9), sharey=True, sharex=True)
        if 'energies' in plot_key:
            axs.flatten()[0].set_yscale('symlog')
            axs.flatten()[1].set_yscale('symlog')
            axs.flatten()[2].set_yscale('symlog')
            axs.flatten()[3].set_yscale('symlog')
        plt.subplots_adjust(wspace=0.0)
        plt.subplots_adjust(hspace=0.07)
        '''

        file = open(pickle_file, 'rb')
        #superplot_dict, Sink_bound_birth, Sink_formation_times = pickle.load(file)
        superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
        file.close()
        del Sink_formation_times
        print('Read means pickle on rank', rank)
        
        file = open(true_birth_con_pickle, 'rb')
        Sink_birth_all = pickle.load(file)
        file.close()

        sys.stdout.flush()
        CW.Barrier()

        SFE_5_ind = np.argmin(abs(np.array(superplot_dict['SFE'])-0.05))
        SFE_5_time = superplot_dict['Times'][SFE_5_ind]
         
        new_sys_id = superplot_dict['N_stars'][-1]
        first_new_sys_id = new_sys_id
        Sub_sys_indices = []

        sim_start_time = np.nan
        plotted_sinks = []
        set_start_time = True
        
        Grad_1e3 = []
        Grad_1e4 = []
        Initial_gradients = [[],[],[],[]]
        Initial_gradients_1000 = [[],[],[],[]]
        Initial_gradients_10000 = [[],[],[],[]]
        Initial_gradients_100000 = [[],[],[],[]]
        Used_sub_sys = []
        
        print('initialised gradient arrays')
        sys.stdout.flush()
        CW.Barrier()
        rit = -1
        for time_key in superplot_dict['System_times'].keys():
            rit = rit + 1
            if rit == size:
                rit = 0
            if rank == rit:
                print("processing system:", time_key, "on rank", rank)
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
                    sep_ind = -1
                    while reduced == False:
                        open_braket_ind = []
                        for char_it in range(len(sys_comps)):
                            if sys_comps[char_it] == '[':
                                open_braket_ind.append(char_it)
                            if sys_comps[char_it] == ']':
                                open_ind = open_braket_ind.pop()
                                sub_sys = eval(sys_comps[open_ind:char_it+1])
                                if str(sub_sys) not in Used_sub_sys:
                                    real_sink_inds = np.where(np.array(sub_sys)<superplot_dict['N_stars'][-1])[0]
                                    real_sinks = np.array(sub_sys)[real_sink_inds]
                                    not_plotted_sinks = list(set(real_sinks).difference(set(plotted_sinks)))
                                    if len(not_plotted_sinks) > 0:
                                        other_sys = np.min(sub_sys)
                                        if Sink_birth_all[str(np.max(not_plotted_sinks))][0] == True and str(other_sys) == str(Sink_birth_all[str(np.max(not_plotted_sinks))][1]):
                                            print("Core_frag | The birth conditions for", np.max(not_plotted_sinks), "is", Sink_birth_all[str(np.max(not_plotted_sinks))], "| full system:", time_key, "sub_sys:", sub_sys)
                                            print("-------------------------------------------------------")
                                            #if birth_conditions[0] == True and birth_conditions[1] in key_inds:
                                            axis_ind = 0
                                            line_style = '-'
                                            color = 'b'
                                        elif str(other_sys) == str(Sink_birth_all[str(np.max(not_plotted_sinks))][1]) and np.sum(np.array(flatten(eval(Sink_birth_all[str(np.max(not_plotted_sinks))][2])))>np.max(not_plotted_sinks))==0:
                                            print("Delayed_core_frag | The birth conditions for", np.max(not_plotted_sinks), "is", Sink_birth_all[str(np.max(not_plotted_sinks))], "| full system:", time_key, "sub_sys:", sub_sys)
                                            print("-------------------------------------------------------")
                                            #elif birth_conditions[0] == False and birth_conditions[1] in key_inds:
                                            axis_ind = 1
                                            line_style = ':'
                                            color='r'
                                        else:
                                            print("Dynamical_capt | The birth conditions for", np.max(not_plotted_sinks), "is", Sink_birth_all[str(np.max(not_plotted_sinks))], "| full system:", time_key, "sub_sys:", sub_sys)
                                            print("-------------------------------------------------------")
                                            axis_ind = 2
                                            line_style = '-'
                                            color='k'
                                    else:
                                        axis_ind = 3
                                        line_style = '-'
                                        color='k'
                                    if axis_ind == 0:
                                        form_path = 'Core_frag'
                                    elif axis_ind == 1:
                                        form_path = 'Delayed_core_frag'
                                    elif axis_ind == 2:
                                        form_path = 'Dynamical_capt'
                                    elif axis_ind == 3:
                                        form_path = 'Other'
                                    sep_ind = sep_ind + 1
                                    Sep_arr = np.array(superplot_dict[plot_key][time_key]).T[sep_ind][:sep_end_ind+1]
                                    
                                    Sep_arr_true = np.array(superplot_dict['System_seps'][time_key]).T[sep_ind]
                                    Time_arr_full = np.array(superplot_dict['System_times'][time_key]) - superplot_dict['System_times'][time_key][0]
                                    peri_inds = np.where((Sep_arr_true[1:-1] < Sep_arr_true[:-2]) & (Sep_arr_true[1:-1] < Sep_arr_true[2:]))[0]
                                    
                                    if len(peri_inds) > 2 and len(Sep_arr_true)>peri_inds[0]:
                                        plt.clf()
                                        plt.figure(figsize=(15, 3))
                                        plt.semilogy(Time_arr_full, Sep_arr_true)
                                        plt.scatter(Time_arr_full[1:-1][peri_inds], Sep_arr_true[1:-1][peri_inds])
                                        plt.xlim([0, 100000])
                                        plt.savefig('Peri_check_'+str(sub_sys).replace(' ', '')+'.png')
                                        
                                        
                                        initial_a = Sep_arr[1:][peri_inds[0]]
                                        initial_t = Time_arr[1:][peri_inds[0]]
                                        end_t = initial_t + 10000
                                        end_t_ind = np.argmin(abs(Time_arr - end_t))
                                        end_a = Sep_arr[end_t_ind]
                                        end_t_data = Time_arr[end_t_ind]
                                        mean_grad = (end_a-initial_a)/(end_t_data-initial_t)
                                        Initial_gradients_10000[axis_ind].append([mean_grad])
 
                                    '''
                                    import pdb
                                    pdb.set_trace()
                                        
                                    dsep = Sep_arr[1:] - Sep_arr[:-1]
                                    dtime = Time_arr[1:] - Time_arr[:-1]
                                    grad = dsep/dtime
                                    mean_x = (Time_arr[1:] + Time_arr[:-1])/2
                                    try:
                                        sub_1000_inds = np.where(Time_arr<1000)[0]
                                        dt = Time_arr[sub_1000_inds[-1]] - Time_arr[sub_1000_inds[0]]
                                        ds = Sep_arr[sub_1000_inds[-1]] - Sep_arr[sub_1000_inds[0]]
                                        mean_grad = ds/dt
                                        if mean_grad in np.array(Initial_gradients_1000[axis_ind]):
                                            import pdb
                                            pdb.set_trace()
                                        Initial_gradients_1000[axis_ind].append([mean_grad])
                                    except:
                                        print('system has not mean times < 1000yr')
                                    try:
                                        sub_10000_inds = np.where(Time_arr<10000)[0]
                                        non_nan_ind = np.where(np.isnan(Sep_arr)==False)[0]
                                        if len(sub_10000_inds)>9 and Time_arr[non_nan_ind][-1]>10000 and np.max(sub_10000_inds[1:] - sub_10000_inds[:-1])==1:
                                            dt = Time_arr[sub_10000_inds[-1]] - Time_arr[sub_10000_inds[0]]
                                            ds = Sep_arr[sub_10000_inds[-1]] - Sep_arr[sub_10000_inds[0]]
                                            mean_grad = ds/dt
                                            if mean_grad in np.array(Initial_gradients_10000[axis_ind]):
                                                import pdb
                                                pdb.set_trace()
                                            Initial_gradients_10000[axis_ind].append([mean_grad])
                                            if mean_grad < -1e4:
                                                Grad_1e4.append(time_key)
                                            if mean_grad < -1e3:
                                                Grad_1e3.append(time_key)
                                            if mean_grad < -1e2:
                                                plt.clf()
                                                fig, axs = plt.subplots(ncols=1, nrows=3, figsize=(two_col_width, single_col_width), sharex=True)
                                                plt.subplots_adjust(wspace=0.0)
                                                plt.subplots_adjust(hspace=0.0)
                                                axs[0].set_title('System:'+time_key+', form_path:'+form_path+', mean_grad:'+str(mean_grad))
                                                axs[0].semilogy((np.array(superplot_dict['System_times'][time_key]) - superplot_dict['System_times'][time_key][0]), np.array(superplot_dict['System_semimajor'][time_key]).T[sep_ind], label='Semimajor axis')
                                                axs[1].semilogy((np.array(superplot_dict['System_times'][time_key]) - superplot_dict['System_times'][time_key][0]), np.array(superplot_dict['System_seps'][time_key]).T[sep_ind], label='Separation')
                                                axs[1].set_ylim([10, 10000])
                                                axs[2].plot((np.array(superplot_dict['System_times'][time_key]) - superplot_dict['System_times'][time_key][0]), np.array(superplot_dict['System_ecc'][time_key]).T[sep_ind], label='Eccentricity')
                                                axs[2].set_ylim([0.0, 1.1])
                                                axs[0].set_ylabel('Semimajor Axis (au)')
                                                axs[1].set_ylabel('Separation (au)')
                                                axs[2].set_ylabel('Eccentricity')
                                                axs[2].set_xlabel('Time (yr)')
                                                plt.savefig('System:'+time_key+'.png')
                                        else:
                                            print('Not enough points to suggest system stays bound for first  10000yr')
                                    except:
                                        print('system has no mean times < 10000yr')
                                    try:
                                        sub_100000_inds = np.where(Time_arr<100000)[0]
                                        non_nan_ind = np.where(np.isnan(Sep_arr)==False)[0]
                                        if len(sub_10000_inds)>9 and Time_arr[non_nan_ind][-1]>100000 and np.max(sub_100000_inds[1:] - sub_100000_inds[:-1])==1:
                                            dt = Time_arr[sub_100000_inds[-1]] - Time_arr[sub_100000_inds[0]]
                                            ds = Sep_arr[sub_100000_inds[-1]] - Sep_arr[sub_100000_inds[0]]
                                            mean_grad = ds/dt
                                            if mean_grad in np.array(Initial_gradients_100000[axis_ind]):
                                                import pdb
                                                pdb.set_trace()
                                            Initial_gradients_100000[axis_ind].append([mean_grad])
                                        else:
                                            print('Not enough points to suggest system stays bound for first  10000yr')
                                    except:
                                        print('system has not mean times < 100000yr')
                                    if len(grad) > 0:
                                        if grad[0] in np.array(Initial_gradients[axis_ind]):
                                            import pdb
                                            pdb.set_trace()
                                        Initial_gradients[axis_ind].append(grad[0])
                                    #lets trying smoothing over a 1000year window, but not centred
                                    
                                    Used_sub_sys.append(str(sub_sys))
                                    
                                    '''
                                if plot_gradient == True:
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
                                if str(sub_sys) not in Sub_sys_indices:
                                    Sub_sys_indices.append(str(sub_sys))
                                    replace_string = str(new_sys_id)
                                    new_sys_id = new_sys_id + 1
                                else:
                                    replace_string = str(Sub_sys_indices.index(str(sub_sys)) + first_new_sys_id)
                                sys_comps = sys_comps[:open_ind] + replace_string + sys_comps[char_it+1:]
                                if '[' not in sys_comps:
                                    reduced = True
                                break
                    
                if plot_gradient == True:
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

            #plt.savefig(savename+'.pdf', format='pdf', bbox_inches='tight')
            #print('Created '+ savename+'.png')
        if size > 1:
            grad_pickle = 'grad_pickle_'+str(rank)+'.pkl'
        else:
            grad_pickle = 'grad_pickle.pkl'
        file = open(grad_pickle, 'wb')
        #pickle.dump((Initial_gradients, Initial_gradients_10000, Grad_1e4, Grad_1e3), file)
        pickle.dump((Initial_gradients, Initial_gradients_1000, Initial_gradients_10000, Initial_gradients_100000, Grad_1e4, Grad_1e3), file)
        file.close()
        print('saved_gradients')
        
        #Compile together pickles
        if size > 1:
            import pdb
            pdb.set_trace()

grad_pickle = 'grad_pickle.pkl'
file = open(grad_pickle, 'rb')
Initial_gradients, Initial_gradients_1000, Initial_gradients_10000, Initial_gradients_100000, Grad_1e4, Grad_1e3 = pickle.load(file)
file.close()

#Defining gradient bins and getting tick labels
grad_bins = np.concatenate((-1*np.logspace(5,-6,12), np.array([0, 1.e10]))) #np.concatenate((-1*np.logspace(5,-3,9), np.array([0, 1.e10])))
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

#Plotting initial gradients
#grad_bins = np.concatenate((-1*np.logspace(4,-3,15), np.array([0, 1.e10])))
grad_hist_core, grad_bins = np.histogram(Initial_gradients[0], bins=grad_bins)
grad_hist_core_delayed, grad_bins = np.histogram(Initial_gradients[1], bins=grad_bins)
grad_hist_capt, grad_bins = np.histogram(Initial_gradients[2], bins=grad_bins)
grad_hist_misc, grad_bins = np.histogram(Initial_gradients[3], bins=grad_bins)

x_range = np.arange(len(grad_hist_core)+1)

grad_hist_core_norm = grad_hist_core/np.sum(grad_hist_core)
grad_hist_core_delayed_norm = grad_hist_core_delayed/np.sum(grad_hist_core_delayed)
grad_hist_capt_norm = grad_hist_capt/np.sum(grad_hist_capt)
grad_hist_misc_norm = grad_hist_misc/np.sum(grad_hist_misc)

grad_hist_core_norm = np.concatenate((grad_hist_core_norm, np.array([grad_hist_core_norm[-1]])))
grad_hist_core_delayed_norm = np.concatenate((grad_hist_core_delayed_norm, np.array([grad_hist_core_delayed_norm[-1]])))
grad_hist_capt_norm = np.concatenate((grad_hist_capt_norm, np.array([grad_hist_capt_norm[-1]])))
grad_hist_misc_norm = np.concatenate((grad_hist_misc_norm, np.array([grad_hist_misc_norm[-1]])))

plt.clf()
fig, ax = plt.subplots()
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
ax.set_xticklabels(ticklabels[::2])
ax.set_xlabel('Inspiral rate (au/yr)')
ax.set_ylabel('#')
ax.set_ylim([0,0.4])
ax.legend(loc='best')
plt.savefig('Initial_grad_hist.png')

mean_grads = [Initial_gradients_1000, Initial_gradients_10000, Initial_gradients_100000]

#Plotting mean gradients
#calculate means
time_means = [1000, 10000, 100000]
time_means_counter = 0
for Initial_mean_grad in mean_grads:
    core_mean = []
    for grit in range(len(Initial_mean_grad[0])):
        core_mean.append(np.mean(Initial_mean_grad[0][grit]))
    grad_hist_core_mean, grad_bins = np.histogram(core_mean, bins=grad_bins)
    grad_hist_core_mean_std = np.sqrt(grad_hist_core_mean)
    grad_hist_core_mean_rel_err = grad_hist_core_mean_std/grad_hist_core_mean
    grad_hist_core_mean_rel_err = np.nan_to_num(grad_hist_core_mean_rel_err)
    #grad_hist_core_mean_err = 1/np.sqrt(len(core_mean))
    #grad_hist_core_mean_rel_err = grad_hist_core_mean_err/grad_hist_core_mean
    #grad_hist_core_mean_rel_err[np.isinf(grad_hist_core_mean_rel_err)] = 0

    core_delayed_mean = []
    for grit in range(len(Initial_mean_grad[1])):
        core_delayed_mean.append(np.mean(Initial_mean_grad[1][grit]))
    grad_hist_core_delayed_mean, grad_bins = np.histogram(core_delayed_mean, bins=grad_bins)
    grad_hist_core_delayed_mean_std = np.sqrt(grad_hist_core_delayed_mean)
    grad_hist_core_delayed_mean_rel_err = grad_hist_core_delayed_mean_std/grad_hist_core_delayed_mean
    grad_hist_core_delayed_mean_rel_err = np.nan_to_num(grad_hist_core_delayed_mean_rel_err)
    #grad_hist_core_delayed_mean_err = 1/np.sqrt(len(core_delayed_mean))
    #grad_hist_core_delayed_mean_rel_err = grad_hist_core_delayed_mean_err/grad_hist_core_delayed_mean
    #grad_hist_core_delayed_mean_rel_err[np.isinf(grad_hist_core_delayed_mean_rel_err)] = 0

    capt_mean = []
    for grit in range(len(Initial_mean_grad[2])):
        capt_mean.append(np.mean(Initial_mean_grad[2][grit]))
    grad_hist_capt_mean, grad_bins = np.histogram(capt_mean, bins=grad_bins)
    grad_hist_capt_mean_std = np.sqrt(grad_hist_capt_mean)
    grad_hist_capt_mean_rel_err = grad_hist_capt_mean_std/grad_hist_capt_mean
    grad_hist_capt_mean_rel_err = np.nan_to_num(grad_hist_capt_mean_rel_err)

    misc_mean = []
    for grit in range(len(Initial_mean_grad[3])):
        misc_mean.append(np.mean(Initial_mean_grad[3][grit]))
    grad_hist_misc_mean, grad_bins = np.histogram(misc_mean, bins=grad_bins)
    grad_hist_misc_mean_std = np.sqrt(grad_hist_misc_mean)
    grad_hist_misc_mean_rel_err = grad_hist_misc_mean_std/grad_hist_misc_mean
    grad_hist_misc_mean_rel_err = np.nan_to_num(grad_hist_misc_mean_rel_err)

    grad_hist_core_mean_norm = grad_hist_core_mean/np.sum(grad_hist_core_mean)
    grad_hist_core_delayed_mean_norm = grad_hist_core_delayed_mean/np.sum(grad_hist_core_delayed_mean)
    grad_hist_capt_mean_norm = grad_hist_capt_mean/np.sum(grad_hist_capt_mean)
    grad_hist_misc_mean_norm = grad_hist_misc_mean/np.sum(grad_hist_misc_mean)

    grad_hist_core_mean_norm = np.concatenate((grad_hist_core_mean_norm, np.array([grad_hist_core_mean_norm[-1]])))
    grad_hist_core_delayed_mean_norm = np.concatenate((grad_hist_core_delayed_mean_norm, np.array([grad_hist_core_delayed_mean_norm[-1]])))
    grad_hist_capt_mean_norm = np.concatenate((grad_hist_capt_mean_norm, np.array([grad_hist_capt_mean_norm[-1]])))
    grad_hist_misc_mean_norm = np.concatenate((grad_hist_misc_mean_norm, np.array([grad_hist_misc_mean_norm[-1]])))

    grad_hist_core_mean_rel_err = np.concatenate((grad_hist_core_mean_rel_err, np.array([grad_hist_core_mean_rel_err[-1]])))
    grad_hist_core_delayed_mean_rel_err = np.concatenate((grad_hist_core_delayed_mean_rel_err, np.array([grad_hist_core_delayed_mean_rel_err[-1]])))
    grad_hist_capt_mean_rel_err = np.concatenate((grad_hist_capt_mean_rel_err, np.array([grad_hist_capt_mean_rel_err[-1]])))
    grad_hist_misc_mean_rel_err = np.concatenate((grad_hist_misc_mean_rel_err, np.array([grad_hist_misc_mean_rel_err[-1]])))

    plt.clf()
    fig, ax = plt.subplots()
    ax.step(x_range, grad_hist_core_mean_norm, where='post', label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
    ax.errorbar(x_range+0.5, grad_hist_core_mean_norm, yerr=(grad_hist_core_mean_rel_err*grad_hist_core_mean_norm), fmt='none', linewidth=2, color='b', alpha=0.5)
    ax.step(x_range, grad_hist_core_delayed_mean_norm, where='post', label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
    ax.errorbar(x_range+0.5, grad_hist_core_delayed_mean_norm, yerr=(grad_hist_core_delayed_mean_rel_err*grad_hist_core_delayed_mean_norm), fmt='none', linewidth=2, color='purple', alpha=0.5)
    ax.step(x_range, grad_hist_capt_mean_norm, where='post', label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
    ax.errorbar(x_range+0.5, grad_hist_capt_mean_norm, yerr=(grad_hist_capt_mean_rel_err*grad_hist_capt_mean_norm), fmt='none', linewidth=2, color='red', alpha=0.5)
    ax.step(x_range, grad_hist_misc_mean_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
    ax.errorbar(x_range+0.5, grad_hist_misc_mean_norm, yerr=(grad_hist_misc_mean_rel_err*grad_hist_misc_mean_norm), fmt='none', linewidth=2, color='orange', alpha=0.5)
    '''
    ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_mean/np.sum(grad_hist_core_mean), label="Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='b')
    ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_core_delayed_mean/np.sum(grad_hist_core_delayed_mean), label="Delayed Core Fragmentation", width=1, color='None', linewidth=2, edgecolor='purple')
    ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_capt_mean/np.sum(grad_hist_capt_mean), label="Dynamical Capture", width=1, color='None', linewidth=2, edgecolor='red')
    ax.bar(np.arange(len(grad_hist_core))+0.5, grad_hist_misc_mean/np.sum(grad_hist_misc_mean), label="Other", width=1, color='None', linewidth=2, edgecolor='orange')
    '''
    ax.set_xlim([x_range[0], x_range[-1]])
    ax.set_xticklabels(ticklabels[::2])
    ax.set_xlabel('Inspiral rate (au/yr)')
    ax.set_ylabel('#')
    ax.set_ylim(bottom=0)
    ax.legend(loc='best')
    plt.savefig('Initial_mean_grad_'+str(time_means[time_means_counter])+'.png')
    time_means_counter = time_means_counter + 1

"""

def Gaussian(x,scale,mean,sigma):
    return scale*stats.norm.pdf(x, mean, sigma)
    
def Skewed_Gaussian(x, scale, mean, sigma, skew, kurt):
    return scale*stats.skewnorm.pdf(x,mean, sigma, skew, kurt)
    
#Fitting
import scipy.stats as stats
from scipy.optimize import curve_fit
x_log = np.log10(-1*grad_bins[:-2])-0.5
x_log_fit = np.linspace(x_log[0], x_log[-1], 100)

y_log_core = grad_hist_core_mean_norm[:-2]#[::-1]
y_log_core_err = (grad_hist_core_mean_rel_err*grad_hist_core_mean_norm)[:-2]#[::-1]
y_log_core_delayed = grad_hist_core_delayed_mean_norm[:-2]#[::-1]
y_log_core_delayed_err = (grad_hist_core_delayed_mean_rel_err*grad_hist_core_delayed_mean_norm)[:-2]#[::-1]
y_log_capt = grad_hist_capt_mean_norm[:-2]#[::-1]
y_log_capt_err = (grad_hist_capt_mean_rel_err*grad_hist_capt_mean_norm)[:-2]#[::-1]
y_log_misc = grad_hist_misc_mean_norm[:-2]#[::-1]
y_log_misc_err = (grad_hist_misc_mean_rel_err*grad_hist_misc_mean_norm)[:-2]#[::-1]

mean_guess = x_log[np.argmax(y_log_core)]#-0.5
mode_guess = x_log[np.argmax(y_log_core)]

sig_guess = np.std(np.log10(-1*np.array(core_mean)[np.array(core_mean)<0]))
scale_guess = np.sum(grad_hist_core_mean[:-2])/np.sum(grad_hist_core_mean)
popt_core, pcov_core = curve_fit(Gaussian, x_log, y_log_core, [scale_guess, mean_guess, sig_guess], sigma=y_log_core_err)#, bounds=((0, , 0, 0), (1, np.inf, 1, 1))))
fit_core = Gaussian(x_log_fit, *popt_core)
print('Core frag mean inspiral:', popt_core[1], '\pm', popt_core[2])

usable_inds = np.where((np.isnan(core_delayed_mean)==False)&(np.array(core_delayed_mean)<0))[0]
core_delayed_power = np.log10(-1*np.array(core_delayed_mean)[usable_inds])
mean_guess = np.mean(core_delayed_power)#-0.5
median_guess = np.median(core_delayed_power)
mode_guess = x_log[np.argmax(y_log_core_delayed)]
sig_guess = np.std(core_delayed_power)
skew_guess = (3*(mean_guess - median_guess))/sig_guess
scale_guess = np.sum(grad_hist_core_delayed_mean[:-2])/np.sum(grad_hist_core_delayed_mean)
popt_core_delayed, pcov_core_delayed = curve_fit(Gaussian, x_log, y_log_core_delayed, [scale_guess, mean_guess, sig_guess], sigma=y_log_core_delayed_err)
fit_core_delayed = Gaussian(x_log_fit, *popt_core_delayed)
print('Delayed core frag mean inspiral:', popt_core_delayed[1], '\pm', popt_core_delayed[2])

plt.clf()
plt.errorbar(x_log, y_log_core_delayed, yerr=y_log_core_delayed_err, marker='o', color='purple', linestyle="")
#y_skew = scale_guess*skewnorm.pdf(x_log_fit, skew_guess, loc=mode, scale=sig_guess)
y_skew = skewnorm.pdf(x_log_fit, -4, loc=2, scale=2)
plt.plot(x_log_fit, y_skew)
plt.xlim([x_log[0], x_log[-1]])
plt.ylim(bottom=0)
plt.savefig('skew_test.png')

mean_guess = x_log[np.argmax(y_log_capt)]#-0.5
sig_guess = np.std(np.log10(-1*np.array(capt_mean)[np.array(capt_mean)<0]))
scale_guess = np.sum(grad_hist_capt_mean[:-2])/np.sum(grad_hist_capt_mean)
popt_capt, pcov_core = curve_fit(Gaussian, x_log, y_log_capt, [scale_guess, mean_guess, sig_guess], sigma=y_log_capt_err)
fit_capt = Gaussian(x_log_fit, *popt_capt)
print('Core frag mean inspiral:', popt_capt[1], '\pm', popt_capt[2])

mean_guess = x_log[np.argmax(y_log_misc)]#-0.5
sig_guess = np.std(np.log10(-1*np.array(misc_mean)[np.array(misc_mean)<0]))
scale_guess = np.sum(grad_hist_misc_mean[:-2])/np.sum(grad_hist_misc_mean)
popt_misc, pcov_misc = curve_fit(Gaussian, x_log, y_log_misc, [scale_guess, mean_guess, sig_guess], sigma=y_log_misc_err)
fit_misc = Gaussian(x_log_fit, *popt_misc)
print('Core frag mean inspiral:', popt_misc[1], '\pm', popt_misc[2])


plt.clf()
plt.errorbar(x_log, y_log_core, yerr=y_log_core_err, marker='o', color='b', linestyle="")
plt.plot(x_log_fit, fit_core, color='b')

plt.errorbar(x_log, y_log_core_delayed, yerr=y_log_core_delayed_err, marker='o', color='purple', linestyle="")
plt.plot(x_log_fit, fit_core_delayed, color='purple')

plt.errorbar(x_log, y_log_capt, yerr=y_log_capt_err, marker='o', color='red', linestyle="")
plt.plot(x_log_fit, fit_capt, color='red')

plt.errorbar(x_log, y_log_misc, yerr=y_log_misc_err, marker='o', color='orange', linestyle="")
plt.plot(x_log_fit, fit_misc, color='orange')
#plt.xscale("log", nonpositive='clip')
#plt.semilogx(x_log, y_log)
plt.xlim([x_log[0], x_log[-1]])
plt.ylim(bottom=0)
plt.savefig('log_normal_fit.png')
"""
