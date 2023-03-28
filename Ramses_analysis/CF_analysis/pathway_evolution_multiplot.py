import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
from scipy import stats
from scipy.optimize import curve_fit

def Gaussian(x,scale,mean,sigma):
    return scale*stats.norm.pdf(x, mean, sigma)
    
def Gaussian_cdf(x,scale,mean,sigma):
    return scale*stats.norm.cdf(x, mean, sigma)

def Skewed_Gaussian(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.pdf(x, skew, loc=mean, scale=sigma)
    
def Skewed_Gaussian_cdf(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.cdf(x, skew, loc=mean, scale=sigma)

grad_pickles = ['/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G50/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G100/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G125/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G150/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G200/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G400/1000_yr/grad_pickle.pkl']

Sim_ids = ['G50', 'G100', 'G125', 'G150', 'G200', 'G400']
#Defining gradient bins and getting tick labels
grad_bins = np.concatenate((-1*np.logspace(1.5,-4,12), np.array([0, 1.e10])))#np.concatenate((-1*np.logspace(0.5,-5,23), np.array([0, 1.e10])))
#inpiral only:
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

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(grad_pickles), figsize=(single_col_width*1.5, single_col_width*3), sharex=True, sharey=True)#, hspace=0.0)
iter_range = range(0, len(grad_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

Mean_grads = [[], [], [], []]
Median_grads = [[], [], [], []]
Std_grads = [[], [], [], []]
All_grads = [[], [], []]
Grad_over_sep_mean = [[], [], [], []]
Grad_over_sep_median = [[], [], [], []]
Grad_over_sep_std = [[], [], [], []]
Grad_over_sep_all = [[], [], [], []]
for grad_it in range(len(grad_pickles)):
    file = open(grad_pickles[grad_it], 'rb')
    Initial_rel_vel, Initial_eccentricity, Initial_gradients, Grad_over_init_sep, Grad_1e4, Grad_1e3 = pickle.load(file)
    file.close()

    core_inspiral_inds = np.where(np.array(Initial_gradients[0]) < 0)[0]
    delayed_core_inspiral_inds = np.where(np.array(Initial_gradients[1]) < 0)[0]
    capt_inspiral_inds = np.where(np.array(Initial_gradients[2]) < 0)[0]
    other_inspiral_inds = np.where(np.array(Initial_gradients[3]) < 0)[0]
    
    All_grads[0] = All_grads[0] + np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds]).tolist()
    All_grads[1] = All_grads[1] + np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds]).tolist()
    All_grads[2] = All_grads[2] + np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds]).tolist()
    
    Mean_grads[0].append(np.mean(np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds])))
    Median_grads[0].append(np.median(np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds])))
    Std_grads[0].append(np.std(np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds])))
    
    Mean_grads[1].append(np.mean(np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds])))
    Median_grads[1].append(np.median(np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds])))
    Std_grads[1].append(np.std(np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds])))
    
    Mean_grads[2].append(np.mean(np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds])))
    Median_grads[2].append(np.median(np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds])))
    Std_grads[2].append(np.std(np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds])))
    
    Mean_grads[3].append(np.mean(np.log10(-1*np.array(Initial_gradients[3]).T[0][other_inspiral_inds])))
    Median_grads[3].append(np.median(np.log10(-1*np.array(Initial_gradients[3]).T[0][other_inspiral_inds])))
    Std_grads[3].append(np.std(np.log10(-1*np.array(Initial_gradients[3]).T[0][other_inspiral_inds])))
    
    #======================================================================================================
    Grad_over_sep_all[0] = Grad_over_sep_all[0] + np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds]).tolist()
    Grad_over_sep_all[1] = Grad_over_sep_all[1] + np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds]).tolist()
    Grad_over_sep_all[2] = Grad_over_sep_all[2] + np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds]).tolist()
    
    Grad_over_sep_mean[0].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds])))
    Grad_over_sep_median[0].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds])))
    Grad_over_sep_std[0].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds])))
    
    Grad_over_sep_mean[1].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds])))
    Grad_over_sep_median[1].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds])))
    Grad_over_sep_std[1].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds])))
    
    Grad_over_sep_mean[2].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds])))
    Grad_over_sep_median[2].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds])))
    Grad_over_sep_std[2].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds])))
    
    Grad_over_sep_mean[3].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[3]).T[0][other_inspiral_inds])))
    Grad_over_sep_median[3].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[3]).T[0][other_inspiral_inds])))
    Grad_over_sep_std[3].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[3]).T[0][other_inspiral_inds])))
    #======================================================================================================
    
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

    axs[grad_it].step(x_range, grad_hist_core_norm, where='post', label="Bound core frag.", linewidth=2, color='b', alpha=0.5, ls='-')
    axs[grad_it].step(x_range, grad_hist_core_delayed_norm, where='post', label="Unbound core frag.", linewidth=2, color='purple', alpha=0.5, ls='--')
    axs[grad_it].step(x_range, grad_hist_capt_norm, where='post', label="Dynamical capture", linewidth=2, color='red', alpha=0.5, ls='-.')
    #axs[grad_it].step(x_range, grad_hist_misc_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
    
    axs[grad_it].set_ylabel('')
    
    axs[grad_it].tick_params(axis='both', which='major', labelsize=font_size)
    axs[grad_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[grad_it].tick_params(axis='x', direction='in')
    axs[grad_it].tick_params(axis='y', direction='in')
    
    axs[grad_it].set_ylim([0, 0.8])
    axs[grad_it].set_xlim([x_range[0], x_range[-1]])

fig.savefig('Initial_grad_hist.png', bbox_inches='tight', pad_inches=0.02)

grad_pickles = ['/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G50/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G100/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G125/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G150/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G200/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Inspiral_rates/G400/grad_pickle.pkl']

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(grad_pickles), figsize=(single_col_width*1.5, single_col_width*3), sharex=True, sharey=True)#, hspace=0.0)
iter_range = range(0, len(grad_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

Mean_grads_10000 = [[], [], [], []]
Median_grads_10000 = [[], [], [], []]
Std_grads_10000 = [[], [], [], []]
All_grads_10000 = [[], [], []]
Grad_over_sep_mean_10000 = [[], [], [], []]
Grad_over_sep_median_10000 = [[], [], [], []]
Grad_over_sep_std_10000 = [[], [], [], []]
Grad_over_sep_all_10000 = [[], [], [], []]
for grad_it in range(len(grad_pickles)):
    file = open(grad_pickles[grad_it], 'rb')
    Initial_rel_vel, Initial_eccentricity, Initial_gradients, Grad_over_init_sep, Grad_1e4, Grad_1e3 = pickle.load(file)
    file.close()

    core_inspiral_inds = np.where(np.array(Initial_gradients[0]) < 0)[0]
    delayed_core_inspiral_inds = np.where(np.array(Initial_gradients[1]) < 0)[0]
    capt_inspiral_inds = np.where(np.array(Initial_gradients[2]) < 0)[0]
    other_inspiral_inds = np.where(np.array(Initial_gradients[3]) < 0)[0]
    
    All_grads_10000[0] = All_grads_10000[0] + np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds]).tolist()
    All_grads_10000[1] = All_grads_10000[1] + np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds]).tolist()
    All_grads_10000[2] = All_grads_10000[2] + np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds]).tolist()
    
    Mean_grads_10000[0].append(np.mean(np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds])))
    Median_grads_10000[0].append(np.median(np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds])))
    Std_grads_10000[0].append(np.std(np.log10(-1*np.array(Initial_gradients[0]).T[0][core_inspiral_inds])))
    
    Mean_grads_10000[1].append(np.mean(np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds])))
    Median_grads_10000[1].append(np.median(np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds])))
    Std_grads_10000[1].append(np.std(np.log10(-1*np.array(Initial_gradients[1]).T[0][delayed_core_inspiral_inds])))
    
    Mean_grads_10000[2].append(np.mean(np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds])))
    Median_grads_10000[2].append(np.median(np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds])))
    Std_grads_10000[2].append(np.std(np.log10(-1*np.array(Initial_gradients[2]).T[0][capt_inspiral_inds])))
    
    Mean_grads_10000[3].append(np.mean(np.log10(-1*np.array(Initial_gradients[3]).T[0][other_inspiral_inds])))
    Median_grads_10000[3].append(np.median(np.log10(-1*np.array(Initial_gradients[3]).T[0][other_inspiral_inds])))
    Std_grads_10000[3].append(np.std(np.log10(-1*np.array(Initial_gradients[3]).T[0][other_inspiral_inds])))
    
    #======================================================================================================
    Grad_over_sep_all_10000[0] = Grad_over_sep_all_10000[0] + np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds]).tolist()
    Grad_over_sep_all_10000[1] = Grad_over_sep_all_10000[1] + np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds]).tolist()
    Grad_over_sep_all_10000[2] = Grad_over_sep_all_10000[2] + np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds]).tolist()
    
    Grad_over_sep_mean_10000[0].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds])))
    Grad_over_sep_median_10000[0].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds])))
    Grad_over_sep_std_10000[0].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[0]).T[0][core_inspiral_inds])))
    
    Grad_over_sep_mean_10000[1].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds])))
    Grad_over_sep_median_10000[1].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds])))
    Grad_over_sep_std_10000[1].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[1]).T[0][delayed_core_inspiral_inds])))
    
    Grad_over_sep_mean_10000[2].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds])))
    Grad_over_sep_median_10000[2].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds])))
    Grad_over_sep_std_10000[2].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[2]).T[0][capt_inspiral_inds])))
    
    Grad_over_sep_mean_10000[3].append(np.mean(np.log10(-1*np.array(Grad_over_init_sep[3]).T[0][other_inspiral_inds])))
    Grad_over_sep_median_10000[3].append(np.median(np.log10(-1*np.array(Grad_over_init_sep[3]).T[0][other_inspiral_inds])))
    Grad_over_sep_std_10000[3].append(np.std(np.log10(-1*np.array(Grad_over_init_sep[3]).T[0][other_inspiral_inds])))
    #======================================================================================================
    
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

    axs[grad_it].step(x_range, grad_hist_core_norm, where='post', label="Bound core frag.", linewidth=2, color='b', alpha=0.5, ls='-')
    axs[grad_it].step(x_range, grad_hist_core_delayed_norm, where='post', label="Unbound core frag.", linewidth=2, color='purple', alpha=0.5, ls='--')
    axs[grad_it].step(x_range, grad_hist_capt_norm, where='post', label="Dynamical capture", linewidth=2, color='red', alpha=0.5, ls='-.')
    axs[grad_it].step(x_range, grad_hist_misc_norm, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
    
    axs[grad_it].set_ylabel('')
    
    axs[grad_it].tick_params(axis='both', which='major', labelsize=font_size)
    axs[grad_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[grad_it].tick_params(axis='x', direction='in')
    axs[grad_it].tick_params(axis='y', direction='in')
    
    axs[grad_it].set_ylim([0, 0.8])
    axs[grad_it].set_xlim([x_range[0], x_range[-1]])

fig.savefig('Initial_grad_hist.png', bbox_inches='tight', pad_inches=0.02)

masses = [1500, 3000, 3750, 4500, 6000, 12000]
plt.clf()
plt.errorbar(masses, Mean_grads[0], yerr=Std_grads[0], label='Core Fragmentation', color='b')
plt.errorbar(masses, Mean_grads[1], yerr=Std_grads[1], label='Delayed Core Fragmentation', color='purple')
plt.errorbar(masses, Mean_grads[2], yerr=Std_grads[2], label='Dynamical Capture', color='r')
plt.errorbar(masses, Mean_grads[3], yerr=Std_grads[3], label='Other', color='orange')
plt.legend()
plt.xlabel('Gas Mass')
plt.ylabel('Inspiral rate (Log$_{10}$(AU/yr))')
plt.savefig('inspiral_rate_comparison_means')

#==============================================================================

plt.clf()
plt.errorbar(masses, Grad_over_sep_mean[0], yerr=Grad_over_sep_std[0], label='Core Fragmentation', color='b')
plt.errorbar(masses, Grad_over_sep_mean[1], yerr=Grad_over_sep_std[1], label='Delayed Core Fragmentation', color='purple')
plt.errorbar(masses, Grad_over_sep_mean[2], yerr=Grad_over_sep_std[2], label='Dynamical Capture', color='r')
plt.errorbar(masses, Grad_over_sep_mean[3], yerr=Grad_over_sep_std[3], label='Other', color='orange')
plt.legend()
plt.xlabel('Gas Mass')
plt.ylabel('Inspiral rate (Log$_{10}$(AU/yr))')
plt.savefig('inspiral_rate_comparison_means')

#==============================================================================

Core_bounds = [np.array(Mean_grads[0])-np.array(Std_grads[0]), np.array(Mean_grads[0])+np.array(Std_grads[0])]
Delayed_core_bounds = [np.array(Mean_grads[1])-np.array(Std_grads[1]), np.array(Mean_grads[1])+np.array(Std_grads[1])]
Capt_bounds = [np.array(Mean_grads[2])-np.array(Std_grads[2]), np.array(Mean_grads[2])+np.array(Std_grads[2])]
Other_bounds = [np.array(Mean_grads[3])-np.array(Std_grads[3]), np.array(Mean_grads[3])+np.array(Std_grads[3])]

Core_err = [np.array(Median_grads[0]) - Core_bounds[0], Core_bounds[1] - np.array(Median_grads[0])]
Delayed_core_err = [np.array(Median_grads[1]) - Delayed_core_bounds[0], Delayed_core_bounds[1] - np.array(Median_grads[1])]
Capt_err = [np.array(Median_grads[2]) - Capt_bounds[0], Capt_bounds[1] - np.array(Median_grads[2])]
Other_err = [np.array(Median_grads[3]) - Other_bounds[0], Other_bounds[1] - np.array(Median_grads[3])]

Core_bounds_10000 = [np.array(Mean_grads_10000[0])-np.array(Std_grads_10000[0]), np.array(Mean_grads_10000[0])+np.array(Std_grads_10000[0])]
Delayed_core_bounds_10000 = [np.array(Mean_grads_10000[1])-np.array(Std_grads_10000[1]), np.array(Mean_grads_10000[1])+np.array(Std_grads_10000[1])]
Capt_bounds_10000 = [np.array(Mean_grads_10000[2])-np.array(Std_grads_10000[2]), np.array(Mean_grads_10000[2])+np.array(Std_grads_10000[2])]
Other_bounds_10000 = [np.array(Mean_grads_10000[3])-np.array(Std_grads_10000[3]), np.array(Mean_grads_10000[3])+np.array(Std_grads_10000[3])]

Core_err_10000 = [np.array(Median_grads_10000[0]) - Core_bounds_10000[0], Core_bounds_10000[1] - np.array(Median_grads_10000[0])]
Delayed_core_err_10000 = [np.array(Median_grads_10000[1]) - Delayed_core_bounds_10000[0], Delayed_core_bounds_10000[1] - np.array(Median_grads_10000[1])]
Capt_err_10000 = [np.array(Median_grads_10000[2]) - Capt_bounds_10000[0], Capt_bounds_10000[1] - np.array(Median_grads_10000[2])]
Other_err_10000 = [np.array(Median_grads_10000[3]) - Other_bounds_10000[0], Other_bounds_10000[1] - np.array(Median_grads_10000[3])]

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(single_col_width, single_col_width*1.5), sharex=True, sharey=True)#, hspace=0.0)
iter_range = range(0, len(grad_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

axs[0].errorbar(np.array(masses)-130, Median_grads[0], yerr=Core_err, label='Bound core frag.', color='b')
axs[0].errorbar(np.array(masses), Median_grads[1], yerr=Delayed_core_err, label='Unbound core frag.', color='purple')
axs[0].errorbar(np.array(masses)+130, Median_grads[2], yerr=Capt_err, label='Dynamical capture', color='r')

axs[1].errorbar(np.array(masses)-130, Median_grads_10000[0], yerr=Core_err_10000, label='Bound core frag.', color='b')
axs[1].errorbar(np.array(masses), Median_grads_10000[1], yerr=Delayed_core_err_10000, label='Unbound core frag.', color='purple')
axs[1].errorbar(np.array(masses)+130, Median_grads_10000[2], yerr=Capt_err_10000, label='Dynamical capture', color='r')

axs[0].tick_params(which='both', direction='in')
axs[0].tick_params(axis='both', which='major', labelsize=font_size, right=True, top=True)
axs[0].tick_params(axis='both', which='minor', labelsize=font_size, right=True, top=True)
axs[0].text(6500, -2.1, "Baseline=$1\,000\,\mathrm{yr}$", zorder=11, size=font_size)

axs[1].tick_params(which='both', direction='in')
axs[1].tick_params(axis='both', which='major', labelsize=font_size, right=True, top=True)
axs[1].tick_params(axis='both', which='minor', labelsize=font_size, right=True, top=True)
axs[1].text(6500, -0.5, "Baseline=$10\,000\,\mathrm{yr}$", zorder=11, size=font_size)

axs[0].legend(loc='upper right', fontsize=font_size)
axs[1].set_xlabel('Molecular cloud mass (M$_\odot$)', size=font_size)
axs[0].set_ylabel('Inspiral rate (Log$_{10}$(AU/yr))', size=font_size)
axs[1].set_ylabel('Inspiral rate (Log$_{10}$(AU/yr))', size=font_size)
#plt.ylim(top=1.5)
plt.savefig('inspiral_rate_comparison_medians.pdf', bbox_inches='tight', pad_inches=0.02)

#=================================================================================================
Core_bounds = [np.array(Grad_over_sep_mean[0])-np.array(Grad_over_sep_std[0]), np.array(Grad_over_sep_mean[0])+np.array(Grad_over_sep_std[0])]
Delayed_core_bounds = [np.array(Grad_over_sep_mean[1])-np.array(Grad_over_sep_std[1]), np.array(Grad_over_sep_mean[1])+np.array(Grad_over_sep_std[1])]
Capt_bounds = [np.array(Grad_over_sep_mean[2])-np.array(Grad_over_sep_std[2]), np.array(Grad_over_sep_mean[2])+np.array(Grad_over_sep_std[2])]
Other_bounds = [np.array(Grad_over_sep_mean[3])-np.array(Grad_over_sep_std[3]), np.array(Grad_over_sep_mean[3])+np.array(Grad_over_sep_std[3])]

Core_err = [np.array(Grad_over_sep_median[0]) - Core_bounds[0], Core_bounds[1] - np.array(Grad_over_sep_median[0])]
Delayed_core_err = [np.array(Grad_over_sep_median[1]) - Delayed_core_bounds[0], Delayed_core_bounds[1] - np.array(Grad_over_sep_median[1])]
Capt_err = [np.array(Grad_over_sep_median[2]) - Capt_bounds[0], Capt_bounds[1] - np.array(Grad_over_sep_median[2])]
Other_err = [np.array(Grad_over_sep_median[3]) - Other_bounds[0], Other_bounds[1] - np.array(Grad_over_sep_median[3])]

Core_bounds_10000 = [np.array(Grad_over_sep_mean_10000[0])-np.array(Grad_over_sep_std_10000[0]), np.array(Grad_over_sep_mean_10000[0])+np.array(Grad_over_sep_std_10000[0])]
Delayed_core_bounds_10000 = [np.array(Grad_over_sep_mean_10000[1])-np.array(Grad_over_sep_std_10000[1]), np.array(Grad_over_sep_mean_10000[1])+np.array(Grad_over_sep_std_10000[1])]
Capt_bounds_10000 = [np.array(Grad_over_sep_mean_10000[2])-np.array(Grad_over_sep_std_10000[2]), np.array(Grad_over_sep_mean_10000[2])+np.array(Grad_over_sep_std_10000[2])]
Other_bounds_10000 = [np.array(Grad_over_sep_mean_10000[3])-np.array(Grad_over_sep_std_10000[3]), np.array(Grad_over_sep_mean_10000[3])+np.array(Grad_over_sep_std_10000[3])]

Core_err_10000 = [np.array(Grad_over_sep_median_10000[0]) - Core_bounds_10000[0], Core_bounds_10000[1] - np.array(Grad_over_sep_median_10000[0])]
Delayed_core_err_10000 = [np.array(Grad_over_sep_median_10000[1]) - Delayed_core_bounds_10000[0], Delayed_core_bounds_10000[1] - np.array(Grad_over_sep_median_10000[1])]
Capt_err_10000 = [np.array(Grad_over_sep_median_10000[2]) - Capt_bounds_10000[0], Capt_bounds_10000[1] - np.array(Grad_over_sep_median_10000[2])]
Other_err_10000 = [np.array(Grad_over_sep_median_10000[3]) - Other_bounds_10000[0], Other_bounds_10000[1] - np.array(Grad_over_sep_median_10000[3])]

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(single_col_width, single_col_width*1.5), sharex=True, sharey=True)#, hspace=0.0)
iter_range = range(0, len(grad_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

axs[0].errorbar(np.array(masses)-130, Grad_over_sep_median[0], yerr=Core_err, label='Bound core frag.', color='b')
axs[0].errorbar(np.array(masses), Grad_over_sep_median[1], yerr=Delayed_core_err, label='Unbound core frag.', color='purple')
axs[0].errorbar(np.array(masses)+130, Grad_over_sep_median[2], yerr=Capt_err, label='Dynamical capture', color='r')

axs[1].errorbar(np.array(masses)-130, Grad_over_sep_median_10000[0], yerr=Core_err_10000, label='Bound core frag.', color='b')
axs[1].errorbar(np.array(masses), Grad_over_sep_median_10000[1], yerr=Delayed_core_err_10000, label='Unbound core frag.', color='purple')
axs[1].errorbar(np.array(masses)+130, Grad_over_sep_median_10000[2], yerr=Capt_err_10000, label='Dynamical capture', color='r')

axs[0].tick_params(which='both', direction='in')
axs[0].tick_params(axis='both', which='major', labelsize=font_size, right=True, top=True)
axs[0].tick_params(axis='both', which='minor', labelsize=font_size, right=True, top=True)
axs[0].text(6500, -4.25, "Baseline=$1\,000\,\mathrm{yr}$", zorder=11, size=font_size)

axs[1].tick_params(which='both', direction='in')
axs[1].tick_params(axis='both', which='major', labelsize=font_size, right=True, top=True)
axs[1].tick_params(axis='both', which='minor', labelsize=font_size, right=True, top=True)
axs[1].text(6500, -3.5, "Baseline=$10\,000\,\mathrm{yr}$", zorder=11, size=font_size)

axs[0].legend(loc='lower right', fontsize=font_size)
axs[1].set_xlabel('Molecular cloud mass (M$_\odot$)', size=font_size)
axs[0].set_ylabel('Inspiral rate (Log$_{10}$(\dot{a}/a))', size=font_size)
axs[1].set_ylabel('Inspiral rate (Log$_{10}$(\dot{a}/a))', size=font_size)
#plt.ylim(top=1.5)
plt.savefig('inspiral_rate_over_sep_comparison_medians.pdf', bbox_inches='tight', pad_inches=0.02)
#=================================================================================================

#calculate characteristic inspiral rates
print("For 1000 year:")
mean_val = np.mean((All_grads[0]+All_grads[1]+All_grads[2]))
std_val = np.std((All_grads[0]+All_grads[1]+All_grads[2]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads[0]+All_grads[1]+All_grads[2]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("All: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))

mean_val = np.mean((All_grads[0]))
std_val = np.std((All_grads[0]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads[0]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("Bound core frag: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))

mean_val = np.mean((All_grads[1]))
std_val = np.std((All_grads[1]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads[1]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("Unbound core frag: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))

mean_val = np.mean((All_grads[2]))
std_val = np.std((All_grads[2]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads[2]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("Dynamical cap: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))
print("For 10000 year:")

mean_val = np.mean((All_grads_10000[0]+All_grads_10000[1]+All_grads_10000[2]))
std_val = np.std((All_grads_10000[0]+All_grads_10000[1]+All_grads_10000[2]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads_10000[0]+All_grads_10000[1]+All_grads_10000[2]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("All: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))

mean_val = np.mean((All_grads_10000[0]))
std_val = np.std((All_grads_10000[0]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads_10000[0]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("Bound core frag: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))

mean_val = np.mean((All_grads_10000[1]))
std_val = np.std((All_grads_10000[1]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads_10000[1]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("Unbound core frag: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))

mean_val = np.mean((All_grads_10000[2]))
std_val = np.std((All_grads_10000[2]))
upp_val = mean_val + std_val
low_val = mean_val - std_val
median_val = np.median((All_grads_10000[2]))
upp_err = upp_val - median_val
low_err = median_val - low_val
print("Dynamical cap: median:", np.round(median_val, decimals=2), "+/-", [np.round(upp_err, decimals=2), np.round(low_err, decimals=2)], " mean:", np.round(mean_val, decimals=2), "+/-", np.round(std_val, decimals=2))
