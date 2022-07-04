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

#grad_pickles = ['/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G50/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G100/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G125/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G150/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G200/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G400/grad_pickle.pkl']

grad_pickles = ['/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G50/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G100/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G125/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G150/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G200/1000_yr/grad_pickle.pkl', '/lustre/astro/rlk/Analysis_plots/Pathway_evolution/G400/1000_yr/grad_pickle.pkl']

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
for grad_it in range(len(grad_pickles)):
    file = open(grad_pickles[grad_it], 'rb')
    Initial_gradients, Initial_gradients_1000, Initial_gradients_10000, Initial_gradients_100000, Grad_1e4, Grad_1e3 = pickle.load(file)
    file.close()

    core_inspiral_inds = np.where(np.array(Initial_gradients_10000[0]) < 0)[0]
    delayed_core_inspiral_inds = np.where(np.array(Initial_gradients_10000[1]) < 0)[0]
    capt_inspiral_inds = np.where(np.array(Initial_gradients_10000[2]) < 0)[0]
    other_inspiral_inds = np.where(np.array(Initial_gradients_10000[3]) < 0)[0]
    
    Mean_grads[0].append(np.mean(np.log10(-1*np.array(Initial_gradients_10000[0]).T[0][core_inspiral_inds])))
    Median_grads[0].append(np.median(np.log10(-1*np.array(Initial_gradients_10000[0]).T[0][core_inspiral_inds])))
    Std_grads[0].append(np.std(np.log10(-1*np.array(Initial_gradients_10000[0]).T[0][core_inspiral_inds])))
    
    Mean_grads[1].append(np.mean(np.log10(-1*np.array(Initial_gradients_10000[1]).T[0][delayed_core_inspiral_inds])))
    Median_grads[1].append(np.median(np.log10(-1*np.array(Initial_gradients_10000[1]).T[0][delayed_core_inspiral_inds])))
    Std_grads[1].append(np.std(np.log10(-1*np.array(Initial_gradients_10000[1]).T[0][delayed_core_inspiral_inds])))
    
    Mean_grads[2].append(np.mean(np.log10(-1*np.array(Initial_gradients_10000[2]).T[0][capt_inspiral_inds])))
    Median_grads[2].append(np.median(np.log10(-1*np.array(Initial_gradients_10000[2]).T[0][capt_inspiral_inds])))
    Std_grads[2].append(np.std(np.log10(-1*np.array(Initial_gradients_10000[2]).T[0][capt_inspiral_inds])))
    
    Mean_grads[3].append(np.mean(np.log10(-1*np.array(Initial_gradients_10000[3]).T[0][other_inspiral_inds])))
    Median_grads[3].append(np.median(np.log10(-1*np.array(Initial_gradients_10000[3]).T[0][other_inspiral_inds])))
    Std_grads[3].append(np.std(np.log10(-1*np.array(Initial_gradients_10000[3]).T[0][other_inspiral_inds])))
    
    #Plotting initial gradients
    #grad_bins = np.concatenate((-1*np.logspace(4,-3,15), np.array([0, 1.e10])))
    grad_hist_core, grad_bins = np.histogram(Initial_gradients_10000[0], bins=grad_bins)
    grad_hist_core_delayed, grad_bins = np.histogram(Initial_gradients_10000[1], bins=grad_bins)
    grad_hist_capt, grad_bins = np.histogram(Initial_gradients_10000[2], bins=grad_bins)
    grad_hist_misc, grad_bins = np.histogram(Initial_gradients_10000[3], bins=grad_bins)

    x_range = np.arange(len(grad_hist_core)+1)

    grad_hist_core_norm = grad_hist_core/np.sum(grad_hist_core)
    grad_hist_core_delayed_norm = grad_hist_core_delayed/np.sum(grad_hist_core_delayed)
    grad_hist_capt_norm = grad_hist_capt/np.sum(grad_hist_capt)
    grad_hist_misc_norm = grad_hist_misc/np.sum(grad_hist_misc)
    
    grad_hist_core_norm = np.concatenate((grad_hist_core_norm, np.array([grad_hist_core_norm[-1]])))
    grad_hist_core_delayed_norm = np.concatenate((grad_hist_core_delayed_norm, np.array([grad_hist_core_delayed_norm[-1]])))
    grad_hist_capt_norm = np.concatenate((grad_hist_capt_norm, np.array([grad_hist_capt_norm[-1]])))
    grad_hist_misc_norm = np.concatenate((grad_hist_misc_norm, np.array([grad_hist_misc_norm[-1]])))

    axs[grad_it].step(x_range, grad_hist_core_norm, where='post', label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
    axs[grad_it].step(x_range, grad_hist_core_delayed_norm, where='post', label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
    axs[grad_it].step(x_range, grad_hist_capt_norm, where='post', label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
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
plt.ylabel('Log Inspiral rate (au/yr)')
plt.savefig('inspiral_rate_comparison_means')

Core_bounds = [np.array(Mean_grads[0])-np.array(Std_grads[0]), np.array(Mean_grads[0])+np.array(Std_grads[0])]
Delayed_core_bounds = [np.array(Mean_grads[1])-np.array(Std_grads[1]), np.array(Mean_grads[1])+np.array(Std_grads[1])]
Capt_bounds = [np.array(Mean_grads[2])-np.array(Std_grads[2]), np.array(Mean_grads[2])+np.array(Std_grads[2])]
Other_bounds = [np.array(Mean_grads[3])-np.array(Std_grads[3]), np.array(Mean_grads[3])+np.array(Std_grads[3])]

Core_err = [np.array(Median_grads[0]) - Core_bounds[0], Core_bounds[1] - np.array(Median_grads[0])]
Delayed_core_err = [np.array(Median_grads[1]) - Delayed_core_bounds[0], Delayed_core_bounds[1] - np.array(Median_grads[1])]
Capt_err = [np.array(Median_grads[2]) - Capt_bounds[0], Capt_bounds[1] - np.array(Median_grads[2])]
Other_err = [np.array(Median_grads[3]) - Other_bounds[0], Other_bounds[1] - np.array(Median_grads[3])]

plt.clf()
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(single_col_width, 1.25*single_col_width)
plt.errorbar(masses, Median_grads[0], yerr=Core_err, label='Core Fragmentation', color='b')
plt.errorbar(np.array(masses)+10, Median_grads[1], yerr=Delayed_core_err, label='Delayed Core Fragmentation', color='purple')
plt.errorbar(np.array(masses)+20, Median_grads[2], yerr=Capt_err, label='Dynamical Capture', color='r')
plt.errorbar(np.array(masses)+30, Median_grads[3], yerr=Other_err, label='Other', color='orange')
plt.legend()
plt.xlabel('Gas Mass')
plt.ylabel('Log Inspiral rate (au/yr)')
plt.savefig('inspiral_rate_comparison_medians', bbox_inches='tight')
