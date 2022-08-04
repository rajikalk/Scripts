import numpy as np
import matplotlib.pyplot as plt
import pickle
import glob
import scipy.stats as stats
from scipy.optimize import curve_fit
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
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

def Gaussian(x,scale,mean,sigma):
    return scale*stats.norm.pdf(x, mean, sigma)
    
def Gaussian_cdf(x,scale,mean,sigma):
    return scale*stats.norm.cdf(x, mean, sigma)

def Skewed_Gaussian(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.pdf(x, skew, loc=mean, scale=sigma)
    
def Skewed_Gaussian_cdf(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.cdf(x, skew, loc=mean, scale=sigma)

pickles = ['/Users/reggie/Documents/Papers/Multiplicity_statistics/formation_pathway_1500.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/formation_pathway_3000.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/formation_pathway_3750.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/formation_pathway_4500.pkl', '/Users/reggie/Documents/Papers/Multiplicity_statistics/formation_pathway_6000.pkl',
    '/Users/reggie/Documents/Papers/Multiplicity_statistics/formation_pathway_12000.pkl']

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
label_height = [2.5, 16, 32, 40, 65, 135]
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

S_bins = np.logspace(1,4,13)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
step_centers = np.append(bin_centers, bin_centers[-1]+0.25)
x_fit = np.linspace(0,6,10000)
fit_params = []
guess_params = []

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickles), figsize=(single_col_width, single_col_width*2), sharex=True)#, sharey=True)
iter_range = range(0, len(pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.01)

for pick_it in range(len(pickles)):
    file = open(pickles[pick_it], 'rb')
    pathway_counters, Initial_Seps, Initial_Seps_100000 = pickle.load(file)
    file.close()
    
    core_sep_hist, bins = np.histogram(Initial_Seps[0], S_bins)
    core_delayed_sep_hist, bins = np.histogram(Initial_Seps[1], S_bins)
    capt_sep_hist, bins = np.histogram(Initial_Seps[2], S_bins)
    misc_sep_hist, bins = np.histogram(Initial_Seps[3], S_bins)
    
    #core_sep_hist = np.concatenate((core_sep_hist, np.array([core_sep_hist[-1]])))
    #core_delayed_sep_hist = np.concatenate((core_delayed_sep_hist, np.array([core_delayed_sep_hist[-1]])))
    #capt_sep_hist = np.concatenate((capt_sep_hist, np.array([capt_sep_hist[-1]])))
    
    #plt.clf()
    p1 = axs.flatten()[pick_it].bar(bin_centers, core_sep_hist, width=0.25, color='b')#, hatch='+')
    p2 = axs.flatten()[pick_it].bar(bin_centers, core_delayed_sep_hist, width=0.25, bottom=core_sep_hist, color='m')#, hatch='x')
    p3 = axs.flatten()[pick_it].bar(bin_centers, capt_sep_hist, width=0.25, bottom=(np.array(core_sep_hist)+np.array(core_delayed_sep_hist)), color='r')#, hatch='O')
    if pick_it == 0:
        axs.flatten()[pick_it].legend((p1[0], p2[0], p3[0]), ('Bound core frag.', 'Unbound core frag.', 'Dynamical capture'), loc='upper left', fontsize=font_size, labelspacing=0.2, handletextpad=0.6, borderaxespad=0.3, borderpad=0.2)
        axs.flatten()[pick_it].text((1.1), label_height[pick_it], subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((1.1), label_height[pick_it], subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 3:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[-1], visible=False)
    if pick_it == 5:
        xticklabels =axs.flatten()[pick_it].get_xticklabels()
        plt.setp(xticklabels[-1], visible=False)
    axs.flatten()[pick_it].set_ylim(bottom=0)
    axs.flatten()[pick_it].set_xlim([1,4])
    
    axs.flatten()[pick_it].step(step_centers, np.append(core_sep_hist, np.array([core_sep_hist[-1]])), 'k', where="mid", linewidth=1)
    axs.flatten()[pick_it].step(step_centers, np.append(core_sep_hist+core_delayed_sep_hist, np.array([(core_sep_hist+core_delayed_sep_hist)[-1]])), 'k', where="mid", linewidth=1)
    axs.flatten()[pick_it].step(step_centers, np.append(core_sep_hist+core_delayed_sep_hist+capt_sep_hist, np.array([(core_sep_hist+core_delayed_sep_hist+capt_sep_hist)[-1]])), 'k', where="mid", linewidth=1)
    '''
    plt.clf()
    fig, ax = plt.subplots()
    ax.step(np.log10(S_bins), core_sep_hist, label="Core Fragmentation", linewidth=2, color='b', alpha=0.5, ls='-')
    ax.step(np.log10(S_bins), core_delayed_sep_hist, label="Delayed Core Fragmentation", linewidth=2, color='purple', alpha=0.5, ls='--')
    ax.step(np.log10(S_bins), capt_sep_hist, label="Dynamical Capture", linewidth=2, color='red', alpha=0.5, ls='-.')
    #ax.step(bin_centers, misc_sep_hist, where='post', label="Other", linewidth=2, color='orange', alpha=0.5, ls=':')
    plt.ylim(bottom=0)
    plt.xlim([1,4])
    plt.legend(loc='best')
    '''
    
    #fitting core fragmentation
    core_frag_total_separations = np.log10(Initial_Seps[0] + Initial_Seps[1])
    poisson_err = np.sqrt(len(core_frag_total_separations))
    mean_guess = np.nanmean(core_frag_total_separations)
    median_guess = np.nanmedian(core_frag_total_separations)
    std_guess = np.nanstd(core_frag_total_separations)
    skew_guess = (3*(mean_guess - median_guess))/std_guess
    scale_guess = np.max(core_sep_hist+core_delayed_sep_hist)
    guess_params.append([mean_guess, median_guess, std_guess])
    #popt, pcov = curve_fit(Skewed_Gaussian, bin_centers, (core_sep_hist+core_delayed_sep_hist), [scale_guess, mean_guess, std_guess, skew_guess], sigma=poisson_err*np.ones(np.shape(bin_centers)), absolute_sigma=True)#, bounds=((0, , 0, 0), (1, np.inf, 1, 1))))
    popt, pcov = curve_fit(Gaussian, bin_centers, (core_sep_hist+core_delayed_sep_hist), [scale_guess, mean_guess, std_guess])#, bounds=((0, , 0, 0), (1, np.inf, 1, 1))))
    fit_core = Gaussian(x_fit, *popt)
    #fit_core = Skewed_Gaussian(x_fit, *popt)
    axs.flatten()[pick_it].plot(x_fit, fit_core, 'k-')
    #cdf_fit = Skewed_Gaussian_cdf(x_fit, *popt)
    cdf_fit = Gaussian_cdf(x_fit, *popt)
    cdf_norm = cdf_fit/cdf_fit[-1]
    mode = x_fit[np.argmax(fit_core)]
    median = x_fit[np.argmin(abs(cdf_norm-0.5))]
    err_lower = x_fit[np.argmin(abs(cdf_norm-0.16))]
    err_higher = x_fit[np.argmin(abs(cdf_norm-0.84))]
    print('mean, std, mode, median, +/-:', popt[1], popt[2], mode, median, err_higher, err_lower)
    #fit_params.append([popt[1], median, err_lower, err_higher])
    fit_params.append(popt)
    if pick_it == len(pickles) - 1:
        axs.flatten()[pick_it].set_xlabel('log Separation (au)', labelpad=-0.5, fontsize=font_size)
    axs.flatten()[pick_it].set_ylabel('# Systems', fontsize=font_size)
    axs.flatten()[pick_it].tick_params(axis='both', which='major', labelsize=font_size)
    axs.flatten()[pick_it].tick_params(axis='both', which='minor', labelsize=font_size)
    axs.flatten()[pick_it].tick_params(axis='x', direction='in')
    axs.flatten()[pick_it].tick_params(axis='y', direction='in')

    save_num = pickles[pick_it].split('_')[-1].split('.')[0]
    plt.savefig('initial_sep_dist_'+save_num+'.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.8*single_col_width), sharex=True)
x_val = [1500, 3000, 3750, 4500, 6000, 12000]
y_mean = np.array(fit_params).T[1]#0]
y_std = np.array(fit_params).T[2]

plt.errorbar(x_val, y_mean, y_std)
plt.xlabel('Initial Gas Mass (M$_\odot$)', fontsize=font_size)
plt.ylabel('Core fragmentation scale (log au)', fontsize=font_size)

axs.tick_params(axis='both', which='major', labelsize=font_size, right=True)
axs.tick_params(axis='both', which='minor', labelsize=font_size)
axs.tick_params(axis='x', direction='in')
axs.tick_params(axis='y', direction='in')

plt.savefig('core_fragmentation_scales_fit.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, single_col_width), sharex=True)
x_val = [1500, 3000, 3750, 4500, 6000]# 12000]
y_mean = np.array(guess_params).T[0]
y_median = np.array(guess_params).T[1]
import pdb
pdb.set_trace()
low_bounds = y_mean - np.array(guess_params).T[2]
high_bounds = y_mean + np.array(guess_params).T[2]
y_err_low = y_median - low_bounds
y_err_high = high_bounds - y_median


plt.errorbar(x_val, y_median, np.array([y_err_low, y_err_high]))
plt.xlabel('GMC mass (M$_\odot$)')
plt.ylabel('Core fragmentation scales (log au)', bbox_inches='tight', pad_inches=0.02)
plt.savefig('core_fragmentation_scales_data.pdf')
