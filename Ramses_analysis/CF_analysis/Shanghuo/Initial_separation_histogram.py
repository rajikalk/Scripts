import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import collections

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

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

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]

#Formation_pathway = [[18, 11, 17], [42, 25, 39], [86, 29, 84], [82, 50, 113], [97, 45, 102], [17, 8, 21]]
Formation_pathway = []
birth_con_pickles = ['/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/formation_pathway_1500.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/formation_pathway_3000.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/formation_pathway_3750.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/formation_pathway_4500.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/formation_pathway_6000.pkl', '/lustre/astro/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/formation_pathway_12000.pkl']

#birth_con_pickles = ['/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/Max_iter_100/formation_pathway_1500.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/Max_iter_100/formation_pathway_3000.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/Max_iter_100/formation_pathway_3750.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/Max_iter_100/formation_pathway_4500.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/Max_iter_100/formation_pathway_6000.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/SEP_vs_SFE_multi/Max_iter_100/formation_pathway_12000.pkl']

Core_frag_fracs = []
Delayed_core_frag_fracs = []
Dynamical_capt_fracs = []
Unbound_core_fracs = []
Unbound_dynamical_fracs = []
'''
for pathway_numbers in Formation_pathway:
    Total_sys_no = np.sum(pathway_numbers)
    Core_frag_frac = pathway_numbers[0]/Total_sys_no
    Delayed_core_frag_frac = pathway_numbers[1]/Total_sys_no
    Dynamical_capt_frac = pathway_numbers[2]/Total_sys_no
    
    Core_frag_fracs.append(Core_frag_frac)
    Delayed_core_frag_fracs.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs.append(Dynamical_capt_frac)
'''
Initial_Seps_all = []
for birth_con_pickle in birth_con_pickles:
    file = open(birth_con_pickle, 'rb')
    pathway_counters, Initial_Seps, Initial_Seps_100000 = pickle.load(file)
    file.close()
    
    Unbound_tot = np.sum(pathway_counters[1:3])
    unbound_core_frac = pathway_counters[1]/Unbound_tot
    unbound_dynamical_fracs = pathway_counters[2]/Unbound_tot
    Unbound_core_fracs.append(unbound_core_frac)
    Unbound_dynamical_fracs.append(unbound_dynamical_fracs)
    
    Total_sys_no = np.sum(pathway_counters[:3])
    Core_frag_frac = pathway_counters[0]/Total_sys_no
    Delayed_core_frag_frac = pathway_counters[1]/Total_sys_no
    Dynamical_capt_frac = pathway_counters[2]/Total_sys_no
    Initial_Seps_all.append(Initial_Seps)
    Formation_pathway.append([pathway_counters[0], pathway_counters[1], pathway_counters[2]])
    
    Core_frag_fracs.append(Core_frag_frac)
    Delayed_core_frag_fracs.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs.append(Dynamical_capt_frac)
    
    '''
    file = open(birth_con_pickle, 'rb')
    Sink_birth_all = pickle.load(file)
    file.close()
    
    Core_frag_seps = []
    Delayed_core_frag_seps = []
    Dynamical_capt_seps = []
    prev_seps = []
    prev_energies = []
    for key in range(len(Sink_birth_all.keys())):
        if Sink_birth_all[str(key)][2] != 'nan':
            if Sink_birth_all[str(key)][0] == True:
                Core_frag_seps.append(Sink_birth_all[str(key)][3])
            elif Sink_birth_all[str(key)][1] in flatten(eval(Sink_birth_all[str(key)][2])):
                Delayed_core_frag_seps.append(Sink_birth_all[str(key)][3])
            else:
                Dynamical_capt_seps.append(Sink_birth_all[str(key)][3])
            
    remove_dyn = []
    for dyn_sep in Dynamical_capt_seps:
        if dyn_sep in Delayed_core_frag_seps or dyn_sep in Core_frag_seps:
            remove_dyn.append(dyn_sep)
    for dyn_sep in remove_dyn:
        Dynamical_capt_seps.remove(dyn_sep)
    Dynamical_capt_seps = list(set(Dynamical_capt_seps).intersection(set(Dynamical_capt_seps)))
    import pdb
    pdb.set_trace()
    Core_frag_no = len(Core_frag_seps)
    Delayed_core_frag_no = len(Delayed_core_frag_seps)
    Dynamical_capt_no = len(Dynamical_capt_seps) #True number should be 14
    Total_sys_no = Core_frag_no + Delayed_core_frag_no + Dynamical_capt_no
    Core_frag_frac = Core_frag_no/Total_sys_no
    Delayed_core_frag_frac = Delayed_core_frag_no/Total_sys_no
    Dynamical_capt_frac = Dynamical_capt_no/Total_sys_no
    
    Core_frag_fracs.append(Core_frag_frac)
    Delayed_core_frag_fracs.append(Delayed_core_frag_frac)
    Dynamical_capt_fracs.append(Dynamical_capt_frac)
    
    Initial_Seps.append([Core_frag_seps, Delayed_core_frag_seps, Dynamical_capt_seps])
    '''

x_labels = ['Bound core frag.', 'Unbound core frag.', 'Dynamical capture']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
fig, ax = plt.subplots(1, 1, figsize=(single_col_width, 0.8*single_col_width))

ind = np.arange(len(subplot_titles))
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10

p1 = plt.bar(ind, Core_frag_fracs, 0.95, color='b', linewidth=1, edgecolor='k')#, hatch='+'
p2 = plt.bar(ind, Delayed_core_frag_fracs, 0.95, bottom=Core_frag_fracs, color='m', linewidth=1, edgecolor='k')#, hatch='x'
p3 = plt.bar(ind, Dynamical_capt_fracs, 0.95, bottom=(np.array(Delayed_core_frag_fracs)+np.array(Core_frag_fracs)), color='r', linewidth=1, edgecolor='k')#, hatch='O'

plt.xlim([-0.6, 5.6])
plt.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=font_size, right=True)
ax.tick_params(axis='both', which='minor', labelsize=font_size, left=True, right=True, top=False, bottom=False)
plt.xticks(ind, ("1500", "3000", "3750", "4500", "6000", "12000"))
ax.tick_params(which='both', direction='in')
plt.xlabel('Molecular cloud mass (M$_\odot$)', fontsize=font_size, labelpad=-0.5)

plt.legend((p3[0], p2[0], p1[0]), ('Dynamical capture', 'Unbound core frag.', 'Bound core frag.'), loc='upper right', fontsize=font_size)
plt.ylabel('Fraction', fontsize=font_size, labelpad=-0.5)
plt.ylim([0,1])
plt.savefig('formation_pathway.pdf', format='pdf', bbox_inches='tight', pad_inches = 0.02)

#------------------------------------------------------------------------------------------------
x_labels = ['Unbound core frag.', 'Dynamical capture']
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
fig, ax = plt.subplots(1, 1, figsize=(single_col_width, 0.8*single_col_width))

ind = np.arange(len(subplot_titles))
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472
font_size = 10

p1 = plt.bar(ind, Unbound_core_fracs, 0.95, color='m', linewidth=1, edgecolor='k')#, hatch='x'
p2 = plt.bar(ind, Unbound_dynamical_fracs, 0.95, bottom=np.array(Unbound_core_fracs), color='r', linewidth=1, edgecolor='k')#, hatch='O'

plt.xlim([-0.6, 5.6])
plt.minorticks_on()
ax.tick_params(axis='both', which='major', labelsize=font_size, right=True)
ax.tick_params(axis='both', which='minor', labelsize=font_size, left=True, right=True, top=False, bottom=False)
plt.xticks(ind, ("1500", "3000", "3750", "4500", "6000", "12000"))
ax.tick_params(which='both', direction='in')
plt.xlabel('Molecular cloud mass (M$_\odot$)', fontsize=font_size, labelpad=-0.5)

plt.legend((p2[0], p1[0]), ('Dynamical capture', 'Unbound core frag.'), loc='upper right', fontsize=font_size)
plt.ylabel('Fraction', fontsize=font_size, labelpad=-0.5)
plt.ylim([0,1])
plt.savefig('formation_pathway_unbound.png', format='png', bbox_inches='tight', pad_inches = 0.02)

### Make histograms
import scipy.stats as stats
from scipy.optimize import curve_fit
def Gaussian(x,scale,mean,sigma):
    return scale*stats.norm.pdf(x, mean, sigma)
    
def Gaussian_cdf(x,scale,mean,sigma):
    return scale*stats.norm.cdf(x, mean, sigma)

def Skewed_Gaussian(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.pdf(x, skew, loc=mean, scale=sigma)
    
def Skewed_Gaussian_cdf(x, scale, mean, sigma, skew):
    return scale*stats.skewnorm.cdf(x, skew, loc=mean, scale=sigma)

S_bins = np.logspace(1,4,13)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2
#step_centers = np.append(bin_centers, bin_centers[-1]+0.25)
step_centers = np.concatenate(([bin_centers[0]-0.25], bin_centers, [bin_centers[-1]+0.25]))
x_fit = np.linspace(0,6,10000)
fit_params = []
guess_params = []

Shnaghuo_sep = [692.1945689148221,737.4803164795995,1213.8896276468101,467.17906905588245,
 610.2815890074502,326.753154960137,1405.8571018041273,479.7568722618951,
 756.8869692206754,530.3654108966807,520.1109734513876,1176.7210281092603,
583.9974996285858]
Shanghuo_hist, bins = np.histogram(Shnaghuo_sep, S_bins)

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=len(birth_con_pickles), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey='row')
iter_range = range(0, len(birth_con_pickles))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pick_it in range(len(Initial_Seps_all)):
    core_sep_hist, bins = np.histogram(Initial_Seps_all[pick_it][0], S_bins)
    core_delayed_sep_hist, bins = np.histogram(Initial_Seps_all[pick_it][1], S_bins)
    capt_sep_hist, bins = np.histogram(Initial_Seps_all[pick_it][2], S_bins)
    
    p1 = axs[pick_it][0].bar(bin_centers, core_sep_hist, width=0.25, color='b', label='Bound core frag.')#, hatch='+')
    p2 = axs[pick_it][0].bar(bin_centers, core_delayed_sep_hist, width=0.25, bottom=core_sep_hist, color='m', label='Unbound core frag.')#, hatch='x')
    p3 = axs[pick_it][1].bar(bin_centers, capt_sep_hist, width=0.25, color='r', label='Dynamical capture')#, hatch='O')
    if pick_it == 0:
        axs[pick_it][0].legend(loc='upper left', fontsize=font_size, labelspacing=0.2, handletextpad=0.6, borderaxespad=0.3, borderpad=0.2)
        axs[pick_it][1].legend(loc='upper left', fontsize=font_size, labelspacing=0.2, handletextpad=0.6, borderaxespad=0.3, borderpad=0.2)
        axs[pick_it][0].text((1.1), np.max(core_sep_hist+core_delayed_sep_hist)-0.5*np.max(core_sep_hist+core_delayed_sep_hist), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs[pick_it][0].text((1.1), np.max(core_sep_hist+core_delayed_sep_hist)-0.15*np.max(core_sep_hist+core_delayed_sep_hist), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels = axs[pick_it][0].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
    if pick_it == 5:
        xticklabels = axs[pick_it][0].get_xticklabels()
        plt.setp(xticklabels[-1], visible=False)
    axs[pick_it][0].set_ylim(bottom=0)
    axs[pick_it][1].set_ylim(bottom=0)
    axs[pick_it][0].set_xlim([1,4])
    axs[pick_it][1].set_xlim([1,4])
    
    if pick_it == 5:
        Core_seps_2D = (1/np.sqrt(2))*np.array(Initial_Seps_all[pick_it][0])
        core_sep_hist, bins = np.histogram(Core_seps_2D, S_bins)
        
        Unbound_core_seps_2D = (1/np.sqrt(2))*np.array(Initial_Seps_all[pick_it][1])
        core_delayed_sep_hist, bins = np.histogram(Unbound_core_seps_2D, S_bins)
        
        total_hist = core_sep_hist + core_delayed_sep_hist
        total_hist_normalised = total_hist/(np.sum(total_hist))
        core_sep_hist_normalised = core_sep_hist/(np.sum(total_hist))
        core_delayed_sep_hist_normalised = core_delayed_sep_hist/(np.sum(total_hist))
        
        Shanghuo_hist_normalised = Shanghuo_hist/np.sum(Shanghuo_hist)
        
        plt.clf()
        plt.cla()
        plt.bar(bin_centers, core_sep_hist_normalised, width=0.25, color='b', label='Bound core frag. (Kuruwita \& Haug{\o}lle, 2023)')
        plt.bar(bin_centers, core_delayed_sep_hist_normalised, width=0.25, bottom=core_sep_hist_normalised, color='m', label='Unbound core frag. (Kuruwita \& Haug{\o}lle, 2023)')
        plt.bar(bin_centers, Shanghuo_hist_normalised, width=0.25, color='b', label='This Paper')
        
        plt.xlabel('Separation (Log$_{10}$(AU))', labelpad=-0.5, fontsize=font_size)
        plt.ylabel('Normalised # System', fontsize=font_size)
        plt.legend(loc='best')
        
        plt.savefig('Shanghuo.pdf', bbox_inches='tight', pad_inches=0.02)
        import pdb
        pdb.set_trace()
    
    axs[pick_it][0].step(step_centers, np.concatenate(([core_sep_hist[0]], core_sep_hist, [core_sep_hist[-1]])), 'k', where="mid", linewidth=1)
    axs[pick_it][0].step(step_centers, np.concatenate(([(core_sep_hist+core_delayed_sep_hist)[0]], core_sep_hist+core_delayed_sep_hist, [(core_sep_hist+core_delayed_sep_hist)[-1]])), 'k', where="mid", linewidth=1)
    axs[pick_it][1].step(step_centers, np.concatenate(([capt_sep_hist[0]], capt_sep_hist, [capt_sep_hist[-1]])), 'k', where="mid", linewidth=1)
    
    #axs[pick_it][0].step(step_centers, np.append(core_sep_hist[0], np.append(core_sep_hist, np.array([core_sep_hist[-1]]))), 'k', where="mid", linewidth=1)
    #axs[pick_it][0].step(step_centers, np.append((core_sep_hist+core_delayed_sep_hist)[0], np.append(core_sep_hist+core_delayed_sep_hist, np.array([(core_sep_hist+core_delayed_sep_hist)[-1]]))), 'k', where="mid", linewidth=1)
    #axs[pick_it][1].step(step_centers, np.append(capt_sep_hist[0], np.append(capt_sep_hist, np.array([(capt_sep_hist)[-1]]))), 'k', where="mid", linewidth=1)
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
    core_frag_total_separations = np.log10(Initial_Seps_all[pick_it][0] + Initial_Seps_all[pick_it][1])
    poisson_err = np.sqrt(len(core_frag_total_separations))
    mean_guess = np.nanmean(core_frag_total_separations)
    median_guess = np.nanmedian(core_frag_total_separations)
    std_guess = np.nanstd(core_frag_total_separations)
    skew_guess = (3*(mean_guess - median_guess))/std_guess
    scale_guess = np.max(core_sep_hist+core_delayed_sep_hist)
    guess_params.append([mean_guess, median_guess, std_guess])
    #popt, pcov = curve_fit(Skewed_Gaussian, bin_centers, (core_sep_hist+core_delayed_sep_hist), [scale_guess, mean_guess, std_guess, skew_guess], sigma=poisson_err*np.ones(np.shape(bin_centers)), absolute_sigma=True)#, bounds=((0, , 0, 0), (1, np.inf, 1, 1))))
    popt, pcov = curve_fit(Gaussian, bin_centers[6:], (core_sep_hist+core_delayed_sep_hist)[6:], [scale_guess, mean_guess, std_guess])#, bounds=((0, , 0, 0), (1, np.inf, 1, 1))))
    fit_core = Gaussian(x_fit, *popt)
    #fit_core = Skewed_Gaussian(x_fit, *popt)
    axs[pick_it][0].plot(x_fit, fit_core, 'k-')
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
    if pick_it == len(Initial_Seps_all) - 1:
        axs[pick_it][0].set_xlabel('Separation (Log$_{10}$(AU))', labelpad=-0.5, fontsize=font_size)
        axs[pick_it][1].set_xlabel('Separation (Log$_{10}$(AU))', labelpad=-0.5, fontsize=font_size)
    axs[pick_it][0].set_ylabel('# Systems', fontsize=font_size)
    axs[pick_it][0].tick_params(axis='both', which='major', labelsize=font_size)
    axs[pick_it][0].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pick_it][1].tick_params(axis='both', which='major', labelsize=font_size)
    axs[pick_it][1].tick_params(axis='both', which='minor', labelsize=font_size)
    axs[pick_it][0].tick_params(axis='x', direction='in')
    axs[pick_it][0].tick_params(axis='y', direction='in')
    axs[pick_it][1].tick_params(axis='x', direction='in')
    axs[pick_it][1].tick_params(axis='y', direction='in')

    save_num = subplot_titles[pick_it].split('M')[0]
    plt.savefig('initial_sep_dist_'+save_num+'.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, 0.8*single_col_width), sharex=True)
x_val = [1500, 3000, 3750, 4500, 6000, 12000]
y_mean = np.array(fit_params).T[1]#0]
y_std = np.array(fit_params).T[2]

y_lin = 10**y_mean
y_lower = 10**(y_mean-y_std)
y_upper = 10**(y_mean+y_std)
y_lower_err = y_lin - y_lower
y_upper_err = y_upper - y_lin
y_std_linear_upper = 10**(y_mean+y_std)
y_std_linear_lower = 10**(y_mean-y_std)
y_err = np.array([y_std_linear_lower, y_std_linear_upper])

#import pdb
#pdb.set_trace()
#plt.errorbar(x_val, 10**y_mean, )

#plt.plot(x_val, 10**(y_mean))
mass_arr = np.linspace(1500, 12000, 1000)
frag_fit = np.log10(1200*100/np.sqrt(mass_arr))
plt.errorbar(x_val, y_mean, y_std)
plt.plot(mass_arr, frag_fit, 'k--')
#plt.yscale('log')
plt.xlabel('Molecular cloud mass (M$_\odot$)', fontsize=font_size)
#plt.ylabel('Core fragmentation scale (Log$_{10}$(AU))', fontsize=font_size)
plt.ylabel('Core fragmentation scale (log$_{10}$AU)', fontsize=font_size)

axs.tick_params(axis='both', which='major', labelsize=font_size, right=True)
axs.tick_params(axis='both', which='minor', labelsize=font_size)
axs.tick_params(axis='x', direction='in')
axs.tick_params(axis='y', direction='in')

plt.savefig('core_fragmentation_scales_mean.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=1, figsize=(single_col_width, single_col_width), sharex=True)
x_val = [1500, 3000, 3750, 4500, 6000, 12000]
y_mean = np.array(guess_params).T[0]
y_median = np.array(guess_params).T[1]

low_bounds = y_mean - np.array(guess_params).T[2]
high_bounds = y_mean + np.array(guess_params).T[2]
y_err_low = y_median - low_bounds
y_err_high = high_bounds - y_median

plt.errorbar(x_val, y_median, np.array([y_err_low, y_err_high]))
plt.xlabel('GMC mass (M$_\odot$)')
plt.ylabel('Core fragmentation scales (Log$_{10}$ AU)')
plt.savefig('core_fragmentation_scales_median.pdf', bbox_inches='tight', pad_inches=0.02)
