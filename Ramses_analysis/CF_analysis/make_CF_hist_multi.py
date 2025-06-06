import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle

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

subplot_titles = ["3D-Full", "3D-Limits", "2D-Bound", "2D-Unbound"]
plot_label = ['$L_{max}=120$L$_\odot$', '$L_{max}=55$L$_\odot$']
plot_colours = ['tab:blue', 'tab:orange']

plot_pickles = [['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/3D_full/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/3D_L_lim/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/3D_L_lim/Tobin_limits/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Bound/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Bound/Tobin_limits/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Unbound/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Unbound/Tobin_limits/Closest_pairs/plot_cf_hist.pkl']]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data_all.pkl", "rb")
CF_per_bin_all, CF_errs_all, Perseus_sep = pickle.load(file_open)
file_open.close()
CF_errs_all[0][np.array([1, 3])] = 0
CF_errs_all[1][np.array([1, 3])] = 0

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data_66.pkl", "rb")
CF_per_bin_66, CF_errs_66, Perseus_sep = pickle.load(file_open)
file_open.close()

CF_errs_66[0][np.array([1, 3])] = 0
CF_errs_66[1][np.array([1, 3])] = 0

S_bins = np.logspace(1,4,13)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(two_col_width, 0.8*two_col_width), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pick_it_top in range(len(plot_pickles)):
    for pick_it_bot in range(len(plot_pickles[pick_it_top])):
        file = open(plot_pickles[pick_it_top][pick_it_bot], 'rb')
        bin_centers, CF_median, CF_err = pickle.load(file)
        file.close()
        
        if pick_it_top == 0:
            ecolor = 'tab:blue'
        else:
            ecolor = plot_colours[pick_it_bot]
        #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, yerr=CF_err, edgecolor=ecolor, ecolor=ecolor, color=ecolor, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(2-pick_it_bot)})
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(3/(1+pick_it_bot))}, ecolor=ecolor, edgecolor=ecolor)
        
        print('plotted', plot_pickles[pick_it_top][pick_it_bot])
        
    #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_Tobin_Per, yerr=CF_errs_Per, width=0.25, edgecolor='black', label="Perseus (Tobin et al.)", fill=None, ls='--')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_66, yerr=CF_errs_66, width=0.25, edgecolor='black', label="Perseus", fill=None, ls='-')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_all, width=0.25, edgecolor='black', fill=None, ls='--')
    #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_Tobin_Per, width=0.25, edgecolor='black', fill=None, ls='--')
    if pick_it_top == 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].legend(loc='upper right', fontsize=font_size)
    if pick_it_top == 2:
        yticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_yticklabels()
        plt.setp(yticklabels[-1], visible=False)
    if pick_it_top == 3:
        xticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_xticklabels()
        plt.setp(xticklabels[0], visible=False)
    if pick_it_top > 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlabel('Separation (Log$_{10}$(AU))', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlim([1, 4])
    if np.remainder(pick_it_top, 2) == 0:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylabel('Companion Frequency', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylim([0, 0.2])
        
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='major', labelsize=font_size)
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='minor', labelsize=font_size)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='x', direction='in')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='y', direction='in', right=True)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].text((1.03), (0.187), subplot_titles[pick_it_top], zorder=11, fontsize=font_size)
plt.suptitle('SFE=4.0$\pm0.1$%', y=0.91)
    
plt.savefig('CF_hist_paper.pdf', bbox_inches='tight', pad_inches=0.02)

plot_pickles = [['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/3D_full/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/3D_L_lim/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/3D_L_lim/Tobin_limits/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/2D_obs_Bound/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/2D_obs_Bound/Tobin_limits/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/2D_obs_Unbound/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.2/2D_obs_Unbound/Tobin_limits/Closest_pairs/plot_cf_hist.pkl']]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

S_bins = np.logspace(1,4,13)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(two_col_width, 0.8*two_col_width), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pick_it_top in range(len(plot_pickles)):
    for pick_it_bot in range(len(plot_pickles[pick_it_top])):
        file = open(plot_pickles[pick_it_top][pick_it_bot], 'rb')
        bin_centers, CF_median, CF_err = pickle.load(file)
        file.close()
        
        if pick_it_top == 0:
            ecolor = 'grey'
            fcolor = 'grey'
        else:
            ecolor = plot_colours[pick_it_bot]
            fcolor = plot_colours[pick_it_bot]
        #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, yerr=CF_err, edgecolor=ecolor, ecolor=ecolor, color=ecolor, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(2-pick_it_bot)})
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(3/(1+pick_it_bot))}, ecolor=ecolor, edgecolor=ecolor, color=fcolor)
        
        print('plotted', plot_pickles[pick_it_top][pick_it_bot])
        
    #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_Tobin_Per, yerr=CF_errs_Per, width=0.25, edgecolor='black', label="Perseus (Tobin et al.)", fill=None, ls='--')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_66, yerr=CF_errs_66, width=0.25, edgecolor='black', label="Perseus", fill=None, ls='-')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_all, width=0.25, edgecolor='black', fill=None, ls='--')
    if pick_it_top == 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].legend(loc='upper right', fontsize=font_size)
    if pick_it_top == 2:
        yticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_yticklabels()
        plt.setp(yticklabels[-1], visible=False)
    if pick_it_top == 3:
        xticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_xticklabels()
        plt.setp(xticklabels[0], visible=False)
    if pick_it_top > 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlabel('Separation (Log$_{10}$(AU))', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlim([1, 4])
    if np.remainder(pick_it_top, 2) == 0:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylabel('Companion Frequency', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylim([0, 0.2])
        
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='major', labelsize=font_size)
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='minor', labelsize=font_size)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='x', direction='in')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='y', direction='in', right=True)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].text((1.03), (0.187), subplot_titles[pick_it_top], zorder=11, fontsize=font_size)
plt.suptitle('SFE=4.2$\pm0.1$%', y=0.91)
    
plt.savefig('CF_hist_paper_42.pdf', bbox_inches='tight', pad_inches=0.02)

plot_pickles = [['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/3D_full/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/3D_L_lim/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/3D_L_lim/Tobin_limits/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/2D_obs_Bound/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/2D_obs_Bound/Tobin_limits/Closest_pairs/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/2D_obs_Unbound/Closest_pairs/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/SFE_4.1/2D_obs_Unbound/Tobin_limits/Closest_pairs/plot_cf_hist.pkl']]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

S_bins = np.logspace(1,4,13)
bin_centers = (np.log10(S_bins[:-1])+np.log10(S_bins[1:]))/2

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(two_col_width, 0.8*two_col_width), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

for pick_it_top in range(len(plot_pickles)):
    for pick_it_bot in range(len(plot_pickles[pick_it_top])):
        file = open(plot_pickles[pick_it_top][pick_it_bot], 'rb')
        bin_centers, CF_median, CF_err = pickle.load(file)
        file.close()
        
        if pick_it_top == 0:
            ecolor = 'grey'
            fcolor = 'grey'
        else:
            ecolor = plot_colours[pick_it_bot]
            fcolor = plot_colours[pick_it_bot]
        #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, yerr=CF_err, edgecolor=ecolor, ecolor=ecolor, color=ecolor, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(2-pick_it_bot)})
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(3/(1+pick_it_bot))}, ecolor=ecolor, edgecolor=ecolor, color=fcolor)
        
        print('plotted', plot_pickles[pick_it_top][pick_it_bot])
        
    #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_Tobin_Per, yerr=CF_errs_Per, width=0.25, edgecolor='black', label="Perseus (Tobin et al.)", fill=None, ls='--')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_66, yerr=CF_errs_66, width=0.25, edgecolor='black', label="Perseus", fill=None, ls='-')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_all, width=0.25, edgecolor='black', fill=None, ls='--')
    if pick_it_top == 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].legend(loc='upper right', fontsize=font_size)
    if pick_it_top == 2:
        yticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_yticklabels()
        plt.setp(yticklabels[-1], visible=False)
    if pick_it_top == 3:
        xticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_xticklabels()
        plt.setp(xticklabels[0], visible=False)
    if pick_it_top > 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlabel('Separation (Log$_{10}$(AU))', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlim([1, 4])
    if np.remainder(pick_it_top, 2) == 0:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylabel('Companion Frequency', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylim([0, 0.2])
        
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='major', labelsize=font_size)
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='minor', labelsize=font_size)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='x', direction='in')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='y', direction='in', right=True)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].text((1.03), (0.187), subplot_titles[pick_it_top], zorder=11, fontsize=font_size)
plt.suptitle('SFE=4.1$\pm0.1$%', y=0.91)
    
plt.savefig('CF_hist_paper_41.pdf', bbox_inches='tight', pad_inches=0.02)
