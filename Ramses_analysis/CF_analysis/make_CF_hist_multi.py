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

subplot_titles = ["3D_all_sinks", "3D_L_limits", "2D_L_limits", "2D_unbound"]
plot_label = ['$L_{max}=120$L$_\odot$', '$L_{max}=55.29$L$_\odot$']
plot_colours = ['b', 'tab:orange']

plot_pickles = [['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/3D_full/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/3D_L_lim/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/3D_L_lim/Tobin_limits/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Bound/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Bound/Tobin_limits/plot_cf_hist.pkl'], ['/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Unbound/plot_cf_hist.pkl', '/groups/astro/rlk/rlk/Analysis_plots/CF_analysis_all/Paper_CF_Hist/2D_obs_Unbound/Tobin_limits/plot_cf_hist.pkl']]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

file_open = open("/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/Perseus_data.pkl", "rb")
bin_centers, CF_per_bin_Tobin_Per, CF_errs_Per = pickle.load(file_open)
file_open.close()

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
            ecolor = 'b'
        else:
            ecolor = plot_colours[pick_it_bot]
        #axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, yerr=CF_err, edgecolor=ecolor, ecolor=ecolor, color=ecolor, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(2-pick_it_bot)})
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_median, yerr=CF_err, width=0.25, alpha=0.5, label=plot_label[pick_it_bot], error_kw = {'elinewidth':(2-pick_it_bot)})
        
        print('plotted', plot_pickles[pick_it_top][pick_it_bot])
        
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].bar(bin_centers, CF_per_bin_Tobin_Per, width=0.25, edgecolor='black', label="Perseus", fill=None, ls='--')
    if pick_it_top == 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].legend(loc='upper right', fontsize=font_size)
    if pick_it_top == 2:
        xticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_xticklabels()
        plt.setp(xticklabels[-1], visible=False)
        yticklabels = axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].get_yticklabels()
        plt.setp(yticklabels[-1], visible=False)
    if pick_it_top > 1:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlabel('Separation (Log$_{10}$(AU))', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_xlim([1, 4])
    if np.remainder(pick_it_top, 2) == 0:
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylabel('Companion Frequency', fontsize=font_size)
        axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].set_ylim([0, 0.2])
        
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='major', labelsize=font_size)
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='both', which='minor', labelsize=font_size)
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='x', direction='in')
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].tick_params(axis='y', direction='in')
    
    axs[int(pick_it_top/2)][np.remainder(pick_it_top, 2)].text((1.25), (0.18), subplot_titles[pick_it_top], zorder=11, fontsize=font_size)
    
plt.savefig('CF_hist_paper.pdf')
