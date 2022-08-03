import numpy as np
import matplotlib.pyplot as plt
import yt
import glob
import sys
import collections
import matplotlib as mpl
import pickle
import os
import matplotlib.gridspec as gridspec
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

fit_pickles = ['/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G50/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G100/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G125/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G150/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G200/chi_squared_fit.pkl', '/groups/astro/rlk/rlk/Analysis_plots/Fit_over_time/Obs_bound/G400/chi_squared_fit.pkl']

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
line_styles = [':', (0, (3, 1, 1, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 1, 1, 1)), '--', '-']

fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(single_col_width,1.5*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(pickle_files)*2)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.07)

for fit_pick in range(len(fit_pickles)):
    file = open(fit_pickles[fit_pick], 'rb')
    SFE, reduced_chi_square_selected, reduced_chi_square_tobin, p_values = pickle.load(file)
    file.close()
    
    smoothed_chi = []
    smoothed_upp = []
    smoothed_low = []
    smooth_window = 0.0001
    for SFE_it in range(len(SFE)):
        low_SFE = SFE[SFE_it] - smooth_window
        if low_SFE < 0:
            low_SFE = 0
        high_SFE = SFE[SFE_it] + smooth_window
        if high_SFE > 0.05:
            high_SFE = 5
        low_SFE_it = np.argmin(abs(SFE-low_SFE))
        high_SFE_it = np.argmin(abs(SFE-high_SFE))
        mean_chi = np.mean(reduced_chi_square_tobin[low_SFE_it:high_SFE_it])
        median_chi = np.median(reduced_chi_square_tobin[low_SFE_it:high_SFE_it])
        std_chi = np.std(reduced_chi_square_tobin[low_SFE_it:high_SFE_it])
        err_upper = mean_chi+std_chi
        err_lower = mean_chi-std_chi
        smoothed_chi.append(median_chi)
        smoothed_upp.append(err_upper)
        smoothed_low.append(err_lower)
    axs[0].semilogy(SFE, smoothed_chi, label=subplot_titles[fit_pick], color=colors[fit_pick], ls=line_styles[fit_pick])
    axs[0].fill_between(SFE, smoothed_low, smoothed_upp, alpha=0.2)
    
    smoothed_chi = []
    smoothed_upp = []
    smoothed_low = []
    smooth_window = 0.0001
    for SFE_it in range(len(SFE)):
        low_SFE = SFE[SFE_it] - smooth_window
        if low_SFE < 0:
            low_SFE = 0
        high_SFE = SFE[SFE_it] + smooth_window
        if high_SFE > 0.05:
            high_SFE = 5
        low_SFE_it = np.argmin(abs(SFE-low_SFE))
        high_SFE_it = np.argmin(abs(SFE-high_SFE))
        mean_chi = np.mean(reduced_chi_square_selected[low_SFE_it:high_SFE_it])
        median_chi = np.median(reduced_chi_square_selected[low_SFE_it:high_SFE_it])
        std_chi = np.std(reduced_chi_square_selected[low_SFE_it:high_SFE_it])
        err_upper = mean_chi+std_chi
        err_lower = mean_chi-std_chi
        smoothed_chi.append(median_chi)
        smoothed_upp.append(err_upper)
        smoothed_low.append(err_lower)
    axs[1].semilogy(SFE, smoothed_chi, label=subplot_titles[fit_pick], color=colors[fit_pick], ls=line_styles[fit_pick])
    axs[1].fill_between(SFE, smoothed_low, smoothed_upp, alpha=0.2)
    
axs[0].legend(loc='upper right', ncol=2, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
axs[1].xlabel('SFE')
axs[0].ylabel("Fit (<$\chi^2$>)", labelpad=-0.2)
axs[1].ylabel("Fit (<$\chi^2$>)", labelpad=-0.2)
axs[0].text((0.02), (0.6), "All bins", zorder=11, size=args.text_font)
axs[1].text((0.02), (0.6), "Selected bins", zorder=11, size=args.text_font)
axs[0].xlim([0.01, 0.05])
axs[0].ylim(top=100)
plt.savefig("fit_multi.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
plt.figure(figsize=(single_col_width,0.7*single_col_width))

for fit_pick in range(len(fit_pickles)):
    file = open(fit_pickles[fit_pick], 'rb')
    SFE, reduced_chi_square_selected, reduced_chi_square_tobin, p_values = pickle.load(file)
    file.close()
    
    smoothed_chi = []
    smoothed_upp = []
    smoothed_low = []
    smooth_window = 0.0001
    for SFE_it in range(len(SFE)):
        low_SFE = SFE[SFE_it] - smooth_window
        if low_SFE < 0:
            low_SFE = 0
        high_SFE = SFE[SFE_it] + smooth_window
        if high_SFE > 0.05:
            high_SFE = 5
        low_SFE_it = np.argmin(abs(SFE-low_SFE))
        high_SFE_it = np.argmin(abs(SFE-high_SFE))
        mean_chi = np.mean(reduced_chi_square_selected[low_SFE_it:high_SFE_it])
        median_chi = np.median(reduced_chi_square_selected[low_SFE_it:high_SFE_it])
        std_chi = np.std(reduced_chi_square_selected[low_SFE_it:high_SFE_it])
        err_upper = mean_chi+std_chi
        err_lower = mean_chi-std_chi
        smoothed_chi.append(median_chi)
        smoothed_upp.append(err_upper)
        smoothed_low.append(err_lower)
    plt.semilogy(SFE, smoothed_chi, label=subplot_titles[fit_pick], color=colors[fit_pick], ls=line_styles[fit_pick])
    plt.fill_between(SFE, smoothed_low, smoothed_upp, alpha=0.2)
    
plt.legend(loc='upper right', ncol=2, fontsize=font_size, labelspacing=0.1, handletextpad=0.2, borderaxespad=0.2, borderpad=0.2, columnspacing=0.3)
plt.xlabel('SFE')
plt.ylabel("Fit (<$\chi^2$>)", labelpad=-0.2)
plt.xlim([0.01, 0.05])
plt.ylim(top=100)
plt.savefig("fit_selected_bins_multi.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02)
