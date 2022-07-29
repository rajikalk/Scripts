import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pickle
import matplotlib.patches
import collections
import matplotlib
import matplotlib.ticker
#from mpi4py.MPI import COMM_WORLD as CW

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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-lifetime_thres", "--sys_lifetime_threshold", help="What lifetime threshold do you want to define a stable system", type=float, default=100000.)
    parser.add_argument("-ts", "--timescale", help="Do you want to plot in terms of the of the actual simulation time, or the the time since the formation of the first sink?", type=str, default="first_sink")
    parser.add_argument("-fig_suffix", "--figure_suffix", help="Do you want add a suffix to the figures?", type=str, default="")
    parser.add_argument("-add_hist", "--add_histograms", help="Do you want to add the histograms at the end of the super plots?", type=str, default="False")
    parser.add_argument("-all_sep_evol", "--plot_all_separation_evolution", help="do you want to plot all separation evolution?", type=str, default="True")
    parser.add_argument("-x_field", "--x_field", help="Default for x-axis in the multiplot is time", type=str, default="Time")
    parser.add_argument("-tf", "--text_font", help="what font do you want the text to have?", type=int, default=10)
    parser.add_argument("-plt_key", "--plot_key", help="What dictionary key from superplot_dict do you want to plot?", type=str, default='System_seps')
    parser.add_argument("-smooth", "--smooth_bool", help="Do you want to smooth what you are plotting?", type=str, default='False')
    parser.add_argument("-smooth_window", "--smooth_window_val", help="What big (in yrs) do you want the smoothing window to be?", type=float, default=1000)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

#rank = CW.Get_rank()
#size = CW.Get_size()

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
plot_booleans = [[False, True], [False, False], [True, True], [True, False]]

args = parse_inputs()

#Server pickles
pickle_files = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Visible_star_count/vis_count.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Visible_star_count/vis_count.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Visible_star_count/vis_count.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Visible_star_count/vis_count.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Visible_star_count/vis_count.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Visible_star_count/vis_count.pkl"]

#birth_con_pickles = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Low_Cadence_birth_con/sink_birth_all_delayed_core_frag_cleaned.pkl"]
#Local pickles

plt.clf()
'''
if plot_booleans[rank][0] == True:
    fig, axs = plt.subplots(len(pickle_files), 2, gridspec_kw={'width_ratios': [3, 1]}, figsize=(12, len(pickle_files)*3), sharex='col')
    iter_range = range(0, len(pickle_files)*2, 2)
    args.figure_suffix = args.figure_suffix + '_hist'
else:
'''
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(two_col_width,page_height))
iter_range = range(0, len(pickle_files)*2)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.07)

#if plot_booleans[rank][1] == True:
#    args.figure_suffix = args.figure_suffix + '_thres'

subplot_titles = ["1500M$_\odot$", "3000M$_\odot$", "3750M$_\odot$", "4500M$_\odot$", "6000M$_\odot$", "12000M$_\odot$"]
GMC_mass_arr = [1500, 3000, 3750, 4500, 6000, 12000]
Shrinkage_hist = []
Shrinkage_bins = [-10000, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
Shrinks = []
Lifetime_bins = [0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]
Lifetime_hist = []
S_bins = np.logspace(1,4,13)
CF_hist = np.zeros((len(pickle_files),12)).tolist()
Times_all = []
M_tot_all = []
M_tot_multi_all = []
N_stars_all = []
N_multi_stars_all = []
SFE_all = []
t_max = np.nan
G50_t_max = 0

#================================================================================================
#plot fractions of core fragmentation and dynamical capture
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Visible stars")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_stars_NUL[:sfe_5_ind], label="Total number of stars")
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=55, ls=':', color='k', label='Number of Class 0 in Perseus')
    axs.flatten()[pick_it].axhline(y=92, ls='--', color='k', label='Number of Class 0/I in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Visible_star_comparison_01.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Visible stars")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_stars_UL[:sfe_5_ind], label="Total number of stars")
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=55, ls=':', color='k', label='Number of Class 0 in Perseus')
    axs.flatten()[pick_it].axhline(y=92, ls='--', color='k', label='Number of Class 0/I in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Visible_star_comparison_01_120.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Visible stars")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_stars_UL[:sfe_5_ind], label="Total number of stars")
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=55, ls=':', color='k', label='Number of Class 0 in Perseus')
    axs.flatten()[pick_it].axhline(y=92, ls='--', color='k', label='Number of Class 0/I in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Visible_star_comparison_01_120.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Visible stars")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_tobin[:sfe_5_ind], label="Total number of stars")
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=55, ls=':', color='k', label='Number of Class 0 in Perseus')
    axs.flatten()[pick_it].axhline(y=92, ls='--', color='k', label='Number of Class 0/I in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Visible_star_comparison_009_5529.pdf', bbox_inches='tight', pad_inches=0.02)
    
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Total")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_tobin[:sfe_5_ind], label="L limis only")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_tobin_C0[:sfe_5_ind], label="Class 0/I")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_tobin_C0I[:sfe_5_ind], label="Class 0")
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=55, ls=':', color='k', label='Number of Class 0 in Perseus')
    axs.flatten()[pick_it].axhline(y=92, ls='--', color='k', label='Number of Class 0/I in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Visible_star_comparison_Tobin.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Total number of stars")
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_tobin_C0[:sfe_5_ind], label="Visible stars")
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=55, ls='--', color='k', label='Number of Class 0 in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Visible_star_comparison_Tobin_C0.pdf', bbox_inches='tight', pad_inches=0.02)
    
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(single_col_width, 2.5*single_col_width), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

plot_label = ['$L_{max}=120$L$_\odot$', '$L_{max}=55$L$_\odot$']

for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    Times, SFE, Class_0, Class_0_I, N_total, N_vis_tobin, N_vis_tobin_C0, N_vis_tobin_C0I, N_vis_stars_UL, N_vis_stars_NUL = pickle.load(file)
    file.close()
    
    sfe_5_ind = np.argmin(abs(SFE-0.05))

    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_total[:sfe_5_ind], label="Total number of sinks", color='k')
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_stars_UL[:sfe_5_ind], label='$L_{max}=120$L$_\odot$')
    axs.flatten()[pick_it].plot(SFE[:sfe_5_ind], N_vis_tobin_C0I[:sfe_5_ind], label='$L_{max}=55$L$_\odot$')
    
    axs.flatten()[pick_it].set_ylabel('# Stars', fontsize=font_size)
    axs.flatten()[pick_it].axhline(y=92, ls='--', color='k', label='Number of Class 0/I in Perseus')
    axs.flatten()[pick_it].set_ylim(bottom=0)
    if pick_it == 0:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)
        axs.flatten()[pick_it].set_ylim(top=500)
        axs.flatten()[pick_it].legend(loc='upper left', fontsize=font_size)
        axs.flatten()[pick_it].text((0.03), np.max(N_total[:sfe_5_ind])-0.75*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    else:
        axs.flatten()[pick_it].text((0.002), np.max(N_total[:sfe_5_ind])-0.15*np.max(N_total[:sfe_5_ind]), subplot_titles[pick_it], zorder=11, fontsize=font_size)
    if pick_it == 2:
        yticklabels =axs.flatten()[pick_it].get_yticklabels()
        plt.setp(yticklabels[0], visible=False)

    axs.flatten()[pick_it].set_xlabel('SFE', fontsize=font_size)
    axs.flatten()[pick_it].set_xlim([0, 0.05])
    plt.savefig('Vis_star_paper.pdf', bbox_inches='tight', pad_inches=0.02)

