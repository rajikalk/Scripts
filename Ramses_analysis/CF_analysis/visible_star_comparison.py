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
pickle_files = ["/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G50/Visible_star_count/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G100/Visible_star_count/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G125/Visible_star_count/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G150/Visible_star_count/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G200/Visible_star_count/means_superplot.pkl", "/groups/astro/rlk/rlk/Analysis_plots/Superplot_pickles_entire_sim/G400/Visible_star_count/means_superplot.pkl"]

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
fig, axs = plt.subplots(ncols=1, nrows=len(pickle_files), figsize=(two_col_width, page_height), sharex=True)#, sharey=True)
iter_range = range(0, len(pickle_files))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.02)

CF_hist = np.zeros((len(pickle_files),12)).tolist()

Formation_pathway = []
Not_plotted_sinks = [[],[],[],[],[],[]]
for pick_it in iter_range:
    file_it = pick_it
    file = open(pickle_files[file_it], 'rb')
    superplot_dict, Sink_bound_birth, Sink_formation_times, means_dict, Lifetimes_sys, Sep_maxs, Sep_mins, Initial_Seps, Final_seps = pickle.load(file)
    file.close()
    
    import pdb
    pdb.set_trace()
    
    axs.flatten()[pick_it].semilogy(Time_adjusted_formation[:t_max_ind], superplot_dict[args.plot_key][time_key][:t_max_ind], alpha=alpha, color='k', rasterized=True)
