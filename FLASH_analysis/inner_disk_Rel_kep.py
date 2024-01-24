import numpy as np
import pickle
import yt
yt.enable_parallelism()
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
from mpi4py.MPI import COMM_WORLD as CW
import my_flash_fields as myf
import my_flash_module as mym
import pickle
import argparse
import os

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
#Ploting parameters
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

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
smooth_window = 40

#---------------------------------------------------
Spin_labels = ['0.20', '0.25', '0.30', '0.35']
Mach_labels = ['0.0', '0.1', '0.2']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
#Define arguments

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, 0.7*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']
plot_it = -1
xmax= 0
ymax = 0
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    axs.flatten()[plot_it].set_title('Mach='+mach_lab, pad=-0.2)
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mean_rad_vel, Min_rad_vel, Mean_rel_vel, Min_rel_vel, Mass_all, Separation = pickle.load(file)
            file.close()
            
            #Calculate smoothed arrays
            T_smoothed = []
            Rel_kep_smoothed = []
            for t_ind in range(len(Time_array)):
                t_start = Time_array[t_ind] - smooth_window/2
                t_end = Time_array[t_ind] + smooth_window/2
                if t_start < 0.0:
                    t_start = 0.0
                if t_end > 10000.0:
                    t_end = 10000.0
                t_start_ind = np.argmin(abs(np.array(Time_array)-t_start))
                t_end_ind = np.argmin(abs(np.array(Time_array)-t_end))
                t_smooth_val = np.mean(np.array(Time_array)[t_start_ind:t_end_ind+1])
                rel_kep_smooth_val = np.mean(np.array(Mean_L)[t_start_ind:t_end_ind+1])
                T_smoothed.append(t_smooth_val)
                Rel_kep_smoothed.append(rel_kep_smooth_val)
            
            axs.flatten()[plot_it].plot(T_smoothed, Rel_kep_smoothed, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
            
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
            if mach_lab == '0.0':
                axs.flatten()[plot_it].set_ylabel('Average $v_\mathrm{kep_relative}$', labelpad=-0.2)
            else:
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels, visible=False)
            if mach_lab != '0.2' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
            if mach_lab == '0.0' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                #plt.setp(xticklabels[-1], visible=False)

axs.flatten()[0].legend(loc='best')
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim([0.58, 1.0])
plt.savefig('Relative_kep_10_au.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, 0.7*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']
plot_it = -1
xmax= 0
ymax = 0
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    axs.flatten()[plot_it].set_title('Mach='+mach_lab, pad=-0.2)
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/5au/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mean_rad_vel, Min_rad_vel, Mean_rel_vel, Min_rel_vel, Mass_all, Separation = pickle.load(file)
            file.close()
            
            #Calculate smoothed arrays
            T_smoothed = []
            Rel_kep_smoothed = []
            for t_ind in range(len(Time_array)):
                t_start = Time_array[t_ind] - smooth_window/2
                t_end = Time_array[t_ind] + smooth_window/2
                if t_start < 0.0:
                    t_start = 0.0
                if t_end > 10000.0:
                    t_end = 10000.0
                t_start_ind = np.argmin(abs(np.array(Time_array)-t_start))
                t_end_ind = np.argmin(abs(np.array(Time_array)-t_end))
                t_smooth_val = np.mean(np.array(Time_array)[t_start_ind:t_end_ind+1])
                rel_kep_smooth_val = np.mean(np.array(Mean_L)[t_start_ind:t_end_ind+1])
                T_smoothed.append(t_smooth_val)
                Rel_kep_smoothed.append(rel_kep_smooth_val)
            
            axs.flatten()[plot_it].plot(T_smoothed, Rel_kep_smoothed, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
            
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
            if mach_lab == '0.0':
                axs.flatten()[plot_it].set_ylabel('Average $v_\mathrm{kep_relative}$', labelpad=-0.2)
            else:
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels, visible=False)
            if mach_lab != '0.2' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
            if mach_lab == '0.0' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                #plt.setp(xticklabels[-1], visible=False)

axs.flatten()[0].legend(loc='best')
axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[2].tick_params(axis='x', direction='in', top=True)
axs.flatten()[2].tick_params(axis='y', direction='in', right=True)
axs.flatten()[2].minorticks_on()
axs.flatten()[2].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim([0.58, 1])
plt.savefig('Relative_kep_pro_5au.pdf', bbox_inches='tight', pad_inches=0.02)

