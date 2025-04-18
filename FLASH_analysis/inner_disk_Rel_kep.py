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
smooth_window = 300

#---------------------------------------------------
Spin_labels = ['0.20', '0.25', '0.30', '0.35']
Mach_labels = ['0.0', '0.1', '0.2']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
              
spin_up_start = [[4500, 4500, 5250, 6750], [np.nan, 4000, np.nan, 5750], [3750, 3500, 4250, 4250]]
spin_up_end = [[5400, 5450, 5750, 8000], [np.nan, 4500, np.nan, 7250], [5750, 5300, 5500, 5750]]
peak_times = [[6288.53350698, 5337.46225949, 5961.44966664, 7499.56279945],
 [np.nan, 4396.31866809, np.nan, 7239.5115915],
 [5684.37286105, 5479.48842751, 5771.10385454, 5758.7326983]]
#Define arguments

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, 0.7*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']
plot_it = -1
xmax= 0
ymin = 0.78
ymax = 1
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    axs.flatten()[plot_it].set_title('Mach='+mach_lab, pad=-0.2)
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
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
                rel_kep_smooth_val = np.mean(np.array(Mean_rel_kep)[t_start_ind:t_end_ind+1])
                T_smoothed.append(t_smooth_val)
                Rel_kep_smoothed.append(rel_kep_smooth_val)
            
            axs.flatten()[plot_it].plot(Time_array, Mean_rel_kep, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.1, linewidth=1)
            axs.flatten()[plot_it].plot(T_smoothed, Rel_kep_smoothed, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75, linewidth=1)
            
            highlight_start_time = spin_up_start[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
            if np.isnan(highlight_start_time) == False:
                highlight_end_time = spin_up_end[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
                highlight_start_ind = np.argmin(abs(np.array(Time_array)-highlight_start_time))
                highlight_end_ind = np.argmin(abs(np.array(Time_array)-highlight_end_time))
                highlight_min = np.min(Rel_kep_smoothed[highlight_start_ind:highlight_end_ind])
                highlight_max = np.max(Rel_kep_smoothed[highlight_start_ind:highlight_end_ind])
                axs.flatten()[plot_it].axvspan(highlight_start_time, highlight_end_time, ymin=(highlight_min-ymin)/(ymax-ymin), ymax=(highlight_max-ymin)/(ymax-ymin), alpha=0.30, facecolor=colors[Spin_labels.index(spin_lab)])
            
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
            if mach_lab == '0.0':
                axs.flatten()[plot_it].set_ylabel('Relative $v_{\mathrm{kep}}$', labelpad=-0.2)
            else:
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels, visible=False)
            if mach_lab != '0.2' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
            if mach_lab == '0.0' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                #plt.setp(xticklabels[-1], visible=False)

axs.flatten()[1].legend(loc='best')
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
axs.flatten()[plot_it-1].set_ylim([ymin, ymax])
plt.savefig('Relative_kep_10_au.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, 0.7*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']
smooth_window = 300
plot_it = -1
xmax= 0
ymin = 0.6
ymax = 0.75
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    axs.flatten()[plot_it].set_title('Mach='+mach_lab, pad=-0.2)
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/5au/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
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
                rel_kep_smooth_val = np.mean(np.array(Mean_rel_kep)[t_start_ind:t_end_ind+1])
                T_smoothed.append(t_smooth_val)
                Rel_kep_smoothed.append(rel_kep_smooth_val)
            
            axs.flatten()[plot_it].plot(Time_array, Mean_rel_kep, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.1, linewidth=1)
            axs.flatten()[plot_it].plot(T_smoothed, Rel_kep_smoothed, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75, linewidth=1)
            
            highlight_start_time = spin_up_start[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
            if np.isnan(highlight_start_time) == False:
                highlight_end_time = spin_up_end[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
                highlight_start_ind = np.argmin(abs(np.array(Time_array)-highlight_start_time))
                highlight_end_ind = np.argmin(abs(np.array(Time_array)-highlight_end_time))
                highlight_min = np.min(Rel_kep_smoothed[highlight_start_ind:highlight_end_ind])
                highlight_max = np.max(Rel_kep_smoothed[highlight_start_ind:highlight_end_ind])
                axs.flatten()[plot_it].axvspan(highlight_start_time, highlight_end_time, ymin=(highlight_min-ymin)/(ymax-ymin), ymax=(highlight_max-ymin)/(ymax-ymin), alpha=0.30, facecolor=colors[Spin_labels.index(spin_lab)])
            
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
            if mach_lab == '0.0':
                axs.flatten()[plot_it].set_ylabel('Relative $v_{\mathrm{kep}}$', labelpad=-0.2)
            else:
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels, visible=False)
            if mach_lab != '0.2' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
            if mach_lab == '0.0' and spin_lab == '0.35':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                #plt.setp(xticklabels[-1], visible=False)

axs.flatten()[1].legend(loc='best')
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
axs.flatten()[plot_it-1].set_ylim([ymin, ymax])
plt.savefig('Relative_kep_pro_5au.pdf', bbox_inches='tight', pad_inches=0.02)

#==============================================================================
plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = -1
xmax= 0
ymax = 1
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        plot_it = plot_it + 1
        axs.flatten()[plot_it].grid()

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
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
                rel_kep_smooth_val = np.mean(np.array(Mean_rel_kep)[t_start_ind:t_end_ind+1])
                T_smoothed.append(t_smooth_val)
                Rel_kep_smoothed.append(rel_kep_smooth_val)
            
            axs.flatten()[plot_it].plot(T_smoothed, Rel_kep_smoothed, label='10au', linestyle='-', color=colors[Spin_labels.index(spin_lab)], alpha=0.75, linewidth=1)
        else:
            print("Couldn't open", pickle_file)
            
        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/5au/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
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
                rel_kep_smooth_val = np.mean(np.array(Mean_rel_kep)[t_start_ind:t_end_ind+1])
                T_smoothed.append(t_smooth_val)
                Rel_kep_smoothed.append(rel_kep_smooth_val)
            
            axs.flatten()[plot_it].plot(T_smoothed, Rel_kep_smoothed, label='5au', linestyle='--', color=colors[Spin_labels.index(spin_lab)], alpha=0.75, linewidth=1)
            axs.flatten()[plot_it].axvspan(spin_up_start[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)], spin_up_end[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)], alpha=0.30, facecolor=colors[Spin_labels.index(spin_lab)], edgecolor=None)
        else:
            print("Couldn't open", pickle_file)
            
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')
        
        if spin_lab == '0.20':
            axs.flatten()[plot_it].set_title('Mach='+mach_lab)
        if spin_lab == '0.20' and mach_lab == '0.0':
            axs.flatten()[plot_it].legend(loc='best')
    if spin_lab == '0.35':
        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='x', direction='in')
    axs.flatten()[plot_it].tick_params(axis='y', direction='in')

axs.flatten()[plot_it].set_xlim([0, 10000])
axs.flatten()[plot_it].set_ylim([0.58,1])
plt.savefig('Rel_kep_10_and_5_au.pdf', bbox_inches='tight')

#======================================================================================================================
#Plotting disk mass and specific angular momentum

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, 0.7*single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

line_styles = ['-', '--', '-.', ':']
plot_it = -1
xmax= 0
ymin = 0
ymax = 0.025
for mach_lab in Mach_labels:
    plot_it = plot_it + 1
    axs.flatten()[plot_it].set_title('Mach='+mach_lab, pad=-0.2)
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/Relative_keplerian_velocity/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
            file.close()
            Disk_mass = np.array(Total_mass)/1.98841586e+33

            axs.flatten()[plot_it].plot(Time_array, Disk_mass, label='$\Omega t_{ff}$='+spin_lab, ls=line_styles[Spin_labels.index(spin_lab)], alpha=0.75, color=colors[Spin_labels.index(spin_lab)], linewidth=1)
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
            
            highlight_start_time = spin_up_start[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
            if np.isnan(highlight_start_time) == False:
                highlight_end_time = spin_up_end[Mach_labels.index(mach_lab)][Spin_labels.index(spin_lab)]
                highlight_start_ind = np.argmin(abs(np.array(Time_array)-highlight_start_time))
                highlight_end_ind = np.argmin(abs(np.array(Time_array)-highlight_end_time))
                highlight_min = np.min(Disk_mass[highlight_start_ind:highlight_end_ind])
                highlight_max = np.max(Disk_mass[highlight_start_ind:highlight_end_ind])
                axs.flatten()[plot_it].axvspan(highlight_start_time, highlight_end_time, ymin=(highlight_min-ymin)/(ymax-ymin), ymax=(highlight_max-ymin)/(ymax-ymin), alpha=0.30, facecolor=colors[Spin_labels.index(spin_lab)])
        else:
            print("Couldn't open", pickle_file)
        
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('Disk mass (M$_\odot$)', labelpad=-0.2)
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
#axs.flatten()[plot_it-1].set_ylim([0.6, 0.75])
plt.savefig('Disk_mass.pdf', bbox_inches='tight', pad_inches=0.02)
import pdb
pdb.set_trace()

#======================================================================================================================
#Plotting disk mass and specific angular momentum

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

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
            file.close()

            axs.flatten()[plot_it].plot(Time_array, Mean_rel_kep, label='$\Omega t_{ff}$='+spin_lab, ls=line_styles[Spin_labels.index(spin_lab)], alpha=0.75, color=colors[Spin_labels.index(spin_lab)], linewidth=1)
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
        else:
            print("Couldn't open", pickle_file)
        
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('h (m$^2$/s)', labelpad=-0.2)
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
#axs.flatten()[plot_it-1].set_ylim([0.6, 0.75])
plt.savefig('disk_h.pdf', bbox_inches='tight', pad_inches=0.02)

#======================================================================================================================
#Plotting disk mass radial velocity

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

        pickle_file = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/10au/gathered_profile.pkl'
        
        if os.path.exists(pickle_file):
            file = open(pickle_file, 'rb')
            R_sink, Time_array, Total_mass, Mean_mass, Mean_rel_kep, Min_rel_kep, Mean_rad_vel, Min_rad_vel, Mass_all, Separation = pickle.load(file)
            file.close()

            axs.flatten()[plot_it].plot(Time_array, Min_rad_vel, label='$\Omega t_{ff}$='+spin_lab, ls=line_styles[Spin_labels.index(spin_lab)], alpha=0.75, color=colors[Spin_labels.index(spin_lab)], linewidth=1)
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)', labelpad=-0.2)
        else:
            print("Couldn't open", pickle_file)
        
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('Radial velocity (m/s)', labelpad=-0.2)
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
#axs.flatten()[plot_it-1].set_ylim([0.6, 0.75])
plt.savefig('disk_rad_vel.pdf', bbox_inches='tight', pad_inches=0.02)
