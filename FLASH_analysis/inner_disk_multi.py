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

#---------------------------------------------------
#Define arguments
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-update", "--update_pickles", help="do you want to read the Flash output and update the pickles?", type=str, default='True')
    parser.add_argument("-lref", "--refinment_level", type=str, default='9')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#---------------------------------------------------
#Get simulation files
#input_dir = sys.argv[1]
args = parse_inputs()

Mach_labels = ['0.0', '0.1', '0.2']
Spin_labels = ['0.20', '0.25', '0.30', '0.35']
linestyles = ['-', '--', ':']

#=========================================================================

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
ymin = np.inf
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            if np.nanmax(Separation) > max_sep:
                max_sep = np.nanmax(Separation)

for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            
            #if spin_lab == '0.20' and mach_lab == '0.1':
            #	import pdb
            #	pdb.set_trace()
            
            if np.min(Total_L) < ymin:
            	ymin = np.min(Total_L)
            if spin_lab != '0.20' and mach_lab != '0.1':
            	if np.max(Total_L) > ymax:
            		ymax = np.max(Total_L)
            
            ax2 = axs.flatten()[plot_it].twinx()
            #axs.flatten()[plot_it].semilogy(Time_array, Total_L, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            axs.flatten()[plot_it].plot(Time_array, Total_L, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)], alpha=0.8)
            ax2.plot(Time_array, Separation, color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
            ax2.axhline(y=20, color='k', linewidth=0.5)
            if mach_lab == '0.2':
                ax2.set_ylabel('Separation (AU)')
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].legend(loc='best')
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('L ($g\,cm^2/s$)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels[-1], visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
    plot_it = plot_it + 1
    axs.flatten()[plot_it-1].set_xlim([0, 10000])
    plt.savefig('Total_Inner_disk_L_mach_comp.pdf', bbox_inches='tight')
    
#axs.flatten()[plot_it-1].set_ylim(top=1.e58)
axs.flatten()[plot_it-1].set_ylim([ymin, ymax])
axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Total_Inner_disk_L_mach_comp.pdf', bbox_inches='tight')
print('saved figure Total_Inner_disk_L_mach_comp.pdf')

#================================================================
#Specific L comparison

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
ymin = np.inf
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            if np.nanmax(Separation) > max_sep:
                max_sep = np.nanmax(Separation)

for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            
            #if spin_lab == '0.20' and mach_lab == '0.1':
            #	import pdb
            #	pdb.set_trace()
            
            if np.min(Total_L_spec) < ymin:
            	ymin = np.min(Total_L_spec)
            if spin_lab != '0.20' and mach_lab != '0.1':
            	if np.max(Total_L_spec) > ymax:
            		ymax = np.max(Total_L_spec)
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(Time_array, Total_L_spec, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)], alpha=0.8)
            #axs.flatten()[plot_it].semilogy(Time_array, Total_L_spec, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.plot(Time_array, Separation, color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
            ax2.axhline(y=20, color='k', linewidth=0.5)
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].legend(loc='upper right')
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('h ($cm^2/s$)')
            ax2.set_ylabel('Separation (AU)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels[-1], visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
    plot_it = plot_it + 1
    axs.flatten()[plot_it-1].set_xlim([0, 10000])
    plt.savefig('Total_Inner_disk_L_mach_comp_spec.pdf', bbox_inches='tight')
    

#axs.flatten()[plot_it-1].set_ylim(top=1.e34)
axs.flatten()[plot_it-1].set_ylim([ymin, ymax])
axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Total_Inner_disk_L_mach_comp_spec.pdf', bbox_inches='tight')
print('saved figure Total_Inner_disk_L_mach_comp_spec.pdf')

#=========================================================================

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
ymin = np.inf
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            if np.nanmax(Separation) > max_sep:
                max_sep = np.nanmax(Separation)

for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            
            if np.min(Mean_L) < ymin:
            	ymin = np.min(Mean_L)
            if spin_lab != '0.20' and mach_lab != '0.1':
            	if np.max(Mean_L) > ymax:
            		ymax = np.max(Mean_L)
            
            ax2 = axs.flatten()[plot_it].twinx()
            #axs.flatten()[plot_it].semilogy(Time_array, Total_L, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            axs.flatten()[plot_it].plot(Time_array, Mean_L, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)], alpha=0.8)
            ax2.plot(Time_array, Separation, color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
            ax2.axhline(y=20, color='k', linewidth=0.5)
            if mach_lab == '0.2':
                ax2.set_ylabel('Separation (AU)')
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].legend(loc='best')
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('L ($g\,cm^2/s$)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels[-1], visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
    plot_it = plot_it + 1
    axs.flatten()[plot_it-1].set_xlim([0, 10000])
    plt.savefig('Mean_Inner_disk_L_mach_comp.pdf', bbox_inches='tight')
    
#axs.flatten()[plot_it-1].set_ylim(top=1.e52)
axs.flatten()[plot_it-1].set_ylim([ymin, ymax])
axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Mean_Inner_disk_L_mach_comp.pdf', bbox_inches='tight')
print('saved figure Mean_Inner_disk_L_mach_comp.pdf')

#================================================================
#Specific L comparison

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
ymin = np.inf
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            if np.nanmax(Separation) > max_sep:
                max_sep = np.nanmax(Separation)

for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Total_inner_disk_values/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        #inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/Total_L/gathered_profile.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Separation = pickle.load(file)
            file.close()
            
            if np.min(Mean_L_spec) < ymin:
            	ymin = np.min(Mean_L_spec)
            if spin_lab != '0.20' and mach_lab != '0.1':
            	if np.max(Mean_L_spec) > ymax:
            		ymax = np.max(Mean_L_spec)
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(Time_array, Mean_L_spec, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)], alpha=0.8)
            #axs.flatten()[plot_it].semilogy(Time_array, Total_L_spec, label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.plot(Time_array, Separation, color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
            ax2.axhline(y=20, color='k', linewidth=0.5)
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].legend(loc='upper right')
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('h ($cm^2/s$)')
            ax2.set_ylabel('Separation (AU)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels[-1], visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
    plot_it = plot_it + 1
    axs.flatten()[plot_it-1].set_xlim([0, 10000])
    plt.savefig('Mean_Inner_disk_L_mach_comp_spec.pdf', bbox_inches='tight')
    

#axs.flatten()[plot_it-1].set_ylim(top=1.e24)
axs.flatten()[plot_it-1].set_ylim([ymin, ymax])
axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Mean_Inner_disk_L_mach_comp_spec.pdf', bbox_inches='tight')
print('saved figure Mean_Inner_disk_L_mach_comp_spec.pdf')

