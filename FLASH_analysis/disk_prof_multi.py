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
linestyles = ['-', '--', '-.', ':']

#=========================================================================

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            if np.nanmax(mean_inner_all.T[1]) > max_sep:
                max_sep = np.nanmax(mean_inner_all.T[1])

for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(mean_inner_all.T[0], mean_inner_all.T[2], label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.plot(mean_inner_all.T[0], mean_inner_all.T[1], color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
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
    plt.savefig('Inner_disk_L_mach_comp.pdf', bbox_inches='tight')
    

axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Inner_disk_L_mach_comp.pdf', bbox_inches='tight')
print('saved figure Inner_disk_L_mach_comp.pdf')

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            if np.nanmax(mean_inner_all.T[1]) > max_sep:
                max_sep = np.nanmax(mean_inner_all.T[1])

for mach_lab in Mach_labels:
    for spin_lab in Spin_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(mean_inner_all.T[0], mean_inner_all.T[2], label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.plot(mean_inner_all.T[0], mean_inner_all.T[1], color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
        else:
            print("Couldn't open", inner_pickle)
            
        if mach_lab == '0.0':
            axs.flatten()[plot_it].legend(loc='best')
            axs.flatten()[plot_it].set_ylabel('L ($g\,cm^2/s$)')
        if mach_lab == '0.20':
            ax2.set_ylabel('Separation (AU)')
        if spin_lab != '0.0':
            yticklabels = axs.flatten()[plot_it].get_yticklabels()
            plt.setp(yticklabels, visible=False)
                
    axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
    plot_it = plot_it + 1
    axs.flatten()[plot_it-1].set_xlim([0, 10000])
    plt.savefig('Inner_disk_L_spin_comp.pdf', bbox_inches='tight')
    

axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Inner_disk_L_spin_comp.pdf', bbox_inches='tight')
print('saved figure Inner_disk_L_spin_comp.pdf')


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
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            if np.nanmax(mean_inner_all.T[1]) > max_sep:
                max_sep = np.nanmax(mean_inner_all.T[1])

for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(mean_inner_all.T[0], mean_inner_all.T[3], label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.plot(mean_inner_all.T[0], mean_inner_all.T[1], color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].legend(loc='best')
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
    plt.savefig('Inner_disk_L_mach_comp_spec.pdf', bbox_inches='tight')
    

axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Inner_disk_L_mach_comp_spec.pdf', bbox_inches='tight')
print('saved figure Inner_disk_L_mach_comp_spec.pdf')

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=1, figsize=(two_col_width, single_col_width), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
max_sep = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            if np.nanmax(mean_inner_all.T[1]) > max_sep:
                max_sep = np.nanmax(mean_inner_all.T[1])

for mach_lab in Mach_labels:
    for spin_lab in Spin_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(mean_inner_all.T[0], mean_inner_all.T[3], label='$\mathcal{M}$='+mach_lab, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.plot(mean_inner_all.T[0], mean_inner_all.T[1], color='k', alpha=0.20, ls=linestyles[Mach_labels.index(mach_lab)])
            ax2.set_ylim([0, max_sep])
        else:
            print("Couldn't open", inner_pickle)
            
        if mach_lab == '0.0':
            axs.flatten()[plot_it].legend(loc='best')
            axs.flatten()[plot_it].set_ylabel('h ($g\,cm^2/s$)')
        if mach_lab == '0.20':
            ax2.set_ylabel('Separation (AU)')
        if spin_lab != '0.0':
            yticklabels = axs.flatten()[plot_it].get_yticklabels()
            plt.setp(yticklabels, visible=False)
                
    axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
    plot_it = plot_it + 1
    axs.flatten()[plot_it-1].set_xlim([0, 10000])
    plt.savefig('Inner_disk_L_spin_comp_spec.pdf', bbox_inches='tight')
    

axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Inner_disk_L_spin_comp_spec.pdf', bbox_inches='tight')
print('saved figure Inner_disk_L_spin_comp_spec.pdf')

#================================================================


plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/mean_inner_L.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            mean_inner_all = pickle.load(file)
            file.close()
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(mean_inner_all.T[0], mean_inner_all.T[3], label='<h>')
            ax2.plot(mean_inner_all.T[0], mean_inner_all.T[1], color='k', alpha=0.20, label='Separation')
            ax2.set_ylim([0, max_sep])
            if mach_lab == '0.2':
                ax2.set_ylabel('Separation (AU)')
            if mach_lab != '0.2':
                yticklabels = ax2.get_yticklabels()
                plt.setp(yticklabels, visible=False)
            
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].set_title('Mach ='+mach_lab)
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('h ($cm^2/s$)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels[-1], visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
            if mach_lab != '0.2':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
        
        plot_it = plot_it + 1
        axs.flatten()[plot_it-1].set_xlim([0, 10000])
        plt.savefig('Inner_disk_L_spec.pdf', bbox_inches='tight')

axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Inner_disk_L_spec.pdf', bbox_inches='tight')
print('saved figure Inner_disk_L_spec.pdf')

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        inner_pickle = '/home/kuruwira/fast/Analysis/Disk_L_profiles/Spin_'+spin_lab+'/Mach_'+mach_lab+'/gathered_profile.pkl'
        axs.flatten()[plot_it].grid()
        #single_pickle

        if os.path.exists(inner_pickle):
            file = open(inner_pickle, 'rb')
            Time_array, Radius_array, All_profiles_array, Total_inner_disk = pickle.load(file)
            file.close()
            
            Radius_array = np.array(Radius_array).astype(float)
            Radius_array[np.where(Radius_array==100)] = np.nan
            
            ax2 = axs.flatten()[plot_it].twinx()
            axs.flatten()[plot_it].plot(Time_array, Total_inner_disk, label='L_{total}')
            ax2.plot(Time_array, Radius_array, color='k', alpha=0.20, label='Separation')
            ax2.set_ylim([0, max_sep])
            if mach_lab == '0.2':
                ax2.set_ylabel('Separation (AU)')
            if mach_lab != '0.2':
                yticklabels = ax2.get_yticklabels()
                plt.setp(yticklabels, visible=False)
            
        else:
            print("Couldn't open", inner_pickle)
            
        if spin_lab == '0.20':
            axs.flatten()[plot_it].set_title('Mach ='+mach_lab)
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('L$_{total}$ ($cm^2/s$)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels, visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
            if mach_lab != '0.2':
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
        
        plot_it = plot_it + 1
        axs.flatten()[plot_it-1].set_xlim([0, 10000])
        plt.savefig('Inner_total_L.pdf', bbox_inches='tight')

axs.flatten()[plot_it-1].set_xlim([0, 10000])
plt.savefig('Inner_total_L.pdf', bbox_inches='tight')
print('saved figure Inner_total_L.pdf')




