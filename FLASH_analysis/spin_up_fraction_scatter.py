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
max_time = [[10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000]]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

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
args = parse_inputs()
Mach_labels = ['0.0', '0.1', '0.2']
Spin_labels = ['0.20', '0.25', '0.30', '0.35']

spin_val = [0.20, 0.25, 0.3, 0.35]
spin_up = [[], [], []]
spin_up_spec = [[], [], []]
spin_up_spec_peak = [[], [], []]


for mach_lab in Mach_labels:
    for spin_lab in Spin_labels:
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan

            if len(sink_data.keys()) == 1:
                spin_up_percentage = np.nan
                spin_up_percentage_spec = np.nan
                spin_up_percentage_spec_peak = np.nan
            else:
                prime_id = list(sink_data.keys())[0]
                second_id = list(sink_data.keys())[1]
                secondary_form_time = sink_data[second_id]['time'][0] - yt.YTQuantity(500, 'yr').in_units('s').value
                #secondary_form_time = yt.YTQuantity(2000, 'yr').in_units('s').value + sink_data[prime_id]['time'][0]
                secondary_form_ind = np.argmin(abs(sink_data[prime_id]['time'] - secondary_form_time))
                
                pre_sec_L = np.sqrt(sink_data[prime_id]['anglx'][secondary_form_ind-1]**2 + sink_data[prime_id]['angly'][secondary_form_ind-1]**2 + sink_data[prime_id]['anglz'][secondary_form_ind-1]**2)
                pre_sec_L_spec = pre_sec_L/sink_data[prime_id]['mass'][secondary_form_ind-1]
                
                post_sec_L = np.sqrt(sink_data[prime_id]['anglx'][secondary_form_ind:]**2 + sink_data[prime_id]['angly'][secondary_form_ind:]**2 + sink_data[prime_id]['anglz'][secondary_form_ind:]**2)
                post_sec_L_spec = post_sec_L/sink_data[prime_id]['mass'][secondary_form_ind:]
                
                post_sec_L_val = post_sec_L[-1]
                post_sec_L_spec_last = post_sec_L_spec[-1]
                post_sec_L_spec_peak = np.max(post_sec_L_spec)
                
                DL = post_sec_L_val - pre_sec_L
                spin_up_percentage = DL/pre_sec_L * 100
                
                DL_spec = post_sec_L_spec_last - pre_sec_L_spec
                spin_up_percentage_spec = DL_spec/pre_sec_L_spec * 100
                
                DL_spec_peak = post_sec_L_spec_peak - pre_sec_L_spec
                spin_up_percentage_spec_peak = DL_spec_peak/pre_sec_L_spec * 100
                print("Spin up: L=", spin_up_percentage, "%, L_spec=", spin_up_percentage_spec, "%")
        else:
            spin_up_percentage = np.nan
            spin_up_percentage_spec = np.nan
            spin_up_percentage_spec_peak = np.nan
        spin_up[int(mach_lab.split('.')[-1])].append(spin_up_percentage)
        spin_up_spec[int(mach_lab.split('.')[-1])].append(spin_up_percentage_spec)
        spin_up_spec_peak[int(mach_lab.split('.')[-1])].append(spin_up_percentage_spec_peak)
        
plt.clf()
fig = plt.figure(figsize=(single_col_width, 0.7*single_col_width))
for mach_lab in Mach_labels:
    plt.plot(spin_val, spin_up[int(mach_lab.split('.')[-1])], label='$\mathcal{M}=$'+mach_lab)
    plt.scatter(spin_val, spin_up[int(mach_lab.split('.')[-1])])
plt.xlabel('Initial cloud spin ($\Omega t_{ff}$)')
plt.ylabel('L spin up percentage (%)')
plt.legend(loc='best')
plt.savefig('spin_up_percentage_L.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig = plt.figure(figsize=(single_col_width, 0.7*single_col_width))
for mach_lab in Mach_labels:
    plt.plot(spin_val, spin_up_spec[int(mach_lab.split('.')[-1])], label='$\mathcal{M}=$'+mach_lab)
    plt.scatter(spin_val, spin_up_spec[int(mach_lab.split('.')[-1])])
plt.xlabel('Initial cloud spin ($\Omega t_{ff}$)')
plt.ylabel('h spin up percentage (%)')
plt.legend(loc='best')
plt.savefig('spin_up_percentage_h.pdf', bbox_inches='tight', pad_inches=0.02)

plt.clf()
fig = plt.figure(figsize=(single_col_width, 0.7*single_col_width))
for mach_lab in Mach_labels:
    plt.plot(spin_val, spin_up_spec_peak[int(mach_lab.split('.')[-1])], label='$\mathcal{M}=$'+mach_lab)
    plt.scatter(spin_val, spin_up_spec_peak[int(mach_lab.split('.')[-1])])
plt.xlabel('Initial cloud spin ($\Omega t_{ff}$)')
plt.ylabel('h spin up percentage (%)')
plt.legend(loc='best')
plt.savefig('spin_up_percentage_h_peak.pdf', bbox_inches='tight', pad_inches=0.02)
            
linestyle = ['-', '--', '-.']
plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=2, figsize=(single_col_width, single_col_width*1.5), sharex=True, sharey='row')
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)
axs.flatten()[0].grid()
for mach_lab in Mach_labels:
    axs.flatten()[0].plot(spin_val, spin_up[int(mach_lab.split('.')[-1])], label='$\mathcal{M}=$'+mach_lab, ls=linestyle[Mach_labels.index(mach_lab)], color='k')
    axs.flatten()[0].scatter(spin_val, spin_up[int(mach_lab.split('.')[-1])], color='k')
axs.flatten()[0].set_ylabel('$\Delta L$ (%)')
axs.flatten()[0].legend(loc='best')

axs.flatten()[0].tick_params(axis='x', direction='in', top=True)
axs.flatten()[0].tick_params(axis='y', direction='in', right=True)
axs.flatten()[0].minorticks_on()
axs.flatten()[0].tick_params(which='both', direction='in', axis='both', right=True, top=True)

axs.flatten()[1].grid()
for mach_lab in Mach_labels:
    axs.flatten()[1].plot(spin_val, spin_up_spec_peak[int(mach_lab.split('.')[-1])], label='$\mathcal{M}=$'+mach_lab, ls=linestyle[Mach_labels.index(mach_lab)], color='k')
    axs.flatten()[1].scatter(spin_val, spin_up_spec_peak[int(mach_lab.split('.')[-1])], color='k')
axs.flatten()[1].set_xlabel('Initial cloud spin ($\Omega t_{ff}$)')
axs.flatten()[1].set_ylabel('$\Delta h$ (%)')
axs.flatten()[1].tick_params(axis='x', direction='in', top=True)
axs.flatten()[1].tick_params(axis='y', direction='in', right=True)
axs.flatten()[1].minorticks_on()
axs.flatten()[1].tick_params(which='both', direction='in', axis='both', right=True, top=True)
plt.savefig('spin_up_percentage.pdf', bbox_inches='tight', pad_inches=0.02)
