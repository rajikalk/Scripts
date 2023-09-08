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
#Get simulation files
#input_dir = sys.argv[1]
args = parse_inputs()

'''
plt.clf()
sink_ind = 0
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    sink_data = pickle.load(file)
    file.close()
    form_time = np.nan
    
    for sink_id in sink_data.keys():
        if np.isnan(form_time):
            form_time = sink_data[sink_id]['time'][0]
        #L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
        #L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
        L_tot = np.sqrt((sink_data[sink_id]['anglx']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['angly']/sink_data[sink_id]['mass'])**2 + (sink_data[sink_id]['anglz']/sink_data[sink_id]['mass'])**2)
        L_tot = yt.YTArray(L_tot, 'cm**2/s')
        time = sink_data[sink_id]['time'] - form_time
        time = yt.YTArray(time, 's')
        plt.plot(time.in_units('yr'), L_tot, label=sink_id + "_" + pickle_file.split('/')[-1])

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.ylim(bottom=0)
plt.legend()
plt.savefig('spin_comp_primary.png')
'''
Mach_labels = ['0.0', '0.1', '0.2']
Spin_labels = ['0.20', '0.25', '0.30', '0.35']
max_time = [[10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000]]
plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*1.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan
            
            #for sink_id in sink_data.keys():
            #sink_id = list(sink_data.keys())[0]
            for sink_id in sink_data.keys():
                if np.isnan(form_time):
                    form_time = sink_data[sink_id]['time'][0]
                mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
                time = sink_data[sink_id]['time'] - form_time
                time = yt.YTArray(time, 's')
                if max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)] == None:
                    plot_time = time.in_units('yr')
                    plot_mass = mass.in_units('msun')
                else:
                    end_time = max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)]
                    end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
                    plot_time = time.in_units('yr')[:end_ind+1]
                    plot_mass = mass.in_units('msun')[:end_ind+1]
                axs.flatten()[plot_it].plot(plot_time, plot_mass)
        else:
            print("Couldn't open", single_pickle)

        if spin_lab == '0.20':
            axs.flatten()[plot_it].set_title('Mach ='+mach_lab)
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: Sink Mass (M$_\odot$)')
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
axs.flatten()[plot_it-1].set_ylim(bottom=0)
plt.savefig('Mass_evol.pdf', bbox_inches='tight')


plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = -1
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    plot_it = plot_it + 1
    for mach_lab in Mach_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan
            
            for sink_id in sink_data.keys():
            #sink_id = list(sink_data.keys())[0]
                if np.isnan(form_time):
                    form_time = sink_data[list(sink_data.keys())[0]]['time'][0]
                mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
                time = sink_data[sink_id]['time'] - form_time
                time = yt.YTArray(time, 's')
                if max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)] == None:
                    plot_time = time.in_units('yr')
                    plot_mass = mass.in_units('msun')
                else:
                    end_time = max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)]
                    end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
                    plot_time = time.in_units('yr')[:end_ind+1]
                    plot_mass = mass.in_units('msun')[:end_ind+1]
                if sink_id == list(sink_data.keys())[0]:
                    axs.flatten()[plot_it].plot(plot_time, plot_mass, label='$\mathcal{M}=$'+mach_lab, color=colors[Mach_labels.index(mach_lab)], alpha=0.75)
                else:
                    axs.flatten()[plot_it].plot(plot_time, plot_mass, color=colors[Mach_labels.index(mach_lab)], alpha=0.75)
        else:
            print("Couldn't open", single_pickle)

        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('Sink Mass (M$_\odot$)')
            if spin_lab != '0.20':
                yticklabels = axs.flatten()[plot_it].get_yticklabels()
                plt.setp(yticklabels[-1], visible=False)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
            if mach_lab != '0.2' and spin_lab == Spin_labels[-1]:
                xticklabels = axs.flatten()[plot_it].get_xticklabels()
                plt.setp(xticklabels[-1], visible=False)
        

axs.flatten()[0].legend(loc='best')
axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim(bottom=0)
plt.savefig('Mass_evol_comp.pdf', bbox_inches='tight')

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
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan
            
            #for sink_id in sink_data.keys():
            #sink_id = list(sink_data.keys())[0]
            if np.isnan(form_time):
                form_time = sink_data[list(sink_data.keys())[0]]['time'][0]
            mass = yt.YTArray(sink_data[list(sink_data.keys())[0]]['mass'], 'g')
            time = sink_data[list(sink_data.keys())[0]]['time'] - form_time
            time = yt.YTArray(time, 's')
            if max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)] == None:
                plot_time = time.in_units('yr')
                plot_mass = mass.in_units('msun')
            else:
                end_time = max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)]
                end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
                plot_time = time.in_units('yr')[:end_ind+1]
                plot_mass = mass.in_units('msun')[:end_ind+1]
            axs.flatten()[plot_it].plot(plot_time, plot_mass, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], alpha=0.75)
        else:
            print("Couldn't open", single_pickle)

        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('Sink Mass (M$_\odot$)')
        else:
            yticklabels = axs.flatten()[plot_it].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        if mach_lab != '0.2' and spin_lab == Spin_labels[-1]:
            xticklabels = axs.flatten()[plot_it].get_xticklabels()
            plt.setp(xticklabels[-1], visible=False)

axs.flatten()[0].legend(loc='best')
axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim(bottom=0)
plt.savefig('Mass_evol_spin_comp.pdf', bbox_inches='tight')

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
    axs.flatten()[plot_it].set_title('Mach='+mach_lab)
    for spin_lab in Spin_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_'+args.refinment_level+'.pkl'
        #binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data, line_counter = pickle.load(file)
            file.close()
            form_time = np.nan
            
            Total_accrete_mass = []
            Total_time_arr = []
            for sink_id in sink_data.keys():
                #sink_id = list(sink_data.keys())[0]
                if np.isnan(form_time):
                    form_time = sink_data[sink_id]['time'][0]
                mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
                time = sink_data[sink_id]['time'] - form_time
                time = yt.YTArray(time, 's')
                if max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)] == None:
                    plot_time = time.in_units('yr')
                    plot_mass = mass.in_units('msun')
                else:
                    end_time = max_time[Spin_labels.index(spin_lab)][Mach_labels.index(mach_lab)]
                    end_ind = np.argmin(abs(time.in_units('yr').value - end_time))
                    plot_time = time.in_units('yr')[:end_ind+1]
                    plot_mass = mass. in_units('msun')[:end_ind+1]
                if sink_id == list(sink_data.keys())[0]:
                    axs.flatten()[plot_it].plot(plot_time, plot_mass, label='$\Omega t_{ff}$='+spin_lab, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
                else:
                    axs.flatten()[plot_it].plot(plot_time, plot_mass, linestyle=line_styles[Spin_labels.index(spin_lab)], color=colors[Spin_labels.index(spin_lab)], alpha=0.75)
                if len(Total_time_arr) == 0:
                    Total_accrete_mass = mass.in_units('msun').value
                    Total_time_arr = time.in_units('yr')
                else:
                    add_mass = np.concatenate((np.zeros(len(Total_accrete_mass) - len(mass.in_units('msun'))), mass.in_units('msun').value))
                    Total_accrete_mass = Total_accrete_mass + add_mass
            axs.flatten()[plot_it].plot(Total_time_arr, Total_accrete_mass, linewidth=3, alpha=0.25, linestyle=line_styles[Spin_labels.index(spin_lab)])
                
        else:
            print("Couldn't open", single_pickle)

        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('Sink Mass (M$_\odot$)')
        else:
            yticklabels = axs.flatten()[plot_it].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        if mach_lab != '0.2' and spin_lab == '0.35':
            xticklabels = axs.flatten()[plot_it].get_xticklabels()
            plt.setp(xticklabels[-1], visible=False)
        if mach_lab == '0.0' and spin_lab == '0.35':
            xticklabels = axs.flatten()[plot_it].get_xticklabels()
            #plt.setp(xticklabels[-1], visible=False)

axs.flatten()[0].legend(loc='upper left')
axs.flatten()[plot_it-1].set_xlim([0, 10000])
axs.flatten()[plot_it-1].set_ylim(bottom=0)
plt.savefig('Mass_evol_spin_comp_all_sinks.pdf', bbox_inches='tight')
