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

Mach_labels = ['0.0', '0.1', '0.2']
mach_ls = ['-', '--', '-.']
Spin_labels = ['0.20', '0.25', '0.30', '0.35']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

#---------------------------------------------------
#Define arguments
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-update", "--update_pickles", help="do you want to read the Flash output and update the pickles?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#---------------------------------------------------
#Get simulation files
args = parse_inputs()
#files = sorted(glob.glob(input_dir + '*plt_cnt*'))

plt.clf()
fig, axs = plt.subplots(ncols=3, nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = -1
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        plot_it = plot_it + 1
        axs.flatten()[plot_it].grid()

        single_pickle = '/home/kuruwira/fast/Analysis/Angular_momentum_budget/Flash_2023/Spin_'+spin_lab+'/Single/Mach_'+mach_lab+'/Lref_9/Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_9_gathered_ang_mom.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            Time_array, L_sink, L_orbit, L_in_gas = pickle.load(file)
            file.close()
            
            for time_it in range(len(L_orbit)):
                try:
                    if len(L_orbit[time_it]) == 3:
                        L_orbit[time_it] = yt.YTQuantity(np.nan, 'cm**2*g/s')
                except:
                    pass
            
            L_sink_tot = np.zeros(np.shape(L_in_gas))
            for sink_id in L_sink.keys():
                L_sink_tot = L_sink_tot + (np.append(np.zeros(len(L_sink_tot)-len(np.array(L_sink[sink_id]).T[1])), np.array(L_sink[sink_id]).T[1]))
            L_tot = L_sink_tot + L_orbit + L_in_gas
            
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orbit/L_tot, linestyle = '--', color='grey', linewidth=3, label='Orbit')
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas/L_tot, linestyle = ':', color='grey', linewidth=3, label='Gas')
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_sink_tot/L_tot, linestyle = '-', color='grey', linewidth=3, label = 'Sink total')
            cit = -1
            for sink_id in L_sink.keys():
                cit = cit + 1
                axs.flatten()[plot_it].semilogy(np.array(L_sink[sink_id]).T[0]-Time_array[0], np.array(L_sink[sink_id]).T[1]/L_tot[-len(np.array(L_sink[sink_id]).T[1]):], linestyle = '-', color=colors[cit], linewidth=1.5)
        else:
            print("Couldn't open", single_pickle)
            
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')
        
        if spin_lab == '0.20':
            axs.flatten()[plot_it].set_title('Mach='+mach_lab)
        if spin_lab == '0.20' and mach_lab == '0.0':
            axs.flatten()[plot_it].legend()
    if spin_lab == '0.35':
        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='x', direction='in')
    axs.flatten()[plot_it].tick_params(axis='y', direction='in')

axs.flatten()[plot_it].set_xlim([0, 10000])
axs.flatten()[plot_it].set_ylim([5.e-5,1.1])
plt.savefig('L_evolution_spin_vs_mach_frac.pdf', bbox_inches='tight')
