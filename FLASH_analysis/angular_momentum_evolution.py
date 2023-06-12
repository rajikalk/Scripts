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
Spin_labels = ['0.20', '0.25', '0.30', '0.35']

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
input_dir = sys.argv[1]
args = parse_inputs()
files = sorted(glob.glob(input_dir + '*plt_cnt*'))

plot_it = 0
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
    
        single_pickle = '/home/kuruwira/fast/Analysis/Angular_momentum_budget/Flash_2023/Spin_'+spin_lab+'/Binary/Mach_'+mach_lab+'/Lref_9/Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_9_gathered_ang_mom.pkl'
    
    file = open('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'gathered_ang_mom.pkl', 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    plt.clf()
    plt.semilogy(Time_array - Time_array[0], L_primary, label='Primary spin')
    plt.semilogy(Time_array - Time_array[0], L_secondary, label='Secondary spin')
    plt.semilogy(Time_array - Time_array[0], L_orbit, label='Orbital L')
    plt.semilogy(Time_array - Time_array[0], L_in_gas, label='L_in_gas')
    plt.xlabel('Time (yr)')
    plt.ylabel('Angular momentum (g cm$^2$/s)')
    plt.legend(loc='best')
    plt.xlim(left=0)
    plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
    plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.savefig('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'L_evolution.png', bbox_inches='tight')
    
    L_tot = L_primary + np.nan_to_num(L_secondary) + L_orbit + L_in_gas
    plt.clf()
    plt.semilogy(Time_array - Time_array[0], L_primary/L_tot, label='Primary spin')
    plt.semilogy(Time_array - Time_array[0], L_secondary/L_tot, label='Secondary spin')
    plt.semilogy(Time_array - Time_array[0], L_orbit/L_tot, label='Orbital L')
    plt.semilogy(Time_array - Time_array[0], L_in_gas/L_tot, label='L_in_gas')
    plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
    plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    plt.tick_params(axis='x', direction='in')
    plt.tick_params(axis='y', direction='in')
    plt.xlabel('Time (yr)')
    plt.ylabel('Angular momentum fraction (%)')
    plt.legend(loc='best')
    plt.xlim(left=0)
    plt.ylim(top=1)
    plt.savefig('_'.join(input_dir.split('Flash_2023/')[-1].split('/'))+'L_evolution_frac.png', bbox_inches='tight')
    
    
    print('saved figure')

sys.stdout.flush()
CW.Barrier()
