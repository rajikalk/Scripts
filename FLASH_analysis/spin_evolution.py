import numpy as np
import pickle
import yt
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib
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

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10

pickle_files = ['/hits/fast/set/kuruwira/Analysis/L_evolution/Spin_0.35_Binary_Mach_0.2_Lref_9_gathered_ang_mom.pkl', '/hits/fast/set/kuruwira/Analysis/L_evolution/Spin_0.35_Single_Mach_0.2_Lref_9_gathered_ang_mom.pkl']

plt.clf()
labels = ['Binary', 'Single']
pit = -1
for pickle_file in pickle_files:
    pit = pit + 1
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    plt.semilogy(Time_array - Time_array[0], L_primary, label=labels[pit] + ' Primary spin')
    plt.semilogy(Time_array - Time_array[0], L_secondary, label=labels[pit] + ' Secondary spin')
    
plt.xlabel('Time (yr)')
plt.ylabel('Angular momentum (g cm$^2$/s)')
plt.legend(loc='best')
plt.xlim(left=0)
plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
plt.tick_params(axis='x', direction='in')
plt.tick_params(axis='y', direction='in')
plt.ylim([1.e-6, 1])
plt.savefig('spin_comp.png', bbox_inches='tight')

plt.clf()
labels = ['Binary', 'Single']
pit = -1
for pickle_file in pickle_files:
    pit = pit + 1
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    L_tot = L_primary + np.nan_to_num(L_secondary) + L_orbit + L_in_gas
    
    plt.semilogy(Time_array - Time_array[0], L_primary/L_tot, label=labels[pit] + ' Primary spin')
    plt.semilogy(Time_array - Time_array[0], L_secondary/L_tot, label=labels[pit] + ' Secondary spin')
    
plt.xlabel('Time (yr)')
plt.ylabel('Angular momentum fraction (%)')
plt.legend(loc='best')
plt.xlim(left=0)
plt.tick_params(axis='both', which='major', labelsize=font_size, right=True)
plt.tick_params(axis='both', which='minor', labelsize=font_size, right=True)
plt.tick_params(axis='x', direction='in')
plt.tick_params(axis='y', direction='in')
plt.ylim([1.e-6, 1])
plt.savefig('spin_comp_frac.png', bbox_inches='tight')
    
    
