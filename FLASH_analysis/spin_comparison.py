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
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#---------------------------------------------------
#Get simulation files
#input_dir = sys.argv[1]
args = parse_inputs()
pickle_files = sorted(sys.argv[1:])

plt.clf()
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    plt.plot(Time_array - Time_array[0], L_primary, label=Spin_label)

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.legend()
plt.savefig('spin_comp_primary.png')


plt.clf()
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    L_tot = L_primary + np.nan_to_num(L_secondary) + L_orbit + L_in_gas
    
    plt.plot(Time_array - Time_array[0], L_primary/L_tot, label=Spin_label)

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.legend()
plt.savefig('spin_comp_primary_frac.png')

plt.clf()
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    plt.plot(Time_array - Time_array[0], L_secondary, label=Spin_label)

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.legend()
plt.savefig('spin_comp_secondary.png')


plt.clf()
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    L_tot = L_primary + np.nan_to_num(L_secondary) + L_orbit + L_in_gas
    
    plt.plot(Time_array - Time_array[0], L_secondary/L_tot, label=Spin_label)

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.legend()
plt.savefig('spin_comp_secondary_frac.png')

plt.clf()
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    plt.plot(Time_array - Time_array[0], L_orbit, label=Spin_label)

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.legend()
plt.savefig('spin_comp_orbit.png')


plt.clf()
for pickle_file in pickle_files:
    Spin_label = " ".join(pickle_file.split('/')[-1].split('_')[:2])
    file = open(pickle_file, 'rb')
    Time_array, L_primary, L_secondary, L_orbit, L_in_gas = pickle.load(file)
    file.close()
    
    L_tot = L_primary + np.nan_to_num(L_secondary) + L_orbit + L_in_gas
    
    plt.plot(Time_array - Time_array[0], L_secondary/L_tot, label=Spin_label)

plt.xlabel('Time since formation')
plt.ylabel('Angular momentum')
plt.xlim(left=0)
plt.legend()
plt.savefig('spin_comp_orbit_frac.png')

