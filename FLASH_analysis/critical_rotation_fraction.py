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

plt.clf()
fig, axs = plt.subplots(ncols=len(Mach_labels), nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0
xmax= 0
ymax = 0
plotted_legend = False
for spin_lab in Spin_labels:
    for mach_lab in Mach_labels:
        axs.flatten()[plot_it].grid()
        #single_pickle
        single_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Single_Mach_'+mach_lab+'_Lref_9.pkl'
        binary_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_lab+'_Binary_Mach_'+mach_lab+'_Lref_9.pkl'
        
        if os.path.exists(single_pickle):
            file = open(single_pickle, 'rb')
            sink_data = pickle.load(file)
            file.close()
            form_time = np.nan
            
            #for sink_id in sink_data.keys():
            #sink_id = list(sink_data.keys())[0]
            for sink_id in sink_data.keys():
                Fast_rotator_rate = (2*np.pi)/yt.YTQuantity(2, 'day').in_units('s')
                Mass = yt.YTArray(sink_data[sink_id]['mass'], 'g')
            
                Radius = yt.YTQuantity(2, 'rsun').in_units('cm')
                break_up_frequency =  np.sqrt(yt.units.gravitational_constant_cgs*Mass/(Radius**3))
                Momentum_of_inertia_sphere = 2/5 * Mass * Radius**2
                L_sphere_break_up = Momentum_of_inertia_sphere * break_up_frequency
                
                
                R_sink = yt.YTQuantity(4.9, 'au').in_units('cm')
                break_up_frequency_sink = np.sqrt(yt.units.gravitational_constant_cgs*Mass/(R_sink**3))
                Momentum_of_inertia_sphere_sink = 2/5 * Mass * R_sink**2
                L_sphere_break_up_sink = Momentum_of_inertia_sphere_sink * break_up_frequency_sink
                
                L_tot = np.sqrt(sink_data[sink_id]['anglx']**2 + sink_data[sink_id]['angly']**2 + sink_data[sink_id]['anglz']**2)
                L_tot = yt.YTArray(L_tot, 'g*cm**2/s')
            
                if np.isnan(form_time):
                    form_time = sink_data[sink_id]['time'][0]
                time = sink_data[sink_id]['time'] - form_time
                time = yt.YTArray(time, 's')
                if time[-1] > xmax:
                    xmax = time[-1]
                if np.max(L_tot) > ymax:
                    ymax = np.max(L_tot)
                axs.flatten()[plot_it].plot(time.in_units('yr'), L_tot, ls='-', label='Sink')
                axs.flatten()[plot_it].plot(time.in_units('yr'), L_sphere_break_up_sink, ls='--', label='Sink breakup')
                axs.flatten()[plot_it].plot(time.in_units('yr'), L_sphere_break_up, ls=':', label='Star breakup')
                if plotted_legend == False:
                    axs.flatten()[plot_it].legend()
                    plotted_legend = True
        else:
            print("Couldn't open", single_pickle)
            
        if plot_it == 0:
            axs.flatten()[plot_it].legend()
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')
            if spin_lab == '0.20':
                axs.flatten()[plot_it].set_title('Mach ='+mach_lab)
        if mach_lab == '0.2':
            if spin_lab == '0.20':
                axs.flatten()[plot_it].set_title('Mach ='+mach_lab)
        if spin_lab == '0.35':
            axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
        
        plot_it = plot_it + 1

axs.flatten()[plot_it-1].set_xlim(left=0)
axs.flatten()[plot_it-1].set_ylim(bottom=0)
plt.savefig('angular_frequency_comp.png')
