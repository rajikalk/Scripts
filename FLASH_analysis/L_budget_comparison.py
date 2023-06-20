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
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
col_title = ['Single', 'Binary']

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

'''
plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = -1
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for col_tit in col_title:
        plot_it = plot_it + 1
        axs.flatten()[plot_it].grid()
        for mach_it in range(len(Mach_labels)):
            if np.remainder(plot_it, 2) == 0:
        
                single_pickle = '/home/kuruwira/fast/Analysis/Angular_momentum_budget/Flash_2023/Spin_'+spin_lab+'/Single/Mach_'+Mach_labels[mach_it]+'/Lref_9/Spin_'+spin_lab+'_Single_Mach_'+Mach_labels[mach_it]+'_Lref_9_gathered_ang_mom.pkl'
                
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
                    
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orbit, label='Orbit', linestyle = mach_ls[mach_it], color=colors[0])
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas, label='Gas', linestyle = mach_ls[mach_it], color=colors[1])
                    #axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_primary, label='Single', linestyle = mach_ls[mach_it], color=colors[2])
                    axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')
                else:
                    print("Couldn't open", single_pickle)
                
            else:
                binary_pickle = '/home/kuruwira/fast/Analysis/Angular_momentum_budget/Flash_2023/Spin_'+spin_lab+'/Binary/Mach_'+Mach_labels[mach_it]+'/Lref_9/Spin_'+spin_lab+'_Binary_Mach_'+Mach_labels[mach_it]+'_Lref_9_gathered_ang_mom.pkl'

                if os.path.exists(binary_pickle):
                    file = open(binary_pickle, 'rb')
                    Time_array, L_sink, L_orbit, L_in_gas = pickle.load(file)
                    file.close()
                    
                    for time_it in range(len(L_orbit)):
                        try:
                            if len(L_orbit[time_it]) == 3:
                                L_orbit[time_it] = yt.YTQuantity(np.nan, 'cm**2*g/s')
                        except:
                            pass
                    
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orbit, label='Orbit', linestyle = mach_ls[mach_it], color=colors[0])
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas, label='Gas', linestyle = mach_ls[mach_it], color=colors[1])
                    #axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_primary, label='Primary', linestyle = mach_ls[mach_it], color=colors[2])
                    #axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_secondary, label='Secondary', linestyle = mach_ls[mach_it], color=colors[3])
                else:
                    print("Couldn't open", binary_pickle)
        
    if spin_lab == '0.20':
        axs.flatten()[plot_it].legend()
        if plot_it == 0:
            axs.flatten()[plot_it].set_title('Single')
        else:
            axs.flatten()[plot_it].set_title('Binary')
    if spin_lab == '0.35':
        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='x', direction='in')
    axs.flatten()[plot_it].tick_params(axis='y', direction='in')

axs.flatten()[plot_it].set_xlim(left=0)
axs.flatten()[plot_it].set_ylim([5.e48, 5.e54])
plt.savefig('L_evolution.png', bbox_inches='tight')

plt.clf()
fig, axs = plt.subplots(ncols=2, nrows=len(Spin_labels), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
iter_range = range(0, len(Spin_labels))
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = -1
xmax= 0
ymax = 0
for spin_lab in Spin_labels:
    for col_tit in col_title:
        plot_it = plot_it + 1
        for mach_it in range(len(Mach_labels)):
            if np.remainder(plot_it, 2) == 0:
                single_pickle = '/home/kuruwira/fast/Analysis/Angular_momentum_budget/Flash_2023/Spin_'+spin_lab+'/Single/Mach_'+Mach_labels[mach_it]+'/Lref_9/Spin_'+spin_lab+'_Single_Mach_'+Mach_labels[mach_it]+'_Lref_9_gathered_ang_mom.pkl'
                
                if os.path.exists(single_pickle):
                    file = open(single_pickle, 'rb')
                    Time_array, L_sink, L_orbit, L_in_gas = pickle.load(file)
                    file.close()
                    
                    L_orb_fixed = []
                    for time_it in range(len(L_orbit)):
                        try:
                            if len(L_orbit[time_it]) == 3:
                                L_orb_fixed.append(np.nan)
                        except:
                            try:
                                L_orb_fixed.append(L_orbit[time_it].value)
                            except:
                                L_orb_fixed.append(L_orbit[time_it])
                    
                    L_orb_fixed = yt.YTArray(L_orb_fixed, 'g*cm**2/s')

                    L_tot = np.nan_to_num(L_primary) + np.nan_to_num(L_secondary) + np.nan_to_num(L_orb_fixed).value + L_in_gas
                    #L_tot = yt.YTArray(np.nan_to_num(L_primary) + np.nan_to_num(L_secondary) + L_in_gas, 'g*cm**2/s')
                    #L_tot = L_tot + L_orbit
                    
                    #if np.min((L_primary/L_tot)[-1*int(len(L_primary)/2):])<1.e-4:
                    #    import pdb
                    #    pdb.set_trace()
                    
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orb_fixed/L_tot, label='Orbit', linestyle = mach_ls[mach_it], color=colors[0])
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas/L_tot, label='Gas', linestyle = mach_ls[mach_it], color=colors[1])
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_primary/L_tot, label='Single', linestyle = mach_ls[mach_it], color=colors[2])
                    axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')
                else:
                    print("Couldn't open", single_pickle)
                
                
            else:
                binary_pickle = '/home/kuruwira/fast/Analysis/Angular_momentum_budget/Flash_2023/Spin_'+spin_lab+'/Binary/Mach_'+Mach_labels[mach_it]+'/Lref_9/Spin_'+spin_lab+'_Binary_Mach_'+Mach_labels[mach_it]+'_Lref_9_gathered_ang_mom.pkl'

                if os.path.exists(binary_pickle):
                    file = open(binary_pickle, 'rb')
                    Time_array, L_sink, L_orbit, L_in_gas = pickle.load(file)
                    file.close()
                    
                    L_orb_fixed = []
                    for time_it in range(len(L_orbit)):
                        try:
                            if len(L_orbit[time_it]) == 3:
                                L_orb_fixed.append(np.nan)
                        except:
                            try:
                                L_orb_fixed.append(L_orbit[time_it].value)
                            except:
                                L_orb_fixed.append(L_orbit[time_it])
                    
                    L_orb_fixed = yt.YTArray(L_orb_fixed, 'g*cm**2/s')
                    
                    L_tot = np.nan_to_num(L_primary) + np.nan_to_num(L_secondary) + np.nan_to_num(L_orb_fixed).value + L_in_gas
                    
                    #if np.min((L_primary/L_tot)[-1*int(len(L_primary)/2):])<1.e-4:
                    #    import pdb
                    #    pdb.set_trace()
                    
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orb_fixed/L_tot, label='Orbit', linestyle = mach_ls[mach_it], color=colors[0])
                    axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas/L_tot, label='Gas', linestyle = mach_ls[mach_it], color=colors[1])
                    #axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_primary/L_tot, label='Primary', linestyle = mach_ls[mach_it], color=colors[2])
                    #axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_secondary/L_tot, label='Secondary', linestyle = mach_ls[mach_it], color=colors[3])
                else:
                    print("Couldn't open", binary_pickle)
    
        axs.flatten()[plot_it].grid()
        
    if spin_lab == '0.20':
        axs.flatten()[plot_it].legend()
        if plot_it == 0:
            axs.flatten()[plot_it].set_title('Single')
        else:
            axs.flatten()[plot_it].set_title('Binary')
    if spin_lab == '0.35':
        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='x', direction='in')
    axs.flatten()[plot_it].tick_params(axis='y', direction='in')
    
axs.flatten()[plot_it].set_ylim([5.e-5,1])
axs.flatten()[plot_it].set_xlim(left=0)
plt.savefig('L_evolution_frac.png', bbox_inches='tight')

sys.stdout.flush()
CW.Barrier()
'''
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
            
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orbit, linestyle = '--', color=colors[0], linewidth=2, label='Orbit')
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas, linestyle = ':', color=colors[1], linewidth=2, label='Gas')
            L_sink_tot = np.zeros(np.shape(L_in_gas))
            for sink_id in L_sink.keys():
                L_sink_tot = L_sink_tot + (np.append(np.zeros(len(L_sink_tot)-len(np.array(L_sink[sink_id]).T[1])), np.array(L_sink[sink_id]).T[1]))
                axs.flatten()[plot_it].semilogy(np.array(L_sink[sink_id]).T[0]-Time_array[0], np.array(L_sink[sink_id]).T[1], linestyle = '-', color=colors[3], linewidth=1)
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_sink_tot, linestyle = '-', color=colors[2], linewidth=2, label='Sink total')
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

axs.flatten()[plot_it].set_xlim(left=0)
axs.flatten()[plot_it].set_ylim([5.e48, 5.e54])
plt.savefig('L_evolution_spin_vs_mach.png', bbox_inches='tight')

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
            
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_orbit/L_tot, linestyle = '--', color=colors[0], linewidth=2, label='Orbit')
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_in_gas/L_tot, linestyle = ':', color=colors[1], linewidth=2, label='Gas')
            for sink_id in L_sink.keys():
                axs.flatten()[plot_it].semilogy(np.array(L_sink[sink_id]).T[0]-Time_array[0], np.array(L_sink[sink_id]).T[1]/L_tot[-len(np.array(L_sink[sink_id]).T[1]):], linestyle = '-', color=colors[3], linewidth=1)
            axs.flatten()[plot_it].semilogy(Time_array - Time_array[0], L_sink_tot/L_tot, linestyle = '-', color=colors[2], linewidth=2, label = 'Sink total')
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

axs.flatten()[plot_it].set_xlim(left=0)
axs.flatten()[plot_it].set_ylim([5.e-5,1])
plt.savefig('L_evolution_spin_vs_mach_frac.png', bbox_inches='tight')

plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2.5), sharex=True, sharey=True)
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
                
            for sink_id in L_sink.keys():
                axs.flatten()[plot_it].plot(np.array(L_sink[sink_id]).T[0]-Time_array[0], np.array(L_sink[sink_id]).T[1]/L_tot[-len(np.array(L_sink[sink_id]).T[1]):], linestyle = mach_ls[Mach_labels.index(mach_lab)], color=colors[Mach_labels.index(mach_lab)], linewidth=1)
            axs.flatten()[plot_it].plot(Time_array - Time_array[0], L_sink_tot/L_tot, linestyle = mach_ls[Mach_labels.index(mach_lab)], color=colors[Mach_labels.index(mach_lab)], linewidth=3, label='Mach '+mach_lab)
        else:
            print("Couldn't open", single_pickle)
            
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')

    if spin_lab == '0.20':
        axs.flatten()[plot_it].legend()
    if spin_lab == '0.35':
        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='x', direction='in')
    axs.flatten()[plot_it].tick_params(axis='y', direction='in')

axs.flatten()[plot_it].set_xlim(left=0)
axs.flatten()[plot_it].set_ylim(bottom=0)
plt.savefig('L_sink_spin_vs_mach_frac.png', bbox_inches='tight')


plt.clf()
fig, axs = plt.subplots(ncols=1, nrows=len(Spin_labels), figsize=(single_col_width, single_col_width*2.5), sharex=True, sharey=True)
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
                
            for sink_id in L_sink.keys():
                axs.flatten()[plot_it].plot(np.array(L_sink[sink_id]).T[0]-Time_array[0], np.array(L_sink[sink_id]).T[1], linestyle = mach_ls[Mach_labels.index(mach_lab)], color=colors[Mach_labels.index(mach_lab)], linewidth=1)
            axs.flatten()[plot_it].plot(Time_array - Time_array[0], L_sink_tot, linestyle = mach_ls[Mach_labels.index(mach_lab)], color=colors[Mach_labels.index(mach_lab)], linewidth=3, label='Mach '+mach_lab)
        else:
            print("Couldn't open", single_pickle)
            
        if mach_lab == '0.0':
            axs.flatten()[plot_it].set_ylabel('$\Omega t_{ff}='+spin_lab+'$: L ($g\,cm^2/s$)')
        
    if spin_lab == '0.20':
        axs.flatten()[plot_it].legend()
    if spin_lab == '0.35':
        axs.flatten()[plot_it].set_xlabel('Time ($yr$)')
    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='both', which='minor', labelsize=font_size, right=True)
    axs.flatten()[plot_it].tick_params(axis='x', direction='in')
    axs.flatten()[plot_it].tick_params(axis='y', direction='in')

axs.flatten()[plot_it].set_xlim(left=0)
axs.flatten()[plot_it].set_ylim(bottom=0)
plt.savefig('L_sink_spin_vs_mach.png', bbox_inches='tight')
