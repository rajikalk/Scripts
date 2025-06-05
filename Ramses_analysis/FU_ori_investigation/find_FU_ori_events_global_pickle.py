import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp
import pickle
import yt

lsun = 3.828e26*1e7 # solar luminosity in erg
FU_temp = np.concatenate((np.zeros(25), np.ones(75)))
time_window = yt.YTQuantity(80, 'yr')

global_pickle = "/groups/astro/rlk/rlk/FU_ori_investigation/Sink_pickles/particle_data_global.pkl"
file_open = open(global_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()

age = yt.YTArray(particle_data['time'])
lacc = yt.YTArray(particle_data['lacc']).T[1]
mass = yt.YTArray(particle_data['mass']).T[1]
mdot = yt.YTArray(particle_data['mdot']).T[1]
lstar = yt.YTArray(0.23*(mass**2.3).value, 'lsun')
facc = 0.5
lstar_offner = 31.3*facc*(mdot/1.e-6)
ltot = lacc + lstar


rank = CW.Get_rank()
size = CW.Get_size()
L_diff_arr = []
time_arr = []
cor_arr = []
for time_it in range(len(age)):
    end_time = age[time_it] + time_window
    end_it = np.argmin(abs(age - end_time))
    useable_times = age[time_it:end_it]
    useable_L = lacc[time_it:end_it]
    if len(useable_L) > 0:
        L_diff = np.max(useable_L)/np.min(useable_L)
        L_diff_arr.append(L_diff)
        time_arr.append(age[time_it])
        scaled_T = useable_times - useable_times[0]
        scaled_L = useable_L - np.min(useable_L)
        scaled_L = scaled_L/np.max(scaled_L)
        cor = np.correlate(scaled_L,FU_temp,'same')
        cor_arr.append(np.nanmax(cor))
        if np.median(cor)>66.6 and L_diff>10: #and mass[time_it] > 0.1
            plt.clf()
            fig, ax1 = plt.subplots()

            ax2 = ax1.twinx()
            ax1.plot(useable_times, scaled_L, label="scaled L_acc", color='b')
            ax1.plot(useable_times, cor[:len(useable_times)]/100., label="correlation", color='r')
            
            ax2.plot(useable_times, useable_L, color='b')

            ax1.set_xlabel('Time (yr)')
            ax1.set_ylabel('scaled L_acc and correlation')
            ax2.set_ylabel('Total log Luminosity')
            
            ax1.set_ylim([0, 1])
            ax2.set_ylim([np.min(useable_L), np.max(useable_L)])
        
            ax1.legend()
            plt.savefig('Sink_' + sink_file.split('sink_')[-1].split('/')[0] + '_time_'+str(age[time_it])+'_mass_'+str(np.round(mass[time_it], decimals=2))+'.png',  bbox_inches='tight')
            print("Found potential match for sink", sink_file.split('sink_')[-1].split('/')[0], "at age", age[time_it])
plt.clf()
plt.plot(time_arr, L_diff_arr)
plt.xlabel('age (yr)')
plt.ylabel('max L diff over 80 yr (log)')
plt.savefig('L_diff.png')
print("plotted L diff history for sink 45 on rank", rank)
