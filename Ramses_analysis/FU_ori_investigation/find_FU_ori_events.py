#Plots MESA luminsoties calculated for all sinks int eh G100 512 run

import matplotlib.pyplot as plt
import numpy as np
import mesaPlot
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
import scipy.interpolate as interp

m=mesaPlot.MESA()
p=mesaPlot.plot()
lsun = 3.828e26*1e7 # solar luminosity in erg
FU_temp = np.concatenate((np.zeros(20), np.ones(80)))
time_window = 100

rank = CW.Get_rank()
size = CW.Get_size()

sink_files = sorted(glob.glob('/data/scratch/troels/IMF_512/mesa/sink_*/LOGS'))
#if rank == 0:
#    print('sink_files:', sink_files)
max_age=150000
rit = -1
for sink_file in sink_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        #try:
        m.log_fold=sink_file
        m.loadHistory()
        mass = m.hist.star_mass
        age = m.hist.star_age
        idx = np.where(age <= max_age)
        age = age[idx]
        lum = 10.**m.hist.log_L
        lacc = m.hist.extra_lum / lsun
        ltot = lum + lacc
        
        #have moving window:
        for time_it in range(len(age)):
            end_time = age[time_it] + time_window
            end_it = np.argmin(abs(age - end_time))
            useable_times = age[time_it:end_it]
            useable_L = ltot[time_it:end_it]
            if len(useable_L) > 0:
                L_diff = np.max(np.log10(useable_L)) - np.min(np.log10(useable_L))
            else:
                L_diff = 0
            scaled_T = useable_times - useable_times[0]
            scaled_L = useable_L - np.min(useable_L)
            scaled_L = scaled_L/np.max(scaled_L)
            cor = np.correlate(scaled_L,FU_temp,'same')
            if np.max(cor)>0.8 and L_diff>1:
                plt.clf()
                plt.plot(useable_times, scaled_L, label="scaled Luminosity")
                plt.plot(useable_times, cor, label="correlation")
                plt.xlabel('Time (yr)')
                plt.ylabel('scaled L and correlation')
                plt.legend()
                plt.savefig('Sink_' + sink_file.split('sink_')[-1].split('/')[0] + '_time_'+str(age[time_it])+'png')
                print("Found potential match for sink", sink_file.split('sink_')[-1].split('/')[0], "at age", age[time_it])
        #except:
        #    print("can't read file")
