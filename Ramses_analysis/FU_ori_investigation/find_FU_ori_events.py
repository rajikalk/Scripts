#Plots MESA luminsoties calculated for all sinks int eh G100 512 run

import matplotlib.pyplot as plt
import numpy as np
import mesaPlot
import os
import glob
from mpi4py.MPI import COMM_WORLD as CW
from scipy.interpolate import interp1d
m=mesaPlot.MESA()
p=mesaPlot.plot()
lsun = 3.828e26*1e7 # solar luminosity in erg
FU_temp = np.concatenate((np.zeros(20), np.ones(80)))
time_window = 100
window_linspace = np.linspace(1, time_window,num=100,endpoint=True)

rank = CW.Get_rank()
size = CW.Get_size()

sink_files = sorted(glob.glob('/data/scratch/troels/IMF_512/mesa/sink_*/LOGS'))
if rank == 0:
    print('sink_files:', sink_files)
max_age=150000
rit = -1
for sink_file in sink_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        try:
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
                end_time = age + time_window
                end_it = np.argmin(abs(age - end_time))
                import pdb
                pdb.set_trace()
                
        except:
            print("can't read file")
