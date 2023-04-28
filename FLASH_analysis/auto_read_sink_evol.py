#!/usr/bin/env python
import os
import glob
import subprocess
from mpi4py.MPI import COMM_WORLD as CW
import sys

#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()
#size=sys.argv[1]

sim_dirs = [x[0] for x in os.walk('/home/kuruwira/fast/Protostellar_spin/Flash_2023/')]

update = True

rit = -1
for sim_dir in sim_dirs:
    if os.path.exists(sim_dir+'/sinks_evol.dat') == True:
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            pickle_name = "_".join(sim_dir.split('Protostellar_spin/')[-1].split('/')) + '.pkl'
            if update == False and os.path.exists(pickle_name) == False:
                proj_run_line = 'python3 /home/kuruwira/Scripts/FLASH_analysis/read_sink_evol.py '+sim_dir+'/sinks_evol.dat /home/kuruwira/fast/Analysis/Sink_evol_pickles/'+pickle_name
        
                subprocess.run(proj_run_line, shell=True)
            elif update == True:
                proj_run_line = 'python3 /home/kuruwira/Scripts/FLASH_analysis/read_sink_evol.py '+sim_dir+'/sinks_evol.dat /home/kuruwira/fast/Analysis/Sink_evol_pickles/'+pickle_name
        
                subprocess.run(proj_run_line, shell=True)
            
