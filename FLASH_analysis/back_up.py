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

run_line = ['rsync -vaz /home/kuruwira/fast/Protostellar_spin/Spin_0.2/Single/Lref_10/Mach_0.2/* backupset:/home/kuruwira/Protostellar_spin/Spin_0.20/Single/.', 'rsync -vaz /home/kuruwira/fast/Protostellar_spin/Spin_0.25/Single/Lref_10/Mach_0.2/* backupset:/home/kuruwira/Protostellar_spin/Spin_0.25/Single/.', 'rsync -vaz /home/kuruwira/fast/Protostellar_spin/Spin_0.2/Binary/Lref_10/Restart/Mach_0.2/From_binary_formation/* backupset:/home/kuruwira/Protostellar_spin/Spin_0.20/Binary/.', 'rsync -vaz /home/kuruwira/fast/Protostellar_spin/Spin_0.25/Binary/Lref_10/Restart/Mach_0.2/From_binary_formation/* backupset:/home/kuruwira/Protostellar_spin/Spin_0.25/Binary/.']

update = True

rit = -1
for sim_dir in run_line:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        subprocess.run(run_line, shell=True)
            
