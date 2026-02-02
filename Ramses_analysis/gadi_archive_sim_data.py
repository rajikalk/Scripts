import subprocess
import glob
import sys
from mpi4py.MPI import COMM_WORLD as CW

rank = CW.Get_rank()
size = CW.Get_size()

data_dirs = sorted(glob.glob('data/output*/'))
rit = -1
for data_dir in data_dir:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        dir_it = int(data_dir.split('output_')[-1][:-1])
        subprocess.run('tar -cvfr output_'+("%05d" % dir_it)+'.tar.gz output_'+("%05d" % dir_it)+'/', shell=True)
        subprocess.run('mdss put output_'+("%05d" % dir_it)+'.tar.gz rlk100/2023_paper/Simulation_data/', shell=True)
        subprocess.run('rm output_'+("%05d" % dir_it)+'.tar.gz', shell=True)
        subprocess.run('rm -rf output_'+("%05d" % dir_it)+'/', shell=True)
    
