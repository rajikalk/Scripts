import subprocess
import glob
from mpi4py.MPI import COMM_WORLD as CW

rank = CW.Get_rank()
size = CW.Get_size()

data_dirs = sorted(glob.glob('data/output*/'))
rit = -1
for data_dir in data_dirs:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        dir_it = int(data_dir.split('output_')[-1][:-1])
        print('Compressing output', dir_it, 'on rank', rit)
        subprocess.run('tar -cvrfr output_'+("%05d" % dir_it)+'.tar.gz data/output_'+("%05d" % dir_it)+'/', shell=True)
        subprocess.run('mdss put output_'+("%05d" % dir_it)+'.tar.gz rlk100/2023_paper/Simulation_data/', shell=True)
        subprocess.run('rm output_'+("%05d" % dir_it)+'.tar.gz', shell=True)
        result = subprocess.run('mdss ls rlk100/2023_paper/Simulation_data/', shell=True, capture_output=True, text=True)
        output = result.stdout.strip()
        output = output.split('\n')
        if 'output_'+("%05d" % dir_it)+'.tar.gz' in output:
            subprocess.run('rm -rf data/output_'+("%05d" % dir_it)+'/*.out0*', shell=True)
    
