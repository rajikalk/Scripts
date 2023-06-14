#!/usr/bin/env python
import os
import glob
#from subprocess import Popen
import subprocess
#from mpi4py.MPI import COMM_WORLD as CW
import sys
import argparse
import shutil

#get mpi size and ranks
#rank = CW.Get_rank()
#size = CW.Get_size()
size=40

sim_dirs = [x[0] for x in os.walk('/hits/fast/set/kuruwira/Protostellar_spin')]

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-make_pickles", "--make_pickles", help="do you want to delete the pickles?", default="False")
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#-------------------------------------------------------
#get input and output directory and arguments
args = parse_inputs()

clean_pickles = eval(args.make_pickles)

for sim_dir in sim_dirs:
    if len(glob.glob(sim_dir + '/sinks_evol.dat')) > 0:
        #check if movie directory exists
        save_dir = '/hits/fast/set/kuruwira/Analysis/Angular_momentum_budget' + sim_dir.split('Protostellar_spin')[-1]
        
        #sys.stdout.flush()
        #CW.Barrier()
        
        if os.path.exists(save_dir) == False:# and rank == 0:
            #make movie directory
            os.makedirs(save_dir)
        
        #sys.stdout.flush()
        #CW.Barrier()
        
        #check if movie is up to date? How... I guess just run the movie line
        #if size > 1:
        run_line = 'mpirun -np ' + str(size) + ' python3 /home/kuruwira/Scripts/FLASH_analysis/angular_momentum_evolution.py ' + sim_dir +'/ '
        #else:
        #    run_line = 'python /home/kuruwira/Scripts/FLASH_analysis/movie_script.py ' + sim_dir +'/ '
            

        if len(glob.glob(sim_dir + '/*plt_cnt*')) == 0:
            shutil.rmtree(save_dir)
        
        if len(glob.glob(sim_dir + '/*plt_cnt*')) > 0:
            if clean_pickles:# and rank == 0:
                for pickle_file in glob.glob(save_dir + '*.pkl'):
                    os.remove(pickle_file)
            
            #sys.stdout.flush()
            #CW.Barrier()
            
            
            proj_run_line = run_line# + save_dir
             
            job_id = ''
            if 'Single' in save_dir:
                job_id = job_id + 'S'
            else:
                job_id = job_id + 'B'
            
            job_id = job_id + 'L'
            job_id = job_id + save_dir.split('Spin_0.')[-1].split('/')[0]
            job_id = job_id + save_dir.split('Mach_0.')[-1].split('/')[0]
            
            job_id = job_id + save_dir.split('Lref_')[-1].split('/')[0]
        
            try:
                os.remove(save_dir+'L_budget.sh')
            except:
                pass
            f = open(save_dir+'/L_budget.sh', 'w')
            
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --job-name='+job_id+'       # shows up in the output of squeue\n')
            f.write('#SBATCH --partition=cascade.p   # specify the partition to run on\n')
            f.write('#SBATCH --time=2:00:00         # specify the requested wall-time\n')
            f.write('#SBATCH --nodes=2               # number of nodes allocated for this job\n')
            f.write('#SBATCH --ntasks-per-node=20    # number of MPI ranks per node\n')
            f.write('#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n')
            f.write('#SBATCH --mail-type=BEGIN,END,FAIL  # NONE, ALL, FAIL, END...\n')
            f.write('#SBATCH --mail-user=rajika.kuruwita@h-its.org\n')
            f.write('#SBATCH --gres=cpuonly\n')

            f.write('source ~/.bashrc\n')
            f.write('chmod a+x /home/kuruwira/Scripts/FLASH_analysis/movie_script.py\n')

            f.write(proj_run_line+ ' 1>budget.out00 2>&1\n')
            f.close()
            
            os.chdir(save_dir)
            subprocess.run('sbatch L_budget.sh', shell=True)
            os.chdir('/hits/fast/set/kuruwira/Movie_frames')
