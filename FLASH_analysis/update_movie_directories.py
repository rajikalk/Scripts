#!/usr/bin/env python
import os
import glob
#from subprocess import Popen
import subprocess
#from mpi4py.MPI import COMM_WORLD as CW
import sys
import argparse

#get mpi size and ranks
#rank = CW.Get_rank()
#size = CW.Get_size()
size=20

sim_dirs = [x[0] for x in os.walk('/hits/fast/set/kuruwira/Protostellar_spin')]

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-make_pickles", "--make_pickles", help="do you want to delete the pickles?", default="False")
    parser.add_argument("-del_images", "--delete_images", type=str, default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#-------------------------------------------------------
#get input and output directory and arguments
args = parse_inputs()

clean_pickles = eval(args.make_pickles)
clean_images = eval(args.delete_images)

proj_dirs = ['/XY/', '/XZ/']
zoom_dirs = ['1000AU/', '250AU/']

for sim_dir in sim_dirs:
    if len(glob.glob(sim_dir + '/sinks_evol.dat')) > 0:
        #check if movie directory exists
        movie_dir = '/hits/fast/set/kuruwira/Movie_frames' + sim_dir.split('Protostellar_spin')[-1]
        
        #sys.stdout.flush()
        #CW.Barrier()
        
        if os.path.exists(movie_dir) == False:# and rank == 0:
            #make movie directory
            os.makedirs(movie_dir)
        
        #sys.stdout.flush()
        #CW.Barrier()
        
        #check if movie is up to date? How... I guess just run the movie line
        #if size > 1:
        run_line = 'mpirun -np ' + str(size) + ' python3 /home/kuruwira/Scripts/FLASH_analysis/movie_script.py ' + sim_dir +'/ '
        #else:
        #    run_line = 'python /home/kuruwira/Scripts/FLASH_analysis/movie_script.py ' + sim_dir +'/ '
            
        
        for proj_dir in proj_dirs:
            for zoom_dir in zoom_dirs:
                save_dir = movie_dir + proj_dir
                if zoom_dir == '250AU/':
                    save_dir = save_dir + zoom_dir
                if os.path.exists(save_dir) == False:
                    os.makedirs(save_dir)
                    
                if clean_pickles:# and rank == 0:
                    for pickle_file in glob.glob(save_dir + '*.pkl'):
                        os.remove(pickle_file)
                
                #sys.stdout.flush()
                #CW.Barrier()
                
                if clean_images:# and rank == 0:
                    for image_file in glob.glob(save_dir + '*.jpg'):
                        os.remove(image_file)
                
                #sys.stdout.flush()
                #CW.Barrier()
                
                proj_run_line = run_line + save_dir
                if proj_dir == '/XZ/':
                    proj_run_line = proj_run_line + " -ax 'y'"
                if zoom_dir == '250AU/':
                    proj_run_line = proj_run_line + " -width 500"
                 
                job_id = ''
                if 'Single' in save_dir:
                    job_id = job_id + 'S'
                else:
                    job_id = job_id + 'B'
                job_id = job_id + save_dir.split('Spin_0.')[-1].split('/')[0]
                job_id = job_id + save_dir.split('Mach_0.')[-1].split('/')[0]
                job_id = job_id + proj_dir[2]
                
                job_id = job_id + save_dir.split('Lref_')[-1].split('/')[0]
                
                job_id = job_id + zoom_dir[:-3]
                
                
                f = open(save_dir+'movie.sh', 'w')
                
                f.write('#!/bin/bash\n')
                f.write('#SBATCH --job-name='+job_id+'       # shows up in the output of squeue\n')
                f.write('#SBATCH --partition=cascade.p   # specify the partition to run on\n')
                f.write('#SBATCH --time=24:00:00         # specify the requested wall-time\n')
                f.write('#SBATCH --nodes=1               # number of nodes allocated for this job\n')
                f.write('#SBATCH --ntasks-per-node=20    # number of MPI ranks per node\n')
                f.write('#SBATCH --cpus-per-task=1       # number of OpenMP threads per MPI rank\n')
                f.write('#SBATCH --mail-type=ALL  # NONE, ALL, FAIL, END...\n')
                f.write('#SBATCH --mail-user=rajika.kuruwita@h-its.org\n')
                f.write('#SBATCH --gres=cpuonly\n')

                f.write('source ~/.bashrc\n')
                f.write('chmod a+x /home/kuruwira/Scripts/FLASH_analysis/movie_script.py\n')

                f.write(proj_run_line+ ' 1>frames.out00 2>&1\n')
                f.close()
                
                os.chdir(save_dir)
                subprocess.run('sbatch movie.sh', shell=True)
                os.chdir('/hits/fast/set/kuruwira/Movie_frames')
                '''
                #proc = Popen(proj_run_line, shell=True)
                #subprocess.run('module list', shell=True)
                subprocess.run('source ~/.bashrc', shell=True)
                #subprocess.run('module list', shell=True)
                
                subprocess.run(proj_run_line, shell=True)
                
                #check all frames were made:
                if len(glob.glob(save_dir + '*.pkl')) != len(glob.glob(save_dir + '*.jpg')):
                    #proc = Popen(proj_run_lines, shell=True)
                    subprocess.run(proj_run_line, shell=True)
                '''
