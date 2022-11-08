#!/usr/bin/env python
import os
import glob
#from subprocess import Popen
import subprocess
#from mpi4py.MPI import COMM_WORLD as CW
import sys

#get mpi size and ranks
#rank = CW.Get_rank()
#size = CW.Get_size()
size=sys.argv[1]

sim_dirs = [x[0] for x in os.walk('/hits/fast/set/kuruwira/Protostellar_spin')]

clean_pickles = False
clean_images = False

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
        run_line = 'mpirun -np ' + str(size) + ' /home/kuruwira/Scripts/FLASH_analysis/movie_script.py ' + sim_dir +'/ '
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
                 
                #proc = Popen(proj_run_line, shell=True)
                #subprocess.run('module list', shell=True)
                subprocess.run('source ~/.bashrc', shell=True)
                #subprocess.run('module list', shell=True)
                
                subprocess.run(proj_run_line, shell=True)
                
                #check all frames were made:
                if len(glob.glob(save_dir + '*.pkl')) != len(glob.glob(save_dir + '*.jpg')):
                    #proc = Popen(proj_run_lines, shell=True)
                    subprocess.run(proj_run_line, shell=True)
