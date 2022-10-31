#!/usr/bin/env python
import os
import glob
import subprocess

sim_dirs = [x[0] for x in os.walk('/hits/fast/set/kuruwira/Protostellar_spin')]
n_processes = 20

for sim_dir in sim_dirs:
    if len(glob.glob(sim_dir + '/sinks_evol.dat')) > 0:
        #check if movie directory exists
        movie_dir = '/hits/fast/set/kuruwira/Movie_frames' + sim_dir.split('Protostellar_spin')[-1]
        if os.path.exists(movie_dir) == False:
            #make movie directory
            os.makedirs(movie_dir)
        
        #check if movie is up to date? How... I guess just run the movie line
        run_line = 'mpirun -np ' + str(n_processes) + '/home/kuruwira/Scripts/FLASH_analysis/movie_script.py ' + sim_dir +'/ '
        
        proj_dirs = ['/XY/', '/XZ/']
        for proj_dir in proj_dirs:
            save_dir = movie_dir + proj_dir
            if os.path.exists(ç) == False:
                os.makedirs(save_dir)
            
            proj_run_line = run_line + save_dir
            if proj_dir == '/XZ/':
                proj_run_line = proj_run_line + " -ax 'y'"
                
            subprocess.run(proj_run_line, shell=True)
            
            #check all frames were made:
            if len(glog.glob(save_dir + '*.pkl') != len(glog.glob(save_dir + '*.jpg'):
                subprocess.run(proj_run_line, shell=True)
