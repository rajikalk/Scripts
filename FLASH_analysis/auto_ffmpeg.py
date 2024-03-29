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

sim_dirs = [x[0] for x in os.walk('/home/kuruwira/fast/Movie_frames')]

update = False

rit = -1
for sim_dir in sim_dirs:
    if len(glob.glob(sim_dir + '/*.jpg')) > 0 and 'Single' in sim_dir and 'Lref_9' in sim_dir and 'XY' in sim_dir and '250' in sim_dir:
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            if len(glob.glob(sim_dir+'/*.jpg')) > 0:
                #make movie!
                movie_name = "_".join(sim_dir.split('Movie_frames/')[-1].split('/')) + '.mp4'
                if update == False and os.path.exists(movie_name) == False:
                    tmp_list = 'jpgtmpfilelist_'+str(rank)+'.txt'
                    proj_run_line = 'python3 /home/kuruwira/Scripts/Automation_Scripts/make_movie.py -o /home/kuruwira/fast/Movies/Protostellar_spin/'+movie_name +' -tmp_list '+tmp_list +' '+sim_dir + '/*.jpg'
            
                    subprocess.run(proj_run_line, shell=True)
                elif update == True:
                    tmp_list = 'jpgtmpfilelist_'+str(rank)+'.txt'
                    proj_run_line = 'python3 /home/kuruwira/Scripts/Automation_Scripts/make_movie.py -o /home/kuruwira/fast/Movies/Protostellar_spin/'+movie_name +' -tmp_list '+tmp_list +' '+sim_dir + '/*.jpg'
            
                    subprocess.run(proj_run_line, shell=True)
                
