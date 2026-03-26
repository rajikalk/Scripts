#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--additions", help="additional inputs", default='', type=str)
parser.add_argument("files", nargs='*')
args = parser.parse_args()

Sim_dirs = [x[0] for x in os.walk('/home/100/rlk100/rlk/RAMSES/Zoom-in/')]
Cleaned_dirs = []
for sim_dir in Sim_dirs:
    if sim_dir[-4:] == 'data':
        Cleaned_dirs.append(sim_dir)
del Sim_dirs
save_dir = '/home/100/rlk100/rlk/RAMSES/Movie_frames/'

dirs = ['XY', 'XZ', 'YZ']

script_lines = ['#!/bin/bash',
                '#PBS -P ek9',
                '#PBS -q expressbw',
                '#PBS -l walltime=24:00:00',
                '#PBS -l ncpus=28',
                '#PBS -l mem=252GB',
                '#PBS -l storage=scratch/ek9+gdata/ek9',
                '#PBS -l wd',
                '#PBS -N ',
                '#PBS -j oe',
                '#PBS -m a',
                '#PBS -M rajika.kuruwita@anu.edu.au',
                '',
                'mpirun -np $PBS_NCPUS ~/Scripts/Ramses_analysis/movie_script.py ',]
#                ' ./ -sink ',
#                ' -sf 0 -dt 50 -cmin 1.e-18 -cmax 1.e-15 -at True -pvl True -ax ',
#                ' -al 1000 -tf 12 -stdv 5 -thickness 2000 -use_gas False -ic 1 -update_alim True -frames_only False -apm True 1>',
#                '.out01 2>&1']

for zoom_dir in Cleaned_dirs:
    Sink_id = int(zoom_dir.split('/Sink_')[-1].split('/')[0])
    job_name = str(Sink_id)
    movie_dir = save_dir + 'Sink_' + str(Sink_id)
    if os.path.exists(movie_dir) == False:
        os.makedirs(movie_dir)
    os.chdir(movie_dir)
    
    if 'Event_restart' in zoom_dir:
        movie_dir = movie_dir + '/Event_Restart'
        job_name = job_name + 'E'
        if os.path.exists(movie_dir) == False:
            os.makedirs(movie_dir)
        os.chdir(movie_dir)
    
    if 'Level' not in zoom_dir:
        movie_dir = movie_dir +'/Level_18'
        job_name = job_name + 'L18'
        if os.path.exists(movie_dir) == False:
            os.makedirs(movie_dir)
    else:
        zoom_level = zoom_dir.split('/Level_')[-1].split('/')[0]
        movie_dir = movie_dir +'/Level_' + zoom_level
        job_name = job_name + 'L'+zoom_level
        if os.path.exists(movie_dir) == False:
            os.makedirs(movie_dir)
    
    os.chdir(movie_dir)
    
    for dir in dirs:
        sub_movie_dir = movie_dir + '/' + dir
        if os.path.exists(sub_movie_dir) == False:
            os.makedirs(sub_movie_dir)
        os.chdir(sub_movie_dir)
        job_name_dir = dir+ job_name
        
        #Check that job script exists
        if os.path.isfile('job.sh') == False:
            f = open('job.sh', 'a')
            f.close()
        f = open('job.sh', 'w')
        for line in script_lines:
            if line == '#PBS -N ':
                write_line = line + job_name_dir
            elif line == 'mpirun -np $PBS_NCPUS ~/Scripts/Ramses_analysis/movie_script.py ':
                write_line = 'mpirun -np $PBS_NCPUS ~/Scripts/Ramses_analysis/movie_script.py '+zoom_dir+'/ ./ -sink '+str(Sink_id)+' -sf 0 -dt 50 -cmin 1.e-18 -cmax 1.e-15 -at True -pvl True -ax '+dir.lower()+' -al 1000 -tf 12 -stdv 5 -thickness 2000 -use_gas False -ic 1 -update_alim True -frames_only False -apm True 1>'+job_name+'.out01 2>&1'
            else:
                write_line = line
            write_line = write_line + '\n'
            f.write(write_line)
        f.close()
        
        #Check if job is already submitted/running
        qstat = subprocess.run(["qstat"], capture_output=True, text=True, check=True)
        if job_name_dir not in qstat.stdout:
            shellcmd = 'qsub job.sh'
            subprocess.call(shellcmd, shell=True)
            print('submitted jobs for', sub_movie_dir)
        os.chdir(movie_dir)
