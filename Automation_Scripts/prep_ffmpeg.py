#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import subprocess
import argparse
import glob

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
                '#PBS -l walltime=1:00:00',
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
        if len(glob.glob(sub_movie_dir+'/*.jpg'))>0:
            shellcmd = 'python ~/Scripts/Automation_Scripts/make_movie.py -o '+job_name+'.mp4 movie_frame_000*.jpg'
            subprocess.call(shellcmd, shell=True)
            print('created movie '+job_name+dir+'.mp4')
        else:
            print("FRAMES DON'T EXIST FOR", job_name+dir)
        os.chdir(movie_dir)
