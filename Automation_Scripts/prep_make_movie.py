#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--additions", help="additional inputs", default='', type=str)
parser.add_argument("files", nargs='*')
args = parser.parse_args()

Sim_dirs = [x[0] for x in os.walk(zoom_directory)]
Cleaned_dirs = []
for sim_dir in Sim_dirs:
    if sim_dir[-4:] == 'data':
        Cleaned_dirs.append(sim_dir)
del Sim_dirs
save_dir = '/home/100/rlk100/rlk/RAMSES/Movie_frames'

dirs = ['XY', 'XZ', 'YZ']

script_lines = ['#!/bin/bash',
                '#PBS -P ek9',
                '#PBS -q express',
                '#PBS -l walltime=24:00:00',
                '#PBS -l ncpus=1',
                '#PBS -l mem=9GB',
                '#PBS -l wd',
                '#PBS -N ',
                '#PBS -j oe',
                '#PBS -m bea',
                '#PBS -M rajika.kuruwita@anu.edu.au',
                '',
                'mpirun -np $PBS_NCPUS ~/Scripts/Ramses_analysis/movie_script.py /home/100/rlk100/rlk/RAMSES/Zoom-in/Sink_17/data/ ./ -sink ',
                ' -sf 0 -dt 50 -cmin 1.e-18 -cmax 1.e-15 -at True -pvl True -ax ',
                ' -al 1000 -tf 12 -stdv 5 -thickness 2000 -use_gas False -ic 1 -update_alim True -frames_only False -apm True 1>',
                '.out01 2>&1']
                
import pdb
pdb.set_trace()

f = open('job.sh', 'r')
reader = f.readlines()
sim_name = reader[7].split(' ')[-1].split('\n')[0]

for dir in dirs:
    job_name = sim_name+dir
    if not os.path.exists(save_dir + dir + '/'):
        os.makedirs(save_dir + dir + '/')
    movie_ls = subprocess.getoutput('ls ' + save_dir+dir + '/').split('\n')
    if len(movie_ls) == 1:
        start_frame = '0'
    else:
        movie_ls = subprocess.getoutput('ls ' + save_dir+dir+'/*jpg').split('\n')
        if "No such file or directory" in movie_ls[0]:
            movie_ls = subprocess.getoutput('ls ' + save_dir+dir+'/*eps').split('\n')
        start_frame = movie_ls[-1].split('_')[-1].split('.')[0]
    script_name = 'movie_'+dir+'.sh'
    if not os.path.isfile(script_name):
        f = open(script_name, 'a')
        f.close()
    f = open(script_name, 'w')
    for line in script_lines:
        if line == '#PBS -N ':
            write_line = line + job_name
        elif line == 'mpirun -np $PBS_NCPUS python ~/Scripts/movie_script_mod.py ':
            if dir == 'xz':
                write_line = line + cur_dir + ' ' + save_dir + dir + '/ -sf ' + start_frame + ' ' + args.additions + '-zt 1.33 -cmin 1.e-16 -cmax 1.e-14 -at True -pvl True -al 1000 -tf 16 ' + '1>' + job_name + '.out 2>&1'
            if dir == 'xy':
                write_line = line + cur_dir + ' ' + save_dir + dir + '/ -sf ' + start_frame + ' -ax "xy" ' + args.additions + ' ' + '1>' + job_name + '.out 2>&1'
            if dir == 'xz_zoom':
                write_line = line + cur_dir + ' ' + save_dir + dir + '/ -sf ' + start_frame + ' -zt 12.3 ' + args.additions + ' ' + '1>' + job_name + '.out 2>&1'
            if dir == 'xy_zoom':
                write_line = line + cur_dir + ' ' + save_dir + dir + '/ -sf ' + start_frame + '-zt 12.3 -cmin 1.e-15 -cmax 1.e-13 -at True -pvl True -ax xy -al 100 -tf 16 ' + args.additions + ' ' + '1>' + job_name + '.out 2>&1'
        else:
            write_line = line
        write_line = write_line + '\n'
        f.write(write_line)




