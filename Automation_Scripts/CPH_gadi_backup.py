#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Christoph Federrath
import subprocess
import csv

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-loc", "--location", help="which server location", default='gadi')
parser.add_argument('files', nargs='*')
args = parser.parse_args()

# ===== MAIN Start =====
try:
    files_done = []
    with open('Copied_dirs_'+args.location+'.txt', 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            files_done.append(row[0])
    file_no = int(files_done[-1].split("_")[-1])
    file_no = file_no + 10
except:
    files_done = []
    file_no = 160

while file_no < 1560:
    curr_file = ("%05d" % file_no)
    if args.location == 'gadi':
        shellcmd = 'scp -r astro03-travel:/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+' gadi:/g/data/ek9/rlk100/RAMSES/Sink_45/data/.'
    elif args.location == 'setonix':
        shellcmd = 'scp -r astro03-travel:/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+' setonix:/home/rkuruwita1/rlk/RAMSES/Zoom-ins/Sink_45/data/.'
    return_code = subprocess.call(shellcmd, shell=True)
    if return_code == 0:
        files_done.append('/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file)
        with open('Copied_dirs_'+args.location+'.txt', 'a') as f:
            f.write(files_done[-1])
        f.close()
        file_no = file_no + 10
    else:
        print("Error while tranferring!")
        break

