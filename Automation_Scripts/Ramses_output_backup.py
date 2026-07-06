#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Christoph Federrath
import subprocess
import csv
import glob
import os

# ===== MAIN Start =====

#Open and read Expiring text file:
with open('/home/100/rlk100/expiring_files.txt', 'r') as f:
    reader = csv.reader(f)
    header = True
    for row in reader:
        if 'output_' in row[0]:
            compress_file = row[0].split(' ')[-1].split('output_')[0] + 'output_' + row[0].split(' ')[-1].split('output_')[1].split('/')[0]
            tar_gz = compress_file+'.tar.gz'
            if
            shellcmd = 'tar -czvf ' + tar_gz + ' compress_file
            result = subprocess.run(shellcmd, stdout=subprocess.PIPE, shell=True, capture_output=True, text=True)
            output_string = result.stdout
            print("Made tar of", compress_file)
            
            mdss_save_dir = 'rlk100/Zoom_in_simulations/Sink_' + compress_file.split('Sink_')[-1].split('data')[0]
            #check is sim dir exists:
            shellcmd = 'mdss ls' + mdss_save_dir
            result = subprocess.run(shellcmd, stdout=subprocess.PIPE, shell=True, capture_output=True, text=True)
            output_string = result.stdout
            import pdb
            pdb.set_trace()
            proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
            shellcmd = 'mdss put '+ compress_file+'.tar.gz '+mdss_save_dir
            proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
            print("put"+compress_file+'.tar.gz '+mdss_save_dir)
            
            backup_done = []

            with open('Backup_dirs.txt', 'a') as f_backup:
                f_backup.write(compress_file +'\n')
            f_backup.close()


f.close()
