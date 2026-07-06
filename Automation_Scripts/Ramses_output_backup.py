#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Christoph Federrath
import subprocess
import csv
import glob
import os
import sys
from mpi4py.MPI import COMM_WORLD as CW

# ===== MAIN Start =====
rank = CW.Get_rank()
size = CW.Get_size()

#Open and read Expiring text file:
rit = -1
with open('/home/100/rlk100/expiring_files.txt', 'r') as f:
    reader = csv.reader(f)
    header = True
    for row in reader:
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            if 'output_' in row[0]:
                compress_file = row[0].split(' ')[-1].split('output_')[0] + 'output_' + row[0].split(' ')[-1].split('output_')[1].split('/')[0]
                tar_gz = compress_file+'.tar.gz'
                #check if mdss save dir is exists
                mdss_save_dir = 'rlk100/Zoom_in_simulations/Sink_' + compress_file.split('Sink_')[-1].split('data')[0]
                shellcmd = 'mdss ls ' + mdss_save_dir
                result = subprocess.run(shellcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                output_string = result.stdout
                if "No such file or directory" in output_string:
                    shellcmd = 'mdss mkdir ' + mdss_save_dir
                    result = subprocess.run(shellcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    output_string = result.stdout
                    print("RANK ", rit, ": Created MDSS directory", mdss_save_dir)
                    sys.stdout.flush()
                #check if tar in mdss
                shellcmd = 'mdss ls ' + mdss_save_dir + tar_gz.split('/')[-1]
                result = subprocess.run(shellcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                output_string = result.stdout
                if "No such file or directory" in output_string:
                    #check if tar exists
                    if os.path.isfile(tar_gz) ==  False:
                        shellcmd = 'tar -czvf ' + tar_gz + ' ' + compress_file
                        result = subprocess.run(shellcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                        output_string = result.stdout
                        print("RANK ", rit, ": Made tar of", compress_file)
                        sys.stdout.flush()
                    #now move to MDSS
                    shellcmd = 'mdss put '+ tar_gz + ' '+mdss_save_dir
                    result = subprocess.run(shellcmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                    output_string = result.stdout
                    print("RANK ", rit, ": put "+tar_gz+' '+mdss_save_dir)
                    sys.stdout.flush()
                    #remove tar file
                    os.remove(tar_gz)
f.close()
