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
            import pdb
            pdb.set_trace()
            compress_file = row[0].split(' ')[-1].split('output_')[0] + 'output_' + row[0].split(' ')[-1].split('output_')[1].split('/')[0]
            shellcmd = 'tar -czvf ' + compress_file+'.tar.gz '+compress_file
            proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
            print("Made tar of", compress_file)
            
            mdss_save_dir = 'rlk100/Zoom_in_simulations/Sink_' + compress_file.split('Sink_')[-1].split('data')[0]
            #check is sim dir exists:
            shellcmd = 'mdss ls '+ mdss_save_dir
            proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
            output_mdss = proc.stdout.read()
            if "No such file or directory" in output_mdss:
                import pdb
                pdb.set_trace()
            
            
            
            backup_done = []
            with open('Backup_dirs.txt', 'r') as f_backup:
                reader = csv.reader(f_backup)
                for row in reader:
                    synced_done.append(row[0])
            if '/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file not in synced_done[-1]:
                shellcmd = 'rsync -vaz astro03-travel:/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+'/* output_'+curr_file+'/.'
                return_code = subprocess.call(shellcmd, shell=True)
                if return_code != 0:
                    print("error on copying file locally")
                    break
                else:
                    #update sync list
                    with open('Backup_dirs.txt', 'a') as f_backup:
                        f_backup.write('/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+'\n')
                    f_backup.close()


f.close()
    
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
    file_no = 150

while file_no < 1560:
    curr_file = ("%05d" % file_no)
    
    #
    shellcmd = "rsync -n -i -a astro03-travel:/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_"+curr_file+" output_"+curr_file+"/"
    proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
    output_Cph_local = proc.stdout.read()
    
    shellcmd = "rsync -n -i -a gadi:/g/data/ek9/rlk100/RAMSES/Sink_45/data/output_"+curr_file+" output_"+curr_file+"/"
    proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
    output_gadi_local = proc.stdout.read()
    #
    if args.location == 'gadi':
        #make empty file and copy to gadi
        if os.path.exists('output_'+curr_file) == False:
            os.makedirs('output_'+curr_file)
            shellcmd = 'rsync -vaz output_'+curr_file+' gadi:/g/data/ek9/rlk100/RAMSES/Sink_45/data/.'
            return_code = subprocess.call(shellcmd, shell=True)
            if return_code != 0:
                print("error with copying empty file")
                break
        
        #Check if CPH file has already been synced with local
        synced_done = []
        with open('Synced_dirs_'+args.location+'.txt', 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                synced_done.append(row[0])
        if '/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file not in synced_done[-1]:
            shellcmd = 'rsync -vaz astro03-travel:/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+'/* output_'+curr_file+'/.'
            return_code = subprocess.call(shellcmd, shell=True)
            if return_code != 0:
                print("error on copying file locally")
                break
            else:
                #update sync list
                with open('Synced_dirs_'+args.location+'.txt', 'a') as f:
                    f.write('/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+'\n')
                f.close()
                
        #move into directory and compare local to Gadi
        os.chdir('output_'+curr_file)
        #Check files on gadi:
        shellcmd='rsync -n -i -a . gadi:/g/data/ek9/rlk100/RAMSES/Sink_45/data/output_'+curr_file+'/'
        proc = subprocess.Popen(shellcmd, stdout=subprocess.PIPE, shell=True)
        output = proc.stdout.read()
        gadi_files = str(output).split('\\n')
        transfer_flag = '<f+++++++ '
        rsync_list = list(filter(lambda x: transfer_flag in x, gadi_files))
        transfer_files = []
        for trans_file in rsync_list:
            transfer_files.append(trans_file.split(transfer_flag)[-1])
        local_files = sorted(glob.glob('*'))
        rm_files = list(set(local_files).symmetric_difference(set(transfer_files)))
        
        #Remove files that are already synced
        
        for rm_file in rm_files:
            os.remove(rm_file)
            
        #sync remaining files
        shellcmd = 'rsync -vaz * gadi:/g/data/ek9/rlk100/RAMSES/Sink_45/data/output_'+curr_file+'/.'
        return_code = subprocess.call(shellcmd, shell=True)
        if return_code != 0:
            print("Error on ryncing to gadi!")
            break
        else:
            files_done.append('/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file)
            with open('Copied_dirs_'+args.location+'.txt', 'a') as f:
                f.write(files_done[-1]+'\n')
            f.close()
            os.chdir('../')
            shellcmd = 'rm -rf output_'+curr_file+'/'
            return_code = subprocess.call(shellcmd, shell=True)
            if return_code == 0:
                print('removed previous file and moving onto the next!')
            file_no = file_no + 10

    elif args.location == 'setonix':
        shellcmd = 'scp -r astro03-travel:/groups/astro/rlk/rlk/FU_ori_investigation/Zoom_in_simulations/Sink_45/data/output_'+curr_file+' setonix:/home/rkuruwita1/rlk/RAMSES/Zoom-ins/Sink_45/data/.'
        return_code = subprocess.call(shellcmd, shell=True)

