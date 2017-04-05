#!/usr/bin/env python
import numpy as np
import h5py
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import glob
import yt
import os
import argparse
import my_module as mym
from mpi4py import MPI

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="What is the input directory?", required=True)
    parser.add_argument("-o", "--output_dir", help="what is the output directory? if this is not defined it will be the same as the input directory")
    parser.add_argument("-of", "--output_file", help="Do you have a file ot write out to?", default="out.csv")
    parser.add_argument("-ff", "--file_frequency", help="every how many files do you want to use?", default=100, type=int)
    parser.add_argument("-dt", "--time_step", help="How frequently do you want to measure the outflow?", default=2, type=int)
    parser.add_argument("-disk", "--measure_disks", help="Do you want to measure disk quantities?", type=bool, default=False)
    #parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#===========================================================================================================================
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

args = parse_inputs()
output_file = args.output_file

#define files to read
dir = args.input_dir
if dir[-1] != '/':
    dir = dir + '/'
if args.output_dir == None:
    split_cur = dir.split('/O')
    save_dir = split_cur[0] + '/YT_O' + split_cur[1]
else:
    save_dir = args.output_dir
if save_dir[-1] != '/':
    save_dir = save_dir + '/'

file_name = 'WIND_hdf5_plt_cnt_*'
files = sorted(glob.glob(dir + file_name))

#define analysis cylinder:
if args.measure_disks != False:
    radius = 20.0
    height = 20.0
    tube_center_1 = None
    tube_center_2 = None
else:
    radius = 500.0
    height = 500.0
    upper_cylinder_bound = 250.0
    lower_cylinder_bound = -250.0
    center_1 = [0.0, 0.0, (upper_cylinder_bound+(height/2.))*yt.units.AU.in_units('cm')]
    center_2 = [0.0, 0.0, (lower_cylinder_bound-(height/2.))*yt.units.AU.in_units('cm')]

#read in sink creation time
part_file = files[-1][:-12] + 'part' + files[-1][-5:]
ds = yt.load(files[-1], particle_filename=part_file)
dd = ds.all_data()
sink_form = np.min(dd['particle_creation_time']/yt.units.yr.in_units('s')).value

#find times:
if rank == 0:
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None)
    print "m_times =", m_times
    usable_files = mym.find_files(m_times, files)
    for other_rank in range(1, size):
        comm.send(usable_files, dest=other_rank, tag=1)
        print "Sent usable files to rank", other_rank
else:
    usable_files = comm.recv(source=0, tag=1)
    print "Received usable files to rank", rank

print "usable_files =", usable_files

#open files to read data out to:
if rank == 0:
    if os.path.isfile(save_dir + output_file) == False:
        f = open(save_dir + output_file, 'a+')
        f.close()
    f = open(save_dir + output_file, 'w')
    if args.measure_disks != False:
        f.write('Lref ' + dir.split('_')[-1] + ': Time, Mass, Momentum, Angular Momentum, Max speed \n')
    else:
        f.write('Lref ' + dir.split('_')[-1] + ': Time, Lz_1, Lz_2 \n')
    f.close()

file_no = 0
rit = 1
for file in usable_files:
    if rank == rit:
        part_file = file[:-12] + 'part' + file[-5:]
        ds = yt.load(file, particle_filename=part_file)
        dd = ds.all_data()
        
        if args.measure_disks != False:
            center_1 = [dd['particle_posx'][0].in_units('cm'), dd['particle_posy'][0].in_units('cm'), dd['particle_posz'][0].in_units('cm')]
            center_2 = [dd['particle_posx'][1].in_units('cm'), dd['particle_posy'][1].in_units('cm'), dd['particle_posz'][1].in_units('cm')]
        
        time_val = ds.current_time.in_units('yr').value - sink_form

        #define cylinders:
        tube_1 = ds.disk(center_1, [0.0, 0.0, 1.0], (radius, 'au'), (height/2., 'au'))
        tube_2 = ds.disk(center_2, [0.0, 0.0, 1.0], (radius, 'au'), (height/2., 'au'))
        
        #calculate mass in cylinders
        # for region 1
        if args.measure_disks != False:
            m1 = np.sum(tube_1['cellmass'].in_units("Msun").value)
            Lz_2 = np.sum(tube_2['angular_momentum_z'].in_units("Msun * km**2 / s").value)
            print "OUTPUT=", time_val, Lz_1, Lz_2, "on rank", rit
                    
            #send data to rank 0 to append to write out.
            write_data = [time_val, Lz_1, Lz_2]
        else:
            pos_pos = np.where(tube_1['velocity_z'] > 0.0)[0]
            if len(pos_pos) == 0:
                mass_1 = 0.0
                speed_1 = 0.0
                max_speed_1 = 0.0
                mom_1 = 0.0
            else:
                mass_1 = tube_1['cell_mass'][pos_pos].in_units('Msun').value
                speed_1 = tube_1['velocity_magnitude'][pos_pos].in_units("km/s").value
                max_speed_1 = np.max(speed_1)
                mom_1 = speed_1*mass_1
            
            #for region 2
            neg_pos = np.where(tube_2['velocity_z'] < 0.0)[0]
            if len(neg_pos) == 0:
                mass_2 = 0.0
                speed_2 = 0.0
                max_speed_2 = 0.0
                mom_2 = 0.0
            else:
                mass_2 = tube_2['cell_mass'][neg_pos].in_units('Msun').value
                speed_2 = tube_2['velocity_magnitude'][neg_pos].in_units("km/s").value
                max_speed_2 = np.max(speed_2)
                mom_2 = speed_2*mass_2
            
            #Calculate outflow mass:
            outflow_mass = np.sum(mass_1) + np.sum(mass_2)
            
            #Calculate max outflow speed:
            if max_speed_1 > max_speed_2:
                max_speed = abs(max_speed_1)
            else:
                max_speed = abs(max_speed_2)

            #Calculate outflow momentum:
            mom = np.sum(mom_1) + np.sum(mom_2)

            #Calculate outflow angular momentum:
            L_x = np.sum(tube_1['angular_momentum_x'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['angular_momentum_x'][neg_pos].in_units("Msun * km**2 / s").value)
            L_y = np.sum(tube_1['angular_momentum_y'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['angular_momentum_y'][neg_pos].in_units("Msun * km**2 / s").value)
            L_z = np.sum(tube_1['angular_momentum_z'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['angular_momentum_z'][neg_pos].in_units("Msun * km**2 / s").value)
            L = np.sqrt((L_x)**2. + (L_y)**2. + (L_z)**2.)
            
            #Append quantities
            if outflow_mass == 0.0:
                outflow_mass = np.nan
                max_speed = np.nan
                mom = np.nan
                L = np.nan
            
            print "OUTPUT=", time_val, outflow_mass, mom, L, max_speed, "on rank", rit

            #send data to rank 0 to append to write out.
            write_data = [time_val, outflow_mass, max_speed, mom, L]
        comm.send(write_data, dest=0, tag=2)
        print "Sent data", write_data, "from rank", rank

    rit = rit + 1
    if rit == size:
        rit = 1

    if rank == 0:
        if (len(usable_files)-file_no) < size:
            end_rank = (len(usable_files)-file_no)+1
        else:
            end_rank = size
        for other_rank in range(1,end_rank):
            write_data = comm.recv(source=other_rank, tag=2)
            print "received data", write_data, "from rank", other_rank
            f = open(save_dir + output_file, 'a')
            write_string = ''
            for datum in write_data:
                if datum != write_data[-1]:
                    write_string = write_string + str(datum) + ','
                else:
                    write_string = write_string + str(datum) + '\n'
            f.write(write_string)
            f.close()

            file_no = file_no + 1
            print "Progress:", file_no, "/", len(usable_files)
        if file_no == len(usable_files):
            print "BREAKING ON RANK", rank
            break

    #print "file =", file
    if file == usable_files[-1]:
        break

    if len(usable_files) < size:
        if rank >= size:
            break

print "Completed job on rank", rank