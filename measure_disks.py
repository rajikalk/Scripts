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
import my_fields as myf

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="What is the input directory?", required=True)
    parser.add_argument("-o", "--output_dir", help="what is the output directory? if this is not defined it will be the same as the input directory")
    parser.add_argument("-of", "--output_file", help="Do you have a file ot write out to?", default="out.csv")
    parser.add_argument("-ff", "--file_frequency", help="every how many files do you want to use?", default=100, type=int)
    parser.add_argument("-dt", "--time_step", help="How frequently do you want to measure the outflow?", default=2, type=int)
    parser.add_argument("-lb", "--lower_bounds", help="lower bounds of box", type=float, default=50.0)
    parser.add_argument("-rad", "--radius", help="What is the radius of the measuring cylindar", type=float, default=6000.0)
    parser.add_argument("-vthres", "--velocity_threshold", help="what velocity threshold do you want to use to outflowing material?", type=float, default=0.0)
    parser.add_argument("-center", "--disk_center", help="do you want two disks centered on the particles or a volume centered on the CoM?", default='CoM', type=str)
    parser.add_argument("-st", "--start_time", help="default is -500yr", default=-500.0, type=float)
    #parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#===========================================================================================================================

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
if args.disk_center == 'CoM':
    radius = 100.0
    height = 100.0
    tube_center_1 = None
else:
    radius = 20.0
    height = 20.0
    tube_center_1 = None
    tube_center_2 = None

#read in sink creation time
part_file = files[-1][:-12] + 'part' + files[-1][-5:]
ds = yt.load(files[-1], particle_filename=part_file)
dd = ds.all_data()
sink_form = np.min(dd['particle_creation_time']/yt.units.yr.in_units('s')).value

#find times:
m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None, start_time=args.start_time)
usable_files = mym.find_files(m_times, files)

print "usable_files =", usable_files

#open files to read data out to:
if os.path.isfile(save_dir + output_file) == False:
    f = open(save_dir + output_file, 'a+')
    f.close()
f = open(save_dir + output_file, 'w')
f.write('Lref ' + dir.split('_')[-1] + ': Time, Disk_mass \n')
f.close()

file_no = 0
for file in usable_files:
    part_file = file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    dd = ds.all_data()
    dist = np.sqrt(dd['CoM'][0]**2. + dd['CoM'][1]**2. + dd['CoM'][2]**2.).in_units('AU').value
    print "CoM dist =", dist, "(AU)"
    

    if args.disk_center == 'CoM':
        center_1 = [dd['CoM'][0].in_units('cm'), dd['CoM'][1].in_units('cm'), dd['CoM'][2].in_units('cm')]
        center_2 = None
    else:
        center_1 = [dd['particle_posx'][0].in_units('cm'), dd['particle_posy'][0].in_units('cm'), dd['particle_posz'][0].in_units('cm')]
        center_2 = [dd['particle_posx'][1].in_units('cm'), dd['particle_posy'][1].in_units('cm'), dd['particle_posz'][1].in_units('cm')]
        
    time_val = ds.current_time.in_units('yr').value - sink_form

    #define cylinders:
    tube_1 = ds.disk(center_1, [0.0, 0.0, 1.0], (radius, 'au'), (height/2., 'au'))
    if center_2 != None:
        tube_2 = ds.disk(center_2, [0.0, 0.0, 1.0], (radius, 'au'), (height/2., 'au'))
        
    if center_2 != None:
        myf.set_center(1)
        #print "doing center 1, with center_vel =", dd['Center_Velocity']
        Lz_1 = np.sum(tube_1['Angular_Momentum_z'].in_units("Msun * km**2 / s").value)
        myf.set_center(2)
        del dd
        dd = ds.all_data()
        #print "doing center 2, with center_vel =", dd['Center_Velocity']
        Lz_2 = np.sum(tube_2['Angular_Momentum_z'].in_units("Msun * km**2 / s").value)
        del dd
        print "OUTPUT=", time_val, Lz_1, Lz_2, "on rank", rit
            
        #send data to rank 0 to append to write out.
        write_data = [time_val, Lz_1, Lz_2]
    else:
        myf.set_center(0)
        myf.set_coordinate_system('sph')
        v_kep = tube_1['Relative_Keplerian_Velocity']
        keplerian_inds = np.where((v_kep>0.8) & (v_kep<1.2))[0]
        disk_mass = np.sum(tube_1['cell_mass'][keplerian_inds].in_units('msun'))
        print "OUTPUT=", time_val, disk_mass
        write_data = [time_val, disk_mass]

    f = open(save_dir + output_file, 'a')
    write_string = ''
    for dit in range(len(write_data)):
        if dit != len(write_data)-1:
            write_string = write_string + str(write_data[dit]) + ','
        else:
            write_string = write_string + str(write_data[dit]) + '\n'
            print "finished creating write line"
    print "write_string =", write_string
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