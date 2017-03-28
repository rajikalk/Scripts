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

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="What is the input directory?", required=True)
    parser.add_argument("-o", "--output_dir", help="what is the output directory? if this is not defined it will be the same as the input directory")
    parser.add_argument("-of", "--output_file", help="Do you have a file ot write out to?", default="out.csv")
    parser.add_argument("-ff", "--file_frequency", help="every how many files do you want to use?", default=100, type=int)
    parser.add_argument("-dt", "--time_step", help="How frequently do you want to measure the outflow?", default=2, type=int)
    #parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

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

#define measured quantities:
time = []
mass = []
maximum_speed = []
momentum = []
ang_momentum = []
#define analysis cylinder:
radius = 500.0
#region:
upper_cylinder_bound = 250.0
lower_cylinder_bound = -250.0
height = 500.0

#read in sink creation time
part_file = files[-1][:-12] + 'part' + files[-1][-5:]
ds = yt.load(files[-1], particle_filename=part_file)
dd = ds.all_data()
sink_form = np.min(dd['particle_creation_time']/yt.units.yr.in_units('s')).value

#find times:
m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None)
usable_files = mym.find_files(m_times, files)

#open files to read data out to:
if os.path.isfile(save_dir + output_file) == False:
    f = open(save_dir + output_file, 'a+')
    f.close()
f = open(save_dir + output_file, 'w')
f.write('Lref ' + dir.split('_')[-1] + ': Time, Mass, Momentum, Angular Momentum, Max speed \n')
f.close()

file_no = 0
for file in usable_files:
    part_file = file[:-12] + 'part' + file[-5:]
    ds = yt.load(file, particle_filename=part_file)
    dd = ds.all_data()
    
    time_val = ds.current_time.in_units('yr').value - sink_form

    #define cylinders:
    tube_1 = ds.disk([0.0, 0.0, (upper_cylinder_bound+(height/2.))*yt.units.AU.in_units('cm')], [0.0, 0.0, 1.0], (radius, 'au'), (height/2., 'au'))
    tube_2 = ds.disk([0.0, 0.0, (lower_cylinder_bound-(height/2.))*yt.units.AU.in_units('cm')], [0.0, 0.0, -1.0], (radius, 'au'), (height/2., 'au'))

    #calculate mass in cylinders
    # for region 1
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
        
    time.append(time_val)
    mass.append(outflow_mass)
    maximum_speed.append(max_speed)
    momentum.append(mom)
    ang_momentum.append(L)
    print "OUTPUT=", time_val, outflow_mass, mom, L, max_speed
    f = open(save_dir + output_file, 'a')
    f.write(str(time_val) + ',' + str(outflow_mass) + ',' + str(mom) + ',' + str(L) + ',' + str(max_speed) + '\n')
    f.close()
    prev_time = time_val

    file_no = file_no + 1
    print "Progress:", file_no, "/", len(usable_files)