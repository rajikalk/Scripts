#!/usr/bin/env python
from mpi4py import MPI
import numpy as np
import h5py
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import glob
import yt
yt.enable_parallelism()
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
    parser.add_argument("-disk", "--measure_disks", help="Do you want to measure disk quantities?", type=bool, default=False)
    parser.add_argument("-lb", "--lower_bounds", help="lower bounds of box", type=float, default=50.0)
    parser.add_argument("-rad", "--radius", help="What is the radius of the measuring cylindar", type=float, default=6000.0)
    parser.add_argument("-vthres", "--velocity_threshold", help="what velocity threshold do you want to use to outflowing material?", type=float, default=0.0)
    parser.add_argument("-center", "--disk_center", help="do you want two disks centered on the particles or a volume centered on the CoM?", default='CoM', type=str)
    parser.add_argument("-st", "--start_time", help="default is -500yr", default=-500.0, type=float)
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
    if args.disk_center == 'CoM':
        radius = 100.0
        height = 100.0
        tube_center_1 = None
    else:
        radius = 20.0
        height = 20.0
        tube_center_1 = None
        tube_center_2 = None
else:
    '''
    radius = 500.0
    height = 500.0
    upper_cylinder_bound = 250.0
    lower_cylinder_bound = -250.0
    '''
    radius = args.radius
    height = 3300.0
    upper_cylinder_bound = args.lower_bounds
    lower_cylinder_bound = -1 * args.lower_bounds
    center_1 = [0.0, 0.0, (upper_cylinder_bound+(height/2.))*yt.units.AU.in_units('cm')]
    center_2 = [0.0, 0.0, (lower_cylinder_bound-(height/2.))*yt.units.AU.in_units('cm')]

#read in sink creation time
part_file = files[-1][:-12] + 'part' + files[-1][-5:]
ds = yt.load(files[-1], particle_filename=part_file)
dd = ds.all_data()
sink_form = np.min(dd['particle_creation_time']/yt.units.yr.in_units('s')).value

#find times:
if rank == 0:
    m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None, start_time=args.start_time)
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
    if args.measure_disks == False:
        f.write('Lref ' + dir.split('_')[-1] + ': Time, Mass, Momentum, Angular Momentum, Max speed, Unbound Mass, CoM dist, Mean speed, F_rad, F_tan \n')
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
        dist = np.sqrt(dd['CoM'][0]**2. + dd['CoM'][1]**2. + dd['CoM'][2]**2.).in_units('AU').value
        print "CoM dist =", dist, "(AU)"
        
        if args.measure_disks != False:
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
        
        #calculate mass in cylinders
        # for region 1
        if args.measure_disks != False:
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
                print "OUTPUT=", time_val, disk_mass, "on rank", rit
                write_data = [time_val, disk_mass]
        else:
            del dd
            myf.set_center(0)
            pos_pos = np.where(tube_1['velz'].in_units('km/s') > args.velocity_threshold)[0]
            #pos_pos = np.where(tube_1['Is_Unbound']==True)[0]
            if len(pos_pos) == 0:
                mass_1 = 0.0
                speed_1 = 0.0
                max_speed_1 = 0.0
                mean_speed_1 = 0.0
                mom_1 = 0.0
                unbound_1 = 0.0
            else:
                mass_1 = tube_1['cell_mass'][pos_pos].in_units('Msun').value
                speed_1 = tube_1['velocity_magnitude'][pos_pos].in_units("km/s").value
                max_speed_1 = np.max(speed_1)
                mean_speed_1 = np.mean(speed_1)
                mom_1 = speed_1*mass_1
                kin = (tube_1['kinetic_energy'][pos_pos] * tube_1['cell_volume'][pos_pos]).in_units('erg')
                pot = (tube_1['gpot'][pos_pos]*tube_1['cell_mass'][pos_pos]).in_units('erg')
                #eint = (tube_1['eint'][pos_pos]*tube_1['cell_mass'][pos_pos]).in_units('erg')
                tot_e = kin + pot# + eint
                un_ints = np.where(tot_e > 0.0)
                unbound_1 =tube_1['cell_mass'][pos_pos][un_ints].in_units('Msun').value
            
            #for region 2
            neg_pos = np.where(tube_2['velz'].in_units('km/s') < -1*args.velocity_threshold)[0]
            #neg_pos = np.where(tube_2['Is_Unbound']==True)[0]
            if len(neg_pos) == 0:
                mass_2 = 0.0
                speed_2 = 0.0
                max_speed_2 = 0.0
                mean_speed_2 = 0.0
                mom_2 = 0.0
                unbound_2 = 0.0
            else:
                mass_2 = tube_2['cell_mass'][neg_pos].in_units('Msun').value
                speed_2 = tube_2['velocity_magnitude'][neg_pos].in_units("km/s").value
                max_speed_2 = np.max(speed_2)
                mean_speed_2 = np.mean(speed_2)
                mom_2 = speed_2*mass_2
                kin = (tube_2['kinetic_energy'][neg_pos] * tube_2['cell_volume'][neg_pos]).in_units('erg')
                pot = (tube_2['gpot'][neg_pos]*tube_2['cell_mass'][neg_pos]).in_units('erg')
                #eint = (tube_2['eint'][neg_pos]*tube_2['cell_mass'][neg_pos]).in_units('erg')
                tot_e = kin + pot# + eint
                un_ints = np.where(tot_e > 0.0)
                unbound_2 =tube_2['cell_mass'][neg_pos][un_ints].in_units('Msun').value
            
            #Calculate outflow mass:
            outflow_mass = np.sum(mass_1) + np.sum(mass_2)
            unbound_mass = np.sum(unbound_1) + np.sum(unbound_2)
            
            #Calculate max outflow speed:
            if max_speed_1 > max_speed_2:
                max_speed = abs(max_speed_1)
            else:
                max_speed = abs(max_speed_2)
            mean_speed = (mean_speed_1 + mean_speed_2)/2.

            #Calculate outflow momentum:
            mom = np.sum(mom_1) + np.sum(mom_2)

            #Calculate outflow angular momentum:
            L_x = np.sum(tube_1['Angular_Momentum_x'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['Angular_Momentum_x'][neg_pos].in_units("Msun * km**2 / s").value)
            L_y = np.sum(tube_1['Angular_Momentum_y'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['Angular_Momentum_y'][neg_pos].in_units("Msun * km**2 / s").value)
            L_z = np.sum(tube_1['Angular_Momentum_z'][pos_pos].in_units("Msun * km**2 / s").value) + np.sum(tube_2['Angular_Momentum_z'][neg_pos].in_units("Msun * km**2 / s").value)
            L = L_z
            #L = np.sqrt((L_x)**2. + (L_y)**2. + (L_z)**2.)
            
            #Append quantities
            if outflow_mass == 0.0:
                outflow_mass = np.nan
                max_speed = np.nan
                mom = np.nan
                L = np.nan
                unbound_mass = np.nan
        
            F_rad = np.sum(dd['Gravitational_Force_on_particles_Rad'].value)
            F_tan = np.sum(dd['Gravitational_Force_on_particles_Tan'].value)
            
            print "OUTPUT=", time_val, outflow_mass, mom, L, max_speed, unbound_mass, dist, mean_speed, F_rad, F_tan, "on rank", rit

            #send data to rank 0 to append to write out.
            write_data = [time_val, outflow_mass, mom, L, max_speed, unbound_mass, dist, mean_speed, F_rad, F_tan]
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

    if len(usable_files) < size:
        if rank >= size:
            break

if rank == 0:
    print "GETTING THE FINAL STUFF TO WRITE OUT"
    for other_rank in range(1,size):
        write_data = comm.recv(source=other_rank, tag=2)
        print "received data", write_data, "from rank", other_rank
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

print "Completed job on rank", rank
