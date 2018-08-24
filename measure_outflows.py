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
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", help="What is the input directory?", required=True)
    parser.add_argument("-o", "--output_dir", help="what is the output directory? if this is not defined it will be the same as the input directory")
    parser.add_argument("-of", "--output_file", help="Do you have a file ot write out to?", default="out.csv")
    parser.add_argument("-ff", "--file_frequency", help="every how many files do you want to use?", default=100, type=int)
    parser.add_argument("-dt", "--time_step", help="How frequently do you want to measure the outflow?", default=2, type=int)
    parser.add_argument("-disk", "--measure_disks", help="Do you want to measure disk quantities?", type=str, default="False")
    parser.add_argument("-lb", "--lower_bounds", help="lower bounds of box", type=float, default=50.0)
    parser.add_argument("-rad", "--radius", help="What is the radius of the measuring cylindar", type=float, default=6000.0)
    parser.add_argument("-vthres", "--velocity_threshold", help="what velocity threshold do you want to use to outflowing material?", type=float, default=0.0)
    parser.add_argument("-center", "--disk_center", help="do you want two disks centered on the particles or a volume centered on the CoM?", default='CoM', type=str)
    parser.add_argument("-st", "--start_time", help="default is -500yr", default=-500.0, type=float)
    parser.add_argument("-m_all", "--measure_all", help="do you want to measure everything tha that is inflowing/outflowing/disk?", default="False")
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
if args.measure_disks != "False":
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
m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None, start_time=args.start_time)
usable_files = mym.find_files(m_times, files)

#open files to read data out to:
if rank == 0:
    if os.path.isfile(save_dir + output_file) == False:
        f = open(save_dir + output_file, 'a+')
        f.close()
    f = open(save_dir + output_file, 'w')
    if args.measure_disks != "False":
        f.write('Lref ' + dir.split('_')[-1] + ': Time, Lz_1, Lz_2 \n')
    elif args.measure_all != "False":
        f.write('Time, Particles, Disk, Infalling, Outflowing \n')
    else:
        f.write('Lref ' + dir.split('_')[-1] + ': Time, Mass, Momentum, Angular Momentum, Max speed, Unbound Mass, CoM dist, Mean speed, F_rad, F_tan \n')
    f.close()

storage = {}
ts = yt.DatasetSeries(usable_files, parallel=True)
for sto, ds in ts.piter(storage=storage):
    if args.measure_disks != "False":
        dd = ds.all_data()
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
    if args.measure_disks != "False":
        if center_2 != None:
            myf.set_center(1)
            #print "doing center 1, with center_vel =", dd['Center_Velocity']
            Lz_1 = np.sum(tube_1['Angular_Momentum_z'].in_units("Msun * km**2 / s").value)
            myf.set_center(2)
            #print "doing center 2, with center_vel =", dd['Center_Velocity']
            Lz_2 = np.sum(tube_2['Angular_Momentum_z'].in_units("Msun * km**2 / s").value)
            
            #send data to rank 0 to append to write out.
            write_data = [time_val, Lz_1, Lz_2]
        else:
            myf.set_center(0)
            myf.set_coordinate_system('sph')
            v_kep = tube_1['Relative_Keplerian_Velocity']
            keplerian_inds = np.where((v_kep>0.8) & (v_kep<1.2))[0]
            disk_mass = np.sum(tube_1['cell_mass'][keplerian_inds].in_units('msun'))
            write_data = [time_val, disk_mass]
    elif args.measure_all != "False":
        disk_volume = ds.disk([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], (radius, 'au'), (upper_cylinder_bound, 'au'))
        ang_disk = np.sum(disk_volume['Angular_Momentum'].value)
        if ('all', u'particle_mass') in ds.field_list:
            ang_part = np.sum(disk_volume['Particle_Angular_Momentum'].value)
        else:
            ang_part = 0.0
        myf.set_center(0)
        pos_pos = np.where(tube_1['velz'].in_units('km/s') > args.velocity_threshold)[0]
        other_pos_pos = np.where(tube_1['velz'].in_units('km/s') < args.velocity_threshold)[0]
        ang_in_1 = np.sum(tube_1['Angular_Momentum'][other_pos_pos].value)
        if len(pos_pos) > 0:
            ang_out_1 = np.sum(tube_1['Angular_Momentum'][pos_pos].value)
        else:
            ang_out_1 = 0.0

        neg_pos = np.where(tube_2['velz'].in_units('km/s') < -1*args.velocity_threshold)[0]
        other_neg_pos = np.where(tube_2['velz'].in_units('km/s') > -1*args.velocity_threshold)[0]
        ang_in_2 = np.sum(tube_2['Angular_Momentum'][other_neg_pos].value)
        if len(pos_pos) > 0:
            ang_out_2 = np.sum(tube_2['Angular_Momentum'][neg_pos].value)
        else:
            ang_out_2 = 0.0
                
        ang_in = ang_in_1 + ang_in_2
        ang_out = ang_out_1 + ang_out_2
        write_data = [time_val, ang_part, ang_disk, ang_in, ang_out]
        print "OUTPUT=", time_val, ang_part, ang_disk, ang_in, ang_out
    else:
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


        print "OUTPUT=", time_val, outflow_mass, mom, L, max_speed, unbound_mass, dist, mean_speed
            
        #send data to rank 0 to append to write out.
        write_data = [time_val, outflow_mass, mom, L, max_speed, unbound_mass, dist, mean_speed]

    write_string = ''
    for dit in range(len(write_data)):
        if dit != len(write_data)-1:
            write_string = write_string + str(write_data[dit]) + ','
        else:
            write_string = write_string + str(write_data[dit]) + '\n'

    sto.result = write_string
    sto.result_id = str(time_val)
    del write_data
    del disk_volume
    del tube_1
    del tube_2
    del pos_pos
    del other_pos_pos
    del neg_pos
    del other_neg_pos

pickle_file = save_dir + "AMB_" + str(float(dir.split('Mach_')[-1].split('/')[0])) + ".pkl"
file = open(pickle_file, 'w+')
pickle.dump((storage), file)
file.close()

if rank == 0:
    del_keys = []
    for key, value in storage.iteritems():
        if value is None:
            del_keys.append(key)
    for key in del_keys:
        del storage[key]
    f = open(save_dir + output_file, 'a')
    for it in sorted(np.array(storage.keys()).astype(np.float)):
        f.write(storage[str(it)])
        print "Printed line:", storage[str(it)]
    f.close()
print "Completed job on rank", rank
