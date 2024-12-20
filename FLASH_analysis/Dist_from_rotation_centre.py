#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
import pickle5 as pickle
import os
import my_flash_module as mym
import my_flash_fields as myf

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#------------------------------------------------------
def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-make_pickles", "--make_movie_pickles", type=str, default='True')
    parser.add_argument("-radius", "--radius_threshold", type=float, default=100)
    
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    parser.add_argument("-start", "--start_time", type=float, default=0)
    parser.add_argument("-end", "--end_time", type=float, default=None)
    
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#-------------------------------------------------------
#get input and output directory and arguments
input_dir = sys.argv[1]
output_dir = sys.argv[2]
args = parse_inputs()
center_pos = [0, 0, 0]

if args.make_movie_pickles == 'True':
    files = sorted(glob.glob(input_dir + '*plt_cnt*'))
    m_times = mym.generate_frame_times(files, args.time_step, start_time=args.start_time, presink_frames=0, end_time=args.end_time)
    print('generated frame times')
    
    Time_array = []
    Distance = []

    if rank == 0:
        pickle_names = 'profile_*.pkl'
        pickle_files = glob.glob(pickle_names)
        
        if len(pickle_files) > 0:
            for pickle_file in pickle_files:
                file = open(pickle_file, 'rb')
                time_val, dist = pickle.load(file)
                file.close()
                
                Time_array = Time_array + time_val
                Distance = Distance + dist
                
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            Distance = list(np.array(Distance)[sorted_inds])
            
            start_time = np.max(Time_array)
        else:
            start_time = m_times[0]
    else:
        start_time = np.nan

    sys.stdout.flush()
    CW.Barrier()

    start_time = CW.bcast(start_time, root=0)

    sys.stdout.flush()
    CW.Barrier()
    
    no_frames = len(m_times)
    try:
        start_frame = m_times.index(start_time)
    except:
        start_frame = np.argmin(abs(np.array(m_times) - start_time))
        
    m_times = m_times[start_frame:]
    usable_files = mym.find_files(m_times, files)
    frames = list(range(start_frame, no_frames))

    sys.stdout.flush()
    CW.Barrier()
    
    usable_files = mym.find_files(m_times, files)
    frames = list(range(len(usable_files)))
    no_frames = len(usable_files)
    print('found usable files for frames')

    #Now let's iterate over the files and get the images we want to plot
    ts = yt.DatasetSeries(usable_files)
    file_int = -1
    rit = -1
    for ds in ts:
        file_int = file_int + 1
        #print(fn, "is going to rank", rank)
        rit = rit + 1
        if rit == size:
            rit = 0
        if rit == rank:

            dd = ds.all_data()
            time_val = int(yt.YTQuantity(ds.current_time.value - np.min(dd['particle_creation_time']).value, 's').in_units('yr').value)
            
            Time_array.append(time_val)
            dd = ds.all_data()

            #Define cylinder!:
            try:
                primary_ind = np.argmin(dd['particle_creation_time'])
            except:
                curr_file = ds.filename
                next_file = curr_file[:-4] + ("%04d"%(int(curr_file[-4:])+1))
                ds = yt.load(next_file)
                dd = ds.all_data()
                primary_ind = np.argmin(dd['particle_creation_time'])
                
            center = yt.YTArray([dd['particle_posx'][primary_ind], dd['particle_posy'][primary_ind], dd['particle_posz'][primary_ind]]).in_units('au')
            normal = yt.YTArray([0, 0, 1], '')
            height = yt.YTQuantity(args.radius_threshold, 'au')
            radius = yt.YTQuantity(args.radius_threshold, 'au')
            disk = ds.disk(center, normal, radius, height)
            com_x = (np.sum(disk['x'].in_units('au')*disk['mass'].in_units('g')))/(np.sum(disk['mass'].in_units('g')))
            com_y = (np.sum(disk['y'].in_units('au')*disk['mass'].in_units('g')))/(np.sum(disk['mass'].in_units('g')))
            com_z = (np.sum(disk['z'].in_units('au')*disk['mass'].in_units('g')))/(np.sum(disk['mass'].in_units('g')))
            com = yt.YTArray([com_x, com_y, com_z])
            
            com_x = (np.sum(disk['velx'].in_units('cm/s')*disk['mass'].in_units('g')))/(np.sum(disk['mass'].in_units('g')))
            com_y = (np.sum(disk['vely'].in_units('cm/s')*disk['mass'].in_units('g')))/(np.sum(disk['mass'].in_units('g')))
            com_z = (np.sum(disk['velz'].in_units('cm/s')*disk['mass'].in_units('g')))/(np.sum(disk['mass'].in_units('g')))
            com_v = yt.YTArray([com_x, com_y, com_z])
            
            dist = np.sqrt(np.sum((center - com)**2))
            Distance.append(dist)

            pickle_file = 'profile_'+str(rank)+'.pkl'
            file = open(pickle_file, 'wb')
            pickle.dump((Time_array, Distance), file)
            file.close()
            print("Calculated angular momentum profile on", rank, "for file", file_int, "of ", no_frames)
    
sys.stdout.flush()
CW.Barrier()

#collect pickles
if rank == 0:
    pickle_names = 'profile_*.pkl'
    pickle_files = glob.glob(pickle_names)
    
    Time_array = []
    Distance = []

    if rank == 0:
        pickle_names = 'profile_*.pkl'
        pickle_files = glob.glob(pickle_names)
        
        if len(pickle_files) > 0:
            for pickle_file in pickle_files:
                file = open(pickle_file, 'rb')
                time_val, dist = pickle.load(file)
                file.close()
                
                Time_array = Time_array + time_val
                Distance = Distance + dist
                
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            Distance = list(np.array(Distance)[sorted_inds])
        
    file = open('gathered_profile.pkl', 'wb')
    pickle.dump((Time_array, Distance), file)
    file.close()
    
sys.stdout.flush()
CW.Barrier()

file = open('gathered_profile.pkl', 'wb')
pickle.dump((Time_array, Distance), file)
file.close()
