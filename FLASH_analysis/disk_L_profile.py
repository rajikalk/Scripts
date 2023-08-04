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
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="z")
    parser.add_argument("-make_pickles", "--make_movie_pickles", type=str, default='True')
    parser.add_argument("-make_frames", "--make_movie_frames", type=str, default='True')
    
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    parser.add_argument("-dt", "--time_step", help="time step between movie frames", default = 10., type=float)
    
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
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None)
    print('generated frame times')
    
    Time_array = []
    Radius_array = []
    All_profiles_array = []

    if rank == 0:
        pickle_names = 'profile_*.pkl'
        pickle_files = glob.glob(pickle_names)
        
        if len(pickle_files) > 0:
            for pickle_file in pickle_files:
                file = open(pickle_file, 'rb')
                time_val, radius, All_profiles = pickle.load(file)
                file.close()
                
                Time_array = Time_array + time_val
                Radius_array = Radius_array + radius
                All_profiles_array = All_profiles_array + All_profiles
                
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            Radius_array = list(np.array(Radius_array)[sorted_inds])
            All_profiles_array = list(np.array(All_profiles_array)[sorted_inds])
            
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
    start_frame = m_times.index(start_time)
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
            if args.plot_time == None:
                time_val = m_times[file_int]#ds.current_time.in_units('yr')
            else:
                dd = ds.all_data()
                time_val = int(yt.YTQuantity(ds.current_time.value - np.min(dd['particle_creation_time']).value, 's').in_units('yr').value)
            
            dd = ds.all_data()

            #Define cylinder!:
            primary_ind = np.argmin(dd['particle_creation_time'])
            center = yt.YTArray([dd['particle_posx'][primary_ind], dd['particle_posy'][primary_ind], dd['particle_posz'][primary_ind]])
            normal = yt.YTArray([0, 0, 1], '')
            height = yt.YTQuantity(50, 'au')
            if len(dd['particle_creation_time']) == 1:
                radius = yt.YTQuantity(100, 'au')
            else:
                part_pos = yt.YTArray([dd['particle_posx'], dd['particle_posy'], dd['particle_posz']])
                nearest_sink_ind = np.argsort(np.sqrt(np.sum((part_pos.T - center)**2, axis=1)).in_units('au'))[1]
                nearest_sink_pos = part_pos.T[nearest_sink_ind]
                radius = np.sqrt(np.sum((center - nearest_sink_pos)**2)).in_units('au')
            disk = ds.disk(center, normal, radius, height)
            Radius_field = disk['radius'].in_units('AU')
            L_disk = disk['L_gas_wrt_primary']
            r_bins = np.arange(0, radius.value+5, 5)
            r_centers = []
            L_means = []
            L_means_spec = []
            for bit in range(1,len(r_bins[1:])):
                usable_inds = np.where((Radius_field>r_bins[bit-1])&(Radius_field<r_bins[bit]))
                weighted_mean = np.sum((disk['L_gas_wrt_primary'][usable_inds]*disk['mass'][usable_inds])/np.sum(disk['mass'][usable_inds])
                means_spec = np.sum((disk['L_gas_wrt_primary'][usable_inds]/disk['mass'][usable_inds])*disk['mass'][usable_inds])/np.sum(disk['mass'][usable_inds])
                r_centers.append(np.mean(r_bins[bit-1:bit+1]))
                L_means.append(weighted_mean)
                L_means_spec.append(means_spec)
            
            Time_array.append(time_val)
            Radius_array.append(radius)
            All_profiles_array.append([r_centers, L_means, L_means_spec])

            pickle_file = 'profile_'+str(rank)+'.pkl'
            file = open(pickle_file, 'wb')
            pickle.dump((Time_array, Radius_array, All_profiles_array), file)
            file.close()
            print("Calculated angular momentum profile on", rank, "for file", file_int, "of ", no_frames)
    
sys.stdout.flush()
CW.Barrier()

#collect pickles
if rank == 0:
    pickle_names = 'profile_*.pkl'
    pickle_files = glob.glob(pickle_names)
    
    if len(pickle_files) > 0:
        for pickle_file in pickle_files:
            file = open(pickle_file, 'rb')
            time_val, radius, All_profiles = pickle.load(file)
            file.close()
            
            Time_array = Time_array + time_val
            Radius_array = Radius_array + radius
            All_profiles_array = All_profiles_array + All_profiles
            
        sorted_inds = np.argsort(Time_array)
        Time_array = list(np.array(Time_array)[sorted_inds])
        Radius_array = list(np.array(Radius_array)[sorted_inds])
        All_profiles_array = np.array(All_profiles_array)[sorted_inds]
        
    file = open('gathered_profile.pkl', 'wb')
    pickle.dump((Time_array, Radius_array, All_profiles_array), file)
    file.close()
    
sys.stdout.flush()
CW.Barrier()

file = open('gathered_profile.pkl', 'rb')
Time_array, Radius_array, All_profiles_array = pickle.load(file)
file.close()

max_val = 0
min_val = np.inf
for profile in All_profiles_array:
    if np.max(profile[1]) > max_val:
        max_val = np.max(profile[1])
    if np.min(profile[1]) < min_val:
        min_val = np.min(profile[1])

sys.stdout.flush()
CW.Barrier()

import matplotlib.pyplot as plt

rit = -1
for time_it in range(len(Time_array)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        inner_inds = np.where(np.array(All_profiles_array[time_it][0])<20)[0]
        mean_L = np.mean(np.array(All_profiles_array[time_it][1])[inner_inds])
        plt.clf()
        plt.semilogy(All_profiles_array[time_it][0], All_profiles_array[time_it][1])
        plt.axvline(x = Radius_array[time_it], color='k')
        plt.axhline(y = mean_L, xmin=0, xmax=0.2, color='r')
        plt.xlim([0,100])
        plt.ylim([min_val, max_val])
        plt.xlabel('Radius (AU)')
        plt.ylabel('L ($g\,cm^2/s$)')
        plt.title('Time:'+str(Time_array[time_it])+'yr')
        
        save_name = 'movie_frame_' + ("%06d" % time_it)
        plt.savefig(save_name+'.jpg', dpi=300, bbox_inches='tight')
        print('Made frame', time_it, 'of', len(Time_array), 'on rank', rank)
