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

    if rank == 0:
        pickle_names = 'profile_*.pkl'
        pickle_files = glob.glob(pickle_names)
        
        Time_array = []
        Radius_array = []
        All_profiles_array = []
        
        if len(pickle_files) > 0:
            for pickle_file in pickle_files:
                file = open(pickle_file, 'rb')
                time_val, radius, All_profiles = pickle.load(file)
                file.close()
                
                Time_array.append(time_val)
                Radius_array.append(radius)
                All_profiles_array = All_profiles_array + All_profiles
                
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            import pdb
            pdb.set_trace()
            
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
    file_int = -1
    for fn in yt.parallel_objects(usable_files, njobs=size):
        if size > 1:
            file_int = usable_files.index(fn)
        else:
            file_int = file_int + 1
            if usable_files[file_int] == usable_files[file_int-1]:
                os.system('cp '+ output_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + output_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        file_counter = usable_files.index(fn)
        pickle_file = output_dir+"movie_frame_"+("%06d" % file_counter)+".pkl"
        make_pickle = True
        #if os.path.isfile(pickle_file) == True:
        #if len(glob.glob(pickle_file)) == 1:
        #    make_pickle = False
        if make_pickle:
            #print(fn, "is going to rank", rank)
            part_file = 'part'.join(fn.split('plt_cnt'))
            ds = yt.load(fn, particle_filename=part_file)
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
            for rit in range(1,len(r_bins[1:])):
                usable_inds = np.where((Radius_field>r_bins[rit-1])&(Radius_field<r_bins[rit]))
                weighted_mean = np.sum(disk['L_gas_wrt_primary'][usable_inds]*disk['mass'][usable_inds])/np.sum(disk['mass'][usable_inds])
                r_centers.append(np.mean(r_bins[rit-1:rit+1]))
                L_means.append(weighted_mean)
            
            All_profiles.append([r_centers, L_means])

            pickle_file = 'profile_'+str(rank)+'.pkl'
            file = open(pickle_file, 'wb')
            pickle.dump((time_val, radius, All_profiles), file)
            file.close()
            print("Calculated angular momentum profile on", rank)
    
#collect pickles

