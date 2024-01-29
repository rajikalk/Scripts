#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import sys
import argparse
from mpi4py.MPI import COMM_WORLD as CW
import numpy as np
try:
    import pickle5 as pickle
except:
    import pickle
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
    parser.add_argument("-field", "--profile_field", type=str, default='L_gas_wrt_primary_spec')
    parser.add_argument("-inner_radius", "--inner_radius_threshold", type=float, default=20)
    parser.add_argument("-height", "--disk_height", type=float, default=20)
    parser.add_argument("-weight", "--weight_field", type=str, default=None)
    
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
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
    if args.plot_time != None:
        m_times = [args.plot_time]
    else:
        m_times = mym.generate_frame_times(files, args.time_step, start_time=args.start_time, presink_frames=0, end_time=args.end_time)
    print('generated frame times')
    
    R_sink = np.nan
    Time_array = []
    Total_L = []
    Total_L_spec = []
    Mean_L = []
    Mean_L_spec = []
    Mean_rad_vel = []
    Min_rad_vel = []
    Mean_rel_vel = []
    Min_rel_vel = []
    Mass_all = []
    Separation = []

    if rank == 0:
        pickle_names = 'profile_*.pkl'
        pickle_files = glob.glob(pickle_names)
        
        if len(pickle_files) > 0:
            for pickle_file in pickle_files:
                file = open(pickle_file, 'rb')
                R_sink, time_val, L_tot, L_spec_tot, L_mean, L_spec_mean, Rad_vel, Rad_vel_min, Mean_rel, Min_rel, Mass, sep = pickle.load(file)
                file.close()
                
                Time_array = Time_array + time_val
                Total_L = Total_L + L_tot
                Total_L_spec = Total_L_spec + L_spec_tot
                Mean_L = Mean_L + L_mean
                Mean_L_spec = Mean_L_spec + L_spec_mean
                Mean_rad_vel = Mean_rad_vel + Rad_vel
                Min_rad_vel = Min_rad_vel + Rad_vel_min
                Mean_rel_vel = Mean_rel_vel + Mean_rel
                Min_rel_vel = Min_rel_vel + Min_rel
                Mass_all = Mass_all + Mass
                Separation = Separation + sep
                
            sorted_inds = np.argsort(Time_array)
            Time_array = list(np.array(Time_array)[sorted_inds])
            Total_L = list(np.array(Total_L)[sorted_inds])
            Total_L_spec = list(np.array(Total_L_spec)[sorted_inds])
            Mean_L = list(np.array(Mean_L)[sorted_inds])
            Mean_L_spec = list(np.array(Mean_L_spec)[sorted_inds])
            Mean_rad_vel = list(np.array(Mean_rad_vel)[sorted_inds])
            Min_rad_vel = list(np.array(Min_rad_vel)[sorted_inds])
            Mean_rel_vel = list(np.array(Mean_rel_vel)[sorted_inds])
            Min_rel_vel = list(np.array(Min_rel_vel)[sorted_inds])
            Mass_all = list(np.array(Mass_all)[sorted_inds])
            Separation = list(np.array(Separation)[sorted_inds])
            
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
                
            center = yt.YTArray([dd['particle_posx'][primary_ind], dd['particle_posy'][primary_ind], dd['particle_posz'][primary_ind]])
            if len(dd['particle_creation_time']) > 1:
                part_pos = yt.YTArray([dd['particle_posx'], dd['particle_posy'], dd['particle_posz']])
                nearest_sink_ind = np.argsort(np.sqrt(np.sum((part_pos.T - center)**2, axis=1)).in_units('au'))[1]
                nearest_sink_pos = part_pos.T[nearest_sink_ind]
                sep = np.sqrt(np.sum((center - nearest_sink_pos)**2)).in_units('au')
            else:
                sep = yt.YTQuantity(np.nan, 'au')
            
            Separation.append(sep)
                
            #normal = yt.YTArray([0, 0, 1], '')
            L_primary = yt.YTArray([dd['particle_x_ang'][primary_ind], dd['particle_y_ang'][primary_ind], dd['particle_z_ang'][primary_ind]])
            normal = L_primary/np.sqrt(np.sum(L_primary**2))
            height = yt.YTQuantity(args.disk_height, 'au')
            radius = yt.YTQuantity(args.inner_radius_threshold, 'au')
            disk = ds.disk(center, normal, radius, height)
            R_sink = 2.5*np.min(disk['dx']).in_units('au')
            
            Radius_field = disk['radius'].in_units('AU')
            Total_L.append(np.sum(disk[args.profile_field]))
            if args.weight_field == None:
                Mean_L.append(np.mean(disk[args.profile_field]))
            else:
                weighted_mean = np.sum(disk[args.profile_field]*disk[args.weight_field])/np.sum(disk[args.weight_field])
                Mean_L.append(weighted_mean)
            #spec_field = args.profile_field.split('_cyl')[0] + '_spec'
            Total_L_spec.append(np.sum(disk['mass']))
            Mean_L_spec.append(np.mean(disk['mass']))
            Mean_rad_vel.append(np.mean(disk[args.profile_field]))
            Min_rad_vel.append(np.min(disk[args.profile_field]))
            Mean_rel_vel.append(np.mean(disk[args.profile_field]))
            Min_rel_vel.append(np.min(disk[args.profile_field]))
            Mass_all.append(np.sum(disk['mass'].in_units('msun')))

            pickle_file = 'profile_'+str(rank)+'.pkl'
            file = open(pickle_file, 'wb')
            pickle.dump((R_sink, Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mean_rad_vel, Min_rad_vel, Mean_rel_vel, Min_rel_vel, Mass_all, Separation), file)
            file.close()
            print("Calculated angular momentum profile on", rank, "for file", file_int, "of ", no_frames)
    
sys.stdout.flush()
CW.Barrier()

#collect pickles
if rank == 0:
    pickle_names = 'profile_*.pkl'
    pickle_files = glob.glob(pickle_names)
    
    R_sink = np.nan
    Time_array = []
    Total_L = []
    Total_L_spec = []
    Mean_L = []
    Mean_L_spec = []
    Mean_rad_vel = []
    Min_rad_vel = []
    Mean_rel_vel = []
    Min_rel_vel = []
    Mass_all = []
    Separation = []
    
    if len(pickle_files) > 0:
        for pickle_file in pickle_files:
            file = open(pickle_file, 'rb')
            R_sink, time_val, L_tot, L_spec_tot, L_mean, L_spec_mean, Rad_vel, Rad_vel_min, Mean_rel, Min_rel, Mass, sep = pickle.load(file)
            file.close()
            
            Time_array = Time_array + time_val
            Total_L = Total_L + L_tot
            Total_L_spec = Total_L_spec + L_spec_tot
            Mean_L = Mean_L + L_mean
            Mean_L_spec = Mean_L_spec + L_spec_mean
            Mean_rad_vel = Mean_rad_vel + Rad_vel
            Min_rad_vel = Min_rad_vel + Rad_vel_min
            Mean_rel_vel = Mean_rel_vel + Mean_rel
            Min_rel_vel = Min_rel_vel + Min_rel
            Mass_all = Mass_all + Mass
            Separation = Separation + sep
            
        sorted_inds = np.argsort(Time_array)
        Time_array = list(np.array(Time_array)[sorted_inds])
        Total_L = list(np.array(Total_L)[sorted_inds])
        Total_L_spec = list(np.array(Total_L_spec)[sorted_inds])
        Mean_L = list(np.array(Mean_L)[sorted_inds])
        Mean_L_spec = list(np.array(Mean_L_spec)[sorted_inds])
        Mean_rad_vel = list(np.array(Mean_rad_vel)[sorted_inds])
        Min_rad_vel = list(np.array(Min_rad_vel)[sorted_inds])
        Mean_rel_vel = list(np.array(Mean_rel_vel)[sorted_inds])
        Min_rel_vel = list(np.array(Min_rel_vel)[sorted_inds])
        Mass_all = list(np.array(Mass_all)[sorted_inds])
        Separation = list(np.array(Separation)[sorted_inds])
        
    file = open('gathered_profile.pkl', 'wb')
    pickle.dump((R_sink, Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mean_rad_vel, Min_rad_vel, Mean_rel_vel, Min_rel_vel, Mass_all, Separation), file)
    file.close()
    
sys.stdout.flush()
CW.Barrier()

file = open('gathered_profile.pkl', 'wb')
pickle.dump((R_sink, Time_array, Total_L, Total_L_spec, Mean_L, Mean_L_spec, Mean_rad_vel, Min_rad_vel,  Mean_rel_vel, Min_rel_vel, Mass_all, Separation), file)
file.close()
