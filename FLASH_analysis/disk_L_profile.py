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
        m_times = mym.generate_frame_times(files, args.time_step, presink_frames=0, end_time=None, )
    print('generated frame times')
    sys.stdout.flush()
    CW.Barrier()
    
    usable_files = mym.find_files(m_times, files)
    frames = list(range(len(usable_files)))
    no_frames = len(usable_files)
    print('found usable files for frames')

    #Now let's iterate over the files and get the images we want to plot
    file_int = -1
    if size > 1:
        njobs = int(size/5)
    else:
        njobs = 1
    for fn in yt.parallel_objects(usable_files, njobs=njobs):
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
            
                        #Get particle data:
            dd = ds.all_data()

            #Define cylinder!:
            import pdb
            pdb.set_trace()



    print("Calculated angular momentum profile on", rank)
