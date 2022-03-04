import pickle
import glob
import sys
import os
from mpi4py.MPI import COMM_WORLD as CW

rank = CW.Get_rank()
size = CW.Get_size()

input_dir = sys.argv[1]
output_dir = sys.argv[2]

file_dirs = sorted(glob.glob(input_dir+"movie_frame_*"))

rit = -1
for file_dir in file_dirs:
    for proj_it in range(8):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
            try:
                pickle_file = file_dir + '/projection_'+str(proj_it)+'.pkl'
                file = open(pickle_file, 'rb')
                X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
                file.close()
                rv_image = vel_rad/100000
                
                save_dir = output_dir + file_dir.split('/')[-1]
                if os.path.exists(save_dir) == False:
                    os.makedirs(save_dir)
                save_pickle = save_dir + '/projection_'+str(proj_it)+'.pkl'
                save_rv_pickle = save_dir + '/projection_'+str(proj_it)+'_rv.pkl'
                file = open(save_pickle, 'wb')
                pickle.dump((X, Y, image.value, ), file)
                file.close()
                file = open(save_rv_pickle, 'wb')
                pickle.dump((X, Y, rv_image.value), file)
                file.close()
                print("saved", save_pickle)
            except:
                print("skipping because", file_dir + '/projection_'+str(proj_it)+'.pkl', "doesn't exist")
