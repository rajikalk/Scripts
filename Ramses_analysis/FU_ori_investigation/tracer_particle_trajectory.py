#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import glob
import numpy as np
import pickle
import gc
    

#=======MAIN=======
tracer_pickle_files = sorted(glob.glob('./movie_frame*.pkl'))

Time_array = []
X_pos = []
Y_pos = []
Z_pos = []

for tracer_file in tracer_pickle_files
    file = open(tracer_file, 'rb')
    tracer_dict = pickle.load(file)
    file.close()
    del tracer_dict['other_positions']
    gc.collect()
    del tracer_dict['burst_velocity']
    gc.collect()
    del tracer_dict['not_accreted_positions']
    gc.collect()
    
    import pdb
    pdb.set_trace()
