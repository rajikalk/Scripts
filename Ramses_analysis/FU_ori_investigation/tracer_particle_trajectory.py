#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import glob
import pickle
import gc
import os
import numpy as np
from my_ramses_fields import projected_vector
    

#=======MAIN=======

if os.path.isfile('tracer_trajectory.pkl') == False:
    tracer_pickle_files = sorted(glob.glob('./movie_frame*.pkl'))

    for tracer_file in tracer_pickle_files:
        print("reading", tracer_file)
        file = open(tracer_file, 'rb')
        tracer_dict = pickle.load(file)
        file.close()
        
        import pdb
        pdb.set_trace()

    #save the data
    #print("saving tracer trajectories")
    #file = open('tracer_trajectory.pkl', 'wb')
    #pickle.dump((Time_array, X_pos, Y_pos, Z_pos, Sink_vec), file)
    #file.close()

#else:
#    print("reading tracer trajectories")
#    file = open('tracer_trajectory.pkl', 'rb')
#    Time_array, X_pos, Y_pos, Z_pos, Sink_vec = pickle.load(file)
#    file.close()
    

#Make plots!
import matplotlib.pyplot as plt
X_pos = np.array(X_pos)
Y_pos = np.array(Y_pos)
alpha_val = np.linspace(0, 1, len(Time_array))

#for time_it in range(len(Time_array)):
    #PLot sink vector

plt.clf()
for part_it in range(np.shape(X_pos)[1]):
    plt.plot(X_pos.T[part_it], Y_pos.T[part_it], alpha=alpha_val)
plt.scatter(0, 0, marker='o', color='magenta', s=3)
plt.xlim([-15, 15])
plt.ylim([-15, 15])
plt.xlabel('X (AU)')
plt.xlabel('Y (AU)')
plt.savefig("XY_tracer_traj.png")
