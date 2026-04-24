#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import glob
import pickle
import gc
import os
import numpy as np
from my_ramses_fields import projected_vector
import yt
    

#=======MAIN=======

if os.path.isfile('tracer_trajectory.pkl') == False:
    tracer_pickle_files = sorted(glob.glob('./movie_frame*.pkl'))

    Time_array = []
    Tracer_parallel = []
    Tracer_perpendicular = []
    Sink_vec = []

    for tracer_file in tracer_pickle_files:
        print("reading", tracer_file)
        file = open(tracer_file, 'rb')
        tracer_dict = pickle.load(file)
        file.close()
        del tracer_dict['other_positions']
        gc.collect()
        del tracer_dict['not_accreted_positions']
        gc.collect()
        del tracer_dict['burst_velocity']
        gc.collect()
        
        Tracer_distance = np.sqrt(np.sum(np.square(tracer_dict['burst_positions']), axis=0))
        Tracer_proj_parallel = projected_vector(yt.YTArray(np.array(tracer_dict['burst_positions']).T, 'au'), tracer_dict['sink_velocity_vector'])
        sign = np.sign(np.dot(Tracer_proj_parallel, tracer_dict['sink_velocity_vector']))
        Tracer_proj_parallel_mag = sign*np.sqrt(np.sum(Tracer_proj_parallel**2, axis=1))
        Tracer_proj_perp = np.sqrt(Tracer_distance**2 - Tracer_proj_parallel_mag.value**2)
        
        Time_array.append(tracer_dict['time'])
        Tracer_parallel.append(Tracer_proj_parallel_mag)
        Tracer_perpendicular.append(Tracer_proj_perp)

    #save the data
    print("saving tracer trajectories")
    file = open('tracer_trajectory.pkl', 'wb')
    pickle.dump((Time_array, Tracer_parallel, Tracer_perpendicular), file)
    file.close()

else:
    print("reading tracer trajectories")
    file = open('tracer_trajectory.pkl', 'rb')
    Time_array, Tracer_parallel, Tracer_perpendicular = pickle.load(file)
    file.close()
    

#Make plots!
import matplotlib.pyplot as plt

plt.clf()
plt.plot(Tracer_parallel, Tracer_perpendicular, alpha=0.25)
plt.scatter(0, 0, marker='o', color='magenta', s=3)
plt.xlim([np.min(Tracer_parallel), 15])
plt.ylim([0, 15])
plt.xlabel('X (AU)')
plt.xlabel('Y (AU)')
plt.savefig("XY_tracer_traj.png")
