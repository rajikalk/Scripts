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
        Tracer_proj_parallel =np.sqrt(np.sum((projected_vector(yt.YTArray(np.array(tracer_dict['burst_positions']).T, 'au'), tracer_dict['sink_velocity_vector']))**2, axis=1))
        import pdb
        pdb.set_trace()
        sign = np.dot(Tracer_proj_parallel, tracer_dict['sink_velocity_vector'].T)
        Tracer_proj_perp = np.sqrt(Tracer_distance**2 - Tracer_proj_parallel.value**2)
        
        Time_array.append(tracer_dict['time'])
        Tracer_parallel.append(Tracer_proj_parallel)
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
plt.xlim([-1, 15])
plt.ylim([0, 15])
plt.xlabel('X (AU)')
plt.xlabel('Y (AU)')
plt.savefig("XY_tracer_traj.png")
