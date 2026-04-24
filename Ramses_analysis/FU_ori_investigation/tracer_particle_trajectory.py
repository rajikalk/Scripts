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
import matplotlib.patches as mpatches
import matplotlib as mpl


Tracer_parallel = -1 * np.array(Tracer_parallel).T
Tracer_perpendicular = np.array(Tracer_perpendicular).T

cmap = mpl.colormaps['plasma']
colors = cmap(np.linspace(0, 1, len(Tracer_parallel[0])-1))
plt.clf()
fig1, ax1 = plt.subplots()
for tracer_it in range(len(Tracer_parallel)):
    start_sep = np.sqrt(Tracer_parallel[tracer_it][0]**2 +Tracer_perpendicular[tracer_it][0]**2)
    end_sep = np.sqrt(Tracer_parallel[tracer_it][-1]**2 +Tracer_perpendicular[tracer_it][-1]**2)
    if end_sep < start_sep:
        for pit in range(1,len(Tracer_parallel[tracer_it])):
            ax1.add_patch(mpatches.FancyArrowPatch((Tracer_parallel[tracer_it][pit-1], Tracer_perpendicular[tracer_it][pit-1]), (Tracer_parallel[tracer_it][pit], Tracer_perpendicular[tracer_it][pit]), color=colors[pit-1], linewidth=0.5, arrowstyle='->', shrinkA=0.0, shrinkB=0.0, alpha=0.5, mutation_scale=5))
ax1.scatter(0, 0, marker='*', color='magenta', s=500)
circle = mpatches.Circle([0, 0], 0.79, fill=False, edgecolor='magenta')
arrow = mpatches.FancyArrowPatch((0, 0), (25, 0), mutation_scale=10, color='magenta')
plt.gca().add_patch(circle)
plt.gca().add_patch(arrow)
plt.xlim([np.min(Tracer_parallel), np.max(Tracer_parallel)])
plt.ylim([np.min(Tracer_perpendicular), np.max(Tracer_perpendicular)])
plt.xlabel('Distance$_\parallel$ (AU)')
plt.ylabel('Distance$_\perp$ (AU)')
plt.gca().set_aspect('equal')
plt.savefig("XY_tracer_traj_full.jpg", bbox_inches='tight', dpi=300)

plt.clf()
fig1, ax1 = plt.subplots()
for tracer_it in range(len(Tracer_parallel)):
    start_sep = np.sqrt(Tracer_parallel[tracer_it][0]**2 +Tracer_perpendicular[tracer_it][0]**2)
    end_sep = np.sqrt(Tracer_parallel[tracer_it][-1]**2 +Tracer_perpendicular[tracer_it][-1]**2)
    if end_sep < start_sep:
        for pit in range(1,len(Tracer_parallel[tracer_it])):
            ax1.add_patch(mpatches.FancyArrowPatch((Tracer_parallel[tracer_it][pit-1], Tracer_perpendicular[tracer_it][pit-1]), (Tracer_parallel[tracer_it][pit], Tracer_perpendicular[tracer_it][pit]), color=colors[pit-1], linewidth=0.5, arrowstyle='->', shrinkA=0.0, shrinkB=0.0, alpha=0.5, mutation_scale=5))
ax1.scatter(0, 0, marker='*', color='magenta', s=500)
circle = mpatches.Circle([0, 0], 0.79, fill=False, edgecolor='magenta')
arrow = mpatches.FancyArrowPatch((0, 0), (5, 0), mutation_scale=10, color='magenta')
plt.gca().add_patch(circle)
plt.gca().add_patch(arrow)
plt.xlim([-15, 15])
plt.ylim([0, 15])
plt.xlabel('Distance$_\parallel$ (AU)')
plt.ylabel('Distance$_\perp$ (AU)')
plt.gca().set_aspect('equal')
plt.savefig("XY_tracer_traj.jpg", bbox_inches='tight', dpi=300)


plt.clf()
fig1, ax1 = plt.subplots()
for tracer_it in range(len(Tracer_parallel)):
    start_sep = np.sqrt(Tracer_parallel[tracer_it][0]**2 +Tracer_perpendicular[tracer_it][0]**2)
    end_sep = np.sqrt(Tracer_parallel[tracer_it][-1]**2 +Tracer_perpendicular[tracer_it][-1]**2)
    if end_sep < start_sep:
        for pit in range(1,len(Tracer_parallel[tracer_it])):
            ax1.add_patch(mpatches.FancyArrowPatch((Tracer_parallel[tracer_it][pit-1], Tracer_perpendicular[tracer_it][pit-1]), (Tracer_parallel[tracer_it][pit], Tracer_perpendicular[tracer_it][pit]), color=colors[pit-1], linewidth=0.5, arrowstyle='->', shrinkA=0.0, shrinkB=0.0, alpha=0.5, mutation_scale=5))
ax1.scatter(0, 0, marker='*', color='magenta', s=500)
circle = mpatches.Circle([0, 0], 0.79, fill=False, edgecolor='magenta')
arrow = mpatches.FancyArrowPatch((0, 0), (1, 0), mutation_scale=10, color='magenta')
plt.gca().add_patch(circle)
plt.gca().add_patch(arrow)
plt.xlim([-5, 5])
plt.ylim([0, 5])
plt.xlabel('Distance$_\parallel$ (AU)')
plt.ylabel('Distance$_\perp$ (AU)')
plt.gca().set_aspect('equal')
plt.savefig("XY_tracer_traj_zoom.jpg", bbox_inches='tight', dpi=300)
