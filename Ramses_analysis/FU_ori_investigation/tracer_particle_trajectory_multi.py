#!/usr/bin/env python
#In this script we are trying to learn how mass is primarily being accreted
import glob
import pickle
import gc
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
mutation_scale = 15
linewidth = 1

Traj_pickles = ['/home/100/rlk100/rlk/RAMSES/Analysis/Tracer_particle_analysis/Event_2/tracer_trajectory.pkl', '/home/100/rlk100/rlk/RAMSES/Analysis/Tracer_particle_analysis/Event_3/tracer_trajectory.pkl', '/home/100/rlk100/rlk/RAMSES/Analysis/Tracer_particle_analysis/Event_4/tracer_trajectory.pkl', '/home/100/rlk100/rlk/RAMSES/Analysis/Tracer_particle_analysis/Event_5/tracer_trajectory.pkl']

fig, axs = plt.subplots(ncols=2, nrows=2, figsize=(two_col_width, 2*single_col_width), sharex=True, sharey=True)
plt.subplots_adjust(hspace=0.0)
plt.subplots_adjust(wspace=0.0)

for Traj_pickle in Traj_pickles:
    if os.path.isfile(Traj_pickle):
        file = open(Traj_pickle, 'rb')
        Time_array, Tracer_parallel, Tracer_perpendicular = pickle.load(file)
        file.close()
        
        Tracer_parallel = -1 * np.array(Tracer_parallel).T
        Tracer_perpendicular = np.array(Tracer_perpendicular).T

        cmap = mpl.colormaps['plasma']
        N = len(Tracer_parallel[0])-1
        colors = cmap(np.linspace(0, 1, N))
        t_event = Time_array[-1] - Time_array[0]
        Time_array = Time_array - Time_array[0]
        Time_norm = Time_array/t_event
        
        ax = axs.flatten()[Traj_pickles.index(Traj_pickle)]

        for tracer_it in range(len(Tracer_parallel)):
            start_sep = np.sqrt(Tracer_parallel[tracer_it][0]**2 +Tracer_perpendicular[tracer_it][0]**2)
            end_sep = np.sqrt(Tracer_parallel[tracer_it][-1]**2 +Tracer_perpendicular[tracer_it][-1]**2)
            if end_sep < 5:
                for pit in range(1,len(Tracer_parallel[tracer_it])):
                    ax.add_patch(mpatches.FancyArrowPatch((Tracer_parallel[tracer_it][pit-1], Tracer_perpendicular[tracer_it][pit-1]), (Tracer_parallel[tracer_it][pit], Tracer_perpendicular[tracer_it][pit]), color=colors[pit-1], linewidth=0.5, arrowstyle='->', shrinkA=0.0, shrinkB=0.0, alpha=0.5, mutation_scale=5))
        ax.scatter(0, 0, marker='*', color='magenta', s=600)
        circle = mpatches.Circle([0, 0], 0.79, fill=False, edgecolor='magenta')
        arrow = mpatches.FancyArrowPatch((0, 0), (2.5, 0), mutation_scale=mutation_scale, color='magenta', linewidth=linewidth)
        plt.gca().add_patch(circle)
        plt.gca().add_patch(arrow)
        ax.set_xlim([-15, 15])
        ax.set_ylim([-1, 15])
        ax.set_xlabel('Distance$_\parallel$ (AU)')
        ax.set_ylabel('Distance$_\perp$ (AU)')
        #plt.gca().set_aspect('equal')
        norm = mpl.colors.Normalize(vmin=0,vmax=1)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        plt.savefig("XY_tracer_traj.pdf", bbox_inches='tight', pad_inches=0.02)
cax = fig.add_axes([0.905, 0.13, 0.03, 0.72])
cbar = plt.colorbar(sm, cax=cax)
cbar.set_label(r"Time Norm (yr)", rotation=270, labelpad=14)
plt.savefig("XY_tracer_traj.pdf", bbox_inches='tight', pad_inches=0.02)
