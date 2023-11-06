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
import subprocess
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects

#------------------------------------------------------
#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()

#Ploting parameters
matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.bf'] = 'Arial:bold'
matplotlib.rcParams['mathtext.it'] = 'Arial:italic'
matplotlib.rcParams['mathtext.rm'] = 'Arial'
matplotlib.rcParams['mathtext.sf'] = 'Arial'
matplotlib.rcParams['mathtext.default'] = 'regular'
matplotlib.rcParams['font.sans-serif'] = 'Arial'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['text.latex.preamble'] = [
       r'\usepackage{siunitx}',   # i need upright \micro symbols, but you need...
       r'\sisetup{detect-all}',   # ...this to force siunitx to actually use your fonts
       r'\usepackage{helvet}',    # set the normal font here
       r'\usepackage{sansmath}',  # load up the sansmath so that math -> helvet
       r'\sansmath'               # <- tricky! -- gotta actually tell tex to use!
]

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
mym.set_global_font_size(font_size)


#------------------------------------------------------
spin_values = ['0.20', '0.25', '0.30', '0.35']
mach_values = ['0.0', '0.1', '0.2']
max_time = [[10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000], [10000, 10000, 10000]]
directory_base = ['/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_','/Single/Mach_', '/Lref_9/']
cmap=plt.cm.gist_heat

plt.clf()
fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(two_col_width, 3), sharex=True, sharey=True)
for ax_it in axs.flatten():
    ax_it.set_aspect('equal')
plt.subplots_adjust(wspace=0.01)
plt.subplots_adjust(hspace=-0.11)

plot_it = -1
width = 500
stdvel = 2

for frame_no in range(1,5):
        pickle_file = 'frame_'+str(frame_no)+'.pkl'
        
        plot_it = plot_it +1
        file = open(pickle_file, 'rb')
        X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
        file.close()
        
        xlim = [-1*width/2, width/2]
        ylim = [-1*width/2, width/2]
        axs.flatten()[plot_it].set_xlim(xlim)
        axs.flatten()[plot_it].set_ylim(ylim)
        
        cbar_lims = [5.e-17, 5.e-15]

        plot = axs.flatten()[plot_it].pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
        axs.flatten()[plot_it].set_aspect('equal')
        
        axs.flatten()[plot_it].streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
        if frame_no == 1:
            plot_velocity_legend = True
        else:
            plot_velocity_legend = False
        mym.my_own_quiver_function(axs.flatten()[plot_it], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
        
        if len(part_info['particle_tag']) > 1:
            if np.min(part_info['particle_form_time'][1:] - part_info['particle_form_time'][:-1]) < 0:
                    sort_inds = np.argsort(part_info['particle_form_time'])
                    part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                    part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                    part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                    part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
        mym.annotate_particles(axs.flatten()[plot_it], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=4)
        
        axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size)
        for line in axs.flatten()[plot_it].xaxis.get_ticklines():
            line.set_color('white')
        for line in axs.flatten()[plot_it].yaxis.get_ticklines():
            line.set_color('white')
            
        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = axs.flatten()[plot_it].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=font_size)
        time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        
        axs.flatten()[plot_it].tick_params(axis='x', which='major', direction='in', color='w', top=True)
        axs.flatten()[plot_it].tick_params(axis='y', which='major', direction='in', color='w', right=True)
        axs.flatten()[plot_it].xaxis.label.set_color('black')
        axs.flatten()[plot_it].yaxis.label.set_color('black')
        axs.flatten()[plot_it].tick_params(axis='both', labelsize=font_size)
            
        
        axs.flatten()[plot_it].set_xlabel('AU', labelpad=-1, fontsize=font_size)
        if frame_no == 1:
            axs.flatten()[plot_it].set_ylabel('AU', fontsize=font_size, labelpad=-20)
        else:
            yticklabels = axs.flatten()[plot_it].get_yticklabels()
            plt.setp(yticklabels, visible=False)
        
        plt.savefig("Fig_1_xz.pdf", format='pdf', bbox_inches='tight')
        
fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.951, 0.123, 0.02, 0.744])
cbar = fig.colorbar(plot, cax=cbar_ax)
cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=0, size=font_size)

plt.savefig("Decra.pdf", format='pdf', bbox_inches='tight', dpi=300)
