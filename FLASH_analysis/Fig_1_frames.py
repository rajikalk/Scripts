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


#------------------------------------------------------
spin_values = ['0.20', '0.25', '0.30', '0.35']
mach_values = ['0.0', '0.1', '0.2']
directory_base = ['/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_','/Single/Mach_', '/Lref_9/']
cmap=plt.cm.gist_heat

plt.clf()
fig, axs = plt.subplots(ncols=len(mach_values), nrows=len(spin_values), figsize=(two_col_width, single_col_width*2.5), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = 0

for spin_val in spin_values:
    for mach_val in mach_values:
        pickle_file = 'Spin_'+spin_val+'_Mach_'+mach_val+'.pkl'
        if os.path.exists(pickle_file) == False:
            runline = "python /home/kuruwira/Scripts/FLASH_analysis/movie_script.py /home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_"+spin_val+"/Single/Mach_"+mach_val+"/Lref_9/ ./ -pt 10000 -width 300"
            cmd = ['python', '/home/kuruwira/Scripts/FLASH_analysis/movie_script.py', '/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_'+spin_val+'/Single/Mach_'+mach_val+'/Lref_9/', './', '-pt', '10000', '-width', '300']
            subprocess.Popen(cmd).wait()
            
            os.rename('movie_frame_000000.pkl', pickle_file)
        
        file = open(pickle_file, 'rb')
        X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
        file.close()
        
        if spin_val == '0.35':
            axs.flatten()[plot_it].set_xlabel('AU', labelpad=-1, fontsize=10)
        if mach_val == '0.0':
            axs.flatten()[plot_it].set_ylabel('AU', fontsize=10) #, labelpad=-20
        xlim = [np.min(X_image).value, np.max(X_image).value]
        ylim = [np.min(Y_image).value, np.max(Y_image).value]
        axs.flatten()[plot_it].set_xlim(xlim)
        axs.flatten()[plot_it].set_ylim(ylim)
        
        cbar_lims = [1.e-15, 1.e-13]
        stdvel = 2

        plot = axs.flatten()[plot_it].pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
        plt.gca().set_aspect('equal')
        
        axs.flatten()[plot_it].streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
        mym.my_own_quiver_function(axs.flatten()[plot_it], X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
        mym.annotate_particles(axs.flatten()[plot_it], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7)
        
        plt.tick_params(axis='both', which='major')# labelsize=16)
        for line in axs.flatten()[plot_it].xaxis.get_ticklines():
            line.set_color('white')
        for line in axs.flatten()[plot_it].yaxis.get_ticklines():
            line.set_color('white')
            
        time_string = "$t$="+str(int(time_val))+"yr"
        time_string_raw = r"{}".format(time_string)
        time_text = axs.flatten()[plot_it].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
        time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
