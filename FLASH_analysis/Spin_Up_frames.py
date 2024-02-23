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

start_times = [[4500, np.nan, 3750], [4500, 4000, 3500], [5250, np.nan, 4250], [6750, 5750, 4250]]
end_times = [[5400, np.nan, 5750], [5450, 4500, 5300], [5750, np.nan, 5500], [8000, 7250, 5750]]

directory_base = ['/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_','/Single/Mach_', '/Lref_9/']
cmap=plt.cm.gist_heat

width = 200
stdvel = 5
n_frames = 5

for spin_val in spin_values:
    for mach_val in mach_values:
        plot_it = -1
        plt.clf()
        fig, axs = plt.subplots(ncols=2, nrows=n_frames, figsize=(0.7*two_col_width, page_height), sharex=True, sharey=True)
        for ax_it in axs.flatten():
            ax_it.set_aspect('equal')
        plt.subplots_adjust(wspace=0.01)
        plt.subplots_adjust(hspace=-0.11)
        
        start_t = start_times[spin_values.index(spin_val)][mach_values.index(mach_val)]
        end_t = end_times[spin_values.index(spin_val)][mach_values.index(mach_val)]
        plot_times = np.linspace(start_t, end_t, n_frames-2)
        plot_times = [plot_times[0]-1000] + list(plot_times) + [plot_times[-1]+1000]
        plot_times = np.array([plot_times, plot_times]).T.flatten()
        
        for plot_time in plot_times:
            plot_it = plot_it + 1
            if np.remainder(plot_it, 2)==0:
                pickle_file = 'Dens_Spin_' + spin_val + '_Mach_' + mach_val + '_time_'+ str(plot_time) +'.pkl'
        
                if os.path.exists(pickle_file) == False:
                    cmd = ['python', '/home/kuruwira/Scripts/FLASH_analysis/movie_script.py', '/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_'+spin_val+'/Single/Mach_'+mach_val+'/Lref_9/', './', '-pt', str(plot_time), '-width', str(width), '-thickness', '200', '-no_quiv', '15', '-stdv', str(stdvel), '-image_center', '1']
                    subprocess.Popen(cmd).wait()
            
                    os.rename('time_'+str(plot_time)+'.pkl', pickle_file)
                
                cbar_lims = [1.e-15, 1.e-13]
            else:
                pickle_file = 'Toomre_Q_Spin_' + spin_val + '_Mach_' + mach_val + '_time_'+ str(plot_time) + '.pkl'
        
                if os.path.exists(pickle_file) == False:
                    cmd = ['python', '/home/kuruwira/Scripts/FLASH_analysis/Toomre_q_frames.py', '/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_'+spin_val+'/Single/Mach_'+mach_val+'/Lref_9/', './', '-pt', str(plot_time), '-width', str(width), '-thickness', '200', '-no_quiv', '15', '-stdv', str(stdvel)]
                    subprocess.Popen(cmd).wait()
            
                    os.rename('time_'+str(plot_time)+'.pkl', pickle_file)
                
                #cbar_lims = [1.e-1, 1.e1]
                cbar_lims = [0, 2]
        
            ax = axs.flatten()[plot_it]
            file = open(pickle_file, 'rb')
            X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
            file.close()
            
            xlim = [-1*width/2, width/2]
            ylim = [-1*width/2, width/2]
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            if np.remainder(plot_it, 2)==0:
                plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
            else:
                plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.RdYlBu, vmin=cbar_lims[0], vmax=cbar_lims[1], rasterized=True, zorder=1)
                #plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.RdYlBu, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
            ax.set_aspect('equal')
            
            if plot_it == len(plot_times)-2:
                #Figure out colorbar
                fig.subplots_adjust(bottom=0.0)
                cbar_ax = fig.add_axes([0.16, 0.05, 0.26, 0.02])
                cbar = fig.colorbar(plot, cax=cbar_ax, orientation='horizontal')
                cbar.set_label(r"Density (g$\,$cm$^{-3}$)", labelpad=0, size=font_size)
            elif plot_it == len(plot_times)-1:
                fig.subplots_adjust(bottom=0.0)
                cbar_ax = fig.add_axes([0.56, 0.05, 0.26, 0.02])
                cbar = fig.colorbar(plot, cax=cbar_ax, orientation='horizontal')
                cbar.set_label(r"Magnetic Toomre Q", labelpad=0, size=font_size)
            
            ax.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            if plot_it == 1:
                plot_velocity_legend = True
            else:
                plot_velocity_legend = False
            mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
            
            if len(part_info['particle_tag']) > 1:
                if np.min(part_info['particle_form_time'][1:] - part_info['particle_form_time'][:-1]) < 0:
                        sort_inds = np.argsort(part_info['particle_form_time'])
                        part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                        part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                        part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                        part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
            if np.remainder(plot_it, 2)==0:
                mym.annotate_particles(axs.flatten()[plot_it], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=4)
            else:
                mym.annotate_particles(axs.flatten()[plot_it], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None, particle_tags=part_info['particle_tag'], zorder=7, split_threshold=4)
            
            ax.tick_params(axis='both', which='major', labelsize=font_size)
            for line in ax.xaxis.get_ticklines():
                line.set_color('white')
            for line in ax.yaxis.get_ticklines():
                line.set_color('white')
                
            if np.remainder(plot_it, 2) == 0:
                time_string = "$t$="+str(int(time_val))+"yr"
                time_string_raw = r"{}".format(time_string)
                time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=font_size)
                time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            ax.tick_params(axis='x', which='major', direction='in', color='w', top=True)
            ax.tick_params(axis='y', which='major', direction='in', color='w', right=True)
            ax.xaxis.label.set_color('black')
            ax.yaxis.label.set_color('black')
            ax.tick_params(axis='both', labelsize=font_size)
                
            xticklabels = ax.get_xticklabels()
            plt.setp(xticklabels, visible=False)
            if np.remainder(plot_it, 2)==0:
                ax.set_ylabel('AU', fontsize=font_size, labelpad=-20)
                if spin_val != '0.20':
                    yticklabels = ax.get_yticklabels()
                    plt.setp(yticklabels[-1], visible=False)
            
            plt.savefig("Spin_"+spin_val+"_Mach_"+mach_val+"_Spin_up.pdf", format='pdf', bbox_inches='tight')
