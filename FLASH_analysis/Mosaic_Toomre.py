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
cmap=plt.cm.gist_heat
sink_pickle = glob.glob("../*.pkl")[0]
Density_frames = sorted(glob.glob("../Density/movie_frame*.pkl"))
Toomre_Q_frames = sorted(glob.glob("../Toomre_Q/movie_frame*.pkl"))
if 'Strong' in os.getcwd():
    time_bounds = [2500, 7000]
else:
    time_bounds = [5250, 9250]

file = open(sink_pickle, 'rb')
sink_data, line_counter = pickle.load(file)
file.close()

primary_ind = list(sink_data.keys())[0]
mass = yt.YTArray(sink_data[primary_ind]['mass'], 'g')
L_tot = np.sqrt(sink_data[primary_ind]['anglx']**2 + sink_data[primary_ind]['angly']**2 + sink_data[primary_ind]['anglz']**2)
L_tot = yt.YTArray(L_tot/sink_data[primary_ind]['mass'], 'cm**2/s')
L_tot = L_tot.in_units('m**2/s')/1.e15
time = sink_data[primary_ind]['time'] - sink_data[primary_ind]['time'][0]
time = yt.YTArray(time, 's')
bound_inds = [np.argmin(abs(time.in_units('yr').value-time_bounds[0])), np.argmin(abs(time.in_units('yr').value-time_bounds[1]))]

if len(Density_frames) != len(Toomre_Q_frames):
    import pdb
    pdb.set_trace()

width = 200
stdvel = 3

rit = -1
for frame_it in range(len(Density_frames)):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        plot_it = -1
        savename = "movie_frame_" + ("%06d" % frame_it) + ".jpg"
        plt.clf()
        fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(two_col_width, 0.33*two_col_width), sharex=False, sharey=False)
        plt.subplots_adjust(wspace=0)
        plt.subplots_adjust(hspace=0)
        axs.flatten()[0].set_aspect('equal')
        axs.flatten()[1].set_aspect('equal')
        if os.path.exists(savename) == False:
            density_pickle = Density_frames[frame_it]
            
            file = open(density_pickle, 'rb')
            X_image, Y_image, image_dens, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
            file.close()
                    
            xlim = [-1*width/2, width/2]
            ylim = [-1*width/2, width/2]
            axs.flatten()[0].set_xlim(xlim)
            axs.flatten()[0].set_ylim(ylim)
                    
            cbar_lims = [1.e-15, 5.e-13]

            plot = axs.flatten()[0].pcolormesh(X_image, Y_image, image_dens, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
            axs.flatten()[0].set_aspect('equal')
            
            axs.flatten()[0].streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            plot_velocity_legend = False
            try:
                mym.my_own_quiver_function(axs.flatten()[0], X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
            except:
                mym.my_own_quiver_function(axs.flatten()[0], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
                    
                    
            if 'particle_tag' in part_info.keys():
                if len(part_info['particle_tag']) > 1:
                    if np.min(part_info['particle_form_time'][1:] - part_info['particle_form_time'][:-1]) < 0:
                            sort_inds = np.argsort(part_info['particle_form_time'])
                            part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                            part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                            part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                            part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
                mym.annotate_particles(axs.flatten()[0], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=4)
                    
            axs.flatten()[0].tick_params(axis='both', which='major', labelsize=font_size)
            for line in axs.flatten()[0].xaxis.get_ticklines():
                line.set_color('white')
            for line in axs.flatten()[0].yaxis.get_ticklines():
                line.set_color('white')
                
            time_string = "$t$="+str(int(time_val))+"yr"
            time_string_raw = r"{}".format(time_string)
            time_text = axs.flatten()[0].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=font_size)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            axs.flatten()[0].tick_params(axis='x', which='major', direction='in', color='w', top=True)
            axs.flatten()[0].tick_params(axis='y', which='major', direction='in', color='w', right=True)
            axs.flatten()[0].set_ylabel('AU', fontsize=font_size, labelpad=-10)
            axs.flatten()[0].xaxis.label.set_color('black')
            axs.flatten()[0].yaxis.label.set_color('black')
            axs.flatten()[0].tick_params(axis='both', labelsize=font_size)
            
            xticklabels = axs.flatten()[0].get_xticklabels()
            plt.setp(xticklabels, visible=False)
            
            cbar_ax = fig.add_axes([0.125, 0.094, 0.26, 0.015])
            cbar = fig.colorbar(plot, cax=cbar_ax, orientation='horizontal')
            cbar.set_label(r"Density (g$\,$cm$^{-3}$)", labelpad=0, size=font_size)
            
            #=======================================================
            
            toomre_q_pickle = Toomre_Q_frames[frame_it]

            file = open(toomre_q_pickle, 'rb')
            X_image, Y_image, image_tq, magx, magy, X_vel_tq, Y_vel_tq, velx_tq, vely_tq, part_info, time_val_tq = pickle.load(file)
            file.close()
            
            if time_val_tq != time_val:
                import pdb
                pdb.set_trace()
                    
            xlim = [-1*width/2, width/2]
            ylim = [-1*width/2, width/2]
            axs.flatten()[1].set_xlim(xlim)
            axs.flatten()[1].set_ylim(ylim)
                    
            cbar_lims = [0, 2]

            plot = axs.flatten()[1].pcolormesh(X_image, Y_image, image_tq, cmap=plt.cm.RdYlBu, vmin=cbar_lims[0], vmax=cbar_lims[1], rasterized=True, zorder=1)
            axs.flatten()[1].set_aspect('equal')
            
            axs.flatten()[1].streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            plot_velocity_legend = True
            try:
                mym.my_own_quiver_function(axs.flatten()[1], X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
            except:
                mym.my_own_quiver_function(axs.flatten()[1], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel)
                    
            if 'particle_tag' in part_info.keys():
                if len(part_info['particle_tag']) > 1:
                    if np.min(part_info['particle_form_time'][1:] - part_info['particle_form_time'][:-1]) < 0:
                            sort_inds = np.argsort(part_info['particle_form_time'])
                            part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                            part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                            part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                            part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
                mym.annotate_particles(axs.flatten()[1], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None, particle_tags=part_info['particle_tag'], zorder=7, split_threshold=4)
                    
            axs.flatten()[1].tick_params(axis='both', which='major', labelsize=font_size)
            for line in axs.flatten()[1].xaxis.get_ticklines():
                line.set_color('white')
            for line in axs.flatten()[1].yaxis.get_ticklines():
                line.set_color('white')
            
            axs.flatten()[1].tick_params(axis='x', which='major', direction='in', color='w', top=True)
            axs.flatten()[1].tick_params(axis='y', which='major', direction='in', color='w', right=True)
            axs.flatten()[1].xaxis.label.set_color('black')
            axs.flatten()[1].yaxis.label.set_color('black')
            axs.flatten()[1].tick_params(axis='both', labelsize=font_size)
            
            xticklabels = axs.flatten()[1].get_xticklabels()
            plt.setp(xticklabels, visible=False)
            yticklabels = axs.flatten()[1].get_yticklabels()
            plt.setp(yticklabels, visible=False)

            #fig.subplots_adjust(bottom=0.05)
            cbar_ax = fig.add_axes([0.385, 0.094, 0.255, 0.015])
            cbar = fig.colorbar(plot, cax=cbar_ax, orientation='horizontal')
            cbar.set_label(r"Magnetic Toomre Q", labelpad=0, size=font_size)
            
            #=========================================
            axs.flatten()[2].plot(time[bound_inds[0]:bound_inds[1]].in_units('yr'), L_tot[bound_inds[0]:bound_inds[1]].in_units('m**2/s'))
            curr_ind = np.argmin(abs(time.in_units('yr').value-time_val))
            axs.flatten()[2].scatter(time[curr_ind].in_units('yr'), L_tot[curr_ind].in_units('m**2/s'))
            axs.flatten()[2].set_xlabel('Time (yr)', labelpad=-0.2)
            axs.flatten()[2].set_xlim([2500, 7000])
            axs.flatten()[2].yaxis.set_label_position("right")
            axs.flatten()[2].yaxis.tick_right()
            axs.flatten()[2].set_ylabel('h ($10^{15}m^2/s$)', labelpad=-0.2, fontsize=font_size)
                    
            plt.savefig("movie_frame_" + ("%06d" % frame_it) + ".jpg", format='jpg', bbox_inches='tight', dpi=300, pad_inches=0.02)
            print("Made frame " + "movie_frame_" + ("%06d" % frame_it) + ".jpg" + " on rank" + str(rank))

