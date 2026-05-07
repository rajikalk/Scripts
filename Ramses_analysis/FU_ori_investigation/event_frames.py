#!/usr/bin/env python
import glob
import sys
import argparse
import numpy as np
import pickle
import my_ramses_module as mym
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import matplotlib.gridspec as gridspec
import os
import subprocess

#------------------------------------------------------

#Ploting parameters
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 7
mym.set_global_font_size(font_size)


#------------------------------------------------------
time_bounds = [[3800, 4900],[5575, 5700], [6580, 6730], [7295, 7340], [7850, 7900]]
burst_bounds = [[], [[5679, 5689]], [[6625, 6635]], [[7309, 7319], [7327, 7337]], [[7858, 7868]]]
cmap=plt.cm.gist_heat

#Start by loading pickel data and then deleting what we don't need
try:
    sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
except:
    sink_pickle = "../particle_data_L20.pkl"
    print("read pickle", sink_pickle)
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    


width = 30
stdvel = 150
n_frames = 5
make_frame = True
event_it = 2
cbar_lims = [1.e-16, 1.e-13]
plot_dt = (time_bounds[event_it -1][1]-time_bounds[event_it -1][0])/4
plot_times = np.arange(time_bounds[event_it -1][0], time_bounds[event_it -1][1]+plot_dt, plot_dt)


plt.clf()
fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
G = gridspec.GridSpec(2, n_frames, height_ratios=[1, 2])
axes_1 = plt.subplot(G[0, :])
plt.subplots_adjust(wspace=0.01)
plt.subplots_adjust(hspace=-0.2)
            
axes_1.set_title("Suppression event "+str(event_it), y=0.8)
start_ind = np.argmin(abs(particle_data['time']-plot_times[0]))
end_ind = np.argmin(abs(particle_data['time']-plot_times[-1]))
#axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='b', ls='-')
axes_1_twin = axes_1.twinx()
axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5)
            
#Plot accretion and separation. This should be loaded from a pickle

axes_1.set_xlabel('Time ($yr$)', labelpad=-0.2)
axes_1.set_ylabel('Accretion rate (M$_\odot/yr$)', labelpad=-0.2, fontsize=font_size)
axes_1_twin.set_ylabel('Separation (au)')
axes_1.tick_params(axis='x', which='major', direction='in', color='k', top=True)
axes_1.tick_params(axis='y', which='major', direction='in', color='k', right=True)
axes_1.xaxis.label.set_color('black')
axes_1.yaxis.label.set_color('black')
axes_1.tick_params(axis='both', labelsize=font_size)
axes_1.set_xlim([plot_times[0], plot_times[-1]])

plt.savefig("Event_"+str(event_it)+"_mosaic.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)

plot_it = -1
for plot_time in plot_times:
    plot_it = plot_it + 1
    plot_time_ind = np.argmin(abs(particle_data['time'] - plot_time))
    axes_1.scatter(particle_data['time'][plot_time_ind], particle_data['mdot'].T[1][plot_time_ind], color='b', marker='o', s=20)
    axes_1_twin.scatter(particle_data['time'][plot_time_ind], particle_data['separation'][plot_time_ind], marker='o', s=20, color='k', alpha=0.5)
    
    movie_plot_pickle = "time_" + str(plot_time) +".pkl"
    if os.path.isfile(movie_plot_pickle) == False:
        #Make movie frame
        cmd = ['python', '/home/100/rlk100/Scripts/Ramses_analysis/movie_script.py', '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/', './', '-sink', '45', '-pt', str(plot_time), '-at', 'True', '-pvl', 'True',  '-ax', 'xy', '-al', '15', '-tf', '12', '-stdv', str(stdvel), '-thickness', '30', '-use_gas', 'False', '-ic', '1', '-update_alim', 'True', '-frames_only', 'False', '-apm', 'True']
        
        subprocess.Popen(cmd).wait()
    
    #load pickle
    file = open(movie_plot_pickle, 'rb')
    X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
    file.close()

    ax = plt.subplot(G[int(plot_it/n_frames)+1, np.remainder(plot_it, n_frames)])
    ax.set_aspect('equal')
    
    
    centre_ind = np.where(part_info['particle_tag']==45)[0]
    X_image = X_image - part_info['particle_position'][0][centre_ind]
    Y_image = Y_image - part_info['particle_position'][1][centre_ind]
    X_vel = X_vel - part_info['particle_position'][0][centre_ind]
    Y_vel = Y_vel - part_info['particle_position'][1][centre_ind]
    part_info['particle_position'][0] = part_info['particle_position'][0] - part_info['particle_position'][0][centre_ind]
    part_info['particle_position'][1] = part_info['particle_position'][1] - part_info['particle_position'][1][centre_ind]

    xlim = [-1*width/2, width/2]
    ylim = [-1*width/2, width/2]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
                
    if plot_it < n_frames:
        plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)

                
    if plot_it == n_frames-1:
        #Figure out colorbar
        #fig.subplots_adjust(bottom=0.0)
        cbar_ax = fig.add_axes([0.90, 0.3, 0.015, 0.257])
        cbar = fig.colorbar(plot, cax=cbar_ax)
        cbar.set_label(r"Density (g$\,$cm$^{-3}$)", labelpad=-8, rotation=270, size=font_size)
        cbar_ticks = cbar.ax.yaxis.get_ticklabels()[2].set_visible(False)
                
    ax.streamplot(X_image, Y_image, magx, magy, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
    if plot_it == 0:
        plot_velocity_legend = True
    else:
        plot_velocity_legend = False
    #mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel, width_ceil = 0.4)
                

    if len(part_info['particle_tag']) > 1:
        sort_inds = np.argsort(part_info['formation_time'])
        part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
        part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
        part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
        part_info['formation_time'] = part_info['formation_time'][sort_inds]
        part_info['particle_velocity'] =  part_info['particle_velocity'][sort_inds]
    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, annotate_velocity=True, standard_vel=stdvel, width_ceil = 1.0, particle_velocity=part_info['particle_velocity'])
                
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    for line in ax.xaxis.get_ticklines():
        line.set_color('white')
    for line in ax.yaxis.get_ticklines():
        line.set_color('white')
                    
    time_string = "$t$="+str(int(plot_time))+"yr"
    time_string_raw = r"{}".format(time_string)
    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.06*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=font_size)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                
    ax.tick_params(axis='x', which='major', direction='in', color='w', top=True)
    ax.tick_params(axis='y', which='major', direction='in', color='w', right=True)
    ax.xaxis.label.set_color('black')
    ax.yaxis.label.set_color('black')
    ax.tick_params(axis='both', labelsize=font_size)
    ax.set_xlabel('AU', fontsize=font_size, labelpad=-1)
                    
    if np.remainder(plot_it, n_frames)==0:
        ax.set_ylabel('AU', fontsize=font_size, labelpad=-5)
        if plot_it == n_frames:
            yticklabels = ax.get_yticklabels()
            plt.setp(yticklabels[-1], visible=False)
    if np.remainder(plot_it, n_frames)!=0:
        yticklabels = ax.get_yticklabels()
        plt.setp(yticklabels, visible=False)
                
    plt.savefig("Event_"+str(event_it)+"_mosaic.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)
