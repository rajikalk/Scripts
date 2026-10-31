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
import yt
import my_ramses_fields_short as myf
import matplotlib.patches as mpatches
#-----------------------------------------------------

parser = argparse.ArgumentParser()
parser.add_argument("-event_id", "--event_identifier", default=2, type=int)
parser.add_argument("-ax", "--axis", default='xy', type=str)
parser.add_argument('files', nargs='*')
args = parser.parse_args()


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "mathtext.fontset": "stixsans"  # Force math to use a sans-serif look
})

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

#Ploting parameters
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 9
mym.set_global_font_size(font_size)


#------------------------------------------------------
cbar_lims_all = [5.e-19, 5.e-16]
cmap=plt.cm.gist_heat

width = 1000
stdvel = 1
n_frames = 5
make_frame = True
event_it = args.event_identifier
cbar_lims = cbar_lims_all[event_it-1]
start_burst = burst_bounds[event_it -1][0]
end_burst = burst_bounds[event_it -1][1]
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
highlight_obs_start = burst_bounds[event_it -1][0] - 15
highlight_obs_end = burst_bounds[event_it -1][0] + 90
highlight_burst_start = burst_bounds[event_it -1][0]
highlight_burst_end = burst_bounds[event_it -1][0] + 20
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_burst = 7340
    end_time = 7340
    
plot_dt = (end_burst-start_burst)/4
plot_times = np.arange(start_burst, end_burst+plot_dt, plot_dt)
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)

'''
plt.clf()
fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
G = gridspec.GridSpec(2, n_frames, height_ratios=[1, 2])
axes_1 = plt.subplot(G[0, :])
plt.subplots_adjust(wspace=0.01)
plt.subplots_adjust(hspace=-0.2)
'''

plt.clf()
fig, axs = plt.subplots(ncols=5, nrows=1, figsize=(two_col_width, single_col_width*0.7), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.01)
cbar_ax = fig.add_axes([0.90, 0.275, 0.015, 0.28])

movie_plot_pickles = sorted(glob.glob("movie_frame*.pkl"))

plot_it = -1
for movie_plot_pickle in movie_plot_pickles:
    plot_it = plot_it + 1
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
    
    part_info_pickle = 'part_info_'+str(plot_time)+'.pkl'
    if os.path.isfile(part_info_pickle):
        file = open(part_info_pickle, 'rb')
        part_info = pickle.load(file)
        file.close()
    #if len(part_info['particle_velocity']) != len(part_info['particle_tag']) and os.path.isfile(part_info_pickle):
    #    os.remove(part_info_pickle)
    #    print("REMOVING PART INFO FILE", part_info_pickle)
    
    #define axis inds
    if args.axis == 'xy':
        ax_x_ind = 0
        ax_y_ind = 1
    elif args.axis == 'xz':
        ax_x_ind = 0
        ax_y_ind = 2
    else:
        ax_x_ind = 1
        ax_y_ind = 2
    
    #get relative velocity
    if os.path.isfile(part_info_pickle) == False:
        sink_id = int(part_info['particle_tag'][-1].value)
        ds = yt.load(files[-1], units_override=units_override)
        sink_form_time = ds.r["sink_particle_form_time"][sink_id]
        usable_files = mym.find_files([plot_time], files, sink_form_time, sink_id, verbatim=True)
        ds = yt.load(usable_files[0], units_override=units_override)
        if len(part_info['particle_tag']) == 1:
            part_info['particle_velocity'][0] = [ds.r['gas', 'sink_particle_vel'+args.axis[0]][45] - ds.r['gas', 'sink_particle_vel'+args.axis[0]][44]]
            part_info['particle_velocity'][1] = [ds.r['gas', 'sink_particle_vel'+args.axis[1]][45] - ds.r['gas', 'sink_particle_vel'+args.axis[1]][44]]
        else:
            part_info['particle_velocity'][0] = [ds.r['gas', 'sink_particle_vel'+args.axis[0]][44] - ds.r['gas', 'sink_particle_vel'+args.axis[0]][45], ds.r['gas', 'sink_particle_vel'+args.axis[0]][45] - ds.r['gas', 'sink_particle_vel'+args.axis[0]][44]]
            part_info['particle_velocity'][1] = [ds.r['gas', 'sink_particle_vel'+args.axis[1]][44] - ds.r['gas', 'sink_particle_vel'+args.axis[1]][45], ds.r['gas', 'sink_particle_vel'+args.axis[1]][45] - ds.r['gas', 'sink_particle_vel'+args.axis[1]][44]]
        file = open(part_info_pickle, 'wb')
        pickle.dump((part_info), file)
        file.close()

    xlim = [-1*width/2, width/2]
    ylim = [-1*width/2, width/2]
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
                
    if plot_it < n_frames:
        plot = axs.flatten()[plot_it].pcolormesh(X_image, Y_image, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)

    if plot_it == n_frames-1:
        cbar = fig.colorbar(plot, cax=cbar_ax)
        cbar.set_label(r"Density (g$\,$cm$^{-3}$)", labelpad=-8, rotation=270, size=font_size)
        cbar_ticks = cbar.ax.yaxis.get_ticklabels()[2].set_visible(False)
                
    if plot_it == 0:
        plot_velocity_legend = True
    else:
        plot_velocity_legend = False
    
    if movie_plot_pickle == movie_plot_pickles[-1]:
        legend_text=str(int(stdvel)) + "km$\,$s$^{-1}$"
        xvel = (0.07*(xlim[1] - xlim[0]))
        yvel = 0
        pos_start = [9, -12.5]
        width_val = 0.8
        width_ceil = 0.8
        axs.flatten()[plot_it].add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=width_val, mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0, alpha=width_val/width_ceil))
        annotate_text = axs.flatten()[plot_it].text(10, -10, legend_text, va="center", ha="center", color='w', fontsize=font_size)
        annotate_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
        
    part_color = None
    if len(part_info['particle_tag']) > 1:
        sort_inds = np.argsort(part_info['formation_time'])
        part_info['particle_position'][0] = part_info['particle_position'][0][sort_inds]
        part_info['particle_position'][1] = part_info['particle_position'][1][sort_inds]
        part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
        part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
        part_info['formation_time'] = part_info['formation_time'][sort_inds]
        part_info['particle_velocity'][0] = part_info['particle_velocity'][0][sort_inds]
        part_info['particle_velocity'][1] = part_info['particle_velocity'][1][sort_inds]
        part_color = ['b', 'cyan']
        if np.max(abs(part_info['particle_position'])) > np.max(xlim):
            part_info['particle_tag'] = [part_info['particle_tag'][-1]]
            part_info['particle_position'] = np.array([[part_info['particle_position'][0][-1]], [part_info['particle_position'][1][-1]]])
            part_info['particle_mass'] = np.array([part_info['particle_mass'][-1]])
            part_info['particle_tag'] = np.array([part_info['particle_tag'][-1]])
            part_info['formation_time'] = np.array([part_info['formation_time'][-1]])
            part_info['particle_velocity'] = np.array([[part_info['particle_velocity'][0][-1]], [part_info['particle_velocity'][1][-1]]])
            part_color = [part_color[-1]]
        
    #Get unit velocity:
    #part_info['particle_velocity'] =part_info['particle_velocity']/np.sqrt(np.sum(part_info['particle_velocity']**2, axis=0))[0]
        
    
    mym.annotate_particles(axs.flatten()[plot_it], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, annotate_velocity=True, standard_vel=1, width_ceil = 1.0, particle_velocity=part_info['particle_velocity'], part_color=part_color, part_tag_split_length=1)

    axs.flatten()[plot_it].tick_params(axis='both', which='major', labelsize=font_size)
    for line in axs.flatten()[plot_it].xaxis.get_ticklines():
        line.set_color('white')
    for line in axs.flatten()[plot_it].yaxis.get_ticklines():
        line.set_color('white')
                    
    time_string = "$t$="+str(int(plot_time))+"yr"
    time_string_raw = r"{}".format(time_string)
    time_text = axs.flatten()[plot_it].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.06*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=font_size)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                
    axs.flatten()[plot_it].tick_params(axis='x', which='major', direction='in', color='w', top=True)
    axs.flatten()[plot_it].tick_params(axis='y', which='major', direction='in', color='w', right=True)
    axs.flatten()[plot_it].xaxis.label.set_color('black')
    axs.flatten()[plot_it].yaxis.label.set_color('black')
    axs.flatten()[plot_it].tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
            
    if np.remainder(plot_it, n_frames)==0:
        axs.flatten()[plot_it].set_ylabel('$'+args.axis[1]+'$ (AU)', fontsize=font_size, labelpad=-5)
        if plot_it == n_frames:
            yticklabels = axs.flatten()[plot_it].get_yticklabels()
            plt.setp(yticklabels[-1], visible=False)
    if np.remainder(plot_it, n_frames)!=0:
        yticklabels = axs.flatten()[plot_it].get_yticklabels()
        plt.setp(yticklabels, visible=False)

    axs.flatten()[plot_it].set_xlabel('$'+args.axis[0]+'$ (AU)', fontsize=font_size, labelpad=-1)
        
    plt.savefig("Environment_"+str(event_it)+"_"+args.axis+"_mosaic.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('saving figure after plotting time', plot_time)
