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
time_bounds = [[3800, 4900],[5575, 5700], [6580, 6720], [7295, 7365], [7850, 7900]]
burst_bounds = [[], [5675, 5700], [6655, 6720], [7325, 7365], [7860, 7900]]
cbar_lims_all = [[], [1.e-15, 1.e-13], [1.e-15, 1.e-13], [1.e-15, 1.e-13], [1.e-15, 1.e-13]]
cmap=plt.cm.gist_heat

#Start by loading pickel data and then deleting what we don't need
try:
    sink_pickle = "/Users/reggie/Documents/Simulation_analysis/FU_ori_analysis/Particle_data_pickles/particle_data_L20.pkl"
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
except:
    sink_pickle = "/scratch/ek9/rlk100/RAMSES/Analysis/Event_plots/particle_data_L20.pkl"
    print("read pickle", sink_pickle)
    file_open = open(sink_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    print("finished reading in pickle")
    


width = 30
stdvel = 1
n_frames = 5
make_frame = True
event_it = args.event_identifier
cbar_lims = cbar_lims_all[event_it-1]
start_burst = burst_bounds[event_it -1][0]
end_burst = burst_bounds[event_it -1][1]
start_time = time_bounds[event_it -1][0]
end_time = time_bounds[event_it -1][1]
if event_it == 4 and os.getcwd().split('/')[-1] == 'End_7340':
    end_burst = 7340
    end_time = 7340
    
plot_dt = (end_burst-start_burst)/4
plot_times = np.arange(start_burst, end_burst+plot_dt, plot_dt)
units_override = {"length_unit":(4.0,"pc"), "velocity_unit":(0.18, "km/s"), "time_unit":(685706129102738.9, "s"), "mass_unit":(2998,"Msun")}
mym.set_units(units_override)


plt.clf()
fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
G = gridspec.GridSpec(2, n_frames, height_ratios=[1, 2])
axes_1 = plt.subplot(G[0, :])
plt.subplots_adjust(wspace=0.01)
plt.subplots_adjust(hspace=-0.2)
            
axes_1.set_title("Burst event "+str(event_it), y=0.8)
start_ind = np.argmin(abs(particle_data['time']-start_time))
end_ind = np.argmin(abs(particle_data['time']-end_time))
#axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[0][start_ind:end_ind], color='b', ls=':')
lns1 = axes_1.semilogy(particle_data['time'][start_ind:end_ind], particle_data['mdot'].T[1][start_ind:end_ind], color='b', ls='-', label="Accretion rate")
axes_1_twin = axes_1.twinx()
lns2 = axes_1_twin.plot(particle_data['time'][start_ind:end_ind], particle_data['separation'][start_ind:end_ind], ls='--', color='k', alpha=0.5, label="Separation")
lns = lns1+lns2
labs = [l.get_label() for l in lns]
axes_1.legend(lns, labs, loc='upper left')
            
#Plot accretion and separation. This should be loaded from a pickle

axes_1.set_xlabel('Time (yr)', labelpad=-0.2, fontsize=font_size) #($yr$)
axes_1.set_ylabel('Accretion rate (M$_\odot$/yr)', labelpad=-0.2, fontsize=font_size)# (M$_\odot/yr$)
axes_1_twin.set_ylabel('Separation (au)', fontsize=font_size)
axes_1.tick_params(axis='x', which='major', direction='in', color='k', top=True)
axes_1.tick_params(axis='y', which='major', direction='in', color='k', right=True)
axes_1.xaxis.label.set_color('black')
axes_1.yaxis.label.set_color('black')
axes_1.tick_params(axis='both', labelsize=font_size)
axes_1.set_xlim([start_time, end_time])
axes_1.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')

plt.savefig("Event_"+str(event_it)+"_"+args.axis+"_mosaic.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)

plot_it = -1
for plot_time in plot_times:
    plot_it = plot_it + 1
    plot_time_ind = np.argmin(abs(particle_data['time'] - plot_time))
    axes_1.scatter(particle_data['time'][plot_time_ind], particle_data['mdot'].T[1][plot_time_ind], color='b', marker='o', s=20)
    axes_1_twin.scatter(particle_data['time'][plot_time_ind], particle_data['separation'][plot_time_ind], marker='o', s=20, color='k', alpha=0.5)
    
    movie_plot_pickle = "time_" + str(plot_time) +".pkl"
    if os.path.isfile(movie_plot_pickle) == False:
        #Make movie frame
        cmd = ['python', '/home/100/rlk100/Scripts/Ramses_analysis/movie_script.py', '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/', './', '-sink', '45', '-pt', str(plot_time), '-at', 'True', '-pvl', 'True',  '-ax', args.axis, '-al', '15', '-tf', '12', '-stdv', str(stdvel), '-thickness', '30', '-use_gas', 'False', '-ic', '1', '-update_alim', 'True', '-frames_only', 'False', '-apm', 'True', '-res', '1028']
        
        subprocess.Popen(cmd).wait()
    tracer_pickle = "tracer_time_" + str(plot_time) + ".pkl"
    if os.path.isfile(tracer_pickle) == False:
        cmd = ['python', '/home/100/rlk100/Scripts/Ramses_analysis/FU_ori_investigation/tracer_particle_analysis.py', '/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/', './', '-pt', str(plot_time)]
        
        subprocess.Popen(cmd).wait()
        
    files = sorted(glob.glob('/home/100/rlk100/gdata/RAMSES/Zoom-in_CPH_sims/Sink_45/Level_19/Level_20/Event_'+str(event_it)+'/data/*/info*'))
    
    #load pickle
    file = open(movie_plot_pickle, 'rb')
    X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
    file.close()
    
    file = open(tracer_pickle, 'rb')
    tracer_data = pickle.load(file)
    file.close()
    
    #if event_it == 2:
    #    usable_inds = np.where(tracer_data['burst_positions'][1]>-5)[0]
    #else:
    usable_inds = np.arange(len(tracer_data['burst_positions'][1]))

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
        plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)

    if plot_it == n_frames-1:
        cbar_ax = fig.add_axes([0.90, 0.267, 0.015, 0.256])
        cbar = fig.colorbar(plot, cax=cbar_ax)
        cbar.set_label(r"Density (g$\,$cm$^{-3}$)", labelpad=-8, rotation=270, size=font_size)
        cbar_ticks = cbar.ax.yaxis.get_ticklabels()[2].set_visible(False)
                
    if plot_it == 0:
        plot_velocity_legend = True
    else:
        plot_velocity_legend = False
        
    if plot_time == plot_times[-1]:
        legend_text=str(int(3)) + "km$\,$s$^{-1}$"
        xvel = (0.07*(xlim[1] - xlim[0]))
        yvel = 0
        pos_start = [9, 12.5]
        width_val = 0.8
        width_ceil = 0.8
        ax.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='magenta', linewidth=width_val, mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0, alpha=width_val/width_ceil))
        annotate_text = ax.text(10, 10, legend_text, va="center", ha="center", color='w', fontsize=font_size)
        annotate_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    
    ax.scatter(tracer_data['other_positions'][ax_x_ind], tracer_data['other_positions'][ax_y_ind], marker='.', s=1, c='orange', edgecolors=None)
    
    ax.scatter(tracer_data['burst_positions'][ax_x_ind][usable_inds], tracer_data['burst_positions'][ax_y_ind][usable_inds], marker='.', s=1, c='magenta', edgecolors=None)

    mym.my_own_quiver_function(ax, tracer_data['burst_positions'][ax_x_ind][usable_inds].value, tracer_data['burst_positions'][ax_y_ind][usable_inds].value, tracer_data['burst_velocity'][ax_x_ind][usable_inds].in_units('cm/s').value, tracer_data['burst_velocity'][ax_y_ind][usable_inds].in_units('cm/s').value, color='magenta', standard_vel=3, plot_velocity_legend=False, pvl_pos=[10, -10])
    
    if plot_time == plot_times[-1]:
        legend_text=str(int(stdvel)) + "km$\,$s$^{-1}$"
        xvel = (0.07*(xlim[1] - xlim[0]))
        yvel = 0
        pos_start = [9, -12.5]
        width_val = 0.8
        width_ceil = 0.8
        ax.add_patch(mpatches.FancyArrowPatch((pos_start[0], pos_start[1]), (pos_start[0]+xvel, pos_start[1]+yvel), arrowstyle='->', color='w', linewidth=width_val, mutation_scale=10.*width_val, shrinkA=0.0, shrinkB=0.0, alpha=width_val/width_ceil))
        annotate_text = ax.text(10, -10, legend_text, va="center", ha="center", color='w', fontsize=font_size)
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
        
    
    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, annotate_velocity=True, standard_vel=1, width_ceil = 1.0, particle_velocity=part_info['particle_velocity'], part_color=part_color, part_tag_split_length=1)

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
    ax.tick_params(axis='both', labelsize=font_size, labelfontfamily='sans-serif')
                    
    if np.remainder(plot_it, n_frames)==0:
        ax.set_ylabel('$'+args.axis[1]+'$ (AU)', fontsize=font_size, labelpad=-5)
        if plot_it == n_frames:
            yticklabels = ax.get_yticklabels()
            plt.setp(yticklabels[-1], visible=False)
    if np.remainder(plot_it, n_frames)!=0:
        yticklabels = ax.get_yticklabels()
        plt.setp(yticklabels, visible=False)

    ax.set_xlabel('$'+args.axis[0]+'$ (AU)', fontsize=font_size, labelpad=-1)
        
    plt.savefig("Event_"+str(event_it)+"_"+args.axis+"_mosaic.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)
    print('saving figure after plotting time', plot_time)
