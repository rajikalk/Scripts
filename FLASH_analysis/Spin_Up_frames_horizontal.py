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
import matplotlib.gridspec as gridspec

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
font_size = 7
mym.set_global_font_size(font_size)


#------------------------------------------------------
spin_values = ['0.20', '0.35']
mach_values = ['0.0', '0.2']

start_times = [[4500, 3750], [6750, 4250]]
end_times = [[5400, 5750], [8000, 5750]]

make_frame = [[False, True], [True, True]]

directory_base = ['/home/kuruwira/fast/Protostellar_spin/Flash_2023/Spin_','/Single/Mach_', '/Lref_9/']
cmap=plt.cm.gist_heat

width = 200
stdvel = 3
n_frames = 5

for spin_val in spin_values:
    for mach_val in mach_values:
        if make_frame[spin_values.index(spin_val)][mach_values.index(mach_val)]:
            plot_it = -1
            
            plt.clf()
            fig = plt.figure(figsize=(two_col_width, 0.6*two_col_width))
            G = gridspec.GridSpec(3, n_frames)
            axes_1 = plt.subplot(G[0, :])
            #fig, axs = plt.subplots(ncols=n_frames, nrows=2, figsize=(two_col_width, 0.4*two_col_width), sharex=True, sharey=True)
            #for ax_it in axs.flatten():
            #    ax_it.set_aspect('equal')
            plt.subplots_adjust(wspace=0.01)
            plt.subplots_adjust(hspace=-0.00)
            
            axes_1.set_title("$\Omega t_{ff}$="+spin_val+", $\mathcal{M}$="+mach_val, y=0.94)
            
            start_t = start_times[spin_values.index(spin_val)][mach_values.index(mach_val)]
            end_t = end_times[spin_values.index(spin_val)][mach_values.index(mach_val)]
            if np.isnan(start_t) or np.isnan(end_t):
                plot_times = []
            else:
                plot_times = np.linspace(start_t, end_t, n_frames-2)
                plot_times = [plot_times[0]-1000] + list(plot_times) + [plot_times[-1]+1000]
                plot_times = plot_times + plot_times
            
            #Plot spin up
            sink_pickle = '/home/kuruwira/fast/Analysis/Sink_evol_pickles/Flash_2023_Spin_'+spin_val+'_Single_Mach_'+mach_val+'_Lref_9.pkl'
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
            plot_inds = []
            for plot_time in plot_times[:n_frames]:
                plot_ind = np.argmin(abs(time.in_units('yr').value - plot_time))
                plot_inds.append(plot_ind)
            axes_1.plot(time[plot_inds[0]:plot_inds[-1]].in_units('yr'), L_tot[plot_inds[0]:plot_inds[-1]])
            axes_1.scatter(time[np.array(plot_inds)].in_units('yr'), L_tot[np.array(plot_inds)])
            axes_1.set_xlim([time[plot_inds[0]].in_units('yr'), time[plot_inds[-1]].in_units('yr')])
            axes_1.set_xlabel('Time ($yr$)', labelpad=-0.2)
            axes_1.set_ylabel('h ($10^{15}m^2/s$)', labelpad=-0.2, fontsize=font_size)
            axes_1.tick_params(axis='x', which='major', direction='in', color='k', top=True)
            axes_1.tick_params(axis='y', which='major', direction='in', color='k', right=True)
            axes_1.xaxis.label.set_color('black')
            axes_1.yaxis.label.set_color('black')
            axes_1.tick_params(axis='both', labelsize=font_size)
            
            for plot_time in plot_times:
                plot_it = plot_it + 1
                if plot_it < n_frames:
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
            
                #ax = axs.flatten()[plot_it]
                ax = plt.subplot(G[int(plot_it/n_frames)+1, np.remainder(plot_it, n_frames)])
                ax.set_aspect('equal')
                file = open(pickle_file, 'rb')
                X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
                file.close()
                
                xlim = [-1*width/2, width/2]
                ylim = [-1*width/2, width/2]
                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                
                if plot_it < n_frames:
                    plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
                else:
                    plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.RdYlBu, vmin=cbar_lims[0], vmax=cbar_lims[1], rasterized=True, zorder=1)
                    #plot = ax.pcolormesh(X_image, Y_image, image, cmap=plt.cm.RdYlBu, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
                ax.set_aspect('equal')
                
                if plot_it == n_frames-1:
                    #Figure out colorbar
                    #fig.subplots_adjust(bottom=0.0)
                    cbar_ax = fig.add_axes([0.90, 0.367, 0.015, 0.257])
                    cbar = fig.colorbar(plot, cax=cbar_ax)
                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", labelpad=-8, rotation=270, size=font_size)
                    cbar_ticks = cbar.ax.yaxis.get_ticklabels()[2].set_visible(False)
                elif plot_it == 2*n_frames-1:
                    #fig.subplots_adjust(bottom=0.05)
                    cbar_ax = fig.add_axes([0.90, 0.11, 0.015, 0.257])
                    cbar = fig.colorbar(plot, cax=cbar_ax)
                    cbar.set_label(r"Magnetic Toomre Q", labelpad=14, rotation=270, size=font_size)
                    cbar_ticks = cbar.ax.yaxis.get_ticklabels()[-1].set_visible(False)
                
                ax.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
                if plot_it == 0:
                    plot_velocity_legend = True
                else:
                    plot_velocity_legend = False
                mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend,limits=[xlim, ylim], Z_val=None, standard_vel=stdvel, width_ceil = 0.4)
                
                if len(part_info['particle_tag']) > 1:
                    if np.min(part_info['particle_form_time'][1:] - part_info['particle_form_time'][:-1]) < 0:
                            sort_inds = np.argsort(part_info['particle_form_time'])
                            part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                            part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                            part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                            part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
                if plot_it >= n_frames:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=3)
                else:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None, particle_tags=part_info['particle_tag'], zorder=7, split_threshold=2)
                
                ax.tick_params(axis='both', which='major', labelsize=font_size)
                for line in ax.xaxis.get_ticklines():
                    line.set_color('white')
                for line in ax.yaxis.get_ticklines():
                    line.set_color('white')
                    
                if plot_it < n_frames:
                    time_string = "$t$="+str(int(time_val))+"yr"
                    time_string_raw = r"{}".format(time_string)
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.06*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=font_size)
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                
                ax.tick_params(axis='x', which='major', direction='in', color='w', top=True)
                ax.tick_params(axis='y', which='major', direction='in', color='w', right=True)
                ax.xaxis.label.set_color('black')
                ax.yaxis.label.set_color('black')
                ax.tick_params(axis='both', labelsize=font_size)
                    
                if plot_it < n_frames:
                    xticklabels = ax.get_xticklabels()
                    plt.setp(xticklabels, visible=False)
                if np.remainder(plot_it, n_frames)==0:
                    ax.set_ylabel('AU', fontsize=font_size, labelpad=-10)
                    if plot_it == n_frames:
                        yticklabels = ax.get_yticklabels()
                        plt.setp(yticklabels[-1], visible=False)
                if plot_it >= n_frames:
                    ax.set_xlabel('AU', fontsize=font_size, labelpad=0)
                    if plot_it > n_frames:
                        xticklabels = ax.get_xticklabels()
                        plt.setp(xticklabels[0], visible=False)
                if np.remainder(plot_it, n_frames)!=0:
                    yticklabels = ax.get_yticklabels()
                    plt.setp(yticklabels, visible=False)
                
                plt.savefig("Spin_"+spin_val+"_Mach_"+mach_val+"_Spin_up_horizontal.pdf", format='pdf', bbox_inches='tight', pad_inches=0.02, dpi=300)
