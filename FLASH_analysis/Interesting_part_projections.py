#!/usr/bin/env python
import os
import glob
import subprocess
from mpi4py.MPI import COMM_WORLD as CW
import sys
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
import my_flash_module as mym
import numpy as np
import pickle

two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
mym.set_global_font_size(font_size)

#get mpi size and ranks
rank = CW.Get_rank()
size = CW.Get_size()
#size=sys.argv[1]
plot_times = np.arange(3800, 5500, 200)

plt.clf()
fig, axs = plt.subplots(ncols=3, nrows=3, figsize=(two_col_width, two_col_width), sharex=True, sharey=True)
plt.subplots_adjust(wspace=0.0)
plt.subplots_adjust(hspace=0.0)

plot_it = -1
for plot_time in plot_times:
    plot_it = plot_it + 1
    pickle_file = 'time_'+str(plot_time)+'.0.pkl'
    if os.path.exists(pickle_file) == False:
        subprocess.run("python3 /home/kuruwira/Scripts/FLASH_analysis/movie_script.py /hits/fast/set/kuruwira/Protostellar_spin/Flash_2023/Spin_0.20/Single/Mach_0.2/Lref_9/ ./ -width 200 -thickness 100 -cmin 5.e-14 -cmax 5.e-12 -make_pickles True -pt " +str(plot_time) +" -pf 0 -stdv 10 -make_frames False -image_center 1", shell=True)
    
    ax = axs.flatten()[plot_it]
    file = open(pickle_file, 'rb')
    X_image, Y_image, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, time_val = pickle.load(file)
    file.close()

    xlim = [-100, 100]
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    
    cbar_lims = [5.e-14, 5.e-12]
    
    stdvel = 10
    
    cmap=plt.cm.gist_heat
    if np.isnan(cbar_lims[0]):
        plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(), rasterized=True, zorder=1)
    else:
        plot = ax.pcolormesh(X_image, Y_image, image, cmap=cmap, norm=LogNorm(vmin=cbar_lims[0], vmax=cbar_lims[1]), rasterized=True, zorder=1)
    plt.gca().set_aspect('equal')

    
    ax.streamplot(X_image.value, Y_image.value, magx.value, magy.value, density=2, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
    
    if plot_it == 0:
        plot_velocity_legend = True
    else:
        plot_velocity_legend = False
    try:
        mym.my_own_quiver_function(ax, X_vel, Y_vel, velx.value, vely.value, plot_velocity_legend=plot_velocity_legend, limits=[xlim, xlim], Z_val=None, standard_vel=stdvel)
    except:
        mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend, limits=[xlim, xlim], Z_val=None, standard_vel=stdvel)

    if len(part_info.keys())>0:
        if 'particle_form_time' in part_info.keys():
            if len(part_info['particle_form_time']) > 1:
                sort_inds = np.argsort(part_info['particle_form_time'])
                part_info['particle_position'] = part_info['particle_position'].T[sort_inds].T
                part_info['particle_mass'] = part_info['particle_mass'][sort_inds]
                part_info['particle_tag'] = part_info['particle_tag'][sort_inds]
                part_info['particle_form_time'] = part_info['particle_form_time'][sort_inds]
        else:
            print("pickle doesn't have sink formation time")
            os.remove(pickle_file)
                
        primary_ind = np.argmax(part_info['particle_mass'])
        mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, xlim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, split_threshold=7)

    plt.tick_params(axis='both', which='major')# labelsize=16)
    for line in ax.xaxis.get_ticklines():
        line.set_color('white')
    for line in ax.yaxis.get_ticklines():
        line.set_color('white')


    '''
    if args.field == 'dens':
        cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=10)
    elif 'Relative_keplerian_velocity_wrt_primary' in args.field:
        cbar.set_label(r"Relative Keplerian Velocity", rotation=270, labelpad=14, size=10)
    elif 'spec' in args.field:
        cbar.set_label(r"Specific angular momentum (g$\,$cm$^{2}/s$)", rotation=270, labelpad=14, size=10)
    else:
        cbar.set_label(r"Angular momentum (g$\,$cm$^{2}/s$)", rotation=270, labelpad=14, size=10)
    '''
    time_string = "$t$="+str(int(time_val))+"yr"
    time_string_raw = r"{}".format(time_string)
    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (xlim[1]-0.07*(xlim[1]-xlim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=10)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    plt.savefig("Interesting_part.pdf", format='pdf', bbox_inches='tight')

fig.subplots_adjust(right=0.95)
cbar_ax = fig.add_axes([0.951, 0.123, 0.02, 0.744])
cbar = fig.colorbar(plot, cax=cbar_ax)
cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=0, size=font_size)
plt.savefig("Interesting_part.pdf", format='pdf', bbox_inches='tight')

