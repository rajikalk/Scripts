#!/usr/bin/env python
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
from mpi4py.MPI import COMM_WORLD as CW
import pickle

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-sf", "--start_frame", help="which frame do yu wantto start making at?", type=int, default=0)
    parser.add_argument("-abs", "--absolute_image", help="do you want to get the absolute value of the image field?", default="False")
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-17')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-15)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=250)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-update_alim", "--update_ax_lim", help="Do you want to update the axes limits by taking away the center position values or not?", type=str, default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()

#Get input and output directories
args = parse_inputs()

#process some input parameters
cbar_max = args.colourbar_max
try:
    cbar_min = float(args.colourbar_min)
except:
    cbar_min = float(args.colourbar_min[1:])
title_parts = args.title.split('_')
title = ''
for part in title_parts:
    if part != title_parts[-1]:
        title = title + part + ' '
    else:
        title = title + part

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

files = sorted(glob.glob(input_dir+"*.pkl"))
if args.start_frame != 0:
    files = files[args.start_frame:]

if rank == 0:
    frame_no = 0
    for file_it in range(len(files)):
        frame_no = int(files[file_it].split("_")[-1].split(".")[0])
        try:
            if files[file_it] == files[file_it-1]:
                os.system('cp '+ input_dir + "movie_frame_" + ("%06d" % frame_no) + ".pkl " + input_dir + "movie_frame_" + ("%06d" % frame_no) + ".pkl ")
        except:
            continue
    print('Finished copying pickles that use the same file for the same frame')
    
files = sorted(glob.glob(input_dir+"*.pkl"))
no_frames = len(files) + args.start_frame

sys.stdout.flush()
CW.Barrier()

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.patheffects as path_effects
matplotlib.use('Agg')
from PIL import Image

sys.stdout.flush()
CW.Barrier()

rit = args.working_rank - 1
for file in files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        frame_no = int(file.split("_")[-1].split(".")[0])
        print("on rank,", rank, "using pickle_file", file)
        file = open(file, 'rb')
        X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo = pickle.load(file)
        #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
        file.close()
        
        xlim = args_dict['xlim']
        ylim = args_dict['ylim']
        if np.round(np.mean(args_dict['xlim'])) != np.round(np.mean(X)):
            shift_x = np.mean(X)
            shift_y = np.mean(Y)
            X = X - shift_x
            Y = Y - shift_y
            X_vel = X_vel - shift_x
            Y_vel = Y_vel - shift_y
            part_info['particle_position'][0] = part_info['particle_position'][0] - shift_x
            part_info['particle_position'][1] = part_info['particle_position'][1] - shift_y
            
        has_particles = args_dict['has_particles']
        time_val = args_dict['time_val']
        xabel = args_dict['xabel']
        yabel = args_dict['yabel']
        plt.clf()
        fig, ax = plt.subplots()
        ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
        ax.set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
      
        if len(files) > 1:
            if args.output_filename == None:
                file_name = save_dir + "movie_frame_" + ("%06d" % frame_no)
            else:
                file_name = args.output_filename + "_" + str(int(time_val))
        else:
            if args.output_filename != None:
                file_name = args.output_filename
            else:
                file_name = save_dir + "time_" + str(args.plot_time)
      
        if 0.0 in (cbar_min, cbar_max):
            plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.brg, rasterized=True, vmin=cbar_min, vmax=cbar_max)
        else:
            plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
        plt.gca().set_aspect('equal')
        if frame_no > 0 or time_val > -1.0:
            plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
        else:
            plt.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5)
        cbar = plt.colorbar(plot, pad=0.0)
        mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)

        if has_particles:
            if args.annotate_particles_mass == True:
                mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'],depth_array=part_info['depth_position'])
            else:
                mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None,depth_array=part_info['depth_position'])

        if args.annotate_time == "True":
            time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), '$t$='+str(int(time_val))+'yr', va="center", ha="left", color='w', fontsize=args.text_font)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
          
        title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
        title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

        if 'dens' in simfo['field']:
            cbar.set_label('Density (gcm$^{-3}$)', rotation=270, labelpad=14, size=args.text_font)
        else:
            cbar.set_label(simfo['field'][1] + ' ($' + simfo['unit_string'] + '$)', rotation=270, labelpad=14, size=args.text_font)

        plt.tick_params(axis='both', which='major')# labelsize=16)
        for line in ax.xaxis.get_ticklines():
            line.set_color('white')
        for line in ax.yaxis.get_ticklines():
            line.set_color('white')

        plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight')
        #Convert to jpeg
        eps_image = Image.open(file_name + ".eps")
        eps_image.load(scale=4)
        eps_image.save(file_name + ".jpg")
  
        #os.system('convert -antialias -quality 100 -density 200 -resize 100% -flatten ' + file_name + '.eps ' + file_name + '.jpg')
        #os.remove(file_name + '.eps')
        print('Created frame', (frame_no+1), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.eps')
        #del image
        #del magx
        #del magy
        #del velx
        #del vely
        
sys.stdout.flush()
CW.Barrier()

print("completed making movie frames on rank", rank)

