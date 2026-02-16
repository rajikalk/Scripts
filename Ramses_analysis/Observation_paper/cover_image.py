#!/usr/bin/env python
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import pickle
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--field", help="What field to you wish to plot?", default="Density")
    parser.add_argument("-f_unit", "--field_unit", help="What units would you like to plot the field?", default="g/cm**3")
    parser.add_argument("-ax", "--axis", help="Along what axis will the plots be made?", default="xz")
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-pzl", "--plot_z_velocities", help="do you want to plot the z velocity?", type=str, default='False')
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-16')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-pd", "--pickle_dump", help="Do you want to dump the plot sata as a pickle? If true, image won't be plotted", default=False)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=250)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-update_alim", "--update_ax_lim", help="Do you want to update the axes limits by taking away the center position values or not?", type=str, default='False')
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#Get input and output directories
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

#Set some plot variables independant on data files
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
mym.set_global_font_size(args.text_font)

#File files
pickle_files = sorted(glob.glob(input_dir+"time*.pkl"))

for pickle_file in pickle_files:
    frame_no = int(pickle_file.split('_')[-1].split('.')[0])
    file = open(pickle_file, 'rb')
    X, Y, image, magx, magy, magz, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()
    
    time_val = args_dict['time_val']
    
    if args.output_filename != None:
        file_name = args.output_filename
    else:
        file_name = pickle_file.split('.pkl')[0]
    
    print(file_name + ".jpg" + " doesn't exist, so plotting image")
    #import pdb
    #pdb.set_trace()
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
    xabel = args_dict['xabel']
    yabel = args_dict['yabel']
    plt.clf()
    fig, ax = plt.subplots()
    ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
    ax.set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    cmap = plt.cm.cividis
    #cmap = plt.cm.gist_heat
    plot = ax.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
    plt.gca().set_aspect('equal')
    #lw = ((np.log10(abs(magz)) - np.min(np.log10(abs(magz)))))/np.max(((np.log10(abs(magz)) - np.min(np.log10(abs(magz))))))
    lw = (abs(magz))/np.max(abs(magz))
    plt.streamplot(X, Y, magx, magy, density=4, arrowstyle='-', minlength=0.5, color='grey', zorder=2, linewidth=1*lw)
    
    
    xy_mag = np.sqrt(velx**2 + velx**2)
    xy_mag = xy_mag/(np.max(xy_mag))
    cbar = plt.colorbar(plot, pad=0.0)
    mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend="False", limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=xy_mag)
        
    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None, zorder=7)
    
    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=300)#, pad_inches=-0.02)
    plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')#, pad_inches=-0.02)
    print('Created:', file_name + '.jpg')
        
print("completed making movie frames on rank")

