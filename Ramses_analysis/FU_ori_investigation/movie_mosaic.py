#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import sys
rank = CW.Get_rank()
size = CW.Get_size()
print("size =", size)
sys.stdout.flush()
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import glob
import pickle
import argparse
import os
import my_ramses_module as mym
import matplotlib.patheffects as path_effects
import gc

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_dir", "--input_dir", help="Path to movie pickles")
    parser.add_argument("-save_dir", "--save_directory", help="do you want define a save directory", type=str, default='./')
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=2.0)
    parser.add_argument("-cmin", "--cbar_min", help="if you don't want to use the limits from the pickles, you can redefine them", type=float, default=None)
    parser.add_argument("-cmax", "--cbar_max", help="if you don't want to use the limits from the pickles, you can redefine them", type=float, default=None)
    args = parser.parse_args()
    return args

#=======MAIN=======
#def main():
args = parse_inputs()
mym.set_global_font_size(args.text_font)

files = sorted(glob.glob(args.input_dir + 'time_*.pkl'))

cmap=plt.cm.gist_heat
two_col_width = 7.20472 #inches
single_col_width = 3.50394 #inches
page_height = 10.62472 #inches
font_size = 10
plot_velocity_legend = False

plt.clf()
#fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(two_col_width, single_col_width), sharex=True, sharey=True)

fig = plt.figure(figsize=(two_col_width, 0.5*two_col_width))
gs = fig.add_gridspec(2, 4, wspace=0, hspace=0)
axs = gs.subplots()


for file_it in range(len(files)):
    file = open(files[file_it], 'rb')
    X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
    #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
    file.close()
    
    cbar_min = args_dict['cbar_min']
    cbar_max = args_dict['cbar_max']
    
    time_val = args_dict['time_val']
    
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
    
    if file_it > 3:
        axs.flatten()[file_it].set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
    if np.remainder(file_it,4) == 0:
        axs.flatten()[file_it].set_ylabel(yabel, fontsize=args.text_font, labelpad=-20)
    axs.flatten()[file_it].set_xlim(xlim)
    axs.flatten()[file_it].set_ylim(ylim)
    
    if 0.0 in (cbar_min, cbar_max) or len(np.where(np.array([cbar_min, cbar_max]) < 0)[0]) > 0:
        plot = axs.flatten()[file_it].pcolormesh(X, Y, image, cmap=plt.cm.bwr, rasterized=True, vmin=cbar_min, vmax=cbar_max, zorder=1)
    else:
        cmap=plt.cm.gist_heat
        plot = axs.flatten()[file_it].pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
    axs.flatten()[file_it].set_aspect('equal')
    #plt.gca().set_aspect('equal')
    
    axs.flatten()[file_it].streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
    if file_it ==  7:
        plot_velocity_legend = True
    mym.my_own_quiver_function(axs.flatten()[file_it], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)

    mym.annotate_particles(axs.flatten()[file_it], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'], zorder=7, ylabel_scale=0.03)

    time_string = "$t$="+str(int(time_val))+"yr"
    time_string_raw = r"{}".format(time_string)
    time_text = axs.flatten()[file_it].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.05*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
    
    if file_it < 4:
        xticklabels = axs.flatten()[file_it].get_xticklabels()
        plt.setp(xticklabels, visible=False)

    if np.remainder(file_it, 4) > 0:
        yticklabels = axs.flatten()[file_it].get_yticklabels()
        plt.setp(yticklabels, visible=False)
            
    axs.flatten()[file_it].tick_params(axis='both', direction='in', color='white', top=True, right=True)

    if file_it == 7:
        plt.savefig(args.save_directory + 'Mosaic.pdf', bbox_inches='tight', dpi=300)
    else:
        plt.savefig(args.save_directory + 'Mosaic.jpg', bbox_inches='tight', dpi=300)
    print("plotting file it", file_it)
