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
    parser.add_argument("-in_pickle", "--input_pickle", help="Path to sink pickle")
    parser.add_argument("-save_dir", "--save_directory", help="do you want define a save directory", type=str, default='./')
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=2.0)
    parser.add_argument("-cmin", "--cbar_min", help="if you don't want to use the limits from the pickles, you can redefine them", type=float, default=None)
    parser.add_argument("-cmax", "--cbar_max", help="if you don't want to use the limits from the pickles, you can redefine them", type=float, default=None)
    parser.add_argument("-end_time", "--end_burst_time", type=float)
    args = parser.parse_args()
    return args

#=======MAIN=======
#def main():
args = parse_inputs()
mym.set_global_font_size(args.text_font)

if rank == 0:
    print("read pickle", args.input_pickle)
    sys.stdout.flush()
    file_open = open(args.input_pickle, 'rb')
    particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
    file_open.close()
    del counter, sink_ind, sink_form_time, particle_data['mass'], particle_data['separation'], particle_data['particle_tag'], args.input_pickle
    gc.collect()
    particle_data['time'] = np.array(particle_data['time'][::5])
    particle_data['mdot'] = np.array(particle_data['mdot'][::5])
    print("finished reading in pickle")
    sys.stdout.flush()
else:
    particle_data = {}

CW.Barrier()
if size > 1:
    print("going to start sending data to other ranks (on rank)", rank)
    sys.stdout.flush()
    particle_data = CW.bcast(particle_data, root=0)
    print("broadcasted particle data on rank", rank)
    sys.stdout.flush()
CW.Barrier()
'''
    print("going to start sending data to other ranks")
    sys.stdout.flush()
    rit = 0
    print("rit =", rit)
    sys.stdout.flush()
    while rit < size:
        rit = rit + 1
        print("rit =", rit)
        sys.stdout.flush()
        CW.send(particle_data, dest=rit, tag=rit)
        print("send particle data to rank", rit)
        sys.stdout.flush()
CW.Barrier()

if size > 1:
    rit = 0
    while rit < size:
        rit = rit + 1
        if rank == rit:
            particle_data = CW.recv(source=0, tag=rit)
            print("particle data received on rank", rank)
            sys.stdout.flush()
'''


no_frames = len(glob.glob(args.input_dir + 'movie_frame*pkl'))
cmap=plt.cm.gist_heat
prev_primary_mass = np.nan

#get start and end time
pickle_file = args.input_dir+'movie_frame_' + ("%06d" % 0) +'.pkl'
file = open(pickle_file, 'rb')
X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
file.close()
time_start = args_dict['time_val']

if args.end_burst_time == None:
    pickle_file = args.input_dir+'movie_frame_' + ("%06d" % (no_frames-1)) +'.pkl'
    file = open(pickle_file, 'rb')
    X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
    file.close()
    time_end = args_dict['time_val']
else:
    time_end = args.end_burst_time


fit = -1
rit = -1
while fit < no_frames:
    fit = fit + 1
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        frame_name = args.save_directory + "movie_frame_" + ("%06d" % fit) + ".jpg"
        if os.path.isfile(frame_name) == False and os.path.isfile(args.input_dir+'movie_frame_' + ("%06d" % fit) +'.pkl'):
    
            fig = plt.figure()
            gs = fig.add_gridspec(1, 2, wspace=0, hspace=0)
            ax1, ax2 = gs.subplots()
            
            #===================YZ proj=====================
        
            pickle_file = args.input_dir+'movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(pickle_file, 'rb')
            X, Y, image, vel_rad, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo, center_vel_rv = pickle.load(file)
            file.close()
            #del pickle_file, simfo, velz, file
            gc.collect()
            
            if len(part_info['particle_tag']) == 2:
                prev_primary_mass = np.max(part_info['particle_mass'])
            elif len(part_info['particle_tag']) == 1:
                part_info['particle_mass'] = np.append(prev_primary_mass, part_info['particle_mass'])
                part_info['particle_position'] = np.array([[0, part_info['particle_position'][0][0]], [0, part_info['particle_position'][1][0]]])
                part_info['particle_tag'] = np.append(44, part_info['particle_tag'])
            
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
                #del shift_x, shift_y
                gc.collect()
            
            yabel = args_dict['yabel']
            if args.cbar_min == None:
                cbar_min = args_dict['cbar_min']
            else:
                cbar_min = args.cbar_min
            
            if args.cbar_max == None:
                cbar_max = args_dict['cbar_max']
            else:
                cbar_max = args.cbar_max
            #del args_dict
            gc.collect()
            
            ax2.set_xlabel(yabel, fontsize=args.text_font, labelpad=-1)#, labelpad=-20)
            ax2.set_xlim(xlim)
            ax2.set_ylim(ylim)
            
            plot = ax2.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            
            #del image
            gc.collect()
            ax2.set_aspect('equal')
            #del X, Y, magx, magy
            gc.collect()
            #del X_vel, Y_vel, velx, vely
            gc.collect()
            mym.annotate_particles(ax2, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
            #del part_info
            gc.collect()
        
            time_string = "$t$="+str(int(time_val))+"yr"
            time_string_raw = r"{}".format(time_string)
            time_text = ax2.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
            #del xlim, ylim
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            yticklabels = ax2.get_yticklabels()
            plt.setp(yticklabels, visible=False)
            #xticklabels = ax2.get_xticklabels()
            #plt.setp(xticklabels[0], visible=False)
            #del xticklabels
            gc.collect()
            
            ax2.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            #===================Accretion profile=====================
            
            ratio =(time_end - time_start)/5
            ax1.set_aspect(ratio)
            ax1.set_xlabel('Time since formation (yr)', labelpad=-1)
            ax1.set_ylabel('Accretion Rate (M$_\odot$/yr)', labelpad=-1)
            ax1.set_xlim([time_start, time_end])
            ax1.set_ylim([1.e-9, 1.e-4])
            
            plot_ind = np.argmin(abs(np.array(particle_data['time']) - time_val))
            ax1.semilogy(particle_data['time'][:plot_ind], particle_data['mdot'].T[0][:plot_ind], color='cyan', linewidth=0.5)
            ax1.semilogy(particle_data['time'][:plot_ind], particle_data['mdot'].T[1][:plot_ind], color='magenta', linewidth=0.5)
            ax1.scatter(particle_data['time'][plot_ind], particle_data['mdot'].T[0][plot_ind], marker='o', color='cyan', s=3)
            ax1.scatter(particle_data['time'][plot_ind], particle_data['mdot'].T[1][plot_ind], marker='o', color='magenta', s=3)
            #del time_val, plot_ind
            gc.collect()
            
            fig.subplots_adjust(right=0.95)
            cbar_ax = fig.add_axes([0.95, 0.22, 0.02, 0.55])
            cbar = fig.colorbar(plot, cax=cbar_ax)
            cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=15, size=args.text_font)
            
            ax1.tick_params(axis='both', direction='in', top=True, right=True)
            
            plt.savefig(frame_name, bbox_inches='tight', dpi=300)
            plt.clf()
            plt.close()
            print("Made frame " + "movie_frame_" + ("%06d" % fit) + ".jpg" + " on rank " + str(rank))
            sys.stdout.flush()
            
print("Finished making frames on rank " + str(rank))
sys.stdout.flush()
