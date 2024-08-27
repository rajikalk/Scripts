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
    particle_data['time'] = particle_data['time'][::5]
    particle_data['mdot'] = particle_data['mdot'][::5]
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


no_frames = np.min([len(glob.glob(args.input_dir + '/XY/movie_frame*pkl')), len(glob.glob(args.input_dir + '/XZ/movie_frame*pkl')), len(glob.glob(args.input_dir + '/YZ/movie_frame*pkl'))])

cmap=plt.cm.gist_heat
prev_primary_mass = np.nan

fit = -1
rit = -1
while fit < no_frames:
    fit = fit + 1
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        frame_name = args.save_directory + "movie_frame_" + ("%06d" % fit) + ".jpg"
        if os.path.isfile(frame_name) == False and os.path.isfile(args.input_dir+'/XY/movie_frame_' + ("%06d" % fit) +'.pkl') and os.path.isfile(args.input_dir+'/XZ/movie_frame_' + ("%06d" % fit) +'.pkl') and os.path.isfile(args.input_dir+'/YZ/movie_frame_' + ("%06d" % fit) +'.pkl'):
    
            plt.clf()
            fig = plt.figure()
            gs = fig.add_gridspec(2, 2, wspace=-0.46, hspace=0)
            (ax1, ax2), (ax3, ax4) = gs.subplots()
            
            #===================YZ proj=====================
        
            pickle_file = args.input_dir+'/YZ/movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(pickle_file, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
            file.close()
            del pickle_file, simfo, velz, file
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
                del shift_x, shift_y
                gc.collect()
            
            yabel = args_dict['yabel']
            cbar_min = args_dict['cbar_min']
            cbar_max = args_dict['cbar_max']
            del args_dict
            gc.collect()
            
            ax1.set_ylabel(yabel, fontsize=args.text_font, labelpad=-5)#, labelpad=-20)
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim)
            
            plot = ax1.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            del image
            gc.collect()
            ax1.set_aspect('equal')
            if fit > 0 or time_val > -1.0:
                ax1.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                ax1.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            del X, Y, magx, magy
            gc.collect()
            mym.my_own_quiver_function(ax1, X_vel, Y_vel, velx, vely, plot_velocity_legend=False, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)
            del X_vel, Y_vel, velx, vely
            gc.collect()
            mym.annotate_particles(ax1, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
            del part_info
            gc.collect()
        
            time_string = "$t$="+str(int(time_val))+"yr"
            time_string_raw = r"{}".format(time_string)
            time_text = ax1.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
            del xlim, ylim
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            xticklabels = ax1.get_xticklabels()
            plt.setp(xticklabels, visible=False)
            del xticklabels
            gc.collect()
            
            ax1.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            #===================XZ proj=====================
            
            pickle_file = args.input_dir+'/XZ/movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(pickle_file, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
            file.close()
            del pickle_file, simfo, velz, file
            gc.collect()
            
            if len(part_info['particle_tag']) == 2:
                prev_primary_mass = np.max(part_info['particle_mass'])
            elif len(part_info['particle_tag']) == 1:
                part_info['particle_mass'] = np.append(prev_primary_mass, part_info['particle_mass'])
                part_info['particle_position'] = np.array([[0, part_info['particle_position'][0][0]], [0, part_info['particle_position'][1][0]]])
                part_info['particle_tag'] = np.append(44, part_info['particle_tag'])
            
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
                del shift_x, shift_y
                gc.collect()
            del args_dict
            gc.collect()
            
            ax2.set_xlim(xlim)
            ax2.set_ylim(ylim)
            
            plot = ax2.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            del image
            gc.collect()
            ax2.set_aspect('equal')
            if fit > 0 or time_val > -1.0:
                ax2.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                ax2.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            del X, Y, magx, magy
            gc.collect()
            mym.my_own_quiver_function(ax2, X_vel, Y_vel, velx, vely, plot_velocity_legend=False, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)
            del X_vel, Y_vel, velx, vely
            gc.collect()
            mym.annotate_particles(ax2, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], particle_tags=part_info['particle_tag'])
            del part_info, xlim, ylim
            gc.collect()
            
            xticklabels = ax2.get_xticklabels()
            plt.setp(xticklabels, visible=False)
            del xticklabels
            gc.collect()
            
            yticklabels = ax2.get_yticklabels()
            plt.setp(yticklabels, visible=False)
            del yticklabels
            gc.collect()
            
            ax2.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            #===================XY proj=====================
            
            pickle_file = args.input_dir+'/XY/movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(pickle_file, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
            file.close()
            del pickle_file, simfo, velz, file
            gc.collect()
            
            if len(part_info['particle_tag']) == 2:
                prev_primary_mass = np.max(part_info['particle_mass'])
            elif len(part_info['particle_tag']) == 1:
                part_info['particle_mass'] = np.append(prev_primary_mass, part_info['particle_mass'])
                part_info['particle_position'] = np.array([[0, part_info['particle_position'][0][0]], [0, part_info['particle_position'][1][0]]])
                part_info['particle_tag'] = np.append(44, part_info['particle_tag'])
            
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
                del shift_x, shift_y
                gc.collect()
            
            xabel = args_dict['xabel']
            del args_dict
            gc.collect()
            
            ax4.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
            ax4.set_xlim(xlim)
            ax4.set_ylim(ylim)
           
            plot = ax4.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            del image
            gc.collect()
            ax4.set_aspect('equal')
            
            if fit > 0 or time_val > -1.0:
                ax4.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                ax4.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            del X, Y, magx, magy
            gc.collect()
            mym.my_own_quiver_function(ax4, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)
            del X_vel, Y_vel, velx, vely
            gc.collect()
            mym.annotate_particles(ax4, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], particle_tags=part_info['particle_tag'])
            del part_info, xlim, ylim
            gc.collect()
            
            yticklabels = ax4.get_yticklabels()
            plt.setp(yticklabels, visible=False)
            del yticklabels
            gc.collect()
            
            ax4.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            #===================Accretion profile=====================
            
            ax3.set_xlabel('Time since formation (yr)', labelpad=-1)
            ax3.set_ylabel('Accretion Rate (M$_\odot$/yr)', labelpad=-1)
            ax3.set_xlim([0, particle_data['time'][-1]])
            ax3.set_ylim([np.min(particle_data['mdot']), np.max(particle_data['mdot'])])
            ax3.set_aspect(1.e3)
            
            plot_ind = np.argmin(abs(np.array(particle_data['time']) - time_val))
            ax3.semilogy(particle_data['time'][:plot_ind], np.array(particle_data['mdot']).T[0][:plot_ind], color='cyan')
            ax3.semilogy(particle_data['time'][:plot_ind], np.array(particle_data['mdot']).T[1][:plot_ind], color='magenta')
            ax3.scatter(particle_data['time'][plot_ind], np.array(particle_data['mdot']).T[0][plot_ind], marker='o', color='cyan', s=5)
            ax3.scatter(particle_data['time'][plot_ind], np.array(particle_data['mdot']).T[1][plot_ind], marker='o', color='magenta', s=5)
            del time_val, plot_ind
            gc.collect()
            
            fig.subplots_adjust(right=0.95)
            cbar_ax = fig.add_axes([0.825, 0.11, 0.02, 0.77])
            cbar = fig.colorbar(plot, cax=cbar_ax)
            cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=0, size=args.text_font)
            
            ax3.tick_params(axis='both', direction='in', top=True, right=True)
            
            plt.savefig(frame_name, format='jpg', bbox_inches='tight', dpi=300)
            print("Made frame " + "movie_frame_" + ("%06d" % fit) + ".jpg" + " on rank " + str(rank))
            sys.stdout.flush()
