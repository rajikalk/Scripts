#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
rank = CW.Get_rank()
size = CW.Get_size()
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv
import glob
import pickle
import argparse
import os
import my_ramses_module as mym
import matplotlib.patheffects as path_effects

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_dir", "--input_dir", help="Path to movie pickles")
    parser.add_argument("-in_pickle", "--input_pickle", help="Path to sink pickle")
    parser.add_argument("-save_dir", "--save_directory", help="do you want define a save directory", type=str, default='./')
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    args = parser.parse_args()
    return args

#=======MAIN=======
#def main():
args = parse_inputs()
'''
print("read pickle", args.input_pickle)
file_open = open(args.input_pickle, 'rb')
particle_data, counter, sink_ind, sink_form_time = pickle.load(file_open)
file_open.close()
del counter, sink_ind, sink_form_time, particle_data['mass']
print("finished reading in pickle")
'''
no_frames = np.min([len(glob.glob(args.input_dir + '/XY/movie_frame*pkl')), len(glob.glob(args.input_dir + '/XZ/movie_frame*pkl')), len(glob.glob(args.input_dir + '/YZ/movie_frame*pkl'))])
cmap=plt.cm.gist_heat

fit = -1
rit = -1
while fit < no_frames:
    fit = fit + 1
    if os.path.isfile(args.input_dir+'/XY/movie_frame_' + ("%06d" % fit) +'.pkl') and os.path.isfile(args.input_dir+'/XZ/movie_frame_' + ("%06d" % fit) +'.pkl') and os.path.isfile(args.input_dir+'/YZ/movie_frame_' + ("%06d" % fit) +'.pkl'):
        rit = rit + 1
        if rit == size:
            rit = 0
        if rank == rit:
        
            fig = plt.figure()
            gs = fig.add_gridspec(2, 2, wspace=-0.4, hspace=0)
            (ax1, ax2), (ax3, ax4) = gs.subplots()
            
            #===================YZ proj=====================
        
            yz_pickle = args.input_dir+'/YZ/movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(yz_pickle, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
            file.close()
            del yz_pickle, simfo, velz
            
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
            
            yabel = args_dict['yabel']
            cbar_min = args_dict['cbar_min']
            cbar_max = args_dict['cbar_max']
            del args_dict
            
            ax1.set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
            ax1.set_xlim(xlim)
            ax1.set_ylim(ylim)
            
            plot = ax1.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            del image
            ax1.set_aspect('equal')
            if fit > 0 or time_val > -1.0:
                ax1.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                ax1.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            del X, Y, magx, magy
            mym.my_own_quiver_function(ax1, X_vel, Y_vel, velx, vely, plot_velocity_legend=False, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)
            del X_vel, Y_vel, velx, vely
            mym.annotate_particles(ax1, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
            del part_info
        
            time_string = "$t$="+str(int(time_val))+"yr"
            time_string_raw = r"{}".format(time_string)
            time_text = ax1.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
            time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
            
            xticklabels = ax1.get_xticklabels()
            plt.setp(xticklabels, visible=False)
            
            ax1.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            #===================XZ proj=====================
            
            xz_pickle = args.input_dir+'/XZ/movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(xz_pickle, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
            file.close()
            del xz_pickle, simfo, velz
            
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
            del args_dict
            
            ax2.set_xlim(xlim)
            ax2.set_ylim(ylim)
            
            plot = ax2.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            del image
            ax2.set_aspect('equal')
            if fit > 0 or time_val > -1.0:
                ax2.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                ax2.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            del X, Y, magx, magy
            mym.my_own_quiver_function(ax2, X_vel, Y_vel, velx, vely, plot_velocity_legend=False, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)
            del X_vel, Y_vel, velx, vely
            mym.annotate_particles(ax1, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], particle_tags=part_info['particle_tag'])
            del part_info
            
            xticklabels = ax2.get_xticklabels()
            plt.setp(xticklabels, visible=False)
            
            yticklabels = ax2.get_yticklabels()
            plt.setp(yticklabels, visible=False)
            
            ax2.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            #===================XY proj=====================
            
            xy_pickle = args.input_dir+'/XY/movie_frame_' + ("%06d" % fit) +'.pkl'
            file = open(xy_pickle, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
            file.close()
            del xy_pickle, simfo, velz
            
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
            
            xabel = args_dict['xabel']
            del args_dict
            
            ax4.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
            ax4.set_xlim(xlim)
            ax4.set_ylim(ylim)
           
            plot = ax4.pcolormesh(X, Y, image, cmap=cmap, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, zorder=1)
            del image
            ax4.set_aspect('equal')
            
            if fit > 0 or time_val > -1.0:
                ax4.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5, color='grey', zorder=2)
            else:
                ax4.streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5, zorder=2)
            del X, Y, magx, magy
            mym.my_own_quiver_function(ax4, X_vel, Y_vel, velx, vely, plot_velocity_legend=True, limits=[xlim, ylim], standard_vel=args.standard_vel, Z_val=None)
            del X_vel, Y_vel, velx, vely
            mym.annotate_particles(ax1, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], particle_tags=part_info['particle_tag'])
            del part_info
            
            yticklabels = ax2.get_yticklabels()
            plt.setp(yticklabels, visible=False)
            
            ax4.tick_params(axis='both', direction='in', color='white', top=True, right=True)
            
            fig.subplots_adjust(right=0.95)
            cbar_ax = fig.add_axes([0.951, 0.111, 0.02, 0.767])
            cbar = fig.colorbar(plot, cax=cbar_ax)
            cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=0, size=args.text_font)
            
            plt.savefig("colourbar_test.jpg", format='jpg', bbox_inches='tight', dpi=300)
            
            import pdb
            pdb.set_trace()
            
            #===================Accretion profile=====================
            
            ax3.set_xlabel('Time since formation (yr)')
            ax3.set_ylabel('Accretion Rate (M$_\odot$/yr)')
            ax3.set_xlim([0, particle_data['time'][-1]])
            ax3.set_ylim([np.min(particle_data['mdot']), np.max(particle_data['mdot'])])
            ax3.set_aspect(1.e3)
            
            plot_ind = np.argmin(abs(np.array(particle_data['time']) - time_val))
            ax3.semilogy(particle_data['time'][:plot_ind], np.array(particle_data['mdot']).T[0][:plot_ind])
            ax3.semilogy(particle_data['time'][:plot_ind], np.array(particle_data['mdot']).T[1][:plot_ind])
            ax3.scatter(particle_data['time'][plot_ind], np.array(particle_data['mdot']).T[0][plot_ind], marker='o')
            ax3.scatter(particle_data['time'][plot_ind], np.array(particle_data['mdot']).T[1][plot_ind], marker='o')
            ax3.axhline(y=2*part_info['accretion_rad'], linestyle='--')
            
            
            plt.savefig("movie_frame_" + ("%06d" % fit) + ".jpg", format='jpg', bbox_inches='tight', dpi=300)
            print("Made frame " + "movie_frame_" + ("%06d" % fit) + ".jpg" + " on rank" + str(rank))
            
        

