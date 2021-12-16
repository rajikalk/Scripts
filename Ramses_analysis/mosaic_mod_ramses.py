#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import os
from subprocess import call
rank = CW.Get_rank()
size = CW.Get_size()
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv
import sys
import glob
import my_ramses_module as mym
import pickle
import matplotlib.patheffects as path_effects
import yt
import my_ramses_fields as myf
import ast
import matplotlib.gridspec as gridspec
import argparse
import traceback
from PIL import Image

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_file", "--input_file", help="file in input data about the mosaic plot")
    parser.add_argument("-save_dir", "--save_directory", help="do you want define a save directory", type=str)
    parser.add_argument("-sx", "--share_x", help="do you want to share the x axis?", default=False)
    parser.add_argument("-sy", "--share_y", help="do you want to share the y axis?", default=True)
    parser.add_argument("-sa", "--share_ax", help="do you want to share axes?", default=True)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-plot_frame", "--plot_this_frame", help="Do you want to plot a particular frame?", type=int, default=None)
    parser.add_argument("-plot_cbar", "--plot_colourbar", help="do you need to plot a colour bar and waste time playing our with grid limits?", type=str, default='False')
    args = parser.parse_args()
    return args

#=======MAIN=======
#def main():
'''
config_dir = mpl.get_configdir()
if os.path.exists(".Process"+str(rank)) == False:
    os.makedirs(".Process"+str(rank))
os.system('export MPLCONFIGDIR=.Process'+str(rank))
call('cp -r '+config_dir+'/* .Process'+str(rank)+'/.', shell=True)
'''
args = parse_inputs()
prev_args = args

# Read in directories:
input_file = args.input_file

# Read in input file
print("Reading in input mosaic file on rank", rank)
positions = []
paths = []
args_dict_all = []
with open(input_file, 'r') as mosaic_file:
    reader = csv.reader(mosaic_file)
    for row in reader:
        if row[0] == 'Grid_inputs:':
            glr = float(row[1])
            grl = float(row[2])
            glw = float(row[3])
            ghspace = float(row[4])
        elif row[0][0] != '#':
            positions.append((int(row[0]), int(row[1])))
            paths.append(row[2])
            if row[3:][0][0] != ' ':
                row[3:][0] = ' '+row[3:][0]
            dict = "{"
            for col in row[3:][0].split(' -')[1:]:
                key_string = col.split(' ')[0]
                value_string = col.split(' ')[1]
                dict = dict + "'"+key_string+"':"+value_string+","
            dict = dict[:-1]+"}"
            dict = ast.literal_eval(dict)
            args_dict_all.append(dict)

positions = np.array(positions)
#import pdb
#pdb.set_trace()

mym.set_global_font_size(args.text_font)
frame_no = len(glob.glob(paths[0]))
save_bool = np.ones(frame_no)
rit = -1
for fit in range(frame_no):
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        file_name = args.save_directory + "movie_frame_" + ("%06d" % fit)
        make_fig = False
        if os.path.exists(file_name+".jpg") == False:
            make_fig = True
            #print("making frame "+file_name+".jpg on rank "+str(rank))
        elif os.path.getsize(file_name+".jpg")/1024 < 160:
            make_fig = True
            #print("making frame "+file_name+".jpg on rank "+str(rank))
        
        if make_fig:
            plt.clf()
            columns = np.max(positions[:,0])
            rows = np.max(positions[:,1])

            width = float(columns)*(14.5/3.)
            height = float(rows)*(17./4.)
            fig =plt.figure(figsize=(width, height))
            
            if args.plot_colourbar != 'False':
                gs_left = gridspec.GridSpec(rows, columns-1)
                gs_right = gridspec.GridSpec(rows, 1)

                gs_left.update(right=glr, wspace=glw, hspace=ghspace)
                gs_right.update(left=grl, hspace=ghspace)
            else:
                gs_left = gridspec.GridSpec(rows, columns)
                gs_left.update(wspace=glw, hspace=ghspace)
                
            axes_dict = {}
            counter = 1

            for pit in range(len(paths)):
                ax_label = 'ax' + str(counter)
                yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                if positions[pit][0] == 1 and positions[pit][1] == 1:
                    if columns > 1:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[0,0])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    else:
                        if args.plot_colourbar != 'False':
                            axes_dict.update({ax_label:fig.add_subplot(gs_right[0,0])})
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[0,0])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                elif positions[pit][0] != columns:
                    if args.share_x and args.share_y:
                        if yit >= len(axes_dict):
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'])})
                            #print "ADDED SUBPLOT:", counter, "on rank", rank
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                            #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_x:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[it][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_y and positions[pit][0]!=1:
                        yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_y:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    else:
                        axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                else:
                    if args.share_x and args.share_y:
                        yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                        if args.plot_colourbar != 'False':
                            axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0], sharex=axes_dict['ax1'], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_x:
                        if args.plot_colourbar != 'False':
                            axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0], sharex=axes_dict['ax1'])})
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharex=axes_dict['ax1'])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    elif args.share_y:
                        yit = np.where(positions[:,1] == positions[pit][1])[0][0]
                        if args.plot_colourbar != 'False':
                            axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1], sharey=axes_dict[list(axes_dict.keys())[yit]])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
                    else:
                        if args.plot_colourbar != 'False':
                            axes_dict.update({ax_label:fig.add_subplot(gs_right[positions[pit][1]-1,0])})
                        else:
                            axes_dict.update({ax_label:fig.add_subplot(gs_left[positions[pit][1]-1,positions[pit][0]-1])})
                        #print "ADDED SUBPLOT:", counter, "on rank", rank
            
                counter = counter + 1
                #axes_dict[ax_label].set(adjustable='box-forced', aspect='equal')
        
                fs = sorted(glob.glob(paths[pit]))[fit]
                try:
                    file = open(fs, 'rb')
                    X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, velz, part_info, args_dict, simfo = pickle.load(file)
                    file.close()
                    del simfo
                except:
                    print("problem with pickle", fs)
                
                
                """
                if np.round(np.mean(args_dict_pit['xlim'])) == np.round(np.mean(X)):
                    xlim = args_dict_pit['xlim']
                    ylim = args_dict_pit['ylim']
                else:
                    xlim = args_dict_pit['xlim'] + np.mean(X)
                    ylim = args_dict_pit['ylim'] + np.mean(Y)
                    
                has_particles = args_dict_pit['has_particles']
                time_val = args_dict_pit['time_val']
                axes_dict[ax_label].set_xlim(xlim)
                axes_dict[ax_label].set_ylim(ylim)
                rv_channel = [fs.split('_')[-2], fs.split('_')[-1].split('.')[0]]
                
                try:
                    image_file = os.getcwd().split('Channel_maps')[0] + 'Time_Series' + os.getcwd().split('Channel_maps')[1] + '/' + fs.split('/')[-3] + '/' + fs.split('/')[-2] + '.pkl'
                    file = open(image_file, 'rb')
                    stuff = pickle.load(file)
                    file.close()
                    background = stuff[2]
                    cbar_min = args_dict_pit['cbar_min']
                    cbar_max = args_dict_pit['cbar_max']
                    if None in (cbar_min, cbar_max):
                        plot = axes_dict[ax_label].pcolormesh(X, Y, background, cmap=plt.cm.Greys, rasterized=True, alpha=0.5)
                    else:
                        plot = axes_dict[ax_label].pcolormesh(X, Y, background, cmap=plt.cm.Greys, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True, alpha=0.1)
                    #print("plotted image")
                except:
                    print("Couldn't get image data")
                    
                non_nan_inds = np.where(np.isnan(image) == False)
                if len(non_nan_inds[0]) > 0:
                    std = np.std(image[non_nan_inds])
                    max = np.max(image[non_nan_inds])
                    mean = np.mean(image[non_nan_inds])
                    #mean = np.mean(image[non_nan_inds])
                    try:
                        level_b = np.arange(1, max.value, std.value*3)
                        level_r = -1*np.arange(1, max.value, std.value*3)[::-1]
                        CS_b = axes_dict[ax_label].contour(X,Y,np.nan_to_num(image*np.sign(sign_array)), levels=level_b, linewidths=0.5, colors='b', alpha=1.0)
                        #CS_b = ax.contour(X,Y,np.nan_to_num(image), levels=level_b, linewidths=0.5, colors='b', alpha=1.0)
                        CS_r = axes_dict[ax_label].contour(X,Y,np.nan_to_num(image*np.sign(sign_array)), levels=level_r, linewidths=0.5, colors='r', alpha=1.0)
                        #ax.clabel(CS,inline=1)
                        #print("Contours plotted")
                    except:
                        print("Couldn't plot contours")
                plt.gca().set_aspect('equal')
                
                mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)
                
                axes_dict[ax_label].set_xlabel('$x$ (AU)', labelpad=-1, fontsize=args.text_font)
                if positions[pit][0] == 1:
                    axes_dict[ax_label].set_ylabel('$y$ (AU)', labelpad=-20, fontsize=args.text_font)
                
                if positions[pit][0] != 1:
                    yticklabels = axes_dict[ax_label].get_yticklabels()
                    plt.setp(yticklabels, visible=False)
                
                if positions[pit][1] == rows:
                    axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
                    if positions[pit][0] != 1:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels[0], visible=False)
                else:
                    xticklabels = axes_dict[ax_label].get_xticklabels()
                    plt.setp(xticklabels, visible=False)
                
                chanel_str = '[' + str(rv_channel[0]) + ',' + str(rv_channel[1]) + ']'
                axes_dict[ax_label].text(0.70*xlim[1], 0.85*ylim[1], chanel_str)
                
                axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font, direction='in')
                axes_dict[ax_label].tick_params(axis='y', which='major', labelsize=args.text_font, direction='in')
                axes_dict[ax_label].xaxis.set_ticks_position('both')
                axes_dict[ax_label].yaxis.set_ticks_position('both')
                """
                
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

                axes_dict[ax_label].set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
                #axes_dict[ax_label].set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
                axes_dict[ax_label].set_xlim(xlim)
                axes_dict[ax_label].set_ylim(ylim)
                
                if 'cmin' in args_dict_all[pit].keys():
                    cbar_min = args_dict_all[pit]['cmin']
                else:
                    cbar_min = args_dict['cbar_min']
                if 'cmax' in args_dict_all[pit].keys():
                    cbar_max = args_dict_all[pit]['cmax']
                else:
                    cbar_max = args_dict['cbar_max']
                
                if 0.0 in (cbar_min, cbar_max) or len(np.where(np.array([cbar_min, cbar_max]) < 0)[0]) > 0 :
                    plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.bwr, rasterized=True, vmin=cbar_min, vmax=cbar_max)
                else:
                    plot = axes_dict[ax_label].pcolormesh(X, Y, image, cmap=plt.cm.gist_heat, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
                plt.gca().set_aspect('equal')
                
                del image
                
                '''
                if frame_no > 0 or time_val > -1.0:
                    axes_dict[ax_label].streamplot(X, Y, magx, magy, density=4, linewidth=0.25, arrowstyle='-', minlength=0.5)
                else:
                    axes_dict[ax_label].streamplot(X, Y, magx, magy, density=4, linewidth=0.25, minlength=0.5)
                '''
                del X
                del Y
                del magx
                del magy
                
                if positions[pit][0] == columns:
                    cbar = plt.colorbar(plot, pad=0.0)
                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
                
                if 'pvl' in args_dict_all[pit].keys():
                    plot_vel_legend = args_dict_all[pit]['pvl']
                    if 'stdv' in args_dict_all[pit].keys():
                        standard_vel = args_dict_all[pit]['stdv']
                    else:
                        standard_vel = 2
                else:
                    plot_vel_legend = False
                    standard_vel = 2
                
                mym.my_own_quiver_function(axes_dict[ax_label], X_vel, Y_vel, velx, vely, plot_velocity_legend=plot_vel_legend, standard_vel=standard_vel, limits=[xlim, ylim])
                
                del X_vel
                del Y_vel
                del velx
                del vely
                
                if 'title' in args_dict_all[pit].keys():
                    title = args_dict_all[pit]['title']
                else:
                    title = args_dict['title']
                
                title_text = axes_dict[ax_label].text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+2))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                
                if 'at' in args_dict_all[pit].keys():
                    annotate_time = args_dict_all[pit]['at']
                else:
                    annotate_time = False
                    
                if annotate_time:
                    time_string = "$t$="+str(int(args_dict['time_val']))+"yr"
                    time_string_raw = r"{}".format(time_string)
                    time_text = axes_dict[ax_label].text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
                    time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

                if 'apm' in args_dict_all[pit].keys():
                    plot_particle_mass = args_dict_all[pit]['apm']
                else:
                    plot_particle_mass = False
                
                if plot_particle_mass:
                    annotate_field = part_info['particle_mass']
                else:
                    annotate_field = None

                if size > 1:
                    try:
                        mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=annotate_field, particle_tags=part_info['particle_tag'])
                    except:
                        save_bool[fit] = 0
                else:
                    mym.annotate_particles(axes_dict[ax_label], part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=annotate_field, particle_tags=part_info['particle_tag'])
                    
                del part_info
                
                
                if positions[pit][0] == 1:
                    axes_dict[ax_label].set_ylabel('$y$ (AU)', labelpad=-20, fontsize=args.text_font)
                
                if positions[pit][0] != 1:
                    yticklabels = axes_dict[ax_label].get_yticklabels()
                    plt.setp(yticklabels, visible=False)
                if positions[pit][1] == rows:
                    axes_dict[ax_label].tick_params(axis='x', which='major', labelsize=args.text_font)
                    if positions[pit][0] != 1:
                        xticklabels = axes_dict[ax_label].get_xticklabels()
                        plt.setp(xticklabels[0], visible=False)
                plt.tick_params(axis='both', which='major', labelsize=args.text_font)
                axes_dict[ax_label].tick_params(direction='inout', color='white')
                for line in axes_dict[ax_label].xaxis.get_ticklines():
                    line.set_color('white')
                for line in axes_dict[ax_label].yaxis.get_ticklines():
                    line.set_color('white')
                    
                if size > 1:
                    if save_bool[fit] == 1 and pit == (len(paths)-1):
                        try:
                            plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=400)
                            #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                            #if pit == len(paths)-1:
                            print('Created frame no '+ str(fit+1) + ' of '+ str(frame_no) +' of projection on rank', rank, 'at time of', str(int(args_dict['time_val'])), 'to save_dir:', file_name + '.jpg')
                        except:
                            print("couldn't save "+ file_name + ".jpg, Try again later")
                    elif save_bool[fit] == 0 and pit == (len(paths)-1):
                        print("Not saving " + file_name +".jpg because of particles" )
                else:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight', dpi=400)
                    #Convert to jpeg
                    #plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    print('Created frame no '+ str(fit+1) + ' of '+ str(frame_no) +' of projection on rank', rank, 'at time of', str(int(args_dict['time_val'])), 'to save_dir:', file_name + '.jpg')
                    
                    


print("completed making movie frames on rank", rank)
