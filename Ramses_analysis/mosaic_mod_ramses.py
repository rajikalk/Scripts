#!/usr/bin/env python
import os
import numpy as np
#from pylab import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import csv
import sys
import pickle
import matplotlib.patheffects as path_effects
import matplotlib.gridspec as gridspec
import argparse
import my_ramses_module as mym
import ast

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
args = parse_inputs()
prev_args = args

# Read in directories:
input_file = args.input_file

# Read in input file
positions = []
plot_type = []
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
            plot_type.append(row[2])
            paths.append(row[3])
            if len(row[4]) > 0:
                dict = "{"
                for col in row[4].split(' -'):
                    key_string = col.split(' ')[0][1:]
                    value_string = col.split(' ')[1]
                    dict = dict + "'"+key_string+"':"+value_string+","
                dict = dict[:-1]+"}"
                dict = eval(dict)
                args_dict_all.append(dict)

positions = np.array(positions)
#import pdb
#pdb.set_trace()

mym.set_global_font_size(args.text_font)
frame_no = len(paths)
for fit in range(frame_no):
    file_name = args.save_directory + "mosaic"
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
    
            if plot_type[counter - 2] == 'time_series':
                import pdb
                pdb.set_trace()
                
            if plot_type[counter - 2] == 'movie_frame':
                import pdb
                pdb.set_trace()
