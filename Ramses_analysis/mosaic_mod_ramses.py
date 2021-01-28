#!/usr/bin/env python
from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD as CW
import h5py
import numpy as np
import numpy.ma as ma
from pylab import *
import matplotlib as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from matplotlib.colors import LogNorm
from subprocess import call
import csv
import sys
import os
from matplotlib import transforms
import glob
import my_module as mym
import pickle
import matplotlib.patheffects as path_effects
import yt
import my_fields as myf
import ast
import matplotlib.gridspec as gridspec
import argparse

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

def define_constants():
    constants = {'year':31557600.0, 'au':1.496e13, 'Msun':1.98841586e+33}
    return constants

def parse_inputs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_file", "--input_file", help="file in input data about the mosaic plot")
    parser.add_argument("-save_dir", "--save_directory", help="do you want define a save directory", type=str)
    parser.add_argument("-sx", "--share_x", help="do you want to share the x axis?", default=False)
    parser.add_argument("-sy", "--share_y", help="do you want to share the y axis?", default=True)
    parser.add_argument("-sa", "--share_ax", help="do you want to share axes?", default=True)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-plot_cbar", "--plot_colourbar", help="do you need to plot a colour bar and waste time playing our with grid limits?", type=str, default='False')
    args = parser.parse_args()
    return args

def get_files(path, args):
    if args.yt_proj == False:
        if args.axis == "xz":
            type = "proj"
        else:
            type = "slice"
        source_directory = sorted(glob.glob(path + 'WIND_' + type + '*'))
    else:
        source_directory = sorted(glob.glob(path + 'WIND_hdf5_plt_cnt*'))
    return source_directory

def has_sinks(file):
    try:
        part_file = file[:-12] + 'part' + file[-5:]
        f = h5py.File(part_file, 'r')
        if f[list(f.keys())[1]][-1][-1] > 0:
            f.close()
            return True
        else:
            f.close()
            return False
    except:
        f = h5py.File(file, 'r')
        if "particlepositions" in list(f.keys()):
            f.close()
            return True
        else:
            f.close()
            return False

def sim_info(path, file, args):
    c = define_constants()
    type = get_files(path, args)[1]
    path_split = path.split('/')
    for p_s in path_split:
        if 'omega' in p_s:
            ang_val = p_s.split('_')[-1]
        if 'lref' in p_s:
            temp_split = p_s.split('_')
            for t_s_it in range(len(temp_split)):
                if 'lref' == temp_split[t_s_it]:
                    lref = temp_split[t_s_it+1]
                if '0.' in temp_split[t_s_it]:
                    den_pert = temp_split[t_s_it]
                else:
                    den_pert = '0.0'
    lref = path.split('lref')[-1].split('/')[0].split('_')[-1]
    if '0.' in  path.split('lref')[0].split('/')[-1]:
        den_pert = path.split('lref')[0].split('/')[-1].split('_')[-2]
    else:
        den_pert = '0.0'
    movie_type = path.split('/')[-2].split('_')[0]

    try:
        f = h5py.File(file, 'r')
        if 'r_accretion' in list(f.keys()):
            racc = f['r_accretion'][0]/yt.units.AU.in_units('cm').value
        else:
            racc = 0.0
        for key in list(f.keys()):
            if args.field in key:
                field = key
        dim = np.shape(f[field])[0]
        if args.zoom_times > 0:
            zoom_cell = int((dim - dim/float(args.zoom_times))/2.)
        else:
            zoom_cell = 0
        xmin = f['minmax_xyz'][0][0]/c['au']
        xmax = f['minmax_xyz'][0][1]/c['au']
        xmin_full = xmin
        cl = (xmax-xmin)/dim
        cell_positions = np.arange(xmin, xmax-1, cl)
        xmin = f['minmax_xyz'][0][0]/c['au'] + zoom_cell*cl
        xmax = f['minmax_xyz'][0][1]/c['au'] - zoom_cell*cl
        if args.axis == "xy":
            ymin = f['minmax_xyz'][1][0]/c['au'] + zoom_cell*cl
            ymax = f['minmax_xyz'][1][1]/c['au'] - zoom_cell*cl
        else:
            ymin = f['minmax_xyz'][2][0]/c['au'] + zoom_cell*cl
            ymax = f['minmax_xyz'][2][1]/c['au'] - zoom_cell*cl
        if args.axis == "xz":
            type = "proj"
        else:
            type = "slice"
        annotate_freq = ((xmax/cl) - (xmin/cl))/31.
    except:
        f = h5py.File(file, 'r')
        if has_sinks(file):
            racc = f[list(f.keys())[18]][109][-1]/yt.units.AU.in_units('cm').value
        else:
            racc = 0.0
        f.close()
        if args.field == 'dens':
            field = ('flash', 'dens')
        else:
            part_file = file[:-12] + 'part' + file[-5:]
            f = yt.load(file, particle_filename=part_file)
            import pdb
            pdb.set_trace()
            field = f.derived_field_list[[x[1] for x in f.derived_field_list].index(args.field)]
        dim = 800
        zoom_cell = 0.0
        if args.ax_lim == None:
            xmin = -1000
            xmax = 1000
            ymin = -1000
            ymax = 1000
        else:
            xmin = -1*args.ax_lim
            xmax = args.ax_lim
            ymin = -1*args.ax_lim
            ymax = args.ax_lim
        cl = (xmax-xmin)/dim
        xmin_full = xmin
        type = "hdf5_plt_cnt"
        annotate_freq = dim/31.
    if args.smooth_cells == None:
        smoothing = annotate_freq/2
    else:
        smoothing = int(args.smooth_cells)
    print("PLOT FIELD IS", field)
    sim_info = {'angular_momentum':ang_val,
                'movie_type':movie_type,
                'field': field,
                'dimension': dim,
                'zoom_cell': zoom_cell,
                'movie_file_type': type,
                'xmin': xmin,
                'xmax': xmax,
                'ymin': ymin,
                'ymax': ymax,
                'cell_length': cl,
                'annotate_freq': annotate_freq,
                'r_accretion': racc,
                'type': type,
                'smoothing': smoothing,
                'refinement_level': lref,
                'den_pert': den_pert,
                'xmin_full':  xmin_full
                }
    f.close()
    return sim_info

def get_image_arrays(f, field, simfo, args, X, Y):
    dim = int(simfo['dimension'])
    image = []
    xpos = int(np.round(np.min(X) - simfo['xmin_full'])/simfo['cell_length'])
    ypos = int(np.round(np.min(Y) - simfo['xmin_full'])/simfo['cell_length'])
    #print "XPOS =", xpos
    for x in range(ypos, ypos+len(X[0])):
        image_val = f[field][x]
        if np.shape(image_val)[0] == 1:
            image_val = image_val.transpose()
        image_val = image_val[xpos:xpos+len(X[0])]
        if simfo['movie_file_type'] == "proj":
            image_val = image_val/(f["minmax_xyz"][1][1]-f["minmax_xyz"][1][0])
        image.append(image_val)
    image = np.array(image)
    if np.shape(image)[-1] == 1:
        image = image[:, :, 0]
    return image


def image_properties(X, Y, args, sim_info):
    if args.axis == "xy":
        ylabel = '$y$ (AU)'
    else:
        ylabel = '$z$ (AU)'
    if args.yt_proj == False:
        xlabel = '$x$ (AU)'
    else:
        xlabel = '$x$ (AU)'
    xlim = [np.min(X), np.max(X)]
    ylim = [np.min(Y), np.max(Y)]
    return xlabel, ylabel, xlim, ylim

def rainbow_text(x,y,ls,lc,**kw):
    """
        Take a list of strings ``ls`` and colors ``lc`` and place them next to each
        other, with text ls[i] being shown in color lc[i].
        
        This example shows how to do both vertical and horizontal text, and will
        pass all keyword arguments to plt.text, so you can set the font size,
        family, etc.
        """
    t = plt.gca().transData
    fig = plt.gcf()
    plt.show()
    
    #horizontal version
    for s,c in zip(ls,lc):
        text = plt.text(x,y," "+s+" ",color=c, transform=t, **kw)
        text.draw(fig.canvas.get_renderer())
        ex = text.get_window_extent()
        t = transforms.offset_copy(text._transform, x=ex.width, units='dots')

#=======MAIN=======
#def main():
rank = CW.Get_rank()
size = CW.Get_size()
args = parse_inputs()
prev_args = args
print("Starting mosaic_mod_script on rank", rank)

# Read in directories:
input_file = args.input_file

# Read in input file
print("Reading in input mosaic file on rank", rank)
positions = []
paths = []
args_dict = []
with open(input_file, 'rU') as mosaic_file:
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
            dict = ""
            for col in row[3:]:
                dict = dict + col
                if col != row[-1]:
                    dict = dict + ','
            dict = ast.literal_eval(dict)
            args_temp = argparse.Namespace(**vars(args))
            for key in list(dict.keys()):
                if key in args:
                    exec("args_temp."+ key + " = " + "str(dict[key])")
            args_dict.append(args_temp)
            del args_temp
            args = prev_args

positions = np.array(positions)

c = define_constants()
mym.set_global_font_size(args.text_font)
frame_no = len(glob.glob(paths[0]))
rit = 0
for fit in range(frame_no):
    if rank == rit:
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
            file = open(fs, 'rb')
            X, Y, image, sign_array, part_info, args_dict_pit, sfo = pickle.load(file)
            file.close()
            
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
            
            file_name = fs.split('channel')[0] + fs.split('channel')[0].split('/')[-2]
            if size > 1:
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    if pit == len(paths)-1:
                        print('Created frame of projection', fs.split('channel')[0].split('/')[-2][-1], 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                except:
                    print("couldn't save for the dviread.py problem. Make frame " + fs.split('channel')[0].split('/')[-2][-1] + " on ipython")
            else:
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    if pit == len(paths)-1:
                        print('Created frame of projection', fs.split('channel')[0].split('/')[-2][-1], 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                except:
                    print("Couldn't save figure " + file_name + ".jpg")
            
            
    rit = rit + 1
    if rit == size:
        rit = 0

print("completed making movie frames on rank", rank)



