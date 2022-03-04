#!/usr/bin/env python
import yt
yt.enable_parallelism()
import glob
import numpy as np
import sys
import os
import my_ramses_module as mym
import my_ramses_fields as myf
from mpi4py.MPI import COMM_WORLD as CW
import pickle
import matplotlib.pyplot as plt
import matplotlib as cm
from matplotlib.colors import LogNorm
from matplotlib import transforms
import matplotlib.patheffects as path_effects
from matplotlib import ticker
from scipy.ndimage import gaussian_filter


def parse_inputs():
    import argparse
    parser = argparse.ArgumentParser()

    #plotting parameters
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-vaf", "--velocity_annotation_frequency", help="how many velocity vectors do you want annotated across one side?", type=float, default=31.)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default=None)#'1.e-16')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=None)#1.e-14)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=150)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("-conv", "--convolve", help="Do you want to convolve the image with a gaussian beam?", type=str, default='True')
    
    #Concerning image output
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-sm", "--skip_made", help="do you want to skip frames already made?", type=str, default='True')
    parser.add_argument("files", nargs='*')
    
    #Are you plotting a particular time? or do you want to start from a particular frame?
    parser.add_argument("-sf", "--start_frame", help="initial frame to start with", default=0, type=int)
    parser.add_argument("-pt", "--plot_time", help="If you want to plot one specific time, specify time in years", type=float)
    args = parser.parse_args()
    return args

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()
if rank == 0:
    print("size =", size)

#Get input and output directories
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)

sys.stdout.flush()
CW.Barrier()

#Set some plot variables
cbar_max = args.colourbar_max
if args.colourbar_min == None:
    cbar_min = args.colourbar_min
else:
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

sys.stdout.flush()
CW.Barrier()

#Section to plot figures:
print("Finished generating projection pickles")
if args.plot_time is None:
    pickle_files = sorted(glob.glob(input_dir+"*/*projection*.pkl"))
    start_pickle = save_dir + "movie_frame_" + ("%06d" % args.start_frame) + "/projection_0.pkl"
else:
    pickle_files = sorted(glob.glob(input_dir+"time_" + str(int(args.plot_time)) + "/projection*.pkl"))
    start_pickle = save_dir+"time_" + str(int(args.plot_time)) + "/projection_0.pkl"
start_ind = pickle_files.index(start_pickle)
pickle_files = pickle_files[start_ind:]
#if rank==0:
#    print("pickle_files =", pickle_files)
rit = -1
for pickle_file in pickle_files:
    rit = rit + 1
    if rit == size:
        rit = 0
    if rank == rit:
        if len(m_times) > 1:
            if args.output_filename == None:
                file_name = pickle_file.split('.pkl')[0]
            else:
                file_name = args.output_filename + "_" + str(int(time_val))
        else:
            if args.output_filename != None:
                file_name = args.output_filename
            else:
                file_name = save_dir + "time_" + str(int(args.plot_time)) + "/time_" + str(args.plot_time) + "_proj_" + pickle_file.split('.pkl')[0][-1]
        if os.path.exists(file_name + ".jpg") and os.path.exists(file_name + "_rv.jpg") and args.skip_made == 'True':
            if os.stat(file_name + ".jpg").st_size == 0 or os.stat(file_name + "_rv.jpg").st_size == 0:
                make_frame = True
            else:
                make_frame = False
                print(file_name + ".jpg already exists, so skipping")
        else:
            make_frame = True
        if make_frame:
            proj_number = pickle_file.split('.pkl')[0][-1]
            print("on rank", rank, "using pickle_file", pickle_file)
            file = open(pickle_file, 'rb')
            image, vel_rad, velx_full, vely_full, part_info, image_dict = pickle.load(file)
            file.close()
            
            #Get position array for the image
            xlim = [-1*args.ax_lim, args.ax_lim]
            ylim = [-1*args.ax_lim, args.ax_lim]
            x = np.linspace(xlim[0], xlim[1], len(image))
            y = np.linspace(ylim[0], ylim[1], len(image))
            X, Y = np.meshgrid(x, y)
            annotate_space = (xlim[1] - xlim[0])/args.velocity_annotation_frequency
            x_ind = []
            y_ind = []
            counter = 0
            while counter < args.velocity_annotation_frequency:
                val = annotate_space*counter + annotate_space/2. + simfo['xmin']
                x_ind.append(int(val))
                y_ind.append(int(val))
                counter = counter + 1
            X_vel, Y_vel = np.meshgrid(x_ind, y_ind)
            velx, vely, velz = mym.get_quiver_arrays(0.0, 0.0, X, velx_full, vely_full, center_vel=image_dict['center_vel_image'])
            
            if args.convolve == 'True':
                res = (xlim[1] - xlim[0])/np.shape(vel_rad)[0]
                beam_rad = np.sqrt(12*27)/2.355
                beam_rad_pixel = beam_rad/res
                image = gaussian_filter(image,sigma=beam_rad_pixel)
                vel_rad = gaussian_filter(vel_rad,sigma=beam_rad_pixel)
                
            time_val = image_dict['time']
            xabel = 'Distance from center (AU)'
            yabel = 'Distance from center (AU)'
            
            plt.clf()
            fig, ax = plt.subplots()
            ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
            ax.set_ylabel(yabel, fontsize=args.text_font)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            #plot image
            if None in (cbar_min, cbar_max):
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.magma, rasterized=True)
            elif 0.0 in (cbar_min, cbar_max):
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.magma, rasterized=True, vmin=cbar_min, vmax=cbar_max)
            else:
                plot = ax.pcolormesh(X, Y, image, cmap=plt.cm.magma, norm=LogNorm(vmin=cbar_min, vmax=cbar_max), rasterized=True)
            plt.gca().set_aspect('equal')
            #plot colour bar
            cbar = plt.colorbar(plot, pad=0.0)
            
            #plot velocity field
            mym.my_own_quiver_function(ax, X_vel, Y_vel, velx, vely, plot_velocity_legend=args.plot_velocity_legend, limits=[xlim, ylim], standard_vel=args.standard_vel)


            if args.annotate_particles_mass == True:
                mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
            else:
                mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)

            #Set colourbar label
            if 'Density' in simfo['field']:
                if args.divide_by_proj_thickness == "True":
                    cbar.set_label(r"Density (g$\,$cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
                else:
                    cbar.set_label(r"Column Density (g$\,$cm$^{-2}$)", rotation=270, labelpad=14, size=args.text_font)
            elif 'Number_Density' in simfo['field']:
                if args.divide_by_proj_thickness == 'True':
                    cbar.set_label(r"Number Density (cm$^{-3}$)", rotation=270, labelpad=14, size=args.text_font)
                else:
                    cbar.set_label(r"Column Density (cm$^{-2}$)", rotation=270, labelpad=14, size=args.text_font)
            else:
                label_string = simfo['field'][1] + ' ($' + args.field_unit + '$)'
                cbar.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=args.text_font)

            if len(title) > 0:
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

            plt.tick_params(axis='both', which='major')# labelsize=16)
            for line in ax.xaxis.get_ticklines():
                line.set_color('white')
            for line in ax.yaxis.get_ticklines():
                line.set_color('white')

            if args.annotate_time == "True":
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    time_string = "$t$="+str(int(time_val))+"yr"
                    time_string_raw = r"{}".format(time_string)
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
                    try:
                        plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                        time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                    except:
                        print("Couldn't outline time string")
                except:
                    print("Couldn't plot time string")
            
            if size > 1:
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    print('Created frame of projection', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                except:
                    print("couldn't save for the dviread.py problem. Make frame " + str(proj_number) + " on ipython")
            else:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame of projection', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                
            #================================================================================================================================
            #Created RV Plot:
            plt.clf()
            fig, ax = plt.subplots()
            ax.set_xlabel(xabel, labelpad=-1, fontsize=args.text_font)
            ax.set_ylabel(yabel, fontsize=args.text_font) #, labelpad=-20
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            
            if len(m_times) > 1:
                if args.output_filename == None:
                    file_name = pickle_file.split('.pkl')[0] + "_rv"
                else:
                    file_name = args.output_filename + "_" + str(int(time_val)) + "_rv"
            else:
                if args.output_filename != None:
                    file_name = args.output_filename +"_rv"
                else:
                    file_name = save_dir + "time_" + str(int(args.plot_time)) + "/time_" + str(args.plot_time) + "_proj_" + proj_number +"_rv"
            
            bool_den_array = image>args.density_threshold
            vel_rad = vel_rad*bool_den_array #bool_den_array*np.nan*vel_rad
            vel_rad[vel_rad == 0] = np.nan
            
            
            v_std = np.std(vel_rad/100000)
            v_cbar_min = -1
            v_cbar_max = 1
            plot = ax.pcolormesh(X, Y, vel_rad/100000, cmap='idl06_r', rasterized=True, vmin=v_cbar_min, vmax=v_cbar_max)
            fmt = ticker.LogFormatterSciNotation()
            fmt.create_dummy_axis()
            if args.density_threshold != 0.0:
                exp_min = np.log10(args.density_threshold)
            else:
                exp_min = np.log10(cbar_min)
            exp_max = np.log10(cbar_max)
            n_level = (exp_max-exp_min)*2 + 1
            contour_levels = np.logspace(exp_min, exp_max, int(n_level))
            CS = ax.contour(X,Y,image, locator=plt.LogLocator(), linewidths=0.5, colors='k', levels=contour_levels)
            #'{:.1e}'.format(your_num)
            def func(x):
                s = "%.0g" % x
                if "e" in s:
                    tup = s.split('e')
                    significand = tup[0].rstrip('0').rstrip('.')
                    sign = tup[1][0].replace('+', '')
                    exponent = tup[1][1:].lstrip('0')
                    s = ('%se%s%s' % (significand, sign, exponent)).rstrip('e')
                return s
            #ax.clabel(CS,CS.levels,fmt=func)

            plt.gca().set_aspect('equal')
            cbar = plt.colorbar(plot, pad=0.0)
            
            if has_particles:
                if args.annotate_particles_mass == True:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=part_info['particle_mass'], particle_tags=part_info['particle_tag'])
                else:
                    mym.annotate_particles(ax, part_info['particle_position'], part_info['accretion_rad'], limits=[xlim, ylim], annotate_field=None)

            label_string = 'Radial Velocity ($km/s$)'
            cbar.set_label(r"{}".format(label_string), rotation=270, labelpad=14, size=args.text_font)

            if len(title) > 0:
                title_text = ax.text((np.mean(xlim)), (ylim[1]-0.03*(ylim[1]-ylim[0])), title, va="center", ha="center", color='w', fontsize=(args.text_font+4))
                title_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])

            plt.tick_params(axis='both', which='major')# labelsize=16)
            for line in ax.xaxis.get_ticklines():
                line.set_color('white')
            for line in ax.yaxis.get_ticklines():
                line.set_color('white')

            if args.annotate_time == "True":
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    time_string = "$t$="+str(int(time_val))+"yr"
                    time_string_raw = r"{}".format(time_string)
                    time_text = ax.text((xlim[0]+0.01*(xlim[1]-xlim[0])), (ylim[1]-0.03*(ylim[1]-ylim[0])), time_string_raw, va="center", ha="left", color='w', fontsize=args.text_font)
                    try:
                        plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                        time_text.set_path_effects([path_effects.Stroke(linewidth=3, foreground='black'), path_effects.Normal()])
                    except:
                        print("Couldn't outline time string")
                except:
                    print("Couldn't plot time string")
                        
            if size > 1:
                try:
                    plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                    plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                    print('Created frame of radial velocity', proj_number, 'of 8 on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.jpg')
                except:
                    print("couldn't save for the dviread.py problem. Make frame " + str(proj_number) + " on ipython")
            else:
                plt.savefig(file_name + ".jpg", format='jpg', bbox_inches='tight')
                plt.savefig(file_name + ".pdf", format='pdf', bbox_inches='tight')
                print('Created frame of radial velocity', proj_number, 'of 8 at time of', str(time_val), 'to save_dir:', file_name + '.jpg')

print("Completed making frames on rank", rank)
