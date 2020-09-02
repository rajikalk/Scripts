#!/usr/bin/env python
import yt
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
    parser.add_argument("-abs", "--absolute_image", help="do you want to get the absolute value of the image field?", default="False")
    parser.add_argument("-o", "--output_filename", help="What will you save your output files as?")
    parser.add_argument("-pvl", "--plot_velocity_legend", help="would you like to annotate the velocity legend?", type=str, default="False")
    parser.add_argument("-wr", "--working_rank", default=0, type=int)
    parser.add_argument("-at", "--annotate_time", help="Would you like to annotate the time that is plotted?", type=str, default="False")
    parser.add_argument("-t", "--title", help="What title would you like the image to have? If left blank it won't show.", default="")
    parser.add_argument("-cmin", "--colourbar_min", help="Input a list with the colour bar ranges", type=str, default='1.e-16')
    parser.add_argument("-cmax", "--colourbar_max", help="Input a list with the colour bar ranges", type=float, default=1.e-14)
    parser.add_argument("-tf", "--text_font", help="What font text do you want to use?", type=int, default=10)
    parser.add_argument("-al", "--ax_lim", help="Want to set the limit of the axis to a nice round number?", type=int, default=250)
    parser.add_argument("-apm", "--annotate_particles_mass", help="Do you want to annotate the particle mass?", default=True)
    parser.add_argument("-stdv", "--standard_vel", help="what is the standard velocity you want to annotate?", type=float, default=5.0)
    parser.add_argument("files", nargs='*')
    args = parser.parse_args()
    return args

#=======MAIN=======

rank = CW.Get_rank()
size = CW.Get_size()

#Get input and output directories
args = parse_inputs()

#Define relevant directories
input_dir = sys.argv[1]
save_dir = sys.argv[2]
if os.path.exists(save_dir) == False:
    os.makedirs(save_dir)



if rank == 0:
    file_int = 0
    for sim_file in usable_files:
        file_int = file_int + 1
        try:
            if usable_files[file_int] == usable_files[file_int-1]:
                os.system('cp '+ save_dir + "movie_frame_" + ("%06d" % frames[file_int-1]) + ".pkl " + save_dir + "movie_frame_" + ("%06d" % frames[file_int]) + ".pkl ")
        except:
            continue
    print('Finished copying pickles that use the same file for the same frame')
    
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
print("Rank =", rank)
for frame_val in range(len(frames)):
    rit = rit + 1
    if rit == size:
        rit = 0
    print("RIT ==", rit, "ON RANK", rank)
    if rank == rit:
        print("creating frame", frames[frame_val], "on rank", rank)
        if args.plot_time != None:
            if args.weight_field == 'None':
                weight_field = None
                pickle_file = save_dir + args.field + "_movie_time_" + (str(args.plot_time)) + "_unweighted.pkl"
            else:
                weight_field = args.weight_field
                #pickle_file = path + args.field + "_movie_time_" + (str(args.plot_time)) + ".pkl"
                pickle_file = save_dir + args.axis + '_' + args.field + '_thickness_' + str(int(args.slice_thickness)) + "_AU_movie_time_" + (str(args.plot_time)) + ".pkl"
        else:
            pickle_file = save_dir + "movie_frame_" + ("%06d" % frames[frame_val]) + ".pkl"
        if os.path.isfile(pickle_file) == False:
            print("NOT USING PICKLE ON RANK", rank)
            time_val = m_times[frame_val]
            print("FILE =", usable_files[frame_val])
            ds = yt.load(usable_files[frame_val], units_override=units_override)
            has_particles = has_sinks(ds)
            if has_particles:
                part_info = mym.get_particle_data(ds, axis=args.axis, sink_id=sink_id)
            else:
                part_info = {}
            if args.ax_lim != None:
                if has_particles and args.image_center != 0:
                    xlim = [-1*args.ax_lim + part_info['particle_position'][0][args.image_center - 1], args.ax_lim + part_info['particle_position'][0][args.image_center - 1]]
                    ylim = [-1*args.ax_lim + part_info['particle_position'][1][args.image_center - 1], args.ax_lim + part_info['particle_position'][1][args.image_center - 1]]
            center_vel = [0.0, 0.0, 0.0]
            if args.image_center != 0 and has_particles:
                original_positions = [X, Y, X_vel, Y_vel]
              
                x_pos = 0#np.round(part_info['particle_position'][0][args.image_center - 1]/cl)*cl
                y_pos = 0#np.round(part_info['particle_position'][1][args.image_center - 1]/cl)*cl
                X = X + x_pos
                Y = Y + y_pos
                X_vel = X_vel + x_pos
                Y_vel = Y_vel + y_pos
              
                part_file = usable_files[frame_val][:-13] + 'part' + usable_files[frame_val][-6:]
                f = h5py.File(part_file, 'r')
                ordered_inds = np.argsort(f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['tag                     '])[0][0]])
                center_vel_x = f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['velx                    '])[0][0]][ordered_inds][args.image_center - 1]
                center_vel_y = f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['vely                    '])[0][0]][ordered_inds][args.image_center - 1]
                center_vel_z = f[list(f.keys())[11]][:,np.where(f[list(f.keys())[5]][:] == ['velz                    '])[0][0]][ordered_inds][args.image_center - 1]
                f.close()
                center_vel = [center_vel_x, center_vel_y, center_vel_z]
                print("CENTER_VEL=", center_vel)

            xabel, yabel, xlim, ylim = image_properties(X, Y, args, simfo)
            #xlim = xlim - center_pos[0]
            #ylim = ylim - center_pos[1]
            if args.axis == 'xy':
                center_vel=center_vel[:2]
            else:
                center_vel=center_vel[::2]
      
        if args.pickle_dump == False:
            print("on rank,", rank, "using pickle_file", pickle_file)
            file = open(pickle_file, 'rb')
            X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo = pickle.load(file)
            #X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, xlim, ylim, has_particles, part_info, simfo, time_val, xabel, yabel = pickle.load(file)
            file.close()
          
            if np.round(np.mean(args_dict['xlim'])) == np.round(np.mean(X)):
                xlim = args_dict['xlim']
                ylim = args_dict['ylim']
            else:
                xlim = args_dict['xlim'] + np.mean(X)
                ylim = args_dict['ylim'] + np.mean(Y)
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
          
            if len(usable_files) > 1:
                if args.output_filename == None:
                    file_name = save_dir + "movie_frame_" + ("%06d" % frames[frame_val])
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
            if frame_val > 0 or time_val > -1.0:
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
                
            '''
            if args.update_ax_lim == 'True':
                new_labels = ax.get_xticks() - np.mean(ax.get_xticks())
                new_labels = new_labels.astype(str).tolist()
                ax.set_xticklabels(new_labels)
                
                new_labels = ax.get_yticks() - np.mean(ax.get_yticks())
                new_labels = new_labels.astype(str).tolist()
                ax.set_yticklabels(new_labels)
            '''

            plt.savefig(file_name + ".eps", format='eps', bbox_inches='tight')
            #Convert to jpeg
            eps_image = Image.open(file_name + ".eps")
            eps_image.load(scale=4)
            eps_image.save(file_name + ".jpg")
      
            #os.system('convert -antialias -quality 100 -density 200 -resize 100% -flatten ' + file_name + '.eps ' + file_name + '.jpg')
            #os.remove(file_name + '.eps')
            print('Created frame', (frames[frame_val]+1), 'of', no_frames, 'on rank', rank, 'at time of', str(time_val), 'to save_dir:', file_name + '.eps')
            #del image
            #del magx
            #del magy
            #del velx
            #del vely
            '''
            if args.image_center != 0 and has_particles:
                X, Y, X_vel, Y_vel = original_positions
            '''

        else:
            print("Creating pickle")
            args_dict = {}
            if args.annotate_time == "True":
                args_dict.update({'annotate_time': '$t$='+str(int(time_val))+'yr'})
            args_dict.update({'field':simfo['field']})
            args_dict.update({'annotate_velocity': args.plot_velocity_legend})
            args_dict.update({'time_val': time_val})
            args_dict.update({'cbar_min': cbar_min})
            args_dict.update({'cbar_max': cbar_max})
            args_dict.update({'title': title})
            args_dict.update({'xabel': xabel})
            args_dict.update({'yabel': yabel})
            args_dict.update({'axlim':args.ax_lim})
            args_dict.update({'xlim':xlim})
            args_dict.update({'ylim':ylim})
            args_dict.update({'has_particles':has_particles})
            print("Built dictionary")
            try:
                pickle_file = save_dir + str(args.axis) +'_' + args_dict['field'].split('_')[0] + "_thickness_" + str(int(args.slice_thickness)) +'_AU_movie_time_' + str(float(args_dict['time_val'])) + '.pkl'
            except:
                pickle_file = save_dir + str(args.axis) +'_' + args_dict['field'][1] + "_thickness_" + str(int(args.slice_thickness)) +'_AU_movie_time_' + str(float(args_dict['time_val'])) + '.pkl'
            print("Got pickle file name:", pickle_file)
            file = open(pickle_file, 'wb')
            print("Opened pickle file")
            #pickle.dump((usable_files[frame_val], X, Y, X_vel, Y_vel, image, velx, vely, part_info, args_dict, simfo, args, magx, magy), file)
            pickle.dump((X, Y, image, magx, magy, X_vel, Y_vel, velx, vely, part_info, args_dict, simfo), file)
            print("Dumped data into pickle")
            file.close()
            print("Created pickle:", pickle_file)
        
sys.stdout.flush()
CW.Barrier()

print("completed making movie frames on rank", rank)

